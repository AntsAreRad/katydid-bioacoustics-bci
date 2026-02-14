# baseline_cooccurrence.R
# Co-occurrence between bird and katydid species:
#   - For each katydid species, test co-occurrence with each bird species
#   - Fisher's exact test on 2x2 presence/absence table across sites
#   - Summary of significant positive and negative associations
#
# Input: integrated_results/ CSV files (presence matrices)
# Output: results/cooccurrence/ CSV and figures

library(tidyverse)
library(ggplot2)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/cooccurrence"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


#' Test pairwise co-occurrence between two groups of species
#'
#' For every pair (one species from group A, one from group B), build a 2x2
#' contingency table of presence/absence across shared sites and run Fisher's
#' exact test. Returns the odds ratio, p-value, and direction of association.
#'
#' @param matrix_a Data frame. Presence/absence matrix for group A (site as
#'   first column).
#' @param matrix_b Data frame. Presence/absence matrix for group B.
#' @param label_a Character. Label for group A (e.g. "katydid").
#' @param label_b Character. Label for group B (e.g. "bird").
#' @param p_adjust_method Character. Method for p.adjust. Default "BH".
#' @return Data frame with columns: species_a, species_b, n_both_present,
#'   n_a_only, n_b_only, n_both_absent, odds_ratio, p_value, p_adjusted,
#'   direction.
#' @export
test_pairwise_cooccurrence <- function(matrix_a, matrix_b,
                                        label_a = "katydid",
                                        label_b = "bird",
                                        p_adjust_method = "BH") {
  # find common sites
  common_sites <- intersect(matrix_a$site, matrix_b$site)
  if (length(common_sites) < 5) {
    warning("Fewer than 5 common sites, co-occurrence test unreliable.")
  }

  ma <- matrix_a %>% filter(site %in% common_sites) %>% arrange(site)
  mb <- matrix_b %>% filter(site %in% common_sites) %>% arrange(site)

  sp_a <- setdiff(colnames(ma), "site")
  sp_b <- setdiff(colnames(mb), "site")

  results <- list()
  idx <- 1

  for (a in sp_a) {
    vec_a <- ma[[a]]
    for (b in sp_b) {
      vec_b <- mb[[b]]

      both_present <- sum(vec_a == 1 & vec_b == 1)
      a_only <- sum(vec_a == 1 & vec_b == 0)
      b_only <- sum(vec_a == 0 & vec_b == 1)
      both_absent <- sum(vec_a == 0 & vec_b == 0)

      tab <- matrix(c(both_present, a_only, b_only, both_absent),
                    nrow = 2, byrow = TRUE)

      ft <- tryCatch(
        fisher.test(tab),
        error = function(e) list(estimate = NA, p.value = NA)
      )

      or <- if (is.numeric(ft$estimate)) ft$estimate else NA
      direction <- case_when(
        is.na(or) ~ "NA",
        or > 1 ~ "positive",
        or < 1 ~ "negative",
        TRUE ~ "none"
      )

      results[[idx]] <- data.frame(
        species_a = a,
        species_b = b,
        group_a = label_a,
        group_b = label_b,
        n_both_present = both_present,
        n_a_only = a_only,
        n_b_only = b_only,
        n_both_absent = both_absent,
        odds_ratio = or,
        p_value = ft$p.value,
        direction = direction,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  out <- bind_rows(results)
  out$p_adjusted <- p.adjust(out$p_value, method = p_adjust_method)
  out %>% arrange(p_adjusted, p_value)
}


#' Summarise co-occurrence results
#'
#' @param cooc_df Data frame from test_pairwise_cooccurrence.
#' @param alpha Numeric. Significance threshold on adjusted p. Default 0.05.
#' @return Data frame summary.
#' @export
summarise_cooccurrence <- function(cooc_df, alpha = 0.05) {
  sig <- cooc_df %>% filter(p_adjusted < alpha)
  data.frame(
    total_pairs_tested = nrow(cooc_df),
    significant_pairs = nrow(sig),
    positive_associations = sum(sig$direction == "positive"),
    negative_associations = sum(sig$direction == "negative"),
    alpha = alpha
  )
}


#' Plot co-occurrence heatmap of significant pairs
#'
#' @param cooc_df Data frame from test_pairwise_cooccurrence.
#' @param alpha Numeric. Significance threshold.
#' @param max_pairs Integer. Max pairs to display.
#' @return ggplot object.
#' @export
plot_cooccurrence_significant <- function(cooc_df, alpha = 0.05,
                                           max_pairs = 40) {
  sig <- cooc_df %>%
    filter(p_adjusted < alpha) %>%
    slice_head(n = max_pairs) %>%
    mutate(pair = paste(species_a, "-", species_b),
           log_or = log2(pmax(odds_ratio, 0.01)))

  if (nrow(sig) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5,
                      label = "No significant co-occurrence pairs found",
                      size = 5) +
             theme_void())
  }

  ggplot(sig, aes(x = reorder(pair, log_or), y = log_or, fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c(positive = "#009E73", negative = "#D55E00")) +
    labs(x = NULL, y = "log2(odds ratio)",
         title = "Significant katydid-bird co-occurrence pairs",
         subtitle = sprintf("Fisher exact test, BH-adjusted p < %.2f", alpha),
         fill = "Association") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 13, face = "bold"))
}


# -- load data --
cat("Loading presence matrices...\n")
katydid_matrix <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                           check.names = FALSE)
bird_matrix <- read.csv(file.path(INPUT_DIR, "bird_presence_matrix.csv"),
                        check.names = FALSE)

# -- pairwise tests --
cat("Testing pairwise co-occurrence (katydid x bird)...\n")
n_katydid <- ncol(katydid_matrix) - 1
n_bird <- ncol(bird_matrix) - 1
cat(sprintf("  %d katydid x %d bird = %d pairs\n",
            n_katydid, n_bird, n_katydid * n_bird))

cooc <- test_pairwise_cooccurrence(katydid_matrix, bird_matrix,
                                    "katydid", "bird")
write.csv(cooc, file.path(OUTPUT_DIR, "katydid_bird_cooccurrence.csv"),
          row.names = FALSE)

# summary
cooc_summary <- summarise_cooccurrence(cooc)
write.csv(cooc_summary, file.path(OUTPUT_DIR, "cooccurrence_summary.csv"),
          row.names = FALSE)

# -- figures --
cat("Creating figures...\n")
sig_pairs <- cooc %>% filter(p_adjusted < 0.05)
if (nrow(sig_pairs) > 0) {
  p1 <- plot_cooccurrence_significant(cooc)
  ggsave(file.path(FIGURE_DIR, "katydid_bird_cooccurrence.jpg"), p1,
         width = 10, height = 8, dpi = 300)
} else {
  cat("  No significant pairs found, skipping figure.\n")
}

# -- summary --
cat("\n-- Co-occurrence summary --\n")
cat(sprintf("  Total pairs tested: %d\n", cooc_summary$total_pairs_tested))
cat(sprintf("  Significant (BH p < 0.05): %d\n", cooc_summary$significant_pairs))
cat(sprintf("    positive: %d, negative: %d\n",
            cooc_summary$positive_associations,
            cooc_summary$negative_associations))
cat("\nCo-occurrence analysis complete.\n")
