# baseline_community_structure.R
# Community structure:
#   - PCoA with Jaccard distance on katydid presence matrix
#   - envfit with vegetation variables (999 permutations)
#   - Katydid-bird richness relationship (GLM Poisson)
#
# Input: integrated_results/ CSV files
# Output: results/community/ CSV and figures

library(tidyverse)
library(ggplot2)
library(vegan)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/community"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

N_PERMUTATIONS <- 999



#' Run PCoA on presence/absence matrix using Jaccard distance
#'
#' @param presence_matrix Data frame with site as first column.
#' @return List with pcoa, dist_matrix, site_scores, eigenvalues, or NULL.
#' @export
run_pcoa_jaccard <- function(presence_matrix) {
  sites <- presence_matrix$site
  sp_cols <- setdiff(colnames(presence_matrix), "site")
  mat <- as.matrix(presence_matrix[, sp_cols])
  rownames(mat) <- sites

  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (nrow(mat) < 3) {
    warning("Fewer than 3 non-empty sites.")
    return(NULL)
  }

  dm <- vegdist(mat, method = "jaccard", binary = TRUE)
  pc <- cmdscale(dm, k = min(nrow(mat) - 1, 5), eig = TRUE)

  scores <- data.frame(site = rownames(mat),
                       Axis1 = pc$points[, 1],
                       Axis2 = pc$points[, 2],
                       stringsAsFactors = FALSE)

  list(pcoa = pc, dist_matrix = dm, site_scores = scores,
       eigenvalues = pc$eig)
}


#' Fit environmental vectors with envfit
#'
#' @param pcoa_result List from run_pcoa_jaccard.
#' @param vegetation_df Data frame with vegetation variables.
#' @param n_perm Integer. Permutations.
#' @param site_col_veg Character. Site column in vegetation data.
#' @return List with envfit result, vectors df, significant df, or NULL.
#' @export
run_envfit <- function(pcoa_result, vegetation_df, n_perm = 999,
                       site_col_veg = "Plot") {
  veg <- standardise_vegetation_sites(vegetation_df)

  common <- intersect(pcoa_result$site_scores$site, veg$site)
  if (length(common) < 4) {
    warning("Fewer than 4 common sites.")
    return(NULL)
  }

  sc <- pcoa_result$site_scores %>% filter(site %in% common) %>% arrange(site)
  vm <- veg %>% filter(site %in% common) %>% arrange(site)

  ord_mat <- as.matrix(sc[, c("Axis1", "Axis2")])
  rownames(ord_mat) <- sc$site

  env_cols <- setdiff(colnames(vm), "site")
  env_cols <- env_cols[sapply(vm[, env_cols, drop = FALSE], is.numeric)]
  env_mat <- vm[, env_cols, drop = FALSE]
  env_mat <- env_mat[, apply(env_mat, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(env_mat) == 0) return(NULL)

  ef <- envfit(ord_mat, env_mat, permutations = n_perm, na.rm = TRUE)

  vectors <- as.data.frame(scores(ef, display = "vectors"))
  vectors$variable <- rownames(vectors)
  vectors$r2 <- ef$vectors$r
  vectors$p_value <- ef$vectors$pvals

  list(envfit = ef, vectors = vectors,
       significant = vectors %>% filter(p_value < 0.05))
}


#' Plot PCoA with optional envfit arrows
#'
#' @param pcoa_result List from run_pcoa_jaccard.
#' @param envfit_result List from run_envfit or NULL.
#' @param richness_df Data frame with site and richness or NULL.
#' @param taxon_label Character.
#' @return ggplot object.
#' @export
plot_pcoa <- function(pcoa_result, envfit_result = NULL,
                      richness_df = NULL, taxon_label = "Katydid") {
  df <- pcoa_result$site_scores
  eig <- pcoa_result$eigenvalues
  ep <- eig[eig > 0]
  vpct <- ep / sum(ep) * 100

  if (!is.null(richness_df))
    df <- left_join(df, richness_df, by = "site")

  p <- ggplot(df, aes(x = Axis1, y = Axis2))
  if ("richness" %in% colnames(df)) {
    p <- p + geom_point(aes(size = richness), alpha = 0.7) +
      scale_size_continuous(name = "Richness", range = c(2, 6))
  } else {
    p <- p + geom_point(size = 3, alpha = 0.7)
  }
  p <- p +
    geom_text(aes(label = site), vjust = -0.8, size = 2.5, alpha = 0.7) +
    labs(x = sprintf("PCoA Axis 1 (%.1f%%)", vpct[1]),
         y = sprintf("PCoA Axis 2 (%.1f%%)", vpct[2]),
         title = paste(taxon_label, "community ordination (PCoA, Jaccard)")) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold"))

  if (!is.null(envfit_result) && nrow(envfit_result$significant) > 0) {
    sc <- max(abs(c(df$Axis1, df$Axis2))) * 0.8
    arr <- envfit_result$significant %>%
      mutate(x_end = Axis1 * sc, y_end = Axis2 * sc)
    p <- p +
      geom_segment(data = arr, aes(x = 0, y = 0, xend = x_end, yend = y_end),
                   arrow = arrow(length = unit(0.2, "cm")),
                   color = "firebrick", linewidth = 0.7, inherit.aes = FALSE) +
      geom_text(data = arr,
                aes(x = x_end * 1.1, y = y_end * 1.1, label = variable),
                color = "firebrick", size = 3, inherit.aes = FALSE)
  }
  p
}


#' Test katydid-bird richness relationship with GLM Poisson
#'
#' @param katydid_richness Data frame with site and richness.
#' @param bird_richness Data frame with site and richness.
#' @return List with data, correlation, glm_model, glm_p, plot, or NULL.
#' @export
test_inter_taxa_richness <- function(katydid_richness, bird_richness) {
  merged <- inner_join(
    katydid_richness %>% rename(katydid_richness = richness),
    bird_richness %>% rename(bird_richness = richness),
    by = "site")
  if (nrow(merged) < 4) return(NULL)

  ct <- cor.test(merged$katydid_richness, merged$bird_richness,
                 method = "spearman", exact = FALSE)
  model <- glm(katydid_richness ~ bird_richness, data = merged, family = poisson)
  ms <- summary(model)
  glm_p <- coef(ms)[2, 4]

  p <- ggplot(merged, aes(x = bird_richness, y = katydid_richness)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "glm", method.args = list(family = "poisson"),
                se = TRUE, color = "steelblue", linewidth = 0.8) +
    labs(x = "Bird species richness", y = "Katydid species richness",
         title = "Katydid vs bird richness (GLM Poisson)",
         subtitle = sprintf("Spearman rho=%.3f (p=%.4f); GLM p=%.4f",
                            ct$estimate, ct$p.value, glm_p)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold"))

  list(data = merged, correlation = ct, glm_model = model,
       glm_p = glm_p, plot = p)
}


# -- load data --
cat("Loading data...\n")
katydid_matrix <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                           check.names = FALSE)
bird_matrix <- read.csv(file.path(INPUT_DIR, "bird_presence_matrix.csv"),
                        check.names = FALSE)
vegetation <- read.csv(file.path(INPUT_DIR, "vegetation_summary.csv"),
                       stringsAsFactors = FALSE)

katydid_richness <- richness_from_matrix(katydid_matrix)
bird_richness <- richness_from_matrix(bird_matrix)

# -- PCoA --
cat("Running PCoA (Jaccard)...\n")
pcoa_result <- run_pcoa_jaccard(katydid_matrix)
if (!is.null(pcoa_result)) {
  write.csv(pcoa_result$site_scores,
            file.path(OUTPUT_DIR, "pcoa_site_scores.csv"), row.names = FALSE)
  eig_df <- data.frame(axis = seq_along(pcoa_result$eigenvalues),
                       eigenvalue = pcoa_result$eigenvalues)
  write.csv(eig_df, file.path(OUTPUT_DIR, "pcoa_eigenvalues.csv"),
            row.names = FALSE)
}

# -- envfit --
cat("Running envfit (%d permutations)...\n", N_PERMUTATIONS)
envfit_result <- NULL
if (!is.null(pcoa_result)) {
  envfit_result <- run_envfit(pcoa_result, vegetation, n_perm = N_PERMUTATIONS)
  if (!is.null(envfit_result))
    write.csv(envfit_result$vectors,
              file.path(OUTPUT_DIR, "envfit_results.csv"), row.names = FALSE)
}

# -- PCoA figure --
if (!is.null(pcoa_result)) {
  p <- plot_pcoa(pcoa_result, envfit_result, katydid_richness, "Katydid")
  ggsave(file.path(FIGURE_DIR, "katydid_pcoa_jaccard.jpg"), p,
         width = 9, height = 7, dpi = 300)
}

# -- katydid-bird relationship --
cat("Testing katydid-bird richness...\n")
taxa_res <- test_inter_taxa_richness(katydid_richness, bird_richness)
if (!is.null(taxa_res)) {
  write.csv(taxa_res$data,
            file.path(OUTPUT_DIR, "katydid_bird_richness.csv"),
            row.names = FALSE)
  ggsave(file.path(FIGURE_DIR, "katydid_bird_richness.jpg"),
         taxa_res$plot, width = 7, height = 6, dpi = 300)
}

# -- summary --
cat("\n-- Community structure summary --\n")
if (!is.null(pcoa_result)) {
  eig <- pcoa_result$eigenvalues
  ep <- eig[eig > 0]
  vpct <- ep / sum(ep) * 100
  cat(sprintf("  PCoA Axis 1: %.1f%%, Axis 2: %.1f%%\n", vpct[1], vpct[2]))
}
if (!is.null(envfit_result) && nrow(envfit_result$significant) > 0) {
  cat("  Significant envfit variables:\n")
  for (i in seq_len(nrow(envfit_result$significant))) {
    r <- envfit_result$significant[i, ]
    cat(sprintf("    %s: r2=%.3f, p=%.4f\n", r$variable, r$r2, r$p_value))
  }
}
if (!is.null(taxa_res)) {
  cat(sprintf("  Katydid-bird: rho=%.3f (p=%.4f), GLM p=%.4f\n",
              taxa_res$correlation$estimate,
              taxa_res$correlation$p.value, taxa_res$glm_p))
}
cat("\nCommunity structure analysis complete.\n")
