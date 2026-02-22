# baseline_site_composition.R
# Site-level species composition similarity analysis.
#
# The mean number of katydid species per day is similar across sites.
# This script answers: "Is that mostly the same species at all sites,
# or are there different species but about the same number?"
#
# Analyses:
#   a) Pairwise Jaccard similarity matrix between all 25 sites
#   b) Beta diversity decomposition (turnover vs nestedness) via betapart
#   c) Species-by-site occurrence heatmap
#   d) Summary statistics for the paper
#
# Input: integrated_results/ CSV files (presence matrices)
# Output: results/site_composition/ CSV and figures

library(tidyverse)
library(ggplot2)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR  <- "integrated_results"
OUTPUT_DIR <- "results/site_composition"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


# --- functions ---------------------------------------------------------------

#' Calculate pairwise Jaccard similarity between sites
#'
#' @param presence_matrix Data frame with site as first column, species as
#'   remaining columns (0/1 values).
#' @return List with similarity_matrix (sites x sites), mean_similarity,
#'   sd_similarity, and long-format data frame for plotting.
#' @export
calculate_jaccard_similarity <- function(presence_matrix) {
  sites <- presence_matrix$site
  sp_cols <- setdiff(colnames(presence_matrix), "site")
  mat <- as.matrix(presence_matrix[, sp_cols])
  rownames(mat) <- sites

  # Remove empty sites/species
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  n <- nrow(mat)

  # Pairwise Jaccard similarity (1 - distance)
  sim_mat <- matrix(NA, n, n, dimnames = list(rownames(mat), rownames(mat)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      shared <- sum(mat[i, ] == 1 & mat[j, ] == 1)
      total  <- sum(mat[i, ] == 1 | mat[j, ] == 1)
      sim_mat[i, j] <- if (total == 0) 0 else shared / total
    }
  }

  # Long format (upper triangle only for stats)
  pairs <- expand.grid(site1 = rownames(mat), site2 = rownames(mat),
                        stringsAsFactors = FALSE) %>%
    mutate(jaccard = as.vector(sim_mat)) %>%
    filter(site1 < site2)

  list(
    similarity_matrix = sim_mat,
    pairs = pairs,
    mean_similarity = mean(pairs$jaccard),
    sd_similarity = sd(pairs$jaccard),
    median_similarity = median(pairs$jaccard),
    n_sites = n,
    n_species = ncol(mat)
  )
}


#' Plot Jaccard similarity heatmap
#'
#' @param sim_result List from calculate_jaccard_similarity.
#' @param taxon_label Character.
#' @return ggplot object.
#' @export
plot_jaccard_heatmap <- function(sim_result, taxon_label = "Katydid") {
  mat <- sim_result$similarity_matrix
  # Order sites by similarity (hierarchical clustering)
  hc <- hclust(as.dist(1 - mat), method = "average")
  site_order <- rownames(mat)[hc$order]

  df <- expand.grid(site1 = rownames(mat), site2 = colnames(mat),
                    stringsAsFactors = FALSE) %>%
    mutate(jaccard = as.vector(mat),
           site1 = factor(site1, levels = site_order),
           site2 = factor(site2, levels = site_order))

  ggplot(df, aes(x = site1, y = site2, fill = jaccard)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.0f", jaccard * 100)),
              size = 2.2, color = ifelse(df$jaccard > 0.6, "white", "black")) +
    scale_fill_gradient2(low = "#f7fbff", mid = "#6baed6", high = "#08306b",
                         midpoint = 0.5, limits = c(0, 1),
                         name = "Jaccard\nsimilarity") +
    labs(x = NULL, y = NULL,
         title = paste(taxon_label, "pairwise site similarity (Jaccard)"),
         subtitle = sprintf("Mean = %.1f%% (SD = %.1f%%), n = %d sites",
                            sim_result$mean_similarity * 100,
                            sim_result$sd_similarity * 100,
                            sim_result$n_sites)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "right") +
    coord_fixed()
}


#' Decompose beta diversity into turnover and nestedness
#'
#' Uses the betapart package. If unavailable, falls back to manual calculation
#' of the Sorensen-based decomposition (Baselga 2010).
#'
#' @param presence_matrix Data frame with site as first column.
#' @return List with total_beta (Sorensen), turnover, nestedness, fractions,
#'   pairwise data frames, and whether betapart was used.
#' @export
decompose_beta_diversity <- function(presence_matrix) {
  sites <- presence_matrix$site
  sp_cols <- setdiff(colnames(presence_matrix), "site")
  mat <- as.matrix(presence_matrix[, sp_cols])
  rownames(mat) <- sites
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  use_betapart <- requireNamespace("betapart", quietly = TRUE)

  if (use_betapart) {
    cat("[>>] Using betapart package for beta diversity decomposition\n")
    bp <- betapart::beta.pair(mat, index.family = "sorensen")
    total   <- as.matrix(bp$beta.sor)
    turnov  <- as.matrix(bp$beta.sim)
    nested  <- as.matrix(bp$beta.sne)

    # Multi-site metrics
    bp_multi <- betapart::beta.multi(mat, index.family = "sorensen")
    multi_total    <- bp_multi$beta.SOR
    multi_turnover <- bp_multi$beta.SIM
    multi_nested   <- bp_multi$beta.SNE
  } else {
    cat("[--] betapart not installed, using manual decomposition\n")
    n <- nrow(mat)
    total <- turnov <- nested <- matrix(0, n, n,
                                         dimnames = list(rownames(mat), rownames(mat)))
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) next
        a <- sum(mat[i, ] == 1 & mat[j, ] == 1)  # shared
        b <- sum(mat[i, ] == 1 & mat[j, ] == 0)  # only in i
        c_val <- sum(mat[i, ] == 0 & mat[j, ] == 1)  # only in j
        # Sorensen dissimilarity
        total[i, j] <- (b + c_val) / (2 * a + b + c_val)
        # Turnover (Simpson)
        turnov[i, j] <- min(b, c_val) / (a + min(b, c_val))
        # Nestedness = total - turnover
        nested[i, j] <- total[i, j] - turnov[i, j]
      }
    }
    # Multi-site (mean of pairwise as approximation)
    idx <- upper.tri(total)
    multi_total    <- mean(total[idx])
    multi_turnover <- mean(turnov[idx])
    multi_nested   <- mean(nested[idx])
  }

  # Pairwise means (upper triangle)
  idx <- upper.tri(total)
  list(
    total_beta = multi_total,
    turnover   = multi_turnover,
    nestedness = multi_nested,
    turnover_fraction  = multi_turnover / multi_total * 100,
    nestedness_fraction = multi_nested / multi_total * 100,
    pairwise_total    = total,
    pairwise_turnover = turnov,
    pairwise_nested   = nested,
    mean_pairwise_total    = mean(total[idx]),
    mean_pairwise_turnover = mean(turnov[idx]),
    mean_pairwise_nested   = mean(nested[idx]),
    n_sites = nrow(mat),
    n_species = ncol(mat),
    used_betapart = use_betapart
  )
}


#' Plot species-by-site occurrence heatmap
#'
#' Shows which species occur at which sites, ordered by prevalence.
#'
#' @param presence_matrix Data frame with site as first column.
#' @param taxon_label Character.
#' @param fill_color Character.
#' @return ggplot object.
#' @export
plot_species_site_heatmap <- function(presence_matrix, taxon_label = "Katydid",
                                      fill_color = "#009E73") {
  sites <- presence_matrix$site
  sp_cols <- setdiff(colnames(presence_matrix), "site")
  mat <- as.matrix(presence_matrix[, sp_cols])
  rownames(mat) <- sites

  # Remove empty
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat) > 0, drop = FALSE]

  # Order species by prevalence (most widespread first)
  sp_order <- names(sort(colSums(mat), decreasing = TRUE))
  # Order sites by richness
  site_order <- names(sort(rowSums(mat), decreasing = TRUE))

  df <- expand.grid(site = rownames(mat), species = colnames(mat),
                    stringsAsFactors = FALSE) %>%
    mutate(present = as.vector(mat),
           site = factor(site, levels = rev(site_order)),
           species = factor(species, levels = sp_order))

  # Compute species prevalence for axis labels
  prevalence <- colSums(mat)
  sp_labels <- paste0(sp_order, " (", prevalence[sp_order], ")")

  ggplot(df, aes(x = species, y = site, fill = factor(present))) +
    geom_tile(color = "grey90", linewidth = 0.2) +
    scale_fill_manual(values = c("0" = "white", "1" = fill_color),
                      labels = c("Absent", "Present"),
                      name = NULL) +
    scale_x_discrete(labels = sp_labels) +
    labs(x = "Species (# sites detected)",
         y = "Site",
         title = paste(taxon_label, "species occurrence across sites"),
         subtitle = sprintf("%d species at %d sites", ncol(mat), nrow(mat))) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "bottom")
}


#' Plot beta diversity decomposition (Sorensen partitioning)
#'
#' @param beta_katydid List from decompose_beta_diversity.
#' @param beta_bird List from decompose_beta_diversity (or NULL).
#' @return ggplot object.
#' @export
plot_beta_decomposition <- function(beta_katydid, beta_bird = NULL) {
  df <- data.frame(
    taxon = "Katydid",
    component = c("Replacement (beta.SIM)", "Richness difference (beta.SNE)"),
    value = c(beta_katydid$turnover, beta_katydid$nestedness),
    fraction = c(beta_katydid$turnover_fraction, beta_katydid$nestedness_fraction)
  )
  if (!is.null(beta_bird)) {
    df <- bind_rows(df, data.frame(
      taxon = "Bird",
      component = c("Replacement (beta.SIM)", "Richness difference (beta.SNE)"),
      value = c(beta_bird$turnover, beta_bird$nestedness),
      fraction = c(beta_bird$turnover_fraction, beta_bird$nestedness_fraction)
    ))
  }
  df$component <- factor(df$component,
                         levels = c("Replacement (beta.SIM)",
                                    "Richness difference (beta.SNE)"))

  ggplot(df, aes(x = taxon, y = value, fill = component)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1f%%", fraction)),
              position = position_stack(vjust = 0.5), size = 4, fontface = "bold") +
    scale_fill_manual(values = c("Replacement (beta.SIM)" = "#D55E00",
                                  "Richness difference (beta.SNE)" = "#56B4E9"),
                      name = "Component") +
    labs(x = NULL, y = "Beta diversity (Sorensen)",
         title = "Sorensen dissimilarity partitioning (Baselga 2010)",
         subtitle = "beta.SOR = beta.SIM (species replacement) + beta.SNE (richness difference)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          legend.position = "bottom")
}


# --- execute -----------------------------------------------------------------

cat("Loading presence matrices...\n")
katydid_pm <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                       stringsAsFactors = FALSE)
bird_pm    <- read.csv(file.path(INPUT_DIR, "bird_presence_matrix.csv"),
                       stringsAsFactors = FALSE)

# -- Jaccard similarity --
cat("[>>] Calculating pairwise Jaccard similarity...\n")
katydid_sim <- calculate_jaccard_similarity(katydid_pm)
bird_sim    <- calculate_jaccard_similarity(bird_pm)
cat(sprintf("[OK] Katydid: mean Jaccard = %.1f%% (SD %.1f%%), %d sites, %d species\n",
            katydid_sim$mean_similarity * 100, katydid_sim$sd_similarity * 100,
            katydid_sim$n_sites, katydid_sim$n_species))
cat(sprintf("[OK] Bird: mean Jaccard = %.1f%% (SD %.1f%%), %d sites, %d species\n",
            bird_sim$mean_similarity * 100, bird_sim$sd_similarity * 100,
            bird_sim$n_sites, bird_sim$n_species))

# Save similarity matrices
write.csv(as.data.frame(katydid_sim$similarity_matrix),
          file.path(OUTPUT_DIR, "katydid_jaccard_similarity_matrix.csv"))
write.csv(as.data.frame(bird_sim$similarity_matrix),
          file.path(OUTPUT_DIR, "bird_jaccard_similarity_matrix.csv"))
write.csv(katydid_sim$pairs,
          file.path(OUTPUT_DIR, "katydid_jaccard_pairwise.csv"),
          row.names = FALSE)
write.csv(bird_sim$pairs,
          file.path(OUTPUT_DIR, "bird_jaccard_pairwise.csv"),
          row.names = FALSE)

# Heatmaps
p1 <- plot_jaccard_heatmap(katydid_sim, "Katydid")
ggsave(file.path(FIGURE_DIR, "katydid_jaccard_heatmap.jpg"), p1,
       width = 10, height = 9, dpi = 300)
p2 <- plot_jaccard_heatmap(bird_sim, "Bird")
ggsave(file.path(FIGURE_DIR, "bird_jaccard_heatmap.jpg"), p2,
       width = 10, height = 9, dpi = 300)

# -- Beta diversity decomposition --
cat("[>>] Decomposing beta diversity (turnover vs nestedness)...\n")
katydid_beta <- decompose_beta_diversity(katydid_pm)
bird_beta    <- decompose_beta_diversity(bird_pm)

cat(sprintf("[OK] Katydid beta: total=%.3f, turnover=%.3f (%.1f%%), nestedness=%.3f (%.1f%%)\n",
            katydid_beta$total_beta, katydid_beta$turnover,
            katydid_beta$turnover_fraction,
            katydid_beta$nestedness, katydid_beta$nestedness_fraction))
cat(sprintf("[OK] Bird beta: total=%.3f, turnover=%.3f (%.1f%%), nestedness=%.3f (%.1f%%)\n",
            bird_beta$total_beta, bird_beta$turnover,
            bird_beta$turnover_fraction,
            bird_beta$nestedness, bird_beta$nestedness_fraction))

# Save beta diversity summary
beta_summary <- data.frame(
  taxon = c("Katydid", "Bird"),
  beta_SOR = c(katydid_beta$total_beta, bird_beta$total_beta),
  beta_SIM = c(katydid_beta$turnover, bird_beta$turnover),
  beta_SNE = c(katydid_beta$nestedness, bird_beta$nestedness),
  beta_SIM_pct = c(katydid_beta$turnover_fraction, bird_beta$turnover_fraction),
  beta_SNE_pct = c(katydid_beta$nestedness_fraction, bird_beta$nestedness_fraction),
  n_sites = c(katydid_beta$n_sites, bird_beta$n_sites),
  n_species = c(katydid_beta$n_species, bird_beta$n_species)
)
write.csv(beta_summary, file.path(OUTPUT_DIR, "beta_diversity_summary.csv"),
          row.names = FALSE)

# Beta decomposition barplot
p3 <- plot_beta_decomposition(katydid_beta, bird_beta)
ggsave(file.path(FIGURE_DIR, "beta_decomposition.jpg"), p3,
       width = 7, height = 6, dpi = 300)

# -- Species-by-site heatmaps --
cat("[>>] Creating species-by-site occurrence heatmaps...\n")
p4 <- plot_species_site_heatmap(katydid_pm, "Katydid", "#009E73")
ggsave(file.path(FIGURE_DIR, "katydid_species_site_heatmap.jpg"), p4,
       width = 14, height = 8, dpi = 300)
p5 <- plot_species_site_heatmap(bird_pm, "Bird", "#0072B2")
ggsave(file.path(FIGURE_DIR, "bird_species_site_heatmap.jpg"), p5,
       width = 18, height = 8, dpi = 300)

# -- Summary for Laurel's question --
cat("\n-- Site composition summary --\n")
cat(sprintf("  Katydid mean Jaccard similarity: %.1f%% (SD %.1f%%)\n",
            katydid_sim$mean_similarity * 100, katydid_sim$sd_similarity * 100))
cat(sprintf("  Sorensen partitioning: beta.SIM = %.1f%%, beta.SNE = %.1f%%\n",
            katydid_beta$turnover_fraction,
            katydid_beta$nestedness_fraction))
cat(sprintf("  Bird mean Jaccard similarity: %.1f%% (SD %.1f%%)\n",
            bird_sim$mean_similarity * 100, bird_sim$sd_similarity * 100))
cat(sprintf("  Sorensen partitioning: beta.SIM = %.1f%%, beta.SNE = %.1f%%\n",
            bird_beta$turnover_fraction,
            bird_beta$nestedness_fraction))

cat("\nSite composition analysis complete.\n")
