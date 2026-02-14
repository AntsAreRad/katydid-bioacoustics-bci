# baseline_detection_richness.R
# Is the number of detections positively or negatively correlated
# with the number of species? Test per site using Spearman.
#
# Also identifies potential outliers (sites far from the trend) and
# reports Spearman results with and without the outlier.
#
# Input: integrated_results/ CSV files
# Output: results/detection_richness/ CSV and figures

library(tidyverse)
library(ggplot2)
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
if (has_ggrepel) library(ggrepel)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/detection_richness"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


#' Calculate total detection count per site
#'
#' @param detections Data frame with at least site and common_name columns.
#' @return Data frame with site and total_detections.
#' @export
calculate_site_detections <- function(detections) {
  detections %>%
    group_by(site) %>%
    summarise(total_detections = n(), .groups = "drop")
}


#' Test correlation between detection count and richness
#'
#' Performs Spearman correlation on all sites. If an outlier is detected
#' (using studentized residuals from a linear model, |resid| > 2.5),
#' also reports the test excluding the outlier.
#'
#' @param detections Data frame of raw detections.
#' @param presence_matrix Data frame. Presence/absence matrix.
#' @param taxon Character. Label.
#' @return List with correlation test, merged data, ggplot, and
#'   optionally outlier info and correlation without outlier.
#' @export
correlate_detection_richness <- function(detections, presence_matrix,
                                         taxon = "unknown") {
  det_counts <- calculate_site_detections(detections)
  richness <- richness_from_matrix(presence_matrix)
  merged <- inner_join(det_counts, richness, by = "site")

  ct <- cor.test(merged$total_detections, merged$richness,
                 method = "spearman", exact = FALSE)

  direction <- if (ct$estimate > 0) "positive" else "negative"
  subtitle <- sprintf("Spearman rho = %.3f (%s), p = %.4f, n = %d sites",
                      ct$estimate, direction, ct$p.value, nrow(merged))

  # Identify outliers via studentized residuals
  lm_fit <- lm(richness ~ total_detections, data = merged)
  merged$resid <- rstudent(lm_fit)
  outlier_sites <- merged %>% filter(abs(resid) > 2.5) %>% pull(site)

  # If outlier found, also run Spearman without it
  ct_no_outlier <- NULL
  subtitle_extra <- ""
  if (length(outlier_sites) > 0) {
    merged_clean <- merged %>% filter(!site %in% outlier_sites)
    ct_no_outlier <- cor.test(merged_clean$total_detections,
                              merged_clean$richness,
                              method = "spearman", exact = FALSE)
    dir_no <- if (ct_no_outlier$estimate > 0) "positive" else "negative"
    subtitle_extra <- sprintf(
      "\nWithout %s: rho = %.3f (%s), p = %.4f, n = %d",
      paste(outlier_sites, collapse = ", "),
      ct_no_outlier$estimate, dir_no, ct_no_outlier$p.value,
      nrow(merged_clean))
  }

  # Site labels
  label_layer <- if (has_ggrepel) {
    ggrepel::geom_text_repel(aes(label = site), size = 3, max.overlaps = 30,
                              segment.color = "grey60", segment.size = 0.3)
  } else {
    geom_text(aes(label = site), size = 3, hjust = -0.2, vjust = -0.5)
  }

  p <- ggplot(merged, aes(x = total_detections, y = richness)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40", linewidth = 0.8) +
    label_layer +
    labs(x = "Total detection count per site",
         y = "Species richness per site",
         title = paste(tools::toTitleCase(taxon),
                       "detection count vs species richness"),
         subtitle = paste0(subtitle, subtitle_extra)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 13, face = "bold"),
          plot.subtitle = element_text(size = 9))

  result <- list(correlation = ct, data = merged, plot = p, taxon = taxon,
                 direction = direction, outlier_sites = outlier_sites)
  if (!is.null(ct_no_outlier)) {
    result$correlation_no_outlier <- ct_no_outlier
  }
  result
}


# -- load data --
cat("Loading data...\n")
katydid_det <- read.csv(file.path(INPUT_DIR, "katydid_detections.csv"),
                        stringsAsFactors = FALSE)
bird_det <- read.csv(file.path(INPUT_DIR, "bird_detections.csv"),
                     stringsAsFactors = FALSE)
katydid_matrix <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                           check.names = FALSE)
bird_matrix <- read.csv(file.path(INPUT_DIR, "bird_presence_matrix.csv"),
                        check.names = FALSE)

# -- analysis --
cat("Testing detection count vs richness...\n")
katydid_res <- correlate_detection_richness(katydid_det, katydid_matrix, "katydid")
bird_res <- correlate_detection_richness(bird_det, bird_matrix, "bird")

# -- save --
write.csv(katydid_res$data %>% select(-resid),
          file.path(OUTPUT_DIR, "katydid_detection_vs_richness.csv"),
          row.names = FALSE)
write.csv(bird_res$data %>% select(-resid),
          file.path(OUTPUT_DIR, "bird_detection_vs_richness.csv"),
          row.names = FALSE)

# Build summary table including outlier info
summary_rows <- list()
for (res in list(katydid_res, bird_res)) {
  row <- data.frame(
    taxon = res$taxon,
    spearman_rho = res$correlation$estimate,
    direction = res$direction,
    p_value = res$correlation$p.value,
    n_sites = nrow(res$data),
    outlier_sites = paste(res$outlier_sites, collapse = "; "),
    stringsAsFactors = FALSE
  )
  if (!is.null(res$correlation_no_outlier)) {
    row$rho_no_outlier <- res$correlation_no_outlier$estimate
    row$p_no_outlier <- res$correlation_no_outlier$p.value
    row$n_no_outlier <- nrow(res$data) - length(res$outlier_sites)
  } else {
    row$rho_no_outlier <- NA
    row$p_no_outlier <- NA
    row$n_no_outlier <- NA
  }
  summary_rows[[length(summary_rows) + 1]] <- row
}
summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(OUTPUT_DIR, "detection_richness_summary.csv"),
          row.names = FALSE)

ggsave(file.path(FIGURE_DIR, "katydid_detection_vs_richness.jpg"),
       katydid_res$plot, width = 8, height = 7, dpi = 300)
ggsave(file.path(FIGURE_DIR, "bird_detection_vs_richness.jpg"),
       bird_res$plot, width = 8, height = 7, dpi = 300)

# -- summary --
cat("\n-- Detection count vs richness --\n")
for (res in list(katydid_res, bird_res)) {
  cat(sprintf("  %s: rho=%.3f (%s), p=%.4f\n",
              tools::toTitleCase(res$taxon),
              res$correlation$estimate, res$direction,
              res$correlation$p.value))
  if (length(res$outlier_sites) > 0) {
    cat(sprintf("    Outlier(s): %s\n",
                paste(res$outlier_sites, collapse = ", ")))
    if (!is.null(res$correlation_no_outlier)) {
      dir_no <- if (res$correlation_no_outlier$estimate > 0) "positive" else "negative"
      cat(sprintf("    Without outlier: rho=%.3f (%s), p=%.4f\n",
                  res$correlation_no_outlier$estimate, dir_no,
                  res$correlation_no_outlier$p.value))
    }
  }
}
cat("\nDetection-richness analysis complete.\n")
