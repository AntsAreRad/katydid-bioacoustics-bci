# baseline_species_distribution.R
# For each bioacoustically detected species:
#   - spatial distribution (which sites, how many sites)
#   - temporal distribution (which months, first/last detection)
#   - average number of detected species per day (across sites)
#
# Input: integrated_results/ CSV files
# Output: results/species_distribution/ CSV and figures

library(tidyverse)
library(lubridate)
library(ggplot2)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/species_distribution"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

UTC_OFFSET_HOURS <- -5


#' Build species-level spatial and temporal summary
#'
#' For each species detected by bioacoustics, compute:
#'   - number of sites where detected
#'   - list of sites
#'   - number of distinct detection days
#'   - first and last detection date
#'   - months with detections
#'
#' @param detections Data frame with file_path, site, common_name.
#' @param species_col Character. Default "common_name".
#' @param utc_offset Numeric. UTC offset.
#' @return Data frame with one row per species.
#' @export
build_species_summary <- function(detections, species_col = "common_name",
                                  utc_offset = -5) {
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))

  det %>%
    group_by(species = !!sym(species_col)) %>%
    summarise(
      n_sites = n_distinct(site),
      sites = paste(sort(unique(site)), collapse = "; "),
      n_detection_days = n_distinct(local_date),
      first_detection = min(local_date),
      last_detection = max(local_date),
      months_detected = paste(sort(unique(format(local_date, "%Y-%m"))),
                              collapse = "; "),
      total_detections = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(n_sites), desc(n_detection_days))
}


#' Calculate average number of species detected per day
#'
#' For each site-day combination, count distinct species. Then summarise
#' across all site-days.
#'
#' @param detections Data frame with file_path, site, common_name.
#' @param species_col Character.
#' @param utc_offset Numeric.
#' @return List with daily_richness (per site-day), site_means (per site),
#'   and global summary.
#' @export
calculate_daily_species_count <- function(detections,
                                          species_col = "common_name",
                                          utc_offset = -5) {
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))

  # species count per site per day
  daily <- det %>%
    group_by(site, local_date) %>%
    summarise(n_species = n_distinct(!!sym(species_col)), .groups = "drop")

  # mean per site
  site_means <- daily %>%
    group_by(site) %>%
    summarise(
      mean_species_per_day = mean(n_species),
      sd_species_per_day = sd(n_species),
      n_days = n(),
      .groups = "drop"
    )

  # global
  global <- daily %>%
    summarise(
      mean_species_per_day = mean(n_species),
      sd = sd(n_species),
      median = median(n_species),
      min = min(n_species),
      max = max(n_species),
      n_site_days = n()
    )

  list(daily_richness = daily, site_means = site_means, global = global)
}


#' Plot average species per day by site
#'
#' @param site_means Data frame from calculate_daily_species_count.
#' @param taxon_label Character.
#' @param fill_color Character.
#' @return ggplot object.
#' @export
plot_daily_species_by_site <- function(site_means, taxon_label = "Katydid",
                                       fill_color = "#009E73") {
  ggplot(site_means, aes(x = reorder(site, mean_species_per_day),
                         y = mean_species_per_day)) +
    geom_col(fill = fill_color, width = 0.7) +
    geom_errorbar(aes(ymin = pmax(mean_species_per_day - sd_species_per_day, 0),
                      ymax = mean_species_per_day + sd_species_per_day),
                  width = 0.3, linewidth = 0.4) +
    coord_flip() +
    labs(x = "Site", y = "Mean species detected per day",
         title = paste(taxon_label,
                       "average number of species detected per day")) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank(),
          plot.title = element_text(size = 13, face = "bold"))
}


#' Plot monthly boxplots of daily species counts
#'
#' Replaces the histogram. Shows a boxplot per month (x = month, y = daily
#' species count), revealing temporal and seasonal variation.
#'
#' @param daily_richness Data frame with local_date and n_species columns.
#' @param taxon_label Character.
#' @param fill_color Character.
#' @return ggplot object.
#' @export
plot_daily_species_monthly_boxplot <- function(daily_richness,
                                               taxon_label = "Katydid",
                                               fill_color = "#009E73") {
  df <- daily_richness %>%
    mutate(year_month = format(local_date, "%Y-%m")) %>%
    arrange(year_month) %>%
    mutate(year_month = factor(year_month, levels = unique(year_month)))

  # Annotate dry/wet season
  dry_months <- c("01", "02", "03", "04", "12")
  df <- df %>%
    mutate(
      month_num = format(local_date, "%m"),
      season = ifelse(month_num %in% dry_months, "Dry", "Wet")
    )

  ggplot(df, aes(x = year_month, y = n_species, fill = season)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
    scale_fill_manual(values = c("Dry" = "#E69F00", "Wet" = fill_color),
                      name = "Season") +
    labs(x = "Month", y = "Number of species detected per day",
         title = paste(taxon_label,
                       "daily species richness by month"),
         subtitle = "BCI dry season: mid-Dec to mid-Apr") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "bottom")
}


# -- load data --
cat("Loading detection data...\n")
katydid_det <- read.csv(file.path(INPUT_DIR, "katydid_detections.csv"),
                        stringsAsFactors = FALSE)
bird_det <- read.csv(file.path(INPUT_DIR, "bird_detections.csv"),
                     stringsAsFactors = FALSE)

# -- species-level summary --
cat("Building species-level summaries...\n")
katydid_species <- build_species_summary(katydid_det, utc_offset = UTC_OFFSET_HOURS)
bird_species <- build_species_summary(bird_det, utc_offset = UTC_OFFSET_HOURS)

write.csv(katydid_species,
          file.path(OUTPUT_DIR, "katydid_species_summary.csv"),
          row.names = FALSE)
write.csv(bird_species,
          file.path(OUTPUT_DIR, "bird_species_summary.csv"),
          row.names = FALSE)
cat("  Saved species summaries\n")

# -- average species per day --
cat("Calculating daily species counts...\n")
katydid_daily <- calculate_daily_species_count(katydid_det,
                                               utc_offset = UTC_OFFSET_HOURS)
bird_daily <- calculate_daily_species_count(bird_det,
                                            utc_offset = UTC_OFFSET_HOURS)

write.csv(katydid_daily$daily_richness,
          file.path(OUTPUT_DIR, "katydid_daily_richness.csv"),
          row.names = FALSE)
write.csv(katydid_daily$site_means,
          file.path(OUTPUT_DIR, "katydid_daily_means_by_site.csv"),
          row.names = FALSE)
write.csv(bird_daily$daily_richness,
          file.path(OUTPUT_DIR, "bird_daily_richness.csv"),
          row.names = FALSE)
write.csv(bird_daily$site_means,
          file.path(OUTPUT_DIR, "bird_daily_means_by_site.csv"),
          row.names = FALSE)

# -- figures --
cat("Creating figures...\n")
# Barplots of mean species per day per site (kept)
p1 <- plot_daily_species_by_site(katydid_daily$site_means, "Katydid", "#009E73")
ggsave(file.path(FIGURE_DIR, "katydid_mean_species_per_day.jpg"), p1,
       width = 8, height = 6, dpi = 300)
p2 <- plot_daily_species_by_site(bird_daily$site_means, "Bird", "#0072B2")
ggsave(file.path(FIGURE_DIR, "bird_mean_species_per_day.jpg"), p2,
       width = 8, height = 6, dpi = 300)

# Monthly boxplots (replace histograms)
p3 <- plot_daily_species_monthly_boxplot(katydid_daily$daily_richness,
                                          "Katydid", "#009E73")
ggsave(file.path(FIGURE_DIR, "katydid_daily_richness_by_month.jpg"), p3,
       width = 9, height = 6, dpi = 300)
p4 <- plot_daily_species_monthly_boxplot(bird_daily$daily_richness,
                                          "Bird", "#0072B2")
ggsave(file.path(FIGURE_DIR, "bird_daily_richness_by_month.jpg"), p4,
       width = 9, height = 6, dpi = 300)

# -- summary --
cat("\n-- Species distribution summary --\n")
cat(sprintf("  Katydid species: %d total\n", nrow(katydid_species)))
cat(sprintf("  Katydid mean species per day: %.1f (sd=%.1f)\n",
            katydid_daily$global$mean_species_per_day,
            katydid_daily$global$sd))
cat(sprintf("  Bird species: %d total\n", nrow(bird_species)))
cat(sprintf("  Bird mean species per day: %.1f (sd=%.1f)\n",
            bird_daily$global$mean_species_per_day,
            bird_daily$global$sd))
cat("\nSpecies distribution analysis complete.\n")
