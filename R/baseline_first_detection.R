# baseline_first_detection.R
# Mean relative first detection day per species (mobility indicator).
#
# For each species, at each site where it is detected within a given month,
# compute the "relative day" of first detection (day 1 = first recording day
# of that site-month). Average across months and sites to obtain a robust
# mean relative first detection day.
#
# Species with a mean close to 1 are detected from the very start of each
# month (likely resident/sedentary). Species with a high mean are detected
# later in average (potentially mobile, arriving on sites during the month).
#
# Input: integrated_results/ CSV files (detections)
# Output: results/first_detection/ CSV only (no figures)

library(tidyverse)
library(lubridate)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/first_detection"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

UTC_OFFSET_HOURS <- -5


#' Calculate mean relative first detection day per species
#'
#' For each month:
#'   - determine day 1 = first recording date at each site that month
#'   - for each species detected at a site that month, compute relative day
#'     of first detection (first detection date - site day 1 + 1)
#' Then average across all month-site combinations per species.
#'
#' @param detections Data frame with file_path, site, common_name.
#' @param species_col Character. Column name for species.
#' @param utc_offset Numeric. UTC offset.
#' @return Data frame with species, mean_relative_day, sd_relative_day,
#'   n_month_sites.
#' @export
calculate_relative_first_detection <- function(detections,
                                                species_col = "common_name",
                                                utc_offset = -5) {
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))
  det$year_month <- format(det$local_date, "%Y-%m")

  # For each site-month, find day 1 (first recording day at that site)
  site_month_start <- det %>%
    group_by(site, year_month) %>%
    summarise(site_month_day1 = min(local_date), .groups = "drop")

  # For each species x site x month, find relative day of first detection
  first_det <- det %>%
    group_by(species = !!sym(species_col), site, year_month) %>%
    summarise(first_detection_date = min(local_date), .groups = "drop") %>%
    left_join(site_month_start, by = c("site", "year_month")) %>%
    mutate(
      relative_day = as.numeric(first_detection_date - site_month_day1) + 1
    )

  # Average across all month-site combos per species
  summary_df <- first_det %>%
    group_by(species) %>%
    summarise(
      mean_relative_day = mean(relative_day),
      sd_relative_day   = sd(relative_day),
      n_month_sites     = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_relative_day)

  list(per_month_site = first_det, summary = summary_df)
}


# -- load data --
cat("Loading detection data...\n")
katydid_det <- read.csv(file.path(INPUT_DIR, "katydid_detections.csv"),
                        stringsAsFactors = FALSE)
bird_det <- read.csv(file.path(INPUT_DIR, "bird_detections.csv"),
                     stringsAsFactors = FALSE)

# -- analysis --
cat("Calculating mean relative first detection day...\n")
katydid_first <- calculate_relative_first_detection(katydid_det,
                                                     utc_offset = UTC_OFFSET_HOURS)
bird_first <- calculate_relative_first_detection(bird_det,
                                                  utc_offset = UTC_OFFSET_HOURS)

# -- save --
write.csv(katydid_first$per_month_site,
          file.path(OUTPUT_DIR, "katydid_relative_first_detection_per_month_site.csv"),
          row.names = FALSE)
write.csv(katydid_first$summary,
          file.path(OUTPUT_DIR, "katydid_mean_relative_first_detection.csv"),
          row.names = FALSE)
write.csv(bird_first$per_month_site,
          file.path(OUTPUT_DIR, "bird_relative_first_detection_per_month_site.csv"),
          row.names = FALSE)
write.csv(bird_first$summary,
          file.path(OUTPUT_DIR, "bird_mean_relative_first_detection.csv"),
          row.names = FALSE)
cat("  Saved first detection CSVs\n")

# -- summary --
cat("\n-- Relative first detection summary --\n")
cat(sprintf("  Katydid: %d species, mean relative day range: %.1f - %.1f\n",
            nrow(katydid_first$summary),
            min(katydid_first$summary$mean_relative_day),
            max(katydid_first$summary$mean_relative_day)))
cat(sprintf("  Bird: %d species, mean relative day range: %.1f - %.1f\n",
            nrow(bird_first$summary),
            min(bird_first$summary$mean_relative_day),
            max(bird_first$summary$mean_relative_day)))
cat("\nFirst detection analysis complete.\n")
