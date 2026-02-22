# baseline_species_activity.R
# Per-species call activity over time.
#
# For each species, shows the number of detections per month with dry/wet
# season coloring. Complements the seasonal richness barplot (which shows
# species counts) by showing call volume instead.
#
# Analyses:
#   a) Monthly total detections per species (katydids) with dry/wet season
#   b) Faceted barplot: one panel per species, ~30 panels
#   c) Heatmap alternative: species x month, color = detection count (log)
#   d) Same for birds (59 species)
#
# Input: integrated_results/ detection CSV files
# Output: results/species_activity/ CSV and figures

library(tidyverse)
library(lubridate)
library(ggplot2)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR  <- "integrated_results"
OUTPUT_DIR <- "results/species_activity"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

UTC_OFFSET_HOURS <- -5


# --- functions ---------------------------------------------------------------

#' Assign BCI season based on date
#'
#' BCI dry season: mid-December to mid-April.
#' Simplified: months 1-4 and 12 = dry, 5-11 = wet.
#'
#' @param dates Date vector.
#' @return Character vector of "Dry" or "Wet".
#' @export
assign_bci_season <- function(dates) {
  m <- month(dates)
  ifelse(m %in% c(1, 2, 3, 4, 12), "Dry", "Wet")
}


#' Calculate monthly detection counts per species
#'
#' @param detections Data frame with file_path, common_name, site.
#' @param species_col Character.
#' @param utc_offset Numeric.
#' @return Data frame with species, year_month, month_label, n_detections,
#'   season, and total_detections_all_months.
#' @export
calculate_monthly_detections <- function(detections, species_col = "common_name",
                                          utc_offset = -5) {
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))

  det$year_month <- floor_date(det$local_date, "month")

  monthly <- det %>%
    group_by(species = !!sym(species_col), year_month) %>%
    summarise(n_detections = n(), .groups = "drop") %>%
    mutate(
      month_label = format(year_month, "%Y-%m"),
      season = assign_bci_season(year_month)
    )

  # Add total detections for ordering
  totals <- monthly %>%
    group_by(species) %>%
    summarise(total_detections = sum(n_detections), .groups = "drop")
  monthly <- monthly %>% left_join(totals, by = "species")

  monthly %>% arrange(desc(total_detections), species, year_month)
}


#' Faceted barplot of monthly detections per species
#'
#' One panel per species, colored by season. Ordered by total detections
#' (most active species first). Uses log10 scale if scale_log = TRUE.
#'
#' @param monthly_det Data frame from calculate_monthly_detections.
#' @param taxon_label Character.
#' @param fill_wet Character. Color for wet season.
#' @param fill_dry Character. Color for dry season.
#' @param scale_log Logical.
#' @param ncol Integer. Number of columns in facet grid.
#' @return ggplot object.
#' @export
plot_species_activity_faceted <- function(monthly_det, taxon_label = "Katydid",
                                          fill_wet = "#009E73",
                                          fill_dry = "#E69F00",
                                          scale_log = FALSE, ncol = 5) {
  # Order species by total detections
  sp_order <- monthly_det %>%
    distinct(species, total_detections) %>%
    arrange(desc(total_detections)) %>%
    pull(species)
  monthly_det <- monthly_det %>%
    mutate(species = factor(species, levels = sp_order))

  # Ensure all months are present for each species
  all_months <- sort(unique(monthly_det$month_label))
  complete_grid <- expand.grid(species = sp_order, month_label = all_months,
                                stringsAsFactors = FALSE) %>%
    mutate(species = factor(species, levels = sp_order))
  monthly_det <- complete_grid %>%
    left_join(monthly_det %>% select(-year_month),
              by = c("species", "month_label")) %>%
    mutate(
      n_detections = replace_na(n_detections, 0),
      season = ifelse(is.na(season),
                      assign_bci_season(as.Date(paste0(month_label, "-15"))),
                      season)
    )

  p <- ggplot(monthly_det, aes(x = month_label, y = n_detections, fill = season)) +
    geom_col(width = 0.8, alpha = 0.85) +
    scale_fill_manual(values = c("Dry" = fill_dry, "Wet" = fill_wet),
                      name = "Season") +
    facet_wrap(~ species, scales = "free_y", ncol = ncol) +
    labs(x = "Month", y = "Number of detections",
         title = paste(taxon_label, "monthly call activity by species"),
         subtitle = "BCI dry season: mid-Dec to mid-Apr") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
          strip.text = element_text(size = 7, face = "bold"),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "bottom",
          panel.spacing = unit(0.3, "lines"))

  if (scale_log) {
    p <- p + scale_y_log10(labels = scales::comma_format()) +
      labs(y = "Number of detections (log scale)")
  }
  p
}


#' Heatmap of species x month detection counts
#'
#' Color represents log10(detections + 1). Species ordered by total detections.
#'
#' @param monthly_det Data frame from calculate_monthly_detections.
#' @param taxon_label Character.
#' @param color_low Character.
#' @param color_high Character.
#' @return ggplot object.
#' @export
plot_species_activity_heatmap <- function(monthly_det, taxon_label = "Katydid",
                                           color_low = "#f7fcf5",
                                           color_high = "#00441b") {
  # Order species by total detections (most active at top)
  sp_order <- monthly_det %>%
    distinct(species, total_detections) %>%
    arrange(total_detections) %>%  # ascending so most active is at top of y
    pull(species)

  all_months <- sort(unique(monthly_det$month_label))

  # Complete grid
  complete_grid <- expand.grid(species = sp_order, month_label = all_months,
                                stringsAsFactors = FALSE)
  heatmap_data <- complete_grid %>%
    left_join(monthly_det %>% select(species, month_label, n_detections, season),
              by = c("species", "month_label")) %>%
    mutate(
      n_detections = replace_na(n_detections, 0),
      log_detections = log10(n_detections + 1),
      species = factor(species, levels = sp_order),
      season = ifelse(is.na(season),
                      assign_bci_season(as.Date(paste0(month_label, "-15"))),
                      season)
    )

  # Season indicator for x-axis
  month_seasons <- heatmap_data %>%
    distinct(month_label, season) %>%
    arrange(month_label)
  x_labels <- paste0(month_seasons$month_label,
                      ifelse(month_seasons$season == "Dry", " *", ""))

  ggplot(heatmap_data, aes(x = month_label, y = species, fill = log_detections)) +
    geom_tile(color = "grey80", linewidth = 0.2) +
    scale_fill_gradient(low = color_low, high = color_high,
                        name = expression(log[10](detections + 1))) +
    scale_x_discrete(labels = x_labels) +
    labs(x = "Month (* = dry season)", y = NULL,
         title = paste(taxon_label, "call activity heatmap (species x month)"),
         subtitle = "Color intensity = log10(detection count + 1)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "right")
}


# --- execute -----------------------------------------------------------------

cat("Loading detection data...\n")
katydid_det <- read.csv(file.path(INPUT_DIR, "katydid_detections.csv"),
                        stringsAsFactors = FALSE)
bird_det    <- read.csv(file.path(INPUT_DIR, "bird_detections.csv"),
                        stringsAsFactors = FALSE)

# -- Katydid monthly detections per species --
cat("[>>] Calculating katydid monthly detections per species...\n")
katydid_monthly <- calculate_monthly_detections(katydid_det,
                                                 utc_offset = UTC_OFFSET_HOURS)
write.csv(katydid_monthly,
          file.path(OUTPUT_DIR, "katydid_monthly_detections_per_species.csv"),
          row.names = FALSE)
n_katydid_sp <- n_distinct(katydid_monthly$species)
cat(sprintf("[OK] %d katydid species, %s total detections\n",
            n_katydid_sp,
            format(sum(katydid_monthly$n_detections), big.mark = ",")))

# Faceted barplot (linear scale)
cat("[>>] Creating katydid faceted activity plot...\n")
p1 <- plot_species_activity_faceted(katydid_monthly, "Katydid", ncol = 5)
height1 <- max(6, ceiling(n_katydid_sp / 5) * 2.5)
ggsave(file.path(FIGURE_DIR, "katydid_species_activity_faceted.jpg"), p1,
       width = 14, height = height1, dpi = 300, limitsize = FALSE)

# Faceted barplot (log scale)
p1_log <- plot_species_activity_faceted(katydid_monthly, "Katydid",
                                         scale_log = TRUE, ncol = 5)
ggsave(file.path(FIGURE_DIR, "katydid_species_activity_faceted_log.jpg"), p1_log,
       width = 14, height = height1, dpi = 300, limitsize = FALSE)

# Heatmap
cat("[>>] Creating katydid activity heatmap...\n")
p2 <- plot_species_activity_heatmap(katydid_monthly, "Katydid")
ggsave(file.path(FIGURE_DIR, "katydid_species_activity_heatmap.jpg"), p2,
       width = 10, height = max(6, n_katydid_sp * 0.35), dpi = 300,
       limitsize = FALSE)

# -- Bird monthly detections per species --
cat("[>>] Calculating bird monthly detections per species...\n")
bird_monthly <- calculate_monthly_detections(bird_det, utc_offset = UTC_OFFSET_HOURS)
write.csv(bird_monthly,
          file.path(OUTPUT_DIR, "bird_monthly_detections_per_species.csv"),
          row.names = FALSE)
n_bird_sp <- n_distinct(bird_monthly$species)
cat(sprintf("[OK] %d bird species, %s total detections\n",
            n_bird_sp,
            format(sum(bird_monthly$n_detections), big.mark = ",")))

# Faceted barplot (linear scale)
cat("[>>] Creating bird faceted activity plot...\n")
p3 <- plot_species_activity_faceted(bird_monthly, "Bird",
                                     fill_wet = "#0072B2", ncol = 6)
height3 <- max(6, ceiling(n_bird_sp / 6) * 2.5)
ggsave(file.path(FIGURE_DIR, "bird_species_activity_faceted.jpg"), p3,
       width = 16, height = height3, dpi = 300, limitsize = FALSE)

# Heatmap
cat("[>>] Creating bird activity heatmap...\n")
p4 <- plot_species_activity_heatmap(bird_monthly, "Bird",
                                     color_low = "#f7fbff", color_high = "#08306b")
ggsave(file.path(FIGURE_DIR, "bird_species_activity_heatmap.jpg"), p4,
       width = 10, height = max(6, n_bird_sp * 0.3), dpi = 300,
       limitsize = FALSE)

# -- Summary --
cat("\n-- Species activity summary --\n")

# Top 5 most active katydid species
top_katydid <- katydid_monthly %>%
  distinct(species, total_detections) %>%
  arrange(desc(total_detections)) %>%
  head(5)
cat("  Top 5 most active katydid species:\n")
for (i in seq_len(nrow(top_katydid))) {
  cat(sprintf("    %d. %s (%s detections)\n", i, top_katydid$species[i],
              format(top_katydid$total_detections[i], big.mark = ",")))
}

# Top 5 most active bird species
top_bird <- bird_monthly %>%
  distinct(species, total_detections) %>%
  arrange(desc(total_detections)) %>%
  head(5)
cat("  Top 5 most active bird species:\n")
for (i in seq_len(nrow(top_bird))) {
  cat(sprintf("    %d. %s (%s detections)\n", i, top_bird$species[i],
              format(top_bird$total_detections[i], big.mark = ",")))
}

cat("\nSpecies activity analysis complete.\n")
