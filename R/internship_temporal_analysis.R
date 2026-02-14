# internship_temporal_analysis.R
# Temporal analysis for M1 IMABEE internship report
# Hourly detection patterns, diel period classification, bird species temporal summary
#
# Reads from: integrated_results/ (pipeline outputs)
# Writes to:  results/internship_m1/temporal/
# Sources:    baseline_helpers.R

# This script is sourced by run_internship_m1_analyses.R which provides `results` and `m1_dir`.

run_internship_temporal <- function(results, m1_dir) {

  cat("\n[INTERNSHIP-M1] Temporal analysis...\n")

  out_dir <- file.path(m1_dir, "temporal")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  UTC_OFFSET_HOURS <- -5

  # ---- 1. Extract temporal info from katydid detections ----
  katydid_det <- results$katydid_data$raw_detections
  if (is.null(katydid_det) || nrow(katydid_det) == 0) {
    katydid_det <- results$katydid_data$detections
  }

  if (is.null(katydid_det) || nrow(katydid_det) == 0) {
    cat("  [!] No katydid detections found, skipping temporal analysis\n")
    return(NULL)
  }

  # Extract timestamp from file_path (YYYYMMDD_HHMMSS pattern)
  pattern <- "(20[0-9]{2}[01][0-9][0-3][0-9]_[0-2][0-9][0-5][0-9][0-5][0-9])"
  katydid_det$timestamp_str <- stringr::str_extract(katydid_det$file_path, pattern)
  katydid_det$utc_time <- lubridate::parse_date_time(
    katydid_det$timestamp_str, "Ymd_HMS", tz = "UTC", quiet = TRUE
  )
  katydid_det$local_time <- katydid_det$utc_time + lubridate::hours(UTC_OFFSET_HOURS)
  katydid_det$hour <- lubridate::hour(katydid_det$local_time)
  katydid_det$local_date <- as.Date(katydid_det$local_time)

  # Remove rows where temporal extraction failed
  temporal_data <- katydid_det %>% dplyr::filter(!is.na(hour))

  if (nrow(temporal_data) == 0) {
    cat("  [!] No temporal information extracted, skipping\n")
    return(NULL)
  }

  cat(sprintf("  Extracted temporal info for %d detections\n", nrow(temporal_data)))

  # ---- 2. Diel period classification ----
  temporal_data <- temporal_data %>%
    dplyr::mutate(
      diel_period = dplyr::case_when(
        hour >= 5  & hour < 9  ~ "Dawn",
        hour >= 9  & hour < 13 ~ "Morning",
        hour >= 13 & hour < 17 ~ "Afternoon",
        hour >= 17 & hour < 21 ~ "Dusk",
        TRUE                   ~ "Night"
      ),
      diel_period = factor(diel_period,
                           levels = c("Dawn", "Morning", "Afternoon", "Dusk", "Night"))
    )

  # ---- 3. Hourly activity statistics (Mean +/- SE per site) ----
  site_hour_counts <- temporal_data %>%
    dplyr::group_by(site, hour) %>%
    dplyr::summarise(detections = dplyr::n(), .groups = "drop")

  hourly_stats <- site_hour_counts %>%
    dplyr::group_by(hour) %>%
    dplyr::summarise(
      n_sites         = dplyr::n(),
      mean_detections = mean(detections),
      sd_detections   = sd(detections),
      se_detections   = sd(detections) / sqrt(dplyr::n()),
      total_detections = sum(detections),
      .groups = "drop"
    ) %>%
    dplyr::arrange(hour)

  # Ensure all 24 hours are present
  all_hours <- data.frame(hour = 0:23)
  hourly_stats <- dplyr::left_join(all_hours, hourly_stats, by = "hour") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                ~ tidyr::replace_na(., 0)))

  # ---- 4. Period summary ----
  period_stats <- temporal_data %>%
    dplyr::group_by(diel_period) %>%
    dplyr::summarise(
      species_count   = dplyr::n_distinct(common_name),
      detection_count = dplyr::n(),
      .groups = "drop"
    )

  # ---- 5. Hourly activity plot (Figure 2 in report) ----
  n_obs <- nrow(site_hour_counts)
  hourly_plot <- ggplot2::ggplot(hourly_stats,
                                 ggplot2::aes(x = hour, y = mean_detections)) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 0.8) +
    ggplot2::geom_point(color = "#0072B2", size = 2.5) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, mean_detections - se_detections),
                   ymax = mean_detections + se_detections),
      width = 0.3, color = "#0072B2", alpha = 0.7
    ) +
    ggplot2::scale_x_continuous(breaks = 0:23) +
    ggplot2::labs(
      title    = expression(paste("Hourly Detection Patterns (Mean ", pm, " SE per site)")),
      subtitle = sprintf("Based on %d site-hour observations", n_obs),
      x = "Hour of Day",
      y = "Mean Number of Detections per Site"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, color = "gray40"),
      axis.text     = ggplot2::element_text(size = 10)
    )

  ggplot2::ggsave(file.path(out_dir, "hourly_activity_pattern.jpg"),
                  hourly_plot, width = 12, height = 6, dpi = 300)
  cat("  Saved: hourly_activity_pattern.jpg\n")

  # ---- 6. Period activity plot ----
  period_plot <- ggplot2::ggplot(period_stats,
                                 ggplot2::aes(x = diel_period, y = detection_count,
                                              fill = diel_period)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::labs(
      title = "Katydid Detections by Diel Period",
      x = "Diel Period (UTC-5)", y = "Total Detections"
    ) +
    ggplot2::theme_minimal()

  ggplot2::ggsave(file.path(out_dir, "diel_period_detections.jpg"),
                  period_plot, width = 8, height = 6, dpi = 300)
  cat("  Saved: diel_period_detections.jpg\n")

  # ---- 7. Bird species temporal summary ----
  bird_det <- results$bird_data$detections
  if (!is.null(bird_det) && nrow(bird_det) > 0 && "file_path" %in% colnames(bird_det)) {
    bird_det$timestamp_str <- stringr::str_extract(bird_det$file_path, pattern)
    bird_det$utc_time <- lubridate::parse_date_time(
      bird_det$timestamp_str, "Ymd_HMS", tz = "UTC", quiet = TRUE
    )
    bird_det$local_time <- bird_det$utc_time + lubridate::hours(UTC_OFFSET_HOURS)
    bird_det$hour <- lubridate::hour(bird_det$local_time)

    bird_hourly <- bird_det %>%
      dplyr::filter(!is.na(hour)) %>%
      dplyr::group_by(hour) %>%
      dplyr::summarise(
        n_species  = dplyr::n_distinct(common_name),
        n_detections = dplyr::n(),
        .groups = "drop"
      )

    write.csv(bird_hourly,
              file.path(out_dir, "bird_hourly_activity.csv"), row.names = FALSE)
    cat("  Saved: bird_hourly_activity.csv\n")
  }

  # ---- 8. Save tables ----
  write.csv(hourly_stats,
            file.path(out_dir, "katydid_hourly_stats.csv"), row.names = FALSE)
  write.csv(period_stats,
            file.path(out_dir, "katydid_diel_period_stats.csv"), row.names = FALSE)

  cat("  Temporal analysis complete\n")

  return(list(
    hourly_activity = hourly_stats,
    period_activity = period_stats,
    hourly_plot     = hourly_plot,
    period_plot     = period_plot
  ))
}
