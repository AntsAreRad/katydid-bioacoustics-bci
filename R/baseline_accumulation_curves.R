# baseline_accumulation_curves.R
# Accumulation curves:
#   a) Global alpha-beta-gamma (single recorder vs composite, full dataset)
#   b) Monthly accumulation curves (one per month, compare saturation speed)
#   c) Seasonal richness barplot (cumulative species at day 7 per month)
#
# Input: integrated_results/ CSV files (detections with timestamps)
# Output: results/accumulation/ CSV and figures

library(tidyverse)
library(lubridate)
library(ggplot2)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/accumulation"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

UTC_OFFSET_HOURS <- -5
MAX_DAYS <- 16
N_PERMUTATIONS <- 100
SEED <- 42


#' Calculate global species accumulation with alpha-beta-gamma decomposition
#'
#' @param detections Data frame with file_path, site, common_name.
#' @param species_col Character. Species column name.
#' @param utc_offset Numeric. UTC offset.
#' @param max_days Integer. Max days for accumulation.
#' @param n_permutations Integer. Permutation count.
#' @param seed Integer. Random seed.
#' @return List with single_recorder, composite, alpha, beta, gamma, n_sites.
#' @export
calculate_accumulation_abg <- function(detections, species_col = "common_name",
                                       utc_offset = -5, max_days = 16,
                                       n_permutations = 100, seed = 42) {
  set.seed(seed)
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))
  sites <- unique(det$site)
  n_sites <- length(sites)

  # single recorder: for each permutation x site, accumulate species
  single_results <- matrix(NA, nrow = n_permutations * n_sites, ncol = max_days)
  row_idx <- 1
  for (perm in seq_len(n_permutations)) {
    for (s in sites) {
      site_det <- det %>% filter(site == s)
      dates <- sort(unique(site_det$local_date))
      if (length(dates) == 0) next
      date_order <- sample(dates)
      cum_sp <- character()
      for (d in seq_len(min(max_days, length(date_order)))) {
        day_sp <- site_det %>%
          filter(local_date == date_order[d]) %>%
          pull(!!sym(species_col)) %>% unique()
        cum_sp <- union(cum_sp, day_sp)
        single_results[row_idx, d] <- length(cum_sp)
      }
      if (length(date_order) < max_days) {
        last_val <- length(cum_sp)
        for (d in (length(date_order) + 1):max_days)
          single_results[row_idx, d] <- last_val
      }
      row_idx <- row_idx + 1
    }
  }
  single_results <- single_results[seq_len(row_idx - 1), , drop = FALSE]
  single_mean <- colMeans(single_results, na.rm = TRUE)
  single_se <- apply(single_results, 2, sd, na.rm = TRUE) /
    sqrt(nrow(single_results))

  # composite: pool all recorders
  composite_results <- matrix(NA, nrow = n_permutations, ncol = max_days)
  for (perm in seq_len(n_permutations)) {
    dates_all <- sort(unique(det$local_date))
    date_order <- sample(dates_all)
    cum_sp <- character()
    for (d in seq_len(min(max_days, length(date_order)))) {
      day_sp <- det %>%
        filter(local_date == date_order[d]) %>%
        pull(!!sym(species_col)) %>% unique()
      cum_sp <- union(cum_sp, day_sp)
      composite_results[perm, d] <- length(cum_sp)
    }
    if (length(date_order) < max_days) {
      last_val <- length(cum_sp)
      for (d in (length(date_order) + 1):max_days)
        composite_results[perm, d] <- last_val
    }
  }
  composite_mean <- colMeans(composite_results, na.rm = TRUE)
  composite_se <- apply(composite_results, 2, sd, na.rm = TRUE) /
    sqrt(n_permutations)

  alpha <- single_mean[1]
  beta_temporal <- if (max_days >= 10) {
    (single_mean[10] - single_mean[1]) / 9
  } else {
    (single_mean[max_days] - single_mean[1]) / (max_days - 1)
  }
  gamma <- composite_mean[max_days] - single_mean[max_days]

  list(
    single_recorder = data.frame(day = seq_len(max_days),
                                 mean_species = single_mean, se = single_se),
    composite = data.frame(day = seq_len(max_days),
                           mean_species = composite_mean, se = composite_se),
    alpha_diversity = alpha, beta_temporal = beta_temporal,
    gamma_diversity = gamma, n_sites = n_sites
  )
}


#' Plot global ABG accumulation with distinct colors
#'
#' Single recorder curve is colored (green for katydids, blue for birds),
#' composite curve is black.  SE envelopes are semi-transparent.
#'
#' @param accum_result List from calculate_accumulation_abg.
#' @param taxon_name Character.
#' @param colors Named vector with single and composite.
#' @return ggplot object.
#' @export
plot_accumulation_abg <- function(accum_result, taxon_name = "katydid",
                                  colors = c(single = "#009E73",
                                             composite = "#000000")) {
  df <- bind_rows(
    accum_result$single_recorder %>% mutate(type = "Single recorder (mean)"),
    accum_result$composite %>% mutate(type = "Composite (25 recorders)")
  )
  annotation <- sprintf(
    "alpha = %.1f spp (day 1)\nbeta temporal = %.2f spp/day\ngamma = +%.1f spp (spatial)",
    accum_result$alpha_diversity, accum_result$beta_temporal,
    accum_result$gamma_diversity)

  ggplot(df, aes(x = day, y = mean_species, color = type, fill = type)) +
    geom_ribbon(aes(ymin = mean_species - se, ymax = mean_species + se),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 1.5) +
    scale_color_manual(values = c("Single recorder (mean)" = unname(colors["single"]),
                                  "Composite (25 recorders)" = unname(colors["composite"]))) +
    scale_fill_manual(values = c("Single recorder (mean)" = unname(colors["single"]),
                                 "Composite (25 recorders)" = unname(colors["composite"]))) +
    annotate("text", x = max(df$day) * 0.55, y = min(df$mean_species) + 1,
             label = annotation, hjust = 0, size = 3.5) +
    labs(x = "Number of recording days",
         y = "Cumulative species detected",
         title = paste0(tools::toTitleCase(taxon_name),
                        " species accumulation (alpha-beta-gamma)"),
         color = NULL, fill = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 13, face = "bold"))
}


#' Calculate monthly accumulation curves (no fill-forward)
#'
#' Curves are truncated to the actual number of recording days available
#' across sites for that month (median number of unique dates per site).
#' This avoids artificial flat lines from fill-forward when only a few
#' days are available.
#'
#' @param detections Data frame with file_path, site, common_name.
#' @param species_col Character.
#' @param utc_offset Numeric.
#' @param max_days Integer. Hard ceiling.
#' @param n_permutations Integer.
#' @param seed Integer.
#' @return Data frame with month_label, day, mean_species, se, effective_days.
#' @export
calculate_monthly_accumulation <- function(detections, species_col = "common_name",
                                            utc_offset = -5, max_days = 16,
                                            n_permutations = 50, seed = 42) {
  set.seed(seed)
  det <- detections
  det$local_date <- extract_local_date(det$file_path, utc_offset)
  det <- det %>% filter(!is.na(local_date))
  det$year_month <- floor_date(det$local_date, "month")
  months_present <- sort(unique(det$year_month))
  all_monthly <- list()

  for (ym in months_present) {
    month_det <- det %>% filter(year_month == ym)
    sites <- unique(month_det$site)
    label <- format(as.Date(ym), "%Y-%m")

    # Determine effective number of days: median unique dates per site
    days_per_site <- month_det %>%
      group_by(site) %>%
      summarise(n_dates = n_distinct(local_date), .groups = "drop")
    effective_days <- min(max_days,
                          as.integer(median(days_per_site$n_dates)))
    if (effective_days < 2) next

    results_mat <- matrix(NA, nrow = n_permutations * length(sites),
                          ncol = effective_days)
    row_idx <- 1
    for (perm in seq_len(n_permutations)) {
      for (s in sites) {
        site_det <- month_det %>% filter(site == s)
        dates <- sort(unique(site_det$local_date))
        if (length(dates) == 0) next
        date_order <- sample(dates)
        cum_sp <- character()
        n_use <- min(effective_days, length(date_order))
        for (d in seq_len(n_use)) {
          day_sp <- site_det %>%
            filter(local_date == date_order[d]) %>%
            pull(!!sym(species_col)) %>% unique()
          cum_sp <- union(cum_sp, day_sp)
          results_mat[row_idx, d] <- length(cum_sp)
        }
        # NO fill-forward: leave remaining columns as NA for sites with
        # fewer dates than effective_days
        row_idx <- row_idx + 1
      }
    }
    results_mat <- results_mat[seq_len(row_idx - 1), , drop = FALSE]
    if (nrow(results_mat) == 0) next
    m_mean <- colMeans(results_mat, na.rm = TRUE)
    m_se <- apply(results_mat, 2, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(0)
      sd(x) / sqrt(length(x))
    })
    all_monthly[[label]] <- data.frame(
      month_label = label, day = seq_len(effective_days),
      mean_species = m_mean, se = m_se, effective_days = effective_days)
  }
  bind_rows(all_monthly)
}


#' Plot monthly accumulation curves (spaghetti plot)
#' @param monthly_df Data frame from calculate_monthly_accumulation.
#' @param taxon_name Character.
#' @return ggplot object.
#' @export
plot_monthly_accumulation <- function(monthly_df, taxon_name = "katydid") {
  ggplot(monthly_df, aes(x = day, y = mean_species,
                         color = month_label, group = month_label)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    labs(x = "Number of recording days",
         y = "Cumulative species detected (single recorder mean)",
         title = paste0(tools::toTitleCase(taxon_name),
                        " monthly accumulation curves"),
         color = "Month") +
    theme_bw(base_size = 12) +
    theme(legend.position = "right",
          plot.title = element_text(size = 13, face = "bold"))
}


#' Plot seasonal richness: cumulative species at a fixed day per month
#'
#' Shows a point+bar plot of mean cumulative species at a reference day
#' for each month, ordered chronologically. Allows visual comparison of
#' seasonal patterns (dry season: mid-Dec to mid-Apr at BCI).
#'
#' @param monthly_df Data frame from calculate_monthly_accumulation.
#' @param taxon_name Character.
#' @param ref_day Integer. Day to extract richness at (default 7).
#' @param fill_color Character. Bar colour.
#' @return ggplot object.
#' @export
plot_seasonal_richness <- function(monthly_df, taxon_name = "katydid",
                                    ref_day = 7, fill_color = "#009E73") {
  # Extract richness at the reference day (or the last available day)
  seasonal <- monthly_df %>%
    group_by(month_label) %>%
    summarise(
      richness_at_ref = {
        if (any(day == ref_day)) {
          mean_species[day == ref_day]
        } else {
          mean_species[day == max(day)]
        }
      },
      se_at_ref = {
        if (any(day == ref_day)) {
          se[day == ref_day]
        } else {
          se[day == max(day)]
        }
      },
      effective_days = max(effective_days),
      day_used = min(ref_day, max(day)),
      .groups = "drop"
    ) %>%
    arrange(month_label) %>%
    mutate(month_label = factor(month_label, levels = month_label))

  # Annotate dry season months (mid-Dec to mid-Apr -> Jan, Feb, Mar, Apr)
  dry_months <- c("01", "02", "03", "04", "12")
  seasonal <- seasonal %>%
    mutate(
      month_num = substr(month_label, 6, 7),
      season = ifelse(month_num %in% dry_months, "Dry", "Wet")
    )

  ggplot(seasonal, aes(x = month_label, y = richness_at_ref, fill = season)) +
    geom_col(width = 0.7, alpha = 0.85) +
    geom_errorbar(aes(ymin = richness_at_ref - se_at_ref,
                      ymax = richness_at_ref + se_at_ref),
                  width = 0.3, linewidth = 0.4) +
    geom_point(size = 2, color = "black") +
    scale_fill_manual(values = c("Dry" = "#E69F00", "Wet" = fill_color),
                      name = "Season") +
    labs(x = "Month",
         y = paste0("Cumulative species at day ", ref_day,
                    " (single recorder mean)"),
         title = paste0(tools::toTitleCase(taxon_name),
                        " seasonal richness pattern"),
         subtitle = paste0("BCI dry season: mid-Dec to mid-Apr")) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 13, face = "bold"),
          legend.position = "bottom")
}


# -- execute --
cat("Loading detection data...\n")
katydid_det <- read.csv(file.path(INPUT_DIR, "katydid_detections.csv"),
                        stringsAsFactors = FALSE)
bird_det <- read.csv(file.path(INPUT_DIR, "bird_detections.csv"),
                     stringsAsFactors = FALSE)

cat("Calculating global katydid accumulation...\n")
katydid_accum <- calculate_accumulation_abg(katydid_det, max_days = MAX_DAYS,
                                            n_permutations = N_PERMUTATIONS, seed = SEED)
cat("Calculating global bird accumulation...\n")
bird_accum <- calculate_accumulation_abg(bird_det, max_days = MAX_DAYS,
                                         n_permutations = N_PERMUTATIONS, seed = SEED)

for (taxon in c("katydid", "bird")) {
  acc <- if (taxon == "katydid") katydid_accum else bird_accum
  df <- data.frame(day = acc$single_recorder$day,
                   single_mean = acc$single_recorder$mean_species,
                   single_se = acc$single_recorder$se,
                   composite_mean = acc$composite$mean_species,
                   composite_se = acc$composite$se)
  write.csv(df, file.path(OUTPUT_DIR, paste0(taxon, "_global_accumulation.csv")),
            row.names = FALSE)
  write.csv(data.frame(taxon = taxon, alpha = acc$alpha_diversity,
                       beta_temporal = acc$beta_temporal,
                       gamma = acc$gamma_diversity, n_sites = acc$n_sites),
            file.path(OUTPUT_DIR, paste0(taxon, "_abg_summary.csv")),
            row.names = FALSE)
}

p1 <- plot_accumulation_abg(katydid_accum, "katydid",
                            c(single = "#009E73", composite = "#000000"))
ggsave(file.path(FIGURE_DIR, "katydid_accumulation_abg.jpg"), p1,
       width = 9, height = 6, dpi = 300)
p2 <- plot_accumulation_abg(bird_accum, "bird",
                            c(single = "#0072B2", composite = "#000000"))
ggsave(file.path(FIGURE_DIR, "bird_accumulation_abg.jpg"), p2,
       width = 9, height = 6, dpi = 300)

cat("Calculating monthly katydid accumulation...\n")
katydid_monthly <- calculate_monthly_accumulation(katydid_det, max_days = MAX_DAYS,
                                                   n_permutations = 50, seed = SEED)
cat("Calculating monthly bird accumulation...\n")
bird_monthly <- calculate_monthly_accumulation(bird_det, max_days = MAX_DAYS,
                                                n_permutations = 50, seed = SEED)

write.csv(katydid_monthly, file.path(OUTPUT_DIR, "katydid_monthly_accumulation.csv"),
          row.names = FALSE)
write.csv(bird_monthly, file.path(OUTPUT_DIR, "bird_monthly_accumulation.csv"),
          row.names = FALSE)

p3 <- plot_monthly_accumulation(katydid_monthly, "katydid")
ggsave(file.path(FIGURE_DIR, "katydid_monthly_accumulation.jpg"), p3,
       width = 10, height = 6, dpi = 300)
p4 <- plot_monthly_accumulation(bird_monthly, "bird")
ggsave(file.path(FIGURE_DIR, "bird_monthly_accumulation.jpg"), p4,
       width = 10, height = 6, dpi = 300)

# Seasonal richness barplots (task 4)
p5 <- plot_seasonal_richness(katydid_monthly, "katydid", ref_day = 7,
                              fill_color = "#009E73")
ggsave(file.path(FIGURE_DIR, "katydid_seasonal_richness.jpg"), p5,
       width = 9, height = 6, dpi = 300)
p6 <- plot_seasonal_richness(bird_monthly, "bird", ref_day = 7,
                              fill_color = "#0072B2")
ggsave(file.path(FIGURE_DIR, "bird_seasonal_richness.jpg"), p6,
       width = 9, height = 6, dpi = 300)

cat("\n-- Accumulation summary --\n")
cat(sprintf("  Katydid: alpha=%.1f, beta=%.2f spp/day, gamma=+%.1f\n",
            katydid_accum$alpha_diversity, katydid_accum$beta_temporal,
            katydid_accum$gamma_diversity))
cat(sprintf("  Bird: alpha=%.1f, beta=%.2f spp/day, gamma=+%.1f\n",
            bird_accum$alpha_diversity, bird_accum$beta_temporal,
            bird_accum$gamma_diversity))
cat(sprintf("  Monthly: %d months (katydid), %d months (bird)\n",
            n_distinct(katydid_monthly$month_label),
            n_distinct(bird_monthly$month_label)))
cat("\nAccumulation analysis complete.\n")
