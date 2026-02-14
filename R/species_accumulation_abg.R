# ==============================================================================
# SPECIES ACCUMULATION WITH ALPHA-BETA-GAMMA DIVERSITY DECOMPOSITION
# ==============================================================================
# Description: Functions for calculating and visualizing species accumulation 
#              curves with temporal (beta) and spatial (gamma) diversity 
#              decomposition, inspired by Laurel Symes' Ecology & Evolution 
#              visualizations.
#
# Author: Leon Brouille
# Project: Katydid Bioacoustics - Biology Letters Submission
# Institution: STRI (Smithsonian Tropical Research Institute) & 
#              Cornell Lab of Ornithology
# Date: February 2026
# ==============================================================================

# Required packages
# library(tidyverse)
# library(ggplot2)

# ==============================================================================
# SECTION 1: SPECIES ACCUMULATION CALCULATION
# ==============================================================================

#' Calculate Species Accumulation with Alpha-Beta-Gamma Decomposition
#'
#' @description
#' Calculates species accumulation curves for single recorders and composite 
#' multi-recorder sampling, enabling decomposition into alpha (point diversity),
#' temporal beta (within-recorder accumulation), and gamma (spatial turnover)
#' diversity components.
#'
#' @details
#' The function implements the following diversity framework:
#' \itemize{
#'   \item \strong{Alpha-diversity}: Expected species detected from one recorder
#'         on one day (initial point on single-recorder curve)
#'   \item \strong{Temporal beta-diversity}: Rate at which new species are 
#'         detected as more days are analyzed from the same recorder
#'   \item \strong{Gamma-diversity}: Additional species detected when using
#'         multiple recorders for the same total recorder-days, representing
#'         spatial turnover between sites
#' }
#'
#' The composite curve is computed by pooling site-days from n recorders and
#' randomly sampling to create an "equivalent effort" comparison. This allows
#' direct visual comparison: at X recorder-days, how many more species do you
#' detect with 10 recorders vs 1?
#'
#' @param detections_data Data frame with acoustic detections. Required columns:
#'   \itemize{
#'     \item \code{file_path}: Character. File path containing timestamp 
#'           (format YYYYMMDD_HHMMSSZ)
#'     \item \code{site}: Character. Site identifier
#'     \item \code{confidence}: Numeric. Detection confidence (0-1), optional
#'   }
#' @param species_col Character. Column name containing species identifiers.
#'   Default: "common_name"
#' @param confidence_threshold Numeric. Minimum confidence for including
#'   detections. Default: 0.9
#' @param n_permutations Integer. Number of permutations for variance estimation.
#'   Default: 100
#' @param n_sites_composite Integer. Number of sites for composite curve.
#'   Default: 10
#' @param max_days Integer. Maximum days to include in accumulation. 
#'   Use 7 for weekly analysis. Default: NULL (use all available)
#' @param seed Integer. Random seed for reproducibility. Default: 42
#'
#' @return A list of class "species_accumulation_abg" containing:
#'   \describe{
#'     \item{single_recorder}{Data frame with single-recorder accumulation:
#'       day, mean_species, se, sd, lower_ci, upper_ci}
#'     \item{composite}{Data frame with composite accumulation, same structure}
#'     \item{site_curves}{List of per-site accumulation curves}
#'     \item{alpha_diversity}{Numeric. Mean species on day 1}
#'     \item{beta_temporal}{Numeric. Slope of single-recorder curve (species/day)}
#'     \item{gamma_diversity}{Numeric. Composite minus single at max days}
#'     \item{parameters}{List of input parameters}
#'     \item{n_sites_used}{Integer. Number of sites in analysis}
#'     \item{max_days_used}{Integer. Number of days analyzed}
#'     \item{total_species}{Integer. Total unique species detected}
#'   }
#'
#' @references
#' Symes, L.B., et al. (2022). Ecology and Evolution.
#' Colwell, R.K., et al. (2012). Models and estimators linking individual-based
#'   and sample-based rarefaction. Journal of Plant Ecology, 5(1), 3-21.
#'
#' @examples
#' \dontrun{
#' # Calculate accumulation for katydids
#' katydid_accum <- calculate_accumulation_abg(
#'   detections_data = katydid_detections,
#'   species_col = "common_name",
#'   confidence_threshold = 0.9,
#'   n_sites_composite = 10,
#'   max_days = 16
#' )
#'
#' # View diversity metrics
#' cat("Alpha:", katydid_accum$alpha_diversity, "\n")
#' cat("Temporal beta:", katydid_accum$beta_temporal, "\n")
#' cat("Gamma:", katydid_accum$gamma_diversity, "\n")
#' }
#'
#' @export
calculate_accumulation_abg <- function(detections_data,
                                       species_col = "common_name",
                                       confidence_threshold = 0.95,
                                       n_permutations = 100,
                                       n_sites_composite = 25,
                                       max_days = NULL,
                                       seed = 42) {
  
  set.seed(seed)
  
  cat("SPECIES ACCUMULATION (Alpha-Beta-Gamma)\n\n")
  
  # Input validation
  required_cols <- c(species_col, "site", "file_path")
  missing_cols <- setdiff(required_cols, names(detections_data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Apply confidence filter
  if ("confidence" %in% names(detections_data)) {
    n_before <- nrow(detections_data)
    detections_data <- detections_data %>%
      filter(confidence >= confidence_threshold)
    cat(sprintf("Confidence filter (>= %.2f): %d -> %d detections\n",
                confidence_threshold, n_before, nrow(detections_data)))
  }
  
  # Extract date from file_path
  detections_data <- detections_data %>%
    mutate(
      timestamp_str = str_extract(file_path, "\\d{8}_\\d{6}Z?"),
      detection_date = as.Date(substr(timestamp_str, 1, 8), format = "%Y%m%d")
    ) %>%
    filter(!is.na(detection_date))
  
  if (nrow(detections_data) == 0) {
    stop("No valid detections with extractable timestamps")
  }
  
  # Get site information
  all_sites <- unique(detections_data$site)
  n_all_sites <- length(all_sites)
  cat(sprintf("Sites available: %d\n", n_all_sites))
  
  # Select sites for composite
  if (n_sites_composite >= n_all_sites) {
    sites_composite <- all_sites
  } else {
    sites_composite <- sample(all_sites, n_sites_composite)
  }
  n_sites_used <- length(sites_composite)
  cat(sprintf("Sites for composite: %d\n", n_sites_used))
  
  # Calculate per-site accumulation curves
  cat("Calculating per-site curves...\n")
  site_curves <- list()
  
  for (site_id in all_sites) {
    site_data <- detections_data %>% filter(site == site_id)
    unique_dates <- sort(unique(site_data$detection_date))
    n_days <- length(unique_dates)
    
    if (n_days < 2) next
    
    # Permutation-based accumulation
    perm_matrix <- matrix(NA, nrow = n_days, ncol = n_permutations)
    
    for (perm in 1:n_permutations) {
      date_order <- sample(unique_dates)
      cumul_species <- character(0)
      
      for (d in 1:n_days) {
        day_species <- site_data %>%
          filter(detection_date == date_order[d]) %>%
          pull(!!sym(species_col)) %>%
          unique()
        cumul_species <- unique(c(cumul_species, day_species))
        perm_matrix[d, perm] <- length(cumul_species)
      }
    }
    
    site_curves[[site_id]] <- data.frame(
      day = 1:n_days,
      mean_species = rowMeans(perm_matrix),
      se = apply(perm_matrix, 1, function(x) sd(x) / sqrt(length(x))),
      sd = apply(perm_matrix, 1, sd)
    )
  }
  
  if (length(site_curves) == 0) {
    stop("No sites with sufficient data (>= 2 days)")
  }
  
  # Determine common max days
  max_days_available <- min(sapply(site_curves, nrow))
  if (!is.null(max_days)) {
    max_days_used <- min(max_days_available, max_days)
  } else {
    max_days_used <- max_days_available
  }
  cat(sprintf("Days for analysis: %d\n", max_days_used))
  
  # Average single-recorder curve
  species_matrix <- matrix(NA, nrow = max_days_used, ncol = length(site_curves))
  for (i in seq_along(site_curves)) {
    species_matrix[, i] <- site_curves[[i]]$mean_species[1:max_days_used]
  }
  
  single_recorder <- data.frame(
    day = 1:max_days_used,
    mean_species = rowMeans(species_matrix, na.rm = TRUE),
    se = apply(species_matrix, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    sd = apply(species_matrix, 1, function(x) sd(x, na.rm = TRUE))
  ) %>%
    mutate(
      lower_ci = mean_species - 1.96 * se,
      upper_ci = mean_species + 1.96 * se
    )
  
  # Composite curve (multiple sites pooled)
  cat("Calculating composite curve...\n")
  composite_data <- detections_data %>% filter(site %in% sites_composite)
  
  site_days <- composite_data %>%
    distinct(site, detection_date) %>%
    arrange(detection_date, site)
  
  n_site_days <- nrow(site_days)
  cat(sprintf("Total site-days for composite: %d\n", n_site_days))
  
  # Composite permutations
  composite_perm <- matrix(NA, nrow = max_days_used, ncol = n_permutations)
  
  for (perm in 1:n_permutations) {
    shuffled_sd <- site_days[sample(1:n_site_days), ]
    cumul_species <- character(0)
    
    for (d in 1:max_days_used) {
      # Sample n_sites_used site-days per "equivalent day"
      start_idx <- (d - 1) * n_sites_used + 1
      end_idx <- min(d * n_sites_used, n_site_days)
      
      if (start_idx > n_site_days) break
      
      selected_sd <- shuffled_sd[start_idx:end_idx, ]
      
      for (j in 1:nrow(selected_sd)) {
        day_species <- composite_data %>%
          filter(site == selected_sd$site[j],
                 detection_date == selected_sd$detection_date[j]) %>%
          pull(!!sym(species_col)) %>%
          unique()
        cumul_species <- unique(c(cumul_species, day_species))
      }
      
      composite_perm[d, perm] <- length(cumul_species)
    }
  }
  
  composite <- data.frame(
    day = 1:max_days_used,
    mean_species = rowMeans(composite_perm, na.rm = TRUE),
    se = apply(composite_perm, 1, function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
    sd = apply(composite_perm, 1, function(x) sd(x, na.rm = TRUE))
  ) %>%
    mutate(
      lower_ci = mean_species - 1.96 * se,
      upper_ci = mean_species + 1.96 * se
    )
  
  # Calculate diversity metrics
  alpha_diversity <- single_recorder$mean_species[1]
  
  # Temporal beta: slope of single-recorder curve (first 5 days)
  fit_days <- min(5, max_days_used)
  if (fit_days >= 2) {
    lm_fit <- lm(mean_species ~ day, data = single_recorder[1:fit_days, ])
    beta_temporal <- as.numeric(coef(lm_fit)[2])
  } else {
    beta_temporal <- NA
  }
  
  # Gamma: difference between composite and single at max day
  gamma_diversity <- composite$mean_species[max_days_used] - 
                     single_recorder$mean_species[max_days_used]
  
  total_species <- n_distinct(detections_data[[species_col]])
  
  cat("\nDIVERSITY METRICS :\n")
  cat(sprintf("Alpha-diversity (day 1): %.1f species\n", alpha_diversity))
  cat(sprintf("Temporal beta (slope): %.2f species/day\n", beta_temporal))
  cat(sprintf("Gamma-diversity (spatial): %.1f species\n", gamma_diversity))
  cat(sprintf("Total species detected: %d\n", total_species))
  
  result <- list(
    single_recorder = single_recorder,
    composite = composite,
    site_curves = site_curves,
    alpha_diversity = alpha_diversity,
    beta_temporal = beta_temporal,
    gamma_diversity = gamma_diversity,
    parameters = list(
      species_col = species_col,
      confidence_threshold = confidence_threshold,
      n_permutations = n_permutations,
      n_sites_composite = n_sites_composite,
      max_days = max_days,
      seed = seed
    ),
    n_sites_used = n_sites_used,
    max_days_used = max_days_used,
    total_species = total_species
  )
  
  class(result) <- c("species_accumulation_abg", "list")
  return(result)
}


# ==============================================================================
# SECTION 2: VISUALIZATION - LAUREL SYMES STYLE
# ==============================================================================

#' Plot Species Accumulation with Alpha-Beta-Gamma Annotations
#'
#' @description
#' Creates a publication-quality species accumulation plot following the style
#' of Laurel Symes' Ecology & Evolution figures. Shows single-recorder curve
#' (solid blue), composite multi-recorder curve (solid black), SE envelopes,
#' and annotated diversity components.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item \strong{Blue solid line}: Average single recorder with dotted SE bounds
#'   \item \strong{Black solid line}: Composite of N recorders with dotted SE bounds
#'   \item \strong{Alpha annotation}: Circle at day 1, labeled
#'   \item \strong{Temporal beta annotation}: Arrow along single-recorder curve
#'   \item \strong{Gamma annotation}: Vertical arrow between curves at right side
#'   \item \strong{Definition box}: Explanatory text in bottom-right corner
#' }
#'
#' Design choices follow Wiley Ecology and Evolution journal standards:
#' clean white background, minimal gridlines, clear legends.
#'
#' @param accum_result Object of class "species_accumulation_abg" from
#'   \code{\link{calculate_accumulation_abg}}
#' @param taxon_name Character. Taxon for title (e.g., "bird", "katydid").
#'   Default: "species"
#' @param show_annotations Logical. Show alpha/beta/gamma labels. Default: TRUE
#' @param show_definitions Logical. Show definition box. Default: TRUE
#' @param colors Named vector with "single" and "composite" colors.
#'   Default: c(single = "#0072B2", composite = "#000000")
#' @param se_line_type Character. Line type for SE bounds. Default: "dotted"
#' @param output_file Character or NULL. Path to save plot. Default: NULL
#' @param width Numeric. Plot width in inches. Default: 10
#' @param height Numeric. Plot height in inches. Default: 7
#' @param dpi Integer. Resolution for saved plot. Default: 300
#'
#' @return A ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' p <- plot_accumulation_abg(katydid_accum, taxon_name = "katydid")
#'
#' # Save to file
#' plot_accumulation_abg(
#'   bird_accum,
#'   taxon_name = "bird",
#'   output_file = "outputs/bird_accumulation_abg.jpg"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon annotate labs
#'   scale_x_continuous scale_y_continuous scale_color_manual theme_minimal
#'   theme element_text element_rect element_line ggsave
#' @importFrom dplyr mutate filter bind_rows
#'
#' @export
plot_accumulation_abg <- function(accum_result,
                                  taxon_name = "species",
                                  show_annotations = TRUE,
                                  show_definitions = TRUE,
                                  colors = c(single = "#0072B2", 
                                             composite = "#000000"),
                                  se_line_type = "dotted",
                                  output_file = NULL,
                                  width = 10,
                                  height = 7,
                                  dpi = 300) {
  
  # Validate input
  if (!inherits(accum_result, "species_accumulation_abg")) {
    warning("Input should be from calculate_accumulation_abg()")
  }
  
  # Extract data
  single_df <- accum_result$single_recorder %>%
    mutate(curve = "single")
  composite_df <- accum_result$composite %>%
    mutate(curve = "composite")
  
  n_sites <- accum_result$n_sites_used
  max_day <- accum_result$max_days_used
  
  # Get y-axis range for positioning
  min_y <- min(c(single_df$lower_ci, composite_df$lower_ci), na.rm = TRUE)
  max_y <- max(c(single_df$upper_ci, composite_df$upper_ci), na.rm = TRUE)
  y_range <- max_y - min_y
  
  # Build base plot
  p <- ggplot() +
    # Composite: main line (draw first so single is on top)
    geom_line(data = composite_df,
              aes(x = day, y = mean_species),
              color = colors["composite"],
              linewidth = 1.2) +
    # Composite: SE bounds (dotted)
    geom_line(data = composite_df,
              aes(x = day, y = lower_ci),
              color = colors["composite"],
              linetype = se_line_type,
              linewidth = 0.6) +
    geom_line(data = composite_df,
              aes(x = day, y = upper_ci),
              color = colors["composite"],
              linetype = se_line_type,
              linewidth = 0.6) +
    # Single recorder: main line
    geom_line(data = single_df, 
              aes(x = day, y = mean_species),
              color = colors["single"], 
              linewidth = 1.2) +
    # Single recorder: SE bounds (dotted)
    geom_line(data = single_df,
              aes(x = day, y = lower_ci),
              color = colors["single"],
              linetype = se_line_type,
              linewidth = 0.6) +
    geom_line(data = single_df,
              aes(x = day, y = upper_ci),
              color = colors["single"],
              linetype = se_line_type,
              linewidth = 0.6) +
    # Axes
    scale_x_continuous(
      breaks = seq(1, max_day, by = 5),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      expand = c(0.05, 0)
    ) +
    # Labels
    labs(
      x = "Number of recorder days analyzed",
      y = sprintf("Total %s species detected", taxon_name)
    ) +
    # Theme: clean, publication-ready
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Add minimal annotations (just the alpha circle marker)
  if (show_annotations) {
    alpha_y <- accum_result$alpha_diversity
    
    # Alpha annotation - just the circle at day 1, no label
    p <- p +
      annotate("point", x = 1, y = alpha_y,
               shape = 1, size = 6, stroke = 1.2,
               color = colors["single"])
  }
  
  # # Add definition box
  # if (show_definitions) {
  #   definitions <- paste0(
  #     "  * alpha-diversity: expected number of species\n",
  #     "      detected from analysis of one recorder on\n",
  #     "      one day.\n",
  #     "  * beta-diversity: rate at which new species are\n",
  #     "      detected as more days are analyzed from the\n",
  #     "      same recorder.\n",
  #     "  * gamma-diversity: increment in species detections\n",
  #     "      from analyzing the same total number of\n",
  #     "      recorder-days but drawn from ", n_sites, " recorders\n",
  #     "      instead of a single recorder."
  #   )
  #   
  #   # Position: bottom right, inside plot area
  #   box_y <- min_y + y_range * 0.25
  #   box_x <- max_day * 0.98
  #   
  #   p <- p +
  #     annotate("label",
  #              x = box_x, y = box_y,
  #              label = definitions,
  #              hjust = 1, vjust = 0.5,
  #              size = 2.8,
  #              fill = alpha("white", 0.9),
  #              color = "grey30",
  #              label.padding = unit(0.5, "lines"),
  #              label.r = unit(0.15, "lines"),
  #              family = "mono")
  # }
  
  # Save if requested
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, dpi = dpi,
           bg = "white")
    cat(sprintf("Plot saved: %s\n", output_file))
  }
  
  return(p)
}


# ==============================================================================
# SECTION 3: DUAL-TAXON ANALYSIS WRAPPER
# ==============================================================================

#' Generate Accumulation Analysis for Birds and Katydids
#'
#' @description
#' Convenience function that calculates and plots alpha-beta-gamma accumulation
#' curves for both bird and katydid datasets.
#'
#' @param bird_detections Data frame with bird detection data
#' @param katydid_detections Data frame with katydid detection data
#' @param bird_confidence Numeric. Confidence threshold for birds. Default: 0.9
#' @param katydid_confidence Numeric. Confidence threshold for katydids. Default: 0.9
#' @param n_sites Integer. Number of sites for composite. Default: 10
#' @param max_days Integer. Maximum days for curves. Default: 16
#' @param output_dir Character. Directory for output files. Default: "outputs"
#' @param save_plots Logical. Save plots to files. Default: TRUE
#'
#' @return List containing:
#'   \describe{
#'     \item{bird}{List with accumulation data and plot for birds}
#'     \item{katydid}{List with accumulation data and plot for katydids}
#'     \item{summary}{Data frame comparing diversity metrics}
#'   }
#'
#' @examples
#' \dontrun{
#' results <- analyze_accumulation_both_taxa(
#'   bird_detections = bird_data$raw_detections,
#'   katydid_detections = katydid_data$raw_detections,
#'   n_sites = 10,
#'   max_days = 16,
#'   output_dir = "outputs/accumulation"
#' )
#' }
#'
#' @export
analyze_accumulation_both_taxa <- function(bird_detections = NULL,
                                           katydid_detections = NULL,
                                           bird_confidence = 0.9,
                                           katydid_confidence = 0.95,
                                           n_sites = 25,
                                           max_days = 16,
                                           output_dir = "outputs",
                                           save_plots = TRUE) {
  
  cat("\n")
  cat("ALPHA-BETA-GAMMA ACCUMULATION ANALYSIS\n")
  cat("\n")
  
  # Create output directory
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list()
  
  # Bird analysis
  if (!is.null(bird_detections) && nrow(bird_detections) > 0) {
    cat("\n--- BIRDS ---\n")
    
    tryCatch({
      results$bird$accumulation <- calculate_accumulation_abg(
        detections_data = bird_detections,
        species_col = "common_name",
        confidence_threshold = bird_confidence,
        n_sites_composite = n_sites,
        max_days = max_days
      )
      
      out_file <- if (save_plots) file.path(output_dir, "bird_accumulation_abg.jpg") else NULL
      
      results$bird$plot <- plot_accumulation_abg(
        results$bird$accumulation,
        taxon_name = "bird",
        output_file = out_file
      )
      
    }, error = function(e) {
      warning(sprintf("Bird analysis failed: %s", e$message))
    })
  }
  
  # Katydid analysis
  if (!is.null(katydid_detections) && nrow(katydid_detections) > 0) {
    cat("\n--- KATYDIDS ---\n")
    
    tryCatch({
      results$katydid$accumulation <- calculate_accumulation_abg(
        detections_data = katydid_detections,
        species_col = "common_name",
        confidence_threshold = katydid_confidence,
        n_sites_composite = n_sites,
        max_days = max_days
      )
      
      out_file <- if (save_plots) file.path(output_dir, "katydid_accumulation_abg.jpg") else NULL
      
      results$katydid$plot <- plot_accumulation_abg(
        results$katydid$accumulation,
        taxon_name = "katydid",
        colors = c(single = "#009E73", composite = "#000000"),
        output_file = out_file
      )
      
    }, error = function(e) {
      warning(sprintf("Katydid analysis failed: %s", e$message))
    })
  }
  
  # Create summary comparison
  summary_data <- data.frame(
    taxon = character(),
    alpha = numeric(),
    beta_temporal = numeric(),
    gamma = numeric(),
    total_species = integer(),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(results$bird$accumulation)) {
    summary_data <- rbind(summary_data, data.frame(
      taxon = "Bird",
      alpha = results$bird$accumulation$alpha_diversity,
      beta_temporal = results$bird$accumulation$beta_temporal,
      gamma = results$bird$accumulation$gamma_diversity,
      total_species = results$bird$accumulation$total_species
    ))
  }
  
  if (!is.null(results$katydid$accumulation)) {
    summary_data <- rbind(summary_data, data.frame(
      taxon = "Katydid",
      alpha = results$katydid$accumulation$alpha_diversity,
      beta_temporal = results$katydid$accumulation$beta_temporal,
      gamma = results$katydid$accumulation$gamma_diversity,
      total_species = results$katydid$accumulation$total_species
    ))
  }
  
  results$summary <- summary_data
  
  # Print summary
  cat("\n")
  cat("SUMMARY: DIVERSITY COMPONENTS\n")
  cat("\n")
  print(summary_data, row.names = FALSE)
  
  return(results)
}


# ==============================================================================
# END OF species_accumulation_abg.R
# ==============================================================================
