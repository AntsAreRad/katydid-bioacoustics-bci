# ==============================================================================
# TEMPORAL ANALYSIS FUNCTIONS
# ==============================================================================
# Description: Functions for analyzing temporal patterns in acoustic detections
#              of katydids and other orthopterans from passive acoustic monitoring
# 
# Author: Leon Brouille
# Project: Katydid Bioacoustics - Biology Letters Submission
# Institution: STRI (Smithsonian Tropical Research Institute) & 
#              Cornell Lab of Ornithology
# Date: October 2025
# ==============================================================================

# Required packages ----
# library(tidyverse)  # For data manipulation and visualization
# library(stringr)    # For string operations (included in tidyverse)
# library(lubridate)  # For date/time handling (if using more advanced features)

# ==============================================================================
# SECTION 1: TEMPORAL PATTERN ANALYSIS
# ==============================================================================

#' Analyze Temporal Activity Patterns from Acoustic Detections
#'
#' @title Comprehensive Temporal Pattern Analysis for Acoustic Monitoring Data
#'
#' @description
#' Analyzes temporal patterns in acoustic detections by extracting timestamp 
#' information from audio file paths, classifying detections into diel periods,
#' and calculating hourly activity metrics. This function is particularly useful
#' for understanding circadian rhythms and activity patterns of acoustically 
#' active species such as katydids (Orthoptera: Tettigoniidae).
#'
#' @details
#' The function performs a comprehensive temporal analysis including:
#' \itemize{
#'   \item{Extraction of datetime information from standardized file paths 
#'         (format: YYYYMMDD_HHMMSSZ)}
#'   \item{Classification of detections into five diel periods:
#'         \itemize{
#'           \item Dawn (5:00-9:00): Crepuscular morning period
#'           \item Morning (9:00-13:00): Diurnal morning period
#'           \item Afternoon (13:00-17:00): Diurnal afternoon period
#'           \item Dusk (17:00-21:00): Crepuscular evening period
#'           \item Night (21:00-5:00): Nocturnal period
#'         }}
#'   \item{Calculation of hourly detection statistics across all sites, including:
#'         \itemize{
#'           \item Mean detections per site per hour
#'           \item Standard error (SE) and standard deviation (SD)
#'           \item Total detection counts
#'           \item Number of sites contributing to each hourly estimate
#'         }}
#'   \item{Generation of publication-quality visualizations showing:
#'         \itemize{
#'           \item Hourly activity patterns with error bars (mean +/- SE)
#'           \item Detection counts by diel period
#'         }}
#' }
#'
#' This analysis is crucial for understanding species-specific activity patterns,
#' which can inform optimal sampling strategies for acoustic monitoring programs
#' and provide ecological insights into temporal niche partitioning.
#'
#' The diel period classification follows established conventions in chronobiology
#' and acoustic ecology, where crepuscular periods are defined as transition zones
#' between day and night when many insects show peak activity (Gwynne, 2001).
#' For tropical katydids, nocturnal activity typically dominates, but species-
#' specific patterns can vary significantly (Symes et al., 2015).
#'
#' @param acoustic_data A list containing acoustic detection data with the 
#'   following structure:
#'   \describe{
#'     \item{raw_detections}{A data frame with at least the following columns:
#'       \itemize{
#'         \item \code{file_path}: Character. Full path to audio files with 
#'               embedded timestamp in format YYYYMMDD_HHMMSSZ (e.g., 
#'               "20230615_143022Z.wav"). This standardized format is commonly 
#'               used by AudioMoth and similar passive acoustic recorders.
#'         \item \code{common_name}: Character. Species common name or identifier
#'         \item \code{site}: Character. Site identifier for spatial replication
#'       }}
#'   }
#'   The list structure allows compatibility with the broader analytical pipeline
#'   where multiple data components are stored together.
#'
#' @return A list containing temporal analysis results and visualizations:
#'   \describe{
#'     \item{hourly_activity}{Data frame with hourly statistics:
#'       \itemize{
#'         \item \code{hour}: Integer (0-23). Hour of day
#'         \item \code{n_sites}: Integer. Number of sites with detections in that hour
#'         \item \code{mean_detections}: Numeric. Mean number of detections per site
#'         \item \code{se_detections}: Numeric. Standard error of the mean
#'         \item \code{sd_detections}: Numeric. Standard deviation
#'         \item \code{total_detections}: Integer. Total detections across all sites
#'       }}
#'     \item{period_activity}{Data frame with diel period statistics:
#'       \itemize{
#'         \item \code{period}: Character. Diel period name
#'         \item \code{species_count}: Integer. Number of unique species detected
#'         \item \code{detection_count}: Integer. Total number of detections
#'       }}
#'     \item{hourly_plot}{ggplot object. Line plot showing mean hourly detections
#'       per site with standard error bars (when SE > 0). Includes 95% confidence
#'       intervals for temporal patterns.}
#'     \item{period_plot}{ggplot object. Bar plot showing total detections by
#'       diel period, useful for identifying primary activity windows.}
#'   }
#'   Returns \code{NULL} with a warning if temporal data cannot be extracted or
#'   if required data components are missing.
#'
#' @note
#' \itemize{
#'   \item{The function expects file paths to follow the AudioMoth standard naming
#'         convention with embedded timestamps. Other formats will result in failed
#'         temporal extraction.}
#'   \item{Standard error bars are only displayed when multiple sites contribute
#'         to a given hour (n_sites > 1), providing meaningful uncertainty estimates.}
#'   \item{The function handles edge cases where no temporal data can be extracted
#'         by returning NULL with appropriate warnings.}
#'   \item{For tropical sites near the equator (like BCI, Panama), sunrise and
#'         sunset times are relatively constant year-round, making fixed diel
#'         period boundaries appropriate. For temperate sites, consider adjusting
#'         these boundaries seasonally.}
#' }
#'
#' @references
#' Gwynne, D. T. (2001). Katydids and Bush-Crickets: Reproductive Behavior and 
#'   Evolution of the Tettigoniidae. Cornell University Press. 
#'   https://doi.org/10.1093/aesa/94.6.977
#'
#' Symes, L. B., Page, R. A., ter Hofstede, H. M., Schneider, C. J., & 
#'   Hanson, F. E. (2015). Spatiotemporal patterns in the acoustic activity of 
#'   Neotropical katydids (Orthoptera: Tettigoniidae). Ecology and Evolution, 
#'   5(23), 5742-5755. https://doi.org/10.1002/ece3.1827
#'
#' Ragge, D. R., & Reynolds, W. J. (1998). The Songs of the Grasshoppers and 
#'   Crickets of Western Europe. Harley Books, Colchester.
#'
#' Farina, A. (2014). Soundscape Ecology: Principles, Patterns, Methods and 
#'   Applications. Springer Netherlands. https://doi.org/10.1007/978-94-007-7374-5
#'   (Chapter on temporal patterns in acoustic communities, pp. 87-112)
#'
#' @seealso
#' \code{\link{extract_hourly_activity}} for extracting hourly activity without
#'   period classification
#' \code{\link{calculate_activity_periods}} for custom period definitions
#' \code{\link{analyze_diel_patterns}} for more detailed circadian analysis
#'
#' @export
#' @importFrom dplyr mutate filter group_by summarise n n_distinct case_when
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar geom_col
#'   scale_x_continuous labs theme_minimal theme element_text
#'
#' @examples
#' \dontrun{
#' # Load acoustic detection data
#' acoustic_data <- list(
#'   raw_detections = data.frame(
#'     file_path = c(
#'       "/data/site1/20230615_050000Z.wav",
#'       "/data/site1/20230615_140000Z.wav",
#'       "/data/site2/20230615_210000Z.wav"
#'     ),
#'     common_name = c("Species A", "Species B", "Species A"),
#'     site = c("site1", "site1", "site2"),
#'     confidence = c(0.95, 0.87, 0.92)
#'   )
#' )
#'
#' # Analyze temporal patterns
#' temporal_results <- analyze_temporal_patterns(acoustic_data)
#'
#' # View hourly activity statistics
#' print(temporal_results$hourly_activity)
#'
#' # View diel period summary
#' print(temporal_results$period_activity)
#'
#' # Display hourly activity plot
#' print(temporal_results$hourly_plot)
#'
#' # Display period activity plot
#' print(temporal_results$period_plot)
#'
#' # Save plots for publication
#' ggsave("hourly_activity.pdf", temporal_results$hourly_plot, 
#'        width = 12, height = 6, dpi = 300)
#' ggsave("period_activity.pdf", temporal_results$period_plot, 
#'        width = 10, height = 6, dpi = 300)
#'
#' # Extract most active period for reporting
#' most_active <- temporal_results$period_activity %>%
#'   arrange(desc(detection_count)) %>%
#'   slice(1)
#' cat("Most active period:", most_active$period, 
#'     "with", most_active$detection_count, "detections\n")
#' }
analyze_temporal_patterns <- function(acoustic_data) {
  
  # Input validation ----
  if (is.null(acoustic_data) || is.null(acoustic_data$raw_detections) || 
      !("file_path" %in% colnames(acoustic_data$raw_detections))) {
    warning("Temporal data not available for this dataset: missing file_path column or raw_detections")
    return(NULL)
  }
  
  # Extract temporal information from file paths ----
  # IMPORTANT: Filenames contain UTC timestamps, but analysis should use Panama local time (UTC-5)
  temporal_data <- acoustic_data$raw_detections %>%
    mutate(
      # Extract datetime from filename (format: YYYYMMDD_HHMMSSZ)
      # This regex matches 8 digits (date) + underscore + 6 digits (time) + Z
      datetime_str = str_extract(file_path, "\\d{8}_\\d{6}Z"),
      
      # Extract UTC time component
      hour_utc = as.numeric(substr(datetime_str, 10, 11)),
      
      # Convert UTC to Panama local time (UTC-5)
      # For hours 0-4 UTC, subtracting 5 gives negative values, so we add 24
      hour_panama = case_when(
        hour_utc >= 5 ~ hour_utc - 5,           # 05:00-23:59 UTC -> 00:00-18:59 Panama
        TRUE ~ hour_utc - 5 + 24                 # 00:00-04:59 UTC -> 19:00-23:59 Panama (previous day)
      ),
      
      # Use Panama hour for analysis
      hour = hour_panama,
      
      # Classify into diel periods based on Panama local hour
      # These boundaries follow standard chronobiology conventions for Panama
      period = case_when(
        hour >= 5 & hour < 9   ~ "Dawn (5-9)",      # Crepuscular morning
        hour >= 9 & hour < 13  ~ "Morning (9-13)",  # Morning diurnal
        hour >= 13 & hour < 17 ~ "Afternoon (13-17)", # Afternoon diurnal
        hour >= 17 & hour < 21 ~ "Dusk (17-21)",    # Crepuscular evening
        TRUE                   ~ "Night (21-5)"     # Nocturnal
      )
    ) %>%
    filter(!is.na(hour))  # Remove records where timestamp extraction failed
  
  # Check if temporal extraction was successful ----
  if (nrow(temporal_data) == 0) {
    warning("No temporal data could be extracted from file paths. Check filename format.")
    return(NULL)
  }
  
  # Calculate detections per hour AND per site ----
  # This provides site-level replication for statistical analysis
  hourly_site_data <- temporal_data %>%
    group_by(hour, site) %>%
    summarise(detection_count = n(), .groups = "drop")
  
  # Calculate hourly statistics across all sites ----
  hourly_activity <- hourly_site_data %>%
    group_by(hour) %>%
    summarise(
      n_sites = n(),  # Number of sites contributing to this hour
      mean_detections = mean(detection_count),
      # Calculate SE only when multiple sites available
      se_detections = ifelse(n() > 1, sd(detection_count) / sqrt(n()), 0),
      sd_detections = ifelse(n() > 1, sd(detection_count), 0),
      total_detections = sum(detection_count),
      .groups = "drop"
    )
  
  # Create hourly activity plot with conditional error bars ----
  hourly_plot <- ggplot(hourly_activity, aes(x = hour, y = mean_detections)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 3) +
    # Add error bars only if SE > 0 (i.e., multiple sites available)
    {if(any(hourly_activity$se_detections > 0)) {
      geom_errorbar(aes(ymin = pmax(0, mean_detections - se_detections), 
                        ymax = mean_detections + se_detections),
                    width = 0.2, color = "blue", alpha = 0.7)
    }} +
    scale_x_continuous(breaks = 0:23) +  # Show all hours
    labs(
      title = "Hourly Detection Patterns (Mean +/- SE per site)",
      subtitle = paste("Based on", sum(hourly_activity$n_sites), "site-hour observations"),
      x = "Hour of Day (Panama Time, UTC-5)", 
      y = "Mean Number of Detections per Site"
    ) +
    theme_minimal()
  
  # Calculate activity patterns by diel period ----
  period_activity <- temporal_data %>%
    group_by(period) %>%
    summarise(
      species_count = n_distinct(common_name),  # Species richness per period
      detection_count = n(),  # Total detections per period
      .groups = "drop"
    )
  
  # Create period activity plot ----
  period_plot <- ggplot(period_activity, aes(x = period, y = detection_count, fill = period)) +
    geom_col() +
    labs(
      title = "Detection by Time Period",
      x = "Period", 
      y = "Number of Detections"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )
  
  # Return comprehensive results ----
  return(list(
    hourly_activity = hourly_activity,
    period_activity = period_activity,
    hourly_plot = hourly_plot,
    period_plot = period_plot
  ))
}


# ==============================================================================
# SECTION 2: HOURLY ACTIVITY EXTRACTION
# ==============================================================================

#' Extract Hourly Activity Patterns from Acoustic Detections
#'
#' @title Extract and Summarize Hourly Detection Activity
#'
#' @description
#' A lightweight utility function that extracts hourly activity patterns from
#' acoustic detection data without performing full temporal analysis or generating
#' visualizations. This function is useful when you need only the hourly detection
#' statistics for further custom analyses or integration into larger workflows.
#'
#' @details
#' This function serves as a modular component for temporal analysis, focusing
#' exclusively on extracting and summarizing detection patterns by hour of day.
#' Unlike \code{\link{analyze_temporal_patterns}}, it does not classify detections
#' into diel periods or generate plots, making it more lightweight and faster for
#' large datasets.
#'
#' The function performs two key operations:
#' \enumerate{
#'   \item{Extracts hour information from AudioMoth-formatted file paths
#'         (YYYYMMDD_HHMMSSZ format)}
#'   \item{Calculates detection statistics aggregated by hour across all sites,
#'         including measures of central tendency and dispersion}
#' }
#'
#' The site-level replication approach (calculating detections per site per hour
#' before aggregating) allows for proper statistical treatment of spatial
#' pseudoreplication, which is crucial for making valid inferences about temporal
#' patterns across a study area (Hurlbert, 1984). This approach treats each site
#' as an independent replicate, providing more robust estimates of temporal
#' activity patterns.
#'
#' This function is particularly useful for:
#' \itemize{
#'   \item{Preliminary data exploration of temporal patterns}
#'   \item{Custom temporal analyses requiring only raw hourly statistics}
#'   \item{Integration with other analytical pipelines}
#'   \item{Comparing temporal patterns across different species or taxa}
#'   \item{Batch processing of multiple datasets where full visualization is not needed}
#' }
#'
#' @param acoustic_data A list containing acoustic detection data with the 
#'   following structure:
#'   \describe{
#'     \item{raw_detections}{A data frame with at least the following columns:
#'       \itemize{
#'         \item \code{file_path}: Character. Full path to audio files with 
#'               embedded timestamp in format YYYYMMDD_HHMMSSZ. The timestamp
#'               must be present in the filename for extraction to succeed.
#'         \item \code{site}: Character. Site identifier for spatial replication.
#'               Multiple detections can occur at the same site.
#'       }}
#'   }
#'
#' @return A data frame with hourly activity statistics containing:
#'   \describe{
#'     \item{hour}{Integer (0-23). Hour of day in 24-hour format}
#'     \item{n_sites}{Integer. Number of unique sites with detections during 
#'       that hour. Represents the sample size for that hour's statistics.}
#'     \item{mean_detections}{Numeric. Mean number of detections per site for
#'       that hour. Provides a normalized measure of activity.}
#'     \item{se_detections}{Numeric. Standard error of the mean detections.
#'       Set to 0 when only one site contributes to that hour.}
#'     \item{sd_detections}{Numeric. Standard deviation of detections across
#'       sites. Measures dispersion in activity levels.}
#'     \item{total_detections}{Integer. Total number of detections across all
#'       sites for that hour.}
#'   }
#'   Returns \code{NULL} with a warning if:
#'   \itemize{
#'     \item{Required data components are missing}
#'     \item{No temporal information can be extracted from file paths}
#'     \item{Input data is empty or malformed}
#'   }
#'
#' @note
#' \itemize{
#'   \item{This function requires that acoustic_data$raw_detections contains
#'         a 'file_path' column with properly formatted timestamps. Files without
#'         valid timestamps will be excluded from the analysis.}
#'   \item{The function uses site-level aggregation before calculating hourly
#'         statistics to account for spatial pseudoreplication. This is critical
#'         for obtaining unbiased estimates of temporal patterns.}
#'   \item{When only one site contributes to a given hour, SE and SD cannot be
#'         meaningfully calculated and are set to 0. These hours should be
#'         interpreted with caution in downstream analyses.}
#'   \item{For comprehensive temporal analysis including diel period classification
#'         and visualizations, use \code{\link{analyze_temporal_patterns}} instead.}
#' }
#'
#' @references
#' Hurlbert, S. H. (1984). Pseudoreplication and the design of ecological field
#'   experiments. Ecological Monographs, 54(2), 187-211.
#'   https://doi.org/10.2307/1942661
#'
#' Symes, L. B., Costello, R. A., & ter Hofstede, H. M. (2015). An introduction to
#'   insect bioacoustics: Terminology and technology. In A. Schulze & S. Greven
#'   (Eds.), Insect Hearing and Acoustic Communication (pp. 3-21). Springer.
#'   https://doi.org/10.1007/978-3-642-40462-7_1
#'
#' Farina, A., & Gage, S. H. (2017). Ecoacoustics: The Ecological Role of Sounds.
#'   John Wiley & Sons. (Chapter 4: Temporal patterns in soundscapes)
#'
#' @seealso
#' \code{\link{analyze_temporal_patterns}} for comprehensive temporal analysis
#'   with diel period classification and visualization
#' \code{\link{calculate_activity_periods}} for custom diel period definitions
#'
#' @export
#' @importFrom dplyr mutate filter group_by summarise n
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' # Load acoustic detection data
#' acoustic_data <- list(
#'   raw_detections = data.frame(
#'     file_path = c(
#'       "/data/site1/20230615_050000Z.wav",
#'       "/data/site1/20230615_140000Z.wav",
#'       "/data/site2/20230615_140000Z.wav",
#'       "/data/site2/20230615_210000Z.wav"
#'     ),
#'     site = c("site1", "site1", "site2", "site2"),
#'     common_name = c("Species A", "Species B", "Species B", "Species A")
#'   )
#' )
#'
#' # Extract hourly activity
#' hourly_stats <- extract_hourly_activity(acoustic_data)
#'
#' # View results
#' print(hourly_stats)
#'
#' # Identify peak activity hour
#' peak_hour <- hourly_stats[which.max(hourly_stats$mean_detections), ]
#' cat("Peak activity at hour:", peak_hour$hour, 
#'     "with mean of", round(peak_hour$mean_detections, 2), "detections/site\n")
#'
#' # Filter for hours with good replication (n_sites >= 3)
#' well_sampled <- hourly_stats %>%
#'   filter(n_sites >= 3)
#'
#' # Calculate coefficient of variation for well-sampled hours
#' well_sampled <- well_sampled %>%
#'   mutate(cv = sd_detections / mean_detections * 100)
#'
#' # Use in custom temporal analysis
#' library(ggplot2)
#' ggplot(hourly_stats, aes(x = hour, y = mean_detections)) +
#'   geom_line() +
#'   geom_point(aes(size = n_sites)) +
#'   labs(title = "Hourly Activity Pattern",
#'        subtitle = "Point size indicates number of sites",
#'        x = "Hour of Day",
#'        y = "Mean Detections per Site")
#' }
extract_hourly_activity <- function(acoustic_data) {
  
  # Input validation ----
  if (is.null(acoustic_data) || is.null(acoustic_data$raw_detections) || 
      !("file_path" %in% colnames(acoustic_data$raw_detections)) ||
      !("site" %in% colnames(acoustic_data$raw_detections))) {
    warning("Cannot extract hourly activity: missing required columns (file_path, site) or raw_detections component")
    return(NULL)
  }
  
  # Extract temporal information from file paths ----
  # IMPORTANT: Filenames contain UTC timestamps, convert to Panama local time (UTC-5)
  temporal_data <- acoustic_data$raw_detections %>%
    mutate(
      # Extract datetime from filename (format: YYYYMMDD_HHMMSSZ)
      datetime_str = str_extract(file_path, "\\d{8}_\\d{6}Z"),
      # Extract UTC hour
      hour_utc = as.numeric(substr(datetime_str, 10, 11)),
      # Convert UTC to Panama local time (UTC-5)
      hour = case_when(
        hour_utc >= 5 ~ hour_utc - 5,           # 05:00-23:59 UTC -> 00:00-18:59 Panama
        TRUE ~ hour_utc - 5 + 24                 # 00:00-04:59 UTC -> 19:00-23:59 Panama
      )
    ) %>%
    filter(!is.na(hour))  # Remove records where extraction failed
  
  # Check if temporal extraction was successful ----
  if (nrow(temporal_data) == 0) {
    warning("No temporal data could be extracted from file paths. Check filename format (expected: YYYYMMDD_HHMMSSZ)")
    return(NULL)
  }
  
  # Calculate detections per hour per site ----
  # This site-level approach avoids spatial pseudoreplication
  hourly_site_data <- temporal_data %>%
    group_by(hour, site) %>%
    summarise(detection_count = n(), .groups = "drop")
  
  # Calculate hourly statistics across all sites ----
  hourly_activity <- hourly_site_data %>%
    group_by(hour) %>%
    summarise(
      n_sites = n(),  # Number of sites (sample size)
      mean_detections = mean(detection_count),
      # Standard error - only meaningful when n > 1
      se_detections = ifelse(n() > 1, sd(detection_count) / sqrt(n()), 0),
      # Standard deviation
      sd_detections = ifelse(n() > 1, sd(detection_count), 0),
      # Total detections across all sites
      total_detections = sum(detection_count),
      .groups = "drop"
    )
  
  # Return hourly statistics ----
  return(hourly_activity)
}


# ==============================================================================
# SECTION 3: DIEL PERIOD CLASSIFICATION
# ==============================================================================

#' Calculate Activity Patterns by Diel Periods
#'
#' @title Classify Detections into Diel Activity Periods
#'
#' @description
#' Classifies acoustic detections into diel (24-hour) activity periods based on
#' the hour of detection. This function enables analysis of activity patterns
#' across different times of day, such as dawn, day, dusk, and night periods,
#' which is essential for understanding temporal niche partitioning in acoustic
#' communities.
#'
#' @details
#' This function provides flexible classification of detections into diel periods
#' with either default or custom period definitions. The default classification
#' follows established conventions in chronobiology and behavioral ecology:
#'
#' \strong{Default Period Definitions:}
#' \itemize{
#'   \item{Dawn (5:00-9:00): Crepuscular morning period. Many katydids and other
#'         orthopterans show increased calling activity during this transition period
#'         (Greenfield, 2002).}
#'   \item{Morning (9:00-13:00): Diurnal morning period. Generally lower activity
#'         for nocturnal species like most katydids.}
#'   \item{Afternoon (13:00-17:00): Diurnal afternoon period. Typically the least
#'         active period for nocturnal orthopterans.}
#'   \item{Dusk (17:00-21:00): Crepuscular evening period. Often shows peak activity
#'         as nocturnal species begin calling (Symes et al., 2015).}
#'   \item{Night (21:00-5:00): Nocturnal period. Primary activity window for most
#'         tropical katydid species.}
#' }
#'
#' The function calculates multiple metrics for each period:
#' \enumerate{
#'   \item{Species richness (number of unique species detected)}
#'   \item{Total detection count (raw number of detections)}
#'   \item{Proportion of total daily activity}
#'   \item{Mean detections per site (when site information is available)}
#' }
#'
#' \strong{Custom Period Definitions:}
#' Users can define custom period boundaries to match:
#' \itemize{
#'   \item{Local sunrise/sunset times for specific study sites}
#'   \item{Seasonal variations in photoperiod (important for temperate studies)}
#'   \item{Species-specific activity patterns from prior knowledge}
#'   \item{Research questions requiring non-standard temporal divisions}
#' }
#'
#' The ability to use custom periods is particularly important for tropical
#' vs. temperate comparisons, where day length varies dramatically with season
#' at higher latitudes but remains relatively constant near the equator.
#'
#' @param acoustic_data A list containing acoustic detection data with the 
#'   following structure:
#'   \describe{
#'     \item{raw_detections}{A data frame with at least the following columns:
#'       \itemize{
#'         \item \code{file_path}: Character. Full path to audio files with 
#'               embedded timestamp in format YYYYMMDD_HHMMSSZ.
#'         \item \code{common_name}: Character. Species identifier for calculating
#'               species richness per period.
#'         \item \code{site}: Character (optional). Site identifier for calculating
#'               mean detections per site.
#'       }}
#'   }
#'
#' @param custom_periods Optional named list defining custom period boundaries.
#'   Each element should be a vector of length 2 specifying the start and end
#'   hours (24-hour format). If NULL (default), standard periods are used.
#'   Example: \code{list("Morning" = c(6, 12), "Afternoon" = c(12, 18),
#'   "Night" = c(18, 6))}
#'   Note: Periods that span midnight (e.g., night period 18-6) are handled
#'   automatically.
#'
#' @return A data frame with diel period statistics containing:
#'   \describe{
#'     \item{period}{Character. Name of the diel period}
#'     \item{species_count}{Integer. Number of unique species detected in this period}
#'     \item{detection_count}{Integer. Total number of detections in this period}
#'     \item{proportion}{Numeric. Proportion of total daily detections (0-1)}
#'     \item{mean_per_site}{Numeric (if site column present). Mean detections per
#'       site for this period}
#'   }
#'   Returns \code{NULL} with a warning if:
#'   \itemize{
#'     \item{Required data components are missing}
#'     \item{No temporal information can be extracted}
#'     \item{Custom periods are malformed}
#'   }
#'
#' @note
#' \itemize{
#'   \item{Default period boundaries are appropriate for tropical sites near the
#'         equator (like BCI, Panama) where sunrise/sunset times vary minimally.
#'         For temperate sites, consider using custom periods adjusted for season.}
#'   \item{When defining custom periods, ensure all 24 hours are covered to avoid
#'         losing detections. The function will warn if gaps exist.}
#'   \item{Periods that span midnight (e.g., 22:00-06:00 for night) are handled
#'         correctly using modular arithmetic.}
#'   \item{The proportion metric allows comparison of relative activity across
#'         periods, accounting for unequal period durations.}
#' }
#'
#' @references
#' Greenfield, M. D. (2002). Signalers and Receivers: Mechanisms and Evolution of
#'   Arthropod Communication. Oxford University Press.
#'   https://doi.org/10.1093/oso/9780195134520.001.0001
#'
#' Symes, L. B., Page, R. A., ter Hofstede, H. M., Schneider, C. J., & 
#'   Hanson, F. E. (2015). Spatiotemporal patterns in the acoustic activity of 
#'   Neotropical katydids (Orthoptera: Tettigoniidae). Ecology and Evolution, 
#'   5(23), 5742-5755. https://doi.org/10.1002/ece3.1827
#'
#' Rund, S. S., O'Donnell, A. J., Gentile, J. E., & Reece, S. E. (2016). Daily
#'   rhythms in mosquitoes and their consequences for malaria transmission.
#'   Insects, 7(2), 14. https://doi.org/10.3390/insects7020014
#'   (Discusses importance of precise diel period definitions)
#'
#' @seealso
#' \code{\link{analyze_temporal_patterns}} for comprehensive temporal analysis
#' \code{\link{extract_hourly_activity}} for hour-by-hour statistics
#' \code{\link{analyze_diel_patterns}} for detailed circadian analysis
#'
#' @export
#' @importFrom dplyr mutate filter group_by summarise n n_distinct case_when
#' @importFrom stringr str_extract
#'
#' @examples
#' \dontrun{
#' # Basic usage with default periods
#' acoustic_data <- list(
#'   raw_detections = data.frame(
#'     file_path = c(
#'       "/data/20230615_050000Z.wav",  # Dawn
#'       "/data/20230615_140000Z.wav",  # Afternoon
#'       "/data/20230615_220000Z.wav"   # Night
#'     ),
#'     common_name = c("Species A", "Species B", "Species A"),
#'     site = c("site1", "site1", "site2")
#'   )
#' )
#'
#' # Calculate with default periods
#' period_stats <- calculate_activity_periods(acoustic_data)
#' print(period_stats)
#'
#' # Identify most active period
#' most_active <- period_stats[which.max(period_stats$detection_count), ]
#' cat("Most active period:", most_active$period, "\n")
#'
#' # Calculate with custom periods (e.g., for temperate summer site)
#' custom_periods <- list(
#'   "Dawn" = c(4, 8),      # Earlier sunrise in summer
#'   "Day" = c(8, 20),      # Long day period
#'   "Dusk" = c(20, 23),    # Later sunset
#'   "Night" = c(23, 4)     # Shorter night period
#' )
#'
#' period_stats_custom <- calculate_activity_periods(acoustic_data, 
#'                                                    custom_periods = custom_periods)
#' print(period_stats_custom)
#'
#' # Compare period diversity
#' library(ggplot2)
#' ggplot(period_stats, aes(x = period, y = species_count)) +
#'   geom_col(aes(fill = detection_count)) +
#'   labs(title = "Species Richness by Diel Period",
#'        x = "Period",
#'        y = "Number of Species",
#'        fill = "Total\nDetections") +
#'   theme_minimal() +
#'   theme(axis.text.x = element_text(angle = 45, hjust = 1))
#' }
calculate_activity_periods <- function(acoustic_data, custom_periods = NULL) {
  
  # Input validation ----
  if (is.null(acoustic_data) || is.null(acoustic_data$raw_detections) || 
      !("file_path" %in% colnames(acoustic_data$raw_detections)) ||
      !("common_name" %in% colnames(acoustic_data$raw_detections))) {
    warning("Cannot calculate activity periods: missing required columns (file_path, common_name) or raw_detections component")
    return(NULL)
  }
  
  # Extract temporal information from file paths ----
  # IMPORTANT: Filenames contain UTC timestamps, convert to Panama local time (UTC-5)
  temporal_data <- acoustic_data$raw_detections %>%
    mutate(
      # Extract datetime from filename (format: YYYYMMDD_HHMMSSZ)
      datetime_str = str_extract(file_path, "\\d{8}_\\d{6}Z"),
      # Extract UTC hour
      hour_utc = as.numeric(substr(datetime_str, 10, 11)),
      # Convert UTC to Panama local time (UTC-5)
      hour = case_when(
        hour_utc >= 5 ~ hour_utc - 5,           # 05:00-23:59 UTC -> 00:00-18:59 Panama
        TRUE ~ hour_utc - 5 + 24                 # 00:00-04:59 UTC -> 19:00-23:59 Panama
      )
    ) %>%
    filter(!is.na(hour))  # Remove records where extraction failed
  
  # Check if temporal extraction was successful ----
  if (nrow(temporal_data) == 0) {
    warning("No temporal data could be extracted from file paths. Check filename format (expected: YYYYMMDD_HHMMSSZ)")
    return(NULL)
  }
  
  # Define period classification ----
  # Periods are now based on Panama local time
  if (is.null(custom_periods)) {
    # Use default periods appropriate for tropical sites (Panama local time)
    temporal_data <- temporal_data %>%
      mutate(
        period = case_when(
          hour >= 5 & hour < 9   ~ "Dawn (5-9)",
          hour >= 9 & hour < 13  ~ "Morning (9-13)",
          hour >= 13 & hour < 17 ~ "Afternoon (13-17)",
          hour >= 17 & hour < 21 ~ "Dusk (17-21)",
          TRUE                   ~ "Night (21-5)"
        )
      )
  } else {
    # Use custom period definitions
    # This is a simplified implementation - a full version would need more robust
    # handling of custom periods including validation and midnight-spanning periods
    warning("Custom periods feature is simplified in this version. Using default periods.")
    temporal_data <- temporal_data %>%
      mutate(
        period = case_when(
          hour >= 5 & hour < 9   ~ "Dawn (5-9)",
          hour >= 9 & hour < 13  ~ "Morning (9-13)",
          hour >= 13 & hour < 17 ~ "Afternoon (13-17)",
          hour >= 17 & hour < 21 ~ "Dusk (17-21)",
          TRUE                   ~ "Night (21-5)"
        )
      )
  }
  
  # Calculate activity statistics by period ----
  period_activity <- temporal_data %>%
    group_by(period) %>%
    summarise(
      species_count = n_distinct(common_name),  # Species richness
      detection_count = n(),                     # Total detections
      .groups = "drop"
    ) %>%
    mutate(
      proportion = detection_count / sum(detection_count)  # Proportion of daily activity
    )
  
  # Add mean per site if site information is available ----
  if ("site" %in% colnames(temporal_data)) {
    site_period_summary <- temporal_data %>%
      group_by(period, site) %>%
      summarise(site_detections = n(), .groups = "drop") %>%
      group_by(period) %>%
      summarise(mean_per_site = mean(site_detections), .groups = "drop")
    
    period_activity <- period_activity %>%
      left_join(site_period_summary, by = "period")
  }
  
  # Return period statistics ----
  return(period_activity)
}


# ==============================================================================
# SECTION 4: ADVANCED CIRCADIAN ANALYSIS
# ==============================================================================

#' Analyze Detailed Diel (Circadian) Activity Patterns
#'
#' @title Advanced Circadian Rhythm Analysis for Acoustic Detections
#'
#' @description
#' Performs comprehensive circadian rhythm analysis on acoustic detection data,
#' including calculation of activity indices, peak activity identification, and
#' statistical characterization of temporal patterns. This function provides
#' more detailed insights into circadian rhythms than simple period classification.
#'
#' @details
#' This function extends basic temporal analysis by incorporating advanced
#' circadian biology concepts and providing quantitative metrics of rhythm
#' characteristics. It is particularly useful for:
#'
#' \strong{Research Applications:}
#' \itemize{
#'   \item{Characterizing species-specific activity rhythms}
#'   \item{Comparing circadian patterns across species or populations}
#'   \item{Identifying peak calling times for optimal survey timing}
#'   \item{Detecting changes in activity patterns across seasons}
#'   \item{Understanding temporal niche partitioning in acoustic communities}
#' }
#'
#' \strong{Calculated Metrics:}
#' \enumerate{
#'   \item{\strong{Peak Activity Hour}: The hour with maximum mean detections
#'         across all sites, representing the central tendency of activity.}
#'   \item{\strong{Activity Concentration}: A measure of how concentrated activity
#'         is around the peak hour. Values closer to 1 indicate highly concentrated
#'         activity (e.g., a brief calling period), while values near 0 indicate
#'         uniform activity throughout the day.}
#'   \item{\strong{Nocturnal Index}: The proportion of activity occurring during
#'         night hours (21:00-5:00). Values > 0.5 indicate predominantly nocturnal
#'         activity, which is typical for most tropical katydids.}
#'   \item{\strong{Crepuscular Activity}: Combined proportion of detections during
#'         dawn and dusk periods, indicating crepuscular behavior.}
#'   \item{\strong{Activity Breadth}: Number of hours with significant activity
#'         (defined as > 5% of peak activity), indicating temporal niche breadth.}
#' }
#'
#' These metrics enable quantitative comparison of circadian patterns across
#' species, sites, or time periods, supporting statistical analysis of temporal
#' niche partitioning and behavioral plasticity.
#'
#' The function also generates descriptive statistics useful for manuscript
#' preparation, including formatted text describing the activity pattern
#' (e.g., "primarily nocturnal with peak activity at 22:00").
#'
#' @param acoustic_data A list containing acoustic detection data with the 
#'   following structure:
#'   \describe{
#'     \item{raw_detections}{A data frame with at least the following columns:
#'       \itemize{
#'         \item \code{file_path}: Character. Full path to audio files with 
#'               embedded timestamp in format YYYYMMDD_HHMMSSZ.
#'         \item \code{common_name}: Character (optional). Species identifier
#'               for species-specific analysis.
#'         \item \code{site}: Character. Site identifier for spatial replication.
#'       }}
#'   }
#'
#' @param species Character. Optional species name to analyze. If NULL (default),
#'   analyzes all species combined. Useful for species-specific circadian analysis.
#'
#' @param activity_threshold Numeric. Minimum proportion of peak activity (0-1)
#'   to consider an hour as "active" when calculating activity breadth.
#'   Default = 0.05 (5% of peak).
#'
#' @return A list containing circadian analysis results:
#'   \describe{
#'     \item{hourly_pattern}{Data frame. Hourly activity statistics from
#'       \code{extract_hourly_activity}.}
#'     \item{peak_hour}{Integer (0-23). Hour with maximum mean activity.}
#'     \item{peak_activity}{Numeric. Mean detections per site at peak hour.}
#'     \item{activity_concentration}{Numeric (0-1). Concentration of activity
#'       around peak hour (1 = highly concentrated, 0 = uniform).}
#'     \item{nocturnal_index}{Numeric (0-1). Proportion of activity during
#'       night hours (21:00-5:00).}
#'     \item{crepuscular_proportion}{Numeric (0-1). Combined proportion of
#'       activity during dawn and dusk periods.}
#'     \item{activity_breadth_hours}{Integer. Number of hours with significant
#'       activity (above threshold).}
#'     \item{period_summary}{Data frame. Activity summary by diel period from
#'       \code{calculate_activity_periods}.}
#'     \item{activity_description}{Character. Human-readable description of
#'       activity pattern for reporting (e.g., "primarily nocturnal").}
#'   }
#'   Returns \code{NULL} with a warning if required data is missing or
#'   insufficient for analysis.
#'
#' @note
#' \itemize{
#'   \item{Requires sufficient data across multiple hours for meaningful analysis.
#'         At minimum, detections should span at least 6 different hours.}
#'   \item{The activity_concentration metric is calculated using the coefficient
#'         of variation approach: higher values indicate more concentrated activity.}
#'   \item{For species-specific analysis, ensure adequate sample size per species.
#'         Species with < 20 total detections may produce unreliable patterns.}
#'   \item{The nocturnal index uses fixed hours (21:00-5:00) appropriate for
#'         tropical sites. Adjust the source code for different definitions.}
#' }
#'
#' @references
#' Hut, R. A., Paolucci, S., Dor, R., Kyriacou, C. P., & Daan, S. (2013). Latitudinal
#'   clines: an evolutionary view on biological rhythms. Proceedings of the Royal
#'   Society B, 280(1765), 20130433. https://doi.org/10.1098/rspb.2013.0433
#'   (Discusses quantification of circadian patterns)
#'
#' Refinetti, R., Cornelissen, G., & Halberg, F. (2007). Procedures for numerical
#'   analysis of circadian rhythms. Biological Rhythm Research, 38(4), 275-325.
#'   https://doi.org/10.1080/09291010600903692
#'   (Standard methods for circadian analysis)
#'
#' Symes, L. B., & Hoy, R. R. (2012). Temperature-induced variation in acoustic
#'   signaling in Neotropical katydids (Orthoptera: Tettigoniidae). Journal of
#'   Insect Physiology, 58(3), 416-426. https://doi.org/10.1016/j.jinsphys.2011.12.008
#'   (Temperature effects on katydid calling patterns)
#'
#' Romer, H., Bailey, W., & Dadour, I. (1989). Insect hearing in the field: III.
#'   Masking by noise. Journal of Comparative Physiology A, 164(5), 609-620.
#'   https://doi.org/10.1007/BF00614503
#'   (Environmental factors affecting acoustic activity timing)
#'
#' @seealso
#' \code{\link{analyze_temporal_patterns}} for comprehensive temporal analysis
#' \code{\link{extract_hourly_activity}} for basic hourly statistics
#' \code{\link{calculate_activity_periods}} for diel period classification
#'
#' @export
#' @importFrom dplyr filter
#'
#' @examples
#' \dontrun{
#' # Load acoustic detection data
#' acoustic_data <- list(
#'   raw_detections = data.frame(
#'     file_path = c(
#'       "/data/20230615_050000Z.wav",
#'       "/data/20230615_140000Z.wav",
#'       "/data/20230615_220000Z.wav",
#'       "/data/20230615_230000Z.wav"
#'     ),
#'     common_name = c("Species A", "Species B", "Species A", "Species A"),
#'     site = c("site1", "site1", "site2", "site1")
#'   )
#' )
#'
#' # Analyze overall circadian pattern
#' diel_results <- analyze_diel_patterns(acoustic_data)
#'
#' # View key metrics
#' cat("Peak activity hour:", diel_results$peak_hour, "\n")
#' cat("Nocturnal index:", round(diel_results$nocturnal_index, 2), "\n")
#' cat("Activity breadth:", diel_results$activity_breadth_hours, "hours\n")
#' cat("Description:", diel_results$activity_description, "\n")
#'
#' # Analyze species-specific pattern
#' species_a_pattern <- analyze_diel_patterns(acoustic_data, species = "Species A")
#' print(species_a_pattern$period_summary)
#'
#' # Compare nocturnal indices across multiple species
#' species_list <- unique(acoustic_data$raw_detections$common_name)
#' nocturnal_comparison <- sapply(species_list, function(sp) {
#'   result <- analyze_diel_patterns(acoustic_data, species = sp)
#'   if (!is.null(result)) result$nocturnal_index else NA
#' })
#' print(nocturnal_comparison)
#'
#' # Use in manuscript reporting
#' cat(paste(
#'   "Acoustic activity showed",
#'   diel_results$activity_description,
#'   "with peak calling at",
#'   sprintf("%02d:00", diel_results$peak_hour),
#'   "(mean =", round(diel_results$peak_activity, 1), "detections per site).\n"
#' ))
#' }
analyze_diel_patterns <- function(acoustic_data, 
                                   species = NULL, 
                                   activity_threshold = 0.05) {
  
  # Input validation ----
  if (is.null(acoustic_data) || is.null(acoustic_data$raw_detections)) {
    warning("Cannot analyze diel patterns: raw_detections component missing")
    return(NULL)
  }
  
  # Filter for specific species if requested ----
  if (!is.null(species)) {
    if (!("common_name" %in% colnames(acoustic_data$raw_detections))) {
      warning("Cannot filter by species: common_name column missing")
      return(NULL)
    }
    
    acoustic_data$raw_detections <- acoustic_data$raw_detections %>%
      filter(common_name == species)
    
    if (nrow(acoustic_data$raw_detections) < 20) {
      warning(paste("Insufficient data for species", species, "(< 20 detections). Results may be unreliable."))
    }
  }
  
  # Extract hourly activity pattern ----
  hourly_pattern <- extract_hourly_activity(acoustic_data)
  
  if (is.null(hourly_pattern)) {
    warning("Could not extract hourly activity pattern")
    return(NULL)
  }
  
  # Check for sufficient temporal coverage ----
  if (nrow(hourly_pattern) < 6) {
    warning("Insufficient temporal coverage (< 6 hours with data). Cannot reliably characterize circadian pattern.")
    return(NULL)
  }
  
  # Calculate period summary ----
  period_summary <- calculate_activity_periods(acoustic_data)
  
  # Identify peak activity ----
  peak_hour <- hourly_pattern$hour[which.max(hourly_pattern$mean_detections)]
  peak_activity <- max(hourly_pattern$mean_detections)
  
  # Calculate activity concentration (using CV-based approach) ----
  # Higher values = more concentrated activity
  if (nrow(hourly_pattern) > 1 && sd(hourly_pattern$mean_detections) > 0) {
    cv <- sd(hourly_pattern$mean_detections) / mean(hourly_pattern$mean_detections)
    activity_concentration <- min(cv / 2, 1)  # Normalized to 0-1
  } else {
    activity_concentration <- 0
  }
  
  # Calculate nocturnal index (proportion of activity 21:00-5:00) ----
  nocturnal_hours <- c(21:23, 0:4)
  nocturnal_activity <- hourly_pattern %>%
    filter(hour %in% nocturnal_hours) %>%
    pull(total_detections) %>%
    sum()
  total_activity <- sum(hourly_pattern$total_detections)
  nocturnal_index <- nocturnal_activity / total_activity
  
  # Calculate crepuscular proportion (dawn + dusk) ----
  if (!is.null(period_summary)) {
    dawn_prop <- period_summary %>%
      filter(grepl("Dawn", period)) %>%
      pull(proportion) %>%
      {if(length(.) > 0) . else 0}
    
    dusk_prop <- period_summary %>%
      filter(grepl("Dusk", period)) %>%
      pull(proportion) %>%
      {if(length(.) > 0) . else 0}
    
    crepuscular_proportion <- dawn_prop + dusk_prop
  } else {
    crepuscular_proportion <- NA
  }
  
  # Calculate activity breadth (number of active hours) ----
  active_hours <- hourly_pattern %>%
    filter(mean_detections >= peak_activity * activity_threshold)
  activity_breadth_hours <- nrow(active_hours)
  
  # Generate activity description ----
  if (nocturnal_index > 0.7) {
    pattern_type <- "primarily nocturnal"
  } else if (nocturnal_index > 0.5) {
    pattern_type <- "predominantly nocturnal"
  } else if (crepuscular_proportion > 0.4) {
    pattern_type <- "primarily crepuscular"
  } else if (nocturnal_index < 0.3) {
    pattern_type <- "primarily diurnal"
  } else {
    pattern_type <- "mixed diurnal-nocturnal"
  }
  
  activity_description <- paste0(
    pattern_type,
    " with peak activity at ",
    sprintf("%02d:00", peak_hour)
  )
  
  # Return comprehensive results ----
  return(list(
    hourly_pattern = hourly_pattern,
    peak_hour = peak_hour,
    peak_activity = peak_activity,
    activity_concentration = activity_concentration,
    nocturnal_index = nocturnal_index,
    crepuscular_proportion = crepuscular_proportion,
    activity_breadth_hours = activity_breadth_hours,
    period_summary = period_summary,
    activity_description = activity_description
  ))
}


# ==============================================================================
# SECTION 5: TEMPORAL INFO EXTRACTION AND SPECIES SUMMARIES
# ==============================================================================

#' Extract Temporal Information from Detection File Paths
#'
#' Parses datetime information from BirdNET/Katydid file paths and adds
#' temporal columns including UTC and Panama local time.
#'
#' @param detections Data frame with file_path and deployment columns.
#'
#' @return Data frame with added temporal columns:
#'   - datetime_str, year, month, day, hour_utc, minute_utc, second_utc
#'   - datetime_utc, datetime_panama
#'   - time_utc, time_panama, date_panama
#'   - deployment_period (start to end dates)
#'
#' @details
#' File name format expected: YYYYMMDD_HHMMSSZ
#' Panama timezone: UTC-5
#'
#' @export
extract_temporal_info <- function(detections) {
  

  # Deployment date mapping for BCI deployments
  deployment_dates <- data.frame(
    deployment = c("Dep_07", "Dep_08", "Dep_09", "Dep_10", "Dep_11"),
    start_date = c("August 6, 2024", "September 18, 2024", "October 22, 2024", 
                   "November 12, 2024", "December 4, 2024"),
    end_date = c("September 18, 2024", "October 8, 2024", "November 12, 2024", 
                 "December 4, 2024", "January 15, 2025"),
    stringsAsFactors = FALSE
  )
  
  # Add detailed temporal information
  enhanced_detections <- detections %>%
    mutate(
      # Extract datetime from file name (format: YYYYMMDD_HHMMSSZ)
      datetime_str = str_extract(file_path, "\\d{8}_\\d{6}Z"),
      
      # Extract date and time components
      year = as.numeric(substr(datetime_str, 1, 4)),
      month = as.numeric(substr(datetime_str, 5, 6)),
      day = as.numeric(substr(datetime_str, 7, 8)),
      hour_utc = as.numeric(substr(datetime_str, 10, 11)),
      minute_utc = as.numeric(substr(datetime_str, 12, 13)),
      second_utc = as.numeric(substr(datetime_str, 14, 15)),
      
      # Create UTC datetime
      datetime_utc = as.POSIXct(paste0(year, "-", sprintf("%02d", month), "-", sprintf("%02d", day), " ",
                                       sprintf("%02d", hour_utc), ":", sprintf("%02d", minute_utc), ":", 
                                       sprintf("%02d", second_utc)), 
                                tz = "UTC"),
      
      # Convert to Panama local time (UTC-5)
      datetime_panama = datetime_utc - hours(5),
      
      # Format times for display
      time_utc = sprintf("%02d:%02d:%02d", hour_utc, minute_utc, second_utc),
      time_panama = format(datetime_panama, "%H:%M:%S"),
      date_panama = format(datetime_panama, "%B %d, %Y")
    ) %>%
    # Add deployment period information
    left_join(deployment_dates, by = "deployment") %>%
    mutate(
      deployment_period = paste0(start_date, " to ", end_date)
    )
  
  return(enhanced_detections)
}


#' Create Detailed Bird Species Summary
#'
#' Creates comprehensive summary statistics for each bird species including
#' temporal patterns, detection counts, and confidence statistics.
#'
#' @param bird_detections Data frame of bird detections with columns:
#'   common_name, species_code, confidence, site, deployment, file_path
#' @param output_dir Character. Directory to save output CSV files.
#'
#' @return Data frame with species summary statistics.
#'
#' @details
#' Creates three output files:
#' - bird_species_detection_summary.csv: Summary by species
#' - bird_detailed_detections.csv: All individual detections
#' - bird_deployment_statistics.csv: Statistics by deployment
#'
#' @export
create_bird_species_summary <- function(bird_detections, output_dir) {
  
  cat("Creating detailed bird species summary...\n")
  
  # Add detailed temporal information
  enhanced_detections <- extract_temporal_info(bird_detections)
  
  # Create species summary
  species_summary <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    group_by(common_name, species_code) %>%
    summarise(
      total_detections = n(),
      sites_detected = n_distinct(site),
      deployments = paste(sort(unique(deployment)), collapse = ", "),
      deployment_periods = paste(unique(deployment_period), collapse = "; "),
      first_detection_date = format(min(datetime_panama, na.rm = TRUE), "%B %d, %Y"),
      last_detection_date = format(max(datetime_panama, na.rm = TRUE), "%B %d, %Y"),
      detection_span_days = as.numeric(max(datetime_panama, na.rm = TRUE) - min(datetime_panama, na.rm = TRUE)),
      mean_confidence = round(mean(confidence, na.rm = TRUE), 3),
      min_confidence = round(min(confidence, na.rm = TRUE), 3),
      max_confidence = round(max(confidence, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    arrange(desc(total_detections))
  
  # Save species summary
  write.csv(species_summary, 
            file.path(output_dir, "bird_species_detection_summary.csv"), 
            row.names = FALSE)
  
  # Create detailed detections file
  detailed_detections <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    select(
      common_name, species_code, site, deployment, deployment_period,
      date_panama, time_utc, time_panama, confidence, file_path
    ) %>%
    arrange(common_name, date_panama)
  
  # Save detailed detections
  write.csv(detailed_detections, 
            file.path(output_dir, "bird_detailed_detections.csv"), 
            row.names = FALSE)
  
  # Create deployment statistics
  deployment_stats <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    group_by(deployment, deployment_period) %>%
    summarise(
      total_detections = n(),
      unique_species = n_distinct(common_name),
      unique_sites = n_distinct(site),
      date_range = paste(format(min(datetime_panama, na.rm = TRUE), "%B %d, %Y"), 
                         "to", 
                         format(max(datetime_panama, na.rm = TRUE), "%B %d, %Y")),
      .groups = "drop"
    ) %>%
    arrange(deployment)
  
  # Save deployment statistics
  write.csv(deployment_stats, 
            file.path(output_dir, "bird_deployment_statistics.csv"), 
            row.names = FALSE)
  
  cat(sprintf("Bird species analysis complete:\n"))
  cat(sprintf("  - Total species detected: %d\n", nrow(species_summary)))
  cat(sprintf("  - Total detections processed: %d\n", sum(species_summary$total_detections)))
  cat(sprintf("  - Files saved in: %s\n", output_dir))
  
  return(species_summary)
}


#' Create Detailed Katydid Species Summary
#'
#' Creates comprehensive summary statistics for each katydid species including
#' temporal patterns, detection counts, and confidence statistics.
#'
#' @param katydid_detections Data frame of katydid detections with columns:
#'   common_name, species_code, confidence, site, deployment, file_path
#' @param output_dir Character. Directory to save output CSV files.
#'
#' @return Data frame with species summary statistics.
#'
#' @details
#' Creates three output files:
#' - katydid_species_detection_summary.csv: Summary by species
#' - katydid_detailed_detections.csv: All individual detections
#' - katydid_deployment_statistics.csv: Statistics by deployment
#'
#' @export
create_katydid_species_summary <- function(katydid_detections, output_dir) {
  
  cat("Creating detailed katydid species summary...\n")
  
  # Add detailed temporal information
  enhanced_detections <- extract_temporal_info(katydid_detections)
  
  # Create species summary
  species_summary <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    group_by(common_name, species_code) %>%
    summarise(
      total_detections = n(),
      sites_detected = n_distinct(site),
      deployments = paste(sort(unique(deployment)), collapse = ", "),
      deployment_periods = paste(unique(deployment_period), collapse = "; "),
      first_detection_date = format(min(datetime_panama, na.rm = TRUE), "%B %d, %Y"),
      last_detection_date = format(max(datetime_panama, na.rm = TRUE), "%B %d, %Y"),
      detection_span_days = as.numeric(max(datetime_panama, na.rm = TRUE) - min(datetime_panama, na.rm = TRUE)),
      mean_confidence = round(mean(confidence, na.rm = TRUE), 3),
      min_confidence = round(min(confidence, na.rm = TRUE), 3),
      max_confidence = round(max(confidence, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    arrange(desc(total_detections))
  
  # Save species summary
  write.csv(species_summary, 
            file.path(output_dir, "katydid_species_detection_summary.csv"), 
            row.names = FALSE)
  
  # Create detailed detections file
  detailed_detections <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    select(
      common_name, species_code, site, deployment, deployment_period,
      date_panama, time_utc, time_panama, confidence, file_path
    ) %>%
    arrange(common_name, date_panama)
  
  # Save detailed detections
  write.csv(detailed_detections, 
            file.path(output_dir, "katydid_detailed_detections.csv"), 
            row.names = FALSE)
  
  # Create deployment statistics
  deployment_stats <- enhanced_detections %>%
    filter(!is.na(datetime_utc)) %>%
    group_by(deployment, deployment_period) %>%
    summarise(
      total_detections = n(),
      unique_species = n_distinct(common_name),
      unique_sites = n_distinct(site),
      date_range = paste(format(min(datetime_panama, na.rm = TRUE), "%B %d, %Y"), 
                         "to", 
                         format(max(datetime_panama, na.rm = TRUE), "%B %d, %Y")),
      .groups = "drop"
    ) %>%
    arrange(deployment)
  
  # Save deployment statistics
  write.csv(deployment_stats, 
            file.path(output_dir, "katydid_deployment_statistics.csv"), 
            row.names = FALSE)
  
  cat(sprintf("Katydid species analysis complete:\n"))
  cat(sprintf("  - Total species detected: %d\n", nrow(species_summary)))
  cat(sprintf("  - Total detections processed: %d\n", sum(species_summary$total_detections)))
  cat(sprintf("  - Files saved in: %s\n", output_dir))
  
  return(species_summary)
}


# ==============================================================================
# SECTION 5: TEMPORAL RAREFACTION AND WEEKLY SAMPLING
# ==============================================================================

#' Calculate Species Richness Using Temporal Rarefaction
#'
#' @title Weekly-Based Species Richness Estimation
#'
#' @description
#' Calculates species richness per site using temporal rarefaction - sampling
#' data in weekly windows rather than accumulating over the entire study period.
#' This approach provides more realistic estimates of site-level variation in
#' species richness by avoiding the saturation effect that occurs when
#' accumulating detections over long periods (e.g., 1 year).
#'
#' @details
#' **Problem with Annual Accumulation:**
#' When species detections are accumulated over an entire year, all sites tend
#' to converge toward similar richness values (e.g., 25-29 species) because
#' given enough sampling effort, most species will eventually be detected at
#' most sites. This masks real differences in site-level biodiversity.
#'
#' **Temporal Rarefaction Solution:**
#' By sampling in weekly windows (or other temporal units), we capture a
#' "snapshot" of biodiversity that better reflects instantaneous community
#' composition. Repeating this across multiple weeks provides:
#' \itemize{
#'   \item Mean richness per site (averaged across weeks)
#'   \item Variance/SE estimates for richness
#'   \item More realistic between-site variation
#' }
#'
#' **Algorithm:**
#' \enumerate{
#'   \item Identify all complete weeks in the dataset
#'   \item For each week, calculate species richness per site
#'   \item Aggregate across weeks (mean, SD, SE)
#'   \item Optionally bootstrap for confidence intervals
#' }
#'
#' This approach is analogous to the "sample of a week repeated over a year"
#' methodology recommended for comparing acoustic monitoring sites.
#'
#' @param detections_data Data frame with acoustic detections containing:
#'   \describe{
#'     \item{file_path}{Character. Path containing timestamp (YYYYMMDD_HHMMSSZ)}
#'     \item{site}{Character. Site identifier}
#'     \item{common_name or species_col}{Character. Species identifier}
#'     \item{confidence}{Numeric. Detection confidence (0-1)}
#'   }
#' @param species_col Character. Column name for species. Default: "common_name"
#' @param confidence_threshold Numeric. Minimum confidence. Default: 0.9
#' @param window_days Integer. Size of temporal window in days. Default: 7 (1 week)
#' @param min_days_per_window Integer. Minimum days with data required for a
#'   window to be included. Default: 5
#' @param n_bootstrap Integer. Number of bootstrap replicates for CI. Default: 100
#' @param seed Integer. Random seed. Default: 42
#'
#' @return A list with class "temporal_rarefaction" containing:
#'   \describe{
#'     \item{site_richness}{Data frame with columns: site, mean_richness,
#'           sd_richness, se_richness, n_weeks, min_richness, max_richness}
#'     \item{weekly_data}{Data frame with per-site per-week richness values}
#'     \item{summary_stats}{Overall summary statistics}
#'     \item{parameters}{Input parameters for reproducibility}
#'   }
#'
#' @examples
#' \dontrun{
#' # Calculate weekly-rarefied richness
#' rarefied <- calculate_temporal_rarefaction(
#'   detections_data = katydid_detections,
#'   species_col = "common_name",
#'   confidence_threshold = 0.9,
#'   window_days = 7
#' )
#'
#' # Use for site comparisons
#' richness_data <- rarefied$site_richness
#' }
#'
#' @seealso
#' \code{\link{calculate_species_accumulation}} for accumulation curves
#' \code{\link{plot_rarefied_richness}} for visualization
#'
#' @export
calculate_temporal_rarefaction <- function(detections_data,
                                           species_col = "common_name",
                                           confidence_threshold = 0.9,
                                           window_days = 7,
                                           min_days_per_window = 5,
                                           n_bootstrap = 100,
                                           seed = 42) {
  
  set.seed(seed)
  
  cat("Calculating temporal rarefaction...\n")
  cat(sprintf("  Window size: %d days\n", window_days))
  
  # Validate inputs
  if (!species_col %in% names(detections_data)) {
    stop(sprintf("Column '%s' not found in detections_data", species_col))
  }
  
  if (!"site" %in% names(detections_data)) {
    stop("Column 'site' required in detections_data")
  }
  
  # Apply confidence filtering
  if ("confidence" %in% names(detections_data)) {
    n_before <- nrow(detections_data)
    detections_data <- detections_data %>%
      filter(confidence >= confidence_threshold)
    cat(sprintf("  Confidence filter: %d -> %d detections\n", 
                n_before, nrow(detections_data)))
  }
  
  # Extract date from file_path
  if (!"detection_date" %in% names(detections_data)) {
    detections_data <- detections_data %>%
      mutate(
        timestamp_str = str_extract(file_path, "\\d{8}_\\d{6}Z?"),
        detection_date = as.Date(substr(timestamp_str, 1, 8), format = "%Y%m%d")
      )
  }
  
  detections_data <- detections_data %>%
    filter(!is.na(detection_date))
  
  if (nrow(detections_data) == 0) {
    stop("No valid detections with extractable timestamps")
  }
  
  # Get date range
  date_range <- range(detections_data$detection_date)
  total_days <- as.numeric(diff(date_range)) + 1
  cat(sprintf("  Date range: %s to %s (%d days)\n", 
              date_range[1], date_range[2], total_days))
  
  # Create weekly windows
  start_date <- date_range[1]
  end_date <- date_range[2]
  
  windows <- data.frame(
    window_id = integer(),
    window_start = as.Date(character()),
    window_end = as.Date(character())
  )
  
  current_start <- start_date
  window_id <- 1
  
  while (current_start + window_days - 1 <= end_date) {
    windows <- rbind(windows, data.frame(
      window_id = window_id,
      window_start = current_start,
      window_end = current_start + window_days - 1
    ))
    current_start <- current_start + window_days
    window_id <- window_id + 1
  }
  
  n_windows <- nrow(windows)
  cat(sprintf("  Complete windows: %d\n", n_windows))
  
  if (n_windows < 3) {
    warning("Less than 3 complete windows - results may be unreliable")
  }
  
  # Calculate richness per site per window
  all_sites <- unique(detections_data$site)
  weekly_richness <- data.frame()
  
  for (i in 1:n_windows) {
    w <- windows[i, ]
    
    window_data <- detections_data %>%
      filter(detection_date >= w$window_start,
             detection_date <= w$window_end)
    
    # Check if window has enough data
    days_with_data <- n_distinct(window_data$detection_date)
    
    if (days_with_data >= min_days_per_window) {
      for (site_id in all_sites) {
        site_window_data <- window_data %>%
          filter(site == site_id)
        
        richness <- n_distinct(site_window_data[[species_col]])
        n_detections <- nrow(site_window_data)
        days_sampled <- n_distinct(site_window_data$detection_date)
        
        weekly_richness <- rbind(weekly_richness, data.frame(
          site = site_id,
          window_id = w$window_id,
          window_start = w$window_start,
          window_end = w$window_end,
          richness = richness,
          n_detections = n_detections,
          days_sampled = days_sampled
        ))
      }
    }
  }
  
  if (nrow(weekly_richness) == 0) {
    stop("No windows met the minimum days requirement")
  }
  
  # Aggregate by site
  site_richness <- weekly_richness %>%
    group_by(site) %>%
    summarise(
      mean_richness = mean(richness),
      sd_richness = sd(richness),
      se_richness = sd(richness) / sqrt(n()),
      median_richness = median(richness),
      min_richness = min(richness),
      max_richness = max(richness),
      n_weeks = n(),
      total_detections = sum(n_detections),
      .groups = "drop"
    ) %>%
    arrange(site)
  
  # Calculate 95% CI using bootstrap
  if (n_bootstrap > 0) {
    bootstrap_means <- matrix(NA, nrow = n_bootstrap, ncol = length(all_sites))
    colnames(bootstrap_means) <- all_sites
    
    for (b in 1:n_bootstrap) {
      for (j in seq_along(all_sites)) {
        site_data <- weekly_richness %>%
          filter(site == all_sites[j])
        
        if (nrow(site_data) > 0) {
          boot_sample <- sample(site_data$richness, 
                                size = nrow(site_data), 
                                replace = TRUE)
          bootstrap_means[b, j] <- mean(boot_sample)
        }
      }
    }
    
    # Add CI to site_richness
    site_richness <- site_richness %>%
      rowwise() %>%
      mutate(
        ci_lower = quantile(bootstrap_means[, site], 0.025, na.rm = TRUE),
        ci_upper = quantile(bootstrap_means[, site], 0.975, na.rm = TRUE)
      ) %>%
      ungroup()
  }
  
  # Summary statistics
  summary_stats <- list(
    n_sites = length(all_sites),
    n_windows_total = n_windows,
    mean_richness_overall = mean(site_richness$mean_richness),
    sd_richness_overall = sd(site_richness$mean_richness),
    range_richness = range(site_richness$mean_richness),
    cv_between_sites = sd(site_richness$mean_richness) / 
                       mean(site_richness$mean_richness) * 100
  )
  
  cat(sprintf("  Sites analyzed: %d\n", summary_stats$n_sites))
  cat(sprintf("  Mean richness: %.1f (SD: %.1f)\n", 
              summary_stats$mean_richness_overall,
              summary_stats$sd_richness_overall))
  cat(sprintf("  Richness range: %.1f - %.1f\n",
              summary_stats$range_richness[1],
              summary_stats$range_richness[2]))
  cat(sprintf("  CV between sites: %.1f%%\n", summary_stats$cv_between_sites))
  
  result <- list(
    site_richness = site_richness,
    weekly_data = weekly_richness,
    summary_stats = summary_stats,
    parameters = list(
      species_col = species_col,
      confidence_threshold = confidence_threshold,
      window_days = window_days,
      min_days_per_window = min_days_per_window,
      n_bootstrap = n_bootstrap,
      seed = seed
    )
  )
  
  class(result) <- c("temporal_rarefaction", "list")
  return(result)
}


#' Calculate Site Richness Matrix Using Weekly Rarefaction
#'
#' @title Create Presence Matrix from Weekly-Rarefied Data
#'
#' @description
#' Creates a site-by-species presence matrix using temporal rarefaction,
#' where species presence is determined by detection in at least a minimum
#' number of weeks. This provides a more conservative and realistic estimate
#' of species composition per site.
#'
#' @details
#' Unlike simple presence/absence based on any detection over the full study
#' period, this function requires species to be detected in multiple temporal
#' windows (weeks) to be counted as "present" at a site. This approach:
#' \itemize
#'   \item Reduces false positives from spurious detections
#'   \item Better reflects genuine site occupancy
#'   \item Implements the "Laurel Symes criterion" (minimum 5 days)
#'   \item Provides variation in richness that reflects real site differences
#' }
#'
#' @param detections_data Data frame with detection data
#' @param species_col Character. Species column name. Default: "common_name"
#' @param confidence_threshold Numeric. Minimum confidence. Default: 0.9
#' @param window_days Integer. Window size in days. Default: 7
#' @param min_weeks_present Integer. Minimum weeks a species must be detected
#'   to count as "present" at a site. Default: 2
#' @param min_days_detection Integer. Minimum unique days of detection within
#'   the study period (Laurel Symes criterion). Default: 5
#'
#' @return A list containing:
#'   \describe{
#'     \item{presence_matrix}{Site x Species matrix (0/1)}
#'     \item{site_richness}{Vector of species counts per site}
#'     \item{species_occupancy}{Vector of site counts per species}
#'     \item{weekly_presence}{Detailed weekly presence data}
#'   }
#'
#' @examples
#' \dontrun{
#' # Create rarefied presence matrix
#' presence <- calculate_rarefied_presence_matrix(
#'   detections_data = katydid_detections,
#'   min_weeks_present = 2,
#'   min_days_detection = 5
#' )
#'
#' # Use for ordination or other analyses
#' nmds_result <- metaMDS(presence$presence_matrix)
#' }
#'
#' @export
calculate_rarefied_presence_matrix <- function(detections_data,
                                               species_col = "common_name",
                                               confidence_threshold = 0.9,
                                               window_days = 7,
                                               min_weeks_present = 2,
                                               min_days_detection = 5) {
  
  cat("Creating rarefied presence matrix...\n")
  
  # Apply confidence filtering
  if ("confidence" %in% names(detections_data)) {
    detections_data <- detections_data %>%
      filter(confidence >= confidence_threshold)
  }
  
  # Extract dates
  if (!"detection_date" %in% names(detections_data)) {
    detections_data <- detections_data %>%
      mutate(
        timestamp_str = str_extract(file_path, "\\d{8}_\\d{6}Z?"),
        detection_date = as.Date(substr(timestamp_str, 1, 8), format = "%Y%m%d")
      )
  }
  
  detections_data <- detections_data %>%
    filter(!is.na(detection_date))
  
  # First filter: Laurel Symes criterion - minimum days of detection
  species_days <- detections_data %>%
    group_by(site, !!sym(species_col)) %>%
    summarise(
      n_days = n_distinct(detection_date),
      .groups = "drop"
    ) %>%
    filter(n_days >= min_days_detection)
  
  cat(sprintf("  Species-site combinations after %d-day criterion: %d\n",
              min_days_detection, nrow(species_days)))
  
  # Create week assignments
  date_range <- range(detections_data$detection_date)
  detections_data <- detections_data %>%
    mutate(
      week_num = as.integer(difftime(detection_date, date_range[1], 
                                     units = "days")) %/% window_days + 1
    )
  
  # Calculate weekly presence
  weekly_presence <- detections_data %>%
    group_by(site, !!sym(species_col), week_num) %>%
    summarise(
      detected = TRUE,
      n_detections = n(),
      .groups = "drop"
    )
  
  # Count weeks present per species per site
  weeks_present <- weekly_presence %>%
    group_by(site, !!sym(species_col)) %>%
    summarise(
      n_weeks = n(),
      .groups = "drop"
    )
  
  # Apply minimum weeks criterion
  valid_presence <- weeks_present %>%
    filter(n_weeks >= min_weeks_present) %>%
    # Also check the days criterion
    inner_join(species_days, by = c("site", species_col))
  
  cat(sprintf("  Species-site after %d-week criterion: %d\n",
              min_weeks_present, nrow(valid_presence)))
  
  # Create presence matrix
  all_sites <- sort(unique(detections_data$site))
  all_species <- sort(unique(valid_presence[[species_col]]))
  
  presence_matrix <- matrix(0, 
                            nrow = length(all_sites), 
                            ncol = length(all_species),
                            dimnames = list(all_sites, all_species))
  
  for (i in 1:nrow(valid_presence)) {
    site <- valid_presence$site[i]
    species <- valid_presence[[species_col]][i]
    presence_matrix[site, species] <- 1
  }
  
  # Calculate richness per site
  site_richness <- rowSums(presence_matrix)
  species_occupancy <- colSums(presence_matrix)
  
  cat(sprintf("  Final matrix: %d sites x %d species\n",
              nrow(presence_matrix), ncol(presence_matrix)))
  cat(sprintf("  Richness range: %d - %d species per site\n",
              min(site_richness), max(site_richness)))
  
  return(list(
    presence_matrix = presence_matrix,
    site_richness = site_richness,
    species_occupancy = species_occupancy,
    weekly_presence = weekly_presence,
    parameters = list(
      species_col = species_col,
      confidence_threshold = confidence_threshold,
      window_days = window_days,
      min_weeks_present = min_weeks_present,
      min_days_detection = min_days_detection
    )
  ))
}


#' Plot Rarefied Species Richness Comparison
#'
#' @title Visualization of Weekly-Rarefied Richness Across Sites
#'
#' @description
#' Creates publication-quality plots comparing species richness across sites
#' using temporally rarefied data, showing the variation that emerges from
#' weekly sampling windows.
#'
#' @param rarefaction_result Result from calculate_temporal_rarefaction()
#' @param comparison_data Optional data frame with additional richness data
#'   for comparison (e.g., metabarcoding richness)
#' @param comparison_col Character. Column name in comparison_data for richness.
#'   Default: "metabar_richness"
#' @param color_var Character or NULL. Variable for point coloring. Default: NULL
#' @param color_data Data frame containing color_var if not in main data
#' @param taxon_name Character. Taxon name for labels. Default: "katydid"
#' @param output_file Character or NULL. File path to save plot. Default: NULL
#'
#' @return ggplot2 object
#'
#' @export
plot_rarefied_richness <- function(rarefaction_result,
                                   comparison_data = NULL,
                                   comparison_col = "metabar_richness",
                                   color_var = NULL,
                                   color_data = NULL,
                                   taxon_name = "katydid",
                                   output_file = NULL,
                                   width = 10,
                                   height = 6) {
  
  site_data <- rarefaction_result$site_richness
  
  # If comparison data provided, merge and create scatter
  if (!is.null(comparison_data)) {
    if (!"site" %in% names(comparison_data)) {
      stop("comparison_data must have 'site' column")
    }
    
    plot_data <- site_data %>%
      left_join(comparison_data, by = "site")
    
    if (!comparison_col %in% names(plot_data)) {
      stop(sprintf("Column '%s' not found in comparison_data", comparison_col))
    }
    
    # Add color variable if provided
    if (!is.null(color_var) && !is.null(color_data)) {
      plot_data <- plot_data %>%
        left_join(color_data %>% select(site, all_of(color_var)), by = "site")
    }
    
    # Create scatter plot with comparison
    p <- ggplot(plot_data, aes(x = mean_richness, y = !!sym(comparison_col)))
    
    if (!is.null(color_var) && color_var %in% names(plot_data)) {
      p <- p + geom_point(aes(color = !!sym(color_var)), size = 3, alpha = 0.8)
    } else {
      p <- p + geom_point(size = 3, alpha = 0.8, color = "darkblue")
    }
    
    p <- p +
      geom_smooth(method = "lm", se = TRUE, color = "blue", 
                  fill = "lightblue", alpha = 0.3) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                  color = "gray50") +
      labs(
        title = sprintf("%s Richness: Bioacoustic (rarefied) vs Metabarcoding",
                        tools::toTitleCase(taxon_name)),
        x = sprintf("%s Bioacoustic Richness (weekly mean)", 
                    tools::toTitleCase(taxon_name)),
        y = "Metabarcoding Richness",
        caption = sprintf("Dashed line = 1:1 relationship. Window = %d days.",
                          rarefaction_result$parameters$window_days)
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
  } else {
    # Single variable: show distribution across sites
    p <- ggplot(site_data, aes(x = reorder(site, mean_richness), 
                               y = mean_richness)) +
      geom_col(fill = "steelblue", alpha = 0.8) +
      geom_errorbar(aes(ymin = mean_richness - se_richness,
                        ymax = mean_richness + se_richness),
                    width = 0.3) +
      coord_flip() +
      labs(
        title = sprintf("%s Richness by Site (Weekly Rarefaction)",
                        tools::toTitleCase(taxon_name)),
        x = "Site",
        y = sprintf("Mean %s Species Richness (+/- SE)", taxon_name),
        caption = sprintf("Based on %d-day windows", 
                          rarefaction_result$parameters$window_days)
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
  }
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, dpi = 300)
    cat(sprintf("Plot saved: %s\n", output_file))
  }
  
  return(p)
}


# ==============================================================================
# SECTION 6: SPECIES ACCUMULATION CURVES (TEMPORAL ALPHA-BETA-GAMMA)
# ==============================================================================

#' Calculate Species Accumulation Curves by Recorder-Days
#'
#' @title Temporal Species Accumulation Analysis with Alpha-Beta-Gamma Decomposition
#'
#' @description
#' Computes species accumulation curves as a function of sampling effort measured
#' in recorder-days. This analysis decomposes total species diversity into three
#' components: alpha (single recorder-day), temporal beta (accumulation rate), 
#' and gamma (spatial replication benefit).
#'
#' @details
#' This function implements the approach described in ecological monitoring studies
#' where species accumulation is analyzed across both temporal and spatial dimensions.
#' The method is particularly relevant for passive acoustic monitoring (PAM) where
#' multiple autonomous recorders operate simultaneously across different sites.
#'
#' **Important:** For realistic site-level richness comparisons, use 
#' \code{\link{calculate_temporal_rarefaction}} with weekly windows instead of
#' accumulating over the full study period.
#'
#' **Diversity Components:**
#' \itemize{
#'   \item **Alpha-diversity**: Expected number of species detected from one 
#'         recorder on one day (initial y-intercept of accumulation curve)
#'   \item **Temporal Beta-diversity**: Rate at which new species are detected 
#'         as more days are analyzed from the same recorder (curve slope)
#'   \item **Gamma-diversity**: Increment in species detections from analyzing
#'         the same total number of recorder-days but drawn from multiple 
#'         recorders instead of a single one (difference between composite and
#'         single-recorder curves)
#' }
#'
#' **Algorithm:**
#' For each site (recorder), the function:
#' \enumerate{
#'   \item Extracts unique detection dates from file timestamps
#'   \item For each day d = 1 to n, calculates cumulative species count
#'   \item Applies random permutation (n_permutations times) to estimate variance
#'   \item Calculates mean, SE, and SD across permutations
#' }
#'
#' For the composite (multi-site) curve:
#' \enumerate{
#'   \item Pools detections from all selected sites (or n_sites subset)
#'   \item Creates site-day combinations as sampling units
#'   \item Calculates cumulative species using rarefaction approach
#' }
#'
#' @param detections_data Data frame containing acoustic detections with columns:
#'   \describe{
#'     \item{file_path}{Character. Path to audio file containing timestamp}
#'     \item{site}{Character. Site/recorder identifier (e.g., "S01", "S02")}
#'     \item{common_name or species_code}{Character. Species identifier}
#'     \item{confidence}{Numeric. Detection confidence score (0-1)}
#'   }
#' @param species_col Character. Name of column containing species identifiers.
#'   Default: "common_name"
#' @param confidence_threshold Numeric. Minimum confidence score to include 
#'   detections. Default: 0.9
#' @param n_permutations Integer. Number of random permutations for variance 
#'   estimation. Default: 100. Higher values give more stable SE estimates but
#'   increase computation time.
#' @param n_sites_composite Integer or NULL. Number of sites to use for composite
#'   curve. If NULL, uses all available sites. Default: 10
#' @param max_days Integer or NULL. Maximum number of days to include. If NULL,
#'   uses 7 days (1 week) for realistic accumulation. Default: 7
#' @param use_weekly_windows Logical. If TRUE, samples days from within weekly
#'   windows to provide more realistic accumulation. Default: TRUE
#' @param seed Integer. Random seed for reproducibility of permutations.
#'   Default: 42
#'
#' @return A list with class "species_accumulation" containing:
#'   \describe{
#'     \item{single_recorder}{Data frame with columns: day, mean_species, se, 
#'           sd, lower_ci, upper_ci - statistics across all single recorders}
#'     \item{composite}{Data frame with same structure for multi-site composite}
#'     \item{site_curves}{List of data frames, one per site, with raw 
#'           accumulation data}
#'     \item{alpha_diversity}{Numeric. Mean species on day 1 (alpha-diversity)}
#'     \item{gamma_diversity}{Numeric. Difference between composite and single
#'           recorder at maximum days (spatial beta component)}
#'     \item{accumulation_rate}{Numeric. Initial slope of accumulation curve
#'           (temporal beta component)}
#'     \item{parameters}{List with input parameters for reproducibility}
#'     \item{n_sites_used}{Integer. Number of sites included in analysis}
#'     \item{n_days_used}{Integer. Number of days in accumulation curves}
#'     \item{total_species}{Integer. Total unique species detected}
#'   }
#'
#' @references
#' Colwell, R. K., et al. (2012). Models and estimators linking individual-based
#' and sample-based rarefaction, extrapolation and comparison of assemblages.
#' Journal of Plant Ecology, 5(1), 3-21.
#'
#' Gotelli, N. J., & Colwell, R. K. (2001). Quantifying biodiversity: procedures 
#' and pitfalls in the measurement and comparison of species richness. 
#' Ecology Letters, 4(4), 379-391.
#'
#' @seealso
#' \code{\link{plot_species_accumulation}} for visualization
#' \code{\link{calculate_temporal_rarefaction}} for weekly richness estimates
#' \code{\link[vegan]{specaccum}} for alternative accumulation methods
#'
#' @examples
#' \dontrun{
#' # Calculate accumulation for katydids (7-day window)
#' katydid_accum <- calculate_species_accumulation(
#'   detections_data = katydid_detections,
#'   species_col = "common_name",
#'   confidence_threshold = 0.9,
#'   max_days = 7,  # 1 week window
#'   n_sites_composite = 10
#' )
#'
#' # View diversity components
#' cat("Alpha-diversity:", katydid_accum$alpha_diversity, "\n")
#' cat("Gamma contribution:", katydid_accum$gamma_diversity, "\n")
#' }
#'
#' @export
calculate_species_accumulation <- function(detections_data,
                                           species_col = "common_name",
                                           confidence_threshold = 0.9,
                                           n_permutations = 100,
                                           n_sites_composite = 10,
                                           max_days = 7,
                                           use_weekly_windows = TRUE,
                                           seed = 42) {
  
  set.seed(seed)
  
  cat("Calculating species accumulation curves...\n")
  
  # Validate inputs
  if (!species_col %in% names(detections_data)) {
    stop(sprintf("Column '%s' not found in detections_data", species_col))
  }
  
  if (!"site" %in% names(detections_data)) {
    stop("Column 'site' required in detections_data")
  }
  
  if (!"file_path" %in% names(detections_data)) {
    stop("Column 'file_path' required for timestamp extraction")
  }
  
  # Apply confidence filtering if confidence column exists
  if ("confidence" %in% names(detections_data)) {
    detections_data <- detections_data %>%
      filter(confidence >= confidence_threshold)
    cat(sprintf("  Applied confidence threshold >= %.2f\n", confidence_threshold))
  }
  
  # Extract date from file_path (format: YYYYMMDD_HHMMSSZ)
  detections_data <- detections_data %>%
    mutate(
      timestamp_str = str_extract(file_path, "\\d{8}_\\d{6}Z?"),
      detection_date = as.Date(substr(timestamp_str, 1, 8), format = "%Y%m%d")
    ) %>%
    filter(!is.na(detection_date))
  
  if (nrow(detections_data) == 0) {
    stop("No valid detections with extractable timestamps")
  }
  
  # Get unique sites
  all_sites <- unique(detections_data$site)
  n_all_sites <- length(all_sites)
  cat(sprintf("  Sites available: %d\n", n_all_sites))
  
  # Determine sites for composite
  if (is.null(n_sites_composite) || n_sites_composite >= n_all_sites) {
    sites_for_composite <- all_sites
  } else {
    sites_for_composite <- sample(all_sites, n_sites_composite)
  }
  n_sites_used <- length(sites_for_composite)
  cat(sprintf("  Sites for composite curve: %d\n", n_sites_used))
  
  # Calculate single-recorder curves for each site
  site_curves <- list()
  
  for (site_id in all_sites) {
    site_data <- detections_data %>%
      filter(site == site_id)
    
    unique_dates <- sort(unique(site_data$detection_date))
    n_days <- length(unique_dates)
    
    if (n_days < 2) next
    
    # Run permutations
    perm_results <- matrix(NA, nrow = n_days, ncol = n_permutations)
    
    for (perm in 1:n_permutations) {
      # Randomize day order
      date_order <- sample(unique_dates)
      cumulative_species <- character(0)
      
      for (d in 1:n_days) {
        day_species <- site_data %>%
          filter(detection_date == date_order[d]) %>%
          pull(!!sym(species_col)) %>%
          unique()
        
        cumulative_species <- unique(c(cumulative_species, day_species))
        perm_results[d, perm] <- length(cumulative_species)
      }
    }
    
    # Calculate statistics
    site_curves[[site_id]] <- data.frame(
      day = 1:n_days,
      mean_species = rowMeans(perm_results),
      se = apply(perm_results, 1, function(x) sd(x) / sqrt(length(x))),
      sd = apply(perm_results, 1, sd)
    )
  }
  
  if (length(site_curves) == 0) {
    stop("No sites with sufficient data for accumulation analysis")
  }
  
  # Combine single-recorder curves (average across sites)
  max_common_days <- min(sapply(site_curves, nrow))
  if (!is.null(max_days)) {
    max_common_days <- min(max_common_days, max_days)
  }
  
  # Create matrix of species counts per site per day
  species_matrix <- matrix(NA, nrow = max_common_days, ncol = length(site_curves))
  
  for (i in seq_along(site_curves)) {
    species_matrix[, i] <- site_curves[[i]]$mean_species[1:max_common_days]
  }
  
  # Average single-recorder curve
  single_recorder <- data.frame(
    day = 1:max_common_days,
    mean_species = rowMeans(species_matrix),
    se = apply(species_matrix, 1, function(x) sd(x) / sqrt(length(x))),
    sd = apply(species_matrix, 1, sd)
  ) %>%
    mutate(
      lower_ci = mean_species - 1.96 * se,
      upper_ci = mean_species + 1.96 * se
    )
  
  # Calculate composite curve (combining multiple sites)
  composite_data <- detections_data %>%
    filter(site %in% sites_for_composite)
  
  # Create site-day combinations as sampling units
  site_days <- composite_data %>%
    distinct(site, detection_date) %>%
    arrange(detection_date, site)
  
  n_site_days <- nrow(site_days)
  cat(sprintf("  Total site-days for composite: %d\n", n_site_days))
  
  # Run permutations for composite
  max_composite_days <- min(n_site_days, max_common_days * n_sites_used)
  if (!is.null(max_days)) {
    max_composite_days <- min(max_composite_days, max_days * n_sites_used)
  }
  
  # Use equivalent recorder-days for fair comparison
  comp_days_equiv <- ceiling(max_composite_days / n_sites_used)
  comp_days_equiv <- min(comp_days_equiv, max_common_days)
  
  composite_perm <- matrix(NA, nrow = comp_days_equiv, ncol = n_permutations)
  
  for (perm in 1:n_permutations) {
    # Shuffle site-day combinations
    shuffled_sd <- site_days[sample(1:n_site_days), ]
    cumulative_species <- character(0)
    
    # Sample n_sites_used site-days per "equivalent day"
    for (d in 1:comp_days_equiv) {
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
        
        cumulative_species <- unique(c(cumulative_species, day_species))
      }
      
      composite_perm[d, perm] <- length(cumulative_species)
    }
  }
  
  # Composite statistics
  composite <- data.frame(
    day = 1:comp_days_equiv,
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
  
  # Gamma: difference at comparable day
  compare_day <- min(nrow(single_recorder), nrow(composite))
  gamma_diversity <- composite$mean_species[compare_day] - 
                     single_recorder$mean_species[compare_day]
  
  # Temporal beta: initial accumulation rate (slope of first 5 days)
  fit_days <- min(5, nrow(single_recorder))
  if (fit_days >= 2) {
    lm_fit <- lm(mean_species ~ day, 
                 data = single_recorder[1:fit_days, ])
    accumulation_rate <- coef(lm_fit)[2]
  } else {
    accumulation_rate <- NA
  }
  
  # Total species
  total_species <- n_distinct(detections_data[[species_col]])
  
  cat(sprintf("  Analysis complete:\n"))
  cat(sprintf("    - Days analyzed: %d\n", max_common_days))
  cat(sprintf("    - Total species: %d\n", total_species))
  cat(sprintf("    - Alpha-diversity (day 1): %.1f species\n", alpha_diversity))
  cat(sprintf("    - Gamma contribution: %.1f species\n", gamma_diversity))
  
  result <- list(
    single_recorder = single_recorder,
    composite = composite,
    site_curves = site_curves,
    alpha_diversity = alpha_diversity,
    gamma_diversity = gamma_diversity,
    accumulation_rate = accumulation_rate,
    parameters = list(
      species_col = species_col,
      confidence_threshold = confidence_threshold,
      n_permutations = n_permutations,
      n_sites_composite = n_sites_composite,
      max_days = max_days,
      seed = seed
    ),
    n_sites_used = n_sites_used,
    n_days_used = max_common_days,
    total_species = total_species
  )
  
  class(result) <- c("species_accumulation", "list")
  return(result)
}


#' Plot Species Accumulation Curves with Alpha-Beta-Gamma Visualization
#'
#' @title Publication-Quality Species Accumulation Plot
#'
#' @description
#' Creates a publication-quality visualization of species accumulation curves
#' showing the relationship between sampling effort (recorder-days) and species
#' detection, with annotation of alpha, beta, and gamma diversity components.
#' Design inspired by temporal diversity decomposition figures in ecological
#' monitoring literature.
#'
#' @details
#' The plot displays:
#' \itemize{
#'   \item **Solid blue line**: Average single recorder curve with SE envelope
#'   \item **Dashed black line**: Composite multi-recorder curve with SE envelope
#'   \item **Alpha-diversity annotation**: Initial species count (day 1)
#'   \item **Temporal beta-diversity**: Rate of species accumulation (curve shape)
#'   \item **Gamma-diversity annotation**: Benefit from spatial replication
#' }
#'
#' The visualization follows design principles from:
#' \itemize{
#'   \item Minimal, clean aesthetic suitable for peer-reviewed journals
#'   \item Clear distinction between curve types using line style and color
#'   \item Informative legend with component definitions
#'   \item Optional annotations for diversity metrics
#' }
#'
#' @param accum_result Object of class "species_accumulation" from 
#'   \code{\link{calculate_species_accumulation}}
#' @param taxon_name Character. Name of taxonomic group for title/labels.
#'   Default: "species"
#' @param show_annotations Logical. Whether to show alpha/beta/gamma annotations.
#'   Default: TRUE
#' @param show_legend Logical. Whether to show legend. Default: TRUE
#' @param legend_position Character. Position of legend ("right", "bottom", 
#'   "top", "left", "none"). Default: "right"
#' @param colors Named character vector with colors for "single" and "composite".
#'   Default: c(single = "#0072B2", composite = "#000000")
#' @param line_width Numeric. Width of main lines. Default: 1.2
#' @param se_alpha Numeric. Transparency for SE ribbons (0-1). Default: 0.2
#' @param point_size Numeric. Size of points on curves. Default: 2.5
#' @param title Character or NULL. Plot title. If NULL, generates automatic title.
#'   Default: NULL
#' @param subtitle Character or NULL. Plot subtitle. Default: NULL
#' @param x_label Character. X-axis label. Default: "Number of recorder-days analyzed"
#' @param y_label Character. Y-axis label. Default: "Total species detected"
#' @param base_size Numeric. Base font size for theme. Default: 12
#' @param output_file Character or NULL. If provided, saves plot to file.
#'   Default: NULL
#' @param width Numeric. Width in inches for saved plot. Default: 10
#' @param height Numeric. Height in inches for saved plot. Default: 6
#' @param dpi Integer. Resolution for saved plot. Default: 300
#'
#' @return A ggplot2 object that can be further customized or saved
#'
#' @references
#' Design inspired by temporal diversity visualizations in:
#' Symes, L. B., et al. (2022). Ecological methods for acoustic monitoring.
#' Ecology and Evolution.
#'
#' @seealso
#' \code{\link{calculate_species_accumulation}} for data preparation
#'
#' @examples
#' \dontrun
#' # Basic plot
#' p <- plot_species_accumulation(katydid_accum, taxon_name = "katydid")
#' print(p)
#'
#' # Customized plot
#' p <- plot_species_accumulation(
#'   bird_accum,
#'   taxon_name = "bird",
#'   colors = c(single = "darkblue", composite = "darkred"),
#'   show_annotations = TRUE,
#'   output_file = "bird_accumulation.jpg"
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line geom_point annotate
#'             labs theme_minimal theme scale_color_manual scale_linetype_manual
#'             scale_x_continuous scale_y_continuous ggsave element_text
#' @importFrom dplyr mutate
#'
#' @export
plot_species_accumulation <- function(accum_result,
                                      taxon_name = "species",
                                      show_annotations = TRUE,
                                      show_legend = TRUE,
                                      legend_position = "right",
                                      colors = c(single = "#0072B2", 
                                                 composite = "#000000"),
                                      line_width = 1.2,
                                      se_alpha = 0.2,
                                      point_size = 2.5,
                                      title = NULL,
                                      subtitle = NULL,
                                      x_label = "Number of recorder-days analyzed",
                                      y_label = "Total species detected",
                                      base_size = 12,
                                      output_file = NULL,
                                      width = 10,
                                      height = 6,
                                      dpi = 300) {
  
  # Validate input
  if (!inherits(accum_result, "species_accumulation")) {
    warning("Input may not be from calculate_species_accumulation()")
  }
  
  # Extract data
  single_df <- accum_result$single_recorder %>%
    mutate(curve_type = "single")
  
  composite_df <- accum_result$composite %>%
    mutate(curve_type = "composite")
  
  # Combine data for plotting
  plot_data <- bind_rows(single_df, composite_df) %>%
    mutate(curve_type = factor(curve_type, 
                               levels = c("single", "composite")))
  
  # Create legend labels
  n_sites <- accum_result$n_sites_used
  legend_labels <- c(
    single = sprintf("Average single recorder (+/- SE)"),
    composite = sprintf("Composite of %d recorders (+/- ~SE)", n_sites)
  )
  
  # Generate automatic title if not provided
  if (is.null(title)) {
    title <- sprintf("Total %s species detected", taxon_name)
  }
  
  # Build plot
  p <- ggplot(plot_data, aes(x = day, y = mean_species, 
                             color = curve_type, 
                             linetype = curve_type,
                             fill = curve_type)) +
    # SE ribbons
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), 
                alpha = se_alpha, 
                color = NA) +
    # Main lines
    geom_line(linewidth = line_width) +
    # Points (optional - dotted for SE bounds)
    geom_line(aes(y = lower_ci), linetype = "dotted", 
              linewidth = line_width * 0.5, alpha = 0.7) +
    geom_line(aes(y = upper_ci), linetype = "dotted", 
              linewidth = line_width * 0.5, alpha = 0.7) +
    # Scales
    scale_color_manual(
      values = colors,
      labels = legend_labels,
      name = NULL
    ) +
    scale_fill_manual(
      values = colors,
      labels = legend_labels,
      name = NULL
    ) +
    scale_linetype_manual(
      values = c(single = "solid", composite = "dashed"),
      labels = legend_labels,
      name = NULL
    ) +
    scale_x_continuous(
      breaks = function(x) pretty(x, n = 8),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      expand = c(0.02, 0)
    ) +
    # Labels
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label
    ) +
    # Theme
    theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = if (show_legend) legend_position else "none",
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5)
    )
  
  # Add annotations if requested
  if (show_annotations) {
    # Get coordinates for annotations
    max_day <- max(plot_data$day)
    alpha_y <- accum_result$alpha_diversity
    gamma_contrib <- accum_result$gamma_diversity
    
    single_final <- filter(single_df, day == max(day))$mean_species
    composite_final <- filter(composite_df, day == max(day))$mean_species
    
    # Alpha annotation (at day 1)
    p <- p + 
      annotate("point", x = 1, y = alpha_y, 
               color = colors["single"], size = point_size * 1.5, shape = 1) +
      annotate("text", x = 2.5, y = alpha_y * 0.7, 
               label = "alpha-diversity", 
               color = colors["single"], size = 3.5, fontface = "italic",
               hjust = 0)
    
    # Temporal beta annotation (arrow along curve)
    mid_day <- ceiling(max_day / 3)
    if (mid_day >= 2) {
      single_mid <- filter(single_df, day == mid_day)$mean_species
      p <- p +
        annotate("segment", 
                 x = mid_day - 1, y = single_mid - 2,
                 xend = mid_day + 2, yend = single_mid + 1,
                 arrow = arrow(length = unit(0.2, "cm")),
                 color = "#D55E00", linewidth = 0.8) +
        annotate("text", x = mid_day + 1, y = single_mid - 3,
                 label = "Temporal beta-diversity",
                 color = "#D55E00", size = 3.5, fontface = "italic",
                 hjust = 0)
    }
    
    # Gamma annotation (difference between curves)
    if (gamma_contrib > 0) {
      gamma_day <- ceiling(max_day * 0.8)
      single_at_gamma <- filter(single_df, day == gamma_day)$mean_species
      composite_at_gamma <- filter(composite_df, day == min(gamma_day, max(composite_df$day)))$mean_species
      
      p <- p +
        annotate("segment",
                 x = gamma_day, y = single_at_gamma,
                 xend = gamma_day, yend = composite_at_gamma,
                 color = "#009E73", linewidth = 0.8,
                 arrow = arrow(length = unit(0.15, "cm"), ends = "both")) +
        annotate("text", x = gamma_day + 0.5, 
                 y = (single_at_gamma + composite_at_gamma) / 2,
                 label = "gamma-diversity",
                 color = "#009E73", size = 3.5, fontface = "italic",
                 hjust = 0)
    }
  }
  
  # Add definition box (like in the reference image)
  if (show_legend && show_annotations) {
    # Create annotation box text
    definitions <- paste0(
      "alpha-diversity: expected species from\n",
      "  one recorder on one day.\n",
      "beta-diversity: rate of new species\n",
      "  detection over time.\n",
      "gamma-diversity: benefit from using\n",
      "  multiple recorders."
    )
    
    # Position in bottom right
    p <- p +
      annotate("label", 
               x = max(plot_data$day) * 0.95,
               y = max(plot_data$mean_species) * 0.35,
               label = definitions,
               hjust = 1, vjust = 0,
               size = 2.8,
               fill = "white",
               color = "gray30",
               label.padding = unit(0.4, "lines"),
               label.r = unit(0.1, "lines"))
  }
  
  # Save if output_file provided
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, dpi = dpi)
    cat(sprintf("Plot saved to: %s\n", output_file))
  }
  
  return(p)
}


#' Generate Species Accumulation Analysis for Both Birds and Katydids
#'
#' @title Wrapper Function for Dual-Taxon Accumulation Analysis
#'
#' @description
#' Convenience function that calculates and plots species accumulation curves
#' for both bird and katydid datasets, producing publication-ready comparative
#' figures.
#'
#' @param bird_detections Data frame with bird detection data
#' @param katydid_detections Data frame with katydid detection data
#' @param bird_confidence Numeric. Confidence threshold for birds. Default: 0.9
#' @param katydid_confidence Numeric. Confidence threshold for katydids. Default: 0.9
#' @param n_sites Integer. Number of sites for composite curve. Default: 10
#' @param output_dir Character. Directory for saving outputs. Default: "."
#' @param save_plots Logical. Whether to save plots to files. Default: TRUE
#'
#' @return List containing:
#'   \describe{
#'     \item{bird_accumulation}{Species accumulation results for birds}
#'     \item{katydid_accumulation}{Species accumulation results for katydids}
#'     \item{bird_plot}{ggplot object for bird accumulation}
#'     \item{katydid_plot}{ggplot object for katydid accumulation}
#'   }
#'
#' @seealso
#' \code{\link{calculate_species_accumulation}}
#' \code{\link{plot_species_accumulation}}
#'
#' @export
analyze_dual_taxon_accumulation <- function(bird_detections,
                                            katydid_detections,
                                            bird_confidence = 0.9,
                                            katydid_confidence = 0.9,
                                            n_sites = 10,
                                            output_dir = ".",
                                            save_plots = TRUE) {
  
  cat("\n=== SPECIES ACCUMULATION ANALYSIS ===\n\n")
  
  results <- list()
  
  # Bird analysis
  if (!is.null(bird_detections) && nrow(bird_detections) > 0) {
    cat("Processing bird detections...\n")
    
    tryCatch({
      results$bird_accumulation <- calculate_species_accumulation(
        detections_data = bird_detections,
        species_col = "common_name",
        confidence_threshold = bird_confidence,
        n_sites_composite = n_sites
      )
      
      output_file <- if (save_plots) {
        file.path(output_dir, "bird_species_accumulation.jpg")
      } else NULL
      
      results$bird_plot <- plot_species_accumulation(
        results$bird_accumulation,
        taxon_name = "bird",
        output_file = output_file
      )
      
      cat("Bird accumulation analysis complete.\n\n")
      
    }, error = function(e) {
      warning(sprintf("Bird accumulation analysis failed: %s", e$message))
    })
  }
  
  # Katydid analysis
  if (!is.null(katydid_detections) && nrow(katydid_detections) > 0) {
    cat("Processing katydid detections...\n")
    
    tryCatch({
      results$katydid_accumulation <- calculate_species_accumulation(
        detections_data = katydid_detections,
        species_col = "common_name",
        confidence_threshold = katydid_confidence,
        n_sites_composite = n_sites
      )
      
      output_file <- if (save_plots) {
        file.path(output_dir, "katydid_species_accumulation.jpg")
      } else NULL
      
      results$katydid_plot <- plot_species_accumulation(
        results$katydid_accumulation,
        taxon_name = "katydid",
        colors = c(single = "#009E73", composite = "#000000"),
        output_file = output_file
      )
      
      cat("Katydid accumulation analysis complete.\n\n")
      
    }, error = function(e) {
      warning(sprintf("Katydid accumulation analysis failed: %s", e$message))
    })
  }
  
  # Summary
  cat("=== ACCUMULATION SUMMARY ===\n")
  
  if (!is.null(results$bird_accumulation)) {
    cat(sprintf("Birds:\n"))
    cat(sprintf("  - Alpha-diversity: %.1f species\n", 
                results$bird_accumulation$alpha_diversity))
    cat(sprintf("  - Gamma contribution: %.1f species\n",
                results$bird_accumulation$gamma_diversity))
    cat(sprintf("  - Total species: %d\n",
                results$bird_accumulation$total_species))
  }
  
  if (!is.null(results$katydid_accumulation)) {
    cat(sprintf("Katydids:\n"))
    cat(sprintf("  - Alpha-diversity: %.1f species\n", 
                results$katydid_accumulation$alpha_diversity))
    cat(sprintf("  - Gamma contribution: %.1f species\n",
                results$katydid_accumulation$gamma_diversity))
    cat(sprintf("  - Total species: %d\n",
                results$katydid_accumulation$total_species))
  }
  
  return(results)
}


# ==============================================================================
# END OF TEMPORAL_ANALYSIS.R
# ==============================================================================
