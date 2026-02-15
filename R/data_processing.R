# Data processing functions
#
# Functions for loading, processing, and transforming bioacoustic and 
# metabarcoding data for katydid biodiversity analysis
#
# Project: Katydid Bioacoustics - BCI 2025-2026
# Context: Comparison of bioacoustic vs DNA metabarcoding methods for 
#          monitoring Orthoptera (katydids) in tropical rainforest
# 
# Location: Barro Colorado Island (BCI), Panama
# Collaboration: STRI (Yves Basset, Greg Lamarre) + Cornell Lab (Laurel Symes)
#
# Author: Leon Brouille (M2 IMABEE)
# Supervisors: Dr. Yves Basset, Dr. Greg Lamarre, Dr. Laurel Symes
# Date: 2025-2026
#
# Required packages:
#   - tidyverse (dplyr, tidyr, readr)
#   - readxl
#   - lubridate
#   - stringr
#

# Required libraries
library(tidyverse)
library(readxl)
library(lubridate)
library(stringr)

# Species-specific confidence thresholds

#' Load species-specific confidence thresholds from CSV
#'
#' Reads a CSV file containing per-species confidence thresholds for BirdNET
#' detections. This allows fine-tuned filtering where some species require
#' higher confidence scores than others (e.g., species frequently confused
#' by the detector, or rare species needing lower thresholds to avoid
#' missing true positives).
#'
#' @param filepath Character. Path to CSV file with species thresholds.
#'   The CSV must contain at minimum:
#'   \itemize{
#'     \item \code{species}: Species identifier (common name matching BirdNET output)
#'     \item \code{confidence_threshold}: Numeric confidence threshold (0-1)
#'   }
#'   Optional columns:
#'   \itemize{
#'     \item \code{species_code}: BirdNET species code (e.g., "bkcchi1")
#'     \item \code{notes}: Justification for the threshold choice
#'     \item \code{source}: Reference for the threshold (e.g., "manual validation")
#'   }
#'
#' @return Data frame with validated species thresholds, or NULL if the file
#'   cannot be loaded. When NULL is returned, the pipeline falls back to
#'   uniform thresholds.
#'
#' @details
#' The CSV is validated upon loading:
#' \itemize{
#'   \item Required columns must be present
#'   \item Thresholds must be numeric and between 0 and 1
#'   \item Duplicate species entries are flagged (last entry kept)
#'   \item Summary statistics are printed for user verification
#' }
#'
#' Species not listed in the CSV will use the default uniform threshold
#' set in \code{CONFIG$bird_confidence}.
#'
#' @examples
#' \dontrun{
#' # Load custom thresholds
#' thresholds <- load_species_thresholds("data/bird_species_thresholds.csv")
#'
#' # Use in pipeline
#' CONFIG$bird_species_thresholds <- thresholds
#' }
#'
#' @seealso \code{\link{apply_species_thresholds}},
#'   \code{\link{generate_threshold_template}}
#'
#' @export

load_species_thresholds <- function(filepath) {
  
  cat(sprintf("\n[*] Loading species-specific thresholds from: %s\n", filepath))
  
  if (!file.exists(filepath)) {
    warning(sprintf("Species thresholds file not found: %s", filepath))
    cat("    Will use uniform confidence threshold for all species\n")
    return(NULL)
  }
  
  tryCatch({
    thresholds <- read.csv(filepath, stringsAsFactors = FALSE, strip.white = TRUE)
    
    # Validate required columns
    if (!"species" %in% colnames(thresholds)) {
      stop("Required column 'species' not found in thresholds file")
    }
    if (!"confidence_threshold" %in% colnames(thresholds)) {
      stop("Required column 'confidence_threshold' not found in thresholds file")
    }
    
    # Coerce and validate threshold values
    thresholds$confidence_threshold <- as.numeric(thresholds$confidence_threshold)
    
    invalid_rows <- is.na(thresholds$confidence_threshold) | 
      thresholds$confidence_threshold < 0 | 
      thresholds$confidence_threshold > 1
    
    if (any(invalid_rows)) {
      n_invalid <- sum(invalid_rows)
      warning(sprintf("%d rows with invalid thresholds (NA or outside 0-1) removed", 
                      n_invalid))
      thresholds <- thresholds[!invalid_rows, ]
    }
    
    # Check for duplicate species (keep last entry with warning)
    if (any(duplicated(thresholds$species))) {
      dup_species <- thresholds$species[duplicated(thresholds$species)]
      warning(sprintf("Duplicate species found: %s. Keeping last entry for each.", 
                      paste(unique(dup_species), collapse = ", ")))
      thresholds <- thresholds %>%
        group_by(species) %>%
        slice_tail(n = 1) %>%
        ungroup()
    }
    
    # Report summary
    cat(sprintf("    [OK] Loaded thresholds for %d species\n", nrow(thresholds)))
    cat(sprintf("    Threshold range: %.2f - %.2f (median: %.2f)\n",
                min(thresholds$confidence_threshold),
                max(thresholds$confidence_threshold),
                median(thresholds$confidence_threshold)))
    
    # Show a few examples
    if (nrow(thresholds) <= 10) {
      for (i in 1:nrow(thresholds)) {
        cat(sprintf("      %s: %.2f\n", 
                    thresholds$species[i], 
                    thresholds$confidence_threshold[i]))
      }
    } else {
      cat(sprintf("    (showing first 5 of %d)\n", nrow(thresholds)))
      for (i in 1:5) {
        cat(sprintf("      %s: %.2f\n", 
                    thresholds$species[i], 
                    thresholds$confidence_threshold[i]))
      }
    }
    
    return(thresholds)
    
  }, error = function(e) {
    warning(sprintf("Error loading species thresholds: %s", e$message))
    cat("    Will use uniform confidence threshold for all species\n")
    return(NULL)
  })
}


#' Apply species-specific confidence thresholds to detection data
#'
#' Filters a detection data frame using per-species confidence thresholds.
#' Species listed in the thresholds table are filtered with their specific
#' threshold; all other species use the default uniform threshold.
#'
#' This function is called internally by \code{process_single_deployment()}
#' and can also be used post-hoc on pre-processed detection data.
#'
#' @param detections Data frame. Detection data with at least a 'confidence'
#'   column and a species identifier column.
#' @param species_thresholds Data frame or NULL. Output of
#'   \code{\link{load_species_thresholds}}. If NULL, falls back to uniform
#'   filtering using \code{default_threshold}.
#' @param default_threshold Numeric. Threshold applied to species not listed
#'   in \code{species_thresholds}, or to all species if
#'   \code{species_thresholds} is NULL. Default: 0.9
#' @param species_col Character. Column name in \code{detections} to match
#'   against species_thresholds$species. Default: "common_name"
#'
#' @return Filtered data frame containing only detections meeting their
#'   respective confidence thresholds. An additional column
#'   \code{threshold_applied} records the threshold used for each detection.
#'
#' @details
#' Filtering logic:
#' \enumerate{
#'   \item If \code{species_thresholds} is NULL: simple uniform filter
#'   \item If provided: join detections with thresholds, use species-specific
#'         where available, default for unlisted species, then filter
#' }
#'
#' @examples
#' \dontrun{
#' # Uniform filtering (no species thresholds)
#' filtered <- apply_species_thresholds(detections, default_threshold = 0.9)
#'
#' # Species-specific filtering
#' thresholds <- load_species_thresholds("bird_thresholds.csv")
#' filtered <- apply_species_thresholds(detections, 
#'                                       species_thresholds = thresholds,
#'                                       default_threshold = 0.9)
#' }
#'
#' @seealso \code{\link{load_species_thresholds}}
#'
#' @export

apply_species_thresholds <- function(detections,
                                     species_thresholds = NULL,
                                     default_threshold = 0.9,
                                     species_col = "common_name") {
  
  if (nrow(detections) == 0) return(detections)
  
  # Case 1: No species-specific thresholds - uniform filtering
  if (is.null(species_thresholds) || nrow(species_thresholds) == 0) {
    filtered <- detections %>%
      filter(confidence >= default_threshold) %>%
      mutate(threshold_applied = default_threshold)
    return(filtered)
  }
  
  # Case 2: Species-specific thresholds
  # Build a lookup: species -> threshold
  threshold_lookup <- setNames(
    species_thresholds$confidence_threshold,
    species_thresholds$species
  )
  
  # Assign threshold per detection: species-specific if available, else default
  detections$threshold_applied <- ifelse(
    detections[[species_col]] %in% names(threshold_lookup),
    threshold_lookup[detections[[species_col]]],
    default_threshold
  )
  
  # Filter each detection against its own threshold
  filtered <- detections %>%
    filter(confidence >= threshold_applied)
  
  return(filtered)
}


#' Generate a template CSV for species-specific thresholds
#'
#' Creates a template CSV file pre-populated with species from existing
#' detection data. Useful as a starting point for defining per-species
#' confidence thresholds based on manual validation.
#'
#' @param detections Data frame. Detection data with species identifiers.
#'   Typically output from a previous pipeline run.
#' @param output_file Character. Path to save the template CSV.
#'   Default: "bird_species_thresholds_template.csv"
#' @param default_threshold Numeric. Default threshold to pre-populate.
#'   Default: 0.9
#' @param species_col Character. Column name for species. Default: "common_name"
#' @param code_col Character. Column name for species codes. Default: "species_code"
#'
#' @return Data frame with the template (also saved as CSV).
#'
#' @details
#' The template includes for each species:
#' \itemize{
#'   \item species: Common name
#'   \item species_code: BirdNET species code
#'   \item confidence_threshold: Pre-filled with default (edit manually)
#'   \item n_detections: Total detections at any confidence level
#'   \item median_confidence: Median confidence score observed
#'   \item notes: Empty column for user annotations
#' }
#'
#' Workflow:
#' 1. Run pipeline once with uniform threshold
#' 2. Generate template from results
#' 3. Review species, adjust thresholds based on validation
#' 4. Re-run pipeline with species-specific thresholds
#'
#' @examples
#' \dontrun{
#' # Generate template from previous run
#' template <- generate_threshold_template(
#'   detections = results$bird_data$detections,
#'   output_file = "data/bird_species_thresholds.csv"
#' )
#' # Then manually edit the CSV and re-run
#' }
#'
#' @seealso \code{\link{load_species_thresholds}}
#'
#' @export

generate_threshold_template <- function(detections,
                                        output_file = "bird_species_thresholds_template.csv",
                                        default_threshold = 0.9,
                                        species_col = "common_name",
                                        code_col = "species_code") {
  
  cat(sprintf("[*] Generating species threshold template...\n"))
  
  if (nrow(detections) == 0) {
    warning("No detections to generate template from")
    return(data.frame())
  }
  
  # Summarize detections per species
  template <- detections %>%
    group_by(across(all_of(c(species_col, code_col)))) %>%
    summarise(
      n_detections = n(),
      median_confidence = round(median(confidence, na.rm = TRUE), 3),
      q25_confidence = round(quantile(confidence, 0.25, na.rm = TRUE), 3),
      q75_confidence = round(quantile(confidence, 0.75, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    rename(species = !!sym(species_col),
           species_code = !!sym(code_col)) %>%
    mutate(
      confidence_threshold = default_threshold,
      notes = "",
      source = "default"
    ) %>%
    arrange(species) %>%
    select(species, species_code, confidence_threshold, 
           n_detections, median_confidence, q25_confidence, q75_confidence,
           notes, source)
  
  # Save template
  write.csv(template, output_file, row.names = FALSE)
  cat(sprintf("    [OK] Template saved: %s (%d species)\n", output_file, nrow(template)))
  cat("    Edit 'confidence_threshold' column, then use with load_species_thresholds()\n")
  
  return(template)
}

#' Process a single BirdNET deployment
#'
#' Processes all detection files (.txt) from a single BirdNET deployment,
#' applying confidence filtering and extracting site information. Uses batch
#' processing to handle large numbers of files efficiently.
#'
#' Supports both uniform confidence thresholds and species-specific thresholds
#' for more nuanced filtering (see \code{\link{load_species_thresholds}}).
#'
#' @param deployment_path Character. Full path to the deployment folder
#' @param deployment_name Character. Name of the deployment (e.g., "Dep_01")
#' @param confidence_threshold Numeric. Default minimum confidence score (0-1).
#'   Used for all species when \code{species_thresholds} is NULL, or for
#'   species not listed in the thresholds file. Default: 0.1
#' @param species_thresholds Data frame or NULL. Species-specific confidence
#'   thresholds loaded via \code{\link{load_species_thresholds}}. If provided,
#'   each species is filtered using its own threshold; species not listed use
#'   the default \code{confidence_threshold}. Default: NULL (uniform threshold).
#' @param max_files_per_batch Integer. Maximum number of files to process in
#'   one batch to avoid memory issues. Default: 3000
#'
#' @return Data frame with columns:
#'   - site: site identifier (extracted from filename)
#'   - deployment: deployment name
#'   - file: source filename
#'   - common_name: species common name
#'   - species_code: species code
#'   - confidence: detection confidence score (0-1)
#'   - start_time: detection start time
#'   - end_time: detection end time
#'
#' @details
#' The function implements:
#' - Batch processing to prevent memory overflow with large datasets
#' - Multiple site name pattern matching (S01_U01, BCIP11_S01, etc.)
#' - Progress reporting for long-running processes
#' - Confidence filtering based on BirdNET score interpretation
#'
#' Site extraction supports patterns:
#' - "S01_U01" format (standard)
#' - "BCIP11_S01" format (BCI prefix)
#'
#' @references
#' Wood, C. M., & Kahl, S. (2024). Guidelines for appropriate use of BirdNET 
#' scores and other detector outputs. Journal of Ornithology, 165(3), 777-782.
#'
#' @seealso \code{\link{combine_all_deployments}}
#'
#' @param species_thresholds Data frame or NULL. Species-specific confidence
#'   thresholds loaded via \code{\link{load_species_thresholds}}. If provided,
#'   each species is filtered using its own threshold; species not listed use
#'   the default \code{confidence_threshold}. Default: NULL (uniform threshold).
#'
#' @export
process_single_deployment <- function(deployment_path, 
                                      deployment_name, 
                                      confidence_threshold = 0.1, 
                                      species_thresholds = NULL,
                                      max_files_per_batch = 3000) {
  
  cat(sprintf("\n - Processing bird deployment: %s \n", deployment_name))
  
  # List all .txt files in deployment
  files <- list.files(path = deployment_path, 
                      pattern = "\\.txt$", 
                      recursive = TRUE, 
                      full.names = TRUE)
  
  cat(sprintf("Number of files found: %d\n", length(files)))
  
  if (length(files) == 0) {
    warning(paste("No files found in", deployment_path))
    return(data.frame())
  }
  
  # Batch processing setup
  total_files <- length(files)
  num_batches <- ceiling(total_files / max_files_per_batch)
  
  cat(sprintf("Processing in %d batch(es) of %d files max\n", 
              num_batches, max_files_per_batch))
  
  all_detections <- data.frame()
  
  # Process each batch
  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * max_files_per_batch + 1
    end_idx <- min(batch_num * max_files_per_batch, total_files)
    batch_files <- files[start_idx:end_idx]
    
    cat(sprintf("  Batch %d/%d: files %d to %d\n", 
                batch_num, num_batches, start_idx, end_idx))
    
    batch_detections <- data.frame()
    
    # Process each file in batch
    for (i in seq_along(batch_files)) {
      file_path <- batch_files[i]
      
      # Progress reporting
      if (i %% 1000 == 0) {
        cat(sprintf("    File %d/%d in batch\n", i, length(batch_files)))
      }
      
      # Extract site name from file path
      site_pattern <- "S(\\d+)_U\\d+"
      site_match <- regexec(site_pattern, file_path)
      
      if (site_match[[1]][1] == -1) {
        # Try alternative pattern
        site_pattern <- "BCIP\\d+_S(\\d+)"
        site_match <- regexec(site_pattern, file_path)
      }
      
      if (site_match[[1]][1] != -1) {
        match_result <- regmatches(file_path, site_match)[[1]]
        site_number <- as.numeric(gsub(".*S(\\d+).*", "\\1", match_result[1]))
        site_id <- sprintf("S%02d", site_number)
        
        # Read detection file
        tryCatch({
          if (file.exists(file_path) && file.size(file_path) > 0) {
            bird_data <- read.delim(file_path, sep = "\t", header = TRUE,
                                    stringsAsFactors = FALSE)
            
            if (nrow(bird_data) > 0) {
              # Check and standardize column names
              required_cols <- c("Species Code", "Common Name", "Confidence")
              alt_cols <- c("Species.Code", "Common.Name", "Confidence")
              
              if (all(required_cols %in% names(bird_data))) {
                names(bird_data)[names(bird_data) %in% required_cols] <- 
                  c("species_code", "common_name", "confidence")
              } else if (all(alt_cols %in% names(bird_data))) {
                names(bird_data)[names(bird_data) %in% alt_cols] <- 
                  c("species_code", "common_name", "confidence")
              } else {
                next
              }
              
              # Standardize time columns
              if ("Begin Time (s)" %in% names(bird_data)) {
                names(bird_data)[names(bird_data) == "Begin Time (s)"] <- "start_time"
              } else if ("Begin.Time..s." %in% names(bird_data)) {
                names(bird_data)[names(bird_data) == "Begin.Time..s."] <- "start_time"
              }
              
              if ("End Time (s)" %in% names(bird_data)) {
                names(bird_data)[names(bird_data) == "End Time (s)"] <- "end_time"
              } else if ("End.Time..s." %in% names(bird_data)) {
                names(bird_data)[names(bird_data) == "End.Time..s."] <- "end_time"
              }
              
              # Select and format data
              detections <- bird_data %>%
                select(species_code, common_name, confidence) %>%
                mutate(site = site_id,
                       file_path = file_path,
                       deployment = deployment_name)
              
              # Apply species-specific or uniform confidence filtering
              detections <- apply_species_thresholds(
                detections, 
                species_thresholds = species_thresholds,
                default_threshold = confidence_threshold,
                species_col = "common_name"
              )
              
              batch_detections <- rbind(batch_detections, detections)
            }
          }
        }, error = function(e) {
          # Ignore errors and continue to next file
        })
      }
    }
    
    all_detections <- rbind(all_detections, batch_detections)
    cat(sprintf("  Batch %d done: %d cumulative detections\n", 
                batch_num, nrow(all_detections)))
    
    # Free memory between batches for large datasets (1 year data)
    gc(verbose = FALSE)
  }
  
  cat(sprintf("\nDeployment %s complete: %d total detections\n", 
              deployment_name, nrow(all_detections)))
  
  return(all_detections)
}


#' Process a single katydid deployment
#'
#' Processes katydid .selections.txt files with
#' proper handling of the Tags/Score format.
#'
#' @param deployment_path Character. Full path to the katydid deployment folder
#' @param deployment_name Character. Name of the deployment
#' @param confidence_threshold Numeric. Minimum confidence score (0-1) for 
#'   including detections. Default: 0.7
#' @param max_files_per_batch Integer. Maximum files per batch. Default: 3000
#'
#' @return Data frame with columns:
#'   - site: site identifier
#'   - deployment: deployment name  
#'   - file: source filename
#'   - common_name: katydid species name (from Tags column)
#'   - species_code: generated species code
#'   - confidence: detection confidence (from Score column)
#'
#' @details
#' The function handles the specific format of katydid annotations:
#' - Reads .selections.txt files (tab-separated)
#' - Extracts species from "Tags" column
#' - Confidence scores from "Score" column
#' - Filters out "background" tags
#' - Site extraction from filename patterns
#'
#' @references
#' Center for Conservation Bioacoustics. (2014). Raven Pro: Interactive Sound 
#' Analysis Software (Version 1.5). The Cornell Lab of Ornithology.
#'
#' @seealso \code{\link{combine_all_katydid_deployments}}
#'
#' @export
process_single_katydid_deployment_FIXED <- function(deployment_path, 
                                                    deployment_name,
                                                    confidence_threshold = 0.7,
                                                    max_files_per_batch = 3000) {
  
  cat(sprintf("\n - Processing katydid deployment: %s \n", deployment_name))
  
  # List all .selections.txt files
  files <- list.files(path = deployment_path,
                      pattern = "\\.selections\\.txt$",
                      recursive = TRUE,
                      full.names = TRUE)
  
  cat(sprintf("Number of .selections.txt files found: %d\n", length(files)))
  
  if (length(files) == 0) {
    warning(paste("No .selections.txt files found in", deployment_path))
    return(data.frame())
  }
  
  # Batch processing setup
  total_files <- length(files)
  num_batches <- ceiling(total_files / max_files_per_batch)
  
  cat(sprintf("Processing in %d batch(es) of %d files max\n", 
              num_batches, max_files_per_batch))
  
  all_detections <- data.frame()
  
  for (batch_num in 1:num_batches) {
    start_idx <- (batch_num - 1) * max_files_per_batch + 1
    end_idx <- min(batch_num * max_files_per_batch, total_files)
    batch_files <- files[start_idx:end_idx]
    
    cat(sprintf("  Batch %d/%d: files %d to %d\n", 
                batch_num, num_batches, start_idx, end_idx))
    
    batch_detections <- data.frame()
    
    for (i in seq_along(batch_files)) {
      file_path <- batch_files[i]
      
      # Progress reporting
      if (i %% 500 == 0) {
        cat(sprintf("    File %d/%d in batch\n", i, length(batch_files)))
      }
      
      # Extract site from filename
      site_pattern <- "S(\\d+)_U\\d+"
      site_match <- regexec(site_pattern, file_path)
      
      if (site_match[[1]][1] != -1) {
        match_result <- regmatches(file_path, site_match)[[1]]
        site_number <- as.numeric(gsub(".*S(\\d+).*", "\\1", match_result[1]))
        site_id <- sprintf("S%02d", site_number)
        
        # Read selections file
        tryCatch({
          if (file.exists(file_path) && file.size(file_path) > 0) {
            katydid_data <- read.delim(file_path, sep = "\t", header = TRUE,
                                       stringsAsFactors = FALSE)
            
            # Check required columns
            required_cols <- c("Selection", "Tags", "Score")
            
            if (all(required_cols %in% colnames(katydid_data))) {
              
              # Filter species detections (exclude background)
              species_detections <- katydid_data %>%
                filter(Tags != "background" & !is.na(Tags) & Tags != "") %>%
                filter(Score >= confidence_threshold)
              
              if (nrow(species_detections) > 0) {
                n_detections <- nrow(species_detections)
                
                # Extract datetime from filename
                # Format: 88260_BCIP12_20240815_143052_FLAC.selections.txt
                datetime_pattern <- "(\\d{8})_(\\d{6})"
                datetime_match <- regexec(datetime_pattern, file_path)
                datetime_result <- regmatches(file_path, datetime_match)[[1]]
                
                if (length(datetime_result) >= 3) {
                  date_str <- datetime_result[2]
                  time_str <- datetime_result[3]
                  
                  datetime_utc <- tryCatch({
                    as.POSIXct(
                      paste0(substr(date_str, 1, 4), "-",
                             substr(date_str, 5, 6), "-",
                             substr(date_str, 7, 8), " ",
                             substr(time_str, 1, 2), ":",
                             substr(time_str, 3, 4), ":",
                             substr(time_str, 5, 6)),
                      tz = "UTC"
                    )
                  }, error = function(e) NA)
                  
                  datetime_panama <- if (!is.na(datetime_utc)) datetime_utc - lubridate::hours(5) else NA
                  date_only <- if (!is.na(datetime_panama)) as.Date(datetime_panama) else NA
                } else {
                  datetime_utc <- NA
                  datetime_panama <- NA
                  date_only <- NA
                }
                
                detections <- data.frame(
                  common_name = species_detections$Tags,
                  species_code = paste0("KATY_", gsub("[^A-Za-z0-9]", "_", 
                                                      species_detections$Tags)),
                  confidence = species_detections$Score,
                  site = rep(site_id, n_detections),
                  file_path = rep(file_path, n_detections),
                  deployment = rep(deployment_name, n_detections),
                  datetime_utc = rep(datetime_utc, n_detections),
                  datetime_panama = rep(datetime_panama, n_detections),
                  date = rep(date_only, n_detections),
                  stringsAsFactors = FALSE
                )
                
                batch_detections <- rbind(batch_detections, detections)
              }
            }
          }
        }, error = function(e) {
          # Continue to next file on error
        })
      }
    }
    
    all_detections <- rbind(all_detections, batch_detections)
    cat(sprintf("  Batch %d done: %d cumulative detections\n", 
                batch_num, nrow(all_detections)))
    
    # Free memory between batches for large datasets (1 year data)
    gc(verbose = FALSE)
  }
  
  # Summary of species found
  if (nrow(all_detections) > 0) {
    unique_species <- unique(all_detections$common_name)
    cat(sprintf("\nSpecies found in %s:\n", deployment_name))
    for (i in seq_along(unique_species)) {
      species_count <- sum(all_detections$common_name == unique_species[i])
      cat(sprintf("  %d. %s: %d detections\n", i, unique_species[i], species_count))
    }
  }
  
  cat(sprintf("Total detections: %d\n", nrow(all_detections)))
  cat(sprintf("Unique species: %d\n", length(unique(all_detections$common_name))))
  
  return(all_detections)
}



# Data combination and aggregation

#' Combine all processed bird deployments
#'
#' Combines results from multiple processed BirdNET deployments into unified
#' datasets. Creates both a complete detection list and a presence/absence
#' matrix across all sites.
#'
#' @param output_dir Character. Directory containing processed deployment results
#'
#' @return List with components:
#'   - detections: combined data frame of all detections
#'   - presence_matrix: site x species presence/absence matrix
#'   - deployment_summary: summary statistics per deployment
#'
#' @details
#' The function:
#' - Reads all *_detections.csv files from deployments/ subdirectory
#' - Combines into single dataset
#' - Creates presence/absence matrix (1 = present, 0 = absent)
#' - Generates summary statistics
#'
#' @seealso \code{\link{process_single_deployment}}
#'
#' @export
combine_all_deployments <- function(output_dir) {
  
  cat("\n=== Combining all bird deployments ===\n")
  
  deployment_dir <- file.path(output_dir, "deployments")
  
  if (!dir.exists(deployment_dir)) {
    stop("Deployments directory not found. Process deployments first.")
  }
  
  # Find all detection files
  detection_files <- list.files(deployment_dir, 
                                pattern = "_detections\\.csv$",
                                full.names = TRUE)
  
  cat(sprintf("Found %d deployment files to combine\n", length(detection_files)))
  
  if (length(detection_files) == 0) {
    stop("No detection files found")
  }
  
  # Combine all detections
  all_detections <- data.frame()
  
  for (file in detection_files) {
    deployment_name <- gsub("_detections\\.csv$", "", basename(file))
    cat(sprintf("Loading %s...\n", deployment_name))
    
    deployment_data <- read.csv(file, stringsAsFactors = FALSE)
    all_detections <- rbind(all_detections, deployment_data)
  }
  
  cat(sprintf("Total detections combined: %d\n", nrow(all_detections)))
  cat(sprintf("Unique species: %d\n", length(unique(all_detections$common_name))))
  cat(sprintf("Unique sites: %d\n", length(unique(all_detections$site))))
  
  # Create presence matrix
  presence_matrix <- all_detections %>%
    select(site, common_name) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = common_name,
                values_from = present,
                values_fill = 0)
  
  # Deployment statistics
  deployment_stats <- all_detections %>%
    group_by(deployment) %>%
    summarise(
      total_detections = n(),
      unique_species = length(unique(common_name)),
      unique_sites = length(unique(site)),
      mean_confidence = mean(confidence, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Save combined results
  write.csv(all_detections, 
            file.path(output_dir, "combined_all_detections.csv"),
            row.names = FALSE)
  write.csv(presence_matrix,
            file.path(output_dir, "combined_presence_matrix.csv"),
            row.names = FALSE)
  write.csv(deployment_stats,
            file.path(output_dir, "bird_deployment_statistics.csv"),
            row.names = FALSE)
  
  cat("\nCombined bird data saved successfully\n")
  
  return(list(
    detections = all_detections,
    presence_matrix = presence_matrix,
    deployment_summary = deployment_stats
  ))
}


#' Combine all processed katydid deployments
#'
#' Combines results from multiple processed katydid deployments into unified
#' datasets. Creates detection lists, presence matrices, and summary statistics
#' specific to katydid acoustic data.
#'
#' @param output_dir Character. Directory containing processed katydid results
#'
#' @return List with components:
#'   - raw_detections: combined data frame of all katydid detections
#'   - presence_matrix: site x species presence/absence matrix
#'   - deployment_summary: summary statistics per deployment
#'   - species_list: vector of all detected species names
#'
#' @seealso \code{\link{process_single_katydid_deployment_FIXED}}
#'
#' @export
combine_all_katydid_deployments <- function(output_dir) {
  
  cat("\n=== Combining all katydid deployments ===\n")
  
  deployment_dir <- file.path(output_dir, "katydid_deployments")
  
  if (!dir.exists(deployment_dir)) {
    stop("Katydid deployments directory not found. Process deployments first.")
  }
  
  # Find all katydid detection files
  detection_files <- list.files(deployment_dir,
                                pattern = "_katydid_detections\\.csv$",
                                full.names = TRUE)
  
  cat(sprintf("Found %d katydid deployment files to combine\n", length(detection_files)))
  
  if (length(detection_files) == 0) {
    stop("No katydid detection files found")
  }
  
  # Combine all detections
  all_detections <- data.frame()
  
  for (file in detection_files) {
    deployment_name <- gsub("_katydid_detections\\.csv$", "", basename(file))
    cat(sprintf("Loading katydid %s...\n", deployment_name))
    
    deployment_data <- read.csv(file, stringsAsFactors = FALSE)
    all_detections <- rbind(all_detections, deployment_data)
  }
  
  cat(sprintf("Total katydid detections combined: %d\n", nrow(all_detections)))
  cat(sprintf("Unique katydid species: %d\n", 
              length(unique(all_detections$common_name))))
  cat(sprintf("Unique sites: %d\n", length(unique(all_detections$site))))
  
  # Create presence matrix
  presence_matrix <- all_detections %>%
    select(site, common_name) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = common_name,
                values_from = present,
                values_fill = 0)
  
  # Deployment statistics
  deployment_stats <- all_detections %>%
    group_by(deployment) %>%
    summarise(
      total_detections = n(),
      unique_species = length(unique(common_name)),
      unique_sites = length(unique(site)),
      mean_confidence = mean(confidence, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Extract species list
  species_list <- sort(unique(all_detections$common_name))
  
  # Save combined results
  write.csv(all_detections,
            file.path(output_dir, "combined_all_katydid_detections.csv"),
            row.names = FALSE)
  write.csv(presence_matrix,
            file.path(output_dir, "combined_katydid_presence_matrix.csv"),
            row.names = FALSE)
  write.csv(deployment_stats,
            file.path(output_dir, "katydid_deployment_statistics.csv"),
            row.names = FALSE)
  
  cat("\nCombined katydid data saved successfully\n")
  
  return(list(
    detections = all_detections,
    presence_matrix = presence_matrix,
    deployment_summary = deployment_stats,
    species_list = species_list
  ))
}



# Katydid data loading

#' Load katydid data with improved structure recognition
#'
#' Flexible function for loading katydid bioacoustic data from either raw
#' .selections.txt files or pre-processed CSV files. Handles directory
#' structure variations and provides fallback to simulated data.
#'
#' @param katydid_dir Character. Path to directory with raw katydid files
#' @param katydid_detection_file Character. Path to pre-processed detections CSV
#' @param katydid_matrix_file Character. Path to pre-processed matrix CSV
#' @param confidence_threshold Numeric. Minimum confidence score. Default: 0.7
#' @param process_katydids_from_scratch Logical. Process raw files? Default: FALSE
#' @param skip_processed_deployments Logical. Skip already processed? Default: TRUE
#' @param output_dir Character. Output directory. Default: "integrated_results"
#'
#' @return List with components:
#'   - raw_detections: detection data frame
#'   - presence_matrix: site x species matrix
#'
#' @details
#' Processing options:
#' 1. From scratch: Process raw .selections.txt files
#' 2. Pre-processed: Load existing CSV files
#' 3. Fallback: Generate simulated data if other options fail
#'
#' @seealso \code{\link{process_single_katydid_deployment_FIXED}}
#'
#' @export
load_katydid_data_improved <- function(katydid_dir = NULL,
                                       katydid_detection_file = NULL, 
                                       katydid_matrix_file = NULL,
                                       confidence_threshold = 0.7, 
                                       min_detection_days = 5,
                                       process_katydids_from_scratch = FALSE,
                                       skip_processed_deployments = TRUE,
                                       output_dir = "integrated_results") {
  
  cat("Loading katydid data with improved structure recognition...\n")
  cat(sprintf("  - katydid_dir: %s\n", ifelse(is.null(katydid_dir), "NULL", katydid_dir)))
  cat(sprintf("  - process_from_scratch: %s\n", process_katydids_from_scratch))
  
  # OPTION 1: Process from scratch
  if (process_katydids_from_scratch && !is.null(katydid_dir)) {
    cat("Processing katydids from raw .selections.txt files...\n")
    
    if (!dir.exists(katydid_dir)) {
      warning(paste("Katydid directory does not exist:", katydid_dir, 
                    "- using simulated data"))
      return(generate_simulated_katydid_data(confidence_threshold))
    }
    
    tryCatch({
      deployments <- identify_katydid_deployments(katydid_dir)
      
      if (nrow(deployments) == 0) {
        warning("No katydid deployments found - using simulated data")
        return(generate_simulated_katydid_data(confidence_threshold))
      }
      
      cat(sprintf("Found %d katydid deployments\n", nrow(deployments)))
      
      # Check processed deployments
      processed_deployments <- check_processed_katydid_deployments(output_dir)
      
      remaining_deployments <- deployments %>%
        filter(!deployment %in% processed_deployments)
      
      if (nrow(remaining_deployments) > 0 || !skip_processed_deployments) {
        
        deployments_to_process <- if (skip_processed_deployments) {
          remaining_deployments
        } else {
          deployments
        }
        
        for (i in 1:nrow(deployments_to_process)) {
          deployment_info <- deployments_to_process[i, ]
          
          cat(sprintf("Processing katydid deployment %d/%d: %s\n", 
                      i, nrow(deployments_to_process), deployment_info$deployment))
          
          detections <- process_single_katydid_deployment_FIXED(
            deployment_path = deployment_info$full_path,
            deployment_name = deployment_info$deployment,
            confidence_threshold = confidence_threshold
          )
          
          save_katydid_deployment_results(detections, deployment_info$deployment, 
                                          output_dir)
        }
      }
      
      # Combine all katydid deployments
      katydid_data <- combine_all_katydid_deployments(output_dir)
      
      # Apply minimum detection days filter (Laurel Symes criterion)
      if (min_detection_days > 0 && "date" %in% colnames(katydid_data$detections)) {
        filtered_detections <- filter_by_detection_days(
          katydid_data$detections, 
          min_days = min_detection_days
        )
        
        # Recalculate presence matrix from filtered detections
        filtered_matrix <- filtered_detections %>%
          select(site, common_name) %>%
          distinct() %>%
          mutate(present = 1) %>%
          pivot_wider(names_from = common_name,
                      values_from = present,
                      values_fill = 0)
        
        return(list(
          raw_detections = filtered_detections,
          presence_matrix = filtered_matrix
        ))
      } else {
        return(list(
          raw_detections = katydid_data$detections,
          presence_matrix = katydid_data$presence_matrix
        ))
      }
      
    }, error = function(e) {
      warning(paste("Error processing katydid directory:", e$message, 
                    "- using simulated data"))
      return(generate_simulated_katydid_data(confidence_threshold))
    })
  }
  
  # OPTION 2: Load from pre-processed files
  if (!is.null(katydid_detection_file) && !is.null(katydid_matrix_file)) {
    cat("Loading pre-processed katydid data files...\n")
    
    if (!file.exists(katydid_detection_file) || !file.exists(katydid_matrix_file)) {
      warning("Katydid files not found - using simulated data")
      return(generate_simulated_katydid_data(confidence_threshold))
    }
    
    tryCatch({
      katydid_detections <- read.csv(katydid_detection_file, stringsAsFactors = FALSE)
      
      # Filter detections by confidence threshold
      filtered_detections <- katydid_detections %>%
        filter(confidence >= confidence_threshold)
      
      # CRITICAL FIX: Recalculate presence matrix from filtered detections
      # The pre-saved matrix contains all species regardless of confidence
      # We must recalculate to ensure only high-confidence detections are included
      filtered_matrix <- filtered_detections %>%
        select(site, common_name) %>%
        distinct() %>%
        mutate(present = 1) %>%
        pivot_wider(names_from = common_name,
                    values_from = present,
                    values_fill = 0)
      
      cat(sprintf("Applied confidence threshold %.2f: %d -> %d detections\n", 
                  confidence_threshold, nrow(katydid_detections), nrow(filtered_detections)))
      
      # Apply minimum detection days filter if date column exists
      if (min_detection_days > 0 && "date" %in% colnames(filtered_detections)) {
        filtered_detections <- filter_by_detection_days(
          filtered_detections, 
          min_days = min_detection_days
        )
        
        # Recalculate presence matrix from filtered detections
        filtered_matrix <- filtered_detections %>%
          select(site, common_name) %>%
          distinct() %>%
          mutate(present = 1) %>%
          pivot_wider(names_from = common_name,
                      values_from = present,
                      values_fill = 0)
      }
      
      cat(sprintf("Species in matrix after all filtering: %d (sites: %d)\n", 
                  ncol(filtered_matrix) - 1, nrow(filtered_matrix)))
      
      return(list(
        raw_detections = filtered_detections,
        presence_matrix = filtered_matrix
      ))
    }, error = function(e) {
      warning(paste("Error reading katydid files:", e$message, "- using simulated data"))
      return(generate_simulated_katydid_data(confidence_threshold))
    })
  }
  
  # OPTION 3: Simulated data (fallback)
  warning("No valid katydid data provided - using simulated data")
  return(generate_simulated_katydid_data(confidence_threshold))
}



# Detection days filtering

#' Filter species by minimum detection days
#'
#' Applies criterion: a species must be detected on at least
#' a minimum number of distinct days at a site to be considered present.
#' This reduces false positives from sporadic/erroneous detections.
#'
#' @param detections Data frame with columns: site, common_name, date
#' @param min_days Integer. Minimum distinct days required. Default: 5
#'
#' @return Filtered data frame containing only detections for species-site
#'   combinations meeting the minimum days criterion
#'
#' @details
#' The function:
#' - Groups detections by site and species
#' - Counts distinct dates for each site-species combination
#' - Filters to keep only combinations with >= min_days distinct dates
#' - Returns the filtered original detections (not aggregated)
#'
#' This criterion was recommended by Dr. Laurel Symes (Cornell Lab) to ensure
#' robust species presence determination in bioacoustic monitoring.
#'
#' @examples
#' \dontrun{
#' filtered <- filter_by_detection_days(katydid_detections, min_days = 5)
#' }
#'
#' @export
filter_by_detection_days <- function(detections, min_days = 5) {
  
  cat(sprintf("\nApplying %d-day minimum detection criterion...\n", min_days))
  
  if (!"date" %in% colnames(detections)) {
    warning("Column 'date' not found in detections - cannot apply minimum days filter")
    warning("Returning unfiltered detections")
    return(detections)
  }
  
  # Count distinct detection days per site-species combination
  detection_days <- detections %>%
    filter(!is.na(date)) %>%
    group_by(site, common_name) %>%
    summarise(
      n_days = n_distinct(date),
      n_detections = n(),
      .groups = "drop"
    )
  
  # Report statistics before filtering
  cat(sprintf("  Site-species combinations before filter: %d\n", nrow(detection_days)))
  cat(sprintf("  Detection days range: %d - %d\n", 
              min(detection_days$n_days), max(detection_days$n_days)))
  
  # Filter to keep only combinations meeting criterion
  valid_combinations <- detection_days %>%
    filter(n_days >= min_days)
  
  cat(sprintf("  Site-species combinations after filter: %d\n", nrow(valid_combinations)))
  
  # Filter original detections to keep only valid combinations
  filtered_detections <- detections %>%
    semi_join(valid_combinations, by = c("site", "common_name"))
  
  # Report filtering results
  n_removed_detections <- nrow(detections) - nrow(filtered_detections)
  pct_removed <- 100 * n_removed_detections / nrow(detections)
  
  species_before <- length(unique(detections$common_name))
  species_after <- length(unique(filtered_detections$common_name))
  
  cat(sprintf("  Detections removed: %d (%.1f%%)\n", n_removed_detections, pct_removed))
  cat(sprintf("  Species before/after: %d / %d\n", species_before, species_after))
  
  # Report per-site richness after filtering
  site_richness <- filtered_detections %>%
    group_by(site) %>%
    summarise(richness = n_distinct(common_name), .groups = "drop")
  
  cat(sprintf("  Site richness range after filter: %d - %d species\n",
              min(site_richness$richness), max(site_richness$richness)))
  
  return(filtered_detections)
}



# Metabarcoding data loading

#' Read metabarcoding data from Excel file
#'
#' Loads DNA metabarcoding data from an Excel file and creates presence
#' matrices for Orthoptera species. Handles site name conversion from
#' Plot format to site format.
#'
#' @param metabarcoding_file Character. Path to Excel file with metabarcoding data
#' @param enrich_with_bold Logical. Whether to enrich with BOLD database. 
#'   Default: FALSE
#'
#' @return List with components:
#'   - raw_data: filtered Orthoptera data frame
#'   - presence_matrix: site x species (BIN) presence matrix
#'   - species_matrix: site x named species presence matrix
#'   - species_list: vector of identified species names
#'   - total_taxa_count: total number of unique BINs
#'   - identified_species_count: number of named species
#'
#' @details
#' The function:
#' - Filters for Orthoptera order or related families
#' - Converts Plot names (P01) to site names (S01)
#' - Creates matrices using BIN URIs as taxa identifiers
#' - Optionally enriches data with BOLD taxonomy
#'
#' @seealso \code{\link{enrich_with_bin_info}}
#'
#' @export

read_metabarcoding_data <- function(metabarcoding_file, enrich_with_bold = FALSE) {
  
  cat("Loading metabarcoding data...\n")
  
  if (!file.exists(metabarcoding_file)) {
    stop(paste("Metabarcoding file does not exist:", metabarcoding_file))
  }
  
  meta_data <- read_excel(metabarcoding_file)
  
  # Filter for Orthoptera
  orthoptera_data <- meta_data %>%
    filter(order == "Orthoptera" | 
             family %in% c("Tettigoniidae", "Gryllidae", "Acrididae", 
                           "Eumastacidae", "Romaleidae", "Gryllacrididae"))
  
  # Optional BOLD enrichment
  if (enrich_with_bold) {
    cat("\nEnriching with BOLD database...\n")
    cat("WARNING: This may take several minutes depending on number of BINs\n")
    orthoptera_data <- enrich_with_bin_info(orthoptera_data)
    species_col <- "enhanced_species"
  } else {
    species_col <- "species"
  }
  
  # Convert Plot to site format
  site_col <- if ("site" %in% colnames(orthoptera_data)) "site" 
  else if ("Plot" %in% colnames(orthoptera_data)) "Plot" 
  else if ("plot" %in% colnames(orthoptera_data)) "plot"
  else NULL
  
  if (is.null(site_col)) {
    stop("Cannot find site/Plot column in metabarcoding data")
  }
  
  if (site_col %in% c("Plot", "plot")) {
    orthoptera_data <- orthoptera_data %>%
      mutate(site = gsub("P", "S", .data[[site_col]]))
  } else {
    orthoptera_data$site <- orthoptera_data[[site_col]]
  }
  
  # Create presence matrices
  # 1. Full matrix (all taxa via bin_uri)
  full_matrix <- orthoptera_data %>%
    group_by(site, bin_uri) %>%
    summarise(present = 1, .groups = "drop") %>%
    pivot_wider(names_from = bin_uri, values_from = present, values_fill = 0)
  
  # 2. Identified species matrix
  identified_species <- orthoptera_data %>%
    filter(!is.na(.data[[species_col]]) & .data[[species_col]] != "") %>%
    group_by(site, !!sym(species_col)) %>%
    summarise(present = 1, .groups = "drop") %>%
    pivot_wider(names_from = !!sym(species_col), values_from = present, values_fill = 0)
  
  # Statistics
  total_taxa <- length(unique(orthoptera_data$bin_uri))
  identified_species_count <- length(unique(orthoptera_data[[species_col]][
    !is.na(orthoptera_data[[species_col]]) & orthoptera_data[[species_col]] != ""]))
  
  cat(sprintf("\n=== Summary ===\n"))
  cat(sprintf("Total taxa (BINs): %d\n", total_taxa))
  cat(sprintf("Identified species: %d\n", identified_species_count))
  
  return(list(
    raw_data = orthoptera_data,
    presence_matrix = full_matrix,
    species_matrix = identified_species,
    species_list = unique(orthoptera_data[[species_col]][
      !is.na(orthoptera_data[[species_col]]) & orthoptera_data[[species_col]] != ""]),
    total_taxa_count = total_taxa,
    identified_species_count = identified_species_count,
    enrichment_used = enrich_with_bold
  ))
}



# Vegetation data loading

#' Read vegetation data from Excel file
#'
#' Loads ForestGEO vegetation census data from an Excel file containing
#' tree and liana information for each plot.
#'
#' @param excel_file Character. Path to Excel file with vegetation data
#'
#' @return List with components:
#'   - summary: data frame with vegetation metrics per plot
#'   - tree_species: site x tree species abundance matrix
#'
#' @details
#' The Excel file must contain:
#' - "Trees" sheet with Genus, SpeciesName, Family, and YB-P* columns
#' - "Liana" sheet with Plot Code, Latin, Rooted Stems, Basal Area
#'
#' Calculated metrics include:
#' - tree_species_richness: number of tree species per plot
#' - tree_abundance: total tree count per plot
#' - tree_total_basal_area: sum of basal areas (or proxy)
#' - liana_species_richness: number of liana species
#' - liana_rooted_stems: total rooted stems
#' - liana_total_basal_area: total liana basal area
#'
#' @export
read_vegetation_data <- function(excel_file) {
  
  cat("Loading vegetation data...\n")
  
  if (!file.exists(excel_file)) {
    stop(paste("Vegetation file does not exist:", excel_file))
  }
  
  # Read tabs with readxl package
  trees_data <- read_excel(excel_file, sheet = "Trees")
  lianas_data <- read_excel(excel_file, sheet = "Liana")
  
  # Identify site columns (YB-P*)
  site_cols <- grep("^YB-P", colnames(trees_data), value = TRUE)
  
  # Transform tree data to long format
  trees_long <- trees_data %>%
    select(Genus, SpeciesName, Family, all_of(site_cols)) %>%
    pivot_longer(
      cols = all_of(site_cols),
      names_to = "Plot",
      values_to = "Count"
    ) %>%
    filter(!is.na(Count) & Count > 0) %>%
    mutate(Latin = paste(Genus, SpeciesName))
  
  # Calculate metrics for each plot
  tree_summary <- trees_long %>%
    group_by(Plot) %>%
    summarise(
      tree_species_richness = n_distinct(Latin),
      tree_abundance = sum(Count),
      tree_total_basal_area = sum(Count * 1000),  # Proxy value
      .groups = "drop"
    )
  
  # Process lianas data
  liana_summary <- lianas_data %>%
    group_by(`Plot Code`) %>%
    summarise(
      liana_species_richness = n_distinct(Latin),
      liana_rooted_stems = sum(`Rooted Stems`, na.rm = TRUE),
      liana_total_basal_area = sum(`Basal Area (mm)`, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Plot = `Plot Code`)
  
  # Merge data as global vegetation summary
  vegetation_summary <- full_join(tree_summary, liana_summary, by = "Plot")
  
  # Create tree species matrix
  tree_species <- trees_long %>%
    group_by(Plot, Latin) %>%
    summarise(abundance = sum(Count), .groups = "drop") %>%
    pivot_wider(names_from = Latin, values_from = abundance, values_fill = 0)
  
  cat(sprintf("Vegetation data loaded for %d plots\n", nrow(vegetation_summary)))
  
  return(list(
    summary = vegetation_summary,
    tree_species = tree_species
  ))
}



# Matrix creation

#' Create katydid presence/absence matrix
#'
#' Converts katydid detection data into a site x species presence/absence matrix
#' suitable for community ecology analyses.
#'
#' @param detections Data frame. Katydid detection data with 'site' and 
#'   'common_name' columns
#'
#' @return Data frame with site as first column and species as remaining columns
#'   (1 = present, 0 = absent)
#'
#' @export
create_katydid_presence_matrix <- function(detections) {
  
  if (nrow(detections) == 0) {
    return(data.frame(site = character()))
  }
  
  presence_matrix <- detections %>%
    select(site, common_name) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = common_name,
                values_from = present,
                values_fill = 0)
  
  return(presence_matrix)
}


#' Create detection count matrix
#'
#' Converts raw detection data into a site x species count matrix.
#'
#' @param raw_detections Data frame. Raw detection data
#' @param species_col Character. Species column name. Default: "common_name"
#' @param site_col Character. Site column name. Default: "site"
#'
#' @return Data frame with site as first column and species counts
#'
#' @export
create_detection_count_matrix <- function(raw_detections, 
                                          species_col = "common_name", 
                                          site_col = "site") {
  
  cat("Creating detection count matrix...\n")
  
  if (nrow(raw_detections) == 0) {
    return(data.frame(site = character()))
  }
  
  detection_counts <- raw_detections %>%
    group_by(across(all_of(c(site_col, species_col)))) %>%
    summarise(detections = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(species_col),
                values_from = detections,
                values_fill = 0)
  
  cat(sprintf("Detection matrix created: %d sites, %d species\n",
              nrow(detection_counts), ncol(detection_counts) - 1))
  cat("Morpheus would be proud\n")
  
  return(detection_counts)
}


# Results saving

#' Save bird deployment results
#'
#' Saves processed bird detection results including raw detections,
#' presence matrix, and summary statistics for a single deployment.
#'
#' @param detections Data frame. Processed bird detections
#' @param deployment_name Character. Name of the deployment
#' @param output_dir Character. Output directory for saving results
#'
#' @return Data frame with summary statistics
#'
#' @export

save_deployment_results <- function(detections, deployment_name, output_dir) {
  
  deployment_output_dir <- file.path(output_dir, "deployments")
  dir.create(deployment_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save raw detections
  detections_file <- file.path(deployment_output_dir, 
                               paste0(deployment_name, "_detections.csv"))
  write.csv(detections, detections_file, row.names = FALSE)
  
  # Create and save presence matrix
  if (nrow(detections) > 0) {
    presence_matrix <- detections %>%
      select(site, common_name) %>%
      distinct() %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = common_name,
                  values_from = present,
                  values_fill = 0)
    
    matrix_file <- file.path(deployment_output_dir,
                             paste0(deployment_name, "_presence_matrix.csv"))
    write.csv(presence_matrix, matrix_file, row.names = FALSE)
    
    # Summary statistics
    summary_stats <- data.frame(
      deployment = deployment_name,
      total_detections = nrow(detections),
      unique_species = length(unique(detections$common_name)),
      unique_sites = length(unique(detections$site)),
      mean_confidence = mean(detections$confidence, na.rm = TRUE),
      processing_time = Sys.time(),
      stringsAsFactors = FALSE
    )
    
    summary_file <- file.path(deployment_output_dir,
                              paste0(deployment_name, "_summary.csv"))
    write.csv(summary_stats, summary_file, row.names = FALSE)
    
    return(summary_stats)
  } else {
    cat(sprintf("No detections for %s - no files saved\n", deployment_name))
    return(data.frame())
  }
}


#' Save katydid deployment results
#'
#' Saves processed katydid detection results including raw detections,
#' presence matrix, and summary statistics for a single deployment.
#'
#' @param detections Data frame. Processed katydid detections
#' @param deployment_name Character. Name of the deployment
#' @param output_dir Character. Output directory for saving results
#'
#' @return Data frame with summary statistics
#'
#' @export

save_katydid_deployment_results <- function(detections, deployment_name, output_dir) {
  
  deployment_output_dir <- file.path(output_dir, "katydid_deployments")
  dir.create(deployment_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Save raw detections
  detections_file <- file.path(deployment_output_dir,
                               paste0(deployment_name, "_katydid_detections.csv"))
  write.csv(detections, detections_file, row.names = FALSE)
  
  # Create and save presence matrix
  if (nrow(detections) > 0) {
    presence_matrix <- detections %>%
      select(site, common_name) %>%
      distinct() %>%
      mutate(present = 1) %>%
      pivot_wider(names_from = common_name,
                  values_from = present,
                  values_fill = 0)
    
    matrix_file <- file.path(deployment_output_dir,
                             paste0(deployment_name, "_katydid_presence_matrix.csv"))
    write.csv(presence_matrix, matrix_file, row.names = FALSE)
    
    # Summary statistics
    summary_stats <- data.frame(
      deployment = deployment_name,
      total_detections = nrow(detections),
      unique_species = length(unique(detections$common_name)),
      unique_sites = length(unique(detections$site)),
      mean_confidence = mean(detections$confidence, na.rm = TRUE),
      processing_time = Sys.time(),
      stringsAsFactors = FALSE
    )
    
    summary_file <- file.path(deployment_output_dir,
                              paste0(deployment_name, "_katydid_summary.csv"))
    write.csv(summary_stats, summary_file, row.names = FALSE)
    
    return(summary_stats)
  } else {
    cat(sprintf("No katydid detections for %s - no files saved\n", deployment_name))
    return(data.frame())
  }
}



# Fallback and utility functions

#' Generate simulated katydid data (fallback function)
#'
#' Creates simulated katydid detection data when real data processing fails (for testing).
#' Used as a fallback to ensure the analysis pipeline continues.
#'
#' @param confidence_threshold Numeric. Minimum confidence threshold. Default: 0.1
#' @param n_sites Integer. Number of sites to simulate. Default: 10
#' @param n_species Integer. Number of species to simulate. Default: 5
#' @param n_detections Integer. Number of detections to generate. Default: 100
#'
#' @return List with raw_detections and presence_matrix components
#'
#' @details
#' WARNING: This is simulated data for testing only.
#'
#' @export

generate_simulated_katydid_data <- function(confidence_threshold = 0.1,
                                            n_sites = 10,
                                            n_species = 5,
                                            n_detections = 100) {
  
  cat("WARNING: Generating SIMULATED katydid data\n")
  cat("   This is fallback data - NOT real detections!\n")
  
  # Common neotropical katydid genus names
  species_pool <- c(
    "Copiphora brevirostris",
    "Anaulacomera spatulata",
    "Thamnobates subfalcata",
    "Euceraia insignis",
    "Dolichocercus latipennis"
  )
  
  # Generate random detections
  simulated_data <- data.frame(
    site = sprintf("S%02d", sample(1:n_sites, n_detections, replace = TRUE)),
    deployment = "SIMULATED",
    file_path = paste0("SIMULATED_", 1:n_detections, ".selections.txt"),
    common_name = sample(species_pool[1:min(n_species, length(species_pool))],
                         n_detections, replace = TRUE),
    species_code = paste0("KATY_SIM_", 1:n_detections),
    confidence = runif(n_detections, confidence_threshold, 1),
    stringsAsFactors = FALSE
  )
  
  # Create presence matrix
  presence_matrix <- simulated_data %>%
    select(site, common_name) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = common_name,
                values_from = present,
                values_fill = 0)
  
  cat(sprintf("Generated %d simulated detections for %d sites and %d species\n",
              nrow(simulated_data), n_sites, n_species))
  
  return(list(
    raw_detections = simulated_data,
    presence_matrix = presence_matrix
  ))
}