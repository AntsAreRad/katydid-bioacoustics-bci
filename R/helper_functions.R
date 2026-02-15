# HELPER FUNCTIONS FOR KATYDID BIOACOUSTICS PROJECT
#
# Project: Comparing Bioacoustics vs DNA Metabarcoding for Katydid Monitoring
# Site: Barro Colorado Island (BCI), Panama
# Internship: M1 IMABEE - Katydid Bioacoustics 2024-2025, M2 IMABEE 2026
# Supervisors: Laurel Symes (Cornell), Yves Basset (STRI), Greg Lamarre (STRI)
#
# This file contains utility functions used across the analysis pipeline:
# - Deployment identification and validation
# - Data quality checks
# - BOLD API integration for taxonomic enrichment
# - File structure testing
#
# Author: Leon Brouille
#

# Required packages
library(tidyverse)  # Data manipulation and piping
library(stringr)    # String operations
library(httr)       # HTTP requests for BOLD API
library(xml2)       # XML parsing for BOLD responses


# Deployments identification function

#' Identify Available BirdNET Deployments
#'
#' Scans a directory to identify valid BirdNET deployment folders. 
#' Valid deployments follow the pattern "Dep_XX" where XX is a number.
#' This function is used to locate processed BirdNET output folders before
#' batch processing.
#'
#' @param directory Character string. Path to the directory containing 
#'   deployment folders.
#'
#' @return A data.frame with three columns:
#'   \itemize{
#'     \item deployment: Character. Standardized deployment name (e.g., "Dep_01")
#'     \item full_path: Character. Complete path to the deployment folder
#'     \item exists: Logical. Whether the folder exists (always TRUE in output)
#'   }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Validates that the input directory exists
#'   \item Searches for folders matching the pattern "Dep_\\d+"
#'   \item Returns only folders that physically exist on disk
#'   \item Prints the number of deployments found to console
#' }
#'
#' @examples
#' \dontrun{
#' # Identify BirdNET deployments
#' bird_deployments <- identify_deployments(
#'   "C:/data/BirdNET_output"
#' )
#' print(bird_deployments)
#' }
#'
#' @seealso \code{\link{identify_katydid_deployments}} for katydid-specific version
#'
#' @export

identify_deployments <- function(directory) {
  # Validate directory existence
  if (!dir.exists(directory)) {
    stop(paste("Directory does not exist:", directory))
  }
  
  # List all subdirectories
  deployment_dirs <- list.dirs(directory, recursive = FALSE, full.names = TRUE)
  deployment_names <- basename(deployment_dirs)
  
  # Filter for valid deployment pattern (Dep_XX)
  valid_deployments <- deployment_names[grepl("^Dep_\\d+$", deployment_names)]
  
  # Create deployment information dataframe
  deployment_info <- data.frame(
    deployment = valid_deployments,
    full_path = file.path(directory, valid_deployments),
    exists = file.exists(file.path(directory, valid_deployments)),
    stringsAsFactors = FALSE
  ) %>%
    filter(exists)
  
  cat(paste("Number of BirdNET deployments found:", nrow(deployment_info), "\n"))
  return(deployment_info)
}


#' Identify Available Katydid Deployments
#'
#' Scans a directory to identify valid katydid acoustic detection folders.
#' Katydid deployments follow a different naming pattern than BirdNET:
#' "88260_BCIPXX_FLAC" where XX is the plot number.
#'
#' @param directory Character string. Path to the directory containing 
#'   katydid deployment folders.
#'
#' @return A data.frame with three columns:
#'   \itemize{
#'     \item deployment: Character. Standardized name (e.g., "Dep_07")
#'     \item full_path: Character. Complete path to the deployment folder
#'     \item exists: Logical. Whether the folder exists (always TRUE in output)
#'   }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Looks for folders matching pattern "88260_BCIP\\d+_FLAC"
#'   \item Extracts plot numbers and converts to standard "Dep_XX" format
#'   \item Returns only existing folders
#'   \item Handles the specific folder structure of Koogu katydid detections
#' }
#'
#' @examples
#' \dontrun{
#' # Identify katydid deployments
#' katydid_deployments <- identify_katydid_deployments(
#'   "C:/data/Koogu_Katydid_detections"
#' )
#' print(katydid_deployments)
#' }
#'
#' @seealso \code{\link{identify_deployments}} for BirdNET version
#'
#' @export

identify_katydid_deployments <- function(directory) {
  # Validate directory existence
  if (!dir.exists(directory)) {
    stop(paste("Directory does not exist:", directory))
  }
  
  # List all subdirectories
  deployment_dirs <- list.dirs(directory, recursive = FALSE, full.names = TRUE)
  deployment_names <- basename(deployment_dirs)
  
  # Filter for katydid deployment pattern: 88260_BCIP##_FLAC
  valid_deployments <- deployment_names[grepl("88260_BCIP\\d+_FLAC", deployment_names)]
  
  # Initialize deployment info dataframe
  deployment_info <- data.frame(
    deployment = character(),
    full_path = character(),
    exists = logical(),
    stringsAsFactors = FALSE
  )
  
  # Extract deployment numbers and standardize names
  for (dep_name in valid_deployments) {
    # Extract the plot number (e.g., "07" from "88260_BCIP07_FLAC")
    dep_number <- str_extract(dep_name, "(?<=BCIP)\\d+")
    standard_name <- paste0("Dep_", dep_number)
    
    deployment_info <- rbind(deployment_info, data.frame(
      deployment = standard_name,
      full_path = file.path(directory, dep_name),
      exists = file.exists(file.path(directory, dep_name)),
      stringsAsFactors = FALSE
    ))
  }
  
  # Filter for existing deployments only
  deployment_info <- deployment_info %>%
    filter(exists)
  
  cat(paste("Number of katydid deployments found:", nrow(deployment_info), "\n"))
  return(deployment_info)
}



# Processed deployment checking function

#' Check Already Processed BirdNET Deployments
#'
#' Identifies which BirdNET deployments have already been processed by checking
#' for summary files in the output directory. This enables resuming interrupted
#' analyses without reprocessing completed deployments.
#'
#' @param output_dir Character string. Path to the output directory where
#'   processed deployment results are stored.
#'
#' @return Character vector of deployment names that have been processed.
#'   Returns empty character vector if no deployments processed or directory
#'   doesn't exist.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Looks for files matching "*_summary.csv" in deployments/ subfolder
#'   \item Extracts deployment names from summary filenames
#'   \item Prints count of processed deployments to console
#'   \item Used to skip already-processed data in batch operations
#' }
#'
#' @examples
#' \dontrun{
#' # Check which deployments are already processed
#' processed <- check_processed_deployments("output/results")
#' print(paste("Already processed:", paste(processed, collapse = ", ")))
#' }
#'
#' @seealso \code{\link{check_processed_katydid_deployments}}
#'
#' @export

check_processed_deployments <- function(output_dir) {
  
  deployment_output_dir <- file.path(output_dir, "deployments")
  
  # Return empty if directory doesn't exist
  if (!dir.exists(deployment_output_dir)) {
    return(character(0))
  }
  
  # Find all summary files
  summary_files <- list.files(deployment_output_dir, 
                               pattern = "_summary\\.csv$", 
                               full.names = FALSE)
  
  # Extract deployment names from filenames
  processed_deployments <- gsub("_summary\\.csv$", "", summary_files)
  
  if (length(processed_deployments) > 0) {
    cat(sprintf("BirdNET deployments already processed: %d\n", 
                length(processed_deployments)))
  }
  
  return(processed_deployments)
}


#' Check Already Processed Katydid Deployments
#'
#' Identifies which katydid deployments have already been processed by checking
#' for katydid summary files in the output directory. Enables resuming 
#' interrupted katydid analyses.
#'
#' @param output_dir Character string. Path to the output directory where
#'   processed katydid deployment results are stored.
#'
#' @return Character vector of deployment names that have been processed.
#'   Returns empty character vector if no deployments processed or directory
#'   doesn't exist.
#'
#' @details
#' Similar to \code{check_processed_deployments} but specifically for katydid data:
#' \itemize{
#'   \item Looks in katydid_deployments/ subfolder
#'   \item Searches for "*_katydid_summary.csv" files
#'   \item Extracts and returns standardized deployment names
#' }
#'
#' @examples
#' \dontrun{
#' # Check processed katydid deployments
#' processed_katydids <- check_processed_katydid_deployments("output/results")
#' print(processed_katydids)
#' }
#'
#' @seealso \code{\link{check_processed_deployments}}
#'
#' @export

check_processed_katydid_deployments <- function(output_dir) {
  
  deployment_output_dir <- file.path(output_dir, "katydid_deployments")
  
  # Return empty if directory doesn't exist
  if (!dir.exists(deployment_output_dir)) {
    return(character(0))
  }
  
  # Find all katydid summary files
  summary_files <- list.files(deployment_output_dir, 
                               pattern = "_katydid_summary\\.csv$", 
                               full.names = FALSE)
  
  # Extract deployment names from filenames
  processed_deployments <- gsub("_katydid_summary\\.csv$", "", summary_files)
  
  if (length(processed_deployments) > 0) {
    cat(sprintf("Katydid deployments already processed: %d\n", 
                length(processed_deployments)))
  }
  
  return(processed_deployments)
}



# File testing functions

#' Test Corrected Processing on Example File
#'
#' Diagnostic function to test katydid file processing on a single example file.
#' Useful for validating detection thresholds and species extraction before
#' running full batch processing.
#'
#' @param example_file Character string. Path to an example .selections.txt file
#'   from Koogu katydid detection output.
#'
#' @return Logical. TRUE if file was processed successfully, FALSE otherwise.
#'   Also prints detailed diagnostic information to console including species
#'   counts at different confidence thresholds.
#'
#' @details
#' The function:
#' \itemize{
#'   \item Validates file existence
#'   \item Reads the tab-delimited detection file
#'   \item Filters out "background" detections
#'   \item Lists unique species with detection counts and score ranges
#'   \item Tests multiple confidence thresholds (0.5, 0.7, 0.9)
#' }
#'
#' @examples
#' \dontrun{
#' # Test a single katydid detection file
#' test_corrected_processing("path/to/example.selections.txt")
#' }
#'
#' @seealso \code{\link{process_single_katydid_deployment_FIXED}}
#'
#' @export
test_corrected_processing <- function(example_file) {
  
  cat("Testing corrected processing on example file...\n")
  
  if (!file.exists(example_file)) {
    cat("[ERROR] Example file not found\n")
    return(FALSE)
  }
  
  # Read and analyze the example file
  katydid_data <- read.delim(example_file, sep = "\t", header = TRUE, 
                             stringsAsFactors = FALSE)
  
  cat(sprintf("File: %s\n", basename(example_file)))
  cat(sprintf("Total selections: %d\n", nrow(katydid_data)))
  
  # Analyze species
  species_tags <- katydid_data$Tags[katydid_data$Tags != "background" & 
                                      !is.na(katydid_data$Tags) & 
                                      katydid_data$Tags != ""]
  
  unique_species <- unique(species_tags)
  cat(sprintf("Species found: %d\n", length(unique_species)))
  
  for (species in unique_species) {
    count <- sum(species_tags == species)
    scores <- katydid_data$Score[katydid_data$Tags == species]
    cat(sprintf("  - %s: %d detections (scores: %.2f - %.2f)\n", 
                species, count, min(scores), max(scores)))
  }
  
  # Test with different confidence thresholds
  for (threshold in c(0.5, 0.7, 0.9)) {
    high_conf_species <- katydid_data %>%
      filter(Tags != "background" & Score >= threshold) %>%
      pull(Tags) %>%
      unique()
    
    cat(sprintf("  At confidence %.1f: %d species\n", threshold, length(high_conf_species)))
  }
  
  return(TRUE)
}


# Data quality check functions

#' Check Data Quality of Analysis Results
#'
#' Performs comprehensive quality checks on integrated analysis results.
#' Validates data dimensions, species counts, site coverage, and method
#' overlap. Essential for ensuring data integrity before publication.
#'
#' @param results List. Complete results object from main_integrated_analysis()
#'   containing bird_data, katydid_data, metabarcoding_data, vegetation_data,
#'   and method_comparison components.
#'
#' @return NULL (invisibly). Prints quality check results to console.
#'
#' @details
#' Quality checks performed:
#' \itemize{
#'   \item Bird data: Number of sites and species
#'   \item Katydid data: Number of sites and species
#'   \item Metabarcoding data: Number of sites and species
#'   \item Vegetation data: Number of plots
#'   \item Method overlap: Species detected by both methods
#'   \item Warnings if overlap is absent or very low (< 10%)
#' }
#'
#' @examples
#' \dontrun{
#' # Run quality check on analysis results
#' results <- main_integrated_analysis(...)
#' check_data_quality(results)
#' }
#'
#' @seealso \code{\link{check_data_quality_with_glm}}
#'
#' @export

check_data_quality <- function(results) {
  
  cat("\n=== DAata quality check ===\n")
  
  # Check bird data
  if (!is.null(results$bird_data)) {
    bird_sites <- length(unique(results$bird_data$presence_matrix$site))
    bird_species <- ncol(results$bird_data$presence_matrix) - 1
    cat(sprintf("[OK] Bird data: %d sites, %d species\n", bird_sites, bird_species))
  }
  
  # Check katydid data
  if (!is.null(results$katydid_data)) {
    katydid_sites <- length(unique(results$katydid_data$presence_matrix$site))
    katydid_species <- ncol(results$katydid_data$presence_matrix) - 1
    cat(sprintf("[OK] Katydid data: %d sites, %d species\n", 
                katydid_sites, katydid_species))
  }
  
  # Check metabarcoding data
  if (!is.null(results$metabarcoding_data)) {
    meta_sites <- length(unique(results$metabarcoding_data$presence_matrix$site))
    meta_species <- ncol(results$metabarcoding_data$presence_matrix) - 1
    cat(sprintf("[OK] Metabarcoding data: %d sites, %d species\n", 
                meta_sites, meta_species))
  }
  
  # Check vegetation data
  if (!is.null(results$vegetation_data)) {
    veg_plots <- nrow(results$vegetation_data$summary)
    cat(sprintf("[OK] Vegetation data: %d plots\n", veg_plots))
  }
  
  # Check method overlap
  if (!is.null(results$method_comparison)) {
    overlap_count <- length(results$method_comparison$both_methods)
    total_species <- length(results$method_comparison$only_acoustic) + 
      length(results$method_comparison$only_metabar) + 
      overlap_count
    
    if (overlap_count == 0) {
      cat("[WARNING] No species overlap between methods\n")
    } else if (overlap_count / total_species < 0.1) {
      cat("[WARNING] Very low species overlap (< 10%)\n")
    } else {
      cat(sprintf("[OK] Species overlap: %d species (%.1f%%)\n", 
                  overlap_count, overlap_count / total_species * 100))
    }
  }
  
  cat("=== Quality check complete! ===\n\n")
}


#' Check Data Quality Including GLM Results
#'
#' Extended quality check that includes validation of GLM (Generalized Linear
#' Model) analysis results. Calls the base quality check function and adds
#' GLM-specific validation.
#'
#' @param results List. Complete results object including glm_analysis component.
#'
#' @return NULL (invisibly). Prints quality check results including GLM
#'   summary to console.
#'
#' @details
#' Performs all checks from \code{check_data_quality} plus:
#' \itemize{
#'   \item Number of significant GLM relationships found
#'   \item Top 3 most significant relationships with p-values
#'   \item Useful for validating statistical analysis completeness
#' }
#'
#' @examples
#' \dontrun{
#' # Run extended quality check with GLM validation
#' results <- main_integrated_analysis(...)
#' results$glm_analysis <- analyze_all_glm_relationships(results, output_dir)
#' check_data_quality_with_glm(results)
#' }
#'
#' @references
#' Wood, S.N. (2017). Generalized Additive Models: An Introduction with R (2nd ed.).
#'   Chapman and Hall/CRC. https://doi.org/10.1201/9781315370279
#'
#' @seealso \code{\link{check_data_quality}}, \code{\link{analyze_all_glm_relationships}}
#'
#' @export

check_data_quality_with_glm <- function(results) {
  
  # Call base quality check function
  check_data_quality(results)
  
  # Add GLM-specific checks
  if (!is.null(results$glm_analysis)) {
    cat(sprintf("[OK] GLM analysis: %d significant relationships found\n", 
                length(results$glm_analysis)))
    
    # Display top 3 most significant relationships
    if (length(results$glm_analysis) > 0) {
      cat("Top 3 most significant relationships:\n")
      glm_summary <- data.frame(
        Relationship = names(results$glm_analysis),
        P_value = sapply(results$glm_analysis, function(x) x$p_value),
        stringsAsFactors = FALSE
      ) %>%
        arrange(P_value) %>%
        head(3)
      
      for (i in 1:nrow(glm_summary)) {
        cat(sprintf("  %d. %s (p = %.4f)\n", 
                    i, glm_summary$Relationship[i], glm_summary$P_value[i]))
      }
    }
  }
  
  cat("=== GLM quality check complete! ===\n\n")
}



# BOLD API functions (Taxonomic Enrichment)

#' Get Taxonomic Information from BOLD Database
#'
#' Retrieves taxonomic information for a given BIN (Barcode Index Number) from
#' the BOLD (Barcode of Life Data) Systems API. Used to enrich DNA metabarcoding
#' results with complete taxonomic classifications.
#'
#' @param bin_id Character string. The BIN identifier (e.g., "BOLD:AAA1234" or
#'   "AAA1234"). The "BOLD:" prefix is optional and will be removed automatically.
#'
#' @return A data.frame with one row and five columns:
#'   \itemize{
#'     \item bin_id: Original BIN identifier
#'     \item species_name: Species name if available (NA otherwise)
#'     \item genus_name: Genus name if available (NA otherwise)
#'     \item family_name: Family name if available (NA otherwise)
#'     \item order_name: Order name if available (NA otherwise)
#'   }
#'
#' @details
#' The function:
#' \itemize{
#'   \item Makes HTTP GET request to BOLD API
#'   \item Parses XML response
#'   \item Handles timeouts (10 second limit)
#'   \item Returns NA values for missing taxonomic levels
#'   \item Catches and logs errors without stopping execution
#' }
#'
#' @examples
#' \dontrun{
#' # Get taxonomic info for a single BIN
#' bin_info <- get_bin_info("BOLD:AAA1234")
#' print(bin_info)
#'
#' # Handle multiple BINs with error checking
#' bins <- c("BOLD:AAA1234", "BOLD:AAB5678")
#' results <- lapply(bins, function(bin) {
#'   Sys.sleep(1)  # Respect API rate limits
#'   get_bin_info(bin)
#' })
#' }
#'
#' @references
#' Ratnasingham, S., & Hebert, P.D.N. (2007). BOLD: The Barcode of Life Data System.
#'   Molecular Ecology Notes 7(3): 355-364. https://doi.org/10.1111/j.1471-8286.2007.01678.x
#'
#' @seealso \code{\link{enrich_with_bin_info}} for batch processing multiple BINs
#'
#' @export

get_bin_info <- function(bin_id) {
  
  # Clean BIN ID (remove "BOLD:" prefix if present)
  clean_bin <- gsub("BOLD:", "", bin_id)
  
  # Construct BOLD API URL
  url <- paste0("http://www.boldsystems.org/index.php/API_Tax/TaxonData?dataTypes=basic&taxId=", 
                clean_bin)
  
  tryCatch({
    # Make API request with timeout
    response <- GET(url, timeout(10))
    
    if (status_code(response) == 200) {
      # Parse XML response
      content_text <- content(response, "text", encoding = "UTF-8")
      xml_data <- read_xml(content_text)
      
      # Extract taxonomic information
      taxon_nodes <- xml_find_all(xml_data, ".//taxon")
      
      if (length(taxon_nodes) > 0) {
        # Extract info from first taxon node
        taxon <- taxon_nodes[1]
        
        species_name <- xml_text(xml_find_first(taxon, ".//species_name"))
        genus_name <- xml_text(xml_find_first(taxon, ".//genus_name"))
        family_name <- xml_text(xml_find_first(taxon, ".//family_name"))
        order_name <- xml_text(xml_find_first(taxon, ".//order_name"))
        
        return(data.frame(
          bin_id = bin_id,
          species_name = ifelse(is.na(species_name) || species_name == "", 
                                NA, species_name),
          genus_name = ifelse(is.na(genus_name) || genus_name == "", 
                              NA, genus_name),
          family_name = ifelse(is.na(family_name) || family_name == "", 
                               NA, family_name),
          order_name = ifelse(is.na(order_name) || order_name == "", 
                              NA, order_name),
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # Return empty record if no data found
    return(data.frame(
      bin_id = bin_id,
      species_name = NA,
      genus_name = NA,
      family_name = NA,
      order_name = NA,
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    # Log error and return empty record
    cat(sprintf("Error retrieving BIN %s: %s\n", bin_id, e$message))
    return(data.frame(
      bin_id = bin_id,
      species_name = NA,
      genus_name = NA,
      family_name = NA,
      order_name = NA,
      stringsAsFactors = FALSE
    ))
  })
}


#' Enrich Orthoptera Data with BOLD Taxonomic Information
#'
#' Batch processes multiple BINs to retrieve complete taxonomic classifications
#' from BOLD. Creates an enhanced species identification column by combining
#' original data with BOLD information. Respects API rate limits.
#'
#' @param orthoptera_data Data.frame containing metabarcoding results with a
#'   'bin_uri' column of BIN identifiers.
#' @param rate_limit_seconds Numeric. Delay in seconds between API requests
#'   to respect BOLD rate limits. Default: 1 second.
#'
#' @return Data.frame. Original data enriched with additional columns:
#'   \itemize{
#'     \item species_name: Species name from BOLD
#'     \item genus_name: Genus name from BOLD
#'     \item family_name: Family name from BOLD
#'     \item order_name: Order name from BOLD
#'     \item enhanced_species: Best available species identification
#'     \item id_source: Source of identification (Original_data, BOLD_species, 
#'       BOLD_genus, or Family_only)
#'   }
#'
#' @details
#' Enhancement strategy (in order of priority):
#' \enumerate{
#'   \item Use original species name if present
#'   \item Use BOLD species name if available
#'   \item Use BOLD genus + "sp." if only genus available
#'   \item Fall back to family level
#' }
#'
#' Progress indicators show:
#' \itemize{
#'   \item Number of unique BINs to process
#'   \item Progress updates every 10 BINs
#'   \item Summary statistics of identification improvements
#'   \item List of newly identified species
#' }
#'
#' @examples
#' \dontrun{
#' # Load metabarcoding data
#' meta_data <- read_excel("metabarcoding_results.xlsx")
#' orthoptera <- meta_data %>% filter(order == "Orthoptera")
#'
#' # Enrich with BOLD data (this takes time!)
#' enriched <- enrich_with_bin_info(orthoptera, rate_limit_seconds = 2)
#'
#' # Compare identification rates
#' table(enriched$id_source)
#' }
#'
#' @references
#' Ratnasingham, S. (2019). mBRAVE: A Barcode Reference library and Analytics 
#'   Visualization Environment. Biodiversity Information Science and Standards 
#'   3: e37313. https://doi.org/10.3897/biss.3.37313
#'
#' @seealso \code{\link{get_bin_info}}
#'
#' @export

enrich_with_bin_info <- function(orthoptera_data, rate_limit_seconds = 1) {
  
  cat("Enriching data with BOLD BIN information...\n")
  
  # Get unique BINs to process
  unique_bins <- unique(orthoptera_data$bin_uri)
  unique_bins <- unique_bins[!is.na(unique_bins) & unique_bins != ""]
  
  cat(sprintf("Searching information for %d unique BINs...\n", length(unique_bins)))
  
  # Initialize results dataframe
  bin_info_df <- data.frame()
  
  # Process each BIN with progress indicators
  for (i in seq_along(unique_bins)) {
    bin_id <- unique_bins[i]
    
    # Show progress every 10 BINs
    if (i %% 10 == 0 || i == length(unique_bins)) {
      cat(sprintf("Progress: %d/%d BINs processed...\n", i, length(unique_bins)))
    }
    
    # Retrieve BIN information
    bin_info <- get_bin_info(bin_id)
    bin_info_df <- rbind(bin_info_df, bin_info)
    
    # Respect API rate limits
    if (i < length(unique_bins)) {
      Sys.sleep(rate_limit_seconds)
    }
  }
  
  # Join BOLD info with original data
  enriched_data <- orthoptera_data %>%
    left_join(bin_info_df, by = c("bin_uri" = "bin_id"))
  
  # Create enhanced species column
  enriched_data <- enriched_data %>%
    mutate(
      enhanced_species = case_when(
        # 1. Original species name if present
        !is.na(species) & species != "" ~ species,
        
        # 2. BOLD species name if available
        !is.na(species_name) & species_name != "" ~ species_name,
        
        # 3. BOLD genus if available
        !is.na(genus_name) & genus_name != "" ~ paste(genus_name, "sp."),
        
        # 4. Fall back to family
        TRUE ~ family_name
      ),
      
      # Track identification source
      id_source = case_when(
        !is.na(species) & species != "" ~ "Original_data",
        !is.na(species_name) & species_name != "" ~ "BOLD_species",
        !is.na(genus_name) & genus_name != "" ~ "BOLD_genus",
        TRUE ~ "Family_only"
      )
    )
  
  # Calculate improvement statistics
  improvement_stats <- enriched_data %>%
    count(id_source) %>%
    mutate(percentage = round(n / nrow(enriched_data) * 100, 1))
  
  cat("\n=== Enrichment results ===\n")
  cat(sprintf("Total records processed: %d\n", nrow(enriched_data)))
  
  for (i in 1:nrow(improvement_stats)) {
    source_type <- improvement_stats$id_source[i]
    count <- improvement_stats$n[i]
    pct <- improvement_stats$percentage[i]
    
    source_desc <- switch(source_type,
                          "Original_data" = "Already identified in original data",
                          "BOLD_species" = "Species identified via BOLD",
                          "BOLD_genus" = "Identified to genus via BOLD",
                          "Family_only" = "Remain at family level")
    
    cat(sprintf("  %s: %d (%.1f%%)\n", source_desc, count, pct))
  }
  
  # Count newly identified species
  newly_identified <- sum(!is.na(enriched_data$species_name) & 
                            (is.na(enriched_data$species) | enriched_data$species == ""))
  
  cat(sprintf("\n[NEW] Species identified via BOLD: %d\n", newly_identified))
  
  # List newly identified species
  if (newly_identified > 0) {
    new_species <- enriched_data %>%
      filter(!is.na(species_name) & (is.na(species) | species == "")) %>%
      select(bin_uri, species_name) %>%
      distinct()
    
    cat("Newly identified species:\n")
    for (i in 1:min(10, nrow(new_species))) {
      cat(sprintf("  %d. %s (%s)\n", i, new_species$species_name[i], 
                  new_species$bin_uri[i]))
    }
    if (nrow(new_species) > 10) {
      cat(sprintf("  ... and %d more\n", nrow(new_species) - 10))
    }
  }
  
  return(enriched_data)
}



# Diagnostic functions

#' Test Katydid Directory Structure
#'
#' Diagnostic function to validate katydid detection folder structure. Helps
#' troubleshoot issues with katydid data loading by examining folder contents
#' and testing deployment identification.
#'
#' @param katydid_dir Character string. Path to the katydid detections directory.
#'
#' @return Logical. TRUE if structure is valid and deployments can be identified,
#'   FALSE otherwise. Also prints detailed diagnostic information to console.
#'
#' @details
#' Diagnostic checks performed:
#' \itemize{
#'   \item Validates directory existence
#'   \item Lists all subdirectories found
#'   \item Counts .selections.txt files in each subdirectory
#'   \item Shows example filenames
#'   \item Tests deployment identification function
#'   \item Reports success or failure with error messages
#' }
#'
#' @examples
#' \dontrun{
#' # Test katydid directory structure
#' is_valid <- test_katydid_structure("C:/data/Koogu_Katydid_detections")
#'
#' if (!is_valid) {
#'   cat("Directory structure issues detected. Check output above.\n")
#' }
#' }
#'
#' @seealso \code{\link{identify_katydid_deployments}}
#'
#' @export

test_katydid_structure <- function(katydid_dir) {
  cat("Testing katydid directory structure...\n")
  
  # Check directory existence
  if (!dir.exists(katydid_dir)) {
    cat("[ERROR] Directory does not exist\n")
    return(FALSE)
  }
  
  # List subdirectories
  subdirs <- list.dirs(katydid_dir, recursive = FALSE)
  cat(sprintf("Found %d subdirectories:\n", length(subdirs)))
  
  for (subdir in subdirs) {
    subdir_name <- basename(subdir)
    cat(sprintf("  - %s\n", subdir_name))
    
    # Count .selections.txt files
    selection_files <- list.files(subdir, 
                                   pattern = "\\.selections\\.txt$", 
                                   recursive = TRUE)
    cat(sprintf("    -> %d .selections.txt files\n", length(selection_files)))
    
    # Show example filenames
    if (length(selection_files) > 0) {
      for (file in selection_files[1:min(3, length(selection_files))]) {
        cat(sprintf("      Example: %s\n", file))
      }
    }
  }
  
  # Test deployment identification
  cat("\nTesting deployment identification:\n")
  tryCatch({
    deployments <- identify_katydid_deployments(katydid_dir)
    cat(sprintf("[OK] Successfully identified %d deployments\n", nrow(deployments)))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("[ERROR] Error identifying deployments: %s\n", e$message))
    return(FALSE)
  })
}
