
# MAIN ANALYSIS SCRIPT - Integrated Bioacoustic vs Metabarcoding Analysis

# Project: Katydid Bioacoustics - BCI Panama
# Author: Leon Brouille (M2 IMABEE)
# Supervisors: Dr. Yves Basset (STRI), Dr. Greg Lamarre (STRI), 
#              Dr. Laurel Symes (Cornell Lab of Ornithology)
# Date: 2025-2026
# Purpose: Main orchestration script for comprehensive comparative analysis
#          between bioacoustic monitoring and DNA metabarcoding methods

# Required packages
# Note: Source all module files before using this script
# source("functions/helper_functions.R")
# source("scripts/data_processing.R")
# source("scripts/statistical_analysis.R")
# source("scripts/temporal_analysis.R")
# source("scripts/plotting_functions.R")
# source("scripts/species_accumulation_abg.R")  # Alpha-beta-gamma curves

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
  library(readxl)
  library(lubridate)
  library(mgcv)
  library(betapart)
  library(ade4)
  library(VennDiagram)
  library(gridExtra)
  library(corrplot)
})



# Main integrated analysis functions

#' Main Integrated Bioacoustic vs Metabarcoding Analysis Pipeline
#'
#' This is the primary orchestration function that executes a comprehensive
#' comparative analysis between bioacoustic monitoring and DNA metabarcoding
#' methods for biodiversity assessment. It processes data from multiple sources,
#' performs various ecological analyses, and generates publication-ready outputs.
#'
#' @description
#' \code{main_integrated_analysis} orchestrates a complete analytical workflow
#' that integrates bioacoustic monitoring data (BirdNET and Katydid recordings),
#' DNA metabarcoding data, and environmental variables (ForestGEO vegetation data)
#' to comprehensively compare detection methods and analyze community patterns.
#' 
#' This function serves as the central entry point for the entire analysis
#' pipeline, calling specialized functions from multiple modules to:
#' \itemize{
#'   \item Process raw acoustic detection data or load pre-processed files
#'   \item Load and standardize metabarcoding and vegetation data
#'   \item Perform methodological comparisons between detection approaches
#'   \item Analyze ecological relationships and community patterns
#'   \item Generate comprehensive visualizations and statistical outputs
#'   \item Save all results in an organized directory structure
#' }
#'
#' @details
#' \strong{WORKFLOW OVERVIEW}
#' 
#' The function executes 8 major analytical steps:
#' 
#' \strong{STEP 1: Bird Data Processing}
#' \itemize{
#'   \item If \code{process_birds_from_scratch = TRUE}: Processes BirdNET 
#'         detection files from raw data using batch processing
#'   \item Otherwise: Loads pre-processed bird detection and presence matrices
#'   \item Creates detailed species summaries and detection statistics
#'   \item Filters detections based on confidence threshold
#' }
#' 
#' \strong{STEP 2: Katydid Data Processing}
#' \itemize{
#'   \item If \code{process_katydids_from_scratch = TRUE}: Processes Katydid
#'         .selections.txt files from raw acoustic data
#'   \item Otherwise: Loads pre-processed katydid detection and presence matrices
#'   \item Creates detailed species summaries and temporal patterns
#'   \item Applies confidence threshold filtering
#' }
#' 
#' \strong{STEP 3: Metabarcoding Data Loading}
#' \itemize{
#'   \item Loads DNA metabarcoding data from Excel file
#'   \item Extracts Orthoptera species list and BIN (Barcode Index Number) information
#'   \item Creates presence/absence matrices for site-level comparisons
#'   \item Optionally enriches data with BOLD (Barcode of Life Database) taxonomy
#' }
#' 
#' \strong{STEP 4: Vegetation Data Loading}
#' \itemize{
#'   \item Loads ForestGEO vegetation census data
#'   \item Processes tree and liana data for each plot
#'   \item Calculates vegetation metrics (basal area, density, diversity)
#'   \item Prepares environmental matrices for constrained ordinations
#' }
#' 
#' \strong{STEP 5: Methodological Comparisons}
#' \itemize{
#'   \item Compares species detection between bioacoustic and metabarcoding methods
#'   \item Calculates detection overlap and complementarity
#'   \item Analyzes site-level similarity between methods (Mantel tests, co-inertia)
#'   \item Examines temporal patterns in acoustic detections
#'   \item Generates Venn diagrams and comparison statistics
#' }
#' 
#' \strong{STEP 6: Ecological Relationships}
#' \itemize{
#'   \item Analyzes vegetation effects on species richness (correlations, GLMs)
#'   \item Tests relationships between taxonomic groups (Katydids-Birds, Katydids-Orthoptera)
#'   \item Examines environmental predictors of community composition
#'   \item Performs diversity partitioning (alpha, beta, gamma)
#' }
#' 
#' \strong{STEP 7: Multivariate Community Analyses}
#' \itemize{
#'   \item Non-metric multidimensional scaling (NMDS) for each taxonomic group
#'   \item Canonical Correspondence Analysis (CCA) with environmental constraints
#'   \item Environmental vector fitting (envfit) to identify significant predictors
#'   \item Generalized Additive Models (GAM) for non-linear relationships
#'   \item Creates ordination biplots and diagnostic plots
#' }
#' 
#' \strong{STEP 8: Visualization and Results Synthesis}
#' \itemize{
#'   \item Generates comprehensive publication-ready figures
#'   \item Creates comparison plots (bar charts, Venn diagrams, scatter plots)
#'   \item Produces temporal activity patterns visualizations
#'   \item Saves all statistical results as CSV files
#'   \item Compiles executive summary with key findings
#' }
#'
#' \strong{DATA PROCESSING OPTIONS}
#' 
#' The function offers flexible data processing workflows:
#' \itemize{
#'   \item \strong{From scratch}: Process raw acoustic files (time-intensive, 
#'         requires computational resources)
#'   \item \strong{Pre-processed}: Load previously processed detection matrices 
#'         (fast, for re-running analyses)
#'   \item \strong{Batch processing}: Automatically handles large datasets by 
#'         processing deployments in batches
#'   \item \strong{Resume capability}: Can skip already-processed deployments 
#'         to resume interrupted analyses
#' }
#'
#' \strong{OUTPUT STRUCTURE}
#' 
#' Results are saved in the \code{output_dir} with the following organization:
#' \preformatted{
#' output_dir/
#' +-- figures/              # All visualization outputs
#'     +-- comparisons/      # Method comparison plots
#'     +-- environmental/    # Vegetation relationship plots
#'     +-- temporal/         # Temporal pattern visualizations
#' +-- tables/               # Statistical results (CSV)
#'     +-- *_comparison.csv
#'     +-- *_correlations.csv
#'     +-- *_summary.csv
#' +-- models/               # Ordination and model outputs
#'     +-- *_ordination.jpg
#'     +-- *_envfit_results.txt
#' }
#'
#' \strong{PERFORMANCE CONSIDERATIONS}
#' 
#' \itemize{
#'   \item Processing from scratch: 30-120 minutes depending on dataset size
#'   \item Loading pre-processed data: 2-5 minutes
#'   \item Memory usage: 2-8 GB RAM depending on dataset size
#'   \item Parallelization: Currently sequential; batch processing prevents memory issues
#' }
#'
#' \strong{ERROR HANDLING}
#' 
#' The function includes robust error handling:
#' \itemize{
#'   \item Validates all input parameters before processing
#'   \item Checks file existence and format compatibility
#'   \item Provides informative error messages for common issues
#'   \item Continues processing even if optional analyses fail
#'   \item Logs warnings for data quality issues
#' }
#'
#' @param bird_dir Character string. Path to directory containing BirdNET 
#'   detection files (*.BirdNET.results.txt). Required if 
#'   \code{process_birds_from_scratch = TRUE}, otherwise optional.
#' @param bird_detections_file Character string. Path to pre-processed bird 
#'   detections CSV file. Required if \code{process_birds_from_scratch = FALSE}.
#' @param bird_matrix_file Character string. Path to pre-processed bird 
#'   presence/absence matrix CSV. Required if \code{process_birds_from_scratch = FALSE}.
#' @param katydid_dir Character string. Path to directory containing Katydid 
#'   detection files (*.selections.txt). Required if 
#'   \code{process_katydids_from_scratch = TRUE}, otherwise optional.
#' @param katydid_detections_file Character string. Path to pre-processed 
#'   katydid detections CSV file. Required if \code{process_katydids_from_scratch = FALSE}.
#' @param katydid_matrix_file Character string. Path to pre-processed katydid 
#'   presence/absence matrix CSV. Required if \code{process_katydids_from_scratch = FALSE}.
#' @param metabarcoding_file Character string. Path to Excel file containing 
#'   DNA metabarcoding data with Orthoptera species and BIN information. 
#'   Required parameter.
#' @param vegetation_file Character string. Path to Excel file containing 
#'   ForestGEO vegetation census data (trees and lianas). Required parameter.
#' @param confidence_threshold_birds Numeric. Minimum confidence threshold for 
#'   filtering bird detections (0-1). Default: 0.9. Higher values increase 
#'   precision but may reduce detection sensitivity.
#' @param confidence_threshold_katydids Numeric. Minimum confidence threshold 
#'   for filtering katydid detections (0-1). Default: 0.7. Lower than birds 
#'   due to higher acoustic complexity in katydid calls.
#' @param output_dir Character string. Path to directory where all results 
#'   will be saved. Will be created if it doesn't exist. Default: "integrated_results".
#' @param process_birds_from_scratch Logical. If TRUE, processes raw BirdNET 
#'   files; if FALSE, loads pre-processed data. Default: FALSE (faster for 
#'   re-running analyses).
#' @param process_katydids_from_scratch Logical. If TRUE, processes raw Katydid 
#'   detection files; if FALSE, loads pre-processed data. Default: FALSE.
#' @param skip_processed_deployments Logical. If TRUE, skips deployments that 
#'   have already been processed (useful for resuming interrupted analyses). 
#'   Only applicable when processing from scratch. Default: TRUE.
#' @param run_temporal_analysis Logical. If TRUE, performs temporal pattern 
#'   analyses (hourly activity, diel patterns). Default: TRUE. Set to FALSE 
#'   to save time if temporal patterns are not of interest.
#' @param create_plots Logical. If TRUE, generates all visualization outputs. 
#'   Default: TRUE. Set to FALSE for faster execution when only statistics 
#'   are needed.
#' @param save_results Logical. If TRUE, saves all results to disk. 
#'   Default: TRUE. Set to FALSE for testing or when only in-memory results 
#'   are needed.
#'
#' @return A named list containing all analysis results with the following components:
#' \describe{
#'   \item{bird_data}{List with bird detections and presence matrix}
#'   \item{katydid_data}{List with katydid detections and presence matrix}
#'   \item{metabarcoding_data}{List with metabarcoding species data and matrices}
#'   \item{vegetation_data}{List with vegetation metrics and plot summaries}
#'   \item{comprehensive_comparison}{Complete method comparison results including:
#'     \itemize{
#'       \item Species lists for each method
#'       \item Detection overlap statistics
#'       \item Venn diagram data
#'       \item Complementarity metrics
#'     }
#'   }
#'   \item{site_similarity}{Site-level similarity analyses including:
#'     \itemize{
#'       \item Mantel test results (correlation between distance matrices)
#'       \item Co-inertia analysis results
#'       \item Site richness comparisons
#'       \item Correlation coefficients and p-values
#'     }
#'   }
#'   \item{temporal_patterns}{Temporal activity patterns including:
#'     \itemize{
#'       \item Hourly detection counts
#'       \item Diel period classifications (day/night/crepuscular)
#'       \item Peak activity times
#'       \item Species-specific activity patterns
#'     }
#'   }
#'   \item{vegetation_effects_katydids}{Vegetation-katydid relationships including:
#'     \itemize{
#'       \item Correlation coefficients with vegetation metrics
#'       \item GLM model results
#'       \item Scatter plots with regression lines
#'       \item Significance tests
#'     }
#'   }
#'   \item{vegetation_effects_metabar}{Vegetation-metabarcoding relationships}
#'   \item{vegetation_effects_birds}{Vegetation-bird relationships}
#'   \item{katydids_birds_relationship}{Cross-taxa relationship analysis}
#'   \item{katydids_ortho_relationship}{Bioacoustic vs metabarcoding relationship}
#'   \item{multivariate_katydids}{NMDS ordination results for katydids including:
#'     \itemize{
#'       \item Ordination coordinates
#'       \item Stress values
#'       \item Environmental vector fitting results
#'       \item Biplot visualization
#'     }
#'   }
#'   \item{multivariate_metabar}{NMDS ordination results for metabarcoding}
#'   \item{multivariate_birds}{NMDS ordination results for birds}
#'   \item{richness_ranges}{Species richness summary statistics per site}
#'   \item{executive_summary}{Text summary of key findings}
#' }
#'
#' @note
#' \strong{Maybe important considerations:}
#' \itemize{
#'   \item Ensure sufficient RAM (>=8 GB recommended) for large datasets
#'   \item Processing from scratch can take a lot of time; plan accordingly
#'   \item Pre-processed files significantly speed up re-analyses
#'   \item The function creates the output directory if it doesn't exist
#'   \item Interrupted analyses can be resumed using \code{skip_processed_deployments = TRUE}
#'   \item Confidence thresholds critically affect detection rates; default values 
#'         may need adjustments
#'   \item Temporal analyses require timestamp information in detection files
#'   \item Some analyses (e.g., co-inertia) require matching sites between methods
#' }
#' 
#' \strong{Data requirements:}
#' \itemize{
#'   \item BirdNET files must be in format: *.BirdNET.results.txt
#'   \item Katydid files must be in format: *.selections.txt
#'   \item Metabarcoding file must contain sheets: "Orthoptera", "BIN_info"
#'   \item Vegetation file must contain: "Trees" and "Lianas" sheets
#'   \item Site names must be consistent across all data sources
#' }
#'
#' @seealso 
#' \strong{Data Processing Functions:}
#' \itemize{
#'   \item \code{\link{identify_deployments}} - Identify available BirdNET deployments
#'   \item \code{\link{process_single_deployment}} - Process individual BirdNET deployment
#'   \item \code{\link{process_single_katydid_deployment_FIXED}} - Process Katydid deployment
#'   \item \code{\link{combine_all_deployments}} - Combine processed bird data
#'   \item \code{\link{read_metabarcoding_data}} - Load metabarcoding data
#'   \item \code{\link{read_vegetation_data}} - Load vegetation data
#' }
#' 
#' \strong{Statistical Analysis Functions:}
#' \itemize{
#'   \item \code{\link{calculate_diversity_metrics}} - Diversity indices
#'   \item \code{\link{analyze_site_similarity}} - Mantel tests and co-inertia
#'   \item \code{\link{analyze_vegetation_effects}} - Environmental relationships
#'   \item \code{\link{analyze_nmds_with_detection_counts}} - NMDS ordination
#'   \item \code{\link{analyze_cca_with_detections}} - Canonical correspondence analysis
#' }
#' 
#' \strong{Visualization Functions:}
#' \itemize{
#'   \item \code{\link{create_comparison_plots}} - Method comparison visualizations
#'   \item \code{\link{save_all_plots}} - Save plots to organized directories
#' }
#'
#' @references
#' Borcard, D., Gillet, F., & Legendre, P. (2018). \emph{Numerical Ecology with R} 
#' (2nd ed.). Springer. \url{https://doi.org/10.1007/978-3-319-71404-2}
#' 
#' Legendre, P., & Legendre, L. (2012). \emph{Numerical Ecology} (3rd ed.). Elsevier.
#' 
#' Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M. (2009). 
#' \emph{Mixed Effects Models and Extensions in Ecology with R}. Springer. 
#' \url{https://doi.org/10.1007/978-0-387-87458-6}
#' 
#' Hampton, S. E., Jones, M. B., Wasser, L. A., Schildhauer, M. P., Supp, S. R., 
#' Brun, J., ... & Aukema, J. E. (2017). Skills and Knowledge for Data-Intensive 
#' Environmental Research. \emph{BioScience}, 67(6), 546-557. 
#' \url{https://doi.org/10.1093/biosci/bix025}
#' 
#' Dray, S., Chessel, D., & Thioulouse, J. (2003). Co-inertia analysis and the 
#' linking of ecological data tables. \emph{Ecology}, 84(11), 3078-3089. 
#' \url{https://doi.org/10.1890/03-0178}
#' 
#' Kahl, S., Wood, C. M., Eibl, M., & Klinck, H. (2021). BirdNET: A deep learning 
#' solution for avian diversity monitoring. \emph{Ecological Informatics}, 61, 101236. 
#' \url{https://doi.org/10.1016/j.ecoinf.2021.101236}
#' 
#' Ji, Y., Ashton, L., Pedley, S. M., Edwards, D. P., Tang, Y., Nakamura, A., ... 
#' & Yu, D. W. (2013). Reliable, verifiable and efficient monitoring of biodiversity 
#' via metabarcoding. \emph{Ecology Letters}, 16(10), 1245-1257. 
#' \url{https://doi.org/10.1111/ele.12162}
#'
#' @export
#' @importFrom tidyverse %>% filter select mutate
#' @importFrom vegan metaMDS cca envfit mantel
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw
#'
#' @examples
#' \dontrun{
# Example 1: Process everything from scratch (comprehensive but time-intensive)
main_integrated_analysis <- function(bird_dir = NULL,
                                     bird_detections_file = NULL,
                                     bird_matrix_file = NULL,
                                     katydid_dir = NULL,
                                     katydid_detections_file = NULL,
                                     katydid_matrix_file = NULL,
                                     metabarcoding_file,
                                     vegetation_file,
                                     confidence_threshold_birds = 0.9,
                                     confidence_threshold_katydids = 0.9,
                                     bird_species_thresholds = NULL,
                                     min_detection_days = 5,
                                     output_dir = "integrated_results",
                                     process_birds_from_scratch = FALSE,
                                     process_katydids_from_scratch = FALSE,
                                     skip_processed_deployments = TRUE,
                                     run_temporal_analysis = TRUE,
                                     run_statistics = TRUE,
                                     create_plots = TRUE,
                                     save_results = TRUE) {
  #' 
  #' # Example 2: Use pre-processed data (fast, for re-running analyses)
  #' results_fast <- main_integrated_analysis(
  #'   bird_detections_file = "results/bird_detections_combined.csv",
  #'   bird_matrix_file = "results/bird_presence_matrix.csv",
  #'   katydid_detections_file = "results/katydid_detections_combined.csv",
  #'   katydid_matrix_file = "results/katydid_presence_matrix.csv",
  #'   metabarcoding_file = "data/Metabarcoding_Orthoptera.xlsx",
  #'   vegetation_file = "data/ForestGEO_vegetation.xlsx",
  #'   output_dir = "results/reanalysis",
  #'   process_birds_from_scratch = FALSE,
  #'   process_katydids_from_scratch = FALSE
  #' )
  #' 
  #' # Example 3: Resume interrupted analysis
  #' results_resumed <- main_integrated_analysis(
  #'   bird_dir = "data/BirdNET_detections/",
  #'   katydid_dir = "data/Katydid_selections/",
  #'   metabarcoding_file = "data/Metabarcoding_Orthoptera.xlsx",
  #'   vegetation_file = "data/ForestGEO_vegetation.xlsx",
  #'   output_dir = "results/interrupted_analysis",
  #'   process_birds_from_scratch = TRUE,
  #'   process_katydids_from_scratch = TRUE,
  #'   skip_processed_deployments = TRUE  # Skip already processed deployments
  #' )
  #' 
  #' # Example 4: Skip temporal analysis and plotting for faster execution
  #' results_stats_only <- main_integrated_analysis(
  #'   bird_detections_file = "results/bird_detections_combined.csv",
  #'   bird_matrix_file = "results/bird_presence_matrix.csv",
  #'   katydid_detections_file = "results/katydid_detections_combined.csv",
  #'   katydid_matrix_file = "results/katydid_presence_matrix.csv",
  #'   metabarcoding_file = "data/Metabarcoding_Orthoptera.xlsx",
  #'   vegetation_file = "data/ForestGEO_vegetation.xlsx",
  #'   output_dir = "results/statistics_only",
  #'   run_temporal_analysis = FALSE,
  #'   create_plots = FALSE
  #' )
  #' 
  #' # Accessing key results
  #' print(results_complete$comprehensive_comparison$stats)
  #' print(results_complete$site_similarity$richness_correlation)
  #' print(results_complete$vegetation_effects_katydids$correlations)
  #' 
  #' # Summary statistics
  #' cat("Total katydid species (bioacoustic):", 
  #'     ncol(results_complete$katydid_data$presence_matrix) - 1, "\n")
  #' cat("Total orthoptera species (metabarcoding):", 
  #'     length(results_complete$metabarcoding_data$species_list), "\n")
  #' cat("Species detected by both methods:", 
  #'     length(results_complete$comprehensive_comparison$katydid_comparison$both_methods), "\n")
  #' }
  

  
  # Initialization

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("     INTEGRATED BIOACOUSTIC vs METABARCODING ANALYSIS\n")
  cat("       ForestGEO BCI - Katydid Monitoring Project 2025-2026\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Create output directory structure
  if (save_results) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "figures/comparisons"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "figures/environmental"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "figures/temporal"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(output_dir, "models"), showWarnings = FALSE, recursive = TRUE)
  }
  
  # Validate required parameters
  if (missing(metabarcoding_file)) {
    stop("Error: metabarcoding_file parameter is required")
  }
  if (missing(vegetation_file)) {
    stop("Error: vegetation_file parameter is required")
  }
  if (!file.exists(metabarcoding_file)) {
    stop("Error: metabarcoding_file not found: ", metabarcoding_file)
  }
  if (!file.exists(vegetation_file)) {
    stop("Error: vegetation_file not found: ", vegetation_file)
  }
  
  # Log analysis parameters
  cat("Analysis parameters:\n")
  cat(sprintf("    - Bird confidence threshold (default): %.2f\n", confidence_threshold_birds))
  if (!is.null(bird_species_thresholds)) {
    cat(sprintf("    - Bird species-specific thresholds: %d species loaded\n", 
                nrow(bird_species_thresholds)))
  } else {
    cat("    - Bird species-specific thresholds: none (uniform)\n")
  }
  cat(sprintf("    - Katydid confidence threshold: %.2f\n", confidence_threshold_katydids))
  cat(sprintf("    - Process birds from scratch: %s\n", process_birds_from_scratch))
  cat(sprintf("    - Process katydids from scratch: %s\n", process_katydids_from_scratch))
  cat(sprintf("    - Run temporal analysis: %s\n", run_temporal_analysis))
  cat(sprintf("    - Run statistics: %s\n", run_statistics))
  cat(sprintf("    - Create plots: %s\n", create_plots))
  cat(sprintf("    - Output directory: %s\n\n", output_dir))
  
  # Check availability of analysis modules
  # On the 'main' branch, only core processing modules are present.
  # Statistical analysis functions are only available on 'analysis/*' branches.
  analysis_modules_available <- all(c(
    exists("comprehensive_method_comparison"),
    exists("analyze_site_similarity"),
    exists("analyze_vegetation_effects"),
    exists("analyze_multivariate_integrated")
  ))
  
  temporal_modules_available <- all(c(
    exists("calculate_temporal_rarefaction"),
    exists("calculate_rarefied_presence_matrix"),
    exists("analyze_temporal_patterns")
  ))
  
  plotting_modules_available <- all(c(
    exists("create_comparison_plots"),
    exists("save_all_plots")
  ))
  
  accumulation_modules_available <- exists("calculate_accumulation_abg") || 
    exists("calculate_species_accumulation")
  
  if (!analysis_modules_available && run_statistics) {
    cat("  [INFO] Statistical analysis modules not loaded.\n")
    cat("         Steps 5-8 (statistics, ordination, plots) will be skipped.\n")
    cat("         To enable: checkout the 'analysis/baseline' branch,\n")
    cat("         or source the analysis R files manually.\n\n")
    run_statistics <- FALSE
  }
  
  if (!temporal_modules_available && run_temporal_analysis) {
    cat("  [INFO] Temporal analysis modules not loaded - temporal analysis skipped.\n\n")
    run_temporal_analysis <- FALSE
  }
  
  if (!plotting_modules_available && create_plots) {
    cat("  [INFO] Plotting modules not loaded - plot creation skipped.\n\n")
    create_plots <- FALSE
  }
  
  

    # STEP 1: Load/process bird data

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 1: Processing bird data...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  bird_data <- tryCatch({
    
    if (process_birds_from_scratch && !is.null(bird_dir)) {
      cat("   Processing birds from BirdNET files (this may take a while...)\n")
      
      # Check if directory exists
      if (!dir.exists(bird_dir)) {
        stop("Bird directory not found: ", bird_dir)
      }
      
      # Identify available deployments
      deployments <- identify_deployments(bird_dir)
      
      if (nrow(deployments) == 0) {
        stop("No BirdNET deployments found in: ", bird_dir)
      }
      
      cat(sprintf("   Found %d deployments to process\n", nrow(deployments)))
      
      # Check for already processed deployments
      processed_deployments <- check_processed_deployments(output_dir)
      
      remaining_deployments <- deployments %>%
        filter(!deployment %in% processed_deployments)
      
      if (nrow(remaining_deployments) > 0 || !skip_processed_deployments) {
        
        deployments_to_process <- if (skip_processed_deployments) {
          remaining_deployments
        } else {
          deployments
        }
        
        cat(sprintf("   Processing %d deployments...\n", nrow(deployments_to_process)))
        
        # Process each deployment
        for (i in 1:nrow(deployments_to_process)) {
          deployment_info <- deployments_to_process[i, ]
          
          cat(sprintf("     [%d/%d] Processing: %s\n", 
                      i, nrow(deployments_to_process), deployment_info$deployment))
          
          detections <- process_single_deployment(
            deployment_path = deployment_info$full_path,
            deployment_name = deployment_info$deployment,
            confidence_threshold = confidence_threshold_birds,
            species_thresholds = bird_species_thresholds
          )
          
          # Save deployment results
          if (save_results) {
            save_deployment_results(detections, deployment_info$deployment, output_dir)
          }
        }
      } else {
        cat("   All deployments already processed, combining results...\n")
      }
      
      # Combine all deployments
      bird_data <- combine_all_deployments(output_dir)
      
    } else if (!is.null(bird_detections_file) && !is.null(bird_matrix_file)) {
      cat("   Loading pre-processed bird data files...\n")
      
      # Validate file existence
      if (!file.exists(bird_detections_file)) {
        stop("Bird detections file not found: ", bird_detections_file)
      }
      if (!file.exists(bird_matrix_file)) {
        stop("Bird matrix file not found: ", bird_matrix_file)
      }
      
      bird_detections <- read.csv(bird_detections_file, stringsAsFactors = FALSE)
      
      # Apply species-specific or uniform confidence filtering
      filtered_bird_detections <- apply_species_thresholds(
        bird_detections,
        species_thresholds = bird_species_thresholds,
        default_threshold = confidence_threshold_birds,
        species_col = "common_name"
      )
      
      # Recalculate presence matrix from filtered detections
      # The pre-saved matrix may contain species at lower confidence levels
      filtered_bird_matrix <- filtered_bird_detections %>%
        select(site, common_name) %>%
        distinct() %>%
        mutate(present = 1) %>%
        pivot_wider(names_from = common_name,
                    values_from = present,
                    values_fill = 0)
      
      cat(sprintf("   Applied confidence threshold %.2f: %d -> %d detections\n", 
                  confidence_threshold_birds, nrow(bird_detections), nrow(filtered_bird_detections)))
      cat(sprintf("   Species in recalculated matrix: %d (sites: %d)\n", 
                  ncol(filtered_bird_matrix) - 1, nrow(filtered_bird_matrix)))
      
      bird_data <- list(
        detections = filtered_bird_detections,
        presence_matrix = filtered_bird_matrix
      )
      
    } else {
      stop("Please provide either bird_dir for processing from scratch, ",
           "or both bird_detections_file and bird_matrix_file")
    }
    
    bird_data
    
  }, error = function(e) {
    cat("   ERROR in bird data processing:", conditionMessage(e), "\n")
    stop(e)
  })
  
  # Report bird data statistics
  n_bird_species <- ncol(bird_data$presence_matrix) - 1
  n_bird_detections <- nrow(bird_data$detections)
  cat(sprintf("   Bird data loaded: %d species, %d detections\n\n", 
              n_bird_species, n_bird_detections))
  
  # Create detailed bird species analysis
  if (process_birds_from_scratch && save_results) {
    cat("   Creating detailed bird species analysis...\n")
    bird_analysis <- create_bird_species_summary(bird_data$detections, output_dir)
    cat("   Bird species summary saved\n\n")
  }
  
  

    # STEP 2: Load/process katydid data

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 2: Loading katydid bioacoustic data...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  katydid_data <- tryCatch({
    
    load_katydid_data_improved(
      katydid_dir = katydid_dir,
      katydid_detection_file = katydid_detections_file,
      katydid_matrix_file = katydid_matrix_file,
      confidence_threshold = confidence_threshold_katydids,
      min_detection_days = min_detection_days,
      process_katydids_from_scratch = process_katydids_from_scratch,
      skip_processed_deployments = skip_processed_deployments,
      output_dir = output_dir
    )
    
  }, error = function(e) {
    cat("   ERROR in katydid data processing:", conditionMessage(e), "\n")
    stop(e)
  })
  
  # Report katydid data statistics
  n_katydid_species <- ncol(katydid_data$presence_matrix) - 1
  n_katydid_detections <- nrow(katydid_data$raw_detections)
  cat(sprintf("   Katydid data loaded: %d species, %d detections\n\n", 
              n_katydid_species, n_katydid_detections))
  
  # Create detailed katydid species analysis
  if (process_katydids_from_scratch && save_results) {
    cat("   Creating detailed katydid species analysis...\n")
    katydid_analysis <- create_katydid_species_summary(katydid_data$raw_detections, output_dir)
    cat("   Katydid species summary saved\n\n")
  }
  
  

  # STEP 3: Load Metabarcoding data

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 3: Loading DNA metabarcoding data...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  metabarcoding_data <- tryCatch({
    
    read_metabarcoding_data(metabarcoding_file)
    
  }, error = function(e) {
    cat("   ERROR in metabarcoding data loading:", conditionMessage(e), "\n")
    stop(e)
  })
  
  n_ortho_species <- length(metabarcoding_data$species_list)
  cat(sprintf("   Metabarcoding data loaded: %d orthoptera species\n\n", 
              n_ortho_species))
  
  

  # STEP 4: Load vegetation data

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 4: Loading vegetation data...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  vegetation_data <- tryCatch({
    
    read_vegetation_data(vegetation_file)
    
  }, error = function(e) {
    cat("   ERROR in vegetation data loading:", conditionMessage(e), "\n")
    stop(e)
  })
  
  n_plots <- nrow(vegetation_data$summary)
  cat(sprintf("   Vegetation data loaded for %d plots\n\n", n_plots))
  
  

  # Save Processed data (always, regardless of analysis mode)
  # This ensures clean matrices are available even in processing-only mode
  
  if (save_results) {
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat("Saving processed data...\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    # Bird presence matrix
    if (!is.null(bird_data$presence_matrix)) {
      write.csv(bird_data$presence_matrix,
                file.path(output_dir, "bird_presence_matrix.csv"),
                row.names = FALSE)
      cat("  [OK] bird_presence_matrix.csv\n")
    }
    
    # Bird detections
    if (!is.null(bird_data$detections)) {
      write.csv(bird_data$detections,
                file.path(output_dir, "bird_detections.csv"),
                row.names = FALSE)
      cat("  [OK] bird_detections.csv\n")
    }
    
    # Katydid presence matrix
    if (!is.null(katydid_data$presence_matrix)) {
      write.csv(katydid_data$presence_matrix,
                file.path(output_dir, "katydid_presence_matrix.csv"),
                row.names = FALSE)
      cat("  [OK] katydid_presence_matrix.csv\n")
    }
    
    # Katydid detections
    if (!is.null(katydid_data$raw_detections)) {
      write.csv(katydid_data$raw_detections,
                file.path(output_dir, "katydid_detections.csv"),
                row.names = FALSE)
      cat("  [OK] katydid_detections.csv\n")
    }
    
    # Metabarcoding data
    if (!is.null(metabarcoding_data$raw_data)) {
      write.csv(metabarcoding_data$raw_data,
                file.path(output_dir, "metabarcoding_orthoptera.csv"),
                row.names = FALSE)
      cat("  [OK] metabarcoding_orthoptera.csv\n")
    }
    if (!is.null(metabarcoding_data$presence_matrix)) {
      write.csv(metabarcoding_data$presence_matrix,
                file.path(output_dir, "metabarcoding_presence_matrix.csv"),
                row.names = FALSE)
      cat("  [OK] metabarcoding_presence_matrix.csv\n")
    }
    
    # Vegetation summary
    if (!is.null(vegetation_data$summary)) {
      write.csv(vegetation_data$summary,
                file.path(output_dir, "vegetation_summary.csv"),
                row.names = FALSE)
      cat("  [OK] vegetation_summary.csv\n")
    }
    
    cat("\n  All processed data saved to: ", output_dir, "\n")
    cat("  These files can be used as input for custom analyses.\n\n")
  }
  
  

  # STEP 5: Methodological comparisons

  comprehensive_results <- NULL
  site_similarity <- NULL
  temporal_patterns_katydids <- NULL
  katydid_rarefied <- NULL
  katydid_rarefied_matrix <- NULL
  bird_rarefied <- NULL
  bird_rarefied_matrix <- NULL
  katydid_accumulation <- NULL
  bird_accumulation <- NULL
  
  if (run_statistics) {
  
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 5: Comparing detection methods...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Katydids: Bioacoustic vs Metabarcoding
  cat("  5a. Katydids: Bioacoustic vs Metabarcoding...\n")
  comprehensive_results <- tryCatch({
    comprehensive_method_comparison(katydid_data, metabarcoding_data, output_dir)
  }, error = function(e) {
    cat("      Warning: Method comparison failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Site similarity analysis
  cat("  5b. Analyzing site-level similarity...\n")
  site_similarity <- tryCatch({
    analyze_site_similarity(katydid_data, metabarcoding_data)
  }, error = function(e) {
    cat("      Warning: Site similarity analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Temporal patterns (katydids)
  if (run_temporal_analysis) {
    cat("  5c. Analyzing temporal patterns...\n")
    temporal_patterns_katydids <- tryCatch({
      analyze_temporal_patterns(katydid_data)
    }, error = function(e) {
      cat("      Warning: Temporal analysis failed:", conditionMessage(e), "\n")
      NULL
    })
  } else {
    temporal_patterns_katydids <- NULL
    cat("  5c. Temporal analysis skipped (run_temporal_analysis = FALSE)\n")
  }
  

  # STEP 5d: Temporal rarefaction for realistic site richness
  # Using weekly windows instead of annual accumulation to get realistic
  # between-site variation in species richness
  
  cat("  5d. Calculating temporal rarefaction (weekly windows)...\n")
  
  # Katydid temporal rarefaction
  katydid_rarefied <- NULL
  katydid_rarefied_matrix <- NULL
  
  if (!is.null(katydid_data) && !is.null(katydid_data$raw_detections)) {
    katydid_rarefied <- tryCatch({
      calculate_temporal_rarefaction(
        detections_data = katydid_data$raw_detections,
        species_col = "common_name",
        confidence_threshold = confidence_threshold_katydids,
        window_days = 7,
        min_days_per_window = 5,
        n_bootstrap = 100
      )
    }, error = function(e) {
      cat("      Warning: Katydid rarefaction failed:", conditionMessage(e), "\n")
      NULL
    })
    
    # Rarefied presence matrix for ordinations
    katydid_rarefied_matrix <- tryCatch({
      calculate_rarefied_presence_matrix(
        detections_data = katydid_data$raw_detections,
        species_col = "common_name",
        confidence_threshold = confidence_threshold_katydids,
        window_days = 7,
        min_weeks_present = 2,
        min_days_detection = 5  # Laurel Symes criterion
      )
    }, error = function(e) {
      cat("      Warning: Katydid rarefied matrix failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(katydid_rarefied)) {
      cat(sprintf("      Katydid rarefied richness: %.1f - %.1f species per site\n",
                  min(katydid_rarefied$site_richness$mean_richness),
                  max(katydid_rarefied$site_richness$mean_richness)))
    }
  }
  
  # Bird temporal rarefaction
  bird_rarefied <- NULL
  bird_rarefied_matrix <- NULL
  
  if (!is.null(bird_data) && !is.null(bird_data$raw_detections)) {
    bird_rarefied <- tryCatch({
      calculate_temporal_rarefaction(
        detections_data = bird_data$raw_detections,
        species_col = "common_name",
        confidence_threshold = confidence_threshold_birds,
        window_days = 7,
        min_days_per_window = 5,
        n_bootstrap = 100
      )
    }, error = function(e) {
      cat("      Warning: Bird rarefaction failed:", conditionMessage(e), "\n")
      NULL
    })
    
    # Rarefied presence matrix
    bird_rarefied_matrix <- tryCatch({
      calculate_rarefied_presence_matrix(
        detections_data = bird_data$raw_detections,
        species_col = "common_name",
        confidence_threshold = confidence_threshold_birds,
        window_days = 7,
        min_weeks_present = 2,
        min_days_detection = 5
      )
    }, error = function(e) {
      cat("      Warning: Bird rarefied matrix failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(bird_rarefied)) {
      cat(sprintf("      Bird rarefied richness: %.1f - %.1f species per site\n",
                  min(bird_rarefied$site_richness$mean_richness),
                  max(bird_rarefied$site_richness$mean_richness)))
    }
  }
  
  cat("   Temporal rarefaction complete\n")
  

  # STEP 5e: Species accumulation curves (ALPHA-BETA-GAMMA)

  cat("  5e. Calculating species accumulation curves (alpha-beta-gamma)...\n")
  
  # Katydid accumulation with alpha-beta-gamma decomposition
  katydid_accumulation <- NULL
  if (!is.null(katydid_data) && !is.null(katydid_data$raw_detections)) {
    katydid_accumulation <- tryCatch({
      # Use new ABG function if available, fall back to original
      if (exists("calculate_accumulation_abg")) {
        calculate_accumulation_abg(
          detections_data = katydid_data$raw_detections,
          species_col = "common_name",
          confidence_threshold = confidence_threshold_katydids,
          n_permutations = 100,
          n_sites_composite = min(25, length(unique(katydid_data$raw_detections$site))),
          max_days = 16,
          seed = 42
        )
      } else {
        calculate_species_accumulation(
          detections_data = katydid_data$raw_detections,
          species_col = "common_name",
          confidence_threshold = confidence_threshold_katydids,
          n_permutations = 100,
          n_sites_composite = min(25, length(unique(katydid_data$raw_detections$site))),
          max_days = 7,
          seed = 42
        )
      }
    }, error = function(e) {
      cat("      Warning: Katydid accumulation failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(katydid_accumulation)) {
      cat(sprintf("      Katydid alpha-diversity: %.1f species (day 1)\n",
                  katydid_accumulation$alpha_diversity))
      if (!is.null(katydid_accumulation$beta_temporal)) {
        cat(sprintf("      Katydid temporal beta: %.2f species/day\n",
                    katydid_accumulation$beta_temporal))
      }
      cat(sprintf("      Katydid gamma contribution: %.1f species\n",
                  katydid_accumulation$gamma_diversity))
    }
  }
  
  # Bird accumulation with alpha-beta-gamma decomposition
  bird_accumulation <- NULL
  if (!is.null(bird_data) && !is.null(bird_data$raw_detections)) {
    bird_accumulation <- tryCatch({
      # Use new ABG function if available, fall back to original
      if (exists("calculate_accumulation_abg")) {
        calculate_accumulation_abg(
          detections_data = bird_data$raw_detections,
          species_col = "common_name",
          confidence_threshold = confidence_threshold_birds,
          n_permutations = 100,
          n_sites_composite = min(25, length(unique(bird_data$raw_detections$site))),
          max_days = 16,
          seed = 42
        )
      } else {
        calculate_species_accumulation(
          detections_data = bird_data$raw_detections,
          species_col = "common_name",
          confidence_threshold = confidence_threshold_birds,
          n_permutations = 100,
          n_sites_composite = min(25, length(unique(bird_data$raw_detections$site))),
          max_days = 7,
          seed = 42
        )
      }
    }, error = function(e) {
      cat("      Warning: Bird accumulation failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (!is.null(bird_accumulation)) {
      cat(sprintf("      Bird alpha-diversity: %.1f species (day 1)\n",
                  bird_accumulation$alpha_diversity))
      if (!is.null(bird_accumulation$beta_temporal)) {
        cat(sprintf("      Bird temporal beta: %.2f species/day\n",
                    bird_accumulation$beta_temporal))
      }
      cat(sprintf("      Bird gamma contribution: %.1f species\n",
                  bird_accumulation$gamma_diversity))
    }
  }
  
  # Save accumulation plots
  if (save_results && create_plots) {
    dir.create(file.path(output_dir, "figures/accumulation"), 
               showWarnings = FALSE, recursive = TRUE)
    
    if (!is.null(katydid_accumulation)) {
      # Use new ABG plot function if available
      if (exists("plot_accumulation_abg")) {
        katydid_accum_plot <- plot_accumulation_abg(
          katydid_accumulation,
          taxon_name = "katydid",
          show_annotations = TRUE,
          show_definitions = TRUE,
          colors = c(single = "#009E73", composite = "#000000"),
          output_file = file.path(output_dir, "figures/accumulation/katydid_accumulation_abg.jpg")
        )
      } else {
        katydid_accum_plot <- plot_species_accumulation(
          katydid_accumulation,
          taxon_name = "katydid",
          colors = c(single = "#009E73", composite = "#000000"),
          output_file = file.path(output_dir, "figures/accumulation/katydid_accumulation.jpg")
        )
      }
    }
    
    if (!is.null(bird_accumulation)) {
      # Use new ABG plot function if available
      if (exists("plot_accumulation_abg")) {
        bird_accum_plot <- plot_accumulation_abg(
          bird_accumulation,
          taxon_name = "bird",
          show_annotations = TRUE,
          show_definitions = TRUE,
          colors = c(single = "#0072B2", composite = "#000000"),
          output_file = file.path(output_dir, "figures/accumulation/bird_accumulation_abg.jpg")
        )
      } else {
        bird_accum_plot <- plot_species_accumulation(
          bird_accumulation,
          taxon_name = "bird",
          colors = c(single = "#0072B2", composite = "#000000"),
          output_file = file.path(output_dir, "figures/accumulation/bird_accumulation.jpg")
        )
      }
    }
  }
  
  cat("   Species accumulation complete\n")
  
  # Save methodological comparison results
  if (save_results && !is.null(comprehensive_results)) {
    write.csv(comprehensive_results$katydid_comparison$stats, 
              file.path(output_dir, "tables/katydid_method_comparison.csv"), 
              row.names = FALSE)
  }
  
  if (save_results && !is.null(site_similarity) && !is.null(site_similarity$site_richness)) {
    write.csv(site_similarity$site_richness, 
              file.path(output_dir, "tables/site_richness_comparison.csv"), 
              row.names = FALSE)
  }
  
  # Save temporal rarefaction results
  if (save_results) {
    # Katydid rarefied richness
    if (!is.null(katydid_rarefied)) {
      write.csv(katydid_rarefied$site_richness,
                file.path(output_dir, "tables/katydid_rarefied_richness.csv"),
                row.names = FALSE)
      write.csv(katydid_rarefied$weekly_data,
                file.path(output_dir, "tables/katydid_weekly_richness.csv"),
                row.names = FALSE)
      cat("      Saved: katydid_rarefied_richness.csv\n")
    }
    
    # Bird rarefied richness
    if (!is.null(bird_rarefied)) {
      write.csv(bird_rarefied$site_richness,
                file.path(output_dir, "tables/bird_rarefied_richness.csv"),
                row.names = FALSE)
      write.csv(bird_rarefied$weekly_data,
                file.path(output_dir, "tables/bird_weekly_richness.csv"),
                row.names = FALSE)
      cat("      Saved: bird_rarefied_richness.csv\n")
    }
    
    # Accumulation curve data
    if (!is.null(katydid_accumulation)) {
      accum_data <- data.frame(
        day = katydid_accumulation$single_recorder$day,
        single_mean = katydid_accumulation$single_recorder$mean_species,
        single_se = katydid_accumulation$single_recorder$se,
        composite_mean = katydid_accumulation$composite$mean_species[1:length(katydid_accumulation$single_recorder$day)],
        composite_se = katydid_accumulation$composite$se[1:length(katydid_accumulation$single_recorder$day)]
      )
      write.csv(accum_data,
                file.path(output_dir, "tables/katydid_accumulation_data.csv"),
                row.names = FALSE)
    }
    
    if (!is.null(bird_accumulation)) {
      accum_data <- data.frame(
        day = bird_accumulation$single_recorder$day,
        single_mean = bird_accumulation$single_recorder$mean_species,
        single_se = bird_accumulation$single_recorder$se,
        composite_mean = bird_accumulation$composite$mean_species[1:length(bird_accumulation$single_recorder$day)],
        composite_se = bird_accumulation$composite$se[1:length(bird_accumulation$single_recorder$day)]
      )
      write.csv(accum_data,
                file.path(output_dir, "tables/bird_accumulation_data.csv"),
                row.names = FALSE)
    }
    
    # Save diversity metrics summary (alpha-beta-gamma)
    diversity_summary <- data.frame(
      taxon = character(),
      alpha_diversity = numeric(),
      beta_temporal = numeric(),
      gamma_diversity = numeric(),
      total_species = integer(),
      n_sites = integer(),
      max_days = integer(),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(katydid_accumulation)) {
      diversity_summary <- rbind(diversity_summary, data.frame(
        taxon = "Katydid",
        alpha_diversity = katydid_accumulation$alpha_diversity,
        beta_temporal = ifelse(is.null(katydid_accumulation$beta_temporal), NA, 
                               katydid_accumulation$beta_temporal),
        gamma_diversity = katydid_accumulation$gamma_diversity,
        total_species = katydid_accumulation$total_species,
        n_sites = katydid_accumulation$n_sites_used,
        max_days = katydid_accumulation$max_days_used
      ))
    }
    
    if (!is.null(bird_accumulation)) {
      diversity_summary <- rbind(diversity_summary, data.frame(
        taxon = "Bird",
        alpha_diversity = bird_accumulation$alpha_diversity,
        beta_temporal = ifelse(is.null(bird_accumulation$beta_temporal), NA,
                               bird_accumulation$beta_temporal),
        gamma_diversity = bird_accumulation$gamma_diversity,
        total_species = bird_accumulation$total_species,
        n_sites = bird_accumulation$n_sites_used,
        max_days = bird_accumulation$max_days_used
      ))
    }
    
    if (nrow(diversity_summary) > 0) {
      write.csv(diversity_summary,
                file.path(output_dir, "tables/diversity_abg_summary.csv"),
                row.names = FALSE)
      cat("      Saved: diversity_abg_summary.csv\n")
    }
  }
  
  cat("   Method comparison complete\n\n")
  
  } else {
    # run_statistics == FALSE: skip Steps 5
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat("STEP 5: Statistical comparisons SKIPPED (analysis modules not loaded)\n")
    cat(paste(rep("-", 70), collapse = ""), "\n\n")
  }
  
  

  # STEP 6: Ecological relationships

  vegetation_effects_katydids <- NULL
  vegetation_effects_metabar <- NULL
  vegetation_effects_birds <- NULL
  katydids_birds_relationship <- NULL
  katydids_ortho_relationship <- NULL
  
  if (run_statistics) {
  
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 6: Analyzing ecological relationships...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Vegetation effects on different groups
  cat("  6a. Vegetation effects on katydids (bioacoustic)...\n")
  vegetation_effects_katydids <- tryCatch({
    analyze_vegetation_effects(katydid_data, vegetation_data, "katydid_bioacoustic")
  }, error = function(e) {
    cat("      Warning: Katydid vegetation analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  cat("  6b. Vegetation effects on orthoptera (metabarcoding)...\n")
  vegetation_effects_metabar <- tryCatch({
    analyze_vegetation_effects(metabarcoding_data, vegetation_data, "orthoptera_metabarcoding")
  }, error = function(e) {
    cat("      Warning: Metabarcoding vegetation analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  cat("  6c. Vegetation effects on birds...\n")
  vegetation_effects_birds <- tryCatch({
    analyze_vegetation_effects(bird_data, vegetation_data, "birds")
  }, error = function(e) {
    cat("      Warning: Bird vegetation analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Relationships between taxonomic groups
  cat("  6d. Relationships between taxonomic groups...\n")
  
  # Katydids vs Birds
  katydids_birds_relationship <- tryCatch({
    analyze_taxonomic_relationships(
      katydid_data, bird_data, vegetation_data,
      group_names = c("Katydids", "Birds")
    )
  }, error = function(e) {
    cat("      Warning: Katydid-bird relationship analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Katydids (bioacoustic) vs Orthoptera (metabarcoding)
  katydids_ortho_relationship <- tryCatch({
    analyze_taxonomic_relationships(
      katydid_data, metabarcoding_data, vegetation_data,
      group_names = c("Katydids Bioacoustic", "Orthoptera Metabarcoding")
    )
  }, error = function(e) {
    cat("      Warning: Katydid-orthoptera relationship analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Save ecological results
  if (save_results) {
    if (!is.null(vegetation_effects_katydids)) {
      write.csv(vegetation_effects_katydids$correlations, 
                file.path(output_dir, "tables/vegetation_correlations_katydids.csv"), 
                row.names = FALSE)
    }
    
    if (!is.null(vegetation_effects_metabar)) {
      write.csv(vegetation_effects_metabar$correlations, 
                file.path(output_dir, "tables/vegetation_correlations_metabarcoding.csv"), 
                row.names = FALSE)
    }
    
    if (!is.null(vegetation_effects_birds)) {
      write.csv(vegetation_effects_birds$correlations, 
                file.path(output_dir, "tables/vegetation_correlations_birds.csv"), 
                row.names = FALSE)
    }
    
    # Create combined correlation summary table (Appendix 4 format)
    correlation_summary <- data.frame()
    
    if (!is.null(vegetation_effects_katydids$correlations)) {
      katydid_corr <- vegetation_effects_katydids$correlations %>%
        mutate(Group = "Katydids", Method = "Bioacoustic")
      correlation_summary <- rbind(correlation_summary, katydid_corr)
    }
    
    if (!is.null(vegetation_effects_metabar$correlations)) {
      meta_corr <- vegetation_effects_metabar$correlations %>%
        mutate(Group = "Orthoptera", Method = "Metabarcoding")
      correlation_summary <- rbind(correlation_summary, meta_corr)
    }
    
    if (!is.null(vegetation_effects_birds$correlations)) {
      bird_corr <- vegetation_effects_birds$correlations %>%
        mutate(Group = "Birds", Method = "Bioacoustic")
      correlation_summary <- rbind(correlation_summary, bird_corr)
    }
    
    if (nrow(correlation_summary) > 0) {
      correlation_summary <- correlation_summary %>%
        select(Group, Method, variable, correlation, p_value) %>%
        arrange(Group, p_value) %>%
        mutate(
          correlation = round(correlation, 3),
          p_value = round(p_value, 4),
          significance = case_when(
            p_value < 0.001 ~ "***",
            p_value < 0.01 ~ "**",
            p_value < 0.05 ~ "*",
            p_value < 0.1 ~ ".",
            TRUE ~ "ns"
          )
        )
      
      write.csv(correlation_summary, 
                file.path(output_dir, "tables/TableA4_Vegetation_Correlations_Summary.csv"), 
                row.names = FALSE)
      cat("      Correlation summary saved (Appendix 4 format)\n")
    }
  }
  
  cat("   Ecological relationships analyzed\n\n")
  
  } else {
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat("STEP 6: Ecological relationships SKIPPED (analysis modules not loaded)\n")
    cat(paste(rep("-", 70), collapse = ""), "\n\n")
  }
  
  

  # STEP 7: Multivariate community analysis

  multivariate_katydids <- NULL
  multivariate_metabar <- NULL
  multivariate_birds <- NULL
  
  if (run_statistics) {
  
  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 7: Multivariate community analyses...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Multivariate analysis for each group
  cat("  7a. Katydid community ordination...\n")
  multivariate_katydids <- tryCatch({
    analyze_multivariate_integrated(
      katydid_data$presence_matrix, vegetation_data, "Katydid"
    )
  }, error = function(e) {
    cat("      Warning: Katydid multivariate analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  cat("  7b. Orthoptera (metabarcoding) community ordination...\n")
  multivariate_metabar <- tryCatch({
    analyze_multivariate_integrated(
      metabarcoding_data$presence_matrix, vegetation_data, "Orthoptera"
    )
  }, error = function(e) {
    cat("      Warning: Metabarcoding multivariate analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  cat("  7c. Bird community ordination...\n")
  multivariate_birds <- tryCatch({
    analyze_multivariate_integrated(
      bird_data$presence_matrix, vegetation_data, "Bird"
    )
  }, error = function(e) {
    cat("      Warning: Bird multivariate analysis failed:", conditionMessage(e), "\n")
    NULL
  })
  
  # Save multivariate results
  if (save_results) {
    if (!is.null(multivariate_katydids)) {
      capture.output(multivariate_katydids$env_fit, 
                     file = file.path(output_dir, "models/katydid_envfit_results.txt"))
      if (create_plots) {
        ggsave(file.path(output_dir, "models/katydid_community_ordination.jpg"), 
               multivariate_katydids$plot, width = 10, height = 7, dpi = 300)
      }
    }
    
    if (!is.null(multivariate_metabar)) {
      capture.output(multivariate_metabar$env_fit, 
                     file = file.path(output_dir, "models/metabarcoding_envfit_results.txt"))
      if (create_plots) {
        ggsave(file.path(output_dir, "models/metabarcoding_community_ordination.jpg"), 
               multivariate_metabar$plot, width = 10, height = 7, dpi = 300)
      }
    }
    
    if (!is.null(multivariate_birds)) {
      capture.output(multivariate_birds$env_fit, 
                     file = file.path(output_dir, "models/bird_envfit_results.txt"))
      if (create_plots) {
        ggsave(file.path(output_dir, "models/bird_community_ordination.jpg"), 
               multivariate_birds$plot, width = 10, height = 7, dpi = 300)
      }
    }
  }
  
  cat("   Multivariate analyses complete\n\n")
  
  } else {
    cat(paste(rep("-", 70), collapse = ""), "\n")
    cat("STEP 7: Multivariate analyses SKIPPED (analysis modules not loaded)\n")
    cat(paste(rep("-", 70), collapse = ""), "\n\n")
  }
  
  

  # STEP 8: Visualization and synthesis

  cat(paste(rep("-", 70), collapse = ""), "\n")
  cat("STEP 8: Creating comprehensive visualizations and synthesis...\n")
  cat(paste(rep("-", 70), collapse = ""), "\n")
  
  # Collect all results into a comprehensive list
  all_results <- list(
    bird_data = bird_data,
    katydid_data = katydid_data,
    metabarcoding_data = metabarcoding_data,
    vegetation_data = vegetation_data,
    comprehensive_comparison = comprehensive_results,
    method_comparison = if (!is.null(comprehensive_results)) {
      comprehensive_results$katydid_comparison
    } else NULL,
    site_similarity = site_similarity,
    temporal_patterns = temporal_patterns_katydids,
    vegetation_effects_katydids = vegetation_effects_katydids,
    vegetation_effects_metabar = vegetation_effects_metabar,
    vegetation_effects_birds = vegetation_effects_birds,
    katydids_birds_relationship = katydids_birds_relationship,
    katydids_ortho_relationship = katydids_ortho_relationship,
    multivariate_katydids = multivariate_katydids,
    multivariate_metabar = multivariate_metabar,
    multivariate_birds = multivariate_birds,
    # NEW: Temporal rarefaction results for realistic site richness
    katydid_rarefied = katydid_rarefied,
    katydid_rarefied_matrix = katydid_rarefied_matrix,
    bird_rarefied = bird_rarefied,
    bird_rarefied_matrix = bird_rarefied_matrix,
    # NEW: Species accumulation curves (alpha-beta-gamma)
    katydid_accumulation = katydid_accumulation,
    bird_accumulation = bird_accumulation
  )
  

  # STEP 8.1: Comprehensive GLM analysis

  glm_results <- NULL
  richness_ranges <- NULL
  
  if (run_statistics && exists("analyze_all_glm_relationships")) {
    cat("   Running comprehensive GLM analysis...\n")
    glm_results <- tryCatch({
      analyze_all_glm_relationships(all_results, output_dir)
    }, error = function(e) {
      cat("      Warning: GLM analysis failed:", conditionMessage(e), "\n")
      NULL
    })
  
  # Add GLM results to all_results
  all_results$glm_analysis <- glm_results
  
  if (!is.null(glm_results) && length(glm_results) > 0) {
    cat(sprintf("      Found %d significant GLM relationships\n", length(glm_results)))
  } else {
    cat("      No significant GLM relationships found (p < 0.05)\n")
  }
  
  # Calculate species richness ranges for the report
  cat("   Calculating species richness ranges...\n")
  richness_ranges <- tryCatch({
    calculate_species_richness_ranges(all_results)
  }, error = function(e) {
    cat("      Warning: Richness calculation failed:", conditionMessage(e), "\n")
    NULL
  })
  all_results$richness_ranges <- richness_ranges
  
  } # end if (run_statistics && exists("analyze_all_glm_relationships"))
  
  # Create comprehensive plots (requires plotting modules)
  if (create_plots && plotting_modules_available) {
    cat("   Creating comprehensive plots...\n")
    
    comparison_plots <- tryCatch({
      create_comparison_plots(all_results)
    }, error = function(e) {
      cat("      Warning: Plot creation failed:", conditionMessage(e), "\n")
      NULL
    })
    
    if (save_results && !is.null(comparison_plots)) {
      cat("   Saving all plots...\n")
      save_results_df <- save_all_plots(
        comparison_plots, 
        output_dir = output_dir  # save_all_plots() cree deja le sous-dossier "figures"
      )
      cat(sprintf("      Saved %d plots successfully\n", 
                  sum(save_results_df$success, na.rm = TRUE)))
    }
  } else {
    cat("   Plot creation skipped (create_plots = FALSE)\n")
  }
  
  # Save additional individual plots if available
  if (save_results && create_plots) {
    
    # Venn diagram
    if (!is.null(comprehensive_results) && 
        !is.null(comprehensive_results$katydid_comparison$venn_plot)) {
      png(file.path(output_dir, "figures/comparisons/species_detection_venn_diagram.png"), 
          width = 800, height = 600)
      grid::grid.draw(comprehensive_results$katydid_comparison$venn_plot)
      dev.off()
    }
    
    # Vegetation relationship plots
    if (!is.null(vegetation_effects_katydids) && 
        !is.null(vegetation_effects_katydids$richness_plot)) {
      ggsave(file.path(output_dir, "figures/environmental/katydid_vegetation_relationship.jpg"), 
             vegetation_effects_katydids$richness_plot, width = 10, height = 7, dpi = 300)
    }
    
    if (!is.null(vegetation_effects_metabar) && 
        !is.null(vegetation_effects_metabar$richness_plot)) {
      ggsave(file.path(output_dir, "figures/environmental/metabarcoding_vegetation_relationship.jpg"), 
             vegetation_effects_metabar$richness_plot, width = 10, height = 7, dpi = 300)
    }
    
    # Taxonomic relationship plots
    if (!is.null(katydids_birds_relationship) && 
        !is.null(katydids_birds_relationship$plot)) {
      ggsave(file.path(output_dir, "figures/comparisons/katydids_birds_relationship.jpg"), 
             katydids_birds_relationship$plot, width = 10, height = 7, dpi = 300)
    }
    
    if (!is.null(katydids_ortho_relationship) && 
        !is.null(katydids_ortho_relationship$plot)) {
      ggsave(file.path(output_dir, "figures/comparisons/katydids_orthoptera_relationship.jpg"), 
             katydids_ortho_relationship$plot, width = 10, height = 7, dpi = 300)
    }
  }
  
  cat("   Visualization and synthesis complete\n\n")
  
  

  # Executive summary

  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("            Analysis Complete - Executive Summary \n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Calculate key metrics
  total_katydid_species_bio <- n_katydid_species
  total_ortho_species_meta <- n_ortho_species
  total_bird_species <- n_bird_species
  
  # Species richness per site
  if (!is.null(richness_ranges)) {
    cat(" Species richness per sites:\n")
    
    if (!is.null(richness_ranges$katydid)) {
      cat(sprintf("     - Katydids: %s species per site (mean: %.1f)\n", 
                  richness_ranges$katydid$range_text, 
                  richness_ranges$katydid$mean))
    }
    
    if (!is.null(richness_ranges$bird)) {
      cat(sprintf("     - Birds: %s species per site (mean: %.1f)\n", 
                  richness_ranges$bird$range_text, 
                  richness_ranges$bird$mean))
    }
    
    if (!is.null(richness_ranges$metabarcoding)) {
      cat(sprintf("     - Metabarcoding: %s species per site (mean: %.1f)\n", 
                  richness_ranges$metabarcoding$range_text, 
                  richness_ranges$metabarcoding$mean))
    }
    cat("\n")
  }
  
  # Species counts
  cat(" Species counts:\n")
  cat(sprintf("     - Katydids (Bioacoustic): %d species\n", total_katydid_species_bio))
  cat(sprintf("     - Orthoptera (Metabarcoding): %d species\n", total_ortho_species_meta))
  cat(sprintf("     - Birds (for comparison): %d species\n\n", total_bird_species))
  
  # Detection overlap
  if (!is.null(comprehensive_results) && !is.null(comprehensive_results$katydid_comparison)) {
    overlap_species <- length(comprehensive_results$katydid_comparison$both_methods)
    only_bio <- length(comprehensive_results$katydid_comparison$only_acoustic)
    only_meta <- length(comprehensive_results$katydid_comparison$only_metabar)
    
    cat("Overlap: detection overlap:\n")
    cat(sprintf("     - Detected by both methods: %d species (%.1f%%)\n", 
                overlap_species, 
                ifelse(total_katydid_species_bio + total_ortho_species_meta > 0,
                       overlap_species / (total_katydid_species_bio + total_ortho_species_meta) * 100, 0)))
    cat(sprintf("     - Only bioacoustic: %d species (%.1f%%)\n", 
                only_bio,
                ifelse(total_katydid_species_bio + total_ortho_species_meta > 0,
                       only_bio / (total_katydid_species_bio + total_ortho_species_meta) * 100, 0)))
    cat(sprintf("     - Only metabarcoding: %d species (%.1f%%)\n\n", 
                only_meta,
                ifelse(total_katydid_species_bio + total_ortho_species_meta > 0,
                       only_meta / (total_katydid_species_bio + total_ortho_species_meta) * 100, 0)))
  }
  
  # Correlation results
  if (!is.null(site_similarity) && !is.null(site_similarity$richness_correlation)) {
    cat("Correlation: Site-levels correlations:\n")
    cat(sprintf("     - Spearman correlation between methods: r = %.3f, p = %.3f\n\n",
                site_similarity$richness_correlation$estimate,
                site_similarity$richness_correlation$p.value))
  }
  
  # Vegetation effects summary
  cat(" Vegetation effects:\n")
  
  if (!is.null(vegetation_effects_katydids) && !is.null(vegetation_effects_katydids$correlations)) {
    sig_correlations_katydids <- sum(vegetation_effects_katydids$correlations$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("     - Katydids: %d significant correlations with vegetation\n", sig_correlations_katydids))
  }
  
  if (!is.null(vegetation_effects_metabar) && !is.null(vegetation_effects_metabar$correlations)) {
    sig_correlations_meta <- sum(vegetation_effects_metabar$correlations$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("     - Metabarcoding: %d significant correlations with vegetation\n", sig_correlations_meta))
  }
  
  if (!is.null(vegetation_effects_birds) && !is.null(vegetation_effects_birds$correlations)) {
    sig_correlations_birds <- sum(vegetation_effects_birds$correlations$p_value < 0.05, na.rm = TRUE)
    cat(sprintf("     - Birds: %d significant correlations with vegetation\n\n", sig_correlations_birds))
  }
  
  # Temporal patterns
  if (!is.null(temporal_patterns_katydids) && !is.null(temporal_patterns_katydids$period_activity)) {
    most_active_period <- temporal_patterns_katydids$period_activity %>%
      filter(detection_count == max(detection_count)) %>%
      pull(period) %>%
      first()
    cat(" Temporal patterns:\n")
    cat(sprintf("     - Most active period for katydids: %s\n\n", most_active_period))
  }
  
  # Output directory information
  cat(sprintf(" ALL RESULTS SAVED IN: %s\n\n", output_dir))
  cat(" Key output files:\n")
  cat("     - tables/katydid_method_comparison.csv - Detection method comparison\n")
  cat("     - tables/site_richness_comparison.csv - Site-level richness data\n")
  cat("     - tables/vegetation_correlations_*.csv - Vegetation relationship analyses\n")
  cat("     - models/*_community_ordination.jpg - Multivariate community plots\n")
  cat("     - figures/comparisons/species_detection_venn_diagram.png - Method overlap\n")
  cat("     - figures/environmental/ - Vegetation-species relationship plots\n")
  cat("     - figures/temporal/ - Temporal activity pattern plots\n\n")
  
  # Main research question answer
  cat(" MAIN RESEARCH QUESTION ANSWER:\n")
  if (total_katydid_species_bio > 0 && total_ortho_species_meta > 0) {
    bio_proportion <- total_katydid_species_bio / (total_katydid_species_bio + total_ortho_species_meta) * 100
    cat("   'Which proportion of sound-producing species can be detected by bioacoustics?'\n")
    cat(sprintf("    %.1f%% of the total orthopteran assemblage was detected by bioacoustics\n", bio_proportion))
    
    if (!is.null(comprehensive_results) && !is.null(comprehensive_results$katydid_comparison)) {
      overlap <- length(comprehensive_results$katydid_comparison$both_methods)
      max_species <- max(total_katydid_species_bio, total_ortho_species_meta)
      cat(sprintf("    The two methods are %s complementary\n", 
                  ifelse(overlap / max_species < 0.5, "highly", "moderately")))
    }
  }
  
  cat("\n")
  cat(" Analysis pipeline completed successfully! \n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Return comprehensive results
  return(all_results)
}


# Resume interrupted analysis function

#' Resume an Interrupted Integrated Analysis
#'
#' Allows resuming an integrated analysis that was interrupted, by detecting
#' already processed deployments and continuing from where it left off.
#'
#' @description
#' This function checks for previously processed bird and katydid deployments
#' in the output directory, combines available results, and continues the 
#' analysis pipeline using pre-processed data. This is particularly useful for:
#' \itemize{
#'   \item Recovering from system crashes or interruptions
#'   \item Incrementally processing large datasets
#'   \item Re-running analyses with updated parameters
#' }
#'
#' @param output_dir Character string. Path to directory containing partially
#'   processed results from a previous run.
#' @param katydid_dir Character string. Path to raw katydid data directory.
#'   Only needed if katydid processing needs to continue from scratch.
#' @param katydid_detections_file Character string. Path to pre-processed
#'   katydid detections CSV (optional fallback).
#' @param katydid_matrix_file Character string. Path to pre-processed
#'   katydid presence matrix CSV (optional fallback).
#' @param metabarcoding_file Character string. Path to metabarcoding Excel file.
#'   Required parameter.
#' @param vegetation_file Character string. Path to vegetation data Excel file.
#'   Required parameter.
#'
#' @return A named list containing all analysis results (same structure as
#'   \code{\link{main_integrated_analysis}}), or NULL if no processed
#'   deployments were found.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks for processed bird deployments in output_dir
#'   \item Checks for processed katydid deployments in output_dir
#'   \item If processed data exists, combines results
#'   \item Calls \code{main_integrated_analysis} with pre-processed data
#'   \item Returns complete analysis results
#' }
#'
#' @note
#' The function assumes that individual deployment results are saved in the
#' standard format produced by \code{save_deployment_results} and
#' \code{save_katydid_deployment_results}.
#'
#' @seealso 
#' \code{\link{main_integrated_analysis}} for the main analysis pipeline
#' \code{\link{check_processed_deployments}} for deployment detection
#' \code{\link{combine_all_deployments}} for combining bird results
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Resume an interrupted analysis
#' results <- resume_integrated_analysis(
#'   output_dir = "integrated_results",
#'   metabarcoding_file = "data/ForLeon.xlsx",
#'   vegetation_file = "data/Vegetation_data_25_plots.xlsx"
#' )
#' 
#' # If no processed data found, returns NULL
#' if (is.null(results)) {
#'   message("No processed data found - starting from scratch")
#' }
#' }

resume_integrated_analysis <- function(output_dir, 
                                       katydid_dir = NULL,
                                       katydid_detections_file = NULL,
                                       katydid_matrix_file = NULL,
                                       metabarcoding_file, 
                                       vegetation_file) {
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("     Resuming integrated analysis\n")
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  # Check what's already been processed for birds
  processed_deployments <- check_processed_deployments(output_dir)
  
  # Check what's already been processed for katydids
  processed_katydid_deployments <- check_processed_katydid_deployments(output_dir)
  
  bird_data <- NULL
  katydid_data <- NULL
  
  # Check for processed bird data
  if (length(processed_deployments) > 0) {
    cat(sprintf("   Found %d processed bird deployments, combining results...\n", 
                length(processed_deployments)))
    
    tryCatch({
      bird_data <- combine_all_deployments(output_dir)
      cat(sprintf("   Bird data combined: %d species, %d detections\n", 
                  ncol(bird_data$presence_matrix) - 1, 
                  nrow(bird_data$detections)))
    }, error = function(e) {
      warning("Failed to combine bird deployments: ", conditionMessage(e))
    })
  } else {
    cat("   No processed bird deployments found\n")
  }
  
  # Check for processed katydid data
  if (length(processed_katydid_deployments) > 0) {
    cat(sprintf("   Found %d processed katydid deployments, combining results...\n", 
                length(processed_katydid_deployments)))
    
    tryCatch({
      katydid_data <- combine_all_katydid_deployments(output_dir)
      cat(sprintf("   Katydid data combined: %d species, %d detections\n", 
                  ncol(katydid_data$presence_matrix) - 1, 
                  nrow(katydid_data$detections)))
    }, error = function(e) {
      warning("Failed to combine katydid deployments: ", conditionMessage(e))
    })
  } else {
    cat("   No processed katydid deployments found\n")
  }
  
  # Continue with analysis if we have data
  if (!is.null(bird_data) || !is.null(katydid_data)) {
    cat("\n   Continuing with analysis using available data...\n\n")
    
    # Prepare file paths for pre-processed data
    bird_det_file <- if (!is.null(bird_data)) {
      file.path(output_dir, "combined_all_detections.csv")
    } else NULL
    
    bird_mat_file <- if (!is.null(bird_data)) {
      file.path(output_dir, "combined_presence_matrix.csv")
    } else NULL
    
    katydid_det_file <- if (!is.null(katydid_data)) {
      file.path(output_dir, "combined_all_katydid_detections.csv")
    } else katydid_detections_file
    
    katydid_mat_file <- if (!is.null(katydid_data)) {
      file.path(output_dir, "combined_katydid_presence_matrix.csv")
    } else katydid_matrix_file
    
    # Run the main analysis with pre-processed data
    results <- main_integrated_analysis(
      bird_detections_file = bird_det_file,
      bird_matrix_file = bird_mat_file,
      katydid_detections_file = katydid_det_file,
      katydid_matrix_file = katydid_mat_file,
      metabarcoding_file = metabarcoding_file,
      vegetation_file = vegetation_file,
      output_dir = output_dir,
      process_birds_from_scratch = FALSE,
      process_katydids_from_scratch = FALSE
    )
    
    return(results)
    
  } else {
    cat("\n   No processed deployments found - cannot resume analysis\n")
    cat("   Please run main_integrated_analysis with process_from_scratch = TRUE\n\n")
    return(NULL)
  }
}
