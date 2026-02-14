
# RUN ANALYSIS - Katydid Bioacoustics Pipeline

#
# Project: Comparing Bioacoustics vs DNA Metabarcoding for Katydid Monitoring
# Site: Barro Colorado Island (BCI), Panama
# Author: Leon Brouille (M1 IMABEE)
# Supervisors: Dr. Yves Basset (STRI), Dr. Greg Lamarre (STRI),
#              Dr. Laurel Symes (Cornell Lab of Ornithology)
# Date: 2025/2026
#
# USAGE:
#   Rscript run_analysis.R
#
#   Or in R:
#   source("run_analysis.R")
#
# REQUIREMENTS:
#   - R >= 4.0.0
#   - All module files in same directory or specified paths
#   - Data files accessible at specified paths
#


cat("
  KATYDID BIOACOUSTICS ANALYSIS PIPELINE
  Barro Colorado Island, Panama - STRI
")



# CONFIGURATION

# Analysis parameters
CONFIG <- list(
  # Confidence thresholds
  bird_confidence = 0.9,          # BirdNET default confidence threshold
  katydid_confidence = 0.95,       # Koogu/Katydid detector threshold
  min_detection_days = 5,         # Minimum days for species to be considered present
  
  # Species-specific thresholds (optional)
  # Set to NULL to use uniform threshold, or provide path to CSV file like data/examples/bird_species_thresholds_example.csv
  # CSV format: species,species_code,confidence_threshold,notes
  # Generate a template with: generate_threshold_template(detections, "my_thresholds.csv")
  bird_species_thresholds_file = NULL,
  
  # Batch processing
  max_files_per_batch = 3000,     # Files per batch to avoid memory issues
  
  # Processing options
  process_birds_from_scratch = TRUE,     # TRUE to process raw BirdNET files
  process_katydids_from_scratch = TRUE,  # TRUE to process raw katydid files
  skip_processed_deployments = TRUE,       # Skip already processed deployments
  
  # Analysis options
  # Legacy flags (kept for compatibility with main_analysis.R, always FALSE)
  # Baseline scripts in R/baseline_*.R replace these old modules
  run_temporal_analysis = FALSE,
  run_statistics = FALSE,
  create_plots = FALSE,
  enrich_with_bold = FALSE,

  run_baseline = TRUE,            # Run baseline paper analyses after processing
  
  # Output
  output_dir = "integrated_results",
  
  # File paths to be modified for your system

  # Option A: Process from scratch (set process_*_from_scratch = TRUE)
  bird_dir = "data/BirdNET_output",
  katydid_dir = "data/Koogu_output",
  
  # Option B: Use pre-processed files (set process_*_from_scratch = FALSE)
  bird_detections_file = "data/preprocessed/combined_all_detections.csv",
  bird_matrix_file = "data/preprocessed/combined_presence_matrix.csv",
  katydid_detections_file = "data/preprocessed/combined_all_katydid_detections.csv",
  katydid_matrix_file = "data/preprocessed/combined_katydid_presence_matrix.csv",
  
  # External data files (required)
  metabarcoding_file = "data/ForLeon.xlsx",
  vegetation_file = "data/Vegetation_data_25_plots.xlsx"
)



# PREREQUISITES CHECK

cat("\n[1/6] Checking prerequisites...\n")

# Check R version
if (getRversion() < "4.0.0") {
  stop("R version 4.0.0 or higher required. Current version: ", getRversion())
}
cat("  [OK] R version:", as.character(getRversion()), "\n")

# Required packages - Core (always needed)
core_packages <- c(
  "tidyverse",   # Data manipulation
  "readxl",      # Excel files
  "lubridate",   # Date handling
  "stringr",     # String operations
  "httr",        # HTTP requests (for BOLD)
  "xml2"         # XML parsing (for BOLD)
)

# Required packages - Analysis (only for analysis/* branches)
analysis_packages <- c(
  "vegan",       # Community ecology
  "ggplot2",     # Visualization
  "betapart",    # Beta diversity
  "ade4",        # Multivariate analysis
  "mgcv",        # GAMs
  "VennDiagram", # Venn diagrams
  "gridExtra",   # Plot arrangement
  "corrplot",    # Correlation plots
  "viridis"      # Color scales
)

# Check and install core packages (required)
missing_core <- core_packages[!sapply(core_packages, requireNamespace, quietly = TRUE)]
if (length(missing_core) > 0) {
  cat("  [WARNING] Missing core packages:", paste(missing_core, collapse = ", "), "\n")
  cat("  Installing missing packages...\n")
  install.packages(missing_core, quiet = TRUE)
}

# Check analysis packages (optional)
missing_analysis <- analysis_packages[!sapply(analysis_packages, requireNamespace, quietly = TRUE)]
if (length(missing_analysis) > 0) {
  cat("  [INFO] Missing analysis packages:", paste(missing_analysis, collapse = ", "), "\n")
  cat("         These are needed for statistical analyses (analysis/* branches)\n")
  cat("         Install with: install.packages(c(\"", 
      paste(missing_analysis, collapse = "\", \""), "\"))\n")
}

# Load core packages
suppressPackageStartupMessages({
  for (pkg in core_packages) {
    library(pkg, character.only = TRUE)
  }
})
cat("  [OK] Core packages loaded\n")

# Load analysis packages if available
analysis_pkgs_loaded <- 0
suppressPackageStartupMessages({
  for (pkg in analysis_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      library(pkg, character.only = TRUE)
      analysis_pkgs_loaded <- analysis_pkgs_loaded + 1
    }
  }
})
if (analysis_pkgs_loaded > 0) {
  cat(sprintf("  [OK] Analysis packages loaded (%d/%d)\n", 
              analysis_pkgs_loaded, length(analysis_packages)))
} else {
  cat("  [INFO] No analysis packages loaded (processing-only mode)\n")
}



# LOAD MODULE FILES

cat("\n[2/6] Loading module files...\n")

# Define module paths (modify if modules are in different location)
# Default: modules in ./R/ subdirectory relative to this script
module_dir <- file.path(dirname(sys.frame(1)$ofile %||% "."), "R")
if (!dir.exists(module_dir)) {
  # Fallback: try R/ in current working directory
  module_dir <- "R"
}
if (!dir.exists(module_dir)) {
  # Fallback: try same directory as this script
  module_dir <- "."
}

# Core modules (required - pipeline will not work without these)
core_modules <- c(
  "helper_functions.R",
  "data_processing.R",
  "main_analysis.R"
)

# Analysis modules (optional - only on analysis/* branches)
analysis_modules <- c(
  "temporal_analysis.R",
  "statistical_analysis.R",
  "plotting_functions.R",
  "species_accumulation_abg.R"
)

# Source core modules (required)
for (module in core_modules) {
  module_path <- file.path(module_dir, module)
  
  if (file.exists(module_path)) {
    tryCatch({
      source(module_path)
      cat(sprintf("  [OK] Loaded: %s\n", module))
    }, error = function(e) {
      cat(sprintf("  [X] Failed to load %s: %s\n", module, e$message))
      stop("Core module loading failed")
    })
  } else {
    cat(sprintf("  [X] File not found: %s\n", module_path))
    stop("Required core module file not found")
  }
}

# Source analysis modules (optional - graceful failure)
analysis_loaded <- c()
for (module in analysis_modules) {
  module_path <- file.path(module_dir, module)
  
  if (file.exists(module_path)) {
    tryCatch({
      source(module_path)
      cat(sprintf("  [OK] Loaded: %s\n", module))
      analysis_loaded <- c(analysis_loaded, module)
    }, error = function(e) {
      cat(sprintf("  [!] Warning: Failed to load %s: %s\n", module, e$message))
      cat("      Analysis functions from this module will not be available\n")
    })
  } else {
    cat(sprintf("  [--] Not found: %s (optional analysis module)\n", module))
  }
}

if (length(analysis_loaded) == 0) {
  cat("\n  [INFO] No analysis modules loaded - running in PROCESSING-ONLY mode\n")
  cat("         Data will be processed and exported as clean matrices.\n")
  cat("         For full analyses, checkout the 'analysis/baseline' branch.\n")
  CONFIG$run_statistics <- FALSE
  CONFIG$run_temporal_analysis <- FALSE
  CONFIG$create_plots <- FALSE
} else {
  cat(sprintf("\n  [OK] Loaded %d/%d analysis modules\n", 
              length(analysis_loaded), length(analysis_modules)))
}



# VERIFY DATA FILES

cat("\n[3/6] Verifying data files...\n")

# Load species-specific thresholds if configured
bird_species_thresholds <- NULL
if (!is.null(CONFIG$bird_species_thresholds_file)) {
  bird_species_thresholds <- load_species_thresholds(CONFIG$bird_species_thresholds_file)
  if (!is.null(bird_species_thresholds)) {
    cat(sprintf("  [OK] Species-specific thresholds loaded for %d species\n", 
                nrow(bird_species_thresholds)))
  }
} else {
  cat(sprintf("  [INFO] Using uniform bird confidence threshold: %.2f\n", 
              CONFIG$bird_confidence))
}

# Check metabarcoding file
if (file.exists(CONFIG$metabarcoding_file)) {
  cat(sprintf("  [OK] Metabarcoding: %s\n", basename(CONFIG$metabarcoding_file)))
} else {
  cat(sprintf("  [X] Metabarcoding file not found: %s\n", CONFIG$metabarcoding_file))
  stop("Metabarcoding file required")
}

# Check vegetation file
if (file.exists(CONFIG$vegetation_file)) {
  cat(sprintf("  [OK] Vegetation: %s\n", basename(CONFIG$vegetation_file)))
} else {
  cat(sprintf("  [X] Vegetation file not found: %s\n", CONFIG$vegetation_file))
  stop("Vegetation file required")
}

# Check bird data
if (CONFIG$process_birds_from_scratch) {
  if (dir.exists(CONFIG$bird_dir)) {
    cat(sprintf("  [OK] Bird directory: %s\n", CONFIG$bird_dir))
  } else {
    cat(sprintf("  [WARNING] Bird directory not found: %s\n", CONFIG$bird_dir))
  }
} else {
  if (file.exists(CONFIG$bird_detections_file) && file.exists(CONFIG$bird_matrix_file)) {
    cat("  [OK] Pre-processed bird files found\n")
  } else {
    cat("  [WARNING] Pre-processed bird files not found\n")
  }
}

# Check katydid data
if (CONFIG$process_katydids_from_scratch) {
  if (dir.exists(CONFIG$katydid_dir)) {
    cat(sprintf("  [OK] Katydid directory: %s\n", CONFIG$katydid_dir))
  } else {
    cat(sprintf("  [WARNING] Katydid directory not found: %s\n", CONFIG$katydid_dir))
  }
} else {
  if (file.exists(CONFIG$katydid_detections_file) && file.exists(CONFIG$katydid_matrix_file)) {
    cat("  [OK] Pre-processed katydid files found\n")
  } else {
    cat("  [WARNING] Pre-processed katydid files not found\n")
  }
}

# Threshold info already reported above


# CREATE OUTPUT DIRECTORY

cat("\n[4/6] Preparing output directory...\n")

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  [OK] Output directory: %s\n", CONFIG$output_dir))


# RUN ANALYSIS

cat("\n[5/6] Running main analysis...\n")
cat(" Grab a coffee this may take a long time\n\n")

start_time <- Sys.time()

# Run the main analysis
results <- tryCatch({
  main_integrated_analysis(
    # Bird data options
    bird_dir = if (CONFIG$process_birds_from_scratch) CONFIG$bird_dir else NULL,
    bird_detections_file = if (!CONFIG$process_birds_from_scratch) CONFIG$bird_detections_file else NULL,
    bird_matrix_file = if (!CONFIG$process_birds_from_scratch) CONFIG$bird_matrix_file else NULL,
    
    # Katydid data options
    katydid_dir = if (CONFIG$process_katydids_from_scratch) CONFIG$katydid_dir else NULL,
    katydid_detections_file = if (!CONFIG$process_katydids_from_scratch) CONFIG$katydid_detections_file else NULL,
    katydid_matrix_file = if (!CONFIG$process_katydids_from_scratch) CONFIG$katydid_matrix_file else NULL,
    
    # External data
    metabarcoding_file = CONFIG$metabarcoding_file,
    vegetation_file = CONFIG$vegetation_file,
    
    # Processing options
    process_birds_from_scratch = CONFIG$process_birds_from_scratch,
    process_katydids_from_scratch = CONFIG$process_katydids_from_scratch,
    skip_processed_deployments = CONFIG$skip_processed_deployments,
    
    # Confidence thresholds
    confidence_threshold_birds = CONFIG$bird_confidence,
    confidence_threshold_katydids = CONFIG$katydid_confidence,
    bird_species_thresholds = bird_species_thresholds,
    min_detection_days = CONFIG$min_detection_days,
    
    # Analysis options
    run_temporal_analysis = CONFIG$run_temporal_analysis,
    run_statistics = CONFIG$run_statistics,
    create_plots = CONFIG$create_plots,
    
    # Output
    output_dir = CONFIG$output_dir
  )
}, error = function(e) {
  cat(sprintf("\n[X] Analysis failed :/ %s\n", e$message))
  cat("\nStack trace:\n")
  traceback()
  return(NULL)
})

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")


# SUMMARY

cat("\n[6/6] Analysis Summary\n")
cat("\n")

if (!is.null(results)) {
  cat("\n[OK] Hurray analysis completed successfully!\n\n")
  
  # Print summary statistics
  cat("Results overview:\n")
  cat(sprintf("  - Duration: %.1f minutes\n", as.numeric(duration)))
  cat(sprintf("  - Output directory: %s\n", CONFIG$output_dir))
  
  # Threshold info
  cat("\nConfidence thresholds:\n")
  if (!is.null(bird_species_thresholds)) {
    cat(sprintf("  - Birds: Species-specific (%d species from %s)\n", 
                nrow(bird_species_thresholds),
                basename(CONFIG$bird_species_thresholds_file)))
    cat(sprintf("           Default for unlisted species: %.2f\n", CONFIG$bird_confidence))
    cat(sprintf("           Range: %.2f - %.2f\n",
                min(bird_species_thresholds$confidence_threshold),
                max(bird_species_thresholds$confidence_threshold)))
  } else {
    cat(sprintf("  - Birds: Uniform threshold = %.2f\n", CONFIG$bird_confidence))
  }
  cat(sprintf("  - Katydids: Uniform threshold = %.2f\n", CONFIG$katydid_confidence))
  
  # Data summary
  cat("\nDataset summary:\n")
  if (!is.null(results$bird_data)) {
    cat(sprintf("  - Bird species: %d\n", 
                ncol(results$bird_data$presence_matrix) - 1))
  }
  
  if (!is.null(results$katydid_data)) {
    cat(sprintf("  - Katydid species: %d\n", 
                ncol(results$katydid_data$presence_matrix) - 1))
  }
  
  if (!is.null(results$metabarcoding_data)) {
    cat(sprintf("  - Metabarcoding taxa: %d\n", 
                results$metabarcoding_data$total_taxa_count))
  }
  
  # Method comparison
  if (!is.null(results$method_comparison)) {
    cat(sprintf("  - Species overlap: %d\n", 
                length(results$method_comparison$both_methods)))
  }
  
  # List output files
  cat("\nGenerated files:\n")
  output_files <- list.files(CONFIG$output_dir, recursive = TRUE)
  for (f in head(output_files, 20)) {
    cat(sprintf("  - %s\n", f))
  }
  if (length(output_files) > 20) {
    cat(sprintf("  ... and %d more files\n", length(output_files) - 20))
  }
  
  # Run quality check
  cat("\nQuick data quality check:\n")
  check_data_quality(results)
  
} else {
  cat("\n[X] Analysis failed :/\n")
  cat("Check error messages above for details\n")
}

cat("\n")
cat("  PROCESSING COMPLETE\n")
cat("\n")

# -- baseline paper analyses --
if (!is.null(results) && isTRUE(CONFIG$run_baseline)) {
  baseline_script <- file.path(module_dir, "run_baseline_analyses.R")
  if (file.exists(baseline_script)) {
    cat("\n[BASELINE] Running baseline paper analyses...\n\n")
    tryCatch({
      source(baseline_script)
    }, error = function(e) {
      cat(sprintf("\n[BASELINE] Failed: %s\n", e$message))
      cat("  You can run it manually: source(\"R/run_baseline_analyses.R\")\n")
    })
  } else {
    cat("\n[BASELINE] run_baseline_analyses.R not found in ", module_dir, "\n")
    cat("  Baseline analyses skipped. Place the baseline scripts in R/ and\n")
    cat("  run: source(\"R/run_baseline_analyses.R\")\n")
  }
} else if (is.null(results)) {
  cat("  Baseline analyses skipped (processing failed)\n")
} else {
  cat("  Baseline analyses skipped (run_baseline = FALSE)\n")
  cat("  Set CONFIG$run_baseline <- TRUE or run manually:\n")
  cat("  source(\"R/run_baseline_analyses.R\")\n")
}

cat("\n")
cat("  DONE\n")
cat("\n")

# Return results (useful when sourcing the script)
invisible(results)
