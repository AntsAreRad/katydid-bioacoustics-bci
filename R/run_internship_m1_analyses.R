# run_internship_m1_analyses.R
# Orchestrator for M1 IMABEE internship analyses
# Runs after the main pipeline (which produces CSVs in integrated_results/)
#
# Usage: source("R/run_internship_m1_analyses.R")
#   Expects `results` object from main_integrated_analysis() to be in scope,
#   or loads data from integrated_results/ if results is not available.
#
# Output: results/internship_m1/
#
# These analyses reproduce the figures and tables from:
#   Leon Brouille (2025) "Passive acoustic monitoring for katydid diversity
#   assessment" - M1 IMABEE internship report

cat("\n")
cat("  INTERNSHIP M1 ANALYSES\n")
cat("  Katydid Bioacoustics - BCI Panama\n")
cat("  Leon Brouille, M1 IMABEE, April-June 2025\n")
cat("\n")

# ---- Configuration ----
M1_OUTPUT_DIR <- file.path(CONFIG$output_dir, "..", "results", "internship_m1")
dir.create(M1_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- Source helper functions ----
# internship_helpers.R contains utility functions for the M1 analyses
helpers_path <- file.path(module_dir, "internship_helpers.R")
if (file.exists(helpers_path)) {
  source(helpers_path)
} else {
  stop("internship_helpers.R not found in ", module_dir,
       "\nMake sure all internship_*.R files are in R/")
}

# ---- Source internship analysis modules ----
m1_modules <- c(
  "internship_temporal_analysis.R",
  "internship_method_comparison.R",
  "internship_glm_analysis.R",
  "internship_multivariate.R",
  "internship_report_annexes.R"
)

for (mod in m1_modules) {
  mod_path <- file.path(module_dir, mod)
  if (file.exists(mod_path)) {
    tryCatch({
      source(mod_path)
      cat(sprintf("  [OK] Loaded: %s\n", mod))
    }, error = function(e) {
      cat(sprintf("  [!] Failed to load %s: %s\n", mod, e$message))
    })
  } else {
    cat(sprintf("  [!] Not found: %s\n", mod_path))
  }
}

# ---- Check that results object exists ----
if (!exists("results") || is.null(results)) {
  cat("\n  [!] No 'results' object found. Attempting to load from pipeline outputs...\n")

  # Try to load from integrated_results/ CSVs
  ir <- CONFIG$output_dir
  if (!dir.exists(ir)) {
    stop("Cannot find pipeline output directory: ", ir,
         "\nRun the main pipeline first (run_analysis.R)")
  }

  cat("  Loading from: ", ir, "\n")

  # Load bird data
  bird_det_file <- file.path(ir, "combined_all_detections.csv")
  bird_mat_file <- file.path(ir, "combined_presence_matrix.csv")
  if (file.exists(bird_det_file) && file.exists(bird_mat_file)) {
    bird_data <- list(
      detections      = read.csv(bird_det_file, stringsAsFactors = FALSE),
      presence_matrix = read.csv(bird_mat_file, stringsAsFactors = FALSE)
    )
  } else {
    stop("Bird data files not found in ", ir)
  }

  # Load katydid data
  kat_det_file <- file.path(ir, "combined_all_katydid_detections.csv")
  kat_mat_file <- file.path(ir, "combined_katydid_presence_matrix.csv")
  if (file.exists(kat_det_file) && file.exists(kat_mat_file)) {
    katydid_data <- list(
      raw_detections  = read.csv(kat_det_file, stringsAsFactors = FALSE),
      detections      = read.csv(kat_det_file, stringsAsFactors = FALSE),
      presence_matrix = read.csv(kat_mat_file, stringsAsFactors = FALSE)
    )
  } else {
    stop("Katydid data files not found in ", ir)
  }

  # Load metabarcoding
  metabarcoding_data <- tryCatch({
    source(file.path(module_dir, "data_processing.R"))
    read_metabarcoding_data(CONFIG$metabarcoding_file)
  }, error = function(e) {
    cat(sprintf("  [!] Could not load metabarcoding: %s\n", e$message))
    NULL
  })

  # Load vegetation
  vegetation_data <- tryCatch({
    read_vegetation_data(CONFIG$vegetation_file)
  }, error = function(e) {
    cat(sprintf("  [!] Could not load vegetation: %s\n", e$message))
    NULL
  })

  results <- list(
    bird_data          = bird_data,
    katydid_data       = katydid_data,
    metabarcoding_data = metabarcoding_data,
    vegetation_data    = vegetation_data,
    output_dir         = ir
  )

  cat("  [OK] Data loaded from pipeline outputs\n")
}

# ---- Run analyses ----
cat("\n  Starting M1 internship analyses...\n")

m1_start <- Sys.time()

# 1. Temporal analysis
temporal_results <- tryCatch({
  run_internship_temporal(results, M1_OUTPUT_DIR)
}, error = function(e) {
  cat(sprintf("  [!] Temporal analysis failed: %s\n", e$message))
  NULL
})

# 2. Method comparison
comparison_results <- tryCatch({
  run_internship_method_comparison(results, M1_OUTPUT_DIR)
}, error = function(e) {
  cat(sprintf("  [!] Method comparison failed: %s\n", e$message))
  NULL
})

# 3. GLM analysis
glm_results <- tryCatch({
  run_internship_glm(results, M1_OUTPUT_DIR)
}, error = function(e) {
  cat(sprintf("  [!] GLM analysis failed: %s\n", e$message))
  NULL
})

# 4. Multivariate (NMDS, CCA, GAM, variance partitioning)
multivariate_results <- tryCatch({
  run_internship_multivariate(results, M1_OUTPUT_DIR)
}, error = function(e) {
  cat(sprintf("  [!] Multivariate analysis failed: %s\n", e$message))
  NULL
})

# 5. Report annexes
annexes_results <- tryCatch({
  run_internship_report_annexes(
    results, M1_OUTPUT_DIR,
    glm_results        = glm_results,
    temporal_results   = temporal_results,
    comparison_results = comparison_results
  )
}, error = function(e) {
  cat(sprintf("  [!] Report annexes failed: %s\n", e$message))
  NULL
})

m1_end <- Sys.time()
m1_duration <- difftime(m1_end, m1_start, units = "mins")

# ---- Summary ----
cat("\n")
cat("  INTERNSHIP M1 ANALYSES COMPLETE\n")
cat(sprintf("  Duration: %.1f minutes\n", as.numeric(m1_duration)))
cat(sprintf("  Output: %s\n", M1_OUTPUT_DIR))
cat("\n")

output_files <- list.files(M1_OUTPUT_DIR, recursive = TRUE)
cat(sprintf("  Generated %d files:\n", length(output_files)))
for (f in head(output_files, 25)) {
  cat(sprintf("    - %s\n", f))
}
if (length(output_files) > 25) {
  cat(sprintf("    ... and %d more\n", length(output_files) - 25))
}

cat("\n")
cat("  DONE\n")
cat("\n")
