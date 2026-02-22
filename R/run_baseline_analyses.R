# run_baseline_analyses.R
# Runs all baseline paper analyses in sequence.
#
# USAGE:
#   Rscript run_baseline_analyses.R
#   or:  source("run_baseline_analyses.R")
#
# PREREQUISITES:
#   - Processing pipeline run first (integrated_results/ CSVs exist)
#   - Packages: tidyverse, ggplot2, lubridate, mgcv, vegan
#   - Optional: VennDiagram
#
# ANALYSES:
#   1. Accumulation curves (global ABG + monthly)
#   2. Species distribution (spatial, temporal, daily species count)
#   3. Detection count vs richness correlation
#   4. Co-occurrence birds x katydids
#   5. Mean first detection date per species
#   6. Venn diagram bioacoustic vs metabarcoding
#   7. Vegetation effects (Spearman, GLM, GAM)
#   8. Community structure (PCoA, envfit, inter-taxa)
#   9. Site composition similarity (Jaccard, beta decomposition)
#  10. Per-species call activity over time

cat("\n")
cat("  BASELINE PAPER ANALYSES\n")
cat("  Katydid Bioacoustics - BCI Panama\n")
cat("\n")

start_time <- Sys.time()

# check inputs
input_dir <- "integrated_results"
required <- c("katydid_presence_matrix.csv", "katydid_detections.csv",
              "bird_presence_matrix.csv", "bird_detections.csv",
              "metabarcoding_presence_matrix.csv", "vegetation_summary.csv")

cat("Checking input files...\n")
missing <- required[!file.exists(file.path(input_dir, required))]
if (length(missing) > 0) {
  for (f in missing) cat(sprintf("  missing: %s\n", f))
  stop("Run the processing pipeline first.")
}
cat("  All input files present.\n\n")

# find script directory
script_dir <- tryCatch(
  file.path(dirname(sys.frame(1)$ofile), "R"),
  error = function(e) "R"
)
if (!dir.exists(script_dir)) script_dir <- "R"
if (!dir.exists(script_dir)) script_dir <- "."

run_script <- function(script_name, label) {
  path <- file.path(script_dir, script_name)
  if (!file.exists(path)) path <- script_name
  if (!file.exists(path)) {
    cat(sprintf("[SKIP] %s not found\n", script_name))
    return(invisible(NULL))
  }
  cat(sprintf("\n[%s]\n", label))
  tryCatch({
    source(path, local = new.env(parent = globalenv()))
    cat(sprintf("[OK] %s\n", label))
  }, error = function(e) {
    cat(sprintf("[FAIL] %s: %s\n", label, e$message))
  })
}

# source shared helpers
helpers_path <- file.path(script_dir, "baseline_helpers.R")
if (file.exists(helpers_path)) {
  source(helpers_path)
  cat("  Loaded baseline_helpers.R\n\n")
} else {
  cat("  [!] baseline_helpers.R not found, scripts may fail\n\n")
}

run_script("baseline_accumulation_curves.R",   "1/10 Accumulation curves")
run_script("baseline_species_distribution.R",  "2/10 Species distribution")
run_script("baseline_detection_richness.R",    "3/10 Detection count vs richness")
run_script("baseline_cooccurrence.R",          "4/10 Co-occurrence birds x katydids")
run_script("baseline_first_detection.R",       "5/10 Mean first detection date")
run_script("baseline_venn_diagram.R",          "6/10 Venn diagram")
run_script("baseline_vegetation_effects.R",    "7/10 Vegetation effects")
run_script("baseline_community_structure.R",   "8/10 Community structure")
run_script("baseline_site_composition.R",     "9/10 Site composition similarity")
run_script("baseline_species_activity.R",    "10/10 Per-species call activity")

end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "mins")

cat("\n")
cat("  ALL BASELINE ANALYSES COMPLETE\n")
cat(sprintf("  Duration: %.1f minutes\n", as.numeric(duration)))
cat("\n")

output_files <- list.files("results", recursive = TRUE)
cat(sprintf("  Total output files: %d\n", length(output_files)))
dirs <- unique(dirname(output_files))
for (d in sort(dirs)) {
  n <- sum(dirname(output_files) == d)
  cat(sprintf("    results/%s/ (%d files)\n", d, n))
}
cat("\n")
