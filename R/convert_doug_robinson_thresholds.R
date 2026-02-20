# CONVERT DOUG ROBINSON BIRDNET THRESHOLDS
#
# Purpose: Convert Doug Robinson's species-by-species BirdNET evaluation CSV
#          into the standardized format expected by the analysis pipeline
#          (load_species_thresholds() in data_processing.R).
#
# Input:  "data/raw/1 - Species key Doug Robinson.csv" (raw evaluation)
# Output: "data/bird_species_thresholds_doug_robinson.csv" (pipeline-ready)
#
# Author: Leon Brouille (M2 IMABEE)
# Expert evaluator: Douglas Robinson
# Date: 2026-02
#
# USAGE:
#   source("convert_doug_robinson_thresholds.R")
#   # Or run from command line:
#   Rscript convert_doug_robinson_thresholds.R
#
# WHAT THIS SCRIPT DOES:
#   Doug evaluated ~469 species detected by BirdNET on BCI. For each species,
#   he assessed whether detections were correct or false positives, and for
#   reliable species, he provided confidence thresholds at 99%, 95%, and 90%
#   precision levels.
#
#   Classification rules applied:
#   - Positive threshold value (e.g. 0.25, 0.84) -> valid species, use threshold
#   - Negative threshold (e.g. -17.7, -82)       -> species misID'd, EXCLUDE
#   - "all correct" / "all positive" in notes     -> reliable, threshold = 0
#   - "all negative" / similar in notes            -> false positives, EXCLUDE
#   - Ambiguous notes without threshold            -> cautious EXCLUDE
#   - No notes, no threshold                       -> unevaluated (not in output)
#
#   The output file has the format expected by load_species_thresholds():
#     species, confidence_threshold, exclude, notes, threshold_99, threshold_95,
#     threshold_90, classification
#
#   Species with exclude=TRUE are removed by the updated load_species_thresholds().
#   Species NOT in the output file use CONFIG$bird_confidence (default 0.9).
#
# THRESHOLD LEVEL:
#   The script stores all threshold levels. Which level is used at runtime
#   is controlled by CONFIG$bird_threshold_level in run_analysis.R (default: "99").
#   Fallback: if chosen level is empty, fall back to a stricter level.
#
# ==============================================================================

library(tidyverse)

# Config - edit paths here if needed

INPUT_FILE  <- "/home/leon/Documents/Katydid/Katydid_fellowship/katydid-bioacoustics-bci/data/1-SpecieskeyDougRobinson.csv"
OUTPUT_FILE <- "data/bird_species_thresholds_doug_robinson.csv"


# Load and parse

cat("=== Converting Doug Robinson BirdNET Thresholds ===\n\n")

if (!file.exists(INPUT_FILE)) {
  stop("Input file not found: ", INPUT_FILE, 
       "\nPlace Doug's CSV in data/raw/ and re-run.")
}

# Read the CSV (UTF-8 with BOM)
raw <- read_csv(INPUT_FILE, show_col_types = FALSE, 
                locale = locale(encoding = "UTF-8"))

cat(sprintf("Loaded %d species from Doug Robinson's evaluation\n", nrow(raw)))

# Standardize column names
raw <- raw %>%
  rename(
    species       = Species,
    complete      = Complete,
    threshold_99  = `Threshold 99`,
    notes         = Notes,
    threshold_95  = `threshold 95`,
    threshold_90  = `threshold 90`,
    threshold_85  = `threshold 85`
  )

# Safely parse numeric columns
# Some threshold columns contain text instead of numbers
# (e.g. threshold_95 for "Black Hawk-Eagle" is "attila")

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

raw <- raw %>%
  mutate(
    threshold_99 = safe_numeric(threshold_99),
    threshold_95 = safe_numeric(threshold_95),
    threshold_90 = safe_numeric(threshold_90),
    threshold_85 = safe_numeric(threshold_85)
  )


# Classify each species

# Notes indicating all detections are false positives
negative_note_patterns <- c(
  "all negative", "all neg", "\\bnope\\b", "\\bnone\\b", "non-bird",
  "probably all negative", "all but one negative",
  "mostly negative", "all incorrect"
)

# Notes indicating species is reliably detected (no threshold needed)
positive_note_patterns <- c(
  "all correct", "all positive", "typically accurate"
)

# Build a single regex for each group
negative_regex <- paste(negative_note_patterns, collapse = "|")
positive_regex <- paste(positive_note_patterns, collapse = "|")

classified <- raw %>%
  mutate(
    notes_lower = tolower(replace_na(notes, "")),
    
    # --- classification logic ---
    classification = case_when(
      # 1. Has a positive numeric threshold at the 99% level
      !is.na(threshold_99) & threshold_99 > 0  ~ "valid",
      
      # 2. Negative threshold -> always exclude
      !is.na(threshold_99) & threshold_99 <= 0  ~ "exclude_negative_threshold",
      
      # 3. No t99, but has positive t95 or t90 -> valid
      is.na(threshold_99) & !is.na(threshold_95) & threshold_95 > 0 ~ "valid",
      is.na(threshold_99) & !is.na(threshold_90) & threshold_90 > 0 ~ "valid",
      
      # 4. Notes say all detections are false
      str_detect(notes_lower, negative_regex) ~ "exclude_all_negative",
      
      # 5. Notes say species is reliably detected
      str_detect(notes_lower, positive_regex)  ~ "valid_all_correct",
      
      # 6. Has notes but ambiguous -> cautious exclude
      nchar(trimws(notes_lower)) > 0           ~ "exclude_ambiguous",
      
      # 7. No threshold, no notes -> unevaluated
      TRUE ~ "unevaluated"
    )
  ) %>%
  select(-notes_lower)


# Build output

# For each threshold level, compute the "effective" value:
#   - positive number -> use it
#   - zero (from "all correct") -> 0
#   - excluded species -> NA (they'll be removed by the exclude column)
#   - negative number -> NA

make_effective <- function(raw_val, classification) {
  case_when(
    classification %in% c("exclude_negative_threshold", 
                           "exclude_all_negative", 
                           "exclude_ambiguous")    ~ NA_real_,
    classification == "valid_all_correct"           ~ 0,
    classification == "unevaluated"                 ~ NA_real_,
    !is.na(raw_val) & raw_val > 0                   ~ raw_val,
    TRUE                                            ~ NA_real_
  )
}

output <- classified %>%
  mutate(
    effective_t99 = make_effective(threshold_99, classification),
    effective_t95 = make_effective(threshold_95, classification),
    effective_t90 = make_effective(threshold_90, classification),
    
    # The "confidence_threshold" column used by load_species_thresholds()
    # will be filled at runtime based on CONFIG$bird_threshold_level.
    # For backwards compatibility, default to threshold_99.
    confidence_threshold = effective_t99,
    
    # Exclude flag
    exclude = classification %in% c("exclude_negative_threshold",
                                     "exclude_all_negative",
                                     "exclude_ambiguous")
  ) %>%
  # Remove unevaluated species (they'll get the default threshold)
  filter(classification != "unevaluated") %>%
  # Select and order output columns
  select(
    species         = species,
    confidence_threshold,
    exclude,
    threshold_99    = effective_t99,
    threshold_95    = effective_t95,
    threshold_90    = effective_t90,
    classification,
    notes           = notes
  ) %>%
  arrange(species)


# Save

dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
write_csv(output, OUTPUT_FILE)


# Summary report

cat("\n=== Classification summary ===\n\n")

summary_tbl <- classified %>%
  count(classification, name = "n_species") %>%
  arrange(desc(n_species))

for (i in 1:nrow(summary_tbl)) {
  cat(sprintf("  %-35s : %3d species\n", 
              summary_tbl$classification[i], summary_tbl$n_species[i]))
}
cat(sprintf("\n  Total evaluated: %d species\n", nrow(classified)))

n_in_output <- nrow(output)
n_excluded  <- sum(output$exclude)
n_valid     <- sum(!output$exclude)

cat(sprintf("\n=== Output files ===\n"))
cat(sprintf("  Species in output: %d\n", n_in_output))
cat(sprintf("    - With custom threshold: %d\n", n_valid))
cat(sprintf("    - Excluded (false positives): %d\n", n_excluded))
cat(sprintf("  Species NOT in output (will use default): %d\n",
            nrow(classified) - n_in_output))

cat(sprintf("\n=== Thresholds availability ===\n"))
cat(sprintf("  threshold_99 available: %d species\n",
            sum(!is.na(output$threshold_99) & !output$exclude)))
cat(sprintf("  threshold_95 available: %d species\n",
            sum(!is.na(output$threshold_95) & !output$exclude)))
cat(sprintf("  threshold_90 available: %d species\n",
            sum(!is.na(output$threshold_90) & !output$exclude)))

cat(sprintf("\nSaved to: %s\n", OUTPUT_FILE))
cat("Ready for the analysis pipeline.\n")
cat("To activate: set CONFIG$bird_species_thresholds_file in run_analysis.R\n")
