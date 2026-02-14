# baseline_venn_diagram.R
# Venn diagram: bioacoustic vs metabarcoding species overlap
#
# The metabarcoding presence matrix uses BOLD BIN URIs as column names,
# while bioacoustic uses latin species names. To compare, we use the raw
# metabarcoding data (metabarcoding_orthoptera.csv) which contains a
# "species" column with resolved latin names for some BINs.
#
# Input: integrated_results/ CSV files
# Output: results/venn/ CSV and figures

library(tidyverse)

venn_available <- requireNamespace("VennDiagram", quietly = TRUE)
if (venn_available) library(VennDiagram)

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/venn"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


#' Compare species lists between bioacoustic and metabarcoding
#'
#' @param acoustic_species Character vector of bioacoustic species names.
#' @param metabar_species Character vector of metabarcoding species names
#'   (resolved latin names, not BINs).
#' @return List with both_methods, only_acoustic, only_metabar, stats.
#' @export
compare_species_lists <- function(acoustic_species, metabar_species) {
  both <- intersect(acoustic_species, metabar_species)
  only_ac <- setdiff(acoustic_species, metabar_species)
  only_mb <- setdiff(metabar_species, acoustic_species)
  total <- length(union(acoustic_species, metabar_species))

  stats <- data.frame(
    category = c("bioacoustic_only", "metabarcoding_only", "both_methods",
                 "total_bioacoustic", "total_metabarcoding", "total_unique"),
    count = c(length(only_ac), length(only_mb), length(both),
              length(acoustic_species), length(metabar_species), total),
    percentage = round(c(length(only_ac), length(only_mb), length(both),
                         length(acoustic_species), length(metabar_species),
                         total) / total * 100, 1),
    stringsAsFactors = FALSE
  )

  list(both_methods = both, only_acoustic = only_ac,
       only_metabar = only_mb, stats = stats)
}


#' Create Venn diagram PNG
#'
#' @param comparison List from compare_species_lists.
#' @param output_file Character.
#' @return Invisible NULL.
#' @export
create_venn_diagram <- function(comparison, output_file) {
  if (!venn_available) {
    cat("  VennDiagram package not available, skipping\n")
    return(invisible(NULL))
  }
  venn <- VennDiagram::venn.diagram(
    x = list(
      Bioacoustic = c(comparison$both_methods, comparison$only_acoustic),
      Metabarcoding = c(comparison$both_methods, comparison$only_metabar)
    ),
    filename = NULL,
    fill = c("#56B4E9", "#E69F00"),
    alpha = 0.5,
    cex = 1.8,
    cat.cex = 1.3,
    cat.fontface = "bold",
    main = "Katydid species detection overlap",
    main.cex = 1.4
  )
  png(output_file, width = 800, height = 600)
  grid::grid.draw(venn)
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(output_file)))
  invisible(NULL)
}


# -- load data --
cat("Loading data for Venn comparison...\n")
katydid_matrix <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                           check.names = FALSE)
acoustic_species <- setdiff(colnames(katydid_matrix), "site")

# Load raw metabarcoding data to get resolved species names
metabar_raw_file <- file.path(INPUT_DIR, "metabarcoding_orthoptera.csv")
if (!file.exists(metabar_raw_file)) {
  cat("  metabarcoding_orthoptera.csv not found, trying presence matrix\n")
  metabar_matrix <- read.csv(file.path(INPUT_DIR, "metabarcoding_presence_matrix.csv"),
                             check.names = FALSE)
  metabar_species <- setdiff(colnames(metabar_matrix), "site")
} else {
  metabar_raw <- read.csv(metabar_raw_file, stringsAsFactors = FALSE)
  # Extract unique resolved species names (non-empty, non-NA)
  if ("species" %in% colnames(metabar_raw)) {
    metabar_species <- metabar_raw %>%
      filter(!is.na(species) & species != "") %>%
      pull(species) %>%
      unique()
    cat(sprintf("  Metabarcoding: %d BINs total, %d with resolved species names\n",
                n_distinct(metabar_raw$bin_uri), length(metabar_species)))
  } else {
    # fallback to BINs
    metabar_species <- unique(metabar_raw$bin_uri)
    cat("  No species column found, using BIN URIs\n")
  }
}

# -- comparison --
cat("Comparing species lists...\n")
comparison <- compare_species_lists(acoustic_species, metabar_species)

# species detail table
all_sp <- union(acoustic_species, metabar_species)
species_detail <- data.frame(species = all_sp, stringsAsFactors = FALSE) %>%
  mutate(
    bioacoustic = species %in% acoustic_species,
    metabarcoding = species %in% metabar_species,
    method = case_when(
      bioacoustic & metabarcoding ~ "both",
      bioacoustic ~ "bioacoustic_only",
      TRUE ~ "metabarcoding_only"
    )
  ) %>%
  arrange(method, species)

write.csv(comparison$stats, file.path(OUTPUT_DIR, "overlap_stats.csv"),
          row.names = FALSE)
write.csv(species_detail, file.path(OUTPUT_DIR, "species_by_method.csv"),
          row.names = FALSE)

# -- Venn diagram --
cat("Creating Venn diagram...\n")
create_venn_diagram(comparison,
                    file.path(FIGURE_DIR, "venn_bioacoustic_metabarcoding.png"))

# -- summary --
cat("\n-- Venn diagram summary --\n")
cat(sprintf("  Bioacoustic species: %d\n", length(acoustic_species)))
cat(sprintf("  Metabarcoding identified species: %d\n", length(metabar_species)))
cat(sprintf("  Both methods: %d\n", length(comparison$both_methods)))
if (length(comparison$both_methods) > 0)
  cat(sprintf("  Shared: %s\n", paste(comparison$both_methods, collapse = ", ")))
cat(sprintf("  Bioacoustic only: %d\n", length(comparison$only_acoustic)))
cat(sprintf("  Metabarcoding only: %d\n", length(comparison$only_metabar)))
cat("\nVenn diagram analysis complete.\n")
