# Katydid Bioacoustics — Barro Colorado Island, Panama

Pipeline for processing and analyzing bioacoustic and DNA metabarcoding data for katydid (Tettigoniidae) diversity assessment on Barro Colorado Island.

**Author:** Léon Brouillé (M2 IMABEE, Université de Rennes / VU Amsterdam)

**Supervisors:** Dr. Yves Basset (STRI), Dr. Greg Lamarre (STRI), Dr. Laurel Symes (Cornell Lab of Ornithology)

**Collaboration:** Smithsonian Tropical Research Institute (STRI) + Cornell Lab of Ornithology

---

## Repository Structure

This repository is organized into **branches**:

| Branch | Contents | Use case |
|--------|----------|----------|
| `main` | Data processing pipeline only | Process BirdNET/Koogu raw outputs into clean matrices for your own analyses |
| `analysis/baseline` | `main` + statistical analyses | Reproduce the baseline paper results (accumulation curves, community structure, seasonal patterns, etc.) |

### Branch `main` — Processing Pipeline

```
run_analysis.R               # Entry point — configure and run
R/
├── helper_functions.R        # Deployment detection, BOLD API, utilities
├── data_processing.R         # BirdNET + Koogu batch processing
└── main_analysis.R           # Integration, matrices, metabarcoding
data/
├── ForLeon.xlsx              # Metabarcoding reference data
├── Vegetation_data_25_plots.xlsx  # Vegetation plot data
└── examples/
    └── bird_species_thresholds_example.csv
```

### Branch `analysis/baseline` — Baseline Paper Analyses

Everything in `main`, plus:

```
R/
├── baseline_helpers.R              # Shared functions (date parsing, richness)
├── baseline_accumulation_curves.R  # Alpha-beta-gamma + monthly + seasonal
├── baseline_species_distribution.R # Daily richness, monthly boxplots
├── baseline_detection_richness.R   # Detection count vs richness (Spearman)
├── baseline_cooccurrence.R         # Bird × katydid co-occurrence
├── baseline_first_detection.R      # Relative first detection day (mobility)
├── baseline_venn_diagram.R         # Bioacoustic vs metabarcoding overlap
├── baseline_vegetation_effects.R   # Vegetation correlations (Spearman, GLM, GAM)
├── baseline_community_structure.R  # PCoA, envfit, inter-taxa correlation
└── run_baseline_analyses.R         # Runs all 8 analyses in sequence
```

---

## Quick Start

### 1. Processing only (`main` branch)

```r
# Edit run_analysis.R to set your paths, then:
source("run_analysis.R")
```

This produces clean CSVs in `integrated_results/`:
- `katydid_detections.csv`, `bird_detections.csv` — all detections
- `katydid_presence_matrix.csv`, `bird_presence_matrix.csv` — site × species matrices
- `metabarcoding_presence_matrix.csv`, `vegetation_summary.csv`

### 2. With baseline analyses (`analysis/baseline` branch)

```bash
git checkout analysis/baseline
```

```r
# Set CONFIG$run_baseline <- TRUE in run_analysis.R, then:
source("run_analysis.R")

# Or run baseline analyses separately after processing:
source("R/run_baseline_analyses.R")
```

Results go to `results/` with CSVs and figures in subdirectories.

---

## Adding Your Own Analyses

After processing, the clean outputs are ready for any custom analysis:

```r
# Load the processed data
katydid_matrix <- read.csv("integrated_results/katydid_presence_matrix.csv", check.names = FALSE)
bird_det <- read.csv("integrated_results/bird_detections.csv")

# Your custom analysis here
```

---

## Study Design

- **25 SwiftOne recorders** deployed across BCI forest plots
- **11 deployments** processed (February 2024 — January 2025)
- **Katydids**: Koogu detector (confidence ≥ 0.95), 5-day minimum detection criterion
- **Birds**: BirdNET (confidence ≥ 0.9), species-specific thresholds supported
- **Timezone**: AudioMoth UTC → Panama local time (UTC−5)
- **Metabarcoding**: DNA-based orthopteran survey for comparison

### Species-Specific BirdNET Thresholds

You can define per-species confidence thresholds:

1. Run the pipeline once with the uniform threshold
2. Generate a template: `generate_threshold_template(results$bird_data$detections, "data/bird_species_thresholds.csv")`
3. Edit the CSV, adjust thresholds per species
4. Set `CONFIG$bird_species_thresholds_file <- "data/bird_species_thresholds.csv"` in `run_analysis.R`

See `data/examples/bird_species_thresholds_example.csv` for format.

---

## Known Issues

- **Site S21**: AudioMoth RTC was not synchronised for some deployments, producing year-2000 timestamps. These are filtered out in `baseline_helpers.R` (`extract_local_date` sets dates < 2024 to NA). The raw detections remain valid.
- **Site S09**: Geographically isolated, consistently low diversity across both taxa. Identified as an outlier in detection–richness analyses.

---

## License

MIT License
