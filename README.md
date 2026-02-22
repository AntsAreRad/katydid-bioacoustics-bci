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
| `analysis/internship-m1` | `main` + M1 internship analyses | Reproduce figures and tables from the M1 IMABEE internship report (temporal patterns, GLM, NMDS, method comparison, etc.) |

### Branch `main` — Processing Pipeline

```
run_analysis.R                           # Entry point — configure and run
convert_doug_robinson_thresholds.R       # Convert Doug Robinson's raw evaluation to pipeline format
R/
├── helper_functions.R                   # Deployment detection, BOLD API, utilities
├── data_processing.R                    # BirdNET + Koogu batch processing, threshold functions
└── main_analysis.R                      # Integration, matrices, metabarcoding
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
├── baseline_site_composition.R     # Jaccard similarity, beta diversity decomposition (Baselga 2010)
├── baseline_species_activity.R     # Per-species monthly call activity, heatmaps
└── run_baseline_analyses.R         # Runs all 10 analyses in sequence
```

### Branch `analysis/internship-m1` — M1 IMABEE Internship Analyses

Everything in `main`, plus refactored analysis modules and dedicated internship scripts:

```
R/
├── internship_helpers.R              # Utility functions (date parsing, richness, site names)
├── temporal_analysis.R               # Temporal pattern functions (hourly, diel, rarefaction)
├── statistical_analysis.R            # Statistical analysis functions (GLM, NMDS, CCA, coinertia)
├── plotting_functions.R              # Plot generation functions
├── species_accumulation_abg.R        # Alpha-beta-gamma accumulation curves
├── internship_temporal_analysis.R    # → Figure 2: hourly activity patterns, diel periods
├── internship_method_comparison.R    # → Figure 1: Venn diagram, coinertia, family analysis
├── internship_glm_analysis.R         # → Figure 3 + Appendix 3: GLM all combinations
├── internship_multivariate.R         # → Figure 4: NMDS detection counts, CCA, GAM, variance partitioning
├── internship_report_annexes.R       # → Appendix 1–5: summary tables
└── run_internship_m1_analyses.R      # Runs all 5 internship analyses in sequence
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

### 3. With M1 internship analyses (`analysis/internship-m1` branch)

```bash
git checkout analysis/internship-m1
```

```r
# Set CONFIG$run_internship_m1 <- TRUE in run_analysis.R, then:
source("run_analysis.R")

# Or run internship analyses separately after processing:
source("R/run_internship_m1_analyses.R")
```

Results go to `results/internship_m1/` with figures and tables matching the M1 report.

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
- **Birds**: BirdNET v2.4 with species-specific confidence thresholds (see below)
- **Timezone**: AudioMoth UTC → Panama local time (UTC-5)
- **Metabarcoding**: DNA-based orthopteran survey for comparison 

### Species-specific BirdNET thresholds

Douglas Robinson (BCI ornithologist) manually evaluated the reliability of BirdNET detections species by species across all 469 species detected on BCI. His evaluation classified each species into one of three categories:

- **Valid with custom threshold** (64 species): BirdNET detects the species reliably above a certain confidence score. Each species gets its own threshold (range: 0.25 to 0.98).
- **Excluded** (405 species): all detections are false positives (misidentified calls, insects, mechanical noise, etc.). These species are removed entirely from the dataset.
- **Not evaluated**: species detected by BirdNET but absent from Doug's file use the default uniform threshold (0.9).

Doug provided thresholds at three precision levels (99%, 95%, 90%). The pipeline uses 99% by default (most conservative) with automatic fallback to 95% or 90% when the stricter level is not available for a given species.

The thresholds are enabled by default. To configure:

```r
# In run_analysis.R:
CONFIG$bird_species_thresholds_file <- "data/bird_species_thresholds_doug_robinson.csv"
CONFIG$bird_threshold_level <- "99"   # "99", "95", or "90"

# To disable and revert to uniform threshold:
CONFIG$bird_species_thresholds_file <- NULL
```

To regenerate the converted thresholds file from Doug's raw evaluation:

```r
source("convert_doug_robinson_thresholds.R")
```

The raw evaluation is stored in `data/raw/1 - Species key Doug Robinson.csv`. See `data/examples/bird_species_thresholds_example.csv` for the simple format if you want to define your own thresholds without the conversion step.

---

## Known Issues

- **Site S21**: Year 2000 timestamps seem to have been produced for some deployments. These are filtered out in `internship_helpers.R` / `baseline_helpers.R` (`extract_local_date` sets dates < 2024 to NA). The raw detections remain valid.
- **Site S09**: Consistently low diversity across both taxa. Identified as an outlier in detection-richness analyses.
- **BirdNET "nocall" entries**: BirdNET outputs a "nocall" label when no species is detected. These are automatically filtered out during threshold application.

## License

MIT License
