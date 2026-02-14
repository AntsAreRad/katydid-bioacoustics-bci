# ==============================================================================
# STATISTICAL ANALYSIS FUNCTIONS - COMPLETE
# ==============================================================================
#
# Functions for statistical analyses of katydid biodiversity data including
# diversity metrics, ordination, multivariate analyses, and vegetation effects
#
# Project: Katydid Bioacoustics - BCI 2025
# Context: Comparison of bioacoustic vs DNA metabarcoding methods for 
#          monitoring Orthoptera (katydids) in tropical rainforest
# 
# Location: Barro Colorado Island (BCI), Panama
# Collaboration: STRI (Yves Basset, Greg Lamarre) + Cornell Lab (Laurel Symes)
#
# Author: Leon Brouille (M1 IMABEE)
# Supervisors: Dr. Yves Basset, Dr. Greg Lamarre, Dr. Laurel Symes
# Date: 2025-10-18
# Last Updated: 2025-10-18 (MERGED ALL SECTIONS)
#
# Required packages:
#   - vegan (ordination, diversity metrics)
#   - betapart (beta diversity)
#   - ade4 (coinertia analysis)
#   - mgcv (GAM models)
#   - tidyverse
#   - ggplot2
#   - viridis (color scales)
#   - VennDiagram (method comparison)
#
# CONTENT:
#   SECTION 1: Site Similarity and Beta Diversity
#   SECTION 2: Vegetation Effects
#   SECTION 3: Variable Selection
#   SECTION 4: NMDS Ordination with Detection Counts
#   SECTION 5: CCA Ordination with Environmental Constraints
#   SECTION 6: Nonlinear Relationships (GAM Analysis)
#   SECTION 7: Comprehensive GLM Analysis
#   SECTION 8: Helper Functions for GLM Analysis
#   SECTION 9: Wrapper Functions for Multivariate Analyses
#   SECTION 10: Method Comparison and Additional Analyses
#
# ==============================================================================

# Required libraries
library(vegan)
library(betapart)
library(ade4)
library(mgcv)
library(tidyverse)
library(ggplot2)
library(viridis)
library(VennDiagram)

# ==============================================================================
# SECTION 1: SITE SIMILARITY AND BETA DIVERSITY
# ==============================================================================

#' Analyze site similarity between acoustic and metabarcoding data
#'
#' Performs comprehensive site similarity analyses comparing species composition
#' detected by bioacoustics versus DNA metabarcoding. Includes distance matrix
#' calculations, Mantel tests, coinertia analysis, and richness comparisons.
#'
#' @param acoustic_data List. Output from combine_all_deployments() or similar,
#'   containing presence_matrix component
#' @param metabarcoding_data List. Metabarcoding dataset with presence_matrix
#'
#' @return List with components:
#'   \itemize{
#'     \item site_richness: data frame with metabar_richness and acoustic_richness
#'       per site
#'     \item richness_correlation: Spearman correlation test result between
#'       richness values
#'     \item mantel_test: Mantel test result comparing distance matrices
#'     \item coinertia: coinertia analysis result (NULL if insufficient data)
#'   }
#'
#' @details
#' The function implements:
#' - Jaccard distance calculations for presence/absence data
#' - Mantel test to assess congruence between methods (Mantel 1967)
#' - Coinertia analysis for multivariate comparison (Doledec & Chessel 1994)
#' - Spearman rank correlation for richness comparison
#'
#' Requires minimum 3 common sites for meaningful analysis. Uses only sites
#' present in both datasets to ensure valid comparisons.
#'
#' The Mantel test evaluates whether sites similar in composition according
#' to one method are also similar according to the other method. A significant
#' positive correlation indicates the methods detect similar biogeographic
#' patterns.
#'
#' @examples
#' \dontrun{
#' # Compare bioacoustic and metabarcoding datasets
#' similarity_results <- analyze_site_similarity(
#'   acoustic_data = katydid_acoustic,
#'   metabarcoding_data = orthoptera_metabar
#' )
#' 
#' # Extract Mantel test results
#' mantel_r <- similarity_results$mantel_test$statistic
#' mantel_p <- similarity_results$mantel_test$signif
#' 
#' # Plot richness comparison
#' library(ggplot2)
#' ggplot(similarity_results$site_richness, 
#'        aes(x = acoustic_richness, y = metabar_richness)) +
#'   geom_point() +
#'   geom_smooth(method = "lm")
#' }
#'
#' @references
#' Mantel, N. (1967). The detection of disease clustering and a generalized
#' regression approach. Cancer Research, 27(2), 209-220.
#'
#' Doledec, S., & Chessel, D. (1994). Co-inertia analysis: an alternative
#' method for studying species-environment relationships. Freshwater Biology,
#' 31(3), 277-294.
#'
#' Legendre, P., & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#'
#' @seealso \code{\link{analyze_vegetation_effects}}, \code{\link[vegan]{mantel}}
#'
#' @export
analyze_site_similarity <- function(acoustic_data, metabarcoding_data) {
  
  # Identify common sites
  acoustic_sites <- acoustic_data$presence_matrix$site
  metabar_sites <- metabarcoding_data$presence_matrix$site
  common_sites <- intersect(acoustic_sites, metabar_sites)
  
  cat(sprintf("Common sites for comparison: %d\n", length(common_sites)))
  cat(sprintf("Acoustic sites: %s\n", paste(acoustic_sites, collapse = ", ")))
  cat(sprintf("Metabarcoding sites: %s\n", paste(metabar_sites, collapse = ", ")))
  
  if (length(common_sites) < 3) {
    cat("Not enough common sites for similarity analysis\n")
    acoustic_dist <- NULL
    metabar_dist <- NULL
    mantel_result <- NULL
  } else {
    # Filter matrices to common sites only
    acoustic_matrix <- acoustic_data$presence_matrix %>%
      filter(site %in% common_sites) %>%
      arrange(site) %>%
      select(-site)
    
    metabar_matrix <- metabarcoding_data$presence_matrix %>%
      filter(site %in% common_sites) %>%
      arrange(site) %>%
      select(-site)
    
    # Calculate distance matrices
    if (ncol(acoustic_matrix) > 0 && ncol(metabar_matrix) > 0) {
      acoustic_dist <- vegdist(acoustic_matrix, method = "jaccard")
      metabar_dist <- vegdist(metabar_matrix, method = "jaccard")
      
      # Mantel test for correlation between distance matrices
      mantel_result <- mantel(acoustic_dist, metabar_dist, method = "spearman")
    } else {
      acoustic_dist <- NULL
      metabar_dist <- NULL
      mantel_result <- NULL
    }
  }
  
  # Coinertia analysis if both datasets available
  coinertia_result <- NULL
  if (!is.null(acoustic_dist) && nrow(acoustic_matrix) > 2 && 
      ncol(acoustic_matrix) > 1 && ncol(metabar_matrix) > 1) {
    
    tryCatch({
      # Prepare data for coinertia
      acoustic_pca <- dudi.pca(as.data.frame(acoustic_matrix), 
                               scannf = FALSE, nf = 2)
      metabar_pca <- dudi.pca(as.data.frame(metabar_matrix), 
                              scannf = FALSE, nf = 2)
      
      # Coinertia analysis
      coinertia_result <- coinertia(acoustic_pca, metabar_pca, 
                                    scannf = FALSE, nf = 2)
    }, error = function(e) {
      cat("Coinertia analysis failed - not enough data or too much similarity\n")
    })
  }
  
  # Species richness per site - use common sites only
  if (length(common_sites) >= 3) {
    site_richness <- data.frame(
      site = common_sites,
      metabar_richness = rowSums(metabar_matrix),
      acoustic_richness = rowSums(acoustic_matrix)
    )
  } else {
    # Fallback if not enough common sites
    site_richness <- data.frame(
      site = metabarcoding_data$presence_matrix$site,
      metabar_richness = rowSums(metabarcoding_data$presence_matrix %>% select(-site))
    ) %>%
      left_join(
        data.frame(
          site = acoustic_data$presence_matrix$site,
          acoustic_richness = rowSums(acoustic_data$presence_matrix %>% select(-site))
        ), by = "site"
      ) %>%
      replace_na(list(acoustic_richness = 0))
  }
  
  # Correlation between richness
  if (sum(site_richness$acoustic_richness > 0) > 3) {
    richness_cor <- cor.test(
      site_richness$acoustic_richness,
      site_richness$metabar_richness,
      method = "spearman"
    )
  } else {
    richness_cor <- NULL
  }
  
  return(list(
    site_richness = site_richness,
    richness_correlation = richness_cor,
    mantel_test = mantel_result,
    coinertia = coinertia_result
  ))
}


# ==============================================================================
# SECTION 2: VEGETATION EFFECTS
# ==============================================================================

#' Analyze vegetation effects on species richness
#'
#' Examines relationships between vegetation structure (trees, lianas) and
#' species richness detected by either acoustic or metabarcoding methods.
#' Includes correlation analyses and linear modeling.
#'
#' @param species_data List. Species detection data with presence_matrix component
#' @param vegetation_data List. Vegetation data with summary component (from
#'   read_vegetation_data())
#' @param data_type Character. Type of data ("acoustic" or "metabarcoding") for
#'   labeling. Default: "acoustic"
#'
#' @return List with components:
#'   \itemize{
#'     \item correlations: data frame with Spearman correlations between species
#'       richness and each vegetation variable
#'     \item model: linear model predicting species richness from vegetation
#'       (NULL if insufficient data)
#'     \item richness_plot: ggplot object showing richness vs vegetation
#'     \item data: merged dataset with species and vegetation data
#'   }
#'
#' @details
#' Vegetation variables analyzed:
#' - tree_species_richness: number of tree species
#' - tree_abundance: total number of trees
#' - tree_total_basal_area: cumulative basal area
#' - liana_species_richness: number of liana species
#' - liana_rooted_stems: number of rooted liana stems
#' - liana_total_basal_area: cumulative liana basal area
#'
#' Uses Spearman correlations (non-parametric) suitable for ecological data
#' which often violate normality assumptions. Linear models include main
#' effects only (no interactions).
#'
#' Site matching: Converts site codes (S01) to Plot IDs (YB-P01) for
#' ForestGEO vegetation plot matching.
#'
#' @examples
#' \dontrun{
#' # Analyze vegetation effects on katydid richness
#' veg_effects <- analyze_vegetation_effects(
#'   species_data = katydid_data,
#'   vegetation_data = forestgeo_data,
#'   data_type = "katydid"
#' )
#' 
#' # View significant correlations
#' sig_cors <- veg_effects$correlations %>%
#'   filter(p_value < 0.05)
#' 
#' # Display plot
#' print(veg_effects$richness_plot)
#' }
#'
#' @references
#' Basset, Y., Lamarre, G. P. A., et al. (2015). Arthropod distribution in a
#' tropical rainforest: Tackling a four dimensional puzzle. PLoS ONE, 10(12),
#' e0144110.
#'
#' Condit, R. (1998). Tropical Forest Census Plots. Springer-Verlag and
#' R. G. Landes Company.
#'
#' @seealso \code{\link{analyze_site_similarity}}, 
#'   \code{\link{analyze_nonlinear_relationships}}
#'
#' @export
analyze_vegetation_effects <- function(species_data, vegetation_data, 
                                       data_type = "acoustic") {
  
  # Merge data, convert site to Plot format for vegetation matching
  if ("site" %in% colnames(species_data$presence_matrix)) {
    species_veg_data <- species_data$presence_matrix %>%
      mutate(
        site_number = as.numeric(gsub("S", "", site)),
        Plot = sprintf("YB-P%02d", site_number)
      ) %>%
      inner_join(vegetation_data$summary, by = "Plot")
  } else {
    warning("No site column found in species data")
    return(NULL)
  }
  
  if (nrow(species_veg_data) == 0) {
    warning(paste("Unable to merge", data_type, "and vegetation data"))
    return(NULL)
  }
  
  # Calculate species richness
  species_cols <- setdiff(colnames(species_data$presence_matrix), c("site"))
  species_richness_data <- species_veg_data %>%
    mutate(species_richness = rowSums(select(., all_of(species_cols)), na.rm = TRUE))
  
  # Analyze correlations with vegetation variables
  veg_vars <- c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area")
  
  correlations <- data.frame(
    variable = character(),
    correlation = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (veg_var in veg_vars) {
    if (veg_var %in% colnames(species_richness_data) && 
        sum(!is.na(species_richness_data[[veg_var]])) > 3) {
      
      cor_test <- cor.test(species_richness_data$species_richness, 
                           species_richness_data[[veg_var]],
                           method = "spearman", 
                           exact = FALSE)
      
      correlations <- rbind(correlations, data.frame(
        variable = veg_var,
        correlation = cor_test$estimate,
        p_value = cor_test$p.value
      ))
    }
  }
  
  # Model species richness as function of vegetation structure
  if (nrow(species_richness_data) > 5) {
    model <- lm(species_richness ~ tree_species_richness + 
                  tree_abundance + 
                  liana_species_richness + 
                  liana_rooted_stems, 
                data = species_richness_data)
  } else {
    model <- NULL
  }
  
  # Visualize relationships
  richness_plot <- ggplot(species_richness_data, 
                          aes(x = liana_species_richness, 
                              y = species_richness)) +
    geom_point(aes(size = liana_rooted_stems, color = tree_species_richness)) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = paste(stringr::str_to_title(data_type), 
                       "Species Richness vs Vegetation Structure"),
         x = "Liana Species Richness", 
         y = paste(stringr::str_to_title(data_type), "Species Richness")) +
    theme_minimal()
  
  return(list(
    correlations = correlations,
    model = model,
    richness_plot = richness_plot,
    data = species_richness_data
  ))
}


# ==============================================================================
# SECTION 3: VARIABLE SELECTION
# ==============================================================================

#' Select optimal environmental variables for ordination
#'
#' Performs forward selection to identify environmental variables that
#' significantly explain species composition patterns. Uses permutation
#' tests within RDA framework.
#'
#' @param species_counts Data frame. Detection count matrix merged with
#'   environmental data
#' @param vegetation_enhanced Data frame. Enhanced vegetation variables from
#'   create_enhanced_environmental_variables()
#' @param method Character. Selection method - currently only "rda" supported.
#'   Default: "rda"
#'
#' @return List with components:
#'   \itemize{
#'     \item selected_variables: character vector of selected variable names
#'     \item model: final RDA model with selected variables
#'     \item species_matrix: species data matrix used
#'     \item env_matrix: environmental data matrix (selected variables only)
#'     \item site_data: complete merged site data
#'   }
#'   Returns NULL if insufficient data for selection.
#'
#' @details
#' The function implements:
#' - Forward selection using ordistep() from vegan package
#' - RDA (Redundancy Analysis) as base model
#' - Permutation tests (199 permutations) for variable significance
#' - Automatic filtering of species with <2 detections
#' - Removal of near-constant environmental variables
#'
#' Variable selection criteria:
#' - Minimum 3 species in matrix
#' - Minimum 2 environmental variables
#' - Minimum 4 sites with complete data
#' - Environmental variables with variance > 0.001
#'
#' The procedure follows Blanchet et al. (2008) recommendations for
#' forward selection in ordination.
#'
#' @examples
#' \dontrun{
#' # Create enhanced variables first
#' veg_enhanced <- create_enhanced_environmental_variables(vegetation_data)
#' 
#' # Select optimal variables
#' selected <- select_optimal_variables(
#'   species_counts = species_detection_matrix,
#'   vegetation_enhanced = veg_enhanced
#' )
#' 
#' # View selected variables
#' print(selected$selected_variables)
#' 
#' # Use in subsequent analyses
#' rda_final <- selected$model
#' }
#'
#' @references
#' Blanchet, F. G., Legendre, P., & Borcard, D. (2008). Forward selection
#' of explanatory variables. Ecology, 89(9), 2623-2632.
#'
#' Legendre, P., & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#'
#' Oksanen, J., et al. (2022). vegan: Community Ecology Package.
#'
#' @seealso \code{\link{create_enhanced_environmental_variables}}, 
#'   \code{\link[vegan]{ordistep}}
#'
#' @export
select_optimal_variables <- function(species_counts, vegetation_enhanced, 
                                     method = "rda") {
  
  cat("Selecting optimal environmental variables...\n")
  
  # Prepare data
  sites_env <- species_counts %>%
    mutate(
      site_number = as.numeric(gsub("S", "", site)),
      Plot = sprintf("YB-P%02d", site_number)
    ) %>%
    left_join(vegetation_enhanced, by = "Plot") %>%
    filter(complete.cases(.))
  
  # Species matrix (filter rare species)
  species_matrix <- sites_env %>%
    select(-c(site, site_number, Plot)) %>%
    select_if(is.numeric) %>%
    select(-starts_with("tree_"), -starts_with("liana_")) %>%
    select_if(~ sum(.) >= 2)  # Keep species with at least 2 detections
  
  # Environmental variables (exclude redundant originals)
  env_vars <- sites_env %>%
    select(total_plant_richness, liana_dominance, structural_complexity,
           log_tree_abundance, log_liana_stems, liana_to_tree_ratio,
           basal_area_ratio, tree_density, liana_density) %>%
    select_if(~ var(., na.rm = TRUE) > 0.001)  # Remove near-constant variables
  
  cat(sprintf("Available variables: %d species, %d environmental variables\n", 
              ncol(species_matrix), ncol(env_vars)))
  
  # Forward selection with ordistep
  if (ncol(species_matrix) >= 3 && ncol(env_vars) >= 2 && nrow(sites_env) >= 4) {
    
    tryCatch({
      # Null and full models
      rda_null <- rda(log1p(species_matrix) ~ 1, data = env_vars)
      rda_full <- rda(log1p(species_matrix) ~ ., data = env_vars)
      
      # Forward selection
      rda_selected <- ordistep(rda_null, scope = formula(rda_full), 
                               direction = "forward", 
                               permutations = 199,
                               trace = FALSE)
      
      selected_vars <- labels(terms(rda_selected))
      cat(sprintf("Selected variables: %s\n", paste(selected_vars, collapse = ", ")))
      
      return(list(
        selected_variables = selected_vars,
        model = rda_selected,
        species_matrix = species_matrix,
        env_matrix = env_vars[, selected_vars, drop = FALSE],
        site_data = sites_env
      ))
      
    }, error = function(e) {
      cat(sprintf("Variable selection failed: %s\n", e$message))
      return(list(
        selected_variables = names(env_vars)[1:min(3, ncol(env_vars))],
        species_matrix = species_matrix,
        env_matrix = env_vars[, 1:min(3, ncol(env_vars)), drop = FALSE],
        site_data = sites_env
      ))
    })
  } else {
    return(NULL)
  }
}


#' Create enhanced environmental variables
#'
#' Generates composite and transformed vegetation variables for improved
#' ecological interpretation and reduced collinearity in statistical models.
#'
#' @param vegetation_data List. Output from read_vegetation_data() with
#'   summary component
#'
#' @return Data frame with original plus enhanced variables:
#'   Composite variables:
#'   \itemize{
#'     \item total_plant_richness: tree + liana species richness
#'     \item liana_dominance: proportion of total richness from lianas
#'     \item structural_complexity: sqrt(liana_stems * tree_abundance)
#'   }
#'   Transformed variables:
#'   \itemize{
#'     \item log_tree_abundance: log(tree_abundance + 1)
#'     \item log_liana_stems: log(liana_stems + 1)
#'     \item sqrt_tree_richness: sqrt(tree_richness)
#'   }
#'   Ratio variables:
#'   \itemize{
#'     \item liana_to_tree_ratio: liana / tree species richness
#'     \item basal_area_ratio: liana / tree basal area
#'     \item tree_density: trees per species
#'     \item liana_density: liana stems per species
#'   }
#'
#' @details
#' Variable transformations serve multiple purposes:
#' - Log transformations: normalize right-skewed abundance data
#' - Square root: moderate transformation for count data
#' - Ratios: capture relative importance and balance
#' - Composite indices: represent multivariate concepts
#'
#' Handles edge cases:
#' - Adds 1 before log to avoid log(0)
#' - Replaces Inf/NaN values with 0
#' - Safe division (adds small constant to denominator)
#'
#' These derived variables often show stronger relationships with biodiversity
#' than raw measurements and reduce multicollinearity among predictors.
#'
#' @examples
#' \dontrun{
#' # Create enhanced variables
#' veg_data <- read_vegetation_data("vegetation.xlsx")
#' veg_enhanced <- create_enhanced_environmental_variables(veg_data)
#' 
#' # Use in modeling
#' model <- lm(species_richness ~ structural_complexity + 
#'             liana_dominance + log_tree_abundance,
#'             data = merged_data)
#' }
#'
#' @references
#' Legendre, P., & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#' (Chapter 1: Variables and transformations)
#'
#' @seealso \code{\link{select_optimal_variables}}, \code{\link{read_vegetation_data}}
#'
#' @export
create_enhanced_environmental_variables <- function(vegetation_data) {
  
  enhanced_env <- vegetation_data$summary %>%
    mutate(
      # Composite variables
      total_plant_richness = tree_species_richness + liana_species_richness,
      liana_dominance = liana_species_richness / 
        (tree_species_richness + liana_species_richness + 1),
      structural_complexity = sqrt(liana_rooted_stems * tree_abundance),
      
      # Transformations for normalization
      log_tree_abundance = log1p(tree_abundance),
      log_liana_stems = log1p(liana_rooted_stems),
      sqrt_tree_richness = sqrt(tree_species_richness),
      
      # Ecologically relevant ratios
      liana_to_tree_ratio = liana_species_richness / tree_species_richness,
      basal_area_ratio = liana_total_basal_area / (tree_total_basal_area + 1),
      
      # Density variables
      tree_density = tree_abundance / tree_species_richness,
      liana_density = liana_rooted_stems / liana_species_richness
    ) %>%
    # Replace Inf/NaN with 0
    mutate(across(everything(), ~ ifelse(is.infinite(.) | is.nan(.), 0, .)))
  
  return(enhanced_env)
}


# ==============================================================================
# SECTION 4: NMDS ORDINATION WITH DETECTION COUNTS
# ==============================================================================

#' Create detection count matrix from raw detection data
#'
#' Transforms raw acoustic detection data into a site-by-species matrix of
#' detection counts (abundance data) suitable for multivariate analyses.
#'
#' @param raw_detections Data frame. Raw detection data with columns for sites,
#'   species, and detection metadata
#' @param species_col Character. Name of column containing species identifiers.
#'   Default: "common_name"
#' @param site_col Character. Name of column containing site identifiers.
#'   Default: "site"
#'
#' @return Data frame with site column plus one column per species containing
#'   detection counts. Returns empty data frame with site column if no detections.
#'
#' @details
#' The function aggregates raw detection records into a summary matrix where:
#' - Rows represent sampling sites
#' - Columns represent species (except first column which is site ID)
#' - Values are integer counts of detections
#'
#' This format is suitable for:
#' - NMDS ordination with Bray-Curtis distance
#' - CCA with abundance data
#' - GAM models of detection rates
#' - Diversity metrics based on counts
#'
#' Multiple detections of the same species at the same site are counted
#' separately, providing information about detection frequency/calling
#' intensity in addition to presence/absence.
#'
#' @examples
#' \dontrun{
#' # Create count matrix from acoustic detections
#' count_matrix <- create_detection_count_matrix(
#'   raw_detections = katydid_detections,
#'   species_col = "common_name",
#'   site_col = "site"
#' )
#' 
#' # View structure
#' str(count_matrix)
#' # > data.frame: 12 obs. of 31 variables (site + 30 species)
#' }
#'
#' @seealso \code{\link{analyze_nmds_with_detection_counts}},
#'   \code{\link{analyze_cca_with_detections}}
#'
#' @export
create_detection_count_matrix <- function(raw_detections, 
                                          species_col = "common_name", 
                                          site_col = "site") {
  
  cat("Creating detection count matrix...\n")
  
  if (nrow(raw_detections) == 0) {
    return(data.frame(site = character()))
  }
  
  # Count detections per site and species
  detection_counts <- raw_detections %>%
    group_by(across(all_of(c(site_col, species_col)))) %>%
    summarise(detections = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(species_col), 
                values_from = detections, 
                values_fill = 0)
  
  cat(sprintf("Detection matrix created: %d sites, %d species\n", 
              nrow(detection_counts), ncol(detection_counts) - 1))
  
  return(detection_counts)
}


#' NMDS ordination with detection count data
#'
#' Performs Non-metric Multidimensional Scaling (NMDS) ordination using
#' detection count matrices from bioacoustic monitoring. Integrates vegetation
#' data and visualizes community patterns with environmental gradients.
#'
#' @param katydid_data List. Output from process/combine katydid deployments,
#'   must contain raw_detections component
#' @param bird_data List. Output from process/combine bird deployments,
#'   must contain detections component
#' @param vegetation_data List. Vegetation data with summary component from
#'   read_vegetation_data()
#' @param output_dir Character. Directory path for saving results and plots
#' @param analysis_name Character. Name prefix for output files.
#'   Default: "Detection_Counts_NMDS"
#'
#' @return List with components:
#'   \itemize{
#'     \item nmds: metaMDS result object with ordination solution
#'     \item scores: data frame with NMDS scores and site metadata
#'     \item env_fit: envfit result showing environmental correlations
#'     \item plot: ggplot object - basic NMDS with detection abundances
#'     \item plot_with_vectors: ggplot object - NMDS with env vectors overlaid
#'   }
#'   Returns NULL if insufficient data or NMDS fails to converge.
#'
#' @details
#' **Ordination method:**
#' - Uses Bray-Curtis dissimilarity on log(x+1) transformed counts
#' - 2 dimensions, up to 100 attempts, max 1000 iterations
#' - Stress < 0.2 considered acceptable (Legendre & Legendre 2012)
#' - No autotransformation (transformation done explicitly)
#'
#' **Data preparation:**
#' - Creates detection count matrices for katydids and birds
#' - Identifies common sites with complete environmental data
#' - Filters rare species (< 2 total detections)
#' - Combines filtered katydid and bird matrices
#' - Log-transforms counts to stabilize variance
#'
#' **Environmental fitting:**
#' - Tests 8 variables: katydid/bird totals + 6 vegetation metrics
#' - Uses envfit() with 999 permutations
#' - Only variables with >3 non-missing values included
#' - Vectors scaled to 0.5-0.6 for visualization clarity
#'
#' **Visualization:**
#' - Point size = katydid detection abundance
#' - Point color = bird detection abundance
#' - Significant environmental vectors (p < 0.05) shown as arrows
#' - Labels show site IDs
#'
#' **Output files:**
#' - nmds_detection_counts.jpg: base ordination plot
#' - nmds_detection_counts_with_vectors.jpg: with env gradients
#' - nmds_detection_envfit_results.txt: full envfit statistics
#'
#' **Minimum requirements:**
#' - At least 3 common sites with complete data
#' - At least 3 species after filtering
#' - NMDS must converge (stress acceptable)
#'
#' @examples
#' \dontrun{
#' # Perform NMDS with detection counts
#' nmds_results <- analyze_nmds_with_detection_counts(
#'   katydid_data = all_katydids,
#'   bird_data = all_birds,
#'   vegetation_data = forestgeo_vegetation,
#'   output_dir = "results/nmds/"
#' )
#' 
#' # Check convergence
#' nmds_results$nmds$stress  # Should be < 0.2
#' 
#' # View significant environmental correlations
#' nmds_results$env_fit
#' 
#' # Display ordination plot
#' print(nmds_results$plot_with_vectors)
#' }
#'
#' @references
#' Legendre, P., & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#'
#' Oksanen, J., et al. (2022). vegan: Community Ecology Package. R package
#' version 2.6-4. https://CRAN.R-project.org/package=vegan
#'
#' Kruskal, J. B. (1964). Nonmetric multidimensional scaling: A numerical
#' method. Psychometrika, 29(2), 115-129.
#'
#' @seealso \code{\link{create_detection_count_matrix}},
#'   \code{\link{analyze_cca_with_detections}},
#'   \code{\link[vegan]{metaMDS}}, \code{\link[vegan]{envfit}}
#'
#' @export
analyze_nmds_with_detection_counts <- function(katydid_data, bird_data, 
                                               vegetation_data, output_dir, 
                                               analysis_name = "Detection_Counts_NMDS") {
  
  cat(sprintf("\n=== %s ANALYSIS ===\n", analysis_name))
  
  # Create detection count matrices
  if (!is.null(katydid_data$raw_detections)) {
    katydid_counts <- create_detection_count_matrix(katydid_data$raw_detections, 
                                                    "common_name", "site")
  } else {
    cat("[X] No katydid raw detections available\n")
    return(NULL)
  }
  
  if (!is.null(bird_data$detections)) {
    bird_counts <- create_detection_count_matrix(bird_data$detections, 
                                                 "common_name", "site")
  } else {
    cat("[X] No bird raw detections available\n")
    return(NULL)
  }
  
  # Prepare data for NMDS - combine katydids and birds with common sites
  common_sites <- intersect(katydid_counts$site, bird_counts$site)
  
  if (length(common_sites) < 3) {
    cat("[X] Not enough common sites for NMDS analysis\n")
    return(NULL)
  }
  
  cat(sprintf("Common sites for analysis: %d\n", length(common_sites)))
  
  # Calculate total detections per site
  katydid_totals <- katydid_counts %>%
    filter(site %in% common_sites) %>%
    mutate(katydid_total_detections = rowSums(select(., -site))) %>%
    select(site, katydid_total_detections)
  
  bird_totals <- bird_counts %>%
    filter(site %in% common_sites) %>%
    mutate(bird_total_detections = rowSums(select(., -site))) %>%
    select(site, bird_total_detections)
  
  # Combine with vegetation data
  detection_env_data <- katydid_totals %>%
    left_join(bird_totals, by = "site") %>%
    mutate(
      site_number = as.numeric(gsub("S", "", site)),
      Plot = sprintf("YB-P%02d", site_number)
    ) %>%
    left_join(vegetation_data$summary, by = "Plot") %>%
    filter(complete.cases(.))
  
  cat(sprintf("Sites with complete data: %d\n", nrow(detection_env_data)))
  
  if (nrow(detection_env_data) < 3) {
    cat("[X] Not enough sites with complete data\n")
    return(NULL)
  }
  
  # Prepare species matrix for NMDS (species with > 2 detections)
  katydid_species_matrix <- katydid_counts %>%
    filter(site %in% detection_env_data$site) %>%
    arrange(site) %>%
    select(-site)
  
  bird_species_matrix <- bird_counts %>%
    filter(site %in% detection_env_data$site) %>%
    arrange(site) %>%
    select(-site)
  
  # Filter rare species (< 2 total detections)
  katydid_species_filtered <- katydid_species_matrix[, 
                                                     colSums(katydid_species_matrix) >= 2, drop = FALSE]
  bird_species_filtered <- bird_species_matrix[, 
                                               colSums(bird_species_matrix) >= 2, drop = FALSE]
  
  # Combine matrices
  combined_species_matrix <- cbind(katydid_species_filtered, bird_species_filtered)
  
  cat(sprintf("Species matrix: %d sites, %d katydid species, %d bird species\n", 
              nrow(combined_species_matrix), ncol(katydid_species_filtered), 
              ncol(bird_species_filtered)))
  
  if (ncol(combined_species_matrix) < 3) {
    cat("[X] Not enough species for NMDS\n")
    return(NULL)
  }
  
  # Execute NMDS
  tryCatch({
    # Log transformation for count data
    log_transformed_matrix <- log1p(combined_species_matrix)
    
    nmds_result <- metaMDS(log_transformed_matrix, 
                           distance = "bray",  # Bray-Curtis for count data
                           k = 2,
                           trymax = 100,
                           maxit = 1000,
                           autotransform = FALSE,
                           trace = 1)
    
    if (!nmds_result$converged) {
      cat("[X] NMDS did not converge\n")
      return(NULL)
    }
    
    cat("[OK] NMDS converged successfully!\n")
    
    # Extract NMDS scores
    nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
    nmds_scores$site <- detection_env_data$site
    nmds_scores$site_number <- detection_env_data$site_number
    nmds_scores$Plot <- detection_env_data$Plot
    
    # Merge with environmental and detection data
    nmds_data <- nmds_scores %>%
      left_join(detection_env_data, by = c("site", "site_number", "Plot"))
    
    # Create NMDS plot
    nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(size = katydid_total_detections, 
                     color = bird_total_detections), alpha = 0.8) +
      geom_text(aes(label = site), vjust = -1.5, size = 3) +
      scale_size_continuous(name = "Katydid\nDetections", range = c(2, 8)) +
      scale_color_viridis_c(name = "Bird\nDetections") +
      labs(
        title = "NMDS: Species Communities with Detection Abundances",
        subtitle = paste("Stress =", round(nmds_result$stress, 3), 
                         "| Size = Katydid detections, Color = Bird detections"),
        x = "NMDS1", 
        y = "NMDS2"
      ) +
      theme_minimal() +
      theme(legend.position = "right")
    
    # Environmental fitting with detection totals as explanatory variables
    env_vars <- detection_env_data %>%
      select(katydid_total_detections, bird_total_detections, 
             tree_species_richness, tree_abundance, tree_total_basal_area,
             liana_species_richness, liana_rooted_stems, liana_total_basal_area) %>%
      select_if(~ sum(!is.na(.)) > 3)  # Keep only variables with enough data
    
    env_fit <- envfit(nmds_result, env_vars, perm = 999)
    
    # Extract significance and filter vectors
    cat("  Environmental fit results:\n")
    print(env_fit)
    
    # Filter to keep only significant vectors (p < 0.05)
    # Use 0.1 threshold for more permissive inclusion if needed
    sig_threshold <- 0.05
    
    nmds_plot_with_vectors <- nmds_plot
    
    if (!is.null(env_fit$vectors) && any(env_fit$vectors$pvals < sig_threshold)) {
      
      # Get significant vectors only
      sig_idx <- env_fit$vectors$pvals < sig_threshold
      sig_vectors <- env_fit$vectors$arrows[sig_idx, , drop = FALSE]
      sig_pvals <- env_fit$vectors$pvals[sig_idx]
      
      cat(sprintf("  [OK] %d significant environmental vectors (p < %.2f):\n", 
                  sum(sig_idx), sig_threshold))
      for (i in seq_along(sig_pvals)) {
        cat(sprintf("       - %s (p = %.3f)\n", rownames(sig_vectors)[i], sig_pvals[i]))
      }
      
      # Add arrows for significant vectors only
      nmds_plot_with_vectors <- nmds_plot_with_vectors +
        geom_segment(data = data.frame(sig_vectors * 0.7, 
                                       variable = rownames(sig_vectors)),
                     aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                     arrow = arrow(length = unit(0.3, "cm")), 
                     color = "red", linewidth = 1.2, inherit.aes = FALSE) +
        geom_text(data = data.frame(sig_vectors * 0.8, 
                                    variable = rownames(sig_vectors)),
                  aes(x = NMDS1, y = NMDS2, label = variable),
                  color = "red", size = 3.5, fontface = "bold", inherit.aes = FALSE)
      
    } else {
      cat(sprintf("  [WARNING] No significant environmental vectors (all p >= %.2f)\n", sig_threshold))
      if (!is.null(env_fit$vectors)) {
        cat("  Vectors p-values:", paste(round(env_fit$vectors$pvals, 3), collapse = ", "), "\n")
      }
    }
    
    # Save results
    ggsave(file.path(output_dir, "nmds_detection_counts.jpg"), 
           nmds_plot, width = 12, height = 8)
    ggsave(file.path(output_dir, "nmds_detection_counts_with_vectors.jpg"), 
           nmds_plot_with_vectors, width = 12, height = 8)
    
    # Save environmental fitting results
    capture.output(env_fit, 
                   file = file.path(output_dir, "nmds_detection_envfit_results.txt"))
    
    return(list(
      nmds = nmds_result,
      scores = nmds_data,
      env_fit = env_fit,
      plot = nmds_plot,
      plot_with_vectors = nmds_plot_with_vectors
    ))
    
  }, error = function(e) {
    cat(sprintf("[X] NMDS failed: %s\n", e$message))
    return(NULL)
  })
}


# ==============================================================================
# SECTION 5: CCA ORDINATION WITH ENVIRONMENTAL CONSTRAINTS
# ==============================================================================

#' Canonical Correspondence Analysis (CCA) with detection counts
#'
#' Performs constrained ordination to analyze relationships between species
#' composition (based on detection counts) and environmental variables. CCA
#' identifies environmental gradients that best explain community patterns.
#'
#' @param katydid_data List. Output from process/combine katydid deployments,
#'   must contain raw_detections component
#' @param bird_data List. Output from process/combine bird deployments
#'   (currently not used in this function, but included for consistency)
#' @param vegetation_data List. Vegetation data with summary component from
#'   read_vegetation_data()
#' @param output_dir Character. Directory path for saving results and plots
#' @param analysis_name Character. Name prefix for output files.
#'   Default: "Detection_CCA"
#'
#' @return List with components:
#'   \itemize{
#'     \item cca: cca result object from vegan package
#'     \item plot: ggplot object - CCA ordination with environmental vectors
#'     \item site_scores: data frame with CCA axis scores and site metadata
#'     \item species_scores: data frame with species positions in CCA space
#'     \item env_scores: data frame with environmental variable vectors (biplot scores)
#'   }
#'   Returns NULL if insufficient data or CCA fails.
#'
#' @details
#' **Method overview:**
#' Canonical Correspondence Analysis (CCA) is a constrained ordination method
#' that directly relates species composition to environmental variables. Unlike
#' NMDS (unconstrained), CCA axes are linear combinations of environmental
#' predictors that maximize the separation of species distributions.
#'
#' **Data preparation:**
#' - Creates detection count matrix for katydids
#' - Filters rare species (< 2 total detections)
#' - Log(x+1) transformation to stabilize variance
#' - Matches sites with complete environmental data
#' - Uses 6 vegetation variables: tree and liana richness, abundance, and basal area
#'
#' **CCA implementation:**
#' - Formula: species ~ all environmental variables
#' - Uses vegan::cca() function
#' - Extracts site, species, and biplot scores
#' - Biplot vectors scaled by factor of 3 for visualization
#'
#' **Interpretation:**
#' - Axis importance shown as % of constrained variance explained
#' - Species scores: center of species distribution along gradients
#' - Site scores: community composition position
#' - Environmental vectors: direction and strength of gradients
#' - Longer arrows = stronger environmental effects
#' - Vector direction = gradient orientation
#' - Sites at arrow tip have high values of that variable
#'
#' **Minimum requirements:**
#' - At least 3 sites with complete environmental data
#' - At least 3 species after filtering
#' - At least 1 environmental variable
#'
#' **Output files:**
#' - cca_detections_environment.jpg: ordination plot with vectors
#' - cca_summary.txt: full CCA results including significance tests
#'
#' @examples
#' \dontrun{
#' # Perform CCA with katydid detections
#' cca_results <- analyze_cca_with_detections(
#'   katydid_data = all_katydids,
#'   bird_data = all_birds,
#'   vegetation_data = forestgeo_vegetation,
#'   output_dir = "results/cca/"
#' )
#' 
#' # Check variance explained
#' summary(cca_results$cca)
#' 
#' # View constrained variance
#' cca_results$cca$CCA$tot.chi / cca_results$cca$tot.chi  # Proportion
#' 
#' # Test significance
#' anova(cca_results$cca, permutations = 999)
#' 
#' # Display ordination plot
#' print(cca_results$plot)
#' }
#'
#' @references
#' ter Braak, C. J. F. (1986). Canonical Correspondence Analysis: A new
#' eigenvector technique for multivariate direct gradient analysis. Ecology,
#' 67(5), 1167-1179.
#'
#' ter Braak, C. J. F., & Smilauer, P. (2002). CANOCO Reference Manual and
#' CanoDraw for Windows User's Guide: Software for Canonical Community
#' Ordination (version 4.5). Microcomputer Power, Ithaca, NY.
#'
#' Legendre, P., & Legendre, L. (2012). Numerical Ecology (3rd ed.). Elsevier.
#' (Chapter 11: Canonical analysis)
#'
#' Oksanen, J., et al. (2022). vegan: Community Ecology Package.
#' https://CRAN.R-project.org/package=vegan
#'
#' @seealso \code{\link{analyze_nmds_with_detection_counts}},
#'   \code{\link{create_detection_count_matrix}},
#'   \code{\link[vegan]{cca}}
#'
#' @export
analyze_cca_with_detections <- function(katydid_data, bird_data, 
                                        vegetation_data, output_dir, 
                                        analysis_name = "Detection_CCA") {
  
  cat(sprintf("\n=== %s ANALYSIS ===\n", analysis_name))
  
  # Create detection count matrix
  if (!is.null(katydid_data$raw_detections)) {
    katydid_counts <- create_detection_count_matrix(katydid_data$raw_detections, 
                                                    "common_name", "site")
  } else {
    cat("[X] No katydid raw detections available\n")
    return(NULL)
  }
  
  # Prepare environmental data
  common_sites <- katydid_counts$site
  
  env_data <- data.frame(
    site = common_sites,
    site_number = as.numeric(gsub("S", "", common_sites)),
    Plot = sprintf("YB-P%02d", as.numeric(gsub("S", "", common_sites)))
  ) %>%
    left_join(vegetation_data$summary, by = "Plot") %>%
    filter(complete.cases(.))
  
  cat(sprintf("Sites with complete data: %d\n", nrow(env_data)))
  
  if (nrow(env_data) < 3) {
    cat("[X] Not enough sites with environmental data\n")
    return(NULL)
  }
  
  # Species matrix (log-transformed)
  species_matrix <- katydid_counts %>%
    filter(site %in% env_data$site) %>%
    arrange(match(site, env_data$site)) %>%
    select(-site)
  
  # Filter rare species
  species_matrix_filtered <- species_matrix[, 
                                            colSums(species_matrix) >= 2, drop = FALSE]
  
  # Log(x+1) transformation
  species_log <- log1p(species_matrix_filtered)
  
  # Environmental variables
  env_matrix <- env_data %>%
    select(tree_species_richness, tree_abundance, tree_total_basal_area,
           liana_species_richness, liana_rooted_stems, liana_total_basal_area) %>%
    select_if(~ sum(!is.na(.)) > 3)
  
  cat(sprintf("CCA input: %d sites, %d species, %d environmental variables\n", 
              nrow(species_log), ncol(species_log), ncol(env_matrix)))
  
  if (ncol(species_log) < 3 || ncol(env_matrix) < 1) {
    cat("[X] Insufficient data for CCA\n")
    return(NULL)
  }
  
  # Execute CCA
  tryCatch({
    cca_result <- cca(species_log ~ ., data = env_matrix)
    
    cat("[OK] CCA completed successfully!\n")
    
    # Extract scores
    site_scores <- as.data.frame(scores(cca_result, display = "sites"))
    species_scores <- as.data.frame(scores(cca_result, display = "species"))
    env_scores <- as.data.frame(scores(cca_result, display = "bp"))
    
    # Add site information
    site_scores$site <- env_data$site
    site_scores <- site_scores %>%
      left_join(env_data, by = "site")
    
    # Create CCA plot
    cca_plot <- ggplot() +
      # Sites
      geom_point(data = site_scores, aes(x = CCA1, y = CCA2), 
                 size = 3, alpha = 0.7, color = "blue") +
      geom_text(data = site_scores, aes(x = CCA1, y = CCA2, label = site), 
                vjust = -1, size = 3) +
      # Environmental variables (biplot vectors)
      geom_segment(data = env_scores, 
                   aes(x = 0, y = 0, xend = CCA1 * 3, yend = CCA2 * 3),
                   arrow = arrow(length = unit(0.3, "cm")), 
                   color = "red", size = 1.2) +
      geom_text(data = env_scores, 
                aes(x = CCA1 * 3.5, y = CCA2 * 3.5, label = rownames(env_scores)),
                color = "red", size = 3, fontface = "bold") +
      labs(
        title = "Canonical Correspondence Analysis (CCA)",
        subtitle = "Detection counts vs Environmental variables",
        x = paste("CCA1 (", round(summary(cca_result)$cont$importance[2,1] * 100, 1), 
                  "%)", sep = ""),
        y = paste("CCA2 (", round(summary(cca_result)$cont$importance[2,2] * 100, 1), 
                  "%)", sep = "")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
    # Save results
    ggsave(file.path(output_dir, "cca_detections_environment.jpg"), 
           cca_plot, width = 12, height = 8)
    
    # Save summary
    capture.output(summary(cca_result), 
                   file = file.path(output_dir, "cca_summary.txt"))
    
    cat("[OK] Results saved to output directory\n")
    
    return(list(
      cca = cca_result,
      plot = cca_plot,
      site_scores = site_scores,
      species_scores = species_scores,
      env_scores = env_scores
    ))
    
  }, error = function(e) {
    cat(sprintf("[X] CCA failed: %s\n", e$message))
    return(NULL)
  })
}


# ==============================================================================
# SECTION 6: NONLINEAR RELATIONSHIPS (GAM ANALYSIS)
# ==============================================================================

#' Analyze nonlinear relationships using Generalized Additive Models (GAM)
#'
#' Tests for nonlinear relationships between species detection rates and
#' environmental variables using GAM with smooth terms. Particularly useful
#' for detecting unimodal, threshold, or complex response patterns.
#'
#' @param species_data List. Species detection data with raw_detections component
#' @param vegetation_data List. Vegetation data with summary component from
#'   read_vegetation_data()
#' @param data_type Character. Type of data for labeling ("katydid", "bird", etc.).
#'   Default: "katydid"
#' @param output_dir Character. Directory path for saving results and plots
#'
#' @return List with one element per significant environmental variable:
#'   \itemize{
#'     \item model: gam object from mgcv package
#'     \item p_value: p-value for smooth term significance
#'     \item plot: ggplot object showing GAM fit with confidence intervals
#'     \item deviance_explained: proportion of deviance explained by model
#'   }
#'   Returns empty list if no significant relationships found, NULL if insufficient data.
#'
#' @details
#' **Method overview:**
#' GAM (Generalized Additive Models) extends GLM by allowing smooth nonlinear
#' functions of predictors. This is ideal for ecological data where responses
#' to environmental gradients are often nonlinear (e.g., optimal ranges, saturation).
#'
#' **Model specification:**
#' - Response: total_detections (count data)
#' - Family: negative binomial (appropriate for overdispersed counts)
#' - Smooth: s(predictor, k=4) with 4 basis functions
#' - One univariate model per environmental variable
#'
#' **Data preparation:**
#' - Creates detection count matrix from raw detections
#' - Calculates species_richness and total_detections per site
#' - Merges with vegetation data using site-to-plot mapping
#' - Filters to complete cases only
#'
#' **Variables tested:**
#' - tree_species_richness: tree diversity
#' - tree_abundance: total tree count
#' - tree_total_basal_area: cumulative tree BA
#' - liana_species_richness: liana diversity
#' - liana_rooted_stems: liana stem count
#' - liana_total_basal_area: cumulative liana BA
#'
#' **Significance testing:**
#' Only relationships with p < 0.05 for the smooth term are retained and plotted.
#' The smooth term test evaluates whether the nonlinear pattern significantly
#' improves fit compared to a null (intercept-only) model.
#'
#' **Visualization:**
#' For significant relationships:
#' - Scatter plot: observed data points
#' - Red line: GAM predicted values
#' - Red ribbon: 95% confidence intervals
#' - Title includes p-value and significance note
#'
#' **Output files:**
#' - gam_[datatype]_[variable].jpg: individual GAM plot for each significant predictor
#' - gam_[datatype]_summary.csv: table of all significant results with p-values
#'   and deviance explained
#'
#' **Minimum requirements:**
#' - At least 5 sites with complete data (GAM needs sufficient df)
#' - raw_detections component in species_data
#'
#' @examples
#' \dontrun{
#' # Analyze nonlinear relationships for katydids
#' gam_results <- analyze_nonlinear_relationships(
#'   species_data = katydid_data,
#'   vegetation_data = forestgeo_vegetation,
#'   data_type = "katydid",
#'   output_dir = "results/gam/"
#' )
#' 
#' # Check which variables show nonlinear patterns
#' names(gam_results)  # Variables with significant smooth terms
#' 
#' # View deviance explained
#' sapply(gam_results, function(x) x$deviance_explained)
#' 
#' # Extract model for detailed inspection
#' tree_richness_gam <- gam_results$tree_species_richness$model
#' summary(tree_richness_gam)
#' gam.check(tree_richness_gam)  # Diagnostic plots
#' 
#' # Display plot
#' print(gam_results$tree_species_richness$plot)
#' }
#'
#' @references
#' Wood, S. N. (2017). Generalized Additive Models: An Introduction with R
#' (2nd ed.). Chapman and Hall/CRC.
#'
#' Hastie, T. J., & Tibshirani, R. J. (1990). Generalized Additive Models.
#' Chapman and Hall/CRC Monographs on Statistics and Applied Probability.
#'
#' Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M.
#' (2009). Mixed Effects Models and Extensions in Ecology with R. Springer.
#' (Chapter 3: GAM for count data)
#'
#' @seealso \code{\link{analyze_vegetation_effects}},
#'   \code{\link{analyze_all_glm_relationships}},
#'   \code{\link[mgcv]{gam}}, \code{\link[mgcv]{nb}}
#'
#' @export
analyze_nonlinear_relationships <- function(species_data, vegetation_data, 
                                            data_type = "katydid", 
                                            output_dir) {
  
  cat(sprintf("\n=== NON-LINEAR ANALYSIS: %s ===\n", toupper(data_type)))
  
  # Load mgcv package for GAM
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    cat(" mgcv package required but not installed\n")
    return(NULL)
  }
  
  # Calculate detection counts per site
  if (!is.null(species_data$raw_detections)) {
    detection_counts <- create_detection_count_matrix(species_data$raw_detections, 
                                                      "common_name", "site")
    
    # Calculate richness and total detections, then merge with vegetation
    species_richness_data <- detection_counts %>%
      mutate(
        species_richness = rowSums(select(., -site)),
        total_detections = rowSums(select(., -site)),
        site_number = as.numeric(gsub("S", "", site)),
        Plot = sprintf("YB-P%02d", site_number)
      ) %>%
      left_join(vegetation_data$summary, by = "Plot") %>%
      filter(complete.cases(.))
  } else {
    cat("[X] No raw detections available\n")
    return(NULL)
  }
  
  if (nrow(species_richness_data) < 5) {
    cat("[X] Not enough data for GAM analysis (need at least 5 sites)\n")
    return(NULL)
  }
  
  cat(sprintf("Sites available for GAM: %d\n", nrow(species_richness_data)))
  
  # Environmental variables to test
  env_vars <- c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area")
  
  gam_results <- list()
  
  for (env_var in env_vars) {
    if (env_var %in% colnames(species_richness_data)) {
      cat(sprintf("  Testing non-linear: %s detections vs %s\n", data_type, env_var))
      
      # GAM for total detections
      tryCatch({
        # GAM with smooth term (k=4 basis functions)
        # Build formula dynamically for variable column name
        gam_formula <- as.formula(paste0("total_detections ~ s(", env_var, ", k = 4)"))
        gam_model <- mgcv::gam(gam_formula, 
                               data = species_richness_data, 
                               family = mgcv::nb())  # Negative binomial for overdispersed counts
        
        # Test smooth term significance
        gam_summary <- summary(gam_model)
        p_value <- gam_summary$s.table[1, "p-value"]
        
        if (p_value < 0.05) {
          cat(sprintf("    [OK] Significant non-linear relationship (p = %.4f)\n", p_value))
          
          # Create predictions for plotting
          env_range <- seq(min(species_richness_data[[env_var]], na.rm = TRUE),
                           max(species_richness_data[[env_var]], na.rm = TRUE),
                           length.out = 100)
          
          pred_data <- data.frame(x = env_range)
          names(pred_data) <- env_var
          predictions <- predict(gam_model, newdata = pred_data, se.fit = TRUE)
          
          plot_data <- data.frame(
            x = env_range,
            fitted = predictions$fit,
            se = predictions$se.fit,
            upper = predictions$fit + 1.96 * predictions$se.fit,
            lower = predictions$fit - 1.96 * predictions$se.fit
          )
          
          # GAM plot with confidence intervals
          gam_plot <- ggplot() +
            geom_point(data = species_richness_data, 
                       aes(x = .data[[env_var]], y = total_detections), 
                       alpha = 0.7, size = 3, color = "darkblue") +
            geom_line(data = plot_data, aes(x = x, y = fitted), 
                      color = "red", size = 1.5) +
            geom_ribbon(data = plot_data, aes(x = x, ymin = lower, ymax = upper), 
                        alpha = 0.3, fill = "red") +
            labs(
              title = paste("Non-linear GAM:", data_type, "detections vs", env_var),
              subtitle = paste("p =", round(p_value, 4), "| Smooth term significant |",
                               "Deviance explained:", round(gam_summary$dev.expl * 100, 1), "%"),
              x = gsub("_", " ", stringr::str_to_title(env_var)),
              y = paste(stringr::str_to_title(data_type), "Total Detections")
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 14, face = "bold", color = "darkgreen"),
              plot.subtitle = element_text(size = 11),
              axis.title = element_text(size = 12)
            )
          
          # Save plot
          plot_filename <- paste0("gam_", data_type, "_", 
                                  gsub("[^A-Za-z0-9]", "_", env_var), ".jpg")
          ggsave(file.path(output_dir, plot_filename), gam_plot, width = 10, height = 7)
          
          # Store results
          gam_results[[env_var]] <- list(
            model = gam_model,
            p_value = p_value,
            plot = gam_plot,
            deviance_explained = gam_summary$dev.expl
          )
          
        } else {
          cat(sprintf("    Non-significant (p = %.4f)\n", p_value))
        }
        
      }, error = function(e) {
        cat(sprintf("     GAM failed for %s: %s\n", env_var, e$message))
      })
    }
  }
  
  # Create summary table of GAM results
  if (length(gam_results) > 0) {
    gam_summary_df <- data.frame(
      variable = names(gam_results),
      p_value = sapply(gam_results, function(x) x$p_value),
      deviance_explained = sapply(gam_results, function(x) x$deviance_explained * 100),
      stringsAsFactors = FALSE
    ) %>%
      arrange(p_value)
    
    # Save summary
    write.csv(gam_summary_df, 
              file.path(output_dir, paste0("gam_", data_type, "_summary.csv")), 
              row.names = FALSE)
    
    cat(sprintf("  [OK] Found %d significant non-linear relationships\n", nrow(gam_summary_df)))
    cat(sprintf("  [OK] Results saved to %s\n", output_dir))
  } else {
    cat("  No significant non-linear relationships detected\n")
  }
  
  return(gam_results)
}


# ==============================================================================
# SECTION 7: COMPREHENSIVE GLM ANALYSIS (MULTI-VARIABLE)
# ==============================================================================

#' Comprehensive GLM analysis of all ecological relationships
#'
#' Performs systematic GLM analyses to test relationships between species richness
#' (katydids, birds) and environmental variables, plus cross-taxa relationships.
#' Generates plots and summary tables for all significant relationships.
#'
#' @param all_results List. Complete analysis results containing:
#'   \itemize{
#'     \item katydid_data: with presence_matrix component
#'     \item bird_data: with presence_matrix component
#'     \item vegetation_data: with summary component
#'   }
#' @param output_dir Character. Directory path for saving GLM results
#'
#' @return Named list of significant GLM results (p < 0.05), where each element contains:
#'   \itemize{
#'     \item model: glm object
#'     \item p_value: model significance p-value
#'     \item aic: Akaike Information Criterion
#'     \item response: name of response variable
#'     \item predictor: name of predictor variable
#'     \item sample_size: number of observations used
#'   }
#'   Returns empty list if no significant relationships found.
#'
#' @details
#' **Analysis workflow:**
#' This comprehensive function tests three categories of relationships:
#' 
#' 1. **Environmental effects on katydids:**
#'    - Response: katydid species richness
#'    - Predictors: 6 vegetation variables (tree/liana richness, abundance, basal area)
#' 
#' 2. **Environmental effects on birds:**
#'    - Response: bird species richness
#'    - Predictors: same 6 vegetation variables
#' 
#' 3. **Cross-taxa relationships:**
#'    - Bird richness ~ katydid richness (do sites with more katydids have more birds?)
#'    - Katydid richness ~ bird richness (inverse relationship)
#'
#' **Data preparation:**
#' - Extracts species richness from presence/absence matrices
#' - Converts site codes (S01-S12) to ForestGEO plot IDs (YB-P01 to YB-P12)
#' - Merges species and vegetation data by plot ID
#' - Filters to complete cases (removes any sites with missing env data)
#' - Uses inner joins for cross-taxa analyses (common sites only)
#'
#' **GLM specifications:**
#' Each relationship tested with separate univariate GLM:
#' - Family: Poisson (for count data - species richness)
#' - Model: richness ~ predictor
#' - Significance threshold: p < 0.05
#' - Model selection: AIC values reported for comparison
#'
#' **Statistical considerations:**
#' - Checks for sufficient sample size (n >= 5 for cross-taxa)
#' - Uses robust GLM wrapper that handles model failures gracefully
#' - Saves only significant relationships (p < 0.05)
#' - Small sample sizes (12 katydid sites) may limit power
#'
#' **Output files:**
#' Created in GLM_analyses/ subdirectory:
#' - Individual plots for each significant relationship (via perform_glm_analysis_robust)
#' - GLM_significant_relationships_summary.csv: table of all significant results
#' - Includes relationship name, response, predictor, p-value, AIC, and sample size
#'
#' **Interpretation guidance:**
#' No significant results may indicate:
#' - Genuine lack of linear relationships
#' - Small sample size (especially for katydids: n=12)
#' - High natural variability in tropical ecosystems
#' - Nonlinear patterns better captured by GAM (see analyze_nonlinear_relationships)
#' - Need for multivariate approaches (CCA, NMDS)
#'
#' @examples
#' \dontrun{
#' # Run comprehensive GLM analysis
#' glm_results <- analyze_all_glm_relationships(
#'   all_results = integrated_results,
#'   output_dir = "results/"
#' )
#' 
#' # Check number of significant relationships
#' length(glm_results)
#' 
#' # View summary of results
#' names(glm_results)  # Names of significant relationships
#' 
#' # Extract specific model
#' tree_rich_model <- glm_results$Katydid_tree_species_richness$model
#' summary(tree_rich_model)
#' 
#' # Get AIC values for model comparison
#' sapply(glm_results, function(x) x$aic)
#' 
#' # Check for cross-taxa relationships
#' if ("Bird_vs_Katydid_Richness" %in% names(glm_results)) {
#'   cross_taxa <- glm_results$Bird_vs_Katydid_Richness
#'   cat("Birds and katydids show correlated richness patterns\n")
#' }
#' }
#'
#' @references
#' McCullagh, P., & Nelder, J. A. (1989). Generalized Linear Models
#' (2nd ed.). Chapman and Hall/CRC.
#'
#' Zuur, A. F., Ieno, E. N., & Elphick, C. S. (2010). A protocol for data
#' exploration to avoid common statistical problems. Methods in Ecology and
#' Evolution, 1(1), 3-14.
#'
#' Burnham, K. P., & Anderson, D. R. (2002). Model Selection and Multimodel
#' Inference: A Practical Information-Theoretic Approach (2nd ed.). Springer.
#'
#' @seealso \code{\link{perform_glm_analysis_robust}},
#'   \code{\link{analyze_nonlinear_relationships}},
#'   \code{\link{analyze_vegetation_effects}}
#'
#' @export
analyze_all_glm_relationships <- function(all_results, output_dir) {
  
  cat("\n=== STARTING COMPREHENSIVE GLM ANALYSIS ===\n")
  
  # Create GLM output directory
  glm_output_dir <- file.path(output_dir, "GLM_analyses")
  dir.create(glm_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  significant_results <- list()
  
  # ============================================================================
  # 1. ENVIRONMENTAL VARIABLES vs KATYDID RICHNESS
  # ============================================================================
  cat("\n1. ENVIRONMENTAL VARIABLES vs KATYDID RICHNESS\n")
  
  if (!is.null(all_results$katydid_data) && !is.null(all_results$vegetation_data)) {
    
    tryCatch({
      # Prepare katydid data
      katydid_matrix <- all_results$katydid_data$presence_matrix
      
      cat(sprintf("Katydid matrix: %d sites  %d columns\n", 
                  nrow(katydid_matrix), ncol(katydid_matrix)))
      
      # Identify species columns (all except 'site')
      species_columns <- colnames(katydid_matrix)[colnames(katydid_matrix) != "site"]
      
      cat(sprintf("  Found %d katydid species columns\n", length(species_columns)))
      
      # Calculate katydid richness
      katydid_species_data <- katydid_matrix[, species_columns, drop = FALSE]
      katydid_richness <- rowSums(katydid_species_data, na.rm = TRUE)
      
      # Create dataframe with site, richness, and Plot conversion
      katydid_richness_df <- data.frame(
        site = katydid_matrix$site,
        katydid_richness = katydid_richness,
        stringsAsFactors = FALSE
      )
      
      # Convert site to Plot (S01 -> YB-P01)
      katydid_richness_df$site_number <- as.numeric(gsub("S", "", katydid_richness_df$site))
      katydid_richness_df$Plot <- sprintf("YB-P%02d", katydid_richness_df$site_number)
      
      cat(sprintf("  Katydid richness: min=%d, max=%d, mean=%.1f\n", 
                  min(katydid_richness_df$katydid_richness), 
                  max(katydid_richness_df$katydid_richness),
                  mean(katydid_richness_df$katydid_richness)))
      
      # Merge with vegetation data
      katydid_env_data <- merge(
        katydid_richness_df[, c("site", "Plot", "katydid_richness")], 
        all_results$vegetation_data$summary, 
        by = "Plot",
        all.x = TRUE
      )
      
      # Remove rows with missing vegetation data
      katydid_env_data <- katydid_env_data[complete.cases(katydid_env_data), ]
      
      cat(sprintf("  Merged katydid-vegetation data: %d rows\n", nrow(katydid_env_data)))
      
      if (nrow(katydid_env_data) > 0) {
        # Environmental variables to test
        env_vars <- c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                      "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area")
        
        # Test each environmental variable
        for (env_var in env_vars) {
          if (env_var %in% colnames(katydid_env_data)) {
            cat(sprintf("    Testing: Katydid richness vs %s\n", env_var))
            
            result <- perform_glm_analysis_robust(
              response_data = katydid_env_data$katydid_richness,
              predictor_data = katydid_env_data[[env_var]],
              response_name = "Katydid Species Richness",
              predictor_name = env_var,
              output_dir = glm_output_dir
            )
            
            if (!is.null(result)) {
              significant_results[[paste("Katydid", env_var, sep = "_")]] <- result
            }
          } else {
            cat(sprintf("[WARNING]     Variable %s not found in vegetation data\n", env_var))
          }
        }
      } else {
        cat("   No katydid-vegetation data after merge\n")
      }
      
    }, error = function(e) {
      cat(sprintf("   Error in katydid analysis: %s\n", e$message))
    })
  }
  
  # ============================================================================
  # 2. ENVIRONMENTAL VARIABLES vs BIRD RICHNESS
  # ============================================================================
  cat("\n2. ENVIRONMENTAL VARIABLES vs BIRD RICHNESS\n")
  
  if (!is.null(all_results$bird_data) && !is.null(all_results$vegetation_data)) {
    
    tryCatch({
      # Prepare bird data
      bird_matrix <- all_results$bird_data$presence_matrix
      
      cat(sprintf("Bird matrix: %d sites  %d columns\n", 
                  nrow(bird_matrix), ncol(bird_matrix)))
      
      # Identify bird species columns
      bird_species_columns <- colnames(bird_matrix)[colnames(bird_matrix) != "site"]
      
      cat(sprintf("  Found %d bird species columns\n", length(bird_species_columns)))
      
      # Calculate bird richness
      bird_species_data <- bird_matrix[, bird_species_columns, drop = FALSE]
      bird_richness <- rowSums(bird_species_data, na.rm = TRUE)
      
      # Create dataframe with site, richness, and Plot conversion
      bird_richness_df <- data.frame(
        site = bird_matrix$site,
        bird_richness = bird_richness,
        stringsAsFactors = FALSE
      )
      
      # Convert site to Plot
      bird_richness_df$site_number <- as.numeric(gsub("S", "", bird_richness_df$site))
      bird_richness_df$Plot <- sprintf("YB-P%02d", bird_richness_df$site_number)
      
      cat(sprintf("  Bird richness: min=%d, max=%d, mean=%.1f\n", 
                  min(bird_richness_df$bird_richness), 
                  max(bird_richness_df$bird_richness),
                  mean(bird_richness_df$bird_richness)))
      
      # Merge with vegetation data
      bird_env_data <- merge(
        bird_richness_df[, c("site", "Plot", "bird_richness")], 
        all_results$vegetation_data$summary, 
        by = "Plot",
        all.x = TRUE
      )
      
      # Remove rows with missing data
      bird_env_data <- bird_env_data[complete.cases(bird_env_data), ]
      
      cat(sprintf("  Merged bird-vegetation data: %d rows\n", nrow(bird_env_data)))
      
      if (nrow(bird_env_data) > 0) {
        # Environmental variables to test
        env_vars <- c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                      "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area")
        
        # Test each environmental variable for birds
        for (env_var in env_vars) {
          if (env_var %in% colnames(bird_env_data)) {
            cat(sprintf("    Testing: Bird richness vs %s\n", env_var))
            
            result <- perform_glm_analysis_robust(
              response_data = bird_env_data$bird_richness,
              predictor_data = bird_env_data[[env_var]],
              response_name = "Bird Species Richness",
              predictor_name = env_var,
              output_dir = glm_output_dir
            )
            
            if (!is.null(result)) {
              significant_results[[paste("Bird", env_var, sep = "_")]] <- result
            }
          }
        }
      } else {
        cat("   No bird-vegetation data after merge\n")
      }
      
    }, error = function(e) {
      cat(sprintf("   Error in bird analysis: %s\n", e$message))
    })
  }
  
  # ============================================================================
  # 3. KATYDID vs BIRD SPECIES RELATIONSHIPS
  # ============================================================================
  cat("\n3. KATYDID vs BIRD SPECIES RELATIONSHIPS\n")
  
  # Check if we have both datasets
  has_katydid_env <- exists("katydid_env_data") && !is.null(katydid_env_data) && 
    nrow(katydid_env_data) > 0
  has_bird_env <- exists("bird_env_data") && !is.null(bird_env_data) && 
    nrow(bird_env_data) > 0
  
  if (has_katydid_env && has_bird_env) {
    
    tryCatch({
      # Combine richness data using Plot as key
      # Only common sites (S01-S12) will be kept
      combined_richness <- merge(
        katydid_env_data[, c("Plot", "site", "katydid_richness")],
        bird_env_data[, c("Plot", "bird_richness")],
        by = "Plot",
        all = FALSE  # Inner join - common sites only
      )
      
      cat(sprintf("  Combined richness data: %d common sites\n", nrow(combined_richness)))
      
      if (nrow(combined_richness) >= 5) {
        
        # GLM: Bird richness vs katydid richness
        result1 <- perform_glm_analysis_robust(
          response_data = combined_richness$bird_richness,
          predictor_data = combined_richness$katydid_richness,
          response_name = "Bird Species Richness",
          predictor_name = "Katydid Species Richness",
          output_dir = glm_output_dir
        )
        
        if (!is.null(result1)) {
          significant_results[["Bird_vs_Katydid_Richness"]] <- result1
        }
        
        # GLM inverse: Katydids vs birds
        result2 <- perform_glm_analysis_robust(
          response_data = combined_richness$katydid_richness,
          predictor_data = combined_richness$bird_richness,
          response_name = "Katydid Species Richness",
          predictor_name = "Bird Species Richness",
          output_dir = glm_output_dir
        )
        
        if (!is.null(result2)) {
          significant_results[["Katydid_vs_Bird_Richness"]] <- result2
        }
        
      } else {
        cat(sprintf("   Not enough combined data for inter-group analysis (n=%d)\n", 
                    nrow(combined_richness)))
      }
      
    }, error = function(e) {
      cat(sprintf("   Error in katydid-bird analysis: %s\n", e$message))
    })
  } else {
    cat("   Missing data for inter-group analysis\n")
    cat(sprintf("    Katydid-env data available: %s\n", has_katydid_env))
    cat(sprintf("    Bird-env data available: %s\n", has_bird_env))
  }
  
  # ============================================================================
  # 4. RESULTS SUMMARY
  # ============================================================================
  cat("\n=== ANALYSIS SUMMARY ===\n")
  
  if (length(significant_results) > 0) {
    
    cat(sprintf("[OK] FOUND %d SIGNIFICANT RELATIONSHIPS (p < 0.05)\n", 
                length(significant_results)))
    
    # Create summary table
    glm_summary <- data.frame(
      Relationship = character(),
      Response = character(),
      Predictor = character(),
      P_value = numeric(),
      AIC = numeric(),
      Sample_size = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (name in names(significant_results)) {
      result <- significant_results[[name]]
      glm_summary <- rbind(glm_summary, data.frame(
        Relationship = name,
        Response = result$response,
        Predictor = result$predictor,
        P_value = result$p_value,
        AIC = result$aic,
        Sample_size = result$sample_size,
        stringsAsFactors = FALSE
      ))
    }
    
    # Sort by p-value
    glm_summary <- glm_summary %>%
      arrange(P_value)
    
    # Save summary
    write.csv(glm_summary, 
              file.path(glm_output_dir, "GLM_significant_relationships_summary.csv"), 
              row.names = FALSE)
    
    # Display results
    cat("\nSIGNIFICANT RELATIONSHIPS:\n")
    for (i in 1:nrow(glm_summary)) {
      cat(sprintf("  %d. %s -> %s (p = %.4f, n = %d)\n", 
                  i, glm_summary$Predictor[i], glm_summary$Response[i], 
                  glm_summary$P_value[i], glm_summary$Sample_size[i]))
    }
    
  } else {
    cat("[X] NO SIGNIFICANT RELATIONSHIPS FOUND (p < 0.05)\n")
    cat("\nPossible explanations:\n")
    cat("   Real lack of linear relationships in this ecosystem\n")
    cat("   Small sample size (12 sites for katydids)\n")
    cat("   High natural variability in tropical forests\n")
    cat("   Nonlinear patterns better captured by GAM\n")
    cat("   Consider multivariate approaches (CCA, NMDS)\n")
  }
  
  cat(sprintf("\n[OK] Results saved in: %s\n", glm_output_dir))
  cat("[COMPLETE] GLM ANALYSIS COMPLETE!\n")
  
  return(significant_results)
}


# ==============================================================================
# SECTION 8: HELPER FUNCTIONS FOR GLM ANALYSIS
# ==============================================================================

#' Perform Robust GLM Analysis with Automatic Family Selection
#'
#' @description
#' Conducts a robust Generalized Linear Model (GLM) analysis between a response 
#' and predictor variable with automatic selection of the appropriate error 
#' distribution family. The function handles data validation, overdispersion 
#' detection, and generates diagnostic plots and detailed outputs for 
#' significant relationships.
#'
#' @details
#' This function implements a comprehensive GLM workflow:
#' 
#' **Data Validation:**
#' - Removes missing values (NA) and infinite values
#' - Checks for sufficient sample size (minimum n=4)
#' - Verifies variance in both response and predictor variables
#' - Reports data summaries and diagnostics
#' 
#' **Family Selection Logic:**
#' - **Count data** (non-negative integers): Initially fits Poisson GLM
#' - **Overdispersion check**: If dispersion ratio > 1.5, switches to quasi-Poisson
#' - **Continuous data**: Uses Gaussian family
#' 
#' **Overdispersion Detection:**
#' Calculates dispersion ratio as residual deviance / residual degrees of freedom.
#' Ratios > 1.5 indicate overdispersion requiring quasi-Poisson correction
#' (Ver Hoef & Boveng, 2007).
#' 
#' **Significance Testing:**
#' Uses Wald test (z-test for Poisson, t-test for Gaussian) to assess predictor
#' significance at  = 0.05 level.
#' 
#' **Output Generation (for significant relationships only):**
#' - Scatter plot with fitted GLM curve and 95% confidence interval
#' - Detailed model summary file (.txt) including coefficients, diagnostics,
#'   and ANOVA table
#' - Both files saved to specified output directory
#'
#' @param response_data Numeric vector. Response variable values
#' @param predictor_data Numeric vector. Predictor variable values (same length as response)
#' @param response_name Character string. Name of response variable for labeling
#' @param predictor_name Character string. Name of predictor variable for labeling
#' @param output_dir Character string. Path to directory for saving plots and model summaries
#'
#' @return 
#' Returns a list with the following elements if relationship is significant (p < 0.05),
#' NULL otherwise:
#' \describe{
#'   \item{response}{Character. Name of response variable}
#'   \item{predictor}{Character. Name of predictor variable}
#'   \item{p_value}{Numeric. P-value from Wald test}
#'   \item{aic}{Numeric. Akaike Information Criterion (NA for quasi-families)}
#'   \item{coefficient}{Numeric. Slope coefficient of predictor}
#'   \item{direction}{Character. "positive" or "negative" relationship}
#'   \item{model}{glm object. Fitted GLM model}
#'   \item{plot_path}{Character. Full path to saved plot file}
#'   \item{model_path}{Character. Full path to saved model summary file}
#'   \item{sample_size}{Integer. Number of observations used in model}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Count data (species richness vs vegetation)
#' richness <- c(12, 15, 8, 22, 18, 14, 10, 20)
#' canopy_cover <- c(75, 82, 60, 90, 88, 78, 65, 85)
#' 
#' result <- perform_glm_analysis_robust(
#'   response_data = richness,
#'   predictor_data = canopy_cover,
#'   response_name = "Species Richness",
#'   predictor_name = "Canopy Cover (%)",
#'   output_dir = "./glm_outputs"
#' )
#' 
#' if (!is.null(result)) {
#'   cat(sprintf("Significant %s relationship (p = %.4f)\n", 
#'               result$direction, result$p_value))
#' }
#' 
#' # Example 2: Continuous data (detection rate vs temperature)
#' detection_rate <- c(0.45, 0.62, 0.38, 0.71, 0.55, 0.49)
#' temperature <- c(24.5, 26.2, 23.1, 27.5, 25.3, 24.8)
#' 
#' result <- perform_glm_analysis_robust(
#'   response_data = detection_rate,
#'   predictor_data = temperature,
#'   response_name = "Detection Rate",
#'   predictor_name = "Temperature (degC)",
#'   output_dir = "./glm_outputs"
#' )
#' }
#'
#' @references
#' McCullagh, P., & Nelder, J. A. (1989). *Generalized Linear Models* (2nd ed.). 
#' Chapman and Hall/CRC. https://doi.org/10.1007/978-1-4899-3242-6
#' 
#' Ver Hoef, J. M., & Boveng, P. L. (2007). Quasi-Poisson vs. negative binomial 
#' regression: How should we model overdispersed count data? *Ecology*, 88(11), 
#' 2766-2772. https://doi.org/10.1890/07-0043.1
#' 
#' Zuur, A. F., Ieno, E. N., Walker, N. J., Saveliev, A. A., & Smith, G. M. (2009). 
#' *Mixed Effects Models and Extensions in Ecology with R*. Springer. 
#' https://doi.org/10.1007/978-0-387-87458-6
#' 
#' Dobson, A. J., & Barnett, A. G. (2018). *An Introduction to Generalized Linear 
#' Models* (4th ed.). CRC Press. https://doi.org/10.1201/9781315182780
#' 
#' Hilbe, J. M. (2011). *Negative Binomial Regression* (2nd ed.). Cambridge 
#' University Press. https://doi.org/10.1017/CBO9780511973420
#'
#' @seealso 
#' \code{\link[stats]{glm}} for GLM fitting
#' \code{\link[stats]{family}} for error distribution families
#' \code{\link{analyze_all_glm_relationships}} which uses this function
#'
#' @export
#' 
#' @importFrom stats glm poisson quasipoisson gaussian deviance df.residual AIC anova
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_minimal theme element_text ggsave
perform_glm_analysis_robust <- function(response_data, predictor_data, response_name, predictor_name, output_dir) {
  
  cat(sprintf("    Analyzing: %s vs %s...\n", response_name, predictor_name))
  
  # Create dataframe for analysis with data validation
  data_for_glm <- data.frame(
    response = as.numeric(response_data),
    predictor = as.numeric(predictor_data)
  )
  
  # Remove missing or infinite values
  data_for_glm <- data_for_glm[complete.cases(data_for_glm), ]
  data_for_glm <- data_for_glm[is.finite(data_for_glm$response) & is.finite(data_for_glm$predictor), ]
  
  cat(sprintf("    Sample size after cleaning: %d\n", nrow(data_for_glm)))
  
  # Check for sufficient sample size
  if (nrow(data_for_glm) < 4) {
    cat(sprintf(" Insufficient data (n=%d)\n", nrow(data_for_glm)))
    return(NULL)
  }
  
  # Check variance (avoid constant variables)
  if (var(data_for_glm$response) < 1e-10) {
    cat(sprintf(" No variance in response variable\n"))
    return(NULL)
  }
  
  if (var(data_for_glm$predictor) < 1e-10) {
    cat(sprintf(" No variance in predictor variable\n"))
    return(NULL)
  }
  
  # Display data summary
  cat(sprintf("    Response range: %.2f - %.2f (mean: %.2f, sd: %.2f)\n", 
              min(data_for_glm$response), max(data_for_glm$response), 
              mean(data_for_glm$response), sd(data_for_glm$response)))
  cat(sprintf("    Predictor range: %.2f - %.2f (mean: %.2f, sd: %.2f)\n", 
              min(data_for_glm$predictor), max(data_for_glm$predictor), 
              mean(data_for_glm$predictor), sd(data_for_glm$predictor)))
  
  # Execute GLM with error handling
  tryCatch({
    
    # Intelligently determine distribution family
    is_count_data <- all(data_for_glm$response >= 0) && all(data_for_glm$response == round(data_for_glm$response))
    
    if (is_count_data && max(data_for_glm$response) > 2) {
      # Count data -> try Poisson
      model <- glm(response ~ predictor, data = data_for_glm, family = poisson())
      family_used <- "poisson"
      
      # Check for overdispersion
      residual_deviance <- deviance(model)
      df_residual <- df.residual(model)
      dispersion_ratio <- residual_deviance / df_residual
      
      if (dispersion_ratio > 1.5) {
        # Use quasi-poisson if overdispersed
        model <- glm(response ~ predictor, data = data_for_glm, family = quasipoisson())
        family_used <- "quasipoisson"
        cat(sprintf("    Using quasi-poisson due to overdispersion (ratio: %.2f)\n", dispersion_ratio))
      }
    } else {
      # Continuous data or low counts -> Gaussian
      model <- glm(response ~ predictor, data = data_for_glm, family = gaussian())
      family_used <- "gaussian"
    }
    
    # Extract results
    summary_model <- summary(model)
    
    # Check that we have a predictor
    if (nrow(summary_model$coefficients) < 2) {
      cat(sprintf("[X] Model failed to fit predictor\n"))
      return(NULL)
    }
    
    p_value <- summary_model$coefficients[2, 4]  # P-value of predictor
    coefficient <- summary_model$coefficients[2, 1]  # Coefficient of predictor
    
    # AIC (handle case where AIC not available for quasi-families)
    if (family_used == "quasipoisson") {
      aic_value <- NA
    } else {
      aic_value <- AIC(model)
    }
    
    cat(sprintf("    Model: %s, Coefficient: %.4f, P-value: %.4f\n", 
                family_used, coefficient, p_value))
    
    # Test for significance
    if (p_value < 0.05) {
      cat(sprintf("[OK] SIGNIFICANT relationship found!\n"))
      
      # Create plot
      plot_filename <- paste0(gsub("[^A-Za-z0-9]", "_", paste(response_name, "vs", predictor_name)), ".jpg")
      plot_path <- file.path(output_dir, plot_filename)
      
      # Determine relationship direction
      direction <- ifelse(coefficient > 0, "positive", "negative")
      
      p <- ggplot(data_for_glm, aes(x = predictor, y = response)) +
        geom_point(size = 3, alpha = 0.8, color = "darkblue") +
        geom_smooth(method = "glm", method.args = list(family = family_used), 
                    se = TRUE, color = "red", fill = "lightcoral", alpha = 0.3) +
        labs(
          title = paste("Significant GLM:", response_name, "vs", predictor_name),
          subtitle = paste0("p = ", round(p_value, 4), " (", direction, " relationship), Family = ", family_used, 
                            ifelse(is.na(aic_value), "", paste0(", AIC = ", round(aic_value, 2)))),
          x = predictor_name,
          y = response_name
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold", color = "darkgreen"),
          plot.subtitle = element_text(size = 11),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)
        )
      
      ggsave(plot_path, p, width = 10, height = 7, dpi = 300)
      cat(sprintf("    Plot saved: %s\n", basename(plot_path)))
      
      # Save detailed model summary
      model_filename <- paste0(gsub("[^A-Za-z0-9]", "_", paste(response_name, "vs", predictor_name)), "_model.txt")
      model_path <- file.path(output_dir, model_filename)
      
      capture.output(
        {
          cat("SIGNIFICANT GLM ANALYSIS RESULTS\n")
          cat("===============================\n")
          cat(sprintf("Response: %s\n", response_name))
          cat(sprintf("Predictor: %s\n", predictor_name))
          cat(sprintf("Sample size: %d\n", nrow(data_for_glm)))
          cat(sprintf("Family: %s\n", family_used))
          cat(sprintf("P-value: %.6f\n", p_value))
          cat(sprintf("Coefficient: %.6f (%s relationship)\n", coefficient, direction))
          if (!is.na(aic_value)) cat(sprintf("AIC: %.2f\n", aic_value))
          cat(sprintf("R-squared (pseudo): %.4f\n", 1 - (deviance(model) / summary_model$null.deviance)))
          cat("\nData summary:\n")
          cat(sprintf("Response: min=%.2f, max=%.2f, mean=%.2f, sd=%.2f\n", 
                      min(data_for_glm$response), max(data_for_glm$response), 
                      mean(data_for_glm$response), sd(data_for_glm$response)))
          cat(sprintf("Predictor: min=%.2f, max=%.2f, mean=%.2f, sd=%.2f\n", 
                      min(data_for_glm$predictor), max(data_for_glm$predictor), 
                      mean(data_for_glm$predictor), sd(data_for_glm$predictor)))
          cat("\nDetailed Model Summary:\n")
          print(summary_model)
          cat("\nANOVA:\n")
          print(anova(model, test = "Chi"))
        },
        file = model_path
      )
      
      # Return structured results list
      return(list(
        response = response_name,
        predictor = predictor_name,
        p_value = p_value,
        aic = ifelse(is.na(aic_value), 0, aic_value),
        coefficient = coefficient,
        direction = direction,
        model = model,
        plot_path = plot_path,
        model_path = model_path,
        sample_size = nrow(data_for_glm)
      ))
      
    } else {
      cat(sprintf(" Not significant (p = %.4f)\n", p_value))
      return(NULL)
    }
    
  }, error = function(e) {
    cat(sprintf("[X] GLM failed: %s\n", e$message))
    return(NULL)
  })
}


# ==============================================================================
# SECTION 9: WRAPPER FUNCTIONS FOR MULTIVARIATE ANALYSES
# ==============================================================================

#' Perform Mantel Test for Matrix Correlation
#'
#' @description
#' Conducts a Mantel test to assess the correlation between two distance or 
#' dissimilarity matrices. The Mantel test uses permutations to test the null 
#' hypothesis that there is no relationship between the two matrices.
#'
#' @details
#' The Mantel test (Mantel, 1967) is a permutation-based procedure for testing 
#' the correlation between two distance matrices. It is commonly used in ecology 
#' to test relationships between community composition and environmental distances, 
#' or to compare community patterns detected by different methods.
#' 
#' **Test Procedure:**
#' 1. Calculates correlation (Pearson, Spearman, or Kendall) between matrices
#' 2. Permutes rows/columns of one matrix many times (default 999)
#' 3. Recalculates correlation for each permutation
#' 4. P-value = proportion of permuted correlations  observed correlation
#' 
#' **Method Selection:**
#' - **Spearman** (default): Robust to outliers, tests monotonic relationships
#' - **Pearson**: Assumes linear relationship, sensitive to outliers
#' - **Kendall**: Most robust but computationally intensive
#' 
#' The function includes error handling for cases where matrices have 
#' incompatible dimensions or insufficient variation.
#'
#' @param matrix1 Distance matrix (dist object or symmetric matrix). 
#'   First dissimilarity/distance matrix
#' @param matrix2 Distance matrix (dist object or symmetric matrix). 
#'   Second dissimilarity/distance matrix (must have same dimensions as matrix1)
#' @param method Character string. Correlation method: "pearson", "spearman" 
#'   (default), or "kendall"
#' @param permutations Integer. Number of permutations for significance test. 
#'   Default = 999. Higher values give more precise p-values but take longer
#' @param na.rm Logical. Should missing values be removed? Default = TRUE
#'
#' @return 
#' Returns a list with the following elements, or NULL if test fails:
#' \describe{
#'   \item{statistic}{Numeric. Mantel correlation coefficient (r)}
#'   \item{signif}{Numeric. P-value based on permutation test}
#'   \item{perm}{Integer. Number of permutations performed}
#'   \item{method}{Character. Correlation method used}
#'   \item{call}{Language. The function call}
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Compare acoustic and metabarcoding community matrices
#' library(vegan)
#' 
#' # Create distance matrices
#' acoustic_dist <- vegdist(acoustic_matrix, method = "jaccard")
#' metabar_dist <- vegdist(metabar_matrix, method = "jaccard")
#' 
#' # Perform Mantel test
#' mantel_result <- perform_mantel_test(
#'   matrix1 = acoustic_dist,
#'   matrix2 = metabar_dist,
#'   method = "spearman",
#'   permutations = 999
#' )
#' 
#' if (!is.null(mantel_result)) {
#'   cat(sprintf("Mantel r = %.3f, p = %.3f\n", 
#'               mantel_result$statistic, mantel_result$signif))
#' }
#' 
#' # Example 2: Test environmental vs community distance
#' env_dist <- dist(scale(environmental_data))
#' comm_dist <- vegdist(community_matrix, method = "bray")
#' 
#' mantel_result <- perform_mantel_test(env_dist, comm_dist, 
#'                                       method = "pearson",
#'                                       permutations = 9999)
#' }
#'
#' @references
#' Mantel, N. (1967). The detection of disease clustering and a generalized 
#' regression approach. *Cancer Research*, 27(2), 209-220.
#' 
#' Legendre, P., & Legendre, L. (2012). *Numerical Ecology* (3rd ed.). 
#' Elsevier. https://doi.org/10.1016/B978-0-444-53868-0.50018-3
#' 
#' Oksanen, J., et al. (2022). vegan: Community Ecology Package. R package 
#' version 2.6-4. https://CRAN.R-project.org/package=vegan
#' 
#' Goslee, S. C., & Urban, D. L. (2007). The ecodist package for dissimilarity-based 
#' analysis of ecological data. *Journal of Statistical Software*, 22(7), 1-19. 
#' https://doi.org/10.18637/jss.v022.i07
#'
#' @seealso 
#' \code{\link[vegan]{mantel}} for the underlying implementation
#' \code{\link[vegan]{vegdist}} for creating distance matrices
#' \code{\link{analyze_site_similarity}} which uses this function
#'
#' @export
#' 
#' @importFrom vegan mantel
perform_mantel_test <- function(matrix1, matrix2, method = "spearman", 
                                permutations = 999, na.rm = TRUE) {
  
  # Load required package
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required for Mantel test. Please install it.")
  }
  
  cat(sprintf("\n[>>] Performing Mantel test (method: %s, permutations: %d)...\n", 
              method, permutations))
  
  # Validate inputs
  if (is.null(matrix1) || is.null(matrix2)) {
    cat(" Error: One or both matrices are NULL\n")
    return(NULL)
  }
  
  # Convert to dist objects if needed
  if (!inherits(matrix1, "dist")) {
    matrix1 <- as.dist(matrix1)
  }
  if (!inherits(matrix2, "dist")) {
    matrix2 <- as.dist(matrix2)
  }
  
  # Check dimensions match
  if (attr(matrix1, "Size") != attr(matrix2, "Size")) {
    cat(sprintf("[X] Error: Matrix dimensions don't match (%d vs %d)\n", 
                attr(matrix1, "Size"), attr(matrix2, "Size")))
    return(NULL)
  }
  
  # Check for sufficient variation
  if (length(unique(as.vector(matrix1))) < 2 || length(unique(as.vector(matrix2))) < 2) {
    cat(" Error: Insufficient variation in one or both matrices\n")
    return(NULL)
  }
  
  # Perform Mantel test with error handling
  tryCatch({
    result <- vegan::mantel(matrix1, matrix2, 
                            method = method, 
                            permutations = permutations,
                            na.rm = na.rm)
    
    cat(sprintf("[OK] Mantel test complete: r = %.4f, p = %.4f\n", 
                result$statistic, result$signif))
    
    # Return results
    return(list(
      statistic = result$statistic,
      signif = result$signif,
      perm = result$perm,
      method = method,
      call = result$call
    ))
    
  }, error = function(e) {
    cat(sprintf(" Mantel test failed: %s\n", e$message))
    return(NULL)
  })
}


#' Perform Coinertia Analysis Between Two Data Tables
#'
#' @description
#' Conducts coinertia analysis (CoIA) to identify and visualize common structures 
#' between two data tables (e.g., species composition from different sampling methods). 
#' CoIA is a symmetric ordination method that finds axes maximizing covariance 
#' between two datasets.
#'
#' @details
#' Coinertia analysis (Doledec & Chessel, 1994) is a multivariate method for 
#' coupling two data tables. Unlike asymmetric methods (e.g., Canonical Correspondence 
#' Analysis), CoIA treats both tables symmetrically, making it ideal for comparing 
#' community patterns detected by different survey methods.
#' 
#' **Analysis Steps:**
#' 1. Performs PCA on each table separately
#' 2. Finds axes that maximize covariance between the two PCAs
#' 3. Projects samples onto common ordination space
#' 4. Calculates RV coefficient (multivariate correlation)
#' 
#' **Key Outputs:**
#' - **RV coefficient**: Overall similarity between tables (0-1 scale)
#' - **Eigenvalues**: Amount of co-structure captured by each axis
#' - **Sample scores**: Coordinates for plotting samples from both methods
#' - **Variable loadings**: Contribution of species/variables to axes
#' 
#' The function automatically standardizes data and handles cases with 
#' insufficient variation or too few observations.
#'
#' @param table1 Data frame or matrix. First data table (samples  variables)
#' @param table2 Data frame or matrix. Second data table (samples  variables, 
#'   same samples as table1)
#' @param nf Integer. Number of coinertia axes to retain. Default = 2
#' @param row.w Numeric vector. Row weights (sample weights). If NULL, uses 
#'   equal weights
#' @param col.w1 Numeric vector. Column weights for table1. If NULL, uses 
#'   equal weights
#' @param col.w2 Numeric vector. Column weights for table2. If NULL, uses 
#'   equal weights
#' @param scannf Logical. Should eigenvalues be plotted for axis selection? 
#'   Default = FALSE
#'
#' @return 
#' Returns a coinertia object (class 'coinertia', 'dudi') with the following 
#' main components, or NULL if analysis fails:
#' \describe{
#'   \item{RV}{Numeric. RV coefficient measuring overall similarity (0-1)}
#'   \item{eig}{Numeric vector. Eigenvalues for each axis}
#'   \item{lw}{Numeric vector. Sample weights}
#'   \item{mX}{Data frame. Sample scores from table1}
#'   \item{mY}{Data frame. Sample scores from table2}
#'   \item{aX}{Data frame. Variable loadings for table1}
#'   \item{aY}{Data frame. Variable loadings for table2}
#'   \item{call}{Language. The function call}
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Compare acoustic and metabarcoding community data
#' library(ade4)
#' library(vegan)
#' 
#' # Prepare data (same sites in both tables)
#' acoustic_data <- acoustic_matrix[common_sites, ]
#' metabar_data <- metabar_matrix[common_sites, ]
#' 
#' # Perform coinertia analysis
#' coia_result <- perform_coinertia_analysis(
#'   table1 = acoustic_data,
#'   table2 = metabar_data,
#'   nf = 2
#' )
#' 
#' if (!is.null(coia_result)) {
#'   cat(sprintf("RV coefficient: %.3f\n", coia_result$RV))
#'   cat(sprintf("Variance explained: %.1f%%\n", 
#'               100 * sum(coia_result$eig[1:2]) / sum(coia_result$eig)))
#'   
#'   # Plot results
#'   plot(coia_result)
#' }
#' }
#'
#' @references
#' Doledec, S., & Chessel, D. (1994). Co-inertia analysis: An alternative method 
#' for studying species-environment relationships. *Freshwater Biology*, 31(3), 
#' 277-294. https://doi.org/10.1111/j.1365-2427.1994.tb01741.x
#' 
#' Dray, S., Chessel, D., & Thioulouse, J. (2003). Co-inertia analysis and the 
#' linking of ecological data tables. *Ecology*, 84(11), 3078-3089. 
#' https://doi.org/10.1890/03-0178
#' 
#' Thioulouse, J., Dray, S., Dufour, A. B., Siberchicot, A., Jombart, T., & 
#' Pavoine, S. (2018). *Multivariate Analysis of Ecological Data with ade4*. 
#' Springer. https://doi.org/10.1007/978-1-4939-8850-1
#' 
#' Dray, S., Dufour, A. B., & Thioulouse, J. (2007). The ade4 package - II: 
#' Two-table and K-table methods. *R News*, 7(2), 47-52.
#'
#' @seealso 
#' \code{\link[ade4]{coinertia}} for the underlying implementation
#' \code{\link[ade4]{dudi.pca}} for PCA implementation
#' \code{\link{analyze_site_similarity}} which uses this function
#'
#' @export
#' 
#' @importFrom ade4 coinertia dudi.pca
perform_coinertia_analysis <- function(table1, table2, nf = 2, 
                                       row.w = NULL, col.w1 = NULL, col.w2 = NULL,
                                       scannf = FALSE) {
  
  # Load required package
  if (!requireNamespace("ade4", quietly = TRUE)) {
    stop("Package 'ade4' is required for coinertia analysis. Please install it.")
  }
  
  cat("\n[>>] Performing coinertia analysis...\n")
  
  # Validate inputs
  if (is.null(table1) || is.null(table2)) {
    cat(" Error: One or both tables are NULL\n")
    return(NULL)
  }
  
  # Convert to data frames if needed
  table1 <- as.data.frame(table1)
  table2 <- as.data.frame(table2)
  
  # Check dimensions
  if (nrow(table1) != nrow(table2)) {
    cat(sprintf("[X] Error: Tables have different numbers of rows (%d vs %d)\n", 
                nrow(table1), nrow(table2)))
    return(NULL)
  }
  
  if (nrow(table1) < 3) {
    cat(sprintf("[X] Error: Insufficient samples (n=%d, need at least 3)\n", nrow(table1)))
    return(NULL)
  }
  
  if (ncol(table1) < 2 || ncol(table2) < 2) {
    cat(sprintf("[X] Error: Insufficient variables (table1: %d, table2: %d, need at least 2 each)\n", 
                ncol(table1), ncol(table2)))
    return(NULL)
  }
  
  cat(sprintf("   Tables: %d samples, %d vs %d variables\n", 
              nrow(table1), ncol(table1), ncol(table2)))
  
  # Perform coinertia analysis with error handling
  tryCatch({
    
    # Step 1: PCA on each table
    cat("   Step 1/3: Performing PCA on table 1...\n")
    pca1 <- ade4::dudi.pca(table1, scannf = FALSE, nf = nf,
                           row.w = row.w, col.w = col.w1)
    
    cat("   Step 2/3: Performing PCA on table 2...\n")
    pca2 <- ade4::dudi.pca(table2, scannf = FALSE, nf = nf,
                           row.w = row.w, col.w = col.w2)
    
    # Step 2: Coinertia analysis
    cat("   Step 3/3: Computing coinertia...\n")
    coia <- ade4::coinertia(pca1, pca2, scannf = scannf, nf = nf)
    
    # Report results
    cat(sprintf("[OK] Coinertia analysis complete\n"))
    cat(sprintf("   RV coefficient: %.4f\n", coia$RV))
    cat(sprintf("   Eigenvalues (first %d axes):\n", min(nf, length(coia$eig))))
    for (i in 1:min(nf, length(coia$eig))) {
      cat(sprintf("     Axis %d: %.4f (%.1f%% variance)\n", 
                  i, coia$eig[i], 100 * coia$eig[i] / sum(coia$eig)))
    }
    
    return(coia)
    
  }, error = function(e) {
    cat(sprintf(" Coinertia analysis failed: %s\n", e$message))
    cat("   This may occur if:\n")
    cat("   - Tables have too much missing data\n")
    cat("   - Variables have zero variance\n")
    cat("   - Tables are too similar (no variation to explain)\n")
    return(NULL)
  })
}


# ==============================================================================
# SECTION 10: METHOD COMPARISON AND ADDITIONAL ANALYSES
# ==============================================================================
# Functions for comparing bioacoustic vs metabarcoding detection methods
# ==============================================================================

compare_katydid_methods <- function(acoustic_data, metabarcoding_data, 
                                    output_dir = "method_comparison_results",
                                    method_names = c("Bioacoustic", "Metabarcoding")) {
  
  cat("\n[>>] COMPREHENSIVE METHOD COMPARISON\n")
  cat("=====================================\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_files <- character()
  
  # -------------------------------------------------------------------------
  # 1. SPECIES OVERLAP ANALYSIS
  # -------------------------------------------------------------------------
  cat("\n[>>] Step 1/4: Analyzing species overlap...\n")
  
  # Extract species lists
  acoustic_species <- colnames(acoustic_data$presence_matrix)[
    !colnames(acoustic_data$presence_matrix) %in% c("site")]
  
  metabar_species <- colnames(metabarcoding_data$presence_matrix)[
    !colnames(metabarcoding_data$presence_matrix) %in% c("site")]
  
  # Calculate overlaps
  only_acoustic <- setdiff(acoustic_species, metabar_species)
  only_metabar <- setdiff(metabar_species, acoustic_species)
  both_methods <- intersect(acoustic_species, metabar_species)
  total_species <- length(union(acoustic_species, metabar_species))
  
  # Detection statistics
  detection_stats <- data.frame(
    Method = c(paste(method_names[1], "only"), 
               paste(method_names[2], "only"), 
               "Both methods"),
    Species_count = c(length(only_acoustic), length(only_metabar), length(both_methods)),
    Percentage = c(
      round(length(only_acoustic) / total_species * 100, 1),
      round(length(only_metabar) / total_species * 100, 1),
      round(length(both_methods) / total_species * 100, 1)
    )
  )
  
  cat(sprintf("   Total species pool: %d\n", total_species))
  cat(sprintf("   %s only: %d (%.1f%%)\n", method_names[1], 
              length(only_acoustic), detection_stats$Percentage[1]))
  cat(sprintf("   %s only: %d (%.1f%%)\n", method_names[2], 
              length(only_metabar), detection_stats$Percentage[2]))
  cat(sprintf("   Both methods: %d (%.1f%%)\n", 
              length(both_methods), detection_stats$Percentage[3]))
  
  # Create Venn diagram
  venn_plot <- NULL
  if (length(acoustic_species) > 0 && length(metabar_species) > 0) {
    tryCatch({
      venn_data <- list()
      venn_data[[method_names[1]]] <- acoustic_species
      venn_data[[method_names[2]]] <- metabar_species
      
      venn_plot <- VennDiagram::venn.diagram(
        venn_data,
        filename = NULL,
        fill = c("lightblue", "lightgreen"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        main = paste("Species Detection:", paste(method_names, collapse = " vs "))
      )
      
      # Save Venn diagram
      venn_file <- file.path(output_dir, "species_venn_diagram.png")
      png(venn_file, width = 800, height = 600)
      grid::grid.draw(venn_plot)
      dev.off()
      output_files <- c(output_files, venn_file)
      cat(sprintf("   [OK] Venn diagram saved: %s\n", basename(venn_file)))
      
    }, error = function(e) {
      cat(sprintf("     Venn diagram creation failed: %s\n", e$message))
    })
  }
  
  # Save detection statistics
  stats_file <- file.path(output_dir, "detection_statistics.csv")
  write.csv(detection_stats, stats_file, row.names = FALSE)
  output_files <- c(output_files, stats_file)
  
  # -------------------------------------------------------------------------
  # 2. FAMILY-LEVEL ANALYSIS
  # -------------------------------------------------------------------------
  cat("\n[>>] Step 2/4: Analyzing family-level patterns...\n")
  
  # Acoustic family summary (all Tettigoniidae for katydids)
  acoustic_summary <- data.frame(
    family = "Tettigoniidae",
    acoustic_species_count = length(acoustic_species),
    stringsAsFactors = FALSE
  )
  
  # Metabarcoding family summary
  if ("family" %in% colnames(metabarcoding_data$raw_data)) {
    metabar_summary <- metabarcoding_data$raw_data %>%
      dplyr::filter(!is.na(species) & species != "") %>%
      dplyr::select(family, species) %>%
      dplyr::distinct() %>%
      dplyr::group_by(family) %>%
      dplyr::summarise(metabar_species_count = dplyr::n(), .groups = "drop")
    
    family_comparison <- dplyr::full_join(acoustic_summary, metabar_summary, by = "family") %>%
      dplyr::mutate(
        acoustic_species_count = ifelse(is.na(acoustic_species_count), 0, acoustic_species_count),
        metabar_species_count = ifelse(is.na(metabar_species_count), 0, metabar_species_count)
      )
  } else {
    family_comparison <- acoustic_summary %>%
      dplyr::mutate(metabar_species_count = length(metabar_species))
  }
  
  cat(sprintf("   Families analyzed: %d\n", nrow(family_comparison)))
  
  # Save family analysis
  family_file <- file.path(output_dir, "family_comparison.csv")
  write.csv(family_comparison, family_file, row.names = FALSE)
  output_files <- c(output_files, family_file)
  
  # -------------------------------------------------------------------------
  # 3. CRYPTIC DIVERSITY ANALYSIS
  # -------------------------------------------------------------------------
  cat("\n[>>] Step 3/4: Analyzing cryptic diversity...\n")
  
  diversity_analysis <- NULL
  diversity_plot <- NULL
  
  if ("bin_uri" %in% colnames(metabarcoding_data$raw_data) && 
      "family" %in% colnames(metabarcoding_data$raw_data)) {
    
    diversity_analysis <- metabarcoding_data$raw_data %>%
      dplyr::group_by(family) %>%
      dplyr::summarise(
        total_bins = dplyr::n_distinct(bin_uri),
        unique_named_species = dplyr::n_distinct(species[!is.na(species) & species != ""]),
        identification_rate = round(unique_named_species / total_bins * 100, 1),
        cryptic_diversity = total_bins - unique_named_species,
        .groups = "drop"
      )
    
    cat(sprintf("   Total BINs (molecular species): %d\n", sum(diversity_analysis$total_bins)))
    cat(sprintf("   Named species: %d\n", sum(diversity_analysis$unique_named_species)))
    cat(sprintf("   Cryptic diversity: %d BINs\n", sum(diversity_analysis$cryptic_diversity)))
    cat(sprintf("   Overall identification rate: %.1f%%\n", 
                mean(diversity_analysis$identification_rate)))
    
    # Create diversity plot
    tryCatch({
      diversity_plot <- diversity_analysis %>%
        dplyr::select(family, unique_named_species, cryptic_diversity) %>%
        tidyr::pivot_longer(
          cols = c(unique_named_species, cryptic_diversity), 
          names_to = "diversity_type", 
          values_to = "count"
        ) %>%
        dplyr::mutate(
          diversity_type = dplyr::case_when(
            diversity_type == "unique_named_species" ~ "Named Species",
            diversity_type == "cryptic_diversity" ~ "Cryptic Diversity (Unnamed BINs)",
            TRUE ~ diversity_type
          )
        ) %>%
        ggplot2::ggplot(ggplot2::aes(x = reorder(family, count), y = count, fill = diversity_type)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
          values = c("Named Species" = "darkgreen", 
                     "Cryptic Diversity (Unnamed BINs)" = "orange")
        ) +
        ggplot2::labs(
          title = "Taxonomic vs Molecular Diversity",
          subtitle = "Named species vs cryptic molecular lineages (BINs)",
          x = "Family", 
          y = "Number of Taxa", 
          fill = "Diversity Type"
        ) +
        ggplot2::theme_minimal()
      
      # Save diversity plot
      plot_file <- file.path(output_dir, "cryptic_diversity_by_family.jpg")
      ggplot2::ggsave(plot_file, diversity_plot, width = 12, height = 8, dpi = 300)
      output_files <- c(output_files, plot_file)
      cat(sprintf("   [OK] Diversity plot saved: %s\n", basename(plot_file)))
      
    }, error = function(e) {
      cat(sprintf("     Diversity plot creation failed: %s\n", e$message))
    })
    
    # Save diversity analysis
    diversity_file <- file.path(output_dir, "cryptic_diversity_analysis.csv")
    write.csv(diversity_analysis, diversity_file, row.names = FALSE)
    output_files <- c(output_files, diversity_file)
  } else {
    cat("     Insufficient data for cryptic diversity analysis\n")
  }
  
  # -------------------------------------------------------------------------
  # 4. SUMMARY STATISTICS
  # -------------------------------------------------------------------------
  cat("\n[>>] Step 4/4: Compiling summary statistics...\n")
  
  summary_stats <- data.frame(
    Metric = c(
      "Total unique species",
      paste(method_names[1], "species"),
      paste(method_names[2], "species"),
      "Species in both methods",
      "Overlap percentage",
      "Method complementarity"
    ),
    Value = c(
      total_species,
      length(acoustic_species),
      length(metabar_species),
      length(both_methods),
      round(length(both_methods) / total_species * 100, 1),
      round((length(only_acoustic) + length(only_metabar)) / total_species * 100, 1)
    )
  )
  
  # Save summary
  summary_file <- file.path(output_dir, "comparison_summary.csv")
  write.csv(summary_stats, summary_file, row.names = FALSE)
  output_files <- c(output_files, summary_file)
  
  cat("\n[OK] COMPARISON COMPLETE\n")
  cat(sprintf("   Output files saved to: %s\n", output_dir))
  cat(sprintf("   Total files created: %d\n", length(output_files)))
  
  # Return comprehensive results
  return(list(
    katydid_comparison = list(
      stats = detection_stats,
      venn_plot = venn_plot,
      only_acoustic = only_acoustic,
      only_metabar = only_metabar,
      both_methods = both_methods,
      total_species = total_species,
      overlap_percentage = round(length(both_methods) / total_species * 100, 1)
    ),
    family_analysis = family_comparison,
    diversity_analysis = list(
      analysis = diversity_analysis,
      plot = diversity_plot
    ),
    summary_stats = summary_stats,
    output_files = output_files
  ))
}


#' Calculate Beta Diversity Components
#'
#' @description
#' Calculates beta diversity and partitions it into turnover and nestedness 
#' components using the betapart framework. Beta diversity measures the 
#' variation in species composition among sites.
#'
#' @details
#' Uses the additive partitioning framework of Baselga (2010) to decompose 
#' Srensen dissimilarity into:
#' - **Turnover** (sim): Species replacement between sites
#' - **Nestedness** (nes): Differences due to species loss/gain
#' - **Total beta** (sor): Overall compositional dissimilarity
#' 
#' Relationship: sor = sim + nes
#'
#' @param community_matrix Matrix or data frame. Sites (rows)  species (columns) 
#'   with presence/absence (0/1) or abundance data
#' @param index Character. Dissimilarity index: "sorensen" (default) or "jaccard"
#'
#' @return List with beta diversity matrices (dist objects): beta.sim, beta.sne, beta.sor
#'
#' @references
#' Baselga, A. (2010). Partitioning the turnover and nestedness components of 
#' beta diversity. *Global Ecology and Biogeography*, 19(1), 134-143.
#'
#' @export
#' @importFrom betapart beta.pair
calculate_beta_diversity <- function(community_matrix, index = "sorensen") {
  
  if (!requireNamespace("betapart", quietly = TRUE)) {
    stop("Package 'betapart' is required. Please install it.")
  }
  
  cat("\n[>>] Calculating beta diversity components...\n")
  
  # Convert to presence/absence if needed
  community_matrix <- ifelse(community_matrix > 0, 1, 0)
  
  tryCatch({
    beta <- betapart::beta.pair(community_matrix, index.family = index)
    
    cat(sprintf("[OK] Beta diversity calculated (index: %s)\n", index))
    cat(sprintf("   Mean turnover: %.3f\n", mean(beta$beta.sim)))
    cat(sprintf("   Mean nestedness: %.3f\n", mean(beta$beta.sne)))
    cat(sprintf("   Mean total beta: %.3f\n", mean(beta$beta.sor)))
    
    return(beta)
  }, error = function(e) {
    cat(sprintf(" Beta diversity calculation failed: %s\n", e$message))
    return(NULL)
  })
}


#' Perform PERMANOVA Test
#'
#' @description
#' Conducts Permutational Multivariate Analysis of Variance (PERMANOVA) to test 
#' for differences in community composition among groups.
#'
#' @details
#' PERMANOVA (Anderson, 2001) tests the null hypothesis that group centroids 
#' are equivalent in multivariate space using permutation tests. Unlike parametric 
#' MANOVA, it makes no distributional assumptions.
#'
#' @param distance_matrix Distance matrix (dist object or square matrix)
#' @param factors Data frame with grouping factors (same order as distance matrix)
#' @param formula Formula specifying model (e.g., ~ treatment + block)
#' @param permutations Integer. Number of permutations. Default = 999
#' @param method Character. Distance method if recalculating. Default = "bray"
#'
#' @return adonis2 object with PERMANOVA results
#'
#' @references
#' Anderson, M. J. (2001). A new method for non-parametric multivariate analysis 
#' of variance. *Austral Ecology*, 26(1), 32-46.
#'
#' @export
#' @importFrom vegan adonis2
perform_permanova <- function(distance_matrix, factors, formula, 
                              permutations = 999, method = "bray") {
  
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required. Please install it.")
  }
  
  cat("\n[>>] Performing PERMANOVA...\n")
  
  tryCatch({
    result <- vegan::adonis2(formula, data = factors, 
                             permutations = permutations, 
                             method = method)
    
    cat("[OK] PERMANOVA complete\n")
    print(result)
    
    return(result)
  }, error = function(e) {
    cat(sprintf(" PERMANOVA failed: %s\n", e$message))
    return(NULL)
  })
}


#' Fit GLM Models (Simple Wrapper)
#'
#' @description
#' Simple wrapper for fitting Generalized Linear Models with standard families.
#'
#' @param formula Formula. Model formula (e.g., response ~ predictor1 + predictor2)
#' @param data Data frame. Data for model fitting
#' @param family Character or family object. Error distribution: "gaussian", "poisson", 
#'   "binomial", etc. Default = "gaussian"
#'
#' @return glm object with fitted model, or NULL if fitting fails
#'
#' @export
#' @importFrom stats glm gaussian poisson binomial
fit_glm_models <- function(formula, data, family = "gaussian") {
  
  cat(sprintf("\n[>>] Fitting GLM (family: %s)...\n", family))
  
  # Convert family string to family object if needed
  if (is.character(family)) {
    family <- switch(family,
                     "gaussian" = gaussian(),
                     "poisson" = poisson(),
                     "binomial" = binomial(),
                     stop("Unknown family. Use 'gaussian', 'poisson', or 'binomial'"))
  }
  
  tryCatch({
    model <- glm(formula, data = data, family = family)
    
    cat("[OK] GLM fitted successfully\n")
    cat(sprintf("   AIC: %.2f\n", AIC(model)))
    cat(sprintf("   Deviance: %.2f\n", deviance(model)))
    
    return(model)
  }, error = function(e) {
    cat(sprintf(" GLM fitting failed: %s\n", e$message))
    return(NULL)
  })
}


#' Fit GAM Models (Simple Wrapper)
#'
#' @description
#' Simple wrapper for fitting Generalized Additive Models for non-linear relationships.
#'
#' @param formula Formula. Model formula with smooth terms (e.g., y ~ s(x1) + s(x2))
#' @param data Data frame. Data for model fitting
#' @param family Character or family object. Error distribution. Default = "gaussian"
#' @param method Character. Smoothness selection method. Default = "REML"
#'
#' @return gam object with fitted model, or NULL if fitting fails
#'
#' @export
#' @importFrom mgcv gam s
fit_gam_models <- function(formula, data, family = "gaussian", method = "REML") {
  
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is required. Please install it.")
  }
  
  cat(sprintf("\n[>>] Fitting GAM (family: %s, method: %s)...\n", family, method))
  
  tryCatch({
    model <- mgcv::gam(formula, data = data, family = family, method = method)
    
    cat("[OK] GAM fitted successfully\n")
    cat(sprintf("   Deviance explained: %.1f%%\n", 
                summary(model)$dev.expl * 100))
    cat(sprintf("   GCV score: %.4f\n", model$gcv.ubre))
    
    return(model)
  }, error = function(e) {
    cat(sprintf(" GAM fitting failed: %s\n", e$message))
    return(NULL)
  })
}

# ==============================================================================
# END OF SECTION 10
# ==============================================================================

# ==============================================================================
# END OF STATISTICAL_ANALYSIS.R - COMPLETE MERGED VERSION
# ==============================================================================
#
# Successfully merged all 10 sections:
# - Sections 1-7: Comprehensive statistical analyses (from BACKUP)
# - Section 8: GLM helper functions (COMPLETE)
# - Section 9: Multivariate wrappers (COMPLETE)
# - Section 10: Method comparison & additional analyses (COMPLETE)
#
# Date merged: 2025-10-19
# ==============================================================================
# Total lines: 3861


# ==============================================================================
# SECTION 11: MISSING FUNCTIONS FOR MAIN_ANALYSIS.R
# ==============================================================================
# Added: 2025-12-02
# Functions required by main_integrated_analysis() that were missing from 
# the refactored code. Extracted from bin_dont_tell_laurel_but_ill_call_it_final.R
# ==============================================================================

#' Compare Identified Katydids Only Between Methods
#'
#' Compares species detection between bioacoustic and metabarcoding methods,
#' focusing specifically on identified katydid species.
#'
#' @param acoustic_data List containing presence_matrix from acoustic detections
#' @param metabarcoding_data List containing metabarcoding species data
#' @param method_names Character vector of length 2 with method names
#'
#' @return List with stats, venn_plot, overlap_rate, both_methods, only_acoustic,
#'         only_metabar, and total_species
#'
#' @details Creates Venn diagram showing species overlap between methods.
#' Uses hardcoded list of identified katydids from metabarcoding for comparison.
#'
#' @references
#' Chen, H. (2018). VennDiagram: Generate High-Resolution Venn and Euler Plots.
#' R package version 1.6.20.
#'
#' @export
compare_identified_katydids_only <- function(acoustic_data, metabarcoding_data, 
                                             method_names = c("Bioacoustic", "Metabarcoding")) {
  
  cat("COMPARAISON 1: KATYDIDS IDENTIFIES SEULEMENT\n")
  
  # Extraire les especes bioacoustiques
  acoustic_species <- colnames(acoustic_data$presence_matrix)[
    !colnames(acoustic_data$presence_matrix) %in% c("site")]
  
  # Extraire dynamiquement les especes identifiees du metabarcoding
  # Chercher dans raw_data la colonne species (non-NA et non-vide)
  identified_katydids <- character(0)
  
  if (!is.null(metabarcoding_data$raw_data)) {
    # Chercher la colonne species dans raw_data
    species_col <- intersect(c("species", "Species", "species_name"), 
                             colnames(metabarcoding_data$raw_data))
    
    if (length(species_col) > 0) {
      # Extraire les especes uniques non-NA et non-vides
      species_values <- metabarcoding_data$raw_data[[species_col[1]]]
      identified_katydids <- unique(species_values[!is.na(species_values) & 
                                                    species_values != "" &
                                                    !grepl("^\\s*$", species_values)])
      cat(sprintf("  -> Extracted %d identified species from metabarcoding data\n", 
                  length(identified_katydids)))
    }
  }
  
  # Si species_matrix existe, utiliser ses colonnes comme fallback
  if (length(identified_katydids) == 0 && !is.null(metabarcoding_data$species_matrix)) {
    identified_katydids <- colnames(metabarcoding_data$species_matrix)[
      !colnames(metabarcoding_data$species_matrix) %in% c("site")]
    cat(sprintf("  -> Using species_matrix columns: %d species\n", 
                length(identified_katydids)))
  }
  
  # Fallback sur species_list si disponible
  if (length(identified_katydids) == 0 && !is.null(metabarcoding_data$species_list)) {
    identified_katydids <- metabarcoding_data$species_list[
      !is.na(metabarcoding_data$species_list) & 
        metabarcoding_data$species_list != ""]
    cat(sprintf("  -> Using species_list: %d species\n", 
                length(identified_katydids)))
  }
  
  cat(sprintf("  Acoustic species: %d\n", length(acoustic_species)))
  cat(sprintf("  Metabarcoding identified species: %d\n", length(identified_katydids)))
  
  # Analyse des chevauchements
  only_acoustic <- setdiff(acoustic_species, identified_katydids)
  only_metabar <- setdiff(identified_katydids, acoustic_species)
  both_methods <- intersect(acoustic_species, identified_katydids)
  
  # Statistiques
  total_katydid_species <- length(union(acoustic_species, identified_katydids))
  overlap_rate <- length(both_methods) / total_katydid_species * 100
  
  detection_stats <- data.frame(
    Method = c(paste(method_names[1], "only"), paste(method_names[2], "only"), "Both methods"),
    Species_count = c(length(only_acoustic), length(only_metabar), length(both_methods)),
    Percentage = c(
      length(only_acoustic) / total_katydid_species * 100,
      length(only_metabar) / total_katydid_species * 100,
      length(both_methods) / total_katydid_species * 100
    )
  )
  
  # Venn diagram pour katydids identifies
  venn_plot <- NULL
  if (length(acoustic_species) > 0 && length(identified_katydids) > 0) {
    tryCatch({
      venn_data <- list()
      venn_data[[method_names[1]]] <- acoustic_species
      venn_data[["Metabarcoding (Katydids only)"]] <- identified_katydids
      
      venn_plot <- VennDiagram::venn.diagram(
        venn_data,
        filename = NULL,
        fill = c("lightblue", "lightgreen"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        main = "Katydid Species Detection Overlap"
      )
    }, error = function(e) {
      cat("Warning: Could not create Venn diagram:", conditionMessage(e), "\n")
    })
  }
  
  return(list(
    stats = detection_stats,
    venn_plot = venn_plot,
    overlap_rate = overlap_rate,
    both_methods = both_methods,
    only_acoustic = only_acoustic,
    only_metabar = only_metabar,
    total_species = total_katydid_species
  ))
}


#' Analyze Detection by Family
#'
#' Compares species detection counts between methods grouped by taxonomic family.
#'
#' @param acoustic_data List containing presence_matrix
#' @param metabarcoding_data List containing raw_data with family information
#' @param family_mapping Optional data frame mapping species to families
#'
#' @return Data frame with family-level comparison statistics
#'
#' @export
analyze_detection_by_family <- function(acoustic_data, metabarcoding_data, family_mapping = NULL) {
  
  if (is.null(family_mapping)) {
    acoustic_families <- data.frame(
      species = colnames(acoustic_data$presence_matrix)[!colnames(acoustic_data$presence_matrix) %in% c("site")],
      family = "Tettigoniidae",
      stringsAsFactors = FALSE
    )
  } else {
    acoustic_families <- family_mapping
  }
  
  metabar_summary <- metabarcoding_data$raw_data %>%
    filter(!is.na(species) & species != "") %>%
    select(family, species) %>%
    distinct() %>%
    group_by(family) %>%
    summarise(metabar_species_count = n(), .groups = "drop")
  
  acoustic_summary <- acoustic_families %>%
    group_by(family) %>%
    summarise(acoustic_species_count = n(), .groups = "drop")
  
  family_comparison <- full_join(acoustic_summary, metabar_summary, by = "family") %>%
    mutate(
      acoustic_species_count = ifelse(is.na(acoustic_species_count), 0, acoustic_species_count),
      metabar_species_count = ifelse(is.na(metabar_species_count), 0, metabar_species_count)
    )
  
  return(family_comparison)
}


#' Analyze Taxonomic Diversity from Metabarcoding Data
#'
#' Analyzes cryptic diversity by comparing BIN counts to named species counts
#' within each taxonomic family.
#'
#' @param metabarcoding_data List containing raw_data with bin_uri, species, and family
#'
#' @return List with analysis data frame and diversity plot
#'
#' @details Calculates identification rate (named species / total BINs) and
#' cryptic diversity (BINs without species names) for each family.
#'
#' @export
analyze_taxonomic_diversity <- function(metabarcoding_data) {
  
  diversity_analysis <- metabarcoding_data$raw_data %>%
    group_by(family) %>%
    summarise(
      total_bins = n_distinct(bin_uri),
      unique_named_species = n_distinct(species[!is.na(species) & species != ""]),
      identification_rate = round(unique_named_species / total_bins * 100, 1),
      cryptic_diversity = total_bins - unique_named_species,
      .groups = "drop"
    )
  
  diversity_plot <- diversity_analysis %>%
    select(family, unique_named_species, cryptic_diversity) %>%
    pivot_longer(cols = c(unique_named_species, cryptic_diversity), 
                 names_to = "diversity_type", values_to = "count") %>%
    mutate(
      diversity_type = case_when(
        diversity_type == "unique_named_species" ~ "Named Species",
        diversity_type == "cryptic_diversity" ~ "Cryptic Diversity"
      )
    ) %>%
    ggplot(aes(x = reorder(family, count), y = count, fill = diversity_type)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("Named Species" = "darkgreen", "Cryptic Diversity" = "orange")) +
    labs(title = "Taxonomic Diversity: Named vs Cryptic Species",
         x = "Family", y = "Number of Taxa", fill = "Diversity Type") +
    theme_minimal()
  
  return(list(analysis = diversity_analysis, plot = diversity_plot))
}


#' Comprehensive Method Comparison
#'
#' Performs complete comparison between bioacoustic and metabarcoding detection
#' methods including species overlap, family-level analysis, and diversity metrics.
#'
#' @param acoustic_data List with presence_matrix from acoustic detections
#' @param metabarcoding_data List with presence_matrix and raw_data
#' @param output_dir Character string for output directory path
#'
#' @return List containing katydid_comparison, family_analysis, and diversity_analysis
#'
#' @details Orchestrates three complementary comparison analyses:
#' 1. Species-level overlap (Venn diagram)
#' 2. Family-level detection comparison
#' 3. Cryptic diversity analysis
#'
#' Results are saved as CSV files and PNG/JPG visualizations.
#'
#' @export
comprehensive_method_comparison <- function(acoustic_data, metabarcoding_data, 
                                            output_dir = "method_comparison_results") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Les 3 comparaisons
  katydid_comparison <- compare_identified_katydids_only(acoustic_data, metabarcoding_data)
  family_analysis <- analyze_detection_by_family(acoustic_data, metabarcoding_data)
  diversity_analysis <- analyze_taxonomic_diversity(metabarcoding_data)
  
  # Sauvegarder les visualisations
  if (!is.null(katydid_comparison$venn_plot)) {
    tryCatch({
      png(file.path(output_dir, "katydids_venn_diagram.png"), width = 800, height = 600)
      grid::grid.draw(katydid_comparison$venn_plot)
      dev.off()
    }, error = function(e) {
      cat("Warning: Could not save Venn diagram:", conditionMessage(e), "\n")
    })
  }
  
  tryCatch({
    ggsave(file.path(output_dir, "cryptic_diversity_by_family.jpg"), 
           diversity_analysis$plot, width = 12, height = 8)
  }, error = function(e) {
    cat("Warning: Could not save diversity plot:", conditionMessage(e), "\n")
  })
  
  # Sauvegarder les resultats
  write.csv(katydid_comparison$stats, file.path(output_dir, "katydid_comparison_stats.csv"), row.names = FALSE)
  write.csv(family_analysis, file.path(output_dir, "family_analysis.csv"), row.names = FALSE)
  write.csv(diversity_analysis$analysis, file.path(output_dir, "diversity_analysis.csv"), row.names = FALSE)
  
  return(list(
    katydid_comparison = katydid_comparison,
    family_analysis = family_analysis,
    diversity_analysis = diversity_analysis
  ))
}


#' Analyze Taxonomic Relationships Between Groups
#'
#' Analyzes correlation and relationships between species richness of two
#' taxonomic groups across sites, optionally incorporating vegetation data.
#'
#' @param group1_data List with presence_matrix for first group
#' @param group2_data List with presence_matrix for second group
#' @param vegetation_data Optional list with vegetation summary data
#' @param group_names Character vector of length 2 with group names
#'
#' @return List with correlation test, linear model, plot, and combined data
#'
#' @details Calculates site-level species richness for both groups, performs
#' Spearman correlation, fits linear model, and creates visualization.
#'
#' @references
#' Zar, J. H. (2010). Biostatistical Analysis. 5th edition. Pearson.
#'
#' @export
analyze_taxonomic_relationships <- function(group1_data, group2_data, vegetation_data = NULL,
                                            group_names = c("Group1", "Group2")) {
  
  # Prepare data for both groups
  group1_site_data <- group1_data$presence_matrix %>%
    mutate(
      site_number = as.numeric(gsub("S", "", site)),
      Plot = sprintf("YB-P%02d", site_number)
    )
  
  group1_richness <- group1_site_data %>%
    mutate(group1_richness = rowSums(select(., -c(site, site_number, Plot)))) %>%
    select(Plot, group1_richness)
  
  group2_site_data <- group2_data$presence_matrix %>%
    mutate(
      site_number = as.numeric(gsub("S", "", site)),
      Plot = sprintf("YB-P%02d", site_number)
    )
  
  group2_richness <- group2_site_data %>%
    mutate(group2_richness = rowSums(select(., -c(site, site_number, Plot)))) %>%
    select(Plot, group2_richness)
  
  # Merge data
  combined_data <- group1_richness %>%
    full_join(group2_richness, by = "Plot") %>%
    replace_na(list(group1_richness = 0, group2_richness = 0))
  
  # Add vegetation data if available
  if (!is.null(vegetation_data) && !is.null(vegetation_data$summary)) {
    combined_data <- combined_data %>%
      left_join(vegetation_data$summary, by = "Plot")
  }
  
  # Correlation between groups
  richness_correlation <- NULL
  if (sum(combined_data$group1_richness > 0) > 3 && 
      sum(combined_data$group2_richness > 0) > 3) {
    
    richness_correlation <- cor.test(
      combined_data$group1_richness,
      combined_data$group2_richness,
      method = "spearman"
    )
  }
  
  # Linear model
  group_model <- NULL
  if (nrow(combined_data) > 5) {
    group_model <- lm(group2_richness ~ group1_richness, data = combined_data)
  }
  
  # Visualization
  relationship_plot <- ggplot(combined_data, aes(x = group1_richness, y = group2_richness)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(title = paste(group_names[2], "vs", group_names[1], "Species Richness"),
         x = paste(group_names[1], "Species Richness"), 
         y = paste(group_names[2], "Species Richness")) +
    theme_minimal()
  
  # Add vegetation info if available
  if (!is.null(vegetation_data) && "liana_species_richness" %in% colnames(combined_data)) {
    relationship_plot <- relationship_plot +
      aes(size = liana_species_richness, color = tree_species_richness) +
      labs(size = "Liana richness", color = "Tree richness")
  }
  
  return(list(
    correlation = richness_correlation,
    model = group_model,
    plot = relationship_plot,
    data = combined_data
  ))
}


#' Integrated Multivariate Analysis (NMDS/PCoA)
#'
#' Performs NMDS ordination on species presence matrix with environmental
#' vector fitting. Falls back to PCoA if NMDS fails to converge.
#'
#' @param species_matrix Data frame with site column and species presence columns
#' @param vegetation_data List with summary data frame containing vegetation metrics
#' @param analysis_name Character string for labeling output
#'
#' @return List with nmds/pcoa object, plot, env_fit results, and method used
#'
#' @details
#' Performs the following steps:
#' 1. Filters species occurring at >= 2 sites and sites with >= 2 species
#' 2. Attempts NMDS with Jaccard distance (100 tries, stress optimization)
#' 3. If NMDS fails, performs PCoA as fallback
#' 4. Fits environmental vectors using envfit()
#' 5. Creates visualization colored by vegetation metrics
#'
#' @references
#' Oksanen, J., et al. (2020). vegan: Community Ecology Package.
#' 
#' Kruskal, J. B. (1964). Nonmetric multidimensional scaling: a numerical method.
#' Psychometrika, 29(2), 115-129.
#'
#' @export
analyze_multivariate_integrated <- function(species_matrix, vegetation_data, analysis_name = "Multivariate") {
  
  # Extract only species presence columns
  species_presence <- species_matrix %>%
    select(-site)
  
  # Filter species and sites with sufficient occurrences
  species_counts <- colSums(species_presence)
  species_presence <- species_presence[, species_counts >= 2, drop = FALSE]
  
  site_counts <- rowSums(species_presence)
  valid_sites <- site_counts >= 2
  
  cat(sprintf("Sites before filtering: %d\n", nrow(species_presence)))
  cat(sprintf("Sites after filtering: %d\n", sum(valid_sites)))
  cat(sprintf("Species before filtering: %d\n", ncol(species_matrix) - 1))
  cat(sprintf("Species after filtering: %d\n", ncol(species_presence)))
  
  if(sum(valid_sites) < 3) {
    warning(paste("Not enough sites with sufficient species for", analysis_name, "NMDS"))
    return(NULL)
  }
  
  # Filter matrices
  species_presence_filtered <- species_presence[valid_sites, , drop = FALSE]
  valid_site_ids <- species_matrix$site[valid_sites]
  
  # Only proceed if we have sufficient data
  if(ncol(species_presence_filtered) > 2 && nrow(species_presence_filtered) > 2) {
    
    # Try NMDS with better parameters
    nmds_result <- tryCatch({
      cat(paste("Attempting NMDS for", analysis_name, "with", nrow(species_presence_filtered), 
                "sites and", ncol(species_presence_filtered), "species\n"))
      
      # NMDS avec parametres ameliores
      species_nmds <- metaMDS(species_presence_filtered, 
                              distance = "jaccard", 
                              k = 2,
                              trymax = 100,
                              maxit = 1000,
                              autotransform = FALSE,
                              wascores = TRUE,
                              trace = 0)
      
      # Verifier si NMDS a converge
      if(species_nmds$converged && species_nmds$stress > 0.001) {
        cat("NMDS converged successfully with stress =", species_nmds$stress, "\n")
        
        # Extract NMDS scores
        nmds_scores <- as.data.frame(scores(species_nmds, display = "sites"))
        
        if(!"NMDS1" %in% colnames(nmds_scores)) {
          colnames(nmds_scores)[1:2] <- c("NMDS1", "NMDS2")
        }
        
        nmds_scores$site <- valid_site_ids
        
        # Convert site to Plot format
        nmds_scores <- nmds_scores %>%
          mutate(
            site_number = as.numeric(gsub("S", "", site)),
            Plot = sprintf("YB-P%02d", site_number)
          )
        
        # Merge with vegetation data
        nmds_data <- nmds_scores %>%
          left_join(vegetation_data$summary, by = "Plot") %>%
          filter(!is.na(tree_species_richness))
        
        if(nrow(nmds_data) < 3) {
          cat("Not enough sites with vegetation data after merge\n")
          return(NULL)
        }
        
        # Create NMDS plot
        nmds_plot <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
          geom_point(aes(size = tree_species_richness, color = liana_species_richness), alpha = 0.7) +
          geom_text(aes(label = site), vjust = -1, size = 3) +
          scale_color_viridis_c() +
          labs(title = paste("NMDS of", analysis_name, "Communities"),
               subtitle = paste("Stress =", round(species_nmds$stress, 3)),
               size = "Tree richness", color = "Liana richness") +
          theme_minimal()
        
        # Environmental fit
        env_data_filtered <- nmds_data %>% 
          select(any_of(c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                          "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area"))) %>%
          select(where(~ sum(!is.na(.)) >= 3))
        
        env_fit <- NULL
        if(ncol(env_data_filtered) > 0) {
          env_fit <- envfit(species_nmds, env_data_filtered, perm = 999)
        }
        
        return(list(nmds = species_nmds, plot = nmds_plot, env_fit = env_fit, method = "NMDS"))
        
      } else {
        cat("NMDS convergence issue, trying PCoA...\n")
        NULL
      }
      
    }, error = function(e) {
      cat(paste("NMDS failed for", analysis_name, ":", e$message, "- trying PCoA...\n"))
      NULL
    })
    
    if (!is.null(nmds_result)) {
      return(nmds_result)
    }
    
    # Fallback to PCoA if NMDS fails
    pcoa_result <- tryCatch({
      cat("Attempting PCoA as fallback...\n")
      
      jac_dist <- vegdist(species_presence_filtered, method = "jaccard")
      
      if(length(unique(as.vector(jac_dist))) < 3) {
        cat("Distance matrix is degenerate\n")
        return(NULL)
      }
      
      pcoa <- cmdscale(jac_dist, k = 2, eig = TRUE)
      
      pcoa_scores <- as.data.frame(pcoa$points)
      names(pcoa_scores) <- c("NMDS1", "NMDS2")
      pcoa_scores$site <- valid_site_ids
      
      pcoa_scores <- pcoa_scores %>%
        mutate(
          site_number = as.numeric(gsub("S", "", site)),
          Plot = sprintf("YB-P%02d", site_number)
        )
      
      pcoa_data <- pcoa_scores %>%
        left_join(vegetation_data$summary, by = "Plot") %>%
        filter(!is.na(tree_species_richness))
      
      if(nrow(pcoa_data) < 3) {
        cat("Not enough sites with vegetation data for PCoA\n")
        return(NULL)
      }
      
      # Variance explained
      eig_values <- pcoa$eig[pcoa$eig > 0]
      var_explained <- eig_values / sum(eig_values) * 100
      
      pcoa_plot <- ggplot(pcoa_data, aes(x = NMDS1, y = NMDS2)) +
        geom_point(aes(size = tree_species_richness, color = liana_species_richness), alpha = 0.7) +
        geom_text(aes(label = site), vjust = -1, size = 3) +
        scale_color_viridis_c() +
        labs(title = paste("PCoA of", analysis_name, "Communities"),
             subtitle = paste("PC1:", round(var_explained[1], 1), "% - PC2:", round(var_explained[2], 1), "%"),
             x = paste("PCoA1 (", round(var_explained[1], 1), "%)", sep = ""),
             y = paste("PCoA2 (", round(var_explained[2], 1), "%)", sep = ""),
             size = "Tree richness", color = "Liana richness") +
        theme_minimal()
      
      # Environmental fit for PCoA
      env_data_filtered <- pcoa_data %>% 
        select(any_of(c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                        "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area"))) %>%
        select(where(~ sum(!is.na(.)) >= 3))
      
      env_fit2 <- NULL
      if(ncol(env_data_filtered) > 0) {
        env_fit2 <- envfit(pcoa$points, env_data_filtered, perm = 999)
      }
      
      return(list(nmds = NULL, pcoa = pcoa, plot = pcoa_plot, env_fit = env_fit2, method = "PCoA"))
      
    }, error = function(e) {
      warning(paste("Both NMDS and PCoA failed for", analysis_name, ":", e$message))
      NULL
    })
    
    return(pcoa_result)
    
  } else {
    warning(paste("Insufficient data for multivariate analysis of", analysis_name))
    return(NULL)
  }
}


# ==============================================================================
# END OF SECTION 11 - MISSING FUNCTIONS ADDED
# ==============================================================================

# ==============================================================================
# ADDITIONAL FUNCTION: SPECIES RICHNESS RANGES
# ==============================================================================
# Added from reference script - required by main_analysis.R
# ==============================================================================

#' Calculate Species Richness Ranges Across Detection Methods
#'
#' Computes summary statistics for species richness per site across different
#' detection methods (bioacoustic, metabarcoding, birds).
#'
#' @param all_results List containing data from all detection methods
#'
#' @return List with richness statistics per method including min, max, mean,
#'   and formatted range text
#'
#' @export
calculate_species_richness_ranges <- function(all_results) {
  
  cat("\nCALCULATING SPECIES RICHNESS RANGES\n")
  
  richness_summary <- list()
  
  # Katydid richness
  if (!is.null(all_results$katydid_data) && !is.null(all_results$katydid_data$presence_matrix)) {
    katydid_matrix <- all_results$katydid_data$presence_matrix
    katydid_species_cols <- colnames(katydid_matrix)[colnames(katydid_matrix) != "site"]
    katydid_richness <- rowSums(katydid_matrix[, katydid_species_cols, drop = FALSE])
    
    richness_summary$katydid <- list(
      min = min(katydid_richness),
      max = max(katydid_richness),
      mean = round(mean(katydid_richness), 1),
      range_text = paste0(min(katydid_richness), "-", max(katydid_richness))
    )
  }
  
  # Bird richness
  if (!is.null(all_results$bird_data) && !is.null(all_results$bird_data$presence_matrix)) {
    bird_matrix <- all_results$bird_data$presence_matrix
    bird_species_cols <- colnames(bird_matrix)[colnames(bird_matrix) != "site"]
    bird_richness <- rowSums(bird_matrix[, bird_species_cols, drop = FALSE])
    
    richness_summary$bird <- list(
      min = min(bird_richness),
      max = max(bird_richness),
      mean = round(mean(bird_richness), 1),
      range_text = paste0(min(bird_richness), "-", max(bird_richness))
    )
  }
  
  # Metabarcoding richness
  if (!is.null(all_results$metabarcoding_data) && !is.null(all_results$metabarcoding_data$presence_matrix)) {
    meta_matrix <- all_results$metabarcoding_data$presence_matrix
    meta_species_cols <- colnames(meta_matrix)[colnames(meta_matrix) != "site"]
    meta_richness <- rowSums(meta_matrix[, meta_species_cols, drop = FALSE])
    
    richness_summary$metabarcoding <- list(
      min = min(meta_richness),
      max = max(meta_richness),
      mean = round(mean(meta_richness), 1),
      range_text = paste0(min(meta_richness), "-", max(meta_richness))
    )
  }
  
  # Print summary
  cat("SPECIES RICHNESS RANGES:\n")
  if (!is.null(richness_summary$katydid)) {
    cat(sprintf("  Katydids: %s species per site (mean: %.1f)\n", 
                richness_summary$katydid$range_text, richness_summary$katydid$mean))
  }
  if (!is.null(richness_summary$bird)) {
    cat(sprintf("  Birds: %s species per site (mean: %.1f)\n", 
                richness_summary$bird$range_text, richness_summary$bird$mean))
  }
  if (!is.null(richness_summary$metabarcoding)) {
    cat(sprintf("  Metabarcoding: %s species per site (mean: %.1f)\n", 
                richness_summary$metabarcoding$range_text, richness_summary$metabarcoding$mean))
  }
  
  return(richness_summary)
}


# ==============================================================================
# SECTION 12: ADDITIONAL MISSING FUNCTIONS
# ==============================================================================
# Added: 2025-12-18
# Functions required by main_integrated_analysis() that were missing
# Source: bin_dont_tell_laurel_but_ill_call_it_final.R
# ==============================================================================


#' Compare Detection Methods (Basic Comparison)
#'
#' Performs basic comparison between bioacoustic and metabarcoding detection 
#' methods, calculating species overlap and creating Venn diagram.
#'
#' @param acoustic_data List containing presence_matrix from acoustic detections
#' @param metabarcoding_data List containing presence_matrix from metabarcoding
#' @param method_names Character vector of length 2 with method names. 
#'   Default = c("Bioacoustic", "Metabarcoding")
#'
#' @return List with:
#' \describe{
#'   \item{stats}{Data frame with detection statistics by method}
#'   \item{venn_data}{List with species sets for Venn diagram}
#'   \item{venn_plot}{Venn diagram grid object or NULL}
#'   \item{only_acoustic}{Character vector of species only detected acoustically}
#'   \item{only_metabar}{Character vector of species only detected by metabarcoding}
#'   \item{both_methods}{Character vector of species detected by both methods}
#' }
#'
#' @details
#' This is a simpler version compared to compare_katydid_methods. Use this for
#' basic overlap analysis without the comprehensive family and diversity analyses.
#'
#' @export
#' @importFrom VennDiagram venn.diagram
compare_detection_methods <- function(acoustic_data, metabarcoding_data, 
                                      method_names = c("Bioacoustic", "Metabarcoding")) {
  
  cat("\n[>>] Comparing detection methods...\n")
  
  # Extract species lists from each method
  acoustic_species <- colnames(acoustic_data$presence_matrix)[
    !colnames(acoustic_data$presence_matrix) %in% c("site")]
  
  metabar_species <- colnames(metabarcoding_data$presence_matrix)[
    !colnames(metabarcoding_data$presence_matrix) %in% c("site")]
  
  # Analyze overlaps and exclusivities
  only_acoustic <- setdiff(acoustic_species, metabar_species)
  only_metabar <- setdiff(metabar_species, acoustic_species)
  both_methods <- intersect(acoustic_species, metabar_species)
  total_species <- length(union(acoustic_species, metabar_species))
  
  # Calculate detection statistics
  detection_stats <- data.frame(
    Method = c(paste(method_names[1], "only"), 
               paste(method_names[2], "only"), 
               "Both methods"),
    Species_count = c(length(only_acoustic), 
                      length(only_metabar), 
                      length(both_methods)),
    Percentage = c(
      length(only_acoustic) / total_species * 100,
      length(only_metabar) / total_species * 100,
      length(both_methods) / total_species * 100
    )
  )
  
  cat(sprintf("   %s species: %d\n", method_names[1], length(acoustic_species)))
  cat(sprintf("   %s species: %d\n", method_names[2], length(metabar_species)))
  cat(sprintf("   Species in both: %d (%.1f%%)\n", 
              length(both_methods), 
              length(both_methods) / total_species * 100))
  
  # Create Venn diagram
  venn_data <- list()
  venn_data[[method_names[1]]] <- acoustic_species
  venn_data[[method_names[2]]] <- metabar_species
  
  venn_plot <- NULL
  if (length(acoustic_species) > 0 && length(metabar_species) > 0) {
    tryCatch({
      venn_plot <- VennDiagram::venn.diagram(
        venn_data,
        filename = NULL,
        fill = c("lightblue", "lightgreen"),
        alpha = 0.5,
        cex = 1.5,
        cat.cex = 1.5,
        main = paste("Species Detection:", method_names[1], "vs", method_names[2])
      )
      cat("[OK] Venn diagram created\n")
    }, error = function(e) {
      cat(sprintf("   Warning: Venn diagram creation failed: %s\n", e$message))
    })
  } else {
    warning("Not enough data for Venn diagram - one method has no species detected")
  }
  
  return(list(
    stats = detection_stats,
    venn_data = venn_data,
    venn_plot = venn_plot,
    only_acoustic = only_acoustic,
    only_metabar = only_metabar,
    both_methods = both_methods
  ))
}


#' Run Additional Detection-Based Analyses
#'
#' Orchestrates additional multivariate analyses including NMDS with detection 
#' counts, CCA, and GAM for non-linear relationships.
#'
#' @param all_results List containing katydid_data, bird_data, and vegetation_data
#' @param output_dir Character string. Base output directory path
#'
#' @return List with results from nmds_detections, cca_detections, gam_katydid, 
#'   and gam_bird analyses
#'
#' @details
#' Calls the following analysis functions:
#' 1. analyze_nmds_with_detection_counts - NMDS ordination
#' 2. analyze_cca_with_detections - Canonical Correspondence Analysis
#' 3. analyze_nonlinear_relationships - GAM for katydids and birds
#'
#' Results are saved in a subdirectory "additional_detection_analyses".
#'
#' @export
run_additional_analyses <- function(all_results, output_dir) {
  
  cat("\n[>>] RUNNING ADDITIONAL DETECTION-BASED ANALYSES\n")
  cat("===============================================\n")
  
  # Create output subdirectory
  additional_output_dir <- file.path(output_dir, "additional_detection_analyses")
  dir.create(additional_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. NMDS with detection counts
  cat("\n1. NMDS with Detection Counts...\n")
  nmds_results <- tryCatch({
    analyze_nmds_with_detection_counts(
      all_results$katydid_data,
      all_results$bird_data,
      all_results$vegetation_data,
      additional_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   NMDS analysis failed: %s\n", e$message))
    NULL
  })
  
  # 2. CCA
  cat("\n2. Canonical Correspondence Analysis...\n")
  cca_results <- tryCatch({
    analyze_cca_with_detections(
      all_results$katydid_data,
      all_results$bird_data,
      all_results$vegetation_data,
      additional_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   CCA analysis failed: %s\n", e$message))
    NULL
  })
  
  # 3. Non-linear analyses (GAM)
  cat("\n3. Non-linear Relationship Analysis (GAM)...\n")
  gam_katydid_results <- tryCatch({
    analyze_nonlinear_relationships(
      all_results$katydid_data,
      all_results$vegetation_data,
      "katydid",
      additional_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   Katydid GAM analysis failed: %s\n", e$message))
    NULL
  })
  
  gam_bird_results <- tryCatch({
    analyze_nonlinear_relationships(
      all_results$bird_data,
      all_results$vegetation_data,
      "bird",
      additional_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   Bird GAM analysis failed: %s\n", e$message))
    NULL
  })
  
  # Compile results
  additional_results <- list(
    nmds_detections = nmds_results,
    cca_detections = cca_results,
    gam_katydid = gam_katydid_results,
    gam_bird = gam_bird_results
  )
  
  cat(sprintf("\n[OK] Additional analyses complete! Results saved in: %s\n", 
              additional_output_dir))
  
  return(additional_results)
}


#' Robust NMDS with Bootstrap Confidence Intervals
#'
#' Performs NMDS ordination with bootstrap resampling to calculate confidence 
#' ellipses for site positions, providing a measure of ordination stability.
#'
#' @param katydid_data List containing raw_detections data frame
#' @param bird_data List containing bird detection data (currently unused but 
#'   included for API consistency)
#' @param vegetation_enhanced Data frame with enhanced environmental variables 
#'   (from create_enhanced_environmental_variables)
#' @param output_dir Character string. Output directory for saving results
#'
#' @return List with:
#' \describe{
#'   \item{nmds}{metaMDS object from main NMDS}
#'   \item{bootstrap_results}{3D array of bootstrap site scores}
#'   \item{confidence_ellipses}{Data frame with ellipse parameters per site}
#'   \item{env_fit}{envfit object with environmental vector fitting}
#'   \item{plot}{ggplot object with NMDS and confidence ellipses}
#'   \item{selected_vars}{Character vector of selected environmental variables}
#' }
#'
#' @details
#' This function:
#' 1. Selects optimal environmental variables via forward selection
#' 2. Performs main NMDS on log-transformed detection counts
#' 3. Runs 99 bootstrap iterations with Procrustes alignment
#' 4. Calculates 95% confidence ellipses from bootstrap results
#' 5. Fits environmental vectors to ordination
#' 6. Creates visualization with ellipses and significant vectors
#'
#' @references
#' Oksanen, J. et al. (2020). vegan: Community Ecology Package.
#' 
#' Peres-Neto, P. R., & Jackson, D. A. (2001). How well do multivariate data 
#' sets match? The advantages of a Procrustean superimposition approach over 
#' the Mantel test. Oecologia, 129(2), 169-178.
#'
#' @export
#' @importFrom vegan metaMDS scores envfit procrustes
#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_segment geom_path 
#'   scale_size_continuous labs theme_minimal theme element_text arrow unit
analyze_robust_nmds <- function(katydid_data, bird_data, vegetation_enhanced, output_dir) {
  
  cat("\n[>>] ROBUST NMDS WITH BOOTSTRAP\n")
  cat("================================\n")
  
  # Select optimal variables
  optimal_vars <- select_optimal_variables(
    create_detection_count_matrix(katydid_data$raw_detections, "common_name", "site"),
    vegetation_enhanced
  )
  
  if (is.null(optimal_vars)) {
    cat("[X] Insufficient data for robust NMDS\n")
    return(NULL)
  }
  
  # NMDS with bootstrap for confidence intervals
  species_log <- log1p(optimal_vars$species_matrix)
  
  tryCatch({
    # Main NMDS
    cat("   Fitting main NMDS...\n")
    nmds_main <- metaMDS(species_log, distance = "bray", k = 2, 
                         trymax = 200, maxit = 2000, trace = 0)
    
    if (!nmds_main$converged) {
      cat("[X] Main NMDS did not converge\n")
      return(NULL)
    }
    
    cat(sprintf("   Main NMDS stress: %.4f\n", nmds_main$stress))
    
    # Bootstrap for confidence intervals
    cat("   Computing bootstrap confidence intervals (99 iterations)...\n")
    n_boot <- 99
    boot_results <- array(NA, dim = c(nrow(species_log), 2, n_boot))
    
    for (i in 1:n_boot) {
      if (i %% 20 == 0) cat(sprintf("     Bootstrap %d/%d...\n", i, n_boot))
      
      # Bootstrap sampling
      boot_indices <- sample(nrow(species_log), replace = TRUE)
      boot_data <- species_log[boot_indices, ]
      
      # Bootstrap NMDS
      nmds_i <- tryCatch({
        metaMDS(boot_data, distance = "bray", k = 2, 
                trymax = 50, trace = 0, previous.best = nmds_main)
      }, error = function(e) NULL)
      
      if (!is.null(nmds_i) && nmds_i$converged) {
        # Procrustes alignment to main NMDS
        proc_result <- procrustes(nmds_main, nmds_i)
        boot_results[boot_indices, , i] <- proc_result$Yrot
      }
    }
    
    # Calculate confidence ellipses
    site_scores <- scores(nmds_main, display = "sites")
    conf_ellipses <- data.frame()
    
    for (i in 1:nrow(site_scores)) {
      valid_boots <- boot_results[i, , !is.na(boot_results[i, 1, ])]
      if (is.matrix(valid_boots) && ncol(valid_boots) > 10) {
        ellipse_data <- data.frame(
          site = optimal_vars$site_data$site[i],
          NMDS1_mean = mean(valid_boots[1, ]),
          NMDS2_mean = mean(valid_boots[2, ]),
          NMDS1_sd = sd(valid_boots[1, ]),
          NMDS2_sd = sd(valid_boots[2, ])
        )
        conf_ellipses <- rbind(conf_ellipses, ellipse_data)
      }
    }
    
    # Environmental fitting with selected variables
    cat("   Fitting environmental vectors...\n")
    env_fit_robust <- envfit(nmds_main, optimal_vars$env_matrix, perm = 999)
    
    # Create enhanced plot
    nmds_scores_df <- as.data.frame(site_scores)
    nmds_scores_df$site <- optimal_vars$site_data$site
    nmds_scores_df <- nmds_scores_df %>%
      left_join(optimal_vars$site_data, by = "site")
    
    # Plot with confidence ellipses
    robust_plot <- ggplot(nmds_scores_df, aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(size = rowSums(optimal_vars$species_matrix)), 
                 alpha = 0.8, color = "darkblue") +
      geom_text(aes(label = site), vjust = -1.5, size = 3, fontface = "bold") +
      scale_size_continuous(name = "Total\nDetections", range = c(3, 8)) +
      labs(
        title = "Robust NMDS with Bootstrap Confidence",
        subtitle = paste("Stress =", round(nmds_main$stress, 3), 
                         "| Selected variables:", 
                         paste(optimal_vars$selected_variables, collapse = ", ")),
        x = "NMDS1", y = "NMDS2"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"))
    
    # Add ellipses if available
    if (nrow(conf_ellipses) > 0) {
      ellipse_data <- data.frame()
      
      for (i in 1:nrow(conf_ellipses)) {
        # Create 95% confidence ellipse
        theta <- seq(0, 2*pi, length.out = 100)
        ellipse_x <- conf_ellipses$NMDS1_mean[i] + 1.96 * conf_ellipses$NMDS1_sd[i] * cos(theta)
        ellipse_y <- conf_ellipses$NMDS2_mean[i] + 1.96 * conf_ellipses$NMDS2_sd[i] * sin(theta)
        
        site_ellipse <- data.frame(
          NMDS1 = ellipse_x,
          NMDS2 = ellipse_y,
          site = conf_ellipses$site[i]
        )
        ellipse_data <- rbind(ellipse_data, site_ellipse)
      }
      
      robust_plot <- robust_plot +
        geom_path(data = ellipse_data, 
                  aes(x = NMDS1, y = NMDS2, group = site),
                  alpha = 0.4, color = "red", linewidth = 0.8, inherit.aes = FALSE)
    }
    
    # Add significant environmental vectors
    if (any(env_fit_robust$vectors$pvals < 0.1)) {
      sig_vectors <- env_fit_robust$vectors$arrows[env_fit_robust$vectors$pvals < 0.1, , drop = FALSE]
      
      if (nrow(sig_vectors) > 0) {
        vector_df <- data.frame(sig_vectors * 0.7)
        vector_df$variable <- rownames(sig_vectors)
        
        robust_plot <- robust_plot +
          geom_segment(data = vector_df,
                       aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                       arrow = arrow(length = unit(0.3, "cm")), 
                       color = "red", linewidth = 1.2, inherit.aes = FALSE) +
          geom_text(data = data.frame(sig_vectors * 0.8, 
                                      variable = rownames(sig_vectors)),
                    aes(x = NMDS1, y = NMDS2, label = variable),
                    color = "red", size = 3, fontface = "bold", inherit.aes = FALSE)
      }
    }
    
    # Save plot
    ggsave(file.path(output_dir, "robust_nmds_bootstrap.jpg"), 
           robust_plot, width = 14, height = 10, dpi = 300)
    
    cat("[OK] Robust NMDS complete\n")
    cat(sprintf("   Plot saved: robust_nmds_bootstrap.jpg\n"))
    
    return(list(
      nmds = nmds_main,
      bootstrap_results = boot_results,
      confidence_ellipses = conf_ellipses,
      env_fit = env_fit_robust,
      plot = robust_plot,
      selected_vars = optimal_vars$selected_variables
    ))
    
  }, error = function(e) {
    cat(sprintf("[X] Robust NMDS failed: %s\n", e$message))
    return(NULL)
  })
}


#' Analyze Variance Partitioning
#'
#' Partitions variance in katydid community composition between structural 
#' and compositional vegetation variables using the varpart function.
#'
#' @param katydid_data List containing raw_detections data frame
#' @param vegetation_enhanced Data frame with enhanced environmental variables
#' @param output_dir Character string. Output directory for saving results
#'
#' @return List with:
#' \describe{
#'   \item{varpart}{varpart object with partitioning results}
#'   \item{plot}{ggplot bar chart of variance components}
#'   \item{partition_data}{Data frame with component names and variance values}
#' }
#'
#' @details
#' Partitions variance into:
#' - Structure: structural_complexity, log_tree_abundance, log_liana_stems
#' - Composition: total_plant_richness, liana_dominance, liana_to_tree_ratio
#' - Shared: Variance explained by both sets
#' - Residual: Unexplained variance
#'
#' @references
#' Borcard, D., Legendre, P., & Drapeau, P. (1992). Partialling out the spatial 
#' component of ecological variation. Ecology, 73(3), 1045-1055.
#'
#' @export
#' @importFrom vegan varpart
analyze_variance_partitioning <- function(katydid_data, vegetation_enhanced, output_dir) {
  
  cat("\n[>>] VARIANCE PARTITIONING\n")
  cat("==========================\n")
  
  # Prepare data
  detection_counts <- create_detection_count_matrix(katydid_data$raw_detections, "common_name", "site")
  
  sites_env <- detection_counts %>%
    mutate(
      site_number = as.numeric(gsub("S", "", site)),
      Plot = sprintf("YB-P%02d", site_number)
    ) %>%
    left_join(vegetation_enhanced, by = "Plot") %>%
    filter(complete.cases(.))
  
  species_matrix <- sites_env %>%
    select(-c(site, site_number, Plot)) %>%
    select_if(is.numeric) %>%
    select(-starts_with("tree_"), -starts_with("liana_"),
           -contains("richness"), -contains("abundance"), 
           -contains("basal")) %>%
    select_if(~ sum(.) >= 2)
  
  if (ncol(species_matrix) < 3 || nrow(sites_env) < 6) {
    cat("[X] Insufficient data for variance partitioning\n")
    return(NULL)
  }
  
  cat(sprintf("   Species: %d, Sites: %d\n", ncol(species_matrix), nrow(sites_env)))
  
  # Groups of environmental variables
  structure_vars <- sites_env %>%
    select(any_of(c("structural_complexity", "log_tree_abundance", "log_liana_stems")))
  
  composition_vars <- sites_env %>%
    select(any_of(c("total_plant_richness", "liana_dominance", "liana_to_tree_ratio")))
  
  if (ncol(structure_vars) < 2 || ncol(composition_vars) < 2) {
    cat("[X] Insufficient environmental variables\n")
    return(NULL)
  }
  
  tryCatch({
    # Variance partitioning
    varpart_result <- varpart(log1p(species_matrix), 
                              structure_vars, 
                              composition_vars)
    
    # Create partition data frame
    partition_df <- data.frame(
      Component = c("Structure", "Composition", "Structure AND Composition", "Residual"),
      Variance = c(
        varpart_result$part$indfract$Adj.R.squared[1],
        varpart_result$part$indfract$Adj.R.squared[2],
        varpart_result$part$indfract$Adj.R.squared[3],
        1 - sum(varpart_result$part$indfract$Adj.R.squared[1:3], na.rm = TRUE)
      )
    ) %>%
      mutate(Variance = pmax(0, Variance))  # Remove negative values
    
    cat("   Variance explained:\n")
    for (i in 1:nrow(partition_df)) {
      cat(sprintf("     %s: %.1f%%\n", partition_df$Component[i], partition_df$Variance[i] * 100))
    }
    
    # Create plot
    partition_plot <- ggplot(partition_df, aes(x = Component, y = Variance, fill = Component)) +
      geom_col() +
      geom_text(aes(label = paste0(round(Variance * 100, 1), "%")), 
                vjust = -0.5, fontface = "bold") +
      scale_fill_viridis_d() +
      labs(
        title = "Variance Partitioning: Katydid Communities",
        subtitle = "Structure vs Composition effects",
        y = "Proportion of Variance Explained",
        x = ""
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 14, face = "bold")
      )
    
    # Save plot
    ggsave(file.path(output_dir, "variance_partitioning.jpg"), 
           partition_plot, width = 10, height = 7, dpi = 300)
    
    # Save detailed results
    capture.output(varpart_result, 
                   file = file.path(output_dir, "variance_partitioning_results.txt"))
    
    cat("[OK] Variance partitioning complete\n")
    cat(sprintf("   Plot saved: variance_partitioning.jpg\n"))
    
    return(list(
      varpart = varpart_result,
      plot = partition_plot,
      partition_data = partition_df
    ))
    
  }, error = function(e) {
    cat(sprintf("[X] Variance partitioning failed: %s\n", e$message))
    return(NULL)
  })
}


#' Run Improved Multivariate Analyses
#'
#' Orchestrates improved multivariate analyses including enhanced environmental 
#' variables, robust NMDS with bootstrap, and variance partitioning.
#'
#' @param all_results List containing katydid_data, bird_data, and vegetation_data
#' @param output_dir Character string. Base output directory path
#'
#' @return List with:
#' \describe{
#'   \item{vegetation_enhanced}{Data frame with enhanced environmental variables}
#'   \item{robust_nmds}{Results from analyze_robust_nmds}
#'   \item{variance_partitioning}{Results from analyze_variance_partitioning}
#' }
#'
#' @details
#' This wrapper function:
#' 1. Creates enhanced environmental variables (composites, transformations)
#' 2. Performs robust NMDS with bootstrap confidence intervals
#' 3. Conducts variance partitioning analysis
#'
#' Results are saved in a subdirectory "improved_multivariate_analyses".
#'
#' @export
run_improved_analyses <- function(all_results, output_dir) {
  
  cat("\n[>>] RUNNING IMPROVED MULTIVARIATE ANALYSES\n")
  cat("============================================\n")
  
  improved_output_dir <- file.path(output_dir, "improved_multivariate_analyses")
  dir.create(improved_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1. Enhanced environmental variables
  cat("\n1. Creating enhanced environmental variables...\n")
  vegetation_enhanced <- tryCatch({
    create_enhanced_environmental_variables(all_results$vegetation_data)
  }, error = function(e) {
    cat(sprintf("   Failed: %s\n", e$message))
    NULL
  })
  
  if (is.null(vegetation_enhanced)) {
    cat("[X] Cannot proceed without enhanced environmental variables\n")
    return(NULL)
  }
  
  # 2. Robust NMDS with bootstrap
  cat("\n2. Robust NMDS with bootstrap...\n")
  robust_nmds <- tryCatch({
    analyze_robust_nmds(
      all_results$katydid_data,
      all_results$bird_data,
      vegetation_enhanced,
      improved_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   Failed: %s\n", e$message))
    NULL
  })
  
  # 3. Variance partitioning
  cat("\n3. Variance partitioning...\n")
  variance_part <- tryCatch({
    analyze_variance_partitioning(
      all_results$katydid_data,
      vegetation_enhanced,
      improved_output_dir
    )
  }, error = function(e) {
    cat(sprintf("   Failed: %s\n", e$message))
    NULL
  })
  
  improved_results <- list(
    vegetation_enhanced = vegetation_enhanced,
    robust_nmds = robust_nmds,
    variance_partitioning = variance_part
  )
  
  cat(sprintf("\n[OK] Improved analyses complete! Results in: %s\n", improved_output_dir))
  
  return(improved_results)
}


# ==============================================================================
# END OF SECTION 12 - ADDITIONAL MISSING FUNCTIONS
# ==============================================================================


# ==============================================================================
# END OF STATISTICAL_ANALYSIS.R
# ==============================================================================
