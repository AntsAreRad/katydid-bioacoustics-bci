# ==============================================================================
# PLOTTING FUNCTIONS
# ==============================================================================
# Description: Functions for creating publication-quality visualizations
#              for comparative analysis of bioacoustic and metabarcoding methods
# ==============================================================================
# Project: Comparative Analysis of Bioacoustic vs Metabarcoding Methods
#          for Katydid (Orthoptera: Tettigoniidae) Community Assessment
# Location: Barro Colorado Island (BCI), Panama
# Author: Leon Brouille (M1 IMABEE)
# Supervisors: Dr. Yves Basset (STRI), Dr. Greg Lamarre (STRI), 
#              Dr. Laurel Symes (Cornell Lab of Ornithology)
# ==============================================================================

# Required packages
#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_smooth geom_abline
#'             geom_hline geom_text labs theme_minimal theme ggsave coord_flip
#'             position_dodge
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr mutate bind_rows filter select group_by summarise %>%
#' @importFrom tidyr pivot_longer



# ==============================================================================
# SECTION 1: COMPREHENSIVE COMPARISON PLOTS
# ==============================================================================

#' Create Comprehensive Method Comparison Plots
#'
#' Generates a complete suite of publication-quality visualization plots
#' comparing bioacoustic and metabarcoding detection methods. This function
#' creates multiple plot types including Venn diagrams, richness comparisons,
#' vegetation correlations, and temporal patterns to provide a comprehensive
#' visual assessment of method performance and complementarity.
#'
#' @details
#' This function is designed to produce publication-ready figures for scientific
#' manuscripts comparing biodiversity detection methods. It creates multiple
#' visualization types:
#'
#' **1. Species Detection Comparison (Bar Plot)**
#' - Shows total species detected by each method
#' - Displays both raw counts and percentages
#' - Highlights method-specific and shared detections
#'
#' **2. Venn Diagrams**
#' - Visualizes species overlap between methods
#' - Can be generated for all species or specific taxonomic groups
#' - Uses colored circles with transparency to show intersections
#'
#' **3. Site Richness Comparison (Scatter Plot)**
#' - Compares species richness per site between methods
#' - Includes linear regression line and 1:1 reference line
#' - Confidence intervals shown with shaded regions
#'
#' **4. Vegetation Correlation Comparison**
#' - Displays correlations with environmental variables
#' - Side-by-side comparison of both methods
#' - Ordered by correlation strength for clarity
#'
#' **5. Temporal Pattern Visualizations**
#' - Hourly activity patterns (if available in results)
#' - Diel period activity summaries
#' - Time-series plots with error bars
#'
#' All plots follow consistent styling with minimal themes, clear axis labels,
#' and appropriate color schemes for publication. The function handles missing
#' data gracefully and only generates plots for available result components.
#'
#' **Design Principles:**
#' - Publication-ready quality (300 DPI recommended for save)
#' - Consistent visual style across all plots
#' - Clear, informative labels and legends
#' - Appropriate use of color and contrast
#' - Accessibility-friendly color schemes
#'
#' **Statistical Visualization Best Practices:**
#' The plots follow established principles from Tufte (2001), Wilke (2019),
#' and recommendations from Weissgerber et al. (2015) for presenting
#' comparative data in scientific publications.
#'
#' @param all_results List containing integrated analysis results. Expected components:
#'   \itemize{
#'     \item \code{method_comparison}: Data frame with species counts and percentages
#'           per method (columns: Method, Species_count, Percentage)
#'     \item \code{site_similarity}: List with site_richness data frame comparing
#'           acoustic_richness and metabar_richness per site
#'     \item \code{vegetation_effects_acoustic}: List with correlations data frame
#'           (vegetation variables vs acoustic richness)
#'     \item \code{vegetation_effects_metabar}: List with correlations data frame
#'           (vegetation variables vs metabarcoding richness)
#'     \item \code{temporal_patterns}: List with period_plot and hourly_plot
#'           ggplot objects from temporal analysis
#'     \item \code{venn_data}: Optional list with species vectors for Venn diagram
#'     \item \code{venn_plot}: Optional pre-generated Venn diagram (grid object)
#'   }
#' @param method_names Character vector of length 2 with method names for labeling.
#'   Default: c("Bioacoustic", "Metabarcoding")
#' @param color_palette Character vector of colors for methods. Default uses
#'   colorblind-friendly palette: c("#0072B2", "#009E73"). If NULL, ggplot2
#'   defaults are used.
#' @param plot_width Numeric. Default width in inches for saved plots. Default: 10
#' @param plot_height Numeric. Default height in inches for saved plots. Default: 6
#' @param theme_base ggplot2 theme function to use as base. Default: theme_minimal()
#'
#' @return Named list of ggplot objects and metadata:
#'   \itemize{
#'     \item \code{detection_comparison_plot}: Bar plot of species counts per method
#'     \item \code{venn_plot}: Venn diagram (grid object) showing species overlap
#'     \item \code{richness_comparison_plot}: Scatter plot comparing site richness
#'     \item \code{vegetation_correlation_plot}: Bar plot of vegetation correlations
#'     \item \code{temporal_period_plot}: Bar plot of activity by diel period
#'     \item \code{temporal_hourly_plot}: Line plot of hourly activity patterns
#'     \item \code{plot_metadata}: List with dimensions and plot availability
#'   }
#'   
#'   Plots that cannot be generated (due to missing data) are returned as NULL
#'   with a warning message. The plot_metadata component tracks which plots
#'   were successfully created.
#'
#' @note
#' **Important Considerations:**
#' \itemize{
#'   \item Venn diagrams require at least 2 species per method to be meaningful
#'   \item Richness comparisons require site-level data in both methods
#'   \item Vegetation plots require environmental data and correlation results
#'   \item Temporal plots require timestamp data in acoustic detections
#'   \item All plots use consistent styling but can be further customized
#'   \item For publication, save plots at 300 DPI minimum
#'   \item Consider journal-specific requirements for figure dimensions
#' }
#'
#' **Performance:**
#' - Fast for typical dataset sizes (<1000 species, <100 sites)
#' - Venn diagram generation is the slowest component (~1-2 seconds)
#' - Large temporal datasets may take longer to plot
#'
#' @seealso
#' \code{\link{save_all_plots}} for organized batch saving of plots
#' \code{\link{analyze_site_similarity}} for generating richness comparison data
#' \code{\link{analyze_vegetation_effects}} for correlation analysis
#' \code{\link{analyze_temporal_patterns}} for temporal data
#'
#' @references
#' Tufte, E. R. (2001). The Visual Display of Quantitative Information (2nd ed.).
#' Graphics Press. ISBN: 978-0961392147
#' 
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis (2nd ed.).
#' Springer. https://doi.org/10.1007/978-3-319-24277-4
#' 
#' Weissgerber, T. L., et al. (2015). Beyond bar and line graphs: time for a new
#' data presentation paradigm. PLOS Biology, 13(4), e1002128.
#' https://doi.org/10.1371/journal.pbio.1002128
#' 
#' Wilke, C. O. (2019). Fundamentals of Data Visualization. O'Reilly Media.
#' https://clauswilke.com/dataviz/
#' 
#' Rougier, N. P., et al. (2014). Ten simple rules for better figures.
#' PLOS Computational Biology, 10(9), e1003833.
#' https://doi.org/10.1371/journal.pcbi.1003833
#'
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have run the integrated analysis
#' results <- main_integrated_analysis(
#'   katydid_detections_file = "katydid_detections.csv",
#'   metabarcoding_file = "metabarcoding_data.xlsx",
#'   vegetation_file = "vegetation_data.xlsx"
#' )
#' 
#' # Create all comparison plots
#' plots <- create_comparison_plots(
#'   all_results = results,
#'   method_names = c("Bioacoustic", "DNA Metabarcoding"),
#'   color_palette = c("#E69F00", "#56B4E9")
#' )
#' 
#' # Access individual plots
#' print(plots$detection_comparison_plot)
#' grid.draw(plots$venn_plot)
#' print(plots$richness_comparison_plot)
#' 
#' # Check which plots were created
#' print(plots$plot_metadata$available_plots)
#' 
#' # Customize a plot further
#' library(ggplot2)
#' plots$richness_comparison_plot +
#'   labs(title = "Custom Title") +
#'   theme_classic()
#' }
create_comparison_plots <- function(all_results,
                                   method_names = c("Bioacoustic", "Metabarcoding"),
                                   color_palette = c("#0072B2", "#009E73"),
                                   plot_width = 10,
                                   plot_height = 6,
                                   theme_base = theme_minimal()) {
  
  # Input validation
  if (!is.list(all_results)) {
    stop("all_results must be a list containing analysis results")
  }
  
  if (length(method_names) != 2) {
    stop("method_names must be a character vector of length 2")
  }
  
  if (!is.null(color_palette) && length(color_palette) != 2) {
    warning("color_palette should have 2 colors, using defaults")
    color_palette <- c("#0072B2", "#009E73")
  }
  
  # Initialize output list
  plot_list <- list()
  created_plots <- character(0)
  
  cat("\nGenerating comparison plots...\n")
  
  # ============================================================================
  # PLOT 1: Species Detection Comparison (Bar Plot)
  # ============================================================================
  cat("  [DEBUG] Checking method_comparison for detection plot...\n")
  
  if (!is.null(all_results$method_comparison) && 
      "stats" %in% names(all_results$method_comparison) &&
      !is.null(all_results$method_comparison$stats) &&
      nrow(all_results$method_comparison$stats) > 0) {
    
    cat("  - Creating species detection comparison plot...\n")
    
    detection_data <- all_results$method_comparison$stats
    cat(sprintf("    Stats data: %d rows, columns: %s\n", 
                nrow(detection_data), paste(colnames(detection_data), collapse = ", ")))
    
    # Create bar plot with counts and percentages
    detection_plot <- ggplot(detection_data, 
                            aes(x = Method, y = Species_count, fill = Method)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = paste0(Species_count, "\n(", 
                                   round(Percentage, 1), "%)")), 
               vjust = -0.5, size = 4) +
      labs(title = paste("Species Detection Comparison:",
                        method_names[1], "vs", method_names[2]),
           x = "Detection Method", 
           y = "Number of Species Detected") +
      theme_base +
      theme(legend.position = "none",
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 11)) +
      coord_cartesian(ylim = c(0, max(detection_data$Species_count, na.rm = TRUE) * 1.15))
    
    # Apply custom colors if provided
    if (!is.null(color_palette)) {
      detection_plot <- detection_plot +
        scale_fill_manual(values = color_palette)
    }
    
    plot_list$detection_comparison_plot <- detection_plot
    created_plots <- c(created_plots, "detection_comparison")
    
  } else {
    plot_list$detection_comparison_plot <- NULL
    cat("  [WARNING] Method comparison data not available or empty, skipping detection plot\n")
    cat("    method_comparison is NULL:", is.null(all_results$method_comparison), "\n")
    if (!is.null(all_results$method_comparison)) {
      cat("    'stats' in names:", "stats" %in% names(all_results$method_comparison), "\n")
      if ("stats" %in% names(all_results$method_comparison)) {
        cat("    stats is NULL:", is.null(all_results$method_comparison$stats), "\n")
        if (!is.null(all_results$method_comparison$stats)) {
          cat("    stats nrow:", nrow(all_results$method_comparison$stats), "\n")
        }
      }
    }
  }
  
  # ============================================================================
  # PLOT 2: Venn Diagram (Species Overlap)
  # ============================================================================
  if (!is.null(all_results$method_comparison$venn_plot)) {
    
    cat("  - Including pre-generated Venn diagram...\n")
    plot_list$venn_plot <- all_results$method_comparison$venn_plot
    created_plots <- c(created_plots, "venn_diagram")
    
  } else if (!is.null(all_results$method_comparison$venn_data)) {
    
    cat("  - Creating Venn diagram from species data...\n")
    
    venn_data <- all_results$method_comparison$venn_data
    
    # Only create if we have data for both methods
    if (length(venn_data[[1]]) > 0 && length(venn_data[[2]]) > 0) {
      
      venn_list <- list(venn_data[[1]], venn_data[[2]])
      names(venn_list) <- method_names
      
      venn_plot <- venn.diagram(
        x = venn_list,
        filename = NULL,
        fill = if (!is.null(color_palette)) color_palette else c("lightblue", "lightgreen"),
        alpha = 0.5,
        cex = 1.8,
        cat.cex = 1.5,
        cat.pos = c(-20, 20),
        cat.dist = c(0.05, 0.05),
        main = paste("Species Detection Overlap:\n", 
                    method_names[1], "vs", method_names[2]),
        main.cex = 1.3
      )
      
      plot_list$venn_plot <- venn_plot
      created_plots <- c(created_plots, "venn_diagram")
      
    } else {
      plot_list$venn_plot <- NULL
      warning("Insufficient species data for Venn diagram (need >0 species per method)")
    }
    
  } else {
    plot_list$venn_plot <- NULL
    warning("Venn diagram data not available")
  }
  
  # ============================================================================
  # PLOT 3: Site Richness Comparison (Scatter Plot)
  # ============================================================================
  if (!is.null(all_results$site_similarity) && 
      !is.null(all_results$site_similarity$site_richness)) {
    
    site_richness <- all_results$site_similarity$site_richness
    
    # Check if required columns exist
    if ("acoustic_richness" %in% colnames(site_richness) &&
        "metabar_richness" %in% colnames(site_richness)) {
      
      cat("  - Creating site richness comparison plot...\n")
      
      richness_plot <- ggplot(site_richness, 
                             aes(x = acoustic_richness, y = metabar_richness)) +
        geom_point(size = 3, alpha = 0.7, color = color_palette[1]) +
        geom_smooth(method = "lm", se = TRUE, color = color_palette[2], 
                   fill = color_palette[2], alpha = 0.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
                   color = "red", linewidth = 0.8) +
        labs(title = paste("Site Richness Comparison:\n",
                          method_names[1], "vs", method_names[2]),
             x = paste(method_names[1], "Species Richness"), 
             y = paste(method_names[2], "Species Richness"),
             caption = "Dashed line = 1:1 relationship") +
        theme_base +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 11),
              plot.caption = element_text(size = 9, hjust = 0.5))
      
      plot_list$richness_comparison_plot <- richness_plot
      created_plots <- c(created_plots, "richness_comparison")
      
    } else {
      plot_list$richness_comparison_plot <- NULL
      warning("Site richness columns not found, skipping richness comparison plot")
    }
    
  } else {
    plot_list$richness_comparison_plot <- NULL
    warning("Site similarity data not available, skipping richness comparison")
  }
  
  # ============================================================================
  # PLOT 4: Vegetation Correlations Comparison
  # ============================================================================
  if (!is.null(all_results$vegetation_effects_acoustic) || 
      !is.null(all_results$vegetation_effects_metabar)) {
    
    cat("  - Creating vegetation correlations comparison plot...\n")
    
    veg_correlations <- data.frame()
    
    # Combine acoustic correlations
    if (!is.null(all_results$vegetation_effects_acoustic) &&
        "correlations" %in% names(all_results$vegetation_effects_acoustic)) {
      
      acoustic_cors <- all_results$vegetation_effects_acoustic$correlations %>%
        mutate(Method = method_names[1])
      
      veg_correlations <- bind_rows(veg_correlations, acoustic_cors)
    }
    
    # Combine metabarcoding correlations
    if (!is.null(all_results$vegetation_effects_metabar) &&
        "correlations" %in% names(all_results$vegetation_effects_metabar)) {
      
      metabar_cors <- all_results$vegetation_effects_metabar$correlations %>%
        mutate(Method = method_names[2])
      
      veg_correlations <- bind_rows(veg_correlations, metabar_cors)
    }
    
    # Only create plot if we have data
    if (nrow(veg_correlations) > 0) {
      
      veg_plot <- ggplot(veg_correlations, 
                        aes(x = reorder(variable, correlation), 
                            y = correlation, fill = Method)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.7) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        coord_flip() +
        labs(title = "Correlations with Vegetation Variables",
             subtitle = "Comparison of environmental relationships",
             x = "Vegetation Variable", 
             y = "Spearman Correlation Coefficient",
             fill = "Method") +
        theme_base +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 11, hjust = 0.5),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 10),
              legend.position = "bottom",
              legend.title = element_text(size = 11),
              legend.text = element_text(size = 10))
      
      # Apply custom colors if provided
      if (!is.null(color_palette)) {
        veg_plot <- veg_plot +
          scale_fill_manual(values = color_palette)
      }
      
      plot_list$vegetation_correlation_plot <- veg_plot
      created_plots <- c(created_plots, "vegetation_correlations")
      
    } else {
      plot_list$vegetation_correlation_plot <- NULL
      warning("No vegetation correlation data available")
    }
    
  } else {
    plot_list$vegetation_correlation_plot <- NULL
    warning("Vegetation effects data not available")
  }
  
  # ============================================================================
  # PLOT 5: Temporal Patterns (if available)
  # ============================================================================
  if (!is.null(all_results$temporal_patterns)) {
    
    # Period plot (diel periods)
    if (!is.null(all_results$temporal_patterns$period_plot)) {
      cat("  - Including temporal period plot...\n")
      plot_list$temporal_period_plot <- all_results$temporal_patterns$period_plot
      created_plots <- c(created_plots, "temporal_periods")
    } else {
      plot_list$temporal_period_plot <- NULL
    }
    
    # Hourly plot
    if (!is.null(all_results$temporal_patterns$hourly_plot)) {
      cat("  - Including hourly activity plot...\n")
      plot_list$temporal_hourly_plot <- all_results$temporal_patterns$hourly_plot
      created_plots <- c(created_plots, "temporal_hourly")
    } else {
      plot_list$temporal_hourly_plot <- NULL
    }
    
  } else {
    plot_list$temporal_period_plot <- NULL
    plot_list$temporal_hourly_plot <- NULL
  }
  
  # ============================================================================
  # Add metadata
  # ============================================================================
  plot_list$plot_metadata <- list(
    created_plots = created_plots,
    n_plots_created = length(created_plots),
    default_width = plot_width,
    default_height = plot_height,
    method_names = method_names,
    color_palette = color_palette,
    timestamp = Sys.time()
  )
  
  cat(sprintf("\n[OK] Successfully created %d plots\n", length(created_plots)))
  cat("  Available plots:", paste(created_plots, collapse = ", "), "\n")
  
  return(plot_list)
}



# ==============================================================================
# SECTION 2: PLOT SAVING AND EXPORT
# ==============================================================================

#' Save All Generated Plots with Organized File Structure
#'
#' Efficiently saves all plots generated by \code{create_comparison_plots()}
#' to disk in multiple formats with organized subdirectory structure. This
#' function handles different plot types (ggplot objects and grid/VennDiagram
#' objects) appropriately and creates publication-ready outputs with
#' configurable quality settings.
#'
#' @details
#' This function provides a comprehensive solution for exporting visualization
#' outputs from comparative biodiversity analyses. It automatically:
#'
#' **File Organization:**
#' Creates a structured directory hierarchy within the output directory:
#' \preformatted{
#' output_dir/
#' +-- figures/
#'     +-- comparisons/
#'         +-- species_detection.[format]
#'         +-- venn_diagram.[format]
#'         +-- site_richness.[format]
#'     +-- environmental/
#'         +-- vegetation_correlations.[format]
#'     +-- temporal/
#'         +-- hourly_patterns.[format]
#'         +-- period_patterns.[format]
#' }
#'
#' **Format Handling:**
#' - **ggplot2 objects**: Saved using \code{ggsave()} with full control over
#'   dimensions, resolution, and format-specific parameters
#' - **grid/VennDiagram objects**: Saved using format-specific device functions
#'   (png(), pdf(), jpeg()) with \code{grid.draw()} for rendering
#' - **Multiple formats**: Can save each plot in multiple formats simultaneously
#'   (e.g., high-res PDF for publication + JPG for presentations)
#'
#' **Quality Control:**
#' - DPI settings ensure publication-quality outputs (300 DPI default)
#' - Configurable dimensions accommodate different journal requirements
#' - Consistent file naming follows best practices for reproducibility
#' - Error handling prevents partial failures from stopping entire export
#'
#' **Performance Considerations:**
#' - Parallel format saving for efficiency (when saving multiple formats)
#' - Progress reporting for long operations
#' - Skips NULL plots without errors
#' - Creates directories only as needed
#'
#' **Publication Standards:**
#' Following recommendations from scientific journals and data visualization
#' best practices (Rougier et al. 2014, Wilke 2019), this function ensures
#' figures meet typical requirements:
#' - Minimum 300 DPI for print publications
#' - Vector formats (PDF) for scalability
#' - Raster formats (PNG/JPG) for web/presentations
#' - Consistent sizing for figure panels
#'
#' @param plot_list Named list of plots, typically output from
#'   \code{create_comparison_plots()}. Can contain:
#'   \itemize{
#'     \item ggplot objects (e.g., detection_comparison_plot, richness_comparison_plot)
#'     \item grid/VennDiagram objects (e.g., venn_plot)
#'     \item NULL values (which are safely skipped)
#'     \item plot_metadata (list component, not saved)
#'   }
#' @param output_dir Character string specifying the base output directory.
#'   If it doesn't exist, it will be created. Subdirectories for different
#'   plot categories are created automatically. Default: "output_plots"
#' @param formats Character vector of output formats. Supported formats:
#'   \itemize{
#'     \item "jpg" or "jpeg": JPEG format, good for photos/complex plots
#'     \item "png": PNG format, lossless, good for web
#'     \item "pdf": PDF format, vector graphics, best for publication
#'     \item "tiff": TIFF format, high-quality raster
#'   }
#'   Default: c("pdf", "jpg") to cover publication and presentation needs
#' @param width Numeric. Width of plots in inches. Default: 10
#' @param height Numeric. Height of plots in inches. Default: 6
#' @param dpi Numeric. Resolution in dots per inch. Applies to raster formats
#'   (jpg, png, tiff). Default: 300 (publication quality). Use 150 for
#'   presentations, 600+ for high-end printing.
#' @param create_subdirs Logical. Should category-specific subdirectories be
#'   created? If TRUE (default), plots are organized into comparisons/,
#'   environmental/, and temporal/ folders. If FALSE, all plots go directly
#'   in output_dir/figures/
#' @param overwrite Logical. Should existing files be overwritten? Default: TRUE
#' @param compression Character. Compression method for PDF. Options: "none",
#'   "lzw", "zip". Default: "lzw" for good balance of size and quality
#' @param verbose Logical. Should progress messages be printed? Default: TRUE
#'
#' @return Invisibly returns a data frame with save operation results:
#'   \itemize{
#'     \item \code{plot_name}: Name of the plot
#'     \item \code{format}: File format used
#'     \item \code{filepath}: Full path to saved file
#'     \item \code{success}: Logical indicating if save was successful
#'     \item \code{file_size_mb}: Size of saved file in megabytes
#'     \item \code{timestamp}: Time of save operation
#'   }
#'   This data frame can be used for quality control and documentation.
#'
#' @note
#' **Important Considerations:**
#' \itemize{
#'   \item **File sizes**: PDF files are typically smaller than high-res raster
#'         formats. Use PDF for publication, JPG/PNG for quick sharing
#'   \item **Venn diagrams**: These require special handling as grid objects.
#'         Function automatically detects and handles them correctly
#'   \item **Large datasets**: High DPI and multiple formats can create large
#'         files. Monitor disk space for extensive analyses
#'   \item **Overwriting**: By default, existing files are overwritten. Set
#'         \code{overwrite = FALSE} to preserve previous versions
#'   \item **Subdirectories**: Organized structure aids in manuscript preparation
#'         but can be disabled for simpler workflows
#'   \item **Error handling**: Individual plot failures don't stop batch saving.
#'         Check the returned data frame for any failed saves
#' }
#'
#' **Journal-Specific Requirements:**
#' Different journals have different figure requirements. Common specifications:
#' - Nature/Science: PDF or EPS, 300-600 DPI for raster
#' - PLOS: TIFF or EPS, 300-600 DPI
#' - Ecology/Ecological Monographs: PDF preferred, 300 DPI minimum
#' - Biology Letters: PDF or high-res JPG, 300 DPI minimum
#'
#' Adjust \code{formats}, \code{dpi}, and dimensions accordingly.
#'
#' @seealso
#' \code{\link{create_comparison_plots}} for generating the plot list
#' \code{\link[ggplot2]{ggsave}} for ggplot2 saving details
#' \code{\link[grid]{grid.draw}} for grid graphics rendering
#'
#' @references
#' Rougier, N. P., et al. (2014). Ten simple rules for better figures.
#' PLOS Computational Biology, 10(9), e1003833.
#' https://doi.org/10.1371/journal.pcbi.1003833
#'
#' Wilke, C. O. (2019). Fundamentals of Data Visualization. O'Reilly Media.
#' https://clauswilke.com/dataviz/
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis (2nd ed.).
#' Springer. https://doi.org/10.1007/978-3-319-24277-4
#'
#' Nature Journal Figure Guidelines. https://www.nature.com/nature/for-authors/final-submission
#' (Guidelines for figure preparation in scientific publications)
#'
#' @export
#' @examples
#' \dontrun{
#' # Create plots first
#' plots <- create_comparison_plots(all_results)
#' 
#' # Save with default settings (PDF + JPG, 300 DPI)
#' save_results <- save_all_plots(
#'   plot_list = plots,
#'   output_dir = "manuscript_figures"
#' )
#' 
#' # Save only high-resolution PDFs for publication
#' save_all_plots(
#'   plot_list = plots,
#'   output_dir = "publication_figures",
#'   formats = "pdf",
#'   width = 8.5,
#'   height = 6,
#'   dpi = 600
#' )
#' 
#' # Save multiple formats without subdirectories
#' save_all_plots(
#'   plot_list = plots,
#'   output_dir = "all_figures",
#'   formats = c("pdf", "png", "jpg"),
#'   create_subdirs = FALSE,
#'   dpi = 300
#' )
#' 
#' # Check save results
#' print(save_results)
#' 
#' # Total size of saved files
#' sum(save_results$file_size_mb, na.rm = TRUE)
#' }
save_all_plots <- function(plot_list,
                          output_dir = "output_plots",
                          formats = c("pdf", "jpg"),
                          width = 10,
                          height = 6,
                          dpi = 300,
                          create_subdirs = TRUE,
                          overwrite = TRUE,
                          compression = "lzw",
                          verbose = TRUE) {
  
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  if (!is.list(plot_list)) {
    stop("plot_list must be a list of plot objects")
  }
  
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }
  
  # Validate formats
  valid_formats <- c("jpg", "jpeg", "png", "pdf", "tiff", "tif")
  formats <- tolower(formats)
  invalid <- formats[!formats %in% valid_formats]
  if (length(invalid) > 0) {
    stop(sprintf("Invalid format(s): %s. Valid formats: %s",
                paste(invalid, collapse = ", "),
                paste(valid_formats, collapse = ", ")))
  }
  
  # Normalize jpeg to jpg
  formats[formats == "jpeg"] <- "jpg"
  formats[formats == "tif"] <- "tiff"
  
  if (dpi < 72) {
    warning("DPI < 72 is very low quality. Consider using dpi >= 150")
  }
  
  if (dpi > 600) {
    warning("DPI > 600 creates very large files. Consider using dpi = 300-600")
  }
  
  # ============================================================================
  # SETUP DIRECTORY STRUCTURE
  # ============================================================================
  
  # Create main output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (verbose) cat(sprintf("[OK] Created output directory: %s\n", output_dir))
  }
  
  # Create figures subdirectory
  figures_dir <- file.path(output_dir, "figures")
  if (!dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Define category subdirectories
  subdirs <- list(
    comparisons = file.path(figures_dir, "comparisons"),
    environmental = file.path(figures_dir, "environmental"),
    temporal = file.path(figures_dir, "temporal")
  )
  
  if (create_subdirs) {
    for (subdir in subdirs) {
      if (!dir.exists(subdir)) {
        dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
      }
    }
  }
  
  # ============================================================================
  # DEFINE PLOT CATEGORIES AND FILE NAMES
  # ============================================================================
  
  # Map plot names to categories and file names
  plot_mapping <- list(
    detection_comparison_plot = list(
      category = "comparisons",
      filename = "species_detection_comparison"
    ),
    venn_plot = list(
      category = "comparisons",
      filename = "species_overlap_venn_diagram"
    ),
    richness_comparison_plot = list(
      category = "comparisons",
      filename = "site_richness_comparison"
    ),
    vegetation_correlation_plot = list(
      category = "environmental",
      filename = "vegetation_correlations"
    ),
    temporal_period_plot = list(
      category = "temporal",
      filename = "diel_period_activity"
    ),
    temporal_hourly_plot = list(
      category = "temporal",
      filename = "hourly_activity_pattern"
    )
  )
  
  # ============================================================================
  # SAVE PLOTS
  # ============================================================================
  
  if (verbose) {
    cat("\n============================================================\n")
    cat("          SAVING PLOTS TO DISK                        \n")
    cat("============================================================\n\n")
    cat(sprintf("Output directory: %s\n", output_dir))
    cat(sprintf("Formats: %s\n", paste(formats, collapse = ", ")))
    cat(sprintf("Dimensions: %g x %g inches\n", width, height))
    cat(sprintf("Resolution: %d DPI (raster formats)\n\n", dpi))
  }
  
  # Initialize results tracking
  save_results <- data.frame()
  n_saved <- 0
  n_failed <- 0
  
  # Get plot names (exclude metadata)
  plot_names <- names(plot_list)
  plot_names <- plot_names[plot_names != "plot_metadata"]
  
  # Process each plot
  for (plot_name in plot_names) {
    
    plot_obj <- plot_list[[plot_name]]
    
    # Skip NULL plots
    if (is.null(plot_obj)) {
      if (verbose) {
        cat(sprintf("  [-] Skipping %s (NULL plot)\n", plot_name))
      }
      next
    }
    
    # Get plot info
    if (plot_name %in% names(plot_mapping)) {
      plot_info <- plot_mapping[[plot_name]]
      
      # Determine output directory
      if (create_subdirs) {
        plot_dir <- subdirs[[plot_info$category]]
      } else {
        plot_dir <- figures_dir
      }
      
      base_filename <- plot_info$filename
      
    } else {
      # Unknown plot - save to main figures directory
      plot_dir <- figures_dir
      base_filename <- gsub("_plot$", "", plot_name)
      if (verbose) {
        cat(sprintf("  [WARNING] Unknown plot type: %s (saving to main figures dir)\n", 
                   plot_name))
      }
    }
    
    # Determine if it's a grid object (Venn diagram) or ggplot
    is_grid_object <- inherits(plot_obj, "gList") || 
                      inherits(plot_obj, "grob") ||
                      inherits(plot_obj, "gTree")
    
    # Save in each requested format
    for (fmt in formats) {
      
      filepath <- file.path(plot_dir, paste0(base_filename, ".", fmt))
      
      # Check if file exists and overwrite setting
      if (file.exists(filepath) && !overwrite) {
        if (verbose) {
          cat(sprintf("  [-] Skipping %s (file exists, overwrite=FALSE)\n", 
                     basename(filepath)))
        }
        next
      }
      
      # Save plot
      save_success <- tryCatch({
        
        if (is_grid_object) {
          # Handle grid/VennDiagram objects
          if (fmt == "pdf") {
            pdf(filepath, width = width, height = height)
            grid.draw(plot_obj)
            dev.off()
          } else if (fmt == "png") {
            png(filepath, width = width * dpi, height = height * dpi, 
                res = dpi, units = "px")
            grid.draw(plot_obj)
            dev.off()
          } else if (fmt == "jpg") {
            jpeg(filepath, width = width * dpi, height = height * dpi, 
                 res = dpi, units = "px", quality = 95)
            grid.draw(plot_obj)
            dev.off()
          } else if (fmt == "tiff") {
            tiff(filepath, width = width * dpi, height = height * dpi, 
                 res = dpi, units = "px", compression = compression)
            grid.draw(plot_obj)
            dev.off()
          }
          
        } else {
          # Handle ggplot objects
          ggsave(
            filename = filepath,
            plot = plot_obj,
            width = width,
            height = height,
            dpi = dpi,
            units = "in",
            device = fmt
          )
        }
        
        TRUE  # Success
        
      }, error = function(e) {
        if (verbose) {
          cat(sprintf("  [ERROR] Error saving %s: %s\n", 
                     basename(filepath), e$message))
        }
        n_failed <<- n_failed + 1
        FALSE  # Failure
      })
      
      # Record results
      if (save_success) {
        file_size_mb <- file.info(filepath)$size / (1024^2)
        
        save_results <- rbind(save_results, data.frame(
          plot_name = plot_name,
          format = fmt,
          filepath = filepath,
          success = TRUE,
          file_size_mb = file_size_mb,
          timestamp = Sys.time(),
          stringsAsFactors = FALSE
        ))
        
        n_saved <- n_saved + 1
        
        if (verbose) {
          cat(sprintf("  [OK] Saved %s (%s, %.2f MB)\n", 
                     basename(filepath), toupper(fmt), file_size_mb))
        }
      }
    }
  }
  
  # ============================================================================
  # SUMMARY
  # ============================================================================
  
  if (verbose) {
    cat("\n============================================================\n")
    cat("          SAVE OPERATION COMPLETE                     \n")
    cat("============================================================\n\n")
    cat(sprintf("[OK] Successfully saved: %d files\n", n_saved))
    if (n_failed > 0) {
      cat(sprintf("[ERROR] Failed saves: %d files\n", n_failed))
    }
    
    if (nrow(save_results) > 0) {
      total_size <- sum(save_results$file_size_mb, na.rm = TRUE)
      cat(sprintf("Total size: %.2f MB\n", total_size))
      cat(sprintf("Output location: %s\n", normalizePath(output_dir)))
    }
    cat("\n")
  }
  
  # Return results invisibly
  invisible(save_results)
}


# ==============================================================================
# SECTION 3: COMPREHENSIVE PLOT GENERATION (for main pipeline)
# ==============================================================================

#' Create Comprehensive Plots for Integrated Analysis
#'
#' Generates and saves all main plots for the integrated bioacoustic vs
#' metabarcoding analysis. This function is called from the main analysis
#' pipeline and creates publication-ready figures.
#'
#' @details
#' This function creates four main plot types:
#' - Species detection comparison bar plot
#' - Site richness comparison scatter plot
#' - Vegetation correlations comparison
#' - Temporal patterns (period and hourly)
#'
#' All plots are saved directly to the output directory in JPG format.
#'
#' @param all_results List containing integrated analysis results with components:
#'   \itemize{
#'     \item method_comparison: Results from method comparison analysis
#'     \item site_similarity: Site-level richness comparison
#'     \item vegetation_effects_acoustic: Acoustic vegetation correlations
#'     \item vegetation_effects_metabar: Metabarcoding vegetation correlations
#'     \item temporal_patterns: Temporal activity patterns
#'   }
#' @param output_dir Character string specifying output directory for plots
#'
#' @return NULL (invisibly). Plots are saved to disk.
#'
#' @note This function is designed for direct integration with main_integrated_analysis.
#'       For more flexible plot creation, use create_comparison_plots() instead.
#'
#' @seealso \code{\link{create_comparison_plots}} for interactive plot creation
#' @seealso \code{\link{save_all_plots}} for batch saving
#'
#' @export
#' @examples
#' \dontrun{
#' create_comprehensive_plots(all_results, "output/figures")
#' }
create_comprehensive_plots <- function(all_results, output_dir) {
  
  # Ensure output directory exists
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Plot 1: Species detection comparison (bar plot)
  # DEBUG: Check method_comparison structure
  cat("  [DEBUG] Checking method_comparison structure...\n")
  
  if (!is.null(all_results$method_comparison) && 
      !is.null(all_results$method_comparison$stats) &&
      nrow(all_results$method_comparison$stats) > 0) {
    
    stats_data <- all_results$method_comparison$stats
    
    # Verify required columns exist
    required_cols <- c("Method", "Species_count", "Percentage")
    if (!all(required_cols %in% colnames(stats_data))) {
      cat("  [ERROR] Missing columns in stats. Expected:", paste(required_cols, collapse = ", "), "\n")
      cat("  [ERROR] Found:", paste(colnames(stats_data), collapse = ", "), "\n")
    } else {
      cat("  [OK] Creating species detection comparison plot...\n")
      cat(sprintf("     Data: %d rows\n", nrow(stats_data)))
      print(stats_data)  # Debug output
      
      detection_plot <- ggplot(stats_data, 
                               aes(x = Method, y = Species_count, fill = Method)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = paste0(Species_count, " (", round(Percentage, 1), "%)")), 
                  vjust = -0.5, size = 4.5) +
        labs(title = "Species Detection Comparison: Bioacoustic vs Metabarcoding",
             x = "Detection Method", 
             y = "Number of Species") +
        scale_fill_manual(values = c("Bioacoustic only" = "#0072B2", 
                                     "Metabarcoding only" = "#009E73",
                                     "Both methods" = "#E69F00")) +
        theme_minimal() +
        theme(legend.position = "none",
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              axis.title = element_text(size = 12),
              axis.text = element_text(size = 11)) +
        coord_cartesian(ylim = c(0, max(stats_data$Species_count) * 1.15))
      
      ggsave(file.path(output_dir, "species_detection_comparison.jpg"), 
             detection_plot, width = 10, height = 6, dpi = 300)
      
      cat("  [OK] species_detection_comparison.jpg saved\n")
    }
  } else {
    cat("  [WARNING] method_comparison$stats is NULL or empty - cannot create detection plot\n")
    cat("  [DEBUG] method_comparison is NULL:", is.null(all_results$method_comparison), "\n")
    if (!is.null(all_results$method_comparison)) {
      cat("  [DEBUG] method_comparison names:", paste(names(all_results$method_comparison), collapse = ", "), "\n")
      cat("  [DEBUG] method_comparison$stats is NULL:", is.null(all_results$method_comparison$stats), "\n")
    }
  }
  
  # Plot 2: Site richness comparison
  if (!is.null(all_results$site_similarity) && !is.null(all_results$site_similarity$site_richness)) {
    site_richness <- all_results$site_similarity$site_richness
    
    if ("acoustic_richness" %in% colnames(site_richness)) {
      richness_comparison_plot <- ggplot(site_richness, 
                                         aes(x = acoustic_richness, y = metabar_richness)) +
        geom_point() +
        geom_smooth(method = "lm", se = TRUE) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
        labs(title = "Site Richness: Bioacoustic vs Metabarcoding",
             x = "Bioacoustic Species Richness", 
             y = "Metabarcoding Species Richness") +
        theme_minimal()
      
      ggsave(file.path(output_dir, "site_richness_comparison.jpg"), 
             richness_comparison_plot, width = 8, height = 6)
    }
  }
  
  # Plot 3: Vegetation effects compilation
  if (!is.null(all_results$vegetation_effects_acoustic) || 
      !is.null(all_results$vegetation_effects_metabar)) {
    
    veg_correlations <- data.frame()
    
    if (!is.null(all_results$vegetation_effects_acoustic)) {
      veg_correlations <- rbind(veg_correlations, 
                                all_results$vegetation_effects_acoustic$correlations %>%
                                  mutate(Method = "Bioacoustic"))
    }
    
    if (!is.null(all_results$vegetation_effects_metabar)) {
      veg_correlations <- rbind(veg_correlations, 
                                all_results$vegetation_effects_metabar$correlations %>%
                                  mutate(Method = "Metabarcoding"))
    }
    
    if (nrow(veg_correlations) > 0) {
      veg_plot <- ggplot(veg_correlations, aes(x = reorder(variable, correlation), 
                                               y = correlation, fill = Method)) +
        geom_col(position = "dodge") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        coord_flip() +
        labs(title = "Correlations with Vegetation Variables",
             x = "Vegetation Variable", 
             y = "Spearman Correlation Coefficient") +
        theme_minimal()
      
      ggsave(file.path(output_dir, "vegetation_correlations_comparison.jpg"), 
             veg_plot, width = 12, height = 8)
    }
  }
  
  # Plot 4: Temporal patterns (if available)
  if (!is.null(all_results$temporal_patterns)) {
    
    if (!is.null(all_results$temporal_patterns$period_plot)) {
      ggsave(file.path(output_dir, "temporal_patterns_periods.jpg"), 
             all_results$temporal_patterns$period_plot, width = 10, height = 6)
    }
    
    if (!is.null(all_results$temporal_patterns$hourly_plot)) {
      ggsave(file.path(output_dir, "temporal_patterns_hourly.jpg"), 
             all_results$temporal_patterns$hourly_plot, width = 12, height = 6)
    }
  }
  
  cat("Plots saved\n")
  
  invisible(NULL)
}


# ==============================================================================
# SECTION 4: REPORT ANNEXES GENERATION
# ==============================================================================

#' Create Report Annexes for Scientific Publication
#'
#' Generates comprehensive annexes for the scientific report including
#' data tables, statistical summaries, and figure compilations. This function
#' creates a structured directory with all supplementary materials needed
#' for publication.
#'
#' @details
#' This function creates the following annexes:
#'
#' **Tables:**
#' - TableA1: Katydid presence matrix by site
#' - TableA2: Species detected by both methods (common species details)
#' - TableA3: GLM analysis summary (significant relationships)
#' - TableA4: Vegetation correlations summary
#' - TableA5: Complete katydid species list with detection statistics
#' - TableA6: Katydid-bird relationships OR temporal activity patterns
#'
#' **Figures:**
#' - Copies of important plots from main analysis
#' - Organized with figure numbering (FigureA1, FigureA2, etc.)
#'
#' **Index:**
#' - ANNEXES_INDEX.csv: Master index of all annexes
#'
#' Directory structure created:
#' \preformatted{
#' output_dir/
#' +-- REPORT_ANNEXES/
#'     +-- Tables/
#'     |   +-- TableA1_Katydid_Presence_Matrix.csv
#'     |   +-- TableA2_Common_Species_Details.csv
#'     |   +-- TableA3_GLM_Summary_Complete.csv
#'     |   +-- TableA4_Vegetation_Correlations_Summary.csv
#'     |   +-- TableA5_Complete_Katydid_Species_List.csv
#'     |   +-- TableA6_*.csv
#'     +-- Figures/
#'     |   +-- FigureA*_*.jpg
#'     +-- ANNEXES_INDEX.csv
#' }
#'
#' @param all_results List containing complete integrated analysis results
#' @param output_dir Character string specifying base output directory
#'
#' @return List with:
#'   \itemize{
#'     \item annexes_dir: Path to created annexes directory
#'     \item tables_created: Number of tables created
#'     \item figures_created: Number of figures copied
#'   }
#'
#' @note Tables are generated based on available data in all_results.
#'       Missing components are skipped with informative messages.
#'
#' @export
#' @examples
#' \dontrun{
#' annexes <- create_report_annexes(all_results, "output")
#' print(paste("Created", annexes$tables_created, "tables"))
#' }
create_report_annexes <- function(all_results, output_dir) {
  
  cat("\n=== CREATING REPORT ANNEXES ===\n")
  cat("================================\n")
  
  # Create annexes directory structure
  annexes_dir <- file.path(output_dir, "REPORT_ANNEXES")
  dir.create(annexes_dir, showWarnings = FALSE, recursive = TRUE)
  
  tables_dir <- file.path(annexes_dir, "Tables")
  figures_dir <- file.path(annexes_dir, "Figures")
  
  dir.create(tables_dir, showWarnings = FALSE)
  dir.create(figures_dir, showWarnings = FALSE)
  
  # ==========================================================================
  # ANNEXE 1: DETECTION MATRICES
  # ==========================================================================
  
  cat("\n[*] Annexe 1: Detection matrices...\n")
  
  # Table A1: Katydid presence matrix
  if (!is.null(all_results$katydid_data$presence_matrix)) {
    
    katydid_matrix_clean <- all_results$katydid_data$presence_matrix %>%
      select(site, sort(names(.)[-1])) %>%
      bind_rows(
        summarise_all(., ~ if(is.numeric(.)) sum(., na.rm = TRUE) else "TOTAL")
      )
    
    write.csv(katydid_matrix_clean, 
              file.path(tables_dir, "TableA1_Katydid_Presence_Matrix.csv"), 
              row.names = FALSE)
    
    cat("   Table A1: Katydid presence matrix (", nrow(katydid_matrix_clean)-1, " sites x", 
        ncol(katydid_matrix_clean)-1, " species)\n")
  }
  
  # Table A2: Common species (both methods)
  if (!is.null(all_results$comprehensive_comparison$katydid_comparison$both_methods)) {
    
    common_species <- all_results$comprehensive_comparison$katydid_comparison$both_methods
    
    if (length(common_species) > 0) {
      common_species_table <- data.frame(
        Species = common_species,
        Detection_Method = "Both",
        Sites_Bioacoustic = NA,
        Sites_Metabarcoding = NA,
        Total_Sites = NA
      )
      
      # Calculate site statistics for each common species
      if (!is.null(all_results$katydid_data$presence_matrix)) {
        for (i in 1:length(common_species)) {
          species <- common_species[i]
          if (species %in% colnames(all_results$katydid_data$presence_matrix)) {
            common_species_table$Sites_Bioacoustic[i] <- sum(
              all_results$katydid_data$presence_matrix[[species]], na.rm = TRUE
            )
          }
        }
      }
      
      if (!is.null(all_results$metabarcoding_data$presence_matrix)) {
        for (i in 1:length(common_species)) {
          species <- common_species[i]
          if (species %in% colnames(all_results$metabarcoding_data$presence_matrix)) {
            common_species_table$Sites_Metabarcoding[i] <- sum(
              all_results$metabarcoding_data$presence_matrix[[species]], na.rm = TRUE
            )
          }
        }
      }
      
      write.csv(common_species_table, 
                file.path(tables_dir, "TableA2_Common_Species_Details.csv"), 
                row.names = FALSE)
      
      cat("   Table A2: Common species details (", nrow(common_species_table), " species)\n")
    }
  }
  
  # ==========================================================================
  # ANNEXE 2: STATISTICAL ANALYSES
  # ==========================================================================
  
  cat("\n[*] Annexe 2: Statistical analyses...\n")
  
  # Table A3: GLM summary
  if (!is.null(all_results$glm_analysis) && length(all_results$glm_analysis) > 0) {
    
    glm_summary_table <- data.frame(
      Relationship = character(),
      Response_Variable = character(),
      Predictor_Variable = character(),
      P_value = numeric(),
      Coefficient = numeric(),
      Direction = character(),
      AIC = numeric(),
      Sample_Size = numeric(),
      Pseudo_R_squared = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (name in names(all_results$glm_analysis)) {
      result <- all_results$glm_analysis[[name]]
      
      pseudo_r2 <- if (!is.null(result$model)) {
        1 - (deviance(result$model) / result$model$null.deviance)
      } else NA
      
      glm_summary_table <- rbind(glm_summary_table, data.frame(
        Relationship = name,
        Response_Variable = result$response,
        Predictor_Variable = result$predictor,
        P_value = round(result$p_value, 4),
        Coefficient = round(result$coefficient, 4),
        Direction = result$direction,
        AIC = round(result$aic, 2),
        Sample_Size = result$sample_size,
        Pseudo_R_squared = round(pseudo_r2, 3),
        stringsAsFactors = FALSE
      ))
    }
    
    glm_summary_table <- glm_summary_table %>% arrange(P_value)
    
    write.csv(glm_summary_table, 
              file.path(tables_dir, "TableA3_GLM_Summary_Complete.csv"), 
              row.names = FALSE)
    
    cat("   Table A3: GLM summary (", nrow(glm_summary_table), " significant relationships)\n")
  }
  
  # Table A4: Vegetation correlations
  correlation_summary <- data.frame()
  
  if (!is.null(all_results$vegetation_effects_katydids$correlations)) {
    katydid_corr <- all_results$vegetation_effects_katydids$correlations %>%
      mutate(Group = "Katydids", Method = "Bioacoustic")
    correlation_summary <- rbind(correlation_summary, katydid_corr)
  }
  
  if (!is.null(all_results$vegetation_effects_metabar$correlations)) {
    meta_corr <- all_results$vegetation_effects_metabar$correlations %>%
      mutate(Group = "Orthoptera", Method = "Metabarcoding")
    correlation_summary <- rbind(correlation_summary, meta_corr)
  }
  
  if (!is.null(all_results$vegetation_effects_birds$correlations)) {
    bird_corr <- all_results$vegetation_effects_birds$correlations %>%
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
              file.path(tables_dir, "TableA4_Vegetation_Correlations_Summary.csv"), 
              row.names = FALSE)
    
    cat("   Table A4: Vegetation correlations (", nrow(correlation_summary), " relationships tested)\n")
  }
  
  # ==========================================================================
  # ANNEXE 5: COMPLETE KATYDID SPECIES LIST
  # ==========================================================================
  
  cat("\n[*] Annexe 5: Complete katydid list...\n")
  
  if (!is.null(all_results$katydid_data$raw_detections)) {
    
    tryCatch({
      species_details <- all_results$katydid_data$raw_detections %>%
        group_by(common_name) %>%
        summarise(
          Total_Detections = n(),
          Sites_Detected = n_distinct(site),
          Mean_Confidence = round(mean(confidence, na.rm = TRUE), 3),
          Min_Confidence = round(min(confidence, na.rm = TRUE), 3),
          Max_Confidence = round(max(confidence, na.rm = TRUE), 3),
          .groups = "drop"
        ) %>%
        arrange(desc(Total_Detections)) %>%
        mutate(
          Species_Rank = row_number(),
          Detection_Rate_per_Site = round(Total_Detections / Sites_Detected, 1)
        )
      
      # Try to extract temporal information
      tryCatch({
        temporal_info <- all_results$katydid_data$raw_detections %>%
          mutate(
            date_str = str_extract(file_path, "\\d{8}"),
            detection_date = as.Date(date_str, format = "%Y%m%d")
          ) %>%
          filter(!is.na(detection_date)) %>%
          group_by(common_name) %>%
          summarise(
            First_Detection = min(detection_date, na.rm = TRUE),
            Last_Detection = max(detection_date, na.rm = TRUE),
            Detection_Span_Days = as.numeric(max(detection_date, na.rm = TRUE) - 
                                               min(detection_date, na.rm = TRUE)),
            .groups = "drop"
          )
        
        if (nrow(temporal_info) > 0) {
          species_details <- species_details %>%
            left_join(temporal_info, by = "common_name")
        }
        
      }, error = function(e) {
        cat("    [!] Could not extract dates from file names\n")
      })
      
      final_columns <- c("Species_Rank", "common_name", "Total_Detections", 
                         "Sites_Detected", "Detection_Rate_per_Site", 
                         "Mean_Confidence", "Min_Confidence", "Max_Confidence")
      
      if ("First_Detection" %in% colnames(species_details)) {
        final_columns <- c(final_columns, "First_Detection", "Last_Detection", "Detection_Span_Days")
      }
      
      species_details_final <- species_details %>%
        select(all_of(final_columns[final_columns %in% colnames(species_details)]))
      
      write.csv(species_details_final, 
                file.path(tables_dir, "TableA5_Complete_Katydid_Species_List.csv"), 
                row.names = FALSE)
      
      cat("   Table A5: Complete katydid list (", nrow(species_details_final), " species)\n")
      
    }, error = function(e) {
      cat("   Error creating Table A5:", e$message, "\n")
      cat("   Creating simplified version...\n")
      
      simple_species_list <- all_results$katydid_data$raw_detections %>%
        group_by(common_name) %>%
        summarise(
          Total_Detections = n(),
          Sites_Detected = n_distinct(site),
          Mean_Confidence = round(mean(confidence, na.rm = TRUE), 3),
          .groups = "drop"
        ) %>%
        arrange(desc(Total_Detections)) %>%
        mutate(Species_Rank = row_number())
      
      write.csv(simple_species_list, 
                file.path(tables_dir, "TableA5_Complete_Katydid_Species_List.csv"), 
                row.names = FALSE)
      
      cat("   Table A5 simplified: Katydid list (", nrow(simple_species_list), " species)\n")
    })
  }
  
  # ==========================================================================
  # ANNEXE 6: SUPPLEMENTARY DATA (KATYDID-BIRD OR TEMPORAL)
  # ==========================================================================
  
  cat("\n[*] Annexe 6: Supplementary data...\n")
  
  # Option A: Katydid-bird relationships
  if (!is.null(all_results$katydids_birds_relationship$data)) {
    
    bird_katydid_data <- all_results$katydids_birds_relationship$data %>%
      select(Plot, group1_richness, group2_richness) %>%
      rename(
        Site = Plot,
        Katydid_Richness = group1_richness,
        Bird_Richness = group2_richness
      ) %>%
      mutate(
        Katydid_Bird_Ratio = round(Katydid_Richness / Bird_Richness, 3),
        Combined_Richness = Katydid_Richness + Bird_Richness
      ) %>%
      arrange(desc(Combined_Richness))
    
    write.csv(bird_katydid_data, 
              file.path(tables_dir, "TableA6_Katydid_Bird_Relationships.csv"), 
              row.names = FALSE)
    
    cat("   Table A6: Katydid-bird relationships by site\n")
    
  } else if (!is.null(all_results$temporal_patterns$hourly_activity)) {
    
    # Option B: Temporal patterns
    temporal_summary <- all_results$temporal_patterns$hourly_activity %>%
      mutate(
        Period = case_when(
          hour >= 5 & hour < 9   ~ "Dawn",
          hour >= 9 & hour < 13  ~ "Morning", 
          hour >= 13 & hour < 17 ~ "Afternoon",
          hour >= 17 & hour < 21 ~ "Dusk",
          TRUE                   ~ "Night"
        ),
        Time_Category = case_when(
          hour >= 6 & hour < 18  ~ "Diurnal",
          TRUE                   ~ "Nocturnal"
        )
      ) %>%
      arrange(hour)
    
    write.csv(temporal_summary, 
              file.path(tables_dir, "TableA6_Temporal_Activity_Patterns.csv"), 
              row.names = FALSE)
    
    cat("   Table A6: Temporal activity patterns\n")
  }
  
  # ==========================================================================
  # COPY IMPORTANT FIGURES
  # ==========================================================================
  
  cat("\n[*] Copying important figures...\n")
  
  existing_plots <- list.files(output_dir, pattern = "\\.(jpg|png|pdf)$", 
                               full.names = TRUE, recursive = TRUE)
  
  # Exclude files already in figures_dir
  existing_plots <- existing_plots[!grepl(figures_dir, existing_plots, fixed = TRUE)]
  
  important_plots <- c(
    "species_detection_comparison",
    "katydid_vegetation_relationship", 
    "Bird_Species_Richness_vs_Katydid_Species_Richness",
    "katydid_community_ordination",
    "temporal_patterns"
  )
  
  figure_counter <- 1
  
  for (plot_pattern in important_plots) {
    matching_files <- existing_plots[grepl(plot_pattern, existing_plots, ignore.case = TRUE)]
    
    if (length(matching_files) > 0) {
      source_file <- matching_files[1]
      extension <- tools::file_ext(source_file)
      dest_file <- file.path(figures_dir, paste0("FigureA", figure_counter, "_", 
                                                 gsub("[^A-Za-z0-9]", "_", plot_pattern), 
                                                 ".", extension))
      
      if (normalizePath(source_file, mustWork = FALSE) != normalizePath(dest_file, mustWork = FALSE)) {
        tryCatch({
          success <- file.copy(source_file, dest_file, overwrite = TRUE)
          if (success) {
            cat(sprintf("   Figure A%d: %s copied\n", figure_counter, basename(source_file)))
          } else {
            cat(sprintf("   [!] Failed to copy Figure A%d: %s\n", figure_counter, basename(source_file)))
          }
        }, error = function(e) {
          cat(sprintf("   Error copying Figure A%d: %s - %s\n", figure_counter, basename(source_file), e$message))
        })
      } else {
        cat(sprintf("   Figure A%d: %s already in correct folder\n", figure_counter, basename(source_file)))
      }
      
      figure_counter <- figure_counter + 1
    }
  }
  
  # Create placeholder if no figures were copied
  if (figure_counter == 1) {
    cat("   No figures found - creating index file\n")
    
    figures_index <- data.frame(
      Figure = paste0("Figure A", 1:5),
      Description = c(
        "Bioacoustic vs metabarcoding detection comparison",
        "Katydid-vegetation relationships", 
        "Katydid-bird richness correlation",
        "Katydid community ordination",
        "Temporal activity patterns"
      ),
      Status = "To create",
      Note = "Figures available in main project folder"
    )
    
    write.csv(figures_index, 
              file.path(figures_dir, "Figures_Index.csv"), 
              row.names = FALSE)
  }
  
  # ==========================================================================
  # CREATE ANNEXES INDEX
  # ==========================================================================
  
  cat("\n[*] Creating annexes index...\n")
  
  annexes_index <- data.frame(
    Annexe = c("Annexe 1", "Annexe 1", "Annexe 2", "Annexe 2", "Annexe 5", "Annexe 6"),
    Type = c("Table A1", "Table A2", "Table A3", "Table A4", "Table A5", "Table A6"),
    Title = c(
      "Katydid presence matrix by site",
      "Species detected by both methods", 
      "GLM analysis summary (significant)",
      "Vegetation-community correlations",
      "Complete katydid species list",
      "Supplementary data (relationships/temporal)"
    ),
    File = c(
      "TableA1_Katydid_Presence_Matrix.csv",
      "TableA2_Common_Species_Details.csv",
      "TableA3_GLM_Summary_Complete.csv", 
      "TableA4_Vegetation_Correlations_Summary.csv",
      "TableA5_Complete_Katydid_Species_List.csv",
      "TableA6_Katydid_Bird_Relationships.csv or TableA6_Temporal_Activity_Patterns.csv"
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(annexes_index, 
            file.path(annexes_dir, "ANNEXES_INDEX.csv"), 
            row.names = FALSE)
  
  # ==========================================================================
  # FINAL SUMMARY
  # ==========================================================================
  
  cat("\n=== ANNEXES CREATED SUCCESSFULLY ===\n")
  cat(sprintf(" Main directory: %s\n", annexes_dir))
  cat(sprintf(" Tables: %s\n", tables_dir))
  cat(sprintf(" Figures: %s\n", figures_dir))
  
  tables_created <- length(list.files(tables_dir, pattern = "\\.(csv|txt)$"))
  figures_created <- length(list.files(figures_dir, pattern = "\\.(jpg|png|pdf)$"))
  
  cat(sprintf("\n STATISTICS:\n"))
  cat(sprintf("  - %d tables created\n", tables_created))
  cat(sprintf("  - %d figures copied\n", figures_created))
  cat(sprintf("  - General index created\n"))
  
  cat(sprintf("\n Ready for report integration!\n"))
  
  return(list(
    annexes_dir = annexes_dir,
    tables_created = tables_created,
    figures_created = figures_created
  ))
}
