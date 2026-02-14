# internship_method_comparison.R
# Methodological comparison for M1 IMABEE internship report
# Venn diagram, coinertia (ade4), family analysis, taxonomic diversity, site similarity
#
# Reads from: integrated_results/ (pipeline outputs)
# Writes to:  results/internship_m1/method_comparison/
# Sources:    baseline_helpers.R

run_internship_method_comparison <- function(results, m1_dir) {

  cat("\n[INTERNSHIP-M1] Method comparison analysis...\n")

  out_dir <- file.path(m1_dir, "method_comparison")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  katydid_data      <- results$katydid_data
  metabarcoding_data <- results$metabarcoding_data

  if (is.null(katydid_data) || is.null(metabarcoding_data)) {
    cat("  [!] Missing katydid or metabarcoding data, skipping\n")
    return(NULL)
  }

  comparison_results <- list()

  # ---- 1. Species overlap & Venn diagram ----
  cat("  1. Species overlap (Venn diagram)...\n")

  acoustic_species <- setdiff(colnames(katydid_data$presence_matrix), "site")
  metabar_species  <- metabarcoding_data$species_list

  both_methods   <- intersect(acoustic_species, metabar_species)
  only_acoustic  <- setdiff(acoustic_species, metabar_species)
  only_metabar   <- setdiff(metabar_species, acoustic_species)
  total_species  <- length(unique(c(acoustic_species, metabar_species)))

  comparison_results$species_overlap <- list(
    both_methods  = both_methods,
    only_acoustic = only_acoustic,
    only_metabar  = only_metabar,
    n_acoustic    = length(acoustic_species),
    n_metabar     = length(metabar_species),
    n_both        = length(both_methods),
    overlap_pct   = round(length(both_methods) / total_species * 100, 1)
  )

  # Venn diagram
  tryCatch({
    venn <- VennDiagram::venn.diagram(
      x = list(
        Bioacoustic    = acoustic_species,
        Metabarcoding  = metabar_species
      ),
      filename = NULL,
      fill     = c("#E41A1C", "#377EB8"),
      alpha    = 0.5,
      cat.cex  = 1.2,
      cex      = 1.5,
      main     = "Species Detection Comparison"
    )
    png(file.path(out_dir, "katydids_venn_diagram.png"), width = 800, height = 600)
    grid::grid.draw(venn)
    dev.off()
    cat("    Saved: katydids_venn_diagram.png\n")
  }, error = function(e) {
    cat(sprintf("    [!] Venn diagram failed: %s\n", e$message))
  })

  # Species overlap barplot (Figure 1 in report)
  overlap_df <- data.frame(
    Category = c("Bioacoustic only", "Both methods", "Metabarcoding only"),
    Count    = c(length(only_acoustic), length(both_methods), length(only_metabar)),
    Pct      = c(
      round(length(only_acoustic) / total_species * 100, 1),
      round(length(both_methods)  / total_species * 100, 1),
      round(length(only_metabar)  / total_species * 100, 1)
    )
  )
  overlap_df$Category <- factor(overlap_df$Category,
                                levels = c("Bioacoustic only", "Both methods", "Metabarcoding only"))

  overlap_plot <- ggplot2::ggplot(overlap_df,
                                  ggplot2::aes(x = Category, y = Count, fill = Category)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d (%.1f%%)", Count, Pct)),
                       vjust = -0.5, size = 4) +
    ggplot2::scale_fill_manual(values = c("#E41A1C", "#4DAF4A", "#377EB8")) +
    ggplot2::labs(
      title = "Species Detection Comparison: Bioacoustic vs Metabarcoding",
      x = "Detection Method", y = "Number of Species"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, max(overlap_df$Count) * 1.15)

  ggplot2::ggsave(file.path(out_dir, "species_detection_comparison.jpg"),
                  overlap_plot, width = 10, height = 7, dpi = 300)
  cat("    Saved: species_detection_comparison.jpg\n")

  # ---- 2. Site-level richness comparison ----
  cat("  2. Site-level richness comparison...\n")

  acoustic_richness <- richness_from_matrix(katydid_data$presence_matrix)
  metabar_richness  <- richness_from_matrix(metabarcoding_data$presence_matrix)

  site_richness <- dplyr::inner_join(
    acoustic_richness %>% dplyr::rename(acoustic_richness = richness),
    metabar_richness  %>% dplyr::rename(metabar_richness  = richness),
    by = "site"
  )

  if (nrow(site_richness) >= 3) {
    cor_test <- cor.test(site_richness$acoustic_richness,
                         site_richness$metabar_richness,
                         method = "spearman")
    comparison_results$richness_correlation <- cor_test
    cat(sprintf("    Spearman rho = %.3f, p = %.4f\n",
                cor_test$estimate, cor_test$p.value))
  }

  write.csv(site_richness,
            file.path(out_dir, "site_richness_comparison.csv"), row.names = FALSE)

  # ---- 3. Coinertia analysis (ade4) ----
  cat("  3. Coinertia analysis...\n")

  tryCatch({
    common_sites <- intersect(katydid_data$presence_matrix$site,
                              metabarcoding_data$presence_matrix$site)

    if (length(common_sites) >= 5) {
      acoustic_mat <- katydid_data$presence_matrix %>%
        dplyr::filter(site %in% common_sites) %>%
        dplyr::arrange(site)
      metabar_mat <- metabarcoding_data$presence_matrix %>%
        dplyr::filter(site %in% common_sites) %>%
        dplyr::arrange(site)

      # Remove zero-variance columns
      ac_num <- acoustic_mat %>% dplyr::select(-site) %>% as.data.frame()
      mb_num <- metabar_mat  %>% dplyr::select(-site) %>% as.data.frame()
      ac_num <- ac_num[, apply(ac_num, 2, var) > 0, drop = FALSE]
      mb_num <- mb_num[, apply(mb_num, 2, var) > 0, drop = FALSE]

      if (ncol(ac_num) >= 2 && ncol(mb_num) >= 2) {
        pca_ac <- ade4::dudi.pca(ac_num, scannf = FALSE, nf = 2)
        pca_mb <- ade4::dudi.pca(mb_num, scannf = FALSE, nf = 2)

        coin <- ade4::coinertia(pca_ac, pca_mb, scannf = FALSE, nf = 2)
        coin_test <- ade4::randtest(coin, nrepet = 999)

        comparison_results$coinertia <- list(
          analysis = coin,
          test     = coin_test,
          rv       = coin$RV
        )

        cat(sprintf("    RV coefficient = %.3f, p = %.4f\n",
                    coin$RV, coin_test$pvalue))

        # Save coinertia results
        capture.output(summary(coin), file = file.path(out_dir, "coinertia_summary.txt"))
        capture.output(coin_test, file = file.path(out_dir, "coinertia_test.txt"))
      }
    } else {
      cat("    [!] Insufficient common sites for coinertia\n")
    }
  }, error = function(e) {
    cat(sprintf("    [!] Coinertia failed: %s\n", e$message))
  })

  # ---- 4. Family-level analysis ----
  cat("  4. Family-level detection analysis...\n")

  tryCatch({
    # Try to load trait matrix for family info
    trait_file <- "data/Species_level_trait_matrix_July_2018.csv"
    if (file.exists(trait_file)) {
      traits <- read.csv(trait_file, stringsAsFactors = FALSE)
      family_mapping <- stats::setNames(traits$Family, traits$Species)
    } else {
      family_mapping <- NULL
    }

    # Count by family if mapping available
    if (!is.null(family_mapping)) {
      acoustic_families <- table(family_mapping[acoustic_species])
      metabar_families  <- table(family_mapping[metabar_species])

      family_df <- data.frame(
        Family    = union(names(acoustic_families), names(metabar_families)),
        Bioacoustic    = 0,
        Metabarcoding  = 0
      )
      for (i in seq_len(nrow(family_df))) {
        fam <- family_df$Family[i]
        if (fam %in% names(acoustic_families)) family_df$Bioacoustic[i] <- acoustic_families[[fam]]
        if (fam %in% names(metabar_families))  family_df$Metabarcoding[i] <- metabar_families[[fam]]
      }

      write.csv(family_df,
                file.path(out_dir, "family_analysis.csv"), row.names = FALSE)
      cat("    Saved: family_analysis.csv\n")
    }
  }, error = function(e) {
    cat(sprintf("    [!] Family analysis failed: %s\n", e$message))
  })

  # ---- 5. Taxonomic diversity summary ----
  cat("  5. Taxonomic diversity...\n")

  if (!is.null(metabarcoding_data$presence_matrix)) {
    meta_sp <- setdiff(colnames(metabarcoding_data$presence_matrix), "site")
    diversity_summary <- data.frame(
      Method             = c("Bioacoustic", "Metabarcoding", "Combined"),
      Total_Species      = c(length(acoustic_species), length(meta_sp), total_species),
      Species_Shared     = c(length(both_methods), length(both_methods), length(both_methods)),
      Species_Exclusive  = c(length(only_acoustic), length(only_metabar), 0)
    )
    write.csv(diversity_summary,
              file.path(out_dir, "taxonomic_diversity_summary.csv"), row.names = FALSE)
    cat("    Saved: taxonomic_diversity_summary.csv\n")
  }

  # ---- 6. Save overall stats ----
  stats_df <- data.frame(
    Metric = c("Bioacoustic species", "Metabarcoding species", "Both methods",
               "Only bioacoustic", "Only metabarcoding", "Total species", "Overlap %"),
    Value  = c(length(acoustic_species), length(metabar_species), length(both_methods),
               length(only_acoustic), length(only_metabar), total_species,
               comparison_results$species_overlap$overlap_pct)
  )
  write.csv(stats_df,
            file.path(out_dir, "method_comparison_stats.csv"), row.names = FALSE)

  cat("  Method comparison complete\n")
  return(comparison_results)
}
