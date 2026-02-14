# internship_report_annexes.R
# Report annexes generation for M1 IMABEE internship report
# Produces Appendix 1-5 tables and figure index
#
# Reads from: integrated_results/ + results/internship_m1/
# Writes to:  results/internship_m1/REPORT_ANNEXES/
# Sources:    baseline_helpers.R

run_internship_report_annexes <- function(results, m1_dir,
                                          glm_results = NULL,
                                          temporal_results = NULL,
                                          comparison_results = NULL) {

  cat("\n[INTERNSHIP-M1] Generating report annexes...\n")

  annexes_dir <- file.path(m1_dir, "REPORT_ANNEXES")
  tables_dir  <- file.path(annexes_dir, "Tables")
  figures_dir <- file.path(annexes_dir, "Figures")
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

  # ===========================================================================
  # APPENDIX 1 - Katydid species list with detection stats
  # ===========================================================================
  cat("  Appendix 1: Katydid species list...\n")

  katydid_det <- results$katydid_data$raw_detections
  if (is.null(katydid_det)) katydid_det <- results$katydid_data$detections

  if (!is.null(katydid_det) && nrow(katydid_det) > 0) {
    species_list <- katydid_det %>%
      dplyr::group_by(common_name) %>%
      dplyr::summarise(
        Total_Detections   = dplyr::n(),
        Sites_Detected     = dplyr::n_distinct(site),
        Detection_Rate     = round(dplyr::n() / dplyr::n_distinct(site), 1),
        Mean_Confidence    = round(mean(confidence, na.rm = TRUE), 3),
        Min_Confidence     = round(min(confidence, na.rm = TRUE), 3),
        First_Detection    = min(as.Date(substr(
          stringr::str_extract(file_path, "20[0-9]{6}"), 1, 8), format = "%Y%m%d"),
          na.rm = TRUE),
        Last_Detection     = max(as.Date(substr(
          stringr::str_extract(file_path, "20[0-9]{6}"), 1, 8), format = "%Y%m%d"),
          na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(Total_Detections)) %>%
      dplyr::mutate(Species_Rank = dplyr::row_number())

    # Add detection span
    species_list$Detection_Span_Days <- as.numeric(
      species_list$Last_Detection - species_list$First_Detection
    )

    write.csv(species_list,
              file.path(tables_dir, "Appendix1_Katydid_Species_List.csv"),
              row.names = FALSE)
    cat(sprintf("    %d species listed\n", nrow(species_list)))
  }

  # ===========================================================================
  # APPENDIX 2 - Species detected by both methods
  # ===========================================================================
  cat("  Appendix 2: Species detected by both methods...\n")

  if (!is.null(comparison_results) && !is.null(comparison_results$species_overlap)) {
    both <- comparison_results$species_overlap$both_methods

    if (length(both) > 0) {
      common_table <- data.frame(Species = both, Detection_Method = "Both",
                                  Sites_Bioacoustic = NA, stringsAsFactors = FALSE)

      for (i in seq_along(both)) {
        sp <- both[i]
        if (sp %in% colnames(results$katydid_data$presence_matrix)) {
          common_table$Sites_Bioacoustic[i] <- sum(
            results$katydid_data$presence_matrix[[sp]], na.rm = TRUE
          )
        }
      }

      write.csv(common_table,
                file.path(tables_dir, "Appendix2_Species_Both_Methods.csv"),
                row.names = FALSE)
      cat(sprintf("    %d species in common\n", nrow(common_table)))
    }
  }

  # ===========================================================================
  # APPENDIX 3 - GLM summary (significant only)
  # ===========================================================================
  cat("  Appendix 3: GLM summary...\n")

  if (!is.null(glm_results) && length(glm_results) > 0) {
    glm_table <- do.call(rbind, lapply(names(glm_results), function(nm) {
      r <- glm_results[[nm]]
      data.frame(
        Relationship    = nm,
        Response        = r$response,
        Predictor       = r$predictor,
        P_value         = round(r$p_value, 4),
        Coefficient     = round(r$coefficient, 4),
        Direction       = r$direction,
        AIC             = round(r$aic, 2),
        Sample_Size     = r$sample_size,
        Pseudo_R_squared = round(r$pseudo_r2, 3),
        stringsAsFactors = FALSE
      )
    })) %>% dplyr::arrange(P_value)

    write.csv(glm_table,
              file.path(tables_dir, "Appendix3_GLM_Summary.csv"), row.names = FALSE)
    cat(sprintf("    %d significant GLM relationships\n", nrow(glm_table)))
  }

  # ===========================================================================
  # APPENDIX 4 - Vegetation correlations
  # ===========================================================================
  cat("  Appendix 4: Vegetation correlations...\n")

  corr_summary <- data.frame()
  veg <- results$vegetation_data$summary
  if ("Plot" %in% colnames(veg)) veg <- standardise_vegetation_sites(veg)

  veg_num_vars <- colnames(veg)[sapply(veg, is.numeric)]
  veg_num_vars <- setdiff(veg_num_vars, "site")

  for (group_info in list(
    list(name = "Katydids", method = "Bioacoustic", data = results$katydid_data),
    list(name = "Birds",    method = "Bioacoustic", data = results$bird_data),
    list(name = "Orthoptera", method = "Metabarcoding", data = results$metabarcoding_data)
  )) {
    rich <- richness_from_matrix(group_info$data$presence_matrix)
    merged <- dplyr::inner_join(rich, veg, by = "site")

    for (vv in veg_num_vars) {
      if (!vv %in% colnames(merged)) next
      valid <- !is.na(merged$richness) & !is.na(merged[[vv]])
      if (sum(valid) < 5) next

      ct <- cor.test(merged$richness[valid], merged[[vv]][valid], method = "spearman")

      corr_summary <- rbind(corr_summary, data.frame(
        Group       = group_info$name,
        Method      = group_info$method,
        Variable    = vv,
        Correlation = round(ct$estimate, 3),
        P_value     = round(ct$p.value, 4),
        Significance = dplyr::case_when(
          ct$p.value < 0.001 ~ "***",
          ct$p.value < 0.01  ~ "**",
          ct$p.value < 0.05  ~ "*",
          ct$p.value < 0.1   ~ ".",
          TRUE               ~ "ns"
        ),
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(corr_summary) > 0) {
    write.csv(corr_summary,
              file.path(tables_dir, "Appendix4_Vegetation_Correlations.csv"),
              row.names = FALSE)
    cat(sprintf("    %d correlations computed\n", nrow(corr_summary)))
  }

  # ===========================================================================
  # APPENDIX 5 - Katydid and bird richness per site
  # ===========================================================================
  cat("  Appendix 5: Richness per site...\n")

  katydid_rich <- richness_from_matrix(results$katydid_data$presence_matrix)
  bird_rich    <- richness_from_matrix(results$bird_data$presence_matrix)

  richness_table <- dplyr::inner_join(
    katydid_rich %>% dplyr::rename(Katydid_Richness = richness),
    bird_rich    %>% dplyr::rename(Bird_Richness    = richness),
    by = "site"
  ) %>%
    dplyr::mutate(
      Katydid_Bird_Ratio = round(Katydid_Richness / Bird_Richness, 3),
      Combined_Richness  = Katydid_Richness + Bird_Richness
    ) %>%
    dplyr::arrange(dplyr::desc(Combined_Richness))

  # Convert site names to YB-P format for the report
  richness_table$Site <- gsub("^S(\\d+)$", "YB-P\\1", richness_table$site)

  write.csv(richness_table %>% dplyr::select(Site, Katydid_Richness, Bird_Richness,
                                              Katydid_Bird_Ratio, Combined_Richness),
            file.path(tables_dir, "Appendix5_Richness_Per_Site.csv"),
            row.names = FALSE)
  cat(sprintf("    %d sites\n", nrow(richness_table)))

  # ===========================================================================
  # Copy important figures
  # ===========================================================================
  cat("  Copying important figures...\n")

  fig_patterns <- c(
    "species_detection_comparison", "katydids_venn_diagram",
    "hourly_activity_pattern", "Bird_vs_Katydid",
    "nmds_detection_counts", "katydid_community_ordination"
  )

  # Search in m1_dir and integrated_results
  search_dirs <- c(m1_dir, results$output_dir %||% "integrated_results")
  fig_counter <- 1

  for (pat in fig_patterns) {
    for (sdir in search_dirs) {
      if (!dir.exists(sdir)) next
      files <- list.files(sdir, pattern = pat, recursive = TRUE, full.names = TRUE)
      files <- files[!grepl("REPORT_ANNEXES", files)]
      if (length(files) > 0) {
        ext <- tools::file_ext(files[1])
        dest <- file.path(figures_dir, sprintf("FigureA%d_%s.%s", fig_counter, pat, ext))
        file.copy(files[1], dest, overwrite = TRUE)
        fig_counter <- fig_counter + 1
        break
      }
    }
  }
  cat(sprintf("    %d figures copied\n", fig_counter - 1))

  # ===========================================================================
  # Annexes index
  # ===========================================================================
  index <- data.frame(
    Appendix = paste0("Appendix ", 1:5),
    Description = c(
      "List of katydid species detected by bioacoustics",
      "Species detected by both methods",
      "Complete GLM analysis summary (significant)",
      "Vegetation correlation summary",
      "Total katydid and bird richness per site"
    ),
    stringsAsFactors = FALSE
  )
  write.csv(index, file.path(annexes_dir, "ANNEXES_INDEX.csv"), row.names = FALSE)

  tables_n  <- length(list.files(tables_dir, pattern = "\\.csv$"))
  figures_n <- length(list.files(figures_dir, pattern = "\\.(jpg|png|pdf)$"))
  cat(sprintf("  Report annexes complete: %d tables, %d figures\n", tables_n, figures_n))

  return(list(annexes_dir = annexes_dir, tables = tables_n, figures = figures_n))
}
