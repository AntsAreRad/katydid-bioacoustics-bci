# internship_glm_analysis.R
# Comprehensive GLM analysis for M1 IMABEE internship report
# All combinations: bird~katydid, bird~veg, katydid~veg richness
# Poisson family, diagnostic plots, pseudo-R2
#
# Reads from: integrated_results/ (pipeline outputs)
# Writes to:  results/internship_m1/glm/
# Sources:    baseline_helpers.R

run_internship_glm <- function(results, m1_dir) {

  cat("\n[INTERNSHIP-M1] GLM analysis...\n")

  out_dir <- file.path(m1_dir, "glm")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Prepare richness data ----
  katydid_rich <- richness_from_matrix(results$katydid_data$presence_matrix)
  bird_rich    <- richness_from_matrix(results$bird_data$presence_matrix)

  # Vegetation data
  veg <- results$vegetation_data$summary
  if ("Plot" %in% colnames(veg)) {
    veg <- standardise_vegetation_sites(veg)
  }

  # Merge all into one data frame
  site_data <- dplyr::inner_join(
    katydid_rich %>% dplyr::rename(katydid_richness = richness),
    bird_rich    %>% dplyr::rename(bird_richness    = richness),
    by = "site"
  )

  if (!is.null(veg) && "site" %in% colnames(veg)) {
    site_data <- dplyr::left_join(site_data, veg, by = "site")
  }

  cat(sprintf("  %d sites with complete data\n", nrow(site_data)))

  # ---- Define all GLM combinations ----
  veg_vars <- c("tree_species_richness", "tree_abundance", "tree_total_basal_area",
                "liana_species_richness", "liana_rooted_stems", "liana_total_basal_area")
  # Keep only variables actually present in data
 veg_vars <- veg_vars[veg_vars %in% colnames(site_data)]

  combinations <- list()

  # Bird ~ Katydid
  combinations[["Bird_vs_Katydid"]] <- list(
    response = "bird_richness", predictor = "katydid_richness"
  )

  # Bird ~ vegetation variables
  for (v in veg_vars) {
    name <- paste0("Bird_vs_", v)
    combinations[[name]] <- list(response = "bird_richness", predictor = v)
  }

  # Katydid ~ vegetation variables
  for (v in veg_vars) {
    name <- paste0("Katydid_vs_", v)
    combinations[[name]] <- list(response = "katydid_richness", predictor = v)
  }

  # ---- Run all GLMs ----
  glm_results <- list()
  summary_rows <- list()

  for (name in names(combinations)) {
    combo <- combinations[[name]]
    resp <- combo$response
    pred <- combo$predictor

    if (!all(c(resp, pred) %in% colnames(site_data))) next

    resp_data <- site_data[[resp]]
    pred_data <- site_data[[pred]]

    # Remove NAs
    valid <- !is.na(resp_data) & !is.na(pred_data)
    if (sum(valid) < 5) next

    tryCatch({
      model <- glm(resp_data[valid] ~ pred_data[valid], family = poisson)
      s <- summary(model)

      p_val   <- coef(s)[2, 4]
      coeff   <- coef(s)[2, 1]
      aic_val <- AIC(model)
      pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
      direction <- ifelse(coeff > 0, "positive", "negative")

      glm_results[[name]] <- list(
        model      = model,
        response   = resp,
        predictor  = pred,
        p_value    = p_val,
        coefficient = coeff,
        direction  = direction,
        aic        = aic_val,
        pseudo_r2  = pseudo_r2,
        sample_size = sum(valid)
      )

      summary_rows[[name]] <- data.frame(
        Relationship    = name,
        Response        = resp,
        Predictor       = pred,
        P_value         = round(p_val, 4),
        Coefficient     = round(coeff, 4),
        Direction       = direction,
        AIC             = round(aic_val, 2),
        Sample_Size     = sum(valid),
        Pseudo_R_squared = round(pseudo_r2, 3),
        Significant     = ifelse(p_val < 0.05, "yes", "no"),
        stringsAsFactors = FALSE
      )

      # ---- Plot significant GLMs (Figure 3 in report) ----
      if (p_val < 0.05) {
        plot_df <- data.frame(x = pred_data[valid], y = resp_data[valid])

        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(size = 3, alpha = 0.7, color = "#0072B2") +
          ggplot2::geom_smooth(method = "glm", method.args = list(family = "poisson"),
                               se = TRUE, color = "#D55E00", fill = "#D55E00",
                               alpha = 0.2) +
          ggplot2::labs(
            title    = sprintf("Significant GLM: %s vs %s", resp, pred),
            subtitle = sprintf("p = %.3f (%s relationship), Family = poisson, AIC = %.2f",
                               p_val, direction, aic_val),
            x = gsub("_", " ", pred),
            y = gsub("_", " ", resp)
          ) +
          ggplot2::theme_minimal()

        plot_name <- paste0(gsub("[^A-Za-z0-9]", "_", name), ".jpg")
        ggplot2::ggsave(file.path(out_dir, plot_name), p, width = 10, height = 7, dpi = 300)
        cat(sprintf("    Significant: %s (p=%.4f, pseudo-R2=%.3f) -> %s\n",
                    name, p_val, pseudo_r2, plot_name))

        # Diagnostic plot
        diag_name <- paste0("diagnostic_", gsub("[^A-Za-z0-9]", "_", name), ".jpg")
        jpeg(file.path(out_dir, diag_name), width = 1000, height = 800, quality = 90)
        par(mfrow = c(2, 2))
        plot(model)
        dev.off()
      }

    }, error = function(e) {
      cat(sprintf("    [!] GLM failed for %s: %s\n", name, e$message))
    })
  }

  # ---- Save summary table ----
  if (length(summary_rows) > 0) {
    glm_summary <- dplyr::bind_rows(summary_rows) %>%
      dplyr::arrange(P_value)

    write.csv(glm_summary,
              file.path(out_dir, "GLM_complete_summary.csv"), row.names = FALSE)

    sig_summary <- glm_summary %>% dplyr::filter(Significant == "yes")
    if (nrow(sig_summary) > 0) {
      write.csv(sig_summary,
                file.path(out_dir, "GLM_significant_relationships_summary.csv"),
                row.names = FALSE)
    }

    cat(sprintf("  GLM analysis complete: %d tested, %d significant\n",
                nrow(glm_summary), nrow(sig_summary)))
  }

  # Return only significant results (compatible with create_report_annexes)
  sig_results <- glm_results[sapply(glm_results, function(x) x$p_value < 0.05)]
  return(sig_results)
}
