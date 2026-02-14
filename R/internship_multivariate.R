# internship_multivariate.R
# Multivariate analyses for M1 IMABEE internship report
# NMDS with detection counts (Bray-Curtis), CCA, GAM, variance partitioning
#
# These analyses complement the PCoA (Jaccard presence/absence) already run
# in the main pipeline. Here we use detection COUNTS (abundance) instead.
#
# Reads from: integrated_results/ (pipeline outputs)
# Writes to:  results/internship_m1/multivariate/
# Sources:    baseline_helpers.R

run_internship_multivariate <- function(results, m1_dir) {

  cat("\n[INTERNSHIP-M1] Multivariate analyses (NMDS, CCA, GAM, variance partitioning)...\n")

  out_dir <- file.path(m1_dir, "multivariate")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  UTC_OFFSET_HOURS <- -5
  multivariate_results <- list()

  # ---- 0. Build detection count matrices from raw detections ----
  cat("  Building detection count matrices...\n")

  build_count_matrix <- function(raw_det, species_col = "common_name", site_col = "site") {
    if (is.null(raw_det) || nrow(raw_det) == 0) return(NULL)
    counts <- raw_det %>%
      dplyr::group_by(.data[[site_col]], .data[[species_col]]) %>%
      dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(species_col),
                         values_from = count, values_fill = 0) %>%
      dplyr::rename(site = dplyr::all_of(site_col))
    counts
  }

  katydid_det <- results$katydid_data$raw_detections
  if (is.null(katydid_det)) katydid_det <- results$katydid_data$detections
  bird_det <- results$bird_data$detections

  katydid_counts <- build_count_matrix(katydid_det)
  bird_counts    <- build_count_matrix(bird_det)

  if (is.null(katydid_counts) || nrow(katydid_counts) < 3) {
    cat("  [!] Insufficient katydid detection data for multivariate analyses\n")
    return(NULL)
  }

  # Vegetation
  veg <- results$vegetation_data$summary
  if ("Plot" %in% colnames(veg)) veg <- standardise_vegetation_sites(veg)

  # Richness values for point aesthetics
  katydid_rich <- richness_from_matrix(results$katydid_data$presence_matrix)
  bird_rich    <- richness_from_matrix(results$bird_data$presence_matrix)

  # ---- 1. NMDS with Bray-Curtis on detection counts (Figure 4 in report) ----
  cat("  1. NMDS ordination (Bray-Curtis, detection counts)...\n")

  tryCatch({
    # Prepare community matrix
    common_sites <- Reduce(intersect, list(
      katydid_counts$site,
      if (!is.null(bird_counts)) bird_counts$site else katydid_counts$site,
      if (!is.null(veg)) veg$site else katydid_counts$site
    ))

    if (length(common_sites) < 5) {
      cat("    [!] Too few common sites for NMDS\n")
    } else {
      kat_mat <- katydid_counts %>%
        dplyr::filter(site %in% common_sites) %>%
        dplyr::arrange(site)
      sites <- kat_mat$site
      kat_num <- kat_mat %>% dplyr::select(-site) %>% as.matrix()

      # Log(x+1) transform
      kat_log <- log1p(kat_num)

      # NMDS
      set.seed(42)
      nmds <- tryCatch({
        vegan::metaMDS(kat_log, distance = "bray", k = 2,
                       trymax = 100, trace = 0)
      }, error = function(e) {
        cat(sprintf("    NMDS failed, using PCoA fallback: %s\n", e$message))
        dist_mat <- vegan::vegdist(kat_log, method = "bray")
        pcoa <- stats::cmdscale(dist_mat, k = 2, eig = TRUE)
        list(points = pcoa$points, stress = NA, converged = FALSE)
      })

      scores_df <- data.frame(
        site  = sites,
        NMDS1 = if (!is.null(nmds$points)) nmds$points[, 1] else nmds$points[, 1],
        NMDS2 = if (!is.null(nmds$points)) nmds$points[, 2] else nmds$points[, 2]
      )

      # Add richness info
      scores_df <- scores_df %>%
        dplyr::left_join(katydid_rich %>% dplyr::rename(katydid_richness = richness), by = "site") %>%
        dplyr::left_join(bird_rich %>% dplyr::rename(bird_richness = richness), by = "site")

      # Add total detection counts
      kat_totals <- data.frame(site = sites, katydid_total_detection = rowSums(kat_num))
      scores_df <- dplyr::left_join(scores_df, kat_totals, by = "site")

      if (!is.null(bird_counts)) {
        bird_site_totals <- bird_counts %>%
          dplyr::filter(site %in% common_sites) %>%
          dplyr::arrange(site) %>%
          dplyr::mutate(bird_total_detection = rowSums(dplyr::select(., -site))) %>%
          dplyr::select(site, bird_total_detection)
        scores_df <- dplyr::left_join(scores_df, bird_site_totals, by = "site")
      }

      # Envfit with vegetation
      if (!is.null(veg)) {
        env_data <- veg %>%
          dplyr::filter(site %in% common_sites) %>%
          dplyr::arrange(site)
        env_vars <- env_data %>%
          dplyr::select(dplyr::where(is.numeric)) %>%
          dplyr::select(-dplyr::matches("^site$"))

        # Add total detections as environmental vectors
        if ("katydid_total_detection" %in% colnames(scores_df))
          env_vars$katydid_total_detection <- scores_df$katydid_total_detection
        if ("bird_total_detection" %in% colnames(scores_df))
          env_vars$bird_total_detection <- scores_df$bird_total_detection

        if (!is.null(nmds$points) && inherits(nmds, "metaMDS")) {
          ef <- vegan::envfit(nmds, env_vars, permutations = 999)
          capture.output(ef, file = file.path(out_dir, "nmds_envfit_results.txt"))
        }
      }

      # NMDS plot (Figure 4 style)
      stress_label <- if (!is.null(nmds$stress)) sprintf("Stress = %.3f", nmds$stress) else "PCoA fallback"

      p_nmds <- ggplot2::ggplot(scores_df, ggplot2::aes(x = NMDS1, y = NMDS2)) +
        ggplot2::geom_point(
          ggplot2::aes(size = katydid_total_detection,
                       color = if ("bird_total_detection" %in% colnames(scores_df))
                         bird_total_detection else bird_richness),
          alpha = 0.8
        ) +
        ggplot2::geom_text(ggplot2::aes(label = site), vjust = -0.8, size = 3) +
        ggplot2::scale_size_continuous(name = "Katydid\nDetections", range = c(2, 10)) +
        ggplot2::scale_color_viridis_c(name = "Bird\nDetections", option = "D") +
        ggplot2::labs(
          title    = "NMDS: Species Communities with Detection Abundances",
          subtitle = sprintf("%s | Size = Katydid detections, Color = Bird detections", stress_label)
        ) +
        ggplot2::theme_minimal()

      ggplot2::ggsave(file.path(out_dir, "nmds_detection_counts.jpg"),
                      p_nmds, width = 12, height = 9, dpi = 300)
      cat(sprintf("    NMDS complete (%s). Saved: nmds_detection_counts.jpg\n", stress_label))

      multivariate_results$nmds <- list(
        ordination = nmds,
        scores     = scores_df,
        plot       = p_nmds
      )
    }
  }, error = function(e) {
    cat(sprintf("    [!] NMDS failed: %s\n", e$message))
  })

  # ---- 2. CCA (Canonical Correspondence Analysis) ----
  cat("  2. CCA analysis...\n")

  tryCatch({
    if (!is.null(veg) && exists("kat_log") && exists("common_sites")) {
      env_cca <- veg %>%
        dplyr::filter(site %in% common_sites) %>%
        dplyr::arrange(site) %>%
        dplyr::select(dplyr::where(is.numeric))

      # Remove constant columns
      env_cca <- env_cca[, apply(env_cca, 2, var, na.rm = TRUE) > 0, drop = FALSE]

      if (ncol(env_cca) >= 1 && nrow(kat_log) >= 5) {
        cca_result <- vegan::cca(kat_log ~ ., data = env_cca)
        cca_anova  <- vegan::anova.cca(cca_result, permutations = 999)

        capture.output(summary(cca_result), file = file.path(out_dir, "cca_summary.txt"))
        capture.output(cca_anova, file = file.path(out_dir, "cca_anova.txt"))

        # CCA plot
        jpeg(file.path(out_dir, "cca_detections_environment.jpg"),
             width = 1000, height = 800, quality = 90)
        plot(cca_result, display = c("sites", "bp"),
             main = "CCA: Katydid Detection Counts ~ Environment")
        dev.off()

        cat("    CCA complete. Saved: cca_detections_environment.jpg\n")
        multivariate_results$cca <- cca_result
      }
    }
  }, error = function(e) {
    cat(sprintf("    [!] CCA failed: %s\n", e$message))
  })

  # ---- 3. GAM (non-linear relationships) ----
  cat("  3. GAM analysis (negative binomial)...\n")

  tryCatch({
    gam_dir <- file.path(out_dir, "gam")
    dir.create(gam_dir, showWarnings = FALSE)

    gam_results <- list()

    for (data_type in c("katydid", "bird")) {
      if (data_type == "katydid") {
        rich_data <- katydid_rich
      } else {
        rich_data <- bird_rich
      }

      merged <- dplyr::inner_join(rich_data, veg, by = "site")

      for (vv in veg_vars_available(veg)) {
        if (!vv %in% colnames(merged)) next
        valid <- !is.na(merged$richness) & !is.na(merged[[vv]])
        if (sum(valid) < 8) next

        gam_model <- tryCatch({
          mgcv::gam(richness ~ s(get(vv), k = 4),
                     data = merged[valid, ],
                     family = mgcv::nb())
        }, error = function(e) NULL)

        if (!is.null(gam_model)) {
          s_gam <- summary(gam_model)
          p_smooth <- s_gam$s.table[1, "p-value"]

          gam_results[[paste0(data_type, "_", vv)]] <- list(
            model   = gam_model,
            p_value = p_smooth,
            edf     = s_gam$s.table[1, "edf"],
            r2      = s_gam$r.sq
          )

          if (p_smooth < 0.1) {
            gam_plot <- ggplot2::ggplot(merged[valid, ],
                                         ggplot2::aes(x = .data[[vv]], y = richness)) +
              ggplot2::geom_point(size = 3, alpha = 0.7) +
              ggplot2::geom_smooth(method = "gam",
                                   method.args = list(family = mgcv::nb(), k = 4),
                                   se = TRUE, color = "#E41A1C") +
              ggplot2::labs(
                title    = sprintf("GAM: %s richness ~ %s", data_type, vv),
                subtitle = sprintf("p = %.4f, edf = %.2f, R2 = %.3f",
                                   p_smooth, s_gam$s.table[1, "edf"], s_gam$r.sq),
                x = gsub("_", " ", vv),
                y = sprintf("%s species richness", data_type)
              ) +
              ggplot2::theme_minimal()

            ggplot2::ggsave(
              file.path(gam_dir, sprintf("gam_%s_%s.jpg", data_type, vv)),
              gam_plot, width = 10, height = 7, dpi = 300
            )
          }
        }
      }
    }

    # Save GAM summary
    if (length(gam_results) > 0) {
      gam_summary <- do.call(rbind, lapply(names(gam_results), function(nm) {
        g <- gam_results[[nm]]
        data.frame(Model = nm, P_value = round(g$p_value, 4),
                   EDF = round(g$edf, 2), R_squared = round(g$r2, 3),
                   stringsAsFactors = FALSE)
      }))
      write.csv(gam_summary, file.path(gam_dir, "gam_summary.csv"), row.names = FALSE)
      cat(sprintf("    GAM: %d models fitted, %d with p < 0.1\n",
                  nrow(gam_summary), sum(gam_summary$P_value < 0.1)))
    }

    multivariate_results$gam <- gam_results
  }, error = function(e) {
    cat(sprintf("    [!] GAM analysis failed: %s\n", e$message))
  })

  # ---- 4. Variance partitioning ----
  cat("  4. Variance partitioning...\n")

  tryCatch({
    if (exists("kat_log") && !is.null(veg) && exists("common_sites")) {
      env_vp <- veg %>%
        dplyr::filter(site %in% common_sites) %>%
        dplyr::arrange(site) %>%
        dplyr::select(dplyr::where(is.numeric))

      # Split into tree and liana variable groups
      tree_vars  <- env_vp %>% dplyr::select(dplyr::starts_with("tree"))
      liana_vars <- env_vp %>% dplyr::select(dplyr::starts_with("liana"))

      if (ncol(tree_vars) >= 1 && ncol(liana_vars) >= 1) {
        vp <- vegan::varpart(kat_log, tree_vars, liana_vars)

        capture.output(vp, file = file.path(out_dir, "variance_partitioning.txt"))

        jpeg(file.path(out_dir, "variance_partitioning.jpg"),
             width = 800, height = 600, quality = 90)
        plot(vp, Xnames = c("Trees", "Lianas"),
             main = "Variance Partitioning: Katydid Community")
        dev.off()

        cat("    Variance partitioning complete. Saved: variance_partitioning.jpg\n")
        multivariate_results$variance_partitioning <- vp
      }
    }
  }, error = function(e) {
    cat(sprintf("    [!] Variance partitioning failed: %s\n", e$message))
  })

  # Save detection count matrices for reference
  if (!is.null(katydid_counts)) {
    write.csv(katydid_counts,
              file.path(out_dir, "katydid_detection_count_matrix.csv"), row.names = FALSE)
  }
  if (!is.null(bird_counts)) {
    write.csv(bird_counts,
              file.path(out_dir, "bird_detection_count_matrix.csv"), row.names = FALSE)
  }

  cat("  Multivariate analyses complete\n")
  return(multivariate_results)
}


# Helper: get available vegetation numeric variable names
veg_vars_available <- function(veg) {
  num_cols <- colnames(veg)[sapply(veg, is.numeric)]
  num_cols[!num_cols %in% c("site")]
}
