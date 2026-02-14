# baseline_vegetation_effects.R
# Vegetation effects on species richness:
#   - Spearman correlations (katydids, birds, metabarcoding vs vegetation)
#   - GLM Poisson for significant relationships
#   - GAM negative binomial (k=4) for non-linear detection
#
# Input: integrated_results/ CSV files
# Output: results/vegetation/ CSV and figures

library(tidyverse)
library(ggplot2)
library(mgcv)

if (!exists("extract_local_date")) {
  for (h in c("R/baseline_helpers.R", "baseline_helpers.R")) {
    if (file.exists(h)) { source(h); break }
  }
}

INPUT_DIR <- "integrated_results"
OUTPUT_DIR <- "results/vegetation"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)


#' Compute Spearman correlations between richness and vegetation variables
#'
#' @param richness_df Data frame with site and richness columns.
#' @param vegetation_df Data frame with site column (named Plot) and numeric
#'   vegetation variables.
#' @param site_col_veg Character. Site column name in vegetation. Default "Plot".
#' @return Data frame with variable, rho, p_value, n.
#' @export
correlate_richness_vegetation <- function(richness_df, vegetation_df,
                                          site_col_veg = "Plot") {
  veg <- standardise_vegetation_sites(vegetation_df)

  merged <- inner_join(richness_df, veg, by = "site")
  if (nrow(merged) < 4) {
    warning("Fewer than 4 common sites.")
    return(data.frame(variable = character(), rho = numeric(),
                      p_value = numeric(), n = integer()))
  }

  veg_vars <- setdiff(colnames(veg), "site")
  veg_vars <- veg_vars[sapply(merged[, veg_vars, drop = FALSE], is.numeric)]

  map_dfr(veg_vars, function(v) {
    vals <- merged[[v]]
    if (all(is.na(vals)) || sd(vals, na.rm = TRUE) == 0)
      return(data.frame(variable = v, rho = NA, p_value = NA,
                        n = sum(!is.na(vals))))
    ct <- cor.test(merged$richness, vals, method = "spearman", exact = FALSE)
    data.frame(variable = v, rho = ct$estimate, p_value = ct$p.value,
               n = sum(!is.na(vals)))
  }) %>% arrange(p_value)
}


#' Fit Poisson GLM for richness ~ vegetation predictor
#'
#' @param richness_df Data frame with site and richness.
#' @param vegetation_df Data frame.
#' @param predictor Character. Vegetation variable name.
#' @param site_col_veg Character.
#' @return List with model, summary, aic, and ggplot.
#' @export
fit_glm_poisson <- function(richness_df, vegetation_df, predictor,
                             site_col_veg = "Plot") {
  veg <- standardise_vegetation_sites(vegetation_df)
  merged <- inner_join(richness_df, veg, by = "site")

  model <- glm(as.formula(paste("richness ~", predictor)),
               data = merged, family = poisson)
  ms <- summary(model)

  pred_range <- range(merged[[predictor]], na.rm = TRUE)
  newdf <- data.frame(seq(pred_range[1], pred_range[2], length.out = 100))
  colnames(newdf) <- predictor
  newdf$predicted <- predict(model, newdata = newdf, type = "response")

  p <- ggplot(merged, aes_string(x = predictor, y = "richness")) +
    geom_point(size = 3, alpha = 0.7) +
    geom_line(data = newdf, aes_string(x = predictor, y = "predicted"),
              color = "firebrick", linewidth = 1) +
    labs(x = gsub("_", " ", predictor), y = "Species richness",
         title = paste("GLM Poisson: richness ~", predictor),
         subtitle = sprintf("p = %.4f, AIC = %.1f",
                            coef(ms)[2, 4], AIC(model))) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  list(model = model, summary = ms, aic = AIC(model), plot = p)
}


#' Fit GAM (negative binomial, k=4) for non-linear relationships
#'
#' @param richness_df Data frame with site and richness.
#' @param vegetation_df Data frame.
#' @param predictor Character.
#' @param k Integer. Basis dimension. Default 4.
#' @param site_col_veg Character.
#' @return List with model, edf, p_value, and ggplot, or NULL on failure.
#' @export
fit_gam_nb <- function(richness_df, vegetation_df, predictor,
                        k = 4, site_col_veg = "Plot") {
  veg <- standardise_vegetation_sites(vegetation_df)
  merged <- inner_join(richness_df, veg, by = "site")

  model <- tryCatch(
    gam(as.formula(paste("richness ~ s(", predictor, ", k =", k, ")")),
        data = merged, family = nb(), method = "REML"),
    error = function(e) { warning("GAM failed: ", e$message); NULL }
  )
  if (is.null(model)) return(NULL)

  ms <- summary(model)
  pred_range <- range(merged[[predictor]], na.rm = TRUE)
  newdf <- data.frame(seq(pred_range[1], pred_range[2], length.out = 100))
  colnames(newdf) <- predictor
  pv <- predict(model, newdata = newdf, type = "response", se.fit = TRUE)
  newdf$predicted <- pv$fit
  newdf$se <- pv$se.fit

  edf <- ms$s.table[1, "edf"]
  pval <- ms$s.table[1, "p-value"]

  p <- ggplot(merged, aes_string(x = predictor, y = "richness")) +
    geom_point(size = 3, alpha = 0.7) +
    geom_ribbon(data = newdf,
                aes_string(x = predictor,
                           ymin = "predicted - 1.96 * se",
                           ymax = "predicted + 1.96 * se"),
                alpha = 0.2, inherit.aes = FALSE) +
    geom_line(data = newdf,
              aes_string(x = predictor, y = "predicted"),
              color = "steelblue", linewidth = 1, inherit.aes = FALSE) +
    labs(x = gsub("_", " ", predictor), y = "Species richness",
         title = paste("GAM NB: richness ~ s(", predictor, ")"),
         subtitle = sprintf("edf = %.2f, p = %.4f", edf, pval)) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 12, face = "bold"))

  list(model = model, edf = edf, p_value = pval, plot = p)
}


# -- load data --
cat("Loading data...\n")
katydid_matrix <- read.csv(file.path(INPUT_DIR, "katydid_presence_matrix.csv"),
                           check.names = FALSE)
bird_matrix <- read.csv(file.path(INPUT_DIR, "bird_presence_matrix.csv"),
                        check.names = FALSE)
metabar_matrix <- tryCatch(
  read.csv(file.path(INPUT_DIR, "metabarcoding_presence_matrix.csv"),
           check.names = FALSE),
  error = function(e) NULL)
vegetation <- read.csv(file.path(INPUT_DIR, "vegetation_summary.csv"),
                       stringsAsFactors = FALSE)

katydid_richness <- richness_from_matrix(katydid_matrix)
bird_richness <- richness_from_matrix(bird_matrix)
metabar_richness <- if (!is.null(metabar_matrix)) richness_from_matrix(metabar_matrix) else NULL

# -- Spearman correlations --
cat("Computing Spearman correlations...\n")
katydid_cor <- correlate_richness_vegetation(katydid_richness, vegetation)
bird_cor <- correlate_richness_vegetation(bird_richness, vegetation)
metabar_cor <- if (!is.null(metabar_richness))
  correlate_richness_vegetation(metabar_richness, vegetation) else data.frame()

katydid_cor$taxon <- "katydid"
bird_cor$taxon <- "bird"
if (nrow(metabar_cor) > 0) metabar_cor$taxon <- "metabarcoding"
all_cor <- bind_rows(katydid_cor, bird_cor, metabar_cor)
write.csv(all_cor, file.path(OUTPUT_DIR, "vegetation_correlations.csv"),
          row.names = FALSE)

# -- GLM and GAM for significant correlations --
cat("Fitting GLM and GAM for significant relationships...\n")
sig_vars <- all_cor %>% filter(p_value < 0.05)
glm_results <- list()
gam_results <- list()

if (nrow(sig_vars) > 0) {
  for (i in seq_len(nrow(sig_vars))) {
    row <- sig_vars[i, ]
    rdf <- switch(row$taxon,
                  katydid = katydid_richness,
                  bird = bird_richness,
                  metabarcoding = metabar_richness)
    if (is.null(rdf)) next
    label <- paste0(row$taxon, "_", row$variable)

    glm_res <- tryCatch(fit_glm_poisson(rdf, vegetation, row$variable),
                         error = function(e) NULL)
    if (!is.null(glm_res)) {
      glm_results[[label]] <- glm_res
      ggsave(file.path(FIGURE_DIR, paste0("glm_", label, ".jpg")),
             glm_res$plot, width = 7, height = 6, dpi = 300)
    }

    gam_res <- tryCatch(fit_gam_nb(rdf, vegetation, row$variable),
                         error = function(e) NULL)
    if (!is.null(gam_res)) {
      gam_results[[label]] <- gam_res
      ggsave(file.path(FIGURE_DIR, paste0("gam_", label, ".jpg")),
             gam_res$plot, width = 7, height = 6, dpi = 300)
    }
  }
}

if (length(glm_results) > 0) {
  glm_df <- map_dfr(names(glm_results), function(l) {
    r <- glm_results[[l]]
    coefs <- coef(r$summary)
    data.frame(model = l, estimate = coefs[2, 1], se = coefs[2, 2],
               p = coefs[2, 4], aic = r$aic)
  })
  write.csv(glm_df, file.path(OUTPUT_DIR, "glm_results.csv"), row.names = FALSE)
}
if (length(gam_results) > 0) {
  gam_df <- map_dfr(names(gam_results), function(l) {
    r <- gam_results[[l]]
    data.frame(model = l, edf = r$edf, p_value = r$p_value)
  })
  write.csv(gam_df, file.path(OUTPUT_DIR, "gam_results.csv"), row.names = FALSE)
}

# -- summary --
cat("\n-- Vegetation effects summary --\n")
for (tx in c("katydid", "bird", "metabarcoding")) {
  sig <- all_cor %>% filter(taxon == tx, p_value < 0.05)
  cat(sprintf("  %s significant correlations: %d\n", tx, nrow(sig)))
  if (nrow(sig) > 0) {
    for (j in seq_len(nrow(sig)))
      cat(sprintf("    %s: rho=%.3f, p=%.4f\n",
                  sig$variable[j], sig$rho[j], sig$p_value[j]))
  }
}
cat(sprintf("  GLM models fitted: %d\n", length(glm_results)))
cat(sprintf("  GAM models fitted: %d\n", length(gam_results)))
cat("\nVegetation effects analysis complete.\n")
