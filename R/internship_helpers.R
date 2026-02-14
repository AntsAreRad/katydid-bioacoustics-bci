# internship_helpers.R
# Shared utility functions for M1 IMABEE internship analyses.
# Sourced by all internship_*.R scripts.
# These functions were originally written during the M1 internship.

#' Extract local date from AudioMoth file paths
#'
#' Searches for YYYYMMDD_HHMMSS pattern in file paths where the year starts
#' with 20 (2000-2099). This avoids false matches on numeric deployment IDs
#' like "88260" in the path.
#'
#' Dates before 2024 are set to NA (known AudioMoth clock issue at S21
#' where the RTC was not synchronised, producing year-2000 timestamps).
#'
#' @param file_paths Character vector.
#' @param utc_offset Numeric. Hours offset from UTC. Default -5 (Panama).
#' @return Date vector (local dates).
#' @export
extract_local_date <- function(file_paths, utc_offset = -5) {
  pattern <- "(20[0-9]{2}[01][0-9][0-3][0-9]_[0-2][0-9][0-5][0-9][0-5][0-9])"
  matches <- stringr::str_extract(file_paths, pattern)
  utc <- lubridate::parse_date_time(matches, "Ymd_HMS", tz = "UTC", quiet = TRUE)
  local <- as.Date(utc + lubridate::hours(utc_offset))
  # Filter out dates before 2024 (AudioMoth clock not synchronised)
  local[!is.na(local) & lubridate::year(local) < 2024] <- NA
  local
}


#' Calculate species richness per site from presence/absence matrix
#'
#' @param presence_matrix Data frame with site as first column, species as
#'   remaining columns (0/1 values).
#' @return Data frame with site and richness.
#' @export
richness_from_matrix <- function(presence_matrix) {
  sp <- setdiff(colnames(presence_matrix), "site")
  data.frame(site = presence_matrix$site,
             richness = rowSums(presence_matrix[, sp, drop = FALSE]),
             stringsAsFactors = FALSE)
}


#' Convert vegetation site names to match bioacoustic site names
#'
#' The vegetation data uses "YB-P01" format, while bioacoustic data uses "S01".
#' This converts "YB-P01" -> "S01", "YB-P02" -> "S02", etc.
#'
#' @param veg_df Data frame with a Plot or site column.
#' @return Data frame with standardised site column.
#' @export
standardise_vegetation_sites <- function(veg_df) {
  if ("Plot" %in% colnames(veg_df)) {
    veg_df$site <- gsub("YB-P", "S", veg_df$Plot)
  } else if ("site" %in% colnames(veg_df)) {
    veg_df$site <- gsub("YB-P", "S", veg_df$site)
  }
  veg_df
}
