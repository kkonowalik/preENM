# Recording effort from GBIF ----------------------------------------------------

#' SQL query for a GBIF recording-effort layer
#'
#' Builds the GBIF SQL download that counts, per grid cell of the study
#' extent: all records of a reference group (`n_records`, default vascular
#' plants), the number of recorded species (`n_species`, less inflated by a
#' few huge atlas datasets: the recommended effort proxy) and the records of
#' a target group of companion species (`n_target`, for a target-group
#' background). GBIF aggregates on its side, so the result is one small
#' table, not millions of records. Run it with [download_gbif_effort()] and
#' turn the table into rasters with [gbif_effort_rasters()]; use a layer as
#' `effort` in [sample_background()].
#'
#' @param extent Study extent `c(xmin, xmax, ymin, ymax)` in degrees.
#' @param cell Cell size in degrees (default 0.1).
#' @param taxon_filter SQL condition defining the reference group (default
#'   vascular plants: `"phylum = 'Tracheophyta'"`).
#' @param target Optional species names of the target group (e.g. companions
#'   with the same habitat and recorders, without your species itself).
#' @param max_uncertainty_m Records with a larger stated uncertainty are
#'   left out.
#' @param basis Accepted `basisOfRecord` values.
#' @return The SQL query (a single string).
#' @export
#' @examples
#' cat(gbif_effort_sql(c(-32, 75, 34, 82),
#'                     target = c("Silene acaulis", "Salix herbacea")))
gbif_effort_sql <- function(extent, cell = 0.1,
                            taxon_filter = "phylum = 'Tracheophyta'",
                            target = NULL, max_uncertainty_m = 10000,
                            basis = c("HUMAN_OBSERVATION",
                                      "PRESERVED_SPECIMEN", "OBSERVATION",
                                      "OCCURRENCE", "MATERIAL_SAMPLE")) {
  e <- unname(extent)
  if (length(e) != 4 || e[1] >= e[2] || e[3] >= e[4]) {
    stop("`extent` must be c(xmin, xmax, ymin, ymax).", call. = FALSE)
  }
  k <- 1 / cell
  if (abs(k - round(k)) > 1e-9) stop("`cell` must divide 1 degree, e.g. 0.1, ",
                                     "0.25 or 0.5.", call. = FALSE)
  k <- round(k)
  quote_sql <- function(x) paste0("'", gsub("'", "''", x), "'",
                                  collapse = ", ")
  tgt <- if (length(target)) {
    sprintf("SUM(CASE WHEN species IN (%s) THEN 1 ELSE 0 END) AS n_target",
            quote_sql(target))
  } else "0 AS n_target"
  paste(
    sprintf(paste("SELECT FLOOR(decimalLatitude * %d) AS lat_bin,",
                  "FLOOR(decimalLongitude * %d) AS lon_bin,"), k, k),
    "COUNT(*) AS n_records, COUNT(DISTINCT speciesKey) AS n_species,", tgt,
    "FROM occurrence",
    "WHERE", taxon_filter,
    "AND hasCoordinate = TRUE AND hasGeospatialIssues = FALSE",
    "AND occurrenceStatus = 'PRESENT'",
    sprintf("AND basisOfRecord IN (%s)", quote_sql(basis)),
    sprintf(paste("AND (coordinateUncertaintyInMeters IS NULL OR",
                  "coordinateUncertaintyInMeters <= %d)"),
            as.integer(max_uncertainty_m)),
    sprintf(paste("AND decimalLatitude >= %s AND decimalLatitude < %s",
                  "AND decimalLongitude >= %s AND decimalLongitude < %s"),
            num_arg(e[3]), num_arg(e[4]), num_arg(e[1]), num_arg(e[2])),
    sprintf(paste("GROUP BY FLOOR(decimalLatitude * %d),",
                  "FLOOR(decimalLongitude * %d)"), k, k))
}

#' Download a GBIF recording-effort table
#'
#' Runs the SQL query of [gbif_effort_sql()] as a GBIF download (package
#' `rgbif` >= 3.8). **A GBIF account is needed**, and its credentials must be
#' stored as environment variables, never in your scripts: run
#' `file.edit("~/.Renviron")` (or `usethis::edit_r_environ()`), add the three
#' lines
#' ```
#' GBIF_USER=your_user_name
#' GBIF_PWD=your_password
#' GBIF_EMAIL=your_email
#' ```
#' and restart R. SQL downloads are marked experimental by GBIF; if your
#' account is refused, download the target-group species normally and count
#' records per cell instead. Large extents can take an hour or more. **Cite
#' the DOI** written to `effort_download_DOI.txt`.
#'
#' @param sql Query from [gbif_effort_sql()].
#' @param path Folder for the download and the table (re-used: an existing
#'   `effort_table.csv` is returned without downloading again).
#' @param wait Wait for GBIF to finish (status checked every minute).
#' @return The path of `effort_table.csv` (or, with `wait = FALSE`, the
#'   download key).
#' @export
#' @examples
#' \dontrun{
#' q <- gbif_effort_sql(c(-32, 75, 34, 82), target = c("Silene acaulis"))
#' tab <- download_gbif_effort(q, "gbif_effort")
#' eff <- gbif_effort_rasters(tab, c(-32, 75, 34, 82))
#' }
download_gbif_effort <- function(sql, path, wait = TRUE) {
  f_tab <- file.path(path, "effort_table.csv")
  if (file.exists(f_tab)) return(f_tab)
  if (!requireNamespace("rgbif", quietly = TRUE) ||
      utils::packageVersion("rgbif") < "3.8") {
    stop("Install rgbif >= 3.8: install.packages(\"rgbif\")", call. = FALSE)
  }
  miss <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")[
    !nzchar(Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")))]
  if (length(miss)) {
    stop("GBIF credentials missing: ", paste(miss, collapse = ", "),
         ". Add them to ~/.Renviron (see ?download_gbif_effort) and restart R.",
         call. = FALSE)
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  rgbif::occ_download_sql_prep(sql)        # validation: stops if GBIF rejects it
  dl <- rgbif::occ_download_sql(sql)
  message("GBIF download key: ", as.character(dl))
  if (!wait) return(as.character(dl))
  rgbif::occ_download_wait(dl, status_ping = 60)
  z <- rgbif::occ_download_get(dl, path = path, overwrite = TRUE)
  tab <- as.data.frame(rgbif::occ_download_import(z))
  names(tab) <- tolower(names(tab))
  utils::write.csv(tab, f_tab, row.names = FALSE)
  doi <- tryCatch(rgbif::occ_download_meta(dl)$doi, error = function(e) NA)
  writeLines(paste0("https://doi.org/", doi),
             file.path(path, "effort_download_DOI.txt"))
  message("Cite this download: https://doi.org/", doi)
  f_tab
}

#' Recording-effort rasters from a GBIF effort table
#'
#' @param table The effort table (data.frame or path to
#'   `effort_table.csv`) with `lat_bin`, `lon_bin` and count columns.
#' @param extent,cell The extent and cell size used in [gbif_effort_sql()].
#' @param path Optional folder to write `effort_<layer>.tif`.
#' @return A [terra::SpatRaster] with layers `n_records`, `n_species`,
#'   `n_target` (those present); cells without records are 0.
#' @export
#' @examples
#' tab <- data.frame(lat_bin = c(450, 451), lon_bin = c(100, 101),
#'                   n_records = c(10, 3), n_species = c(8, 2),
#'                   n_target = c(1, 0))
#' eff <- gbif_effort_rasters(tab, c(9.5, 10.5, 44.5, 45.5))
#' eff
gbif_effort_rasters <- function(table, extent, cell = 0.1, path = NULL) {
  if (is.character(table)) table <- utils::read.csv(table)
  names(table) <- tolower(names(table))
  if (!all(c("lat_bin", "lon_bin") %in% names(table))) {
    stop("The table needs lat_bin and lon_bin columns.", call. = FALSE)
  }
  e <- unname(extent)
  k <- round(1 / cell)
  r0 <- terra::rast(xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4],
                    resolution = cell, crs = "EPSG:4326", vals = 0)
  cl <- terra::cellFromXY(r0, cbind((table$lon_bin + 0.5) / k,
                                    (table$lat_bin + 0.5) / k))
  ok <- !is.na(cl)
  lays <- intersect(c("n_records", "n_species", "n_target"), names(table))
  out <- lapply(lays, function(v) {
    r <- r0
    vv <- terra::values(r, mat = FALSE)
    vv[cl[ok]] <- table[[v]][ok]
    terra::values(r) <- vv
    names(r) <- v
    if (!is.null(path)) {
      if (!dir.exists(path)) dir.create(path, recursive = TRUE)
      terra::writeRaster(r, file.path(path, paste0("effort_", v, ".tif")),
                         overwrite = TRUE)
    }
    r
  })
  terra::rast(out)
}
