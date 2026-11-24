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

# Biogeographical ignorance ---------------------------------------------------------

#' Map of biogeographical ignorance from recording effort
#'
#' Where few records of any species have been made, the absence of records
#' of your species says little, and a model's prediction there rests on
#' extrapolation from elsewhere. A map of biogeographical ignorance shows
#' this uncertainty of the data next to the prediction (Tessarolo et al.
#' 2021), so that low-confidence areas of a projection can be flagged.
#'
#' **Ignorance score.** The effort of each cell (records of a reference
#' group, e.g. from [gbif_effort_rasters()]) is turned into a score from 0
#' (well recorded) to 1 (nothing recorded), with one of the algorithms of
#' Ruete (2015), on which maps of biogeographical ignorance build:
#'
#' * `"half"` (default): `I = O / (N + O)`, where `N` is the effort and `O`
#'   (`half`) the effort that halves ignorance; it does not depend on the
#'   most recorded cell. With `O = 1` one record already halves it; for a
#'   whole plant or bird group, a larger `O` (e.g. 10 or 50) is more
#'   realistic.
#' * `"log"`: `I = 1 - log(N + 1) / log(Nmax + 1)`; relative to the most
#'   recorded cell, separating "few" from "enough" records.
#' * `"linear"`: `I = 1 - N / Nmax`; only for effort without a long tail.
#'
#' With `index = TRUE` the effort is the species observation index
#' `n_records / n_species` (records per species; Ruete 2015), less
#' inflated where a few species are recorded very often. `radius_km`
#' spreads the effort of each cell to its neighbours within that distance
#' before scoring (information decays with distance; Tessarolo et al.
#' 2021 used the range of spatial autocorrelation of the predictors, which
#' [suggest_block_size()] estimates).
#'
#' **Limits.** The full maps of Tessarolo et al. (2021) also weigh each
#' record by its taxonomic reliability, its age and the completeness of the
#' inventory of the cell, and interpolate between cells in geographic and
#' environmental space; that needs record-level data. This function uses
#' the effort counts only.
#'
#' **With a suitability raster** (`suitability`), the ignorance is brought
#' to its grid and the result adds `suitability_known` (the suitability
#' where ignorance is at most `threshold`, `NA` elsewhere) and `class`: the
#' four combinations of suitable / unsuitable (at `suitability_threshold`)
#' and known / ignorant (at `threshold`), as in the joint maps of
#' prediction and ignorance of Tessarolo et al. (2021).
#'
#' @param effort Effort raster(s), e.g. from [gbif_effort_rasters()], or a
#'   path.
#' @param layer The layer with the effort (default `"n_records"`, else the
#'   first).
#' @param method `"half"`, `"log"` or `"linear"`.
#' @param half The effort `O` that halves ignorance (`method = "half"`).
#' @param index Use the species observation index `n_records / n_species`
#'   (needs both layers).
#' @param radius_km Optional distance (km) within which effort counts for
#'   neighbouring cells.
#' @param suitability Optional suitability raster (or path) to annotate.
#' @param threshold Ignorance above which a cell counts as not known.
#' @param suitability_threshold Suitability at or above which a cell counts
#'   as suitable (for `class`).
#' @return A list of class `ignorance_map`: `raster` (layer `ignorance`
#'   and, with `suitability`, `suitability`, `suitability_known` and
#'   `class`), `summary` (share of cells and of the suitable area above the
#'   ignorance threshold) and `settings`.
#' @references Tessarolo, G., Ladle, R.J., Lobo, J.M., Rangel, T.F.,
#'   Hortal, J. (2021) Using maps of biogeographical ignorance to reveal the
#'   uncertainty in distributional data hidden in species distribution
#'   models. *Ecography* 44, 1743-1755. \doi{10.1111/ecog.05793}
#'
#'   Ruete, A. (2015) Displaying bias in sampling effort of data accessed
#'   from biodiversity databases using ignorance maps. *Biodiversity Data
#'   Journal* 3, e5361. \doi{10.3897/BDJ.3.e5361}
#' @seealso [gbif_effort_rasters()], [plot_ignorance()]
#' @export
#' @examples
#' tab <- data.frame(lat_bin = rep(445:454, each = 10),
#'                   lon_bin = rep(95:104, 10),
#'                   n_records = rpois(100, 20) * rep(0:1, 50),
#'                   n_species = 5, n_target = 0)
#' eff <- gbif_effort_rasters(tab, c(9.5, 10.5, 44.5, 45.5))
#' ig <- ignorance_map(eff, half = 10)
#' ig
#' suit <- terra::init(eff[[1]], "x") - 9.5
#' ig2 <- ignorance_map(eff, half = 10, suitability = suit)
#' plot_ignorance(ig2)
ignorance_map <- function(effort, layer = "n_records",
                          method = c("half", "log", "linear"), half = 1,
                          index = FALSE, radius_km = NULL,
                          suitability = NULL, threshold = 0.5,
                          suitability_threshold = 0.5) {
  method <- match.arg(method)
  if (is.character(effort)) effort <- terra::rast(effort)
  if (!inherits(effort, "SpatRaster")) {
    stop("`effort` must be a SpatRaster or a raster path.", call. = FALSE)
  }
  if (!is.numeric(half) || half <= 0) stop("`half` must be positive.",
                                           call. = FALSE)
  if (index) {
    if (!all(c("n_records", "n_species") %in% names(effort))) {
      stop("index = TRUE needs layers n_records and n_species.",
           call. = FALSE)
    }
    N <- terra::ifel(effort[["n_species"]] > 0,
                     effort[["n_records"]] / effort[["n_species"]], 0)
  } else {
    N <- if (layer %in% names(effort)) effort[[layer]] else effort[[1]]
  }
  if (isTRUE(terra::global(N, "min", na.rm = TRUE)[1, 1] < 0)) {
    stop("Effort cannot be negative.", call. = FALSE)
  }
  if (!is.null(radius_km) && radius_km > 0) {
    N <- terra::focal(N, radius_kernel(N, radius_km), fun = "sum",
                      na.rm = TRUE)
  }
  ign <- switch(method,
    half = half / (N + half),
    log = {
      nmax <- terra::global(N, "max", na.rm = TRUE)[1, 1]
      if (!is.finite(nmax) || nmax <= 0) N * 0 + 1 else
        1 - log(N + 1) / log(nmax + 1)
    },
    linear = {
      nmax <- terra::global(N, "max", na.rm = TRUE)[1, 1]
      if (!is.finite(nmax) || nmax <= 0) N * 0 + 1 else 1 - N / nmax
    })
  names(ign) <- "ignorance"
  out <- ign
  suit_summary <- NULL
  if (!is.null(suitability)) {
    if (is.character(suitability)) suitability <- terra::rast(suitability)
    s <- suitability[[1]]
    g <- if (isTRUE(terra::compareGeom(ign, s, stopOnError = FALSE))) ign
      else terra::project(ign, s, method = "near")
    known <- terra::ifel(g <= threshold, s, NA)
    cls <- (s >= suitability_threshold) * 2 + (g > threshold) + 1
    levels(cls) <- data.frame(id = 1:4, class = c(
      "unsuitable, known", "unsuitable, ignorant", "suitable, known",
      "suitable, ignorant"))
    out <- c(g, s, known, cls)
    names(out) <- c("ignorance", "suitability", "suitability_known",
                    "class")
    fr <- terra::freq(cls)
    cnt <- stats::setNames(rep(0, 4), c("unsuitable, known",
                                        "unsuitable, ignorant",
                                        "suitable, known",
                                        "suitable, ignorant"))
    if (!is.null(fr) && nrow(fr)) cnt[as.character(fr$value)] <- fr$count
    suit_summary <- data.frame(class = names(cnt), cells = as.vector(cnt),
                               share = as.vector(cnt) / max(1, sum(cnt)))
  }
  v <- terra::freq(ign > threshold)
  n_ig <- if (!is.null(v) && any(v$value == 1)) v$count[v$value == 1] else 0
  n_all <- if (!is.null(v)) sum(v$count) else 0
  structure(list(raster = out,
                 summary = list(share_ignorant = n_ig / max(1, n_all),
                                cells = n_all, classes = suit_summary),
                 settings = list(method = method, half = half,
                                 index = index, radius_km = radius_km,
                                 threshold = threshold,
                                 suitability_threshold =
                                   suitability_threshold)),
            class = "ignorance_map")
}

# Circular moving window of `radius_km` as a 0/1 matrix of cells (internal);
# for longitude/latitude the cell size is taken at the middle latitude
radius_kernel <- function(r, radius_km) {
  res <- terra::res(r)
  if (terra::is.lonlat(r)) {
    lat <- mean(c(terra::ymin(r), terra::ymax(r))) * pi / 180
    dx <- res[1] * 111.32 * cos(lat)
    dy <- res[2] * 110.57
  } else {
    u <- terra::linearUnits(r)
    if (!is.finite(u) || u <= 0) u <- 1
    dx <- res[1] * u / 1000
    dy <- res[2] * u / 1000
  }
  nx <- floor(radius_km / dx)
  ny <- floor(radius_km / dy)
  if (2 * nx + 1 > 1001 || 2 * ny + 1 > 1001) {
    stop("`radius_km` spans more than 500 cells; aggregate the effort ",
         "raster first.", call. = FALSE)
  }
  g <- outer((-ny:ny) * dy, (-nx:nx) * dx, function(a, b) sqrt(a^2 + b^2))
  (g <= radius_km) * 1
}

#' @export
print.ignorance_map <- function(x, ...) {
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  how <- switch(s$method, half = paste0("half-ignorance, O = ", s$half),
                log = "log-normalised", linear = "normalised")
  cat("Biogeographical ignorance (", how, if (s$index)
    "; records per species", if (!is.null(s$radius_km))
      paste0("; effort within ", s$radius_km, " km"), ")\n\nIN SHORT\n",
    sep = "")
  sh <- x$summary$share_ignorant
  if (sh == 0) {
    say(paste0("No cell has an ignorance above ", s$threshold, ": with this ",
               "effort and `half`, the whole area counts as recorded. Check ",
               "that `half` suits the group (a larger value is stricter)."))
  } else if (sh == 1) {
    say(paste0("Every cell has an ignorance above ", s$threshold, ": the ",
               "effort layer may be empty, or `half` too large for the ",
               "group."))
  } else {
    say(paste0(share_label(sh), " of the cells have an ignorance above ",
               s$threshold, ": there, few records of the reference group ",
               "were made, so missing records of the species say little."))
  }
  cl <- x$summary$classes
  if (!is.null(cl)) {
    suitable <- cl$cells[3] + cl$cells[4]
    if (suitable == 0) {
      say("No cell is predicted suitable at this threshold.")
    } else if (cl$cells[4] == 0) {
      say("Every cell predicted suitable lies where ignorance is low.")
    } else {
      say(paste0("Of the cells predicted suitable (suitability >= ",
                 s$suitability_threshold, "), ",
                 share_label(cl$cells[4] / suitable), " lie where ignorance ",
                 "is high: those predictions rest on few records and should ",
                 "be flagged as low confidence."))
    }
  }
  cat("", strwrap(paste0(
    "Ignorance 0 = well recorded, 1 = nothing recorded (Ruete 2015). The ",
    "threshold of ", s$threshold, " (half of the ignorance removed) is a ",
    "rule of thumb. Effort counts only: taxonomic reliability, age of ",
    "records and completeness, used by Tessarolo et al. (2021), are not ",
    "included."), width = 78), sep = "\n")
  invisible(x)
}

#' Plot a map of biogeographical ignorance
#'
#' The ignorance map and, when a suitability raster was given to
#' [ignorance_map()], the joint map of suitability and ignorance (suitable
#' or not, known or ignorant). Legend entries give the number and share of
#' cells in each class.
#'
#' @param x A result of [ignorance_map()].
#' @param type `"both"` (side by side, when there is a joint map),
#'   `"ignorance"` or `"joint"`.
#' @param max_cells Largest number of cells drawn.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot, or for `type = "both"` invisibly a list of two.
#' @references Tessarolo, G., Ladle, R.J., Lobo, J.M., Rangel, T.F.,
#'   Hortal, J. (2021) Using maps of biogeographical ignorance to reveal the
#'   uncertainty in distributional data hidden in species distribution
#'   models. *Ecography* 44, 1743-1755. \doi{10.1111/ecog.05793}
#' @export
#' @examples
#' tab <- data.frame(lat_bin = rep(445:454, each = 10),
#'                   lon_bin = rep(95:104, 10),
#'                   n_records = rpois(100, 20) * rep(0:1, 50),
#'                   n_species = 5, n_target = 0)
#' eff <- gbif_effort_rasters(tab, c(9.5, 10.5, 44.5, 45.5))
#' plot_ignorance(ignorance_map(eff, half = 10))
plot_ignorance <- function(x, type = c("both", "ignorance", "joint"),
                           max_cells = 5e5, base_size = 11, file = NULL,
                           width = NULL, height = NULL) {
  if (!inherits(x, "ignorance_map")) {
    stop("`x` must come from ignorance_map().", call. = FALSE)
  }
  type <- match.arg(type)
  joint <- "class" %in% names(x$raster)
  if (type == "joint" && !joint) {
    stop("No suitability was given to ignorance_map().", call. = FALSE)
  }
  if (type == "both" && !joint) type <- "ignorance"
  if (is.null(width)) width <- if (type == "both") 12.5 else 6.5
  if (is.null(height)) height <- 5.5
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  r <- x$raster
  if (terra::ncell(r) > max_cells) {
    r <- terra::spatSample(r, max_cells, method = "regular",
                           as.raster = TRUE)
  }
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  crd <- if (terra::is.lonlat(r)) ggplot2::coord_quickmap() else
    ggplot2::coord_equal()
  s <- x$settings
  p_ig <- ggplot2::ggplot(d[!is.na(d$ignorance), ],
                          ggplot2::aes(.data$x, .data$y,
                                       fill = .data$ignorance)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_distiller(palette = "YlOrBr", direction = 1,
                                  limits = c(0, 1), name = "Ignorance") +
    crd +
    ggplot2::labs(x = NULL, y = NULL, title = "Biogeographical ignorance",
                  subtitle = paste0(switch(s$method,
                    half = paste0("Half-ignorance score, O = ", s$half),
                    log = "Log-normalised effort",
                    linear = "Normalised effort"),
                    "; 0 = well recorded, 1 = nothing recorded.")) +
    theme_enm(base_size, ink, muted)
  p_j <- NULL
  if (joint) {
    dj <- d[!is.na(d$class), ]
    lv <- c("unsuitable, known", "unsuitable, ignorant", "suitable, known",
            "suitable, ignorant")
    dj$class <- factor(as.character(dj$class), lv)
    cl <- x$summary$classes
    lab <- paste0(lv, ": ", cl$cells, " (", share_label(cl$share), ")")
    pal <- c("#d9d9d9", "#fdd0a2", "#2171b5", "#a50f15")
    p_j <- ggplot2::ggplot(dj, ggplot2::aes(.data$x, .data$y,
                                            fill = .data$class)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_manual(values = stats::setNames(pal, lv),
                                 labels = lab, drop = FALSE, name = NULL) +
      crd +
      ggplot2::labs(x = NULL, y = NULL,
                    title = "Suitability and ignorance",
                    subtitle = paste0("Suitable: suitability >= ",
                                      s$suitability_threshold,
                                      "; ignorant: ignorance > ",
                                      s$threshold, ".")) +
      theme_enm(base_size, ink, muted)
  }
  p <- switch(type, ignorance = p_ig, joint = p_j,
              list(ignorance = p_ig, joint = p_j))
  if (!is.null(file)) save_figure(p, file, width, height)
  if (type == "both") {
    print_side_by_side(p)
    return(invisible(p))
  }
  p
}
