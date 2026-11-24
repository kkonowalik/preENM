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
  ignorance_object(ign, suitability, threshold, suitability_threshold,
                   list(method = method, half = half, index = index,
                        radius_km = radius_km))
}

# Ignorance raster, optionally with a suitability raster, as an
# `ignorance_map` object (internal; shared by ignorance_map() and
# record_ignorance())
ignorance_object <- function(ign, suitability, threshold,
                             suitability_threshold, settings) {
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
  settings$threshold <- threshold
  settings$suitability_threshold <- suitability_threshold
  structure(list(raster = out,
                 summary = list(share_ignorant = n_ig / max(1, n_all),
                                cells = n_all, classes = suit_summary),
                 settings = settings),
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
                log = "log-normalised", linear = "normalised",
                records = paste0("from ", s$n_records, " records: ",
                                 paste(s$components, collapse = ", ")))
  cat("Biogeographical ignorance (", how, if (isTRUE(s$index))
    "; records per species", if (!is.null(s$radius_km))
      paste0(if (s$method == "records") "; interpolated within " else
        "; effort within ", s$radius_km, " km"), ")\n\nIN SHORT\n",
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
               s$threshold, ": there, the records of the reference group ",
               "are few", if (s$method == "records")
                 ", old, imprecise or from incomplete inventories", ", so ",
               "missing records of the species say little."))
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
    "Ignorance 0 = well recorded, 1 = nothing recorded",
    if (s$method == "records") paste0(
      " (Tessarolo et al. 2021: per record, the components are added and ",
      "rescaled to 0-1 over all records; a cell takes its best record of ",
      "the focal species, or the mean over the other species; empty cells ",
      "are interpolated within ", s$radius_km, " km).") else
      " (Ruete 2015).", " The threshold of ", s$threshold,
    " is a rule of thumb.", if (s$method != "records") paste0(
      " Effort counts only: taxonomic reliability, age of records and ",
      "completeness (Tessarolo et al. 2021) can be added from the records ",
      "with record_ignorance().")), width = 78), sep = "\n")
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
                    linear = "Normalised effort",
                    records = paste0("From records (", paste(s$components,
                                                            collapse = ", "),
                                     ")")),
                    ";\n0 = well recorded, 1 = nothing recorded.")) +
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

#' Map of biogeographical ignorance from occurrence records
#'
#' Builds a map of biogeographical ignorance from the records themselves
#' (e.g. a GBIF download of the reference group), following the components
#' of Tessarolo et al. (2021). Each record gets a quality from the
#' components available:
#'
#' * **Completeness** of the inventory of its cell, from the species
#'   accumulation curve: `1 - f1 / N`, where `N` is the number of records of
#'   the group in the cell and `f1` the number of species recorded only
#'   once. `f1 / N` is the expected slope at the end of a randomised
#'   accumulation curve, the quantity Tessarolo et al. averaged over 100
#'   randomisations. Cells with fewer than `min_records` records have
#'   completeness 0 (Tessarolo et al. used 50 individuals).
#' * **Age** (`year`): a Gaussian decay with the time since the record, 1
#'   for the most recent record; `half_life` is the age at which it falls
#'   to 0.5 (default: half of the period spanned by the records, so the
#'   oldest records are the least informative, as in Tessarolo et al.).
#' * **Taxonomic reliability** (`reliability`): a value from 0 to 1 per
#'   record, e.g. by identifier class (Tessarolo et al. used 1, 0.9 and 0.75
#'   for taxonomists, experts and amateurs). Give a column, or a named
#'   vector translating the values of a column (e.g. `basisOfRecord` or a
#'   class of identifier), or one weight per species from the taxonomic
#'   complexity of its name ([taxon_complexity()]: the number of synonyms
#'   and of subspecies and varieties in World Flora Online or another
#'   backbone). That proxy is not published or validated; if you use it,
#'   say so.
#' * **Precision** (`uncertainty`, e.g. `coordinateUncertaintyInMeters`):
#'   not in Tessarolo et al.; `min(1, (cell size / (2 x uncertainty))^2)`,
#'   roughly the share of the uncertainty circle that falls in the cell.
#'
#' As in Tessarolo et al., the components of each record are added and
#' rescaled to 0-1 over all records, and ignorance is 1 minus that. A cell
#' with records of the focal `species` takes the ignorance of its best one;
#' a cell with records of other species only, the mean over species of
#' their best record (records of the group there make an absence of the
#' focal species informative). Cells without records are interpolated by
#' inverse squared distance from the cells with records within
#' `radius_km` (Tessarolo et al. used the range of spatial autocorrelation
#' of the predictors, see [suggest_block_size()], and a distance combining
#' geography and environment; here it is geographic only); farther cells
#' get ignorance 1.
#'
#' @param records Records of the reference group (including the focal
#'   species): a data.frame or CSV path (e.g. a GBIF download), with
#'   coordinates and the other columns named below.
#' @param template A raster defining the grid (e.g. 10 km cells; the result
#'   has one value per cell). Keep it coarse: every cell is interpolated.
#' @param species Optional name of the focal species (a value of the
#'   species column).
#' @param species_col Column with species names (default: `species`, else
#'   `scientificName`).
#' @param year Column with the year of collection (default `year` if
#'   present; `NULL` or `FALSE` to ignore age).
#' @param half_life Age (years) at which the age component falls to 0.5.
#' @param reliability Optional taxonomic reliability: a column with values
#'   0-1, a named vector translating the values of `reliability_col`, or a
#'   [taxon_complexity()] result (one weight per species, from its number of
#'   synonyms and infraspecific taxa).
#' @param reliability_col Column translated by a named `reliability`.
#' @param uncertainty Column with the coordinate uncertainty in metres
#'   (default `coordinateUncertaintyInMeters` if present; `FALSE` to
#'   ignore).
#' @param min_records Fewest records in a cell for its completeness.
#' @param radius_km Distance within which empty cells are interpolated.
#' @inheritParams ignorance_map
#' @param coords,crs Coordinate columns and their CRS (detected by default).
#' @return An `ignorance_map` (see [ignorance_map()]); `settings` lists the
#'   components used, and `records` the quality of each record.
#' @references Tessarolo, G., Ladle, R.J., Lobo, J.M., Rangel, T.F.,
#'   Hortal, J. (2021) Using maps of biogeographical ignorance to reveal the
#'   uncertainty in distributional data hidden in species distribution
#'   models. *Ecography* 44, 1743-1755. \doi{10.1111/ecog.05793}
#' @seealso [ignorance_map()], [plot_ignorance()]
#' @export
#' @examples
#' set.seed(1)
#' rec <- data.frame(lon = runif(600, 10, 12), lat = runif(600, 50, 51),
#'                   species = sample(paste0("sp", 1:30), 600, TRUE),
#'                   year = sample(1950:2020, 600, TRUE))
#' rec <- rec[rec$lon < 11.5, ]                  # the east is unrecorded
#' grid <- terra::rast(xmin = 10, xmax = 12, ymin = 50, ymax = 51,
#'                     resolution = 0.1, crs = "EPSG:4326")
#' ig <- record_ignorance(rec, grid, species = "sp1", min_records = 5,
#'                        radius_km = 20)
#' ig
#' plot_ignorance(ig)
record_ignorance <- function(records, template, species = NULL,
                             species_col = NULL, year = NULL,
                             half_life = NULL, reliability = NULL,
                             reliability_col = NULL, uncertainty = NULL,
                             min_records = 50, radius_km = 50,
                             suitability = NULL, threshold = 0.5,
                             suitability_threshold = 0.5, coords = NULL,
                             crs = "EPSG:4326") {
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  df <- as.data.frame(records, stringsAsFactors = FALSE)
  r <- as_grid(template)[[1]]
  if (terra::ncell(r) > 5e6) {
    warning("The template has ", terra::ncell(r), " cells; ignorance is ",
            "interpolated cell by cell, so a coarser grid (e.g. 10 km) is ",
            "faster and closer to the scale of recording.", call. = FALSE)
  }
  nms <- names(df)
  sp_col <- if (is.null(species_col)) {
    find_column(nms, c("species", "scientificName", "taxon"))
  } else species_col
  if (is.na(sp_col) || !sp_col %in% nms) {
    stop("No species column; give `species_col`.", call. = FALSE)
  }
  xy_cols <- if (is.null(coords)) guess_coords(nms) else coords
  lon <- as_number(df[[xy_cols[1]]])
  lat <- as_number(df[[xy_cols[2]]])
  ok <- !is.na(lon) & !is.na(lat) & !is.na(df[[sp_col]])
  df <- df[ok, , drop = FALSE]
  p <- terra::vect(cbind(lon[ok], lat[ok]), crs = crs)
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  df$.cell <- terra::cellFromXY(r, terra::crds(p))
  df <- df[!is.na(df$.cell), , drop = FALSE]
  if (!nrow(df)) stop("No record falls on the template.", call. = FALSE)
  sp <- as.character(df[[sp_col]])
  comp <- list()

  # completeness of each cell: 1 - f1 / N (expected final slope of the
  # randomised species accumulation curve)
  cell_n <- tabulate(match(df$.cell, unique(df$.cell)))
  ucell <- unique(df$.cell)
  f1 <- vapply(ucell, function(cc) {
    t <- table(sp[df$.cell == cc])
    sum(t == 1)
  }, numeric(1))
  compl <- ifelse(cell_n >= min_records, 1 - f1 / cell_n, 0)
  comp$completeness <- compl[match(df$.cell, ucell)]

  # age
  yr_col <- if (is.null(year)) find_column(nms, "year") else if
    (isFALSE(year)) NA_character_ else year
  if (!is.na(yr_col) && yr_col %in% nms) {
    yr <- as_number(df[[yr_col]])
    if (any(!is.na(yr))) {
      age <- max(yr, na.rm = TRUE) - yr
      hl <- if (is.null(half_life)) max(1, diff(range(yr, na.rm = TRUE)) / 2)
        else half_life
      w <- exp(-log(2) * (age / hl)^2)
      w[is.na(w)] <- 0
      comp$age <- w
    }
  }

  # taxonomic reliability
  if (!is.null(reliability)) {
    if (is.data.frame(reliability)) {
      # taxon_complexity(): one weight per species name
      if (!all(c("name", "reliability") %in% names(reliability))) {
        stop("A data.frame `reliability` needs columns name and ",
             "reliability (see taxon_complexity()).", call. = FALSE)
      }
      rel <- reliability$reliability[match(sp, reliability$name)]
      if (any(is.na(rel))) {
        message(sum(is.na(rel)), " record(s) of species not found in the ",
                "backbone get reliability 1.")
        rel[is.na(rel)] <- 1
      }
    } else if (is.character(reliability) && length(reliability) == 1 &&
        reliability %in% nms) {
      rel <- as_number(df[[reliability]])
    } else if (is.numeric(reliability) && !is.null(names(reliability))) {
      if (is.null(reliability_col) || !reliability_col %in% nms) {
        stop("A named `reliability` needs `reliability_col`.",
             call. = FALSE)
      }
      rel <- unname(reliability[as.character(df[[reliability_col]])])
    } else {
      stop("`reliability` must be a column name or a named vector.",
           call. = FALSE)
    }
    if (any(rel < 0 | rel > 1, na.rm = TRUE)) {
      stop("Reliability values must lie between 0 and 1.", call. = FALSE)
    }
    rel[is.na(rel)] <- 0
    comp$reliability <- rel
  }

  # precision of the coordinates
  un_col <- if (is.null(uncertainty)) {
    find_column(nms, c("coordinateUncertaintyInMeters", "uncertainty"))
  } else if (isFALSE(uncertainty)) NA_character_ else uncertainty
  if (!is.na(un_col) && un_col %in% nms) {
    u <- as_number(df[[un_col]])
    if (any(!is.na(u))) {
      cs <- sqrt(stats::median(terra::extract(terra::cellSize(
        terra::rast(r), unit = "m"), unique(df$.cell))[[1]], na.rm = TRUE))
      prec <- pmin(1, (cs / (2 * u))^2)
      prec[is.na(prec)] <- 1          # no stated uncertainty: not penalised
      comp$precision <- prec
    }
  }

  # record quality: components added, rescaled 0-1; ignorance = 1 - that
  tot <- Reduce(`+`, comp)
  rg <- range(tot)
  q <- if (diff(rg) > 0) (tot - rg[1]) / diff(rg) else rep(1, length(tot))
  ig_rec <- 1 - q

  # one value per cell with records
  if (!is.null(species) && !species %in% sp) {
    warning("No record of `species` (", species, "); every cell uses the ",
            "records of the other species.", call. = FALSE)
  }
  cell_val <- vapply(ucell, function(cc) {
    i <- which(df$.cell == cc)
    if (!is.null(species) && any(sp[i] == species)) {
      return(min(ig_rec[i][sp[i] == species]))
    }
    mean(tapply(ig_rec[i], sp[i], min))
  }, numeric(1))

  # interpolation of empty cells (inverse squared distance within radius)
  v <- rep(NA_real_, terra::ncell(r))
  v[ucell] <- cell_val
  empty <- setdiff(seq_len(terra::ncell(r)), ucell)
  if (length(empty)) {
    pc <- terra::vect(terra::xyFromCell(r, ucell), crs = terra::crs(r))
    pe <- terra::vect(terra::xyFromCell(r, empty), crs = terra::crs(r))
    Xc <- dist_space(pc)
    Xe <- dist_space(pe)
    if (length(ucell) > 20) {
      nn <- dbscan::kNN(Xc$X, k = 20, query = Xe$X)
    } else {
      # few cells with records: all distances
      dd <- outer(rowSums(Xe$X^2), rowSums(Xc$X^2), "+") -
        2 * Xe$X %*% t(Xc$X)
      dd[] <- sqrt(pmax(dd, 0))
      nn <- list(dist = dd, id = matrix(seq_along(ucell), nrow(dd),
                                        length(ucell), byrow = TRUE))
    }
    d <- matrix(nn$dist, nrow = length(empty))
    if (Xc$lonlat) d[] <- 2 * 6371.0088 * asin(pmin(d / (2 * 6371.0088), 1))
    w <- (d <= radius_km) / pmax(d, 1e-9)^2
    val <- matrix(cell_val[nn$id], nrow = length(empty))
    sw <- rowSums(w)
    v[empty] <- ifelse(sw > 0, rowSums(w * val) / sw, 1)
  }
  ign <- terra::rast(r)
  terra::values(ign) <- v
  if (terra::hasValues(r)) ign <- terra::mask(ign, r)    # e.g. land only
  names(ign) <- "ignorance"
  out <- ignorance_object(ign, suitability, threshold, suitability_threshold,
                          list(method = "records", components = names(comp),
                               n_records = nrow(df), radius_km = radius_km,
                               species = species, index = FALSE))
  out$records <- data.frame(cell = df$.cell, species = sp,
                            as.data.frame(comp), ignorance = ig_rec,
                            stringsAsFactors = FALSE)
  out
}
