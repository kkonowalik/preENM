# Selection of predictor variables across periods and regions -----------------

#' Sample environmental values across periods, scenarios and regions
#'
#' The first step of the variable selection in [select_variables()]. Draws
#' random cells from each data set, i.e. the same predictors for different
#' periods (present, Last Glacial Maximum, ...) or climate scenarios (e.g.
#' 2070 SSP5-8.5), and reads their values. Collinearity among predictors can
#' change between periods; a set chosen only on present-day data may be
#' collinear, or extrapolated, in the period a model is projected to.
#'
#' **Cells.** Cells with values in all variables are drawn at random from
#' each data set ("own" cells). For longitude/latitude rasters the sample is
#' area-proportional (cells are accepted with probability proportional to
#' cos(latitude)), so correlations are not dominated by the many small cells
#' of high latitudes. When all data sets share one grid, cells valid in all of
#' them ("common" cells) are drawn too: they compare the same places between
#' periods (gradient stability, agreement between data sources, anomalies).
#'
#' **Reading.** Rasters are read block by block, sequentially, so memory use
#' stays low and compressed files are not decompressed once per cell. Clip
#' them to the accessible area first (e.g. [crop_rasters()]).
#'
#' **Anomaly ("delta") data sets** combine a high-resolution present
#' climate with a climate model's change, at the common cells only (no maps
#' needed): `anomalies = list(LGM_delta = c(current = "CHELSA", past = "LGM",
#' model_present = "TraCE_present"))` (see [anomaly_climate()]).
#'
#' **Regions** are the areas within which correlations are also checked:
#' polygons (e.g. the Alps, the Carpathians, Scandinavia) or a grid of tiles.
#'
#' @param datasets Named list of data sets; each a [terra::SpatRaster] with
#'   one layer per variable, a vector of raster files or a folder. The first is
#'   the reference (the calibration period) unless `reference` says otherwise.
#' @param vars Variables to use. Layers are matched by name: exactly, or as a
#'   whole part of the layer name between `_`, `-` or `.` (so `bio01` finds
#'   `CHELSA_TraCE21k_bio01_-100_V1.0`). Default: the layer names of the
#'   reference.
#' @param n Number of cells drawn per data set (and for the common cells).
#' @param reference Name of the reference data set (default: the first).
#' @param regions Optional regions: polygons (`sf`, [terra::SpatVector] or a
#'   path) with names in `region_col`, or a list of tile breaks in the raster
#'   coordinates, e.g. `list(x = c(-30, 5, 30, 90), y = c(30, 60, 90))`.
#' @param region_col Column with the region names (default: the first text
#'   column).
#' @param anomalies Optional named list of anomaly data sets, each
#'   `c(current = , past = , model_present = )` naming data sets.
#' @param ratio_vars,ratio_limits Variables that use a ratio anomaly (e.g.
#'   precipitation) and its limits; see [anomaly_climate()].
#' @param bounds Optional limits per variable applied to anomaly values, e.g.
#'   `list(scd = c(0, 365))` for snow-cover days.
#' @param check_scale Warn about variables whose spread (interquartile
#'   range) differs from the reference by a factor of about 10, i.e. stored
#'   x 10 in one data set (as bio07 in some CHELSA-TraCE21k files).
#' @param fix_scale Divide or multiply such variables by 10 to match the
#'   reference.
#' @param area_weighted Area-proportional sampling for longitude/latitude
#'   rasters.
#' @param seed Optional random seed.
#' @param verbose Print progress.
#' @return An object of class `env_sample` for [select_variables()].
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' now <- c(sim$env, sim$env[[1]] * sim$env[[2]])
#' names(now) <- c("t", "p", "tp")
#' past <- now
#' past[["t"]] <- past[["t"]] - 0.3
#' es <- sample_environment(list(present = now, past = past), n = 2000,
#'                          seed = 1)
#' es
sample_environment <- function(datasets, vars = NULL, n = 1e5,
                               reference = NULL, regions = NULL,
                               region_col = NULL, anomalies = NULL,
                               ratio_vars = NULL, ratio_limits = c(0.1, 10),
                               bounds = NULL, check_scale = TRUE,
                               fix_scale = FALSE, area_weighted = TRUE,
                               seed = NULL, verbose = TRUE) {
  say <- function(...) if (verbose) message(...)
  if (!is.list(datasets) || inherits(datasets, "SpatRaster")) {
    datasets <- list(present = datasets)
  }
  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    stop("`datasets` must be a named list, e.g. list(present = ..., ",
         "LGM = ...).", call. = FALSE)
  }
  ds <- lapply(datasets, as_dataset)
  if (is.null(reference)) reference <- names(ds)[1]
  if (!reference %in% names(ds)) {
    stop("`reference` must be one of the data sets.", call. = FALSE)
  }
  if (is.null(vars)) vars <- names(ds[[reference]])
  ds <- lapply(names(ds), function(k) match_layers(ds[[k]], vars, k))
  names(ds) <- names(datasets)
  same_grid <- all(vapply(ds, function(r) {
    isTRUE(terra::compareGeom(r, ds[[reference]], stopOnError = FALSE))
  }, logical(1)))

  say("Finding cells with values in all ", length(vars), " variables ...")
  valid <- lapply(ds, valid_cells)
  for (k in names(valid)) {
    if (!length(valid[[k]])) stop("No cells with values in all variables in ",
                                  k, ".", call. = FALSE)
  }
  out <- with_seed(seed, {
    own <- lapply(names(ds), function(k) draw_cells(ds[[k]], valid[[k]], n,
                                                    area_weighted))
    names(own) <- names(ds)
    common <- NULL
    if (same_grid && length(ds) > 1) {
      cc <- Reduce(intersect, valid)
      if (length(cc)) {
        common <- draw_cells(ds[[reference]], sort(cc), n, area_weighted)
      } else {
        warning("No cells have values in all data sets: comparisons on the ",
                "same cells are skipped.", call. = FALSE)
      }
    } else if (length(ds) > 1) {
      say("Data sets are on different grids: comparisons on the same cells ",
          "(gradient stability, anomalies) are skipped.")
    }
    list(own = own, common = common)
  })
  say("Reading values ...")
  sets <- list()
  com <- list()
  for (k in names(ds)) {
    cells <- sort(unique(c(out$own[[k]], out$common)))
    v <- read_cells(ds[[k]], cells)
    colnames(v) <- vars
    sets[[k]] <- cell_table(ds[[k]], out$own[[k]], v[match(out$own[[k]], cells),
                                                       , drop = FALSE])
    if (!is.null(out$common)) {
      com[[k]] <- cell_table(ds[[k]], out$common,
                             v[match(out$common, cells), , drop = FALSE])
    }
  }
  if (length(com)) {
    ok <- Reduce(`&`, lapply(com, function(d) stats::complete.cases(d[vars])))
    com <- lapply(com, function(d) d[ok, , drop = FALSE])
  }

  # data sets stored in other units (e.g. x 10)
  if (check_scale && length(ds) > 1) {
    # x 10 storage multiplies the spread by 10, while a change of climate
    # shifts the values but barely changes their spread (a ratio of medians
    # would be meaningless for variables centred near zero)
    ref_iqr <- vapply(vars, function(v) stats::IQR(sets[[reference]][[v]]),
                      numeric(1))
    for (k in setdiff(names(ds), reference)) {
      for (v in vars) {
        rt <- stats::IQR(sets[[k]][[v]]) / ref_iqr[[v]]
        if (!is.finite(rt)) next
        f <- if (rt > 8 && rt < 12) 1 / 10 else if (rt > 1 / 12 && rt < 1 / 8)
          10 else NA
        if (is.na(f)) next
        if (fix_scale) {
          sets[[k]][[v]] <- sets[[k]][[v]] * f
          if (length(com)) com[[k]][[v]] <- com[[k]][[v]] * f
          say(v, " in ", k, ": spread ", round(rt, 1), " x the reference; ",
              "multiplied by ", f)
        } else {
          warning(v, " in ", k, " has a spread (IQR) ", signif(rt, 2),
                  " x that of ",
                  reference, ": stored in other units (x 10)? Use ",
                  "fix_scale = TRUE to rescale.", call. = FALSE)
        }
      }
    }
  }

  # anomaly data sets on the common cells
  for (nm in names(anomalies)) {
    a <- anomalies[[nm]]
    if (!length(com)) stop("Anomalies need data sets on one grid.",
                           call. = FALSE)
    if (!all(c("current", "past", "model_present") %in% names(a)) ||
        !all(a %in% names(ds))) {
      stop("Each anomaly needs current, past and model_present naming data ",
           "sets.", call. = FALSE)
    }
    d <- anomaly_climate(com[[a[["current"]]]][vars], com[[a[["past"]]]][vars],
                         com[[a[["model_present"]]]][vars],
                         ratio_vars = ratio_vars, ratio_limits = ratio_limits)
    for (v in intersect(names(bounds), vars)) {
      d[[v]] <- pmin(pmax(d[[v]], bounds[[v]][1]), bounds[[v]][2])
    }
    keep <- setdiff(names(com[[a[["current"]]]]), vars)
    d <- cbind(d, com[[a[["current"]]]][keep])
    com[[nm]] <- d
    sets[[nm]] <- d
  }

  # regions
  if (!is.null(regions)) {
    for (k in names(sets)) sets[[k]]$.region <- assign_regions(sets[[k]],
                                                              regions,
                                                              region_col,
                                                              ds[[reference]])
    for (k in names(com)) com[[k]]$.region <- assign_regions(com[[k]], regions,
                                                            region_col,
                                                            ds[[reference]])
  }
  res <- list(sets = sets, common = if (length(com)) com else NULL,
              vars = vars, reference = reference,
              anomalies = names(anomalies), same_grid = same_grid,
              lonlat = terra::is.lonlat(ds[[reference]]),
              area_weighted = area_weighted,
              regions = !is.null(regions), n = n)
  class(res) <- "env_sample"
  say("Sample: ", paste(names(sets), vapply(sets, nrow, 1L), sep = " = ",
                        collapse = ", "),
      if (length(com)) paste0("; common cells = ", nrow(com[[1]])))
  res
}

#' @export
print.env_sample <- function(x, ...) {
  cat("Environmental sample: ", length(x$vars), " variables (",
      paste(x$vars, collapse = ", "), ")\n", sep = "")
  cat("Data sets (reference: ", x$reference, "): ",
      paste(names(x$sets), vapply(x$sets, nrow, 1L), sep = " = ",
            collapse = ", "), " cells\n", sep = "")
  cat("Common cells: ", if (is.null(x$common)) "none" else nrow(x$common[[1]]),
      "; area-proportional: ", x$lonlat && x$area_weighted, "\n", sep = "")
  if (x$regions) {
    tab <- table(x$sets[[x$reference]]$.region)
    cat("Regions (cells in the reference): ",
        paste(names(tab), tab, sep = " = ", collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

# --- reading helpers ---------------------------------------------------------------

as_dataset <- function(x) {
  if (inherits(x, "SpatRaster")) return(x)
  if (is.character(x)) {
    if (length(x) == 1 && dir.exists(x)) {
      x <- list.files(x, pattern = "\\.tif{1,2}$", full.names = TRUE,
                      ignore.case = TRUE)
    }
    r <- terra::rast(x)
    names(r) <- sub("\\.tif{1,2}$", "", basename(x), ignore.case = TRUE)
    return(r)
  }
  stop("Each data set must be a SpatRaster, raster files or a folder.",
       call. = FALSE)
}

# Pick and rename the layers holding `vars` (exact name or name token)
match_layers <- function(r, vars, what) {
  nm <- names(r)
  idx <- vapply(vars, function(v) {
    hit <- which(nm == v)
    if (!length(hit)) {
      rx <- paste0("(^|[_.-])", v, "([_.-]|$)")
      hit <- which(grepl(rx, nm, ignore.case = TRUE))
    }
    if (length(hit) != 1) {
      stop("Variable '", v, "' matches ", length(hit), " layers in ", what,
           if (length(hit)) paste0(" (", paste(nm[hit], collapse = ", "), ")"),
           ".", call. = FALSE)
    }
    hit
  }, integer(1))
  out <- r[[idx]]
  names(out) <- vars
  out
}

# Cell numbers with values in every layer (block-wise, sequential)
valid_cells <- function(r) {
  terra::readStart(r)
  on.exit(terra::readStop(r))
  b <- terra::blocks(r)
  nc <- terra::ncol(r)
  out <- vector("list", b$n)
  for (i in seq_len(b$n)) {
    v <- terra::readValues(r, row = b$row[i], nrows = b$nrows[i], col = 1,
                           ncols = nc, mat = TRUE)
    out[[i]] <- (b$row[i] - 1) * nc + which(stats::complete.cases(v))
  }
  unlist(out)
}

# Values of the given (sorted) cells, read block by block
read_cells <- function(r, cells) {
  terra::readStart(r)
  on.exit(terra::readStop(r))
  b <- terra::blocks(r)
  nc <- terra::ncol(r)
  out <- matrix(NA_real_, length(cells), terra::nlyr(r))
  for (i in seq_len(b$n)) {
    lo <- (b$row[i] - 1) * nc + 1
    hi <- (b$row[i] + b$nrows[i] - 1) * nc
    j <- which(cells >= lo & cells <= hi)
    if (!length(j)) next
    v <- terra::readValues(r, row = b$row[i], nrows = b$nrows[i], col = 1,
                           ncols = nc, mat = TRUE)
    out[j, ] <- v[cells[j] - lo + 1, , drop = FALSE]
  }
  out
}

# Random cells, area-proportional for longitude/latitude rasters
draw_cells <- function(r, valid, n, area_weighted) {
  if (!length(valid)) return(integer(0))
  idx <- valid[sample.int(length(valid), min(length(valid), n * 4))]
  if (area_weighted && terra::is.lonlat(r)) {
    w <- cos(terra::yFromCell(r, idx) * pi / 180)
    idx <- idx[stats::runif(length(idx)) < w / max(w)]
  }
  if (length(idx) < n && length(valid) > n) {
    warning("Only ", length(idx), " cells drawn (wanted ", n, "): with ",
            "area-proportional sampling the ", length(valid), " cells give ",
            "fewer; use a smaller `n` or a finer grid.", call. = FALSE)
  }
  sort(idx[seq_len(min(n, length(idx)))])
}

cell_table <- function(r, cells, v) {
  xy <- terra::xyFromCell(r, cells)
  d <- as.data.frame(v)
  d$.x <- xy[, 1]
  d$.y <- xy[, 2]
  d[stats::complete.cases(v), , drop = FALSE]
}

assign_regions <- function(d, regions, region_col, r) {
  if (is.list(regions) && !inherits(regions, c("sf", "sfc")) &&
      all(c("x", "y") %in% names(regions))) {
    ix <- cut(d$.x, regions$x, labels = FALSE)
    iy <- cut(d$.y, regions$y, labels = FALSE)
    lab <- paste0("x", ix, "_y", iy)
    lab[is.na(ix) | is.na(iy)] <- NA
    return(lab)
  }
  g <- if (inherits(regions, "sf")) regions else
    sf::st_as_sf(if (is.character(regions)) terra::vect(regions) else
      if (inherits(regions, "SpatVector")) regions else
        stop("`regions` must be polygons or a list(x = , y = ) of tile breaks.",
             call. = FALSE))
  if (is.null(region_col)) {
    chr <- names(g)[vapply(sf::st_drop_geometry(g), is.character, logical(1))]
    region_col <- if (length(chr)) chr[1] else NULL
  }
  nm <- if (is.null(region_col)) paste0("region", seq_len(nrow(g))) else
    as.character(g[[region_col]])
  g <- sf::st_transform(sf::st_make_valid(g), terra::crs(r))
  p <- sf::st_as_sf(d[c(".x", ".y")], coords = c(".x", ".y"),
                    crs = sf::st_crs(g))
  hit <- sf::st_intersects(p, g)
  vapply(hit, function(h) if (length(h)) nm[h[1]] else NA_character_,
         character(1))
}
