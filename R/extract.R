# Reading and harmonising inputs ------------------------------------------------

#' Read model predictions as a SpatRaster
#'
#' Accepts the many ways a user may hold model output and returns a single
#' multi-layer [terra::SpatRaster] with one named layer per model.
#'
#' @param models One of:
#'   * a [terra::SpatRaster] (each layer is treated as one model),
#'   * a character vector of file paths (e.g. GeoTIFF / ASCII grids exported
#'     from Maxent, biomod2, ENMeval, QGIS, ...),
#'   * a (named) list mixing `SpatRaster` objects and file paths,
#'   * a legacy `raster::RasterLayer`/`RasterStack` (converted with
#'     [terra::rast()]).
#' @return A named list of single-layer `SpatRaster`s.
#' @keywords internal
#' @noRd
as_model_list <- function(models) {
  if (inherits(models, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    models <- terra::rast(models)
  }
  if (inherits(models, "SpatRaster")) {
    out <- lapply(seq_len(terra::nlyr(models)), function(i) models[[i]])
    names(out) <- names(models)
  } else if (is.character(models)) {
    nms <- names(models)
    if (is.null(nms)) nms <- sub("\\.[^.]*$", "", basename(models))
    out <- lapply(models, terra::rast)
    out <- lapply(out, function(r) {
      if (terra::nlyr(r) > 1) {
        warning("A file has several layers; only the first one is used.",
                call. = FALSE)
        r <- r[[1]]
      }
      r
    })
    names(out) <- nms
  } else if (is.list(models)) {
    out <- lapply(models, function(m) {
      if (is.character(m)) m <- terra::rast(m)
      if (inherits(m, c("RasterLayer", "RasterStack", "RasterBrick"))) {
        m <- terra::rast(m)
      }
      if (!inherits(m, "SpatRaster")) {
        stop("Every element of `models` must be a SpatRaster or a file path.",
             call. = FALSE)
      }
      m[[1]]
    })
    nms <- names(models)
    if (is.null(nms) || any(nms == "")) {
      auto <- vapply(out, function(r) names(r)[1], character(1))
      nms <- if (is.null(nms)) auto else ifelse(nms == "", auto, nms)
    }
    names(out) <- nms
  } else {
    stop("`models` must be a SpatRaster, file path(s) or a list of them.",
         call. = FALSE)
  }
  if (anyDuplicated(names(out))) names(out) <- make.unique(names(out))
  out
}

#' Read occurrence points as a SpatVector
#'
#' @param x A data.frame / matrix with coordinates, a path to a CSV file,
#'   an `sf` object or a [terra::SpatVector] of points.
#' @param coords Optional character vector of length 2 giving the names of
#'   the x (longitude) and y (latitude) columns. If `NULL`, common names such
#'   as `lon/lat`, `longitude/latitude`, `x/y` or `decimalLongitude/
#'   decimalLatitude` (GBIF) are detected automatically.
#' @param crs Coordinate reference system of tabular coordinates. Defaults to
#'   WGS84 longitude/latitude (`"EPSG:4326"`).
#' @keywords internal
#' @noRd
as_points <- function(x, coords = NULL, crs = "EPSG:4326") {
  if (inherits(x, "SpatVector")) return(x)
  if (inherits(x, "sf")) return(terra::vect(x))
  if (inherits(x, "SpatialPointsDataFrame") || inherits(x, "SpatialPoints")) {
    return(terra::vect(x))
  }
  if (is.character(x) && length(x) == 1) {
    x <- utils::read.csv(x, check.names = FALSE)
  }
  if (is.matrix(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) {
    stop("Points must be a data.frame, matrix, CSV path, sf or SpatVector.",
         call. = FALSE)
  }
  if (is.null(coords)) coords <- guess_coords(names(x))
  if (!all(coords %in% names(x))) {
    stop("Coordinate columns not found: ",
         paste(setdiff(coords, names(x)), collapse = ", "), call. = FALSE)
  }
  xy <- data.frame(x = as.numeric(x[[coords[1]]]),
                   y = as.numeric(x[[coords[2]]]))
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  terra::vect(xy, geom = c("x", "y"), crs = crs)
}

guess_coords <- function(nms) {
  low <- tolower(nms)
  # standard lon/lat first (e.g. merge_occurrences() output, which may also
  # keep the original GBIF columns), then common alternatives
  candidates <- list(
    c("lon", "lat"),
    c("longitude", "latitude"),
    c("decimallongitude", "decimallatitude"),
    c("long", "lat"),
    c("lng", "lat"),
    c("x_long", "y_lat"),
    c("x", "y")
  )
  for (cand in candidates) {
    if (all(cand %in% low)) return(nms[match(cand, low)])
  }
  stop("Could not detect coordinate columns. Please supply ",
       "`coords = c(\"<x/longitude column>\", \"<y/latitude column>\")`.\n",
       "Available columns: ", paste(nms, collapse = ", "), call. = FALSE)
}

#' Extract model predictions at presence and absence/background points
#'
#' Builds the observation/prediction table that all evaluation metrics are
#' computed from. Points are re-projected to each raster's CRS when needed.
#'
#' @param models Model predictions: a [terra::SpatRaster] (one layer per
#'   model), file path(s) to rasters, or a named list of these.
#' @param presence Presence points (data.frame, CSV path, `sf` or
#'   `SpatVector`).
#' @param absence Absence points or, for presence-only data, background /
#'   pseudo-absence points (same formats as `presence`).
#' @param coords Optional names of the x/longitude and y/latitude columns
#'   when points are tabular. Detected automatically if `NULL`.
#' @param crs CRS of tabular coordinates (default WGS84, `"EPSG:4326"`).
#' @param na_action How to handle points that fall on `NA` cells.
#'   `"common"` (default) drops a point for *all* models if it is `NA` in any
#'   model, so every model is evaluated on exactly the same points (a fair
#'   comparison). `"per_model"` keeps points independently for each model
#'   (this reproduces the original 2021 script).
#' @param rescale Bring predictions to the 0-1 scale. Algorithms write
#'   their output on different scales (0-1, 0-100, 0-1000); the scale of
#'   each model is recognised from its largest value, and models on 0-100
#'   or 0-1000 are divided by 100 or 1000, with a message. `FALSE` keeps
#'   the values and warns when they are not on the 0-1 scale. Models whose
#'   scale cannot be recognised (negative values or above 1000, e.g. raw
#'   Maxent output) are left as they are, with a warning.
#'
#' @return A data.frame with columns `obs` (1 = presence, 0 = absence) and one
#'   column of predicted values per model. Attribute `"n_dropped"` holds the
#'   number of points removed because of `NA` predictions, and `"scales"`
#'   the recognised scale of each model and what it was divided by.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models, sim$presence, sim$absence)
#' head(tab)
extract_predictions <- function(models, presence, absence, coords = NULL,
                                crs = "EPSG:4326",
                                na_action = c("common", "per_model"),
                                rescale = TRUE) {
  na_action <- match.arg(na_action)
  mods <- as_model_list(models)
  pres <- as_points(presence, coords, crs)
  abs <- as_points(absence, coords, crs)
  pts <- rbind(pres, abs)
  obs <- c(rep(1L, nrow(pres)), rep(0L, nrow(abs)))

  vals <- lapply(mods, function(r) {
    p <- pts
    if (!is.na(terra::crs(r)) && terra::crs(r) != "" &&
        !is.na(terra::crs(p)) && terra::crs(p) != "" &&
        !terra::same.crs(r, p)) {
      p <- terra::project(p, r)
    }
    terra::extract(r, p, ID = FALSE)[[1]]
  })
  tab <- data.frame(obs = obs, vals, check.names = FALSE)
  names(tab) <- c("obs", names(mods))

  n_before <- nrow(tab)
  if (na_action == "common") {
    tab <- tab[stats::complete.cases(tab), , drop = FALSE]
  }
  dropped <- n_before - nrow(tab)
  if (dropped > 0) {
    message(dropped, " point(s) fell on NA cells (e.g. outside the raster) ",
            "and were removed.")
  }
  rownames(tab) <- NULL
  tab <- rescale_table(tab, rescale)
  attr(tab, "n_dropped") <- dropped
  tab
}

check_prediction_scale <- function(preds) {
  rng <- range(unlist(preds), na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    warning("Some predictions lie outside [0, 1] (range: ",
            signif(rng[1], 3), " to ", signif(rng[2], 3), "). MAE and bias ",
            "assume probability-like suitability values; consider rescaling ",
            "(e.g. Maxent 'cloglog'/'logistic' output).", call. = FALSE)
  }
  invisible(NULL)
}

# Scale of each model's predictions (internal): 1 (0-1), 100 (0-100) or
# 1000 (0-1000), judged from the largest value; NA when it cannot be told
# (negative values, or larger than 1000, e.g. raw Maxent output).
prediction_scales <- function(preds) {
  do.call(rbind, lapply(names(preds), function(m) {
    v <- preds[[m]]
    v <- v[is.finite(v)]
    lo <- if (length(v)) min(v) else NA_real_
    hi <- if (length(v)) max(v) else NA_real_
    sc <- if (!length(v) || lo < 0) NA_real_ else if (hi <= 1) 1 else
      if (hi <= 100) 100 else if (hi <= 1000) 1000 else NA_real_
    data.frame(model = m, min = lo, max = hi, scale = sc,
               stringsAsFactors = FALSE)
  }))
}

# Check, and with rescale = TRUE harmonise, the scales of several models
# (internal). `preds` is a named list (or data.frame) of numeric vectors, or
# of values sampled from rasters. Returns the divisor per model (1 = keep)
# and the scale table; messages and warnings explain what was found.
harmonise_scales <- function(preds, rescale = TRUE, what = "predictions") {
  sc <- prediction_scales(preds)
  div <- stats::setNames(rep(1, nrow(sc)), sc$model)
  unknown <- sc$model[is.na(sc$scale)]
  wide <- sc[!is.na(sc$scale) & sc$scale > 1, , drop = FALSE]
  if (nrow(wide)) {
    lab <- paste0(wide$model, " (0-", wide$scale, ", max ",
                  signif(wide$max, 3), ")")
    if (rescale) {
      div[wide$model] <- wide$scale
      message("The ", what, " of ", paste(lab, collapse = ", "),
              " are not on the 0-1 scale; they were divided by ",
              paste(unique(wide$scale), collapse = " or "),
              " to the 0-1 scale of the other models. Check that this is ",
              "right (rescale = FALSE keeps the original values).")
    } else {
      warning("The ", what, " of ", paste(lab, collapse = ", "),
              " are not on the 0-1 scale", if (any(sc$scale %in% 1))
                ", unlike those of the other models", ". Metrics such as ",
              "MAE, and ensembles, are then not comparable; use ",
              "rescale = TRUE or rescale them yourself.", call. = FALSE)
    }
  }
  if (length(unknown)) {
    warning("The scale of the ", what, " of ", paste(unknown, collapse = ", "),
            " cannot be recognised (negative values or larger than 1000, ",
            "e.g. raw Maxent output); they were left as they are. Use ",
            "probability-like output (e.g. Maxent cloglog) on the 0-1 scale.",
            call. = FALSE)
  }
  sc$divided_by <- unname(div[sc$model])
  list(divisor = div, table = sc)
}

# Divisor bringing one set of predictions to 0-1 (internal, no message):
# 1, 100 or 1000; 1 also when the scale cannot be recognised
scale_divisor <- function(v) {
  sc <- prediction_scales(list(x = as.vector(v)))$scale
  if (is.na(sc)) 1 else sc
}

# Apply harmonise_scales() to the prediction columns of a table (internal)
rescale_table <- function(tab, rescale = TRUE, what = "predictions") {
  mods <- setdiff(names(tab), "obs")
  h <- harmonise_scales(tab[mods], rescale, what)
  for (m in mods) if (h$divisor[[m]] != 1) tab[[m]] <- tab[[m]] / h$divisor[[m]]
  attr(tab, "scales") <- h$table
  tab
}
