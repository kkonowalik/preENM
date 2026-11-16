# Checks before modelling: sample size, record years, block size -----------------

#' Are there enough records?
#'
#' Compares the number of records with common rules of thumb, after
#' filtering and thinning:
#'
#' * at least 14 records for narrow-ranged and 25 for widespread species,
#'   the minimum numbers found by van Proosdij et al. (2016) for real
#'   landscapes (fewer records gave models no better than chance);
#' * at least `per_predictor` (default 10) records per predictor, so that
#'   responses can be estimated without overfitting;
#' * with `k` folds, at least `per_fold` (default 5) records in a fold on
#'   average, or the metrics of the folds are too noisy.
#'
#' With `env`, the records are counted as occupied cells (several records
#' in a cell add no information to a model fitted on that grid) and the
#' share of the valid cells that is occupied is reported as a lower bound
#' of the prevalence.
#'
#' @param points Records (data.frame, CSV path, `sf` or
#'   [terra::SpatVector]), or their number.
#' @param env Optional predictor rasters you will model with.
#' @param n_predictors Number of predictors (default: the layers of `env`).
#' @param k Optional number of cross-validation folds.
#' @param per_predictor,per_fold Records required per predictor and per
#'   fold.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @return A data.frame with `check`, `value`, `required`, `ok` and `note`.
#' @references van Proosdij, A.S.J., Sosef, M.S.M., Wieringa, J.J., Raes, N.
#'   (2016) Minimum required number of specimen records to develop accurate
#'   species distribution models. *Ecography* 39, 542-552.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 40)
#' check_sample_size(sim$presence, sim$env, k = 5)
#' check_sample_size(18, n_predictors = 6)
check_sample_size <- function(points, env = NULL, n_predictors = NULL,
                              k = NULL, per_predictor = 10, per_fold = 5,
                              coords = NULL, points_crs = "EPSG:4326") {
  prevalence <- NA_real_
  if (is.numeric(points) && length(points) == 1) {
    n <- points
    what <- "records"
  } else {
    if (is.character(points) && length(points) == 1) {
      points <- read_any_csv(points)
    }
    if (is.data.frame(points) && !inherits(points, "sf")) {
      xy <- if (is.null(coords)) guess_coords(names(points)) else coords
      x <- as_number(points[[xy[1]]])
      y <- as_number(points[[xy[2]]])
      ok <- !is.na(x) & !is.na(y)
      p <- terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
    } else {
      p <- as_points(points, coords, points_crs)
    }
    n <- nrow(p)
    what <- "records"
    if (!is.null(env)) {
      r <- as_grid(env)
      if (!terra::same.crs(p, r)) p <- terra::project(p, r)
      cell <- terra::cellFromXY(r, terra::crds(p))
      E <- env_at_cells(r, cell)
      n <- length(unique(cell[stats::complete.cases(E)]))
      what <- "occupied cells with values"
      prevalence <- n / length(valid_cells(r))
      if (is.null(n_predictors)) n_predictors <- terra::nlyr(r)
    }
  }
  rows <- list(
    data.frame(check = paste0("number of ", what), value = n, required = 14,
               ok = n >= 14,
               note = if (n < 14) "too few for a reliable model" else
                 if (n < 25) "enough only for a narrow-ranged species" else
                   "enough also for a widespread species"))
  if (!is.na(prevalence)) {
    rows[[length(rows) + 1]] <- data.frame(
      check = "share of the cells occupied (lower bound of prevalence)",
      value = prevalence, required = NA_real_, ok = NA,
      note = "the more widespread a species, the more records it needs")
  }
  if (!is.null(n_predictors)) {
    v <- n / n_predictors
    rows[[length(rows) + 1]] <- data.frame(
      check = paste0("records per predictor (", n_predictors, ")"),
      value = v, required = per_predictor, ok = v >= per_predictor,
      note = if (v >= per_predictor) "" else paste0(
        "use at most ", max(1, floor(n / per_predictor)),
        " predictors, or simpler models"))
  }
  if (!is.null(k)) {
    v <- n / k
    rows[[length(rows) + 1]] <- data.frame(
      check = paste0("records per fold (k = ", k, ")"), value = v,
      required = per_fold, ok = v >= per_fold,
      note = if (v >= per_fold) "" else paste0(
        "use at most ", max(2, floor(n / per_fold)), " folds"))
  }
  out <- do.call(rbind, rows)
  out$value <- signif(out$value, 3)
  out
}

#' Do the records come from the period of the climate data?
#'
#' Climate layers describe a period (e.g. 1981-2010 for CHELSA V2.1,
#' 1970-2000 for WorldClim 2). Records collected long before or after it
#' may come from places where the climate was different then, or where the
#' species no longer occurs. The function labels every record and prints
#' how many fall inside the period, so that you can decide whether to
#' restrict the records, widen the period with `tolerance`, or choose
#' another climatology.
#'
#' @param points Records (data.frame or CSV path).
#' @param period First and last year of the climate data.
#' @param year Column with the year of collection (default: `year`).
#' @param tolerance Years before and after the period that still count as
#'   inside.
#' @param verbose Print the summary.
#' @return The records with added columns `period_match` (`"within"`,
#'   `"before"`, `"after"` or `"unknown"`) and `years_outside` (0 within;
#'   negative before, positive after). `attr(, "summary")` holds the
#'   counts, shares and median year of each class.
#' @export
#' @examples
#' occ <- data.frame(lon = 1:6, lat = 50,
#'                   year = c(1950, 1985, 1999, 2012, 2021, NA))
#' m <- match_climate_period(occ, c(1981, 2010))
#' m[m$period_match == "within", ]
match_climate_period <- function(points, period = c(1981, 2010), year = NULL,
                                 tolerance = 0, verbose = TRUE) {
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  df <- as.data.frame(points, stringsAsFactors = FALSE)
  if (length(period) != 2 || period[1] > period[2]) {
    stop("`period` must be the first and the last year, e.g. c(1981, 2010).",
         call. = FALSE)
  }
  if (is.null(year)) year <- find_column(names(df), "year")
  if (is.na(year) || !year %in% names(df)) {
    stop("No year column found: name it in `year`.", call. = FALSE)
  }
  y <- as_number(df[[year]])
  lo <- period[1] - tolerance
  hi <- period[2] + tolerance
  off <- ifelse(y < lo, y - lo, ifelse(y > hi, y - hi, 0))
  cls <- ifelse(is.na(y), "unknown",
                ifelse(off < 0, "before", ifelse(off > 0, "after", "within")))
  df$period_match <- cls
  df$years_outside <- off
  lev <- c("within", "before", "after", "unknown")
  s <- data.frame(period_match = lev,
                  n = as.integer(table(factor(cls, lev))))
  s$share <- round(s$n / nrow(df), 3)
  s$median_year <- vapply(lev, function(k) {
    if (any(cls == k & !is.na(y))) stats::median(y[cls == k]) else NA_real_
  }, numeric(1))
  attr(df, "summary") <- s
  if (verbose) {
    message("Climate period ", period[1], "-", period[2],
            if (tolerance > 0) paste0(" (+/- ", tolerance, " years)"), ": ",
            paste0(s$n, " ", s$period_match, collapse = ", "),
            "; median year of the records: ",
            stats::median(y, na.rm = TRUE))
  }
  df
}

#' A block size for spatial cross-validation
#'
#' Suggests `block_size` for [spatial_folds()] from the distance over which
#' the predictors are spatially autocorrelated: blocks at least that large
#' make the held-out fold nearly independent of the training data (Roberts
#' et al. 2017; Valavi et al. 2019). For each layer the empirical
#' semivariogram is computed from a sample of cells, and its *range* is
#' taken as the shortest distance at which the semivariance reaches
#' `sill_share` (default 95%) of the layer's variance. The suggestion is
#' the median range over the layers, as in `blockCV`.
#'
#' Layers that change steadily across the whole area (a trend, such as
#' temperature from north to south) have no patch size: their range comes
#' out large, of the order of the size of the area, or is not reached
#' within the distances examined (`reached = FALSE`; the longest distance
#' is then reported). If the suggested size leaves too few blocks for the
#' folds, blocks cannot be made independent in this area, and the largest
#' size that still gives each fold enough records is the practical choice
#' (see the warnings of [spatial_folds()]).
#'
#' @param env Predictor rasters ([terra::SpatRaster] or file paths).
#' @param n Number of cells sampled (all pairs of them are compared).
#' @param max_km Longest distance examined (default: half the largest
#'   distance between sampled cells).
#' @param n_bins Number of distance classes.
#' @param sill_share Share of the variance that defines the range.
#' @param seed Seed for the sample.
#' @return A `block_size` object: a list with `block_size_km` (the median
#'   range), `ranges` (one row per layer: `range_km`, `reached`) and
#'   `variogram` (the semivariance of each layer per distance class, in
#'   units of its variance).
#' @references Roberts, D.R. et al. (2017) Cross-validation strategies for
#'   data with temporal, spatial, hierarchical, or phylogenetic structure.
#'   *Ecography* 40, 913-929. Valavi, R., Elith, J., Lahoz-Monfort, J.J.,
#'   Guillera-Arroita, G. (2019) blockCV: an R package for generating
#'   spatially or environmentally separated folds for k-fold
#'   cross-validation of species distribution models. *Methods in Ecology
#'   and Evolution* 10, 225-232.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' b <- suggest_block_size(sim$env, n = 1000)
#' b
#' # folds <- spatial_folds(points, k = 5, block_size = b$block_size_km)
suggest_block_size <- function(env, n = 2000, max_km = NULL, n_bins = 30,
                               sill_share = 0.95, seed = 1) {
  r <- as_grid(env)
  vars <- names(r)
  valid <- valid_cells(r)
  if (length(valid) < 30) stop("Fewer than 30 cells with values.",
                               call. = FALSE)
  s <- with_seed(seed, sort(valid[sample.int(length(valid),
                                             min(length(valid), n))]))
  Z <- read_cells(r, s)
  p <- terra::vect(terra::xyFromCell(r, s), crs = terra::crs(r))
  if (terra::is.lonlat(r)) {
    ll <- terra::crds(p)
    p <- terra::project(p, sprintf(paste0(
      "+proj=laea +lat_0=%.4f +lon_0=%.4f +x_0=0 +y_0=0 +datum=WGS84 ",
      "+units=m +no_defs"), mean(ll[, 2]), mean(ll[, 1])))
  }
  d <- as.vector(stats::dist(terra::crds(p))) / 1000
  if (is.null(max_km)) max_km <- max(d) / 2
  use <- d <= max_km
  breaks <- seq(0, max_km, length.out = n_bins + 1)
  bin <- findInterval(d[use], breaks, rightmost.closed = TRUE)
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  npair <- tabulate(bin, n_bins)
  vg <- do.call(rbind, lapply(seq_along(vars), function(j) {
    z <- Z[, j]
    v <- stats::var(z)
    # squared differences of all pairs, in the order of dist()
    g <- as.vector(stats::dist(z))[use]^2 / 2
    gamma <- as.vector(tapply(g, factor(bin, seq_len(n_bins)), mean))
    data.frame(layer = vars[j], distance_km = mid, pairs = npair,
               semivariance = if (v > 0) gamma / v else NA_real_)
  }))
  ranges <- do.call(rbind, lapply(vars, function(v) {
    x <- vg[vg$layer == v & vg$pairs > 0, ]
    hit <- which(x$semivariance >= sill_share)
    data.frame(layer = v,
               range_km = if (length(hit)) x$distance_km[hit[1]] else
                 if (all(is.na(x$semivariance))) NA_real_ else max_km,
               reached = length(hit) > 0)
  }))
  out <- list(block_size_km = signif(stats::median(ranges$range_km,
                                                   na.rm = TRUE), 3),
              ranges = ranges, variogram = vg,
              settings = list(n = length(s), max_km = max_km,
                              sill_share = sill_share))
  class(out) <- "block_size"
  out
}

#' @export
print.block_size <- function(x, ...) {
  cat("Suggested block size: ", x$block_size_km, " km (median range of ",
      nrow(x$ranges), " layers; ", x$settings$n, " cells, distances up to ",
      signif(x$settings$max_km, 3), " km)\n", sep = "")
  if (any(!x$ranges$reached, na.rm = TRUE)) {
    cat("  Range not reached within that distance: ",
        paste(x$ranges$layer[!x$ranges$reached], collapse = ", "), "\n",
        sep = "")
  }
  r <- x$ranges
  r$range_km <- signif(r$range_km, 3)
  print(r, row.names = FALSE)
  invisible(x)
}
