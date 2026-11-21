# How much of a projected change comes from the study design? ----------------------

# Points of any accepted kind as an x/y matrix in the CRS of `r` (internal)
points_xy <- function(points, r, coords = NULL, points_crs = "EPSG:4326") {
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  if (is.data.frame(points) && !inherits(points, "sf")) {
    xy_cols <- if (is.null(coords)) guess_coords(names(points)) else coords
    x <- as_number(points[[xy_cols[1]]])
    y <- as_number(points[[xy_cols[2]]])
    ok <- !is.na(x) & !is.na(y)
    p <- terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else {
    p <- as_points(points, coords, points_crs)
  }
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  xy <- terra::crds(p)
  colnames(xy) <- c("x", "y")
  xy
}

# A weight raster on the grid of `r`, as one value per cell (internal)
weight_values <- function(weight, r) {
  w <- as_grid(weight)[[1]]
  if (!isTRUE(terra::compareGeom(w, r, stopOnError = FALSE))) {
    if (!terra::same.crs(w, r)) w <- terra::project(w, terra::crs(r))
    w <- terra::resample(w, r)
  }
  v <- terra::values(w)[, 1]
  if (any(v < 0, na.rm = TRUE)) stop("Negative weights.", call. = FALSE)
  v
}

#' Random "occurrences" for null models, with the sampling bias of the records
#'
#' Null models of ecological niche models are fitted to points that carry no
#' information about the species: locations drawn from the study area instead
#' of the real records (Raes and ter Steege 2007). How the locations are
#' drawn matters. Real records are concentrated where people collect, and
#' models of uniformly random points are then too easy to beat: Warren et
#' al. (2026) found false positive rates of 65-87% for tests of model fit
#' when the records were collected with bias and the null points were not,
#' and about the nominal 5% when the null points followed the same bias.
#' This function draws the point sets, for any modelling software:
#'
#' * `"uniform"`: every cell of the area is equally likely (the original
#'   test).
#' * `"bias"`: cells are drawn in proportion to `weight`, an estimate of
#'   recording effort, e.g. from [gbif_effort_rasters()] (Warren et al.
#'   2026). The recommended choice when such an estimate exists.
#' * `"pattern"`: the real point pattern is turned by a random angle and
#'   moved to a random place in the area, which keeps the clustering of the
#'   records exactly (Beale et al. 2008, Nunes and Pearson 2017). For when
#'   there is no estimate of effort. Points that land outside the area are
#'   drawn again at random (in proportion to `weight` when given); the mean
#'   share of such points is reported, and a high share means that the area
#'   is too small for this method. As Warren et al. (2026) point out, part
#'   of the clustering of real records is the species following a patchy
#'   environment, so this null is conservative.
#'
#' Fit your models to each set exactly as to the real records (same area,
#' background, variables, settings), then give the predictions to
#' [transfer_bias_test()], or their evaluation to your own comparison.
#' [null_model_test()] is the quick version for models that are already
#' fitted, without refitting.
#'
#' @param area The study area (the area the background is taken from): a
#'   [terra::SpatRaster] or file, of which the cells with a value in the
#'   first layer are used.
#' @param points The real records (data.frame, CSV path, `sf` or
#'   `SpatVector`). They give the number of points (records in distinct
#'   cells of the area) and, for `"pattern"`, the pattern.
#' @param n Number of points per set, instead of `points`.
#' @param method `"uniform"`, `"bias"` or `"pattern"`, see above.
#' @param weight Raster of relative recording effort, for `"bias"`.
#' @param n_sets Number of point sets.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @param seed Seed.
#' @return A list of `n_sets` data.frames with `x` and `y` (cell centres,
#'   in the CRS of `area`; for `"pattern"` the moved coordinates), with
#'   class `null_occurrences` and the attributes `method` and, for
#'   `"pattern"`, `redrawn` (mean share of points drawn again).
#' @references Raes, N., ter Steege, H. (2007) A null-model for
#'   significance testing of presence-only species distribution models.
#'   *Ecography* 30, 727-736.
#'
#'   Warren, D.L., Casadei-Ferreira, A., Economo, E.P., Dinnage, R., Kass,
#'   J.M. (2026) Incorporating sampling bias into ENM/SDM permutation
#'   tests: new methods and a case study on neotropical ants. *Ecography*
#'   2026, e08559. \doi{10.1002/ecog.08559}
#'
#'   Beale, C.M., Lennon, J.J., Gimona, A. (2008) Opening the climate
#'   envelope reveals no macroscale associations with climate in European
#'   birds. *PNAS* 105, 14908-14912.
#'
#'   Nunes, L.A., Pearson, R.G. (2017) A null biogeographical test for
#'   assessing ecological niche evolution. *Journal of Biogeography* 44,
#'   1331-1343.
#' @seealso [transfer_bias_test()], [null_model_test()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' effort <- sim$env[[1]]            # more recording towards the east
#' sets <- null_occurrences(sim$env, sim$presence, method = "bias",
#'                          weight = effort, n_sets = 5)
#' head(sets[[1]])
#' # write.csv(sets[[1]], "null_001.csv", row.names = FALSE)
null_occurrences <- function(area, points = NULL, n = NULL,
                             method = c("uniform", "bias", "pattern"),
                             weight = NULL, n_sets = 99, coords = NULL,
                             points_crs = "EPSG:4326", seed = 1) {
  method <- match.arg(method)
  r <- as_grid(area)[[1]]
  valid <- which(!is.na(terra::values(r)[, 1]))
  if (!length(valid)) stop("`area` has no cells with values.", call. = FALSE)
  whole <- function(v, what) {
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 1) {
      stop("`", what, "` must be a number of at least 1.", call. = FALSE)
    }
    as.integer(round(v))
  }
  n_sets <- whole(n_sets, "n_sets")
  xy <- NULL
  if (!is.null(points)) {
    xy <- points_xy(points, r, coords, points_crs)
    cell <- terra::cellFromXY(r, xy)
    keep <- !is.na(cell) & cell %in% valid & !duplicated(cell)
    xy <- xy[keep, , drop = FALSE]
  }
  if (is.null(n)) {
    if (is.null(xy)) stop("Give `points` or `n`.", call. = FALSE)
    n <- nrow(xy)
  }
  n <- whole(n, "n")
  if (n > length(valid)) {
    stop("`area` has fewer cells (", length(valid), ") than points are ",
         "needed (", n, ").", call. = FALSE)
  }
  w <- rep(1, length(valid))
  if (!is.null(weight)) {
    w <- weight_values(weight, r)[valid]
    w[is.na(w)] <- 0
    if (sum(w > 0) < n) {
      stop("`weight` is positive in fewer cells than points are needed.",
           call. = FALSE)
    }
  } else if (method == "bias") {
    stop("`method = \"bias\"` needs `weight`, a raster of recording effort.",
         call. = FALSE)
  }
  draw <- function(k, not = integer()) {
    ww <- w
    ww[valid %in% not] <- 0
    if (method == "uniform" && is.null(weight)) ww[ww > 0] <- 1
    pos <- which(ww > 0)
    valid[pos[weighted_draw(ww[pos], k)]]
  }
  cells_df <- function(cells) {
    p <- terra::xyFromCell(r, cells)
    data.frame(x = p[, 1], y = p[, 2])
  }
  redrawn <- NULL
  sets <- with_seed(seed, {
    if (method != "pattern") {
      if (method == "uniform") w[] <- 1
      lapply(seq_len(n_sets), function(i) cells_df(draw(n)))
    } else {
      if (is.null(xy) || nrow(xy) < 3) {
        stop("`method = \"pattern\"` needs `points` with at least 3 records ",
             "in the area.", call. = FALSE)
      }
      if (n != nrow(xy)) {
        stop("With `method = \"pattern\"` the number of points is that of ",
             "the records; leave `n` out.", call. = FALSE)
      }
      lonlat <- terra::is.lonlat(r)
      c0 <- colMeans(xy)
      # offsets from the centre in map units of equal length
      k0 <- if (lonlat) cos(c0[2] * pi / 180) else 1
      dx <- (xy[, 1] - c0[1]) * k0
      dy <- xy[, 2] - c0[2]
      share <- numeric(n_sets)
      out <- vector("list", n_sets)
      for (i in seq_len(n_sets)) {
        best <- NULL
        for (try in 1:30) {
          a <- stats::runif(1, 0, 2 * pi)
          c1 <- terra::xyFromCell(r, valid[sample.int(length(valid), 1)])
          k1 <- if (lonlat) max(cos(c1[2] * pi / 180), 0.05) else 1
          nx <- c1[1] + (dx * cos(a) - dy * sin(a)) / k1
          ny <- c1[2] + dx * sin(a) + dy * cos(a)
          cl <- terra::cellFromXY(r, cbind(nx, ny))
          inside <- !is.na(cl) & cl %in% valid & !duplicated(cl)
          if (is.null(best) || sum(inside) > sum(best$inside)) {
            best <- list(x = nx, y = ny, cl = cl, inside = inside)
          }
          if (mean(inside) >= 0.9) break
        }
        miss <- which(!best$inside)
        if (length(miss)) {
          add <- terra::xyFromCell(r, draw(length(miss),
                                           not = best$cl[best$inside]))
          best$x[miss] <- add[, 1]
          best$y[miss] <- add[, 2]
        }
        share[i] <- length(miss) / n
        out[[i]] <- data.frame(x = best$x, y = best$y)
      }
      redrawn <- mean(share)
      out
    }
  })
  names(sets) <- sprintf("null_%03d", seq_len(n_sets))
  if (!is.null(redrawn) && redrawn > 0.25) {
    warning(round(100 * redrawn), "% of the points fell outside the area ",
            "and were drawn again at random: the area is small for ",
            "`method = \"pattern\"`, and the clustering is only partly kept.",
            call. = FALSE)
  }
  structure(sets, class = c("null_occurrences", "list"), method = method,
            redrawn = redrawn, crs = terra::crs(r))
}

#' @export
print.null_occurrences <- function(x, ...) {
  cat(length(x), " sets of ", nrow(x[[1]]), " null points (method \"",
      attr(x, "method"), "\")", sep = "")
  rd <- attr(x, "redrawn")
  if (!is.null(rd)) cat(", ", round(100 * rd), "% of the points drawn again",
                        sep = "")
  cat(".\nFit the models to each set as to the real records.\n")
  invisible(x)
}

# Summaries of change for one set of predictions (internal).
# V: cells x layers matrix (baseline first), sel: list of logical vectors
# (scopes), thr: threshold or NULL
change_metrics <- function(V, sel, thr) {
  base <- V[, 1]
  out <- list()
  for (j in seq_len(ncol(V))[-1]) {
    d <- V[, j] - base
    for (s in names(sel)) {
      k <- sel[[s]] & !is.na(d)
      m <- c(mean_change = mean(d[k]),
             share_declining = mean(d[k] < 0))
      if (!is.null(thr)) {
        a0 <- sum(base[k] >= thr)
        a1 <- sum(V[k, j] >= thr)
        m <- c(m, area_change = if (a0 > 0) (a1 - a0) / a0 else NA_real_)
      }
      out[[length(out) + 1]] <- data.frame(
        projection = colnames(V)[j], scope = s, metric = names(m),
        value = unname(m), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, out)
}

#' Is a projected change driven by the data or by the study design?
#'
#' Models projected to another time (or place) can predict a decline or an
#' increase of suitable conditions that has little to do with the species.
#' The set of environments available in the study area shifts between the
#' periods, every algorithm reacts to that shift in its own way, and the
#' result can be nearly the same whatever points the model is fitted to.
#' Warren et al. (2021) showed this by fitting the same models to random
#' points: the "effect of climate change" on a species was reproduced by
#' models that knew nothing about it.
#'
#' The test repeats their procedure for your study. Fit the same models
#' (same area, background, variables, algorithm and settings) to sets of
#' null points from [null_occurrences()], project them to the same
#' scenarios, and give the predictions here. The function measures the
#' predicted change in each null model and places the change predicted
#' from the real records within that distribution. The verdict has four
#' classes:
#'
#' * **design-driven**: the real model predicts what the random-point
#'   models predict, and those agree on the direction. The same result
#'   would have come out for any species; it says nothing about this one.
#' * **not distinguishable**: the real model is within the spread of the
#'   random-point models, which do not agree on a direction.
#' * **design-dominated**: the real model is outside 95% of the
#'   random-point models, but these agree on the direction and reproduce at
#'   least `dominance` (default half) of the observed change. The records
#'   change the prediction, but most of the predicted change comes from the
#'   design. This includes the case in which random points predict a
#'   *larger* change than the records (`design_share` > 1): the records
#'   weaken it, and the predicted change itself is not evidence about the
#'   species.
#' * **data-driven**: the real model is outside 95% of the random-point
#'   models, and these reproduce less than `dominance` of its change (or
#'   point the other way, or do not agree on a direction).
#'
#' **How large is the difference?** Whether the real model lies inside the
#' null interval is only the first question. Three numbers measure how far
#' it is and what it means:
#'
#' * `ses`, the standardised effect size, (observed - null mean) / null SD
#'   (Gotelli and McCabe 2002): how many standard deviations of the
#'   random-point models separate them from the real model. Unlike the
#'   p-value, which cannot go below 1 / (number of null models + 1), it
#'   keeps growing with the distance, and it can be compared between
#'   projections, metrics and studies. |SES| < 2 is roughly within the
#'   range of random points; the interval and the p-value are exact, the
#'   SES assumes a roughly symmetric distribution of the null values.
#' * `design_share`, (null mean - no change) / (observed - no change): the
#'   part of the predicted change that random points reproduce. 0 = none
#'   of it, 1 = all of it, 0.8 = random points alone give 80% of the
#'   decline; above 1 the records weaken the change the design produces,
#'   below 0 the design points the other way. A result can be clearly
#'   different from random points (large SES) and still be mostly design
#'   (design share near 1): "significant" is not the same as "driven by the
#'   records".
#' * `prob_exceed`: the probability that a model of the real records
#'   predicts a change beyond that of a model of random points, in the
#'   direction of the difference. Give several models of the real records
#'   in `observed` (e.g. the cross-validation runs or models of bootstrap
#'   samples of the records) to include their own variation: it is then
#'   the share of all pairs of a real and a random-point model in which the
#'   real one is further (the common-language effect size of McGraw and
#'   Wong 1992; 0.5 = the two overlap completely, 1 = every real model is
#'   beyond every random-point model). With one model it is the share of
#'   random-point models it exceeds.
#'
#' The null models are estimates of how the study design behaves, not a
#' hypothesis that must be rejected before a model may be used: a
#' design-driven prediction can still be right (if the whole area becomes
#' unsuitable, every model will say so). But a result that is the same for
#' any data needs other lines of evidence (Warren et al. 2021).
#'
#' The measures of change, for the whole predicted area and, with `region`,
#' for a part of it (e.g. the present range):
#' `mean_change` (mean of projected minus present suitability),
#' `share_declining` (share of cells where suitability falls) and, with a
#' `threshold`, `area_change` (proportional change of the area at or above
#' the threshold: -0.5 = half is lost, 1 = doubled).
#'
#' With `maps = TRUE` the result also has maps per projection: the mean and
#' the standard deviation of the change across the null models, and the
#' observed change minus the null mean. Where the null mean is far from
#' zero and the standard deviation small, the direction of the prediction
#' is fixed by the design and the data have little leverage (Warren et al.
#' 2021, Fig. 7).
#'
#' preENM does not fit models. Either give the predictions (`observed`,
#' `null`), made with any software, or give your own function `fit` that
#' fits and projects a model, and it is run for the real records and every
#' null set.
#'
#' @param observed Predictions of the model fitted to the real records: a
#'   [terra::SpatRaster] or file(s) with the present prediction and one
#'   layer per projection (other period, scenario or region on the same
#'   grid). Or a list of such rasters, one per model of the real records
#'   (e.g. cross-validation runs): their mean is the observed value, and
#'   `prob_exceed` includes their variation.
#' @param null Predictions of the null models: a list of rasters (or of
#'   file paths) with the same layers as `observed`, one element per null
#'   set.
#' @param fit Instead of `observed` and `null`: a function that takes a
#'   data.frame of points (`x`, `y`) and returns the predictions as
#'   described for `observed`.
#' @param points The real records, for `fit` and for `threshold = "p10"`.
#' @param null_points The null sets from [null_occurrences()], for `fit`
#'   and for `threshold = "p10"`.
#' @param baseline The layer with the present prediction (name or number).
#' @param region Optional part of the area for a second set of summaries:
#'   polygons (`sf`, `SpatVector`) or a raster (cells that are not `NA` and
#'   not 0).
#' @param threshold Optional threshold of suitability for `area_change`: a
#'   number, or `"p10"` for the 10th percentile of the present prediction
#'   at the points each model was fitted to (needs `points` and
#'   `null_points`).
#' @param dominance Smallest `design_share` at which a result outside the
#'   null interval is called "design-dominated" instead of "data-driven"
#'   (default 0.5: random points reproduce at least half of the change).
#' @param maps Also return the maps.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @param verbose Show progress when `fit` is run.
#' @return A list of class `transfer_bias` with `summary` (one row per
#'   projection, scope and metric: `observed`, `null_mean`, `null_sd`,
#'   `null_lo` and `null_hi` (95% interval), `null_agreement` (share of
#'   null models with the direction of the null mean), `signal`
#'   (`observed - null_mean`), `ses`, `design_share`, `prob_exceed`,
#'   `percentile` (share of null values below the observed one),
#'   `p_value` (two-sided) and `verdict`), `null` (the values of every null
#'   model), `observed` (the values of every model of the records), `maps`,
#'   `n_null` and `n_observed`.
#' @references Gotelli, N.J., McCabe, D.J. (2002) Species co-occurrence: a
#'   meta-analysis of J. M. Diamond's assembly rules model. *Ecology* 83,
#'   2091-2096.
#'
#'   McGraw, K.O., Wong, S.P. (1992) A common language effect size
#'   statistic. *Psychological Bulletin* 111, 361-365.
#'
#'   Warren, D.L., Dornburg, A., Zapfe, K., Iglesias, T.L. (2021)
#'   The effects of climate change on Australia's only endemic Pokemon:
#'   measuring bias in species distribution models. *Methods in Ecology and
#'   Evolution* 12, 985-995. \doi{10.1111/2041-210X.13591}
#'
#'   Warren, D.L., Casadei-Ferreira, A., Economo, E.P., Dinnage, R., Kass,
#'   J.M. (2026) Incorporating sampling bias into ENM/SDM permutation
#'   tests: new methods and a case study on neotropical ants. *Ecography*
#'   2026, e08559. \doi{10.1002/ecog.08559}
#' @seealso [null_occurrences()], [environment_shift()],
#'   [extrapolation_risk()], [plot_transfer_bias()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' present <- sim$env
#' future <- sim$env + 0.15              # everything shifts
#' names(future) <- names(present)
#' # your own model: here a simple climate envelope (range of the records)
#' fit <- function(p) {
#'   e <- terra::extract(present, as.matrix(p[, c("x", "y")]))
#'   lim <- apply(e, 2, range, na.rm = TRUE)
#'   inside <- function(r) {
#'     out <- (r[[1]] >= lim[1, 1]) * (r[[1]] <= lim[2, 1]) *
#'       (r[[2]] >= lim[1, 2]) * (r[[2]] <= lim[2, 2])
#'     out
#'   }
#'   out <- c(inside(present), inside(future))
#'   names(out) <- c("present", "future")
#'   out
#' }
#' sets <- null_occurrences(present, sim$presence, n_sets = 19)
#' tb <- transfer_bias_test(fit = fit, points = sim$presence,
#'                          null_points = sets, verbose = FALSE)
#' tb
transfer_bias_test <- function(observed = NULL, null = NULL, fit = NULL,
                               points = NULL, null_points = NULL,
                               baseline = 1, region = NULL, threshold = NULL,
                               dominance = 0.5, maps = TRUE, coords = NULL,
                               points_crs = "EPSG:4326", verbose = TRUE) {
  if (!is.numeric(dominance) || length(dominance) != 1 || dominance <= 0 ||
      dominance >= 1) {
    stop("`dominance` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(fit)) {
    if (!is.function(fit)) stop("`fit` must be a function.", call. = FALSE)
    if (is.null(points) || is.null(null_points)) {
      stop("With `fit`, give `points` and `null_points` (see ",
           "null_occurrences()).", call. = FALSE)
    }
  } else if (is.null(observed) || is.null(null)) {
    stop("Give `observed` and `null` predictions, or `fit` with `points` ",
         "and `null_points`.", call. = FALSE)
  }
  n_null <- if (is.null(fit)) length(null) else length(null_points)
  if (n_null < 2) stop("At least 2 null models are needed.", call. = FALSE)
  pts <- NULL
  # several models of the real records (e.g. cross-validation runs)
  obs_list <- NULL
  if (is.null(fit) && is.list(observed) && !inherits(observed, "SpatRaster")) {
    if (length(observed) < 1) stop("`observed` is empty.", call. = FALSE)
    obs_list <- observed
    observed <- observed[[1]]
  }
  get_obs <- function() {
    if (is.null(fit)) return(as_grid(observed))
    as_grid(fit(pts))
  }
  # the real records as x/y in the CRS of the null points (that of the area)
  if (!is.null(points)) {
    crs <- attr(null_points, "crs")
    if (is.null(crs) && is.null(fit)) crs <- terra::crs(as_grid(observed))
    if (is.null(crs)) crs <- points_crs
    pts <- as.data.frame(points_xy(points, terra::rast(crs = crs), coords,
                                   points_crs))
  }
  obs <- get_obs()
  if (terra::nlyr(obs) < 2) {
    stop("The predictions need the present layer and at least one ",
         "projection.", call. = FALSE)
  }
  lyr <- names(obs)
  if (anyDuplicated(lyr)) stop("Layer names must be unique.", call. = FALSE)
  b <- if (is.character(baseline)) match(baseline, lyr) else baseline
  if (is.na(b) || b < 1 || b > length(lyr)) {
    stop("`baseline` is not a layer of the predictions.", call. = FALSE)
  }
  ord <- c(b, setdiff(seq_along(lyr), b))
  read <- function(x) {
    x <- as_grid(x)
    if (!isTRUE(terra::compareGeom(x, obs, stopOnError = FALSE)) ||
        terra::nlyr(x) != length(lyr)) {
      stop("Every null prediction needs the grid and the layers of the ",
           "observed one.", call. = FALSE)
    }
    V <- terra::values(x)[, ord, drop = FALSE]
    colnames(V) <- lyr[ord]
    V
  }
  Vobs <- terra::values(obs)[, ord, drop = FALSE]
  colnames(Vobs) <- lyr[ord]
  sel <- list(area = rep(TRUE, nrow(Vobs)))
  if (!is.null(region)) {
    sel$region <- in_mask(obs, seq_len(nrow(Vobs)), region)
    if (!any(sel$region)) stop("`region` covers no cell.", call. = FALSE)
  }
  p10 <- identical(threshold, "p10")
  if (p10 && (is.null(pts) || is.null(null_points))) {
    stop("`threshold = \"p10\"` needs `points` and `null_points`.",
         call. = FALSE)
  }
  if (!is.null(threshold) && !p10 &&
      !(is.numeric(threshold) && length(threshold) == 1)) {
    stop("`threshold` must be one number or \"p10\".", call. = FALSE)
  }
  thr_of <- function(V, p) {
    if (is.null(threshold)) return(NULL)
    if (!p10) return(threshold)
    cl <- terra::cellFromXY(obs, as.matrix(as.data.frame(p)[, c("x", "y")]))
    threshold_p10(V[cl, 1])
  }
  # one column per model of the real records
  if (is.null(obs_list)) {
    om <- change_metrics(Vobs, sel, thr_of(Vobs, pts))
    om_reps <- matrix(om$value, ncol = 1)
  } else {
    Vs <- lapply(obs_list, read)
    reps <- lapply(Vs, function(V) change_metrics(V, sel, thr_of(V, pts)))
    om <- reps[[1]]
    om_reps <- vapply(reps, function(m) m$value, numeric(nrow(om)))
    om_reps <- matrix(om_reps, nrow = nrow(om))
    om$value <- rowMeans(om_reps, na.rm = TRUE)
    Vobs <- Reduce(`+`, Vs) / length(Vs)
  }
  nproj <- length(lyr) - 1
  S <- S2 <- matrix(0, nrow(Vobs), nproj)
  N <- matrix(0L, nrow(Vobs), nproj)
  nm <- matrix(NA_real_, nrow(om), n_null)
  for (i in seq_len(n_null)) {
    if (verbose && !is.null(fit)) {
      message("Null model ", i, " of ", n_null)
    }
    V <- if (is.null(fit)) read(null[[i]]) else
      read(fit(as.data.frame(null_points[[i]])))
    p <- if (p10) null_points[[i]] else NULL
    nm[, i] <- change_metrics(V, sel, thr_of(V, p))$value
    if (maps) {
      d <- V[, -1, drop = FALSE] - V[, 1]
      ok <- !is.na(d)
      d[!ok] <- 0
      S <- S + d
      S2 <- S2 + d^2
      N <- N + ok
    }
  }
  q <- function(p) apply(nm, 1, stats::quantile, probs = p, na.rm = TRUE,
                         names = FALSE)
  mu <- rowMeans(nm, na.rm = TRUE)
  # direction is counted against "no change"; share_declining is centred on
  # one half
  centre <- ifelse(om$metric == "share_declining", 0.5, 0)
  agree <- vapply(seq_len(nrow(nm)), function(i) {
    v <- nm[i, ] - centre[i]
    v <- v[!is.na(v)]
    if (!length(v)) return(NA_real_)
    max(mean(v > 0), mean(v < 0))
  }, numeric(1))
  nn <- rowSums(!is.na(nm))
  pval <- (1 + rowSums(abs(nm - mu) >= abs(om$value - mu), na.rm = TRUE)) /
    (nn + 1)
  lo <- q(0.025)
  hi <- q(0.975)
  inside <- !is.na(om$value) & om$value >= lo & om$value <= hi
  directed <- (lo - centre) * (hi - centre) > 0
  sdn <- apply(nm, 1, stats::sd, na.rm = TRUE)
  signal <- om$value - mu
  # standardised effect size: distance from the null mean in null SDs
  ses <- ifelse(sdn > 0, signal / sdn,
                ifelse(abs(signal) > 0, sign(signal) * Inf, 0))
  # share of the observed change (from "no change") that random points
  # reproduce: 1 = all of it, 0 = none, > 1 = the records weaken it,
  # < 0 = the design points the other way
  oc <- om$value - centre
  design_share <- ifelse(abs(oc) > 1e-12, (mu - centre) / oc, NA_real_)
  # probability that a model of the real records predicts a change beyond
  # that of a model of random points, in the direction of the signal
  # (pairs of models; with one model of the records, the share of null
  # models it exceeds)
  prob_exceed <- vapply(seq_len(nrow(nm)), function(i) {
    a <- om_reps[i, ]
    b <- nm[i, ]
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    if (!length(a) || !length(b) || is.na(signal[i])) return(NA_real_)
    gt <- mean(outer(a, b, ">")) + 0.5 * mean(outer(a, b, "=="))
    if (signal[i] >= 0) gt else 1 - gt
  }, numeric(1))
  dominated <- !is.na(design_share) & design_share >= dominance &
    directed
  summary <- data.frame(
    om[c("projection", "scope", "metric")],
    observed = om$value, null_mean = mu, null_sd = sdn,
    null_lo = lo, null_hi = hi, null_agreement = agree,
    signal = signal, ses = ses, design_share = design_share,
    prob_exceed = prob_exceed,
    percentile = rowMeans(nm < om$value, na.rm = TRUE),
    p_value = pval,
    verdict = ifelse(is.na(om$value), NA_character_,
                     ifelse(inside,
                            ifelse(directed, "design-driven",
                                   "not distinguishable"),
                            ifelse(dominated, "design-dominated",
                                   "data-driven"))),
    stringsAsFactors = FALSE, row.names = NULL)
  colnames(nm) <- sprintf("null_%03d", seq_len(n_null))
  colnames(om_reps) <- sprintf("observed_%03d", seq_len(ncol(om_reps)))
  out <- list(summary = summary,
              null = cbind(om[c("projection", "scope", "metric")], nm),
              observed = cbind(om[c("projection", "scope", "metric")],
                               om_reps),
              maps = NULL, n_null = n_null, n_observed = ncol(om_reps),
              baseline = lyr[b], threshold = threshold,
              dominance = dominance)
  if (maps) {
    mean_map <- S / pmax(N, 1)
    sd_map <- sqrt(pmax(S2 / pmax(N, 1) - mean_map^2, 0) *
                     N / pmax(N - 1, 1))
    mean_map[N == 0] <- NA
    sd_map[N < 2] <- NA
    dobs <- Vobs[, -1, drop = FALSE] - Vobs[, 1]
    out$maps <- lapply(seq_len(nproj), function(j) {
      r <- terra::rast(obs[[1]], nlyrs = 4)
      terra::values(r) <- cbind(dobs[, j], mean_map[, j], sd_map[, j],
                                dobs[, j] - mean_map[, j])
      names(r) <- c("observed_change", "null_mean", "null_sd", "signal")
      r
    })
    names(out$maps) <- lyr[ord][-1]
  }
  structure(out, class = "transfer_bias")
}

#' @export
print.transfer_bias <- function(x, digits = 3, ...) {
  s <- x$summary
  nobs <- if (is.null(x$n_observed)) 1 else x$n_observed
  cat("Design bias of projected change (", x$n_null, " null models, ",
      if (nobs > 1) paste0(nobs, " models of the records, "),
      "present = ", x$baseline, ")\n\n", sep = "")
  show <- s[intersect(c("projection", "scope", "metric", "observed",
                        "null_mean", "null_lo", "null_hi", "ses",
                        "design_share", "p_value", "verdict"), names(s))]
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], signif, digits)
  if (length(unique(show$scope)) == 1) show$scope <- NULL
  print(show, row.names = FALSE)
  main <- s[s$metric == "mean_change" & !is.na(s$verdict), , drop = FALSE]
  cat("\nIN SHORT (mean change of suitability)\n")
  for (i in seq_len(nrow(main))) {
    r <- main[i, ]
    where <- if (length(unique(s$scope)) > 1) paste0(" [", r$scope, "]") else ""
    cat(strwrap(paste0(r$projection, where, ": ", transfer_text(r, digits)),
                width = 78, exdent = 2, prefix = "  ", initial = "- "),
        sep = "\n")
  }
  cat("\nSES: distance from the random-point models in their SDs. ",
      "design share: part of\nthe change that random points reproduce ",
      "(1 = all of it).\n", sep = "")
  if (x$n_null < 99) {
    cat("Fewer than 99 null models: the intervals are rough.\n")
  }
  invisible(x)
}

# One projection in words (internal)
transfer_text <- function(r, digits = 3) {
  dir <- function(v) if (v < 0) "a decline" else "an increase"
  f <- function(v) signif(v, digits)
  ses <- if (is.finite(r$ses)) sprintf("%.1f", r$ses) else
    if (r$ses > 0) "+Inf" else "-Inf"
  pct <- function(v) paste0(round(100 * v), "%")
  ds <- r$design_share
  switch(
    r$verdict,
    "design-driven" = paste0(
      "models of random points predict ", dir(r$null_mean), " too (",
      round(100 * r$null_agreement), "% agree), and the real model is ",
      "within their range (SES = ", ses, "). This result does not ",
      "depend on the records."),
    "not distinguishable" = paste0(
      "the real model is within the range of models of random points ",
      "(SES = ", ses, "), which do not agree on a direction. The records ",
      "add little."),
    "design-dominated" = if (ds > 1) {
      paste0(
        "random points alone predict ", dir(r$null_mean), " of ",
        f(r$null_mean), ", larger than that of the real model (",
        f(r$observed), "; SES = ", ses, "). The records weaken the ",
        "change; the predicted ", sub("^an? ", "", dir(r$observed)),
        " itself comes from the design, not from this species.")
    } else {
      paste0(
        "the real model (", f(r$observed), ") differs from models of ",
        "random points (SES = ", ses, "), but random points alone give ",
        pct(ds), " of this change (", f(r$null_mean), "). Most of the ",
        "predicted change is a product of the design; the records add ",
        sprintf("%+.3g", r$signal), ".")
    },
    "data-driven" = paste0(
      "the real model (", f(r$observed), ") differs from models of random ",
      "points (SES = ", ses, ", p = ", f(r$p_value), "). ",
      if (is.na(ds)) "" else if (ds < 0) {
        paste0("Random points predict the opposite direction (",
               f(r$null_mean), "), so the change comes from the records.")
      } else {
        paste0("Random points alone give ", pct(ds), " of it (",
               f(r$null_mean), "); the records add ",
               sprintf("%+.3g", r$signal), ".")
      }),
    "")
}

#' Plot the design bias of a projected change
#'
#' One row per projection: the distribution of the change predicted by
#' the null models (grey; the line marks their 95% interval), "no change"
#' (dashed) and the change predicted from the real records (the point,
#' coloured by the verdict), labelled with the standardised effect size
#' (SES) and the design share (see [transfer_bias_test()]).
#'
#' @param x Result of [transfer_bias_test()].
#' @param metric `"mean_change"`, `"share_declining"` or `"area_change"`.
#' @param scope `"area"` or `"region"`.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @export
#' @examples
#' # see transfer_bias_test() for how `tb` is made
#' \dontrun{
#' plot_transfer_bias(tb)
#' }
plot_transfer_bias <- function(x, metric = "mean_change", scope = "area",
                               base_size = 11,
                               file = NULL, width = 8.5, height = 4.8) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
  if (!inherits(x, "transfer_bias")) {
    stop("`x` must come from transfer_bias_test().", call. = FALSE)
  }
  s <- x$summary[x$summary$metric == metric & x$summary$scope == scope, ,
                 drop = FALSE]
  if (!nrow(s)) {
    stop("No results for metric \"", metric, "\" and scope \"", scope,
         "\".", call. = FALSE)
  }
  nl <- x$null[x$null$metric == metric & x$null$scope == scope, ,
               drop = FALSE]
  long <- do.call(rbind, lapply(seq_len(nrow(nl)), function(i) {
    data.frame(projection = nl$projection[i],
               value = unlist(nl[i, -(1:3)], use.names = FALSE))
  }))
  long <- long[!is.na(long$value), , drop = FALSE]
  lev <- unique(s$projection)
  long$projection <- factor(long$projection, lev)
  s$projection <- factor(s$projection, lev)
  verdicts <- c("data-driven", "design-dominated", "not distinguishable",
                "design-driven")
  s$verdict <- factor(s$verdict, verdicts)
  fmt_ses <- function(v) ifelse(is.finite(v), sprintf("%.1f", v),
                                ifelse(v > 0, "+Inf", "-Inf"))
  ds <- s$design_share
  s$label <- paste0("SES ", fmt_ses(s$ses),
                    ifelse(is.na(ds), "",
                    ifelse(ds < 0, "  |  design: opposite direction",
                    ifelse(ds > 1, sprintf("  |  design alone %.1fx", ds),
                           paste0("  |  design ", round(100 * ds), "%")))))
  ink <- "#1f2933"
  muted <- "#6b7280"
  zero <- if (metric == "share_declining") 0.5 else 0
  lab <- c(mean_change = "Mean change of suitability",
           share_declining = "Share of cells with falling suitability",
           area_change = "Proportional change of suitable area")[[metric]]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$projection)) +
    ggplot2::geom_vline(xintercept = zero, linetype = "dashed",
                        colour = muted) +
    ggplot2::geom_violin(fill = "#d1d5db", colour = NA, scale = "width",
                         width = 0.7) +
    ggplot2::geom_segment(
      data = s, inherit.aes = FALSE, colour = ink, linewidth = 0.7,
      ggplot2::aes(x = .data$null_lo, xend = .data$null_hi,
                   y = .data$projection, yend = .data$projection)) +
    ggplot2::geom_point(
      data = s, inherit.aes = FALSE, size = 3.2,
      ggplot2::aes(x = .data$observed, y = .data$projection,
                   colour = .data$verdict)) +
    ggplot2::geom_text(
      data = s, inherit.aes = FALSE, colour = ink, size = base_size / 3.6,
      vjust = -1.4,
      ggplot2::aes(x = .data$observed, y = .data$projection,
                   label = .data$label)) +
    ggplot2::scale_colour_manual(
      values = c(`data-driven` = "#1d4ed8", `design-dominated` = "#7c3aed",
                 `not distinguishable` = "#d97706",
                 `design-driven` = "#b91c1c"),
      labels = c(`data-driven` = "data-driven",
                 `design-dominated` = "differs, but mostly design",
                 `not distinguishable` = "not distinguishable",
                 `design-driven` = "design-driven"),
      name = "Real records:", drop = TRUE) +
    ggplot2::labs(
      x = lab, y = NULL,
      title = "Projected change: real records against random points",
      subtitle = paste0("Grey: ", x$n_null, " models of random points, ",
                        "line = their 95% interval; dashed = no change.\n",
                        "SES = distance from them in their SDs.\n",
                        "design = share of the change they reproduce ",
                        "(above 1x: more change than the records).")) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom")
}

#' How much does the available environment change between two periods?
#'
#' A check to run before any model is projected. Between the present and a
#' projection (another period, scenario or region) some conditions of the
#' study area disappear and others appear that are not there today
#' (Williams et al. 2007). This shift of the *available* environment, not
#' the species, is what makes projections biased: a model can only call
#' suitable what it saw during training, so when present conditions
#' disappear most methods predict a decline for any set of points, and
#' where new conditions appear the prediction is extrapolated and depends
#' on the algorithm (Warren et al. 2021).
#'
#' For every projection and variable the function reports the shift of the
#' mean in units of the present standard deviation (`shift_sd`), the share
#' of present cells whose value no longer occurs in the projection
#' (`disappearing`) and the share of projected cells whose value does not
#' occur at present (`novel`); the row `(any variable)` counts cells for
#' which this is true for at least one variable. "Occurs" means within the
#' range of the other period after dropping `trim` of its cells at each
#' end. A large `disappearing` share says that a decline will be predicted
#' whatever the records are; test it with [transfer_bias_test()] and map
#' the new conditions with [extrapolation_risk()].
#'
#' @param reference Present variables of the study area
#'   ([terra::SpatRaster] or files).
#' @param projections The same variables for the projection: a raster, or
#'   a named list of rasters for several projections. Layer names must
#'   match those of `reference`.
#' @param trim Share of cells dropped at each end of a range.
#' @param n Maximum number of cells sampled from each raster.
#' @param seed Seed.
#' @return A data.frame of class `environment_shift` with `projection`,
#'   `variable`, `shift_sd`, `disappearing` and `novel`.
#' @references Williams, J.W., Jackson, S.T., Kutzbach, J.E. (2007)
#'   Projected distributions of novel and disappearing climates by 2100
#'   AD. *PNAS* 104, 5738-5742.
#'
#'   Warren, D.L., Dornburg, A., Zapfe, K., Iglesias, T.L. (2021) The
#'   effects of climate change on Australia's only endemic Pokemon:
#'   measuring bias in species distribution models. *Methods in Ecology and
#'   Evolution* 12, 985-995.
#' @seealso [transfer_bias_test()], [extrapolation_risk()],
#'   [niche_truncation()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' warmer <- sim$env
#' warmer[[1]] <- warmer[[1]] + 0.3
#' environment_shift(sim$env, list(warmer = warmer))
environment_shift <- function(reference, projections, trim = 0.005,
                              n = 50000, seed = 1) {
  ref <- as_grid(reference)
  vars <- names(ref)
  if (!is.list(projections)) projections <- list(projection = projections)
  if (is.null(names(projections)) || any(names(projections) == "")) {
    names(projections) <- paste0("projection_", seq_along(projections))
  }
  if (trim < 0 || trim >= 0.25) {
    stop("`trim` must be between 0 and 0.25.", call. = FALSE)
  }
  smp <- function(x) {
    valid <- valid_cells(x)
    if (!length(valid)) {
      stop("A raster has no cells with values in all layers.", call. = FALSE)
    }
    s <- with_seed(seed, sort(valid[sample.int(length(valid),
                                               min(length(valid), n))]))
    M <- read_cells(x, s)
    colnames(M) <- names(x)
    M
  }
  R <- smp(ref)
  rng <- function(M) apply(M, 2, stats::quantile, probs = c(trim, 1 - trim),
                           names = FALSE)
  outside <- function(M, lim) {
    t(t(M) < lim[1, ]) | t(t(M) > lim[2, ])
  }
  rr <- rng(R)
  out <- lapply(names(projections), function(nm) {
    p <- as_grid(projections[[nm]])
    miss <- setdiff(vars, names(p))
    if (length(miss)) {
      stop("Projection \"", nm, "\" lacks layer(s): ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    P <- smp(p[[vars]])
    lost <- outside(R, rng(P))
    new <- outside(P, rr)
    s <- apply(R, 2, stats::sd)
    data.frame(projection = nm, variable = c(vars, "(any variable)"),
               shift_sd = c((colMeans(P) - colMeans(R)) /
                              ifelse(s > 0, s, NA), NA),
               disappearing = c(colMeans(lost), mean(rowSums(lost) > 0)),
               novel = c(colMeans(new), mean(rowSums(new) > 0)),
               stringsAsFactors = FALSE, row.names = NULL)
  })
  structure(do.call(rbind, out),
            class = c("environment_shift", "data.frame"))
}

#' @export
print.environment_shift <- function(x, digits = 2, ...) {
  cat("Change of the available environment\n\n")
  show <- as.data.frame(x)
  show$shift_sd <- round(show$shift_sd, digits)
  show$disappearing <- round(show$disappearing, 3)
  show$novel <- round(show$novel, 3)
  print(show, row.names = FALSE)
  any_row <- x[x$variable == "(any variable)", , drop = FALSE]
  cat("\n")
  for (i in seq_len(nrow(any_row))) {
    r <- any_row[i, ]
    txt <- if (r$disappearing < 0.05 && r$novel < 0.05) {
      "little change of the available conditions; little design bias expected."
    } else {
      paste0(
        round(100 * r$disappearing), "% of present cells have conditions ",
        "that disappear and ", round(100 * r$novel), "% of projected cells ",
        "have new ones. ",
        if (r$disappearing >= 0.05) paste0(
          "Expect a predicted decline that does not depend on the records",
          if (r$novel >= 0.05) "; " else "."),
        if (r$novel >= 0.05) paste0(
          "predictions in the new conditions are extrapolated and depend ",
          "on the algorithm."),
        " Check with transfer_bias_test().")
    }
    cat(strwrap(paste0(r$projection, ": ", txt), width = 78, exdent = 2,
                prefix = "  ", initial = "- "), sep = "\n")
  }
  invisible(x)
}
