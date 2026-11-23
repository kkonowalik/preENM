# Spatial cross-validation ------------------------------------------------------

#' Assign points to spatial cross-validation folds
#'
#' Random cross-validation of species distribution models is over-optimistic
#' because training and test points lie next to each other and share the same
#' environment (spatial autocorrelation). Spatial blocking splits the study
#' area into square blocks and keeps each block entirely in one fold, so
#' every model is tested on areas it has not seen.
#'
#' Points are projected to an equal-area projection centred on the data (when
#' they are in longitude/latitude), cut into square blocks of `block_size` km,
#' and the blocks are assigned to `k` folds at random. Out of `n_iter` random
#' assignments the one with the most even numbers of presences and of
#' absences/background points per fold is kept.
#'
#' **Choosing `block_size`.** Blocks should be larger than the range of
#' spatial autocorrelation of the predictors, which can be estimated with
#' `blockCV::cv_spatial_autocor()`. As a rule of thumb, values between 50 and
#' 500 km work for regional to continental studies. With few presences, use
#' fewer folds (e.g. `k = 3`) or smaller blocks so that every fold still
#' contains presences.
#'
#' @param points Points to split: a data.frame with coordinates, a CSV path,
#'   an `sf` object or a [terra::SpatVector].
#' @param k Number of folds.
#' @param block_size Block side length in kilometres.
#' @param occ Optional 0/1 vector (1 = presence) of the same length as
#'   `points`. When given, presences and absences are balanced separately.
#' @param n_iter Number of random block-to-fold assignments tried.
#' @param seed Optional random seed.
#' @inheritParams extract_predictions
#' @return An integer vector of fold numbers (1 to `k`), one per point, with
#'   attributes `block` (block id of each point) and `block_size_km`.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
#' pts$fold <- spatial_folds(pts, k = 5, block_size = 100, occ = pts$occ,
#'                           seed = 1)
#' table(pts$fold, pts$occ)
spatial_folds <- function(points, k = 5, block_size, occ = NULL,
                          n_iter = 100, seed = NULL, coords = NULL,
                          crs = "EPSG:4326") {
  if (missing(block_size) || !is.numeric(block_size) || block_size <= 0) {
    stop("Give `block_size` in km (e.g. 100). See ?spatial_folds for how ",
         "to choose it.", call. = FALSE)
  }
  pts <- as_points(points, coords, crs)
  n <- nrow(pts)
  if (is.null(occ)) occ <- rep(1L, n)
  if (length(occ) != n) {
    stop("`occ` must have one value per point (", n, ").", call. = FALSE)
  }
  if (terra::is.lonlat(pts)) {
    xy <- terra::crds(pts)
    laea <- sprintf(paste0("+proj=laea +lat_0=%.4f +lon_0=%.4f +x_0=0 ",
                           "+y_0=0 +datum=WGS84 +units=m +no_defs"),
                    mean(xy[, 2]), mean(xy[, 1]))
    pts <- terra::project(pts, laea)
  }
  xy <- terra::crds(pts)
  size <- block_size * 1000
  block <- paste(floor(xy[, 1] / size), floor(xy[, 2] / size))
  ub <- unique(block)
  bi <- match(block, ub)
  if (length(ub) < k) {
    stop("Only ", length(ub), " blocks contain points, fewer than k = ", k,
         ". Use a smaller `block_size` or fewer folds.", call. = FALSE)
  }
  np <- tabulate(bi[occ == 1], length(ub))
  na <- tabulate(bi[occ != 1], length(ub))
  cv <- function(v) if (mean(v) == 0) 0 else stats::sd(v) / mean(v)

  best <- NULL
  best_score <- Inf
  with_seed(seed, {
    for (i in seq_len(max(1, n_iter))) {
      f <- sample(rep_len(seq_len(k), length(ub)))
      sp <- tabulate(rep(f, np), k)
      sa <- tabulate(rep(f, na), k)
      # a fold without presences cannot be evaluated: heavy penalty
      score <- cv(sp) + cv(sa) + 100 * sum(sp == 0 & sum(np) >= k)
      if (score < best_score) {
        best <- f
        best_score <- score
      }
    }
  })
  fold <- best[bi]
  per_fold <- tabulate(fold[occ == 1], k)
  if (any(per_fold == 0) && any(occ != 1)) {
    warning("Fold(s) ", paste(which(per_fold == 0), collapse = ", "),
            " contain no presences; AUC cannot be computed for them. Use ",
            "fewer folds or smaller blocks.", call. = FALSE)
  } else if (any(per_fold < 5) && any(occ != 1)) {
    warning("Some folds have fewer than 5 presences (", paste(per_fold,
            collapse = ", "), "); their metrics will be noisy.",
            call. = FALSE)
  }
  attr(fold, "block") <- bi
  attr(fold, "block_size_km") <- block_size
  fold
}

#' Evaluate cross-validated models on their held-out folds
#'
#' Summarises spatial (or any) cross-validation results per model or
#' algorithm. You supply, for every cross-validation run, the predictions of
#' the model *on the points it was not trained on*; the function computes
#' AUC, MAE, bias, TSS and the Boyce index for each run and averages them.
#'
#' Runs can also be repeated over several background (pseudo-absence) sets via
#' `set`. The output then shows whether the background sets or the spatial
#' folds cause more variation (`AUC_sd_sets` vs `AUC_sd_folds`): if the sets
#' vary much less than the folds, adding more background sets gains little.
#'
#' @param data A data.frame in long format with one row per held-out point and
#'   model run, containing the columns named by `model`, `fold`, `obs` and
#'   `pred` (and `set`, if used). Predictions from biomod2, ENMeval, Maxent,
#'   tidymodels or your own code can all be brought into this form (see the
#'   README).
#' @param model Column with the model or algorithm name (results are
#'   summarised per value of this column).
#' @param fold Column with the fold that was held out.
#' @param set Optional column with the background set / replicate.
#' @param obs Column with observations (1 = presence, 0 = absence or
#'   background).
#' @param pred Column with predicted suitability (0-1).
#' @param threshold Optional column with a fixed threshold per run (e.g. the
#'   TSS-optimal or 10th-percentile threshold found on that run's
#'   *calibration* data, as biomod2 does for its validation TSS). Without it,
#'   TSS uses the best threshold on the held-out data, which is slightly
#'   optimistic.
#' @param metrics Optional extra metrics: a named list of functions
#'   `function(obs, pred)` and/or names of optional metrics
#'   (`"calibration_slope"`, `"calibration_intercept"`, `"brier"`,
#'   `"Boyce_smooth"`),
#'   computed per run and summarised like the built-in ones (see
#'   [evaluate_models()]).
#' @param coords Optional names of the coordinate columns of `data` (e.g.
#'   `c("lon", "lat")`): the folds are then audited with [fold_report()]
#'   (within each `set`, on the rows of the first model), the audit is
#'   printed with the metrics and a warning is given when leakage is found.
#' @param crs Coordinate reference system of `coords`.
#' @param env,block_size Optional predictor rasters and block size (km) for
#'   the audit; see [fold_report()].
#' @param audit Alternatively, a ready [fold_report()].
#' @return An object of class `enm_cv`: a data.frame with one row per model,
#'   the mean of each metric over runs, their SDs (`AUC_sd`, ...), `n_runs`
#'   and, when available, `AUC_sd_folds` and `AUC_sd_sets`. Per-run metrics
#'   are in `attr(x, "runs")` and mean AUC per fold in `attr(x, "by_fold")`
#'   (a fold that is low for every model points to a region that is hard to
#'   predict, not to a bad model). With `coords` or `audit`, the
#'   [fold_report()] is in `attr(x, "audit")`.
#' @export
#' @examples
#' # see the README for a full example with spatial folds and two GLMs
#' set.seed(1)
#' d <- data.frame(model = rep(c("good", "poor"), each = 200),
#'                 fold = rep(1:4, 100), obs = rep(rbinom(200, 1, 0.3), 2))
#' d$pred <- ifelse(d$model == "good", plogis(3 * d$obs - 1.5 +
#'                  rnorm(400)), runif(400))
#' evaluate_cv(d)
evaluate_cv <- function(data, model = "model", fold = "fold", set = NULL,
                        obs = "obs", pred = "pred", threshold = NULL,
                        metrics = NULL, coords = NULL, crs = "EPSG:4326",
                        env = NULL, block_size = NULL, audit = NULL) {
  extra <- check_extra_metrics(metrics)
  cols <- c(model, fold, set, obs, pred, threshold)
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop("Columns not found in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  if (!all(data[[obs]] %in% c(0, 1))) {
    stop("`", obs, "` must contain only 0 and 1.", call. = FALSE)
  }
  check_prediction_scale(data[pred])
  set_v <- if (is.null(set)) rep(1L, nrow(data)) else data[[set]]
  key <- interaction(data[[model]], set_v, data[[fold]], drop = TRUE,
                     lex.order = TRUE, sep = "\r")
  idx <- split(seq_len(nrow(data)), key)
  runs <- do.call(rbind, lapply(idx, function(i) {
    m <- metric_row(data[[obs]][i], data[[pred]][i], extra)
    if (!is.null(threshold)) {
      thr <- unique(data[[threshold]][i])
      if (length(thr) != 1) {
        stop("`", threshold, "` must be constant within each run.",
             call. = FALSE)
      }
      ok <- !is.na(data[[pred]][i])
      m[["TSS"]] <- tss_score(data[[obs]][i][ok], data[[pred]][i][ok],
                              threshold = thr)[1]
    }
    data.frame(model = as.character(data[[model]][i[1]]), set = set_v[i[1]],
               fold = data[[fold]][i[1]], t(m), check.names = FALSE)
  }))
  rownames(runs) <- NULL
  runs$n_presence <- as.integer(runs$n_presence)
  runs$n_absence <- as.integer(runs$n_absence)
  bad <- runs$n_presence == 0 | runs$n_absence == 0
  if (any(bad)) {
    warning(sum(bad), " run(s) have no presences or no absences in the ",
            "held-out fold; their AUC/TSS/Boyce are NA and are ignored in ",
            "the means.", call. = FALSE)
  }

  mets <- names(metric_funs(extra))
  out <- do.call(rbind, lapply(split(runs, factor(runs$model,
                                                  unique(runs$model))),
                               function(r) {
    mn <- vapply(r[mets], mean, numeric(1), na.rm = TRUE)
    sdv <- vapply(r[mets], stats::sd, numeric(1), na.rm = TRUE)
    d <- data.frame(model = r$model[1], n_runs = sum(!is.na(r$AUC)),
                    t(mn), t(stats::setNames(sdv, paste0(mets, "_sd"))),
                    check.names = FALSE)
    by_set <- split(r$AUC, r$set)
    d$AUC_sd_folds <- mean(vapply(by_set, stats::sd, numeric(1),
                                  na.rm = TRUE), na.rm = TRUE)
    if (!is.null(set)) {
      d$AUC_sd_sets <- stats::sd(vapply(by_set, mean, numeric(1),
                                        na.rm = TRUE))
    }
    d
  }))
  rownames(out) <- NULL
  out[] <- lapply(out, function(v) if (is.numeric(v)) {
    v[is.nan(v)] <- NA
    v
  } else v)

  by_fold <- stats::aggregate(runs$AUC, list(model = runs$model,
                                             fold = runs$fold),
                              mean, na.rm = TRUE)
  by_fold <- stats::reshape(by_fold, idvar = "model", timevar = "fold",
                            direction = "wide")
  names(by_fold) <- sub("^x\\.", "fold_", names(by_fold))
  rownames(by_fold) <- NULL

  if (is.null(audit) && !is.null(coords)) {
    miss <- setdiff(coords, names(data))
    if (length(miss)) {
      stop("Coordinate columns not found in `data`: ",
           paste(miss, collapse = ", "), call. = FALSE)
    }
    first <- data[[model]] == data[[model]][1]
    d1 <- data[first, , drop = FALSE]
    audit <- fold_report(d1[coords], d1[[fold]], occ = d1[[obs]],
                         block_size = block_size, env = env,
                         set = if (is.null(set)) NULL else d1[[set]],
                         coords = coords, crs = crs)
  }
  if (!is.null(audit)) {
    if (!inherits(audit, "fold_report")) {
      stop("`audit` must come from fold_report().", call. = FALSE)
    }
    if (audit$leakage) {
      lk <- audit$leaks
      warning("Leakage between cross-validation folds: ", nrow(lk),
              " records in folds ", paste(sort(unique(lk$fold)),
                                         collapse = ", "),
              " have the same coordinates or raster cell as records in ",
              "other folds, so their metrics are inflated. See ",
              "attr(x, \"audit\").", call. = FALSE)
    }
  }
  attr(out, "audit") <- audit
  attr(out, "runs") <- runs
  attr(out, "by_fold") <- by_fold
  attr(out, "has_sets") <- !is.null(set)
  class(out) <- c("enm_cv", "data.frame")
  out
}

#' @export
print.enm_cv <- function(x, digits = 3, ...) {
  if (is.null(attr(x, "runs"))) return(print_plain(x, digits, ...))
  runs <- attr(x, "runs")
  cat("Cross-validated evaluation: ", nrow(x), " model(s), ",
      length(unique(runs$fold)), " folds",
      if (isTRUE(attr(x, "has_sets"))) {
        paste0(", ", length(unique(runs$set)), " background sets")
      }, "\n\n", sep = "")
  df <- as.data.frame(x)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  bf <- attr(x, "by_fold")
  if (!is.null(bf)) {
    cat("\nMean hold-out AUC per fold:\n")
    bf[-1] <- lapply(bf[-1], round, digits = digits)
    print(bf, row.names = FALSE)
  }
  if (isTRUE(attr(x, "has_sets"))) {
    cat("\nMedian SD of AUC between background sets: ",
        round(stats::median(x$AUC_sd_sets, na.rm = TRUE), 4),
        "; between folds: ",
        round(stats::median(x$AUC_sd_folds, na.rm = TRUE), 4), "\n", sep = "")
  }
  if (!is.null(attr(x, "audit"))) {
    cat("\n")
    print(attr(x, "audit"))
  }
  invisible(x)
}

# Evaluate `expr` with a temporary random seed, restoring the user's RNG state
with_seed <- function(seed, expr) {
  if (is.null(seed)) return(invisible(force(expr)))
  env <- globalenv()
  had <- exists(".Random.seed", envir = env, inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = env)
  on.exit(if (had) assign(".Random.seed", old, envir = env)
          else if (exists(".Random.seed", envir = env, inherits = FALSE))
            rm(".Random.seed", envir = env))
  set.seed(seed)
  invisible(force(expr))
}

# Print a (subsetted) result as a plain rounded data.frame (internal)
print_plain <- function(x, digits = 3, ...) {
  df <- as.data.frame(unclass(x), check.names = FALSE)
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], round, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}

# Auditing cross-validation folds ----------------------------------------------------

# Coordinates for nearest-neighbour distances in km (internal): points on a
# sphere for longitude/latitude (chord lengths, converted to great-circle
# distances by nn_km()), metres to km for projected data.
dist_space <- function(pts) {
  xy <- terra::crds(pts)
  if (terra::is.lonlat(pts)) {
    lon <- xy[, 1] * pi / 180
    lat <- xy[, 2] * pi / 180
    list(X = cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)) *
           6371.0088, lonlat = TRUE)
  } else {
    u <- terra::linearUnits(pts)
    list(X = xy * (if (is.finite(u) && u > 0) u else 1) / 1000,
         lonlat = FALSE)
  }
}

# Distance (km) from each row of Q to its nearest row of X (internal)
nn_km <- function(Q, X, lonlat) {
  Q <- as.matrix(Q)
  X <- as.matrix(X)
  if (!nrow(Q) || !nrow(X)) return(list(dist = rep(NA_real_, nrow(Q)),
                                        id = rep(NA_integer_, nrow(Q))))
  nn <- dbscan::kNN(X, k = 1, query = Q)
  d <- as.vector(nn$dist)
  if (lonlat) d <- 2 * 6371.0088 * asin(pmin(1, d / (2 * 6371.0088)))
  list(dist = d, id = as.vector(nn$id))
}

#' Audit a cross-validation fold assignment for leakage
#'
#' Cross-validation is meant to test a model on data it has not seen. Data
#' leakage, when information about the test data reaches the model during
#' training, inflates the reported performance and is widespread in
#' ecological machine learning (Stock et al. 2023). `fold_report()` checks a
#' fold assignment, for each fold held out for testing against the other
#' folds used for training:
#'
#' * **Leakage**: records with the *same coordinates* in more than one fold,
#'   and records in the *same raster cell* split between folds (the model
#'   sees identical predictor values in training and testing). Either means
#'   the metrics of the affected folds are inflated.
#' * **Spatial separation**: the distance from each test record to the
#'   nearest training record, compared with the block size and with the
#'   typical distance between neighbouring records. Folds whose test records
#'   are as close to training records as records are to each other are not
#'   spatially separated (as in random cross-validation), so spatial
#'   autocorrelation makes their metrics optimistic (Yoo and Koo 2024;
#'   Koldasbayeva and Zaytsev 2025).
#' * **Environmental separation** (with `env`): the overlap of the
#'   distribution of each predictor between training and testing records
#'   (1 = identical, 0 = no overlap), and the share of test records outside
#'   the training range of at least one predictor (extrapolation).
#' * **Spatial sorting bias** (with `occ`): Hijmans' (2012) ratio of the mean
#'   distance from test presences to the nearest training presence to that
#'   from test background points; see [spatial_sorting_bias()].
#'
#' Duplicate coordinates and shared cells are definite leakage; the
#' distances and overlaps describe how hard the test is and are reported
#' with graded wording.
#'
#' @param points The records that were split: a data.frame with coordinates,
#'   a CSV path, an `sf` object or a [terra::SpatVector], one row per value
#'   of `folds`.
#' @param folds Fold of each record (e.g. from [spatial_folds()] or
#'   [environmental_folds()]); `NA` = not used.
#' @param occ Optional 0/1 vector (1 = presence): counts are given per class
#'   and spatial sorting bias is computed.
#' @param block_size Optional block size (km) used to make the folds; taken
#'   from `attr(folds, "block_size_km")` when present.
#' @param env Optional predictor rasters: shared cells and environmental
#'   separation are checked (values are read only at the records).
#' @param set Optional vector of background sets or replicates: the checks
#'   are made within each set.
#' @param n_bins Number of classes per predictor for the overlap.
#' @inheritParams extract_predictions
#' @return A list of class `fold_report`: `folds` (one row per fold and
#'   set: `n_test`, `n_train` and, with `occ`, `presences_test`,
#'   `presences_train`; `duplicates`, `shared_cells`, `min_dist_km`,
#'   `median_dist_km`, `within_block` (share of test records closer than the
#'   block size to a training record), `separation` (median test-to-training
#'   distance divided by the median distance between neighbouring records),
#'   `env_overlap`, `env_overlap_min`, `outside_range`, `ssb`), `leaks` (one
#'   row per leaking record: `index`, `fold`, `type`, `other_folds`),
#'   `leakage` (`TRUE`/`FALSE`), `overall` and `points` (coordinates and
#'   folds, for [plot_folds()]).
#' @references Stock, A., Gregr, E.J., Chan, K.M.A. (2023) Data leakage
#'   jeopardizes ecological applications of machine learning. *Nature
#'   Ecology & Evolution* 7, 1743-1745. \doi{10.1038/s41559-023-02162-1}
#'
#'   Yoo, M., Koo, H. (2024) Exploring the impact of spatial autocorrelation
#'   on optimistic bias in cross-validation and assessing the effectiveness
#'   of spatial cross-validation. *Cartography and Geographic Information
#'   Science* 52, 596-609. \doi{10.1080/15230406.2024.2422593}
#'
#'   Koldasbayeva, D., Zaytsev, A. (2025) Foundation for unbiased
#'   cross-validation of spatio-temporal models for species distribution
#'   modeling. *Ecological Informatics* 92, 103521.
#'   \doi{10.1016/j.ecoinf.2025.103521}
#'
#'   Hijmans, R.J. (2012) Cross-validation of species distribution models:
#'   removing spatial sorting bias and calibration with a null model.
#'   *Ecology* 93, 679-688. \doi{10.1890/11-0826.1}
#' @seealso [spatial_folds()], [environmental_folds()], [evaluate_cv()],
#'   [plot_folds()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
#' f <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)
#' fold_report(pts, f, occ = pts$occ, env = sim$env)
#' # random folds with a duplicated record: leakage
#' pts2 <- rbind(pts, pts[1, ])
#' f2 <- c(sample(1:4, nrow(pts), replace = TRUE), 1 + f[1] %% 4)
#' fold_report(pts2, f2, occ = pts2$occ)
fold_report <- function(points, folds, occ = NULL, block_size = NULL,
                        env = NULL, set = NULL, n_bins = 20, coords = NULL,
                        crs = "EPSG:4326") {
  if (is.null(block_size)) block_size <- attr(folds, "block_size_km")
  pts <- as_points(points, coords, crs)
  n <- nrow(pts)
  if (length(folds) != n) {
    stop("`folds` must have one value per point with coordinates (", n,
         ").", call. = FALSE)
  }
  folds <- as.vector(folds)
  if (!is.null(occ) && length(occ) != n) {
    stop("`occ` must have one value per point (", n, ").", call. = FALSE)
  }
  if (is.null(set)) set <- rep(1L, n)
  if (length(set) != n) stop("`set` must have one value per point.",
                             call. = FALSE)
  xy <- terra::crds(pts)
  key <- paste(signif(xy[, 1], 12), signif(xy[, 2], 12))
  ds <- dist_space(pts)
  cell <- NULL
  E <- NULL
  if (!is.null(env)) {
    r <- as_grid(env)
    pr <- if (terra::same.crs(pts, r)) pts else terra::project(pts, r)
    cell <- terra::cellFromXY(r, terra::crds(pr))
    E <- env_at_cells(r, cell)
  }
  use <- !is.na(folds)

  # leaking records: same coordinates / same cell in more than one fold
  leaks <- list()
  flag <- function(group, type, skip = rep(FALSE, n)) {
    ok <- use & !is.na(group) & !skip
    g <- paste(set, group)[ok]
    f <- folds[ok]
    nf <- tapply(f, g, function(v) length(unique(v)))
    bad <- g %in% names(nf)[nf > 1]
    if (!any(bad)) return(rep(FALSE, n))
    idx <- which(ok)[bad]
    other <- vapply(seq_along(idx), function(j) {
      i <- idx[j]
      same <- which(ok)[g == g[bad][j]]
      paste(sort(setdiff(unique(folds[same]), folds[i])), collapse = ",")
    }, character(1))
    leaks[[length(leaks) + 1]] <<- data.frame(
      index = idx, set = set[idx], fold = folds[idx], type = type,
      other_folds = other, stringsAsFactors = FALSE)
    out <- rep(FALSE, n)
    out[idx] <- TRUE
    out
  }
  dup <- flag(key, "same coordinates")
  shared <- if (is.null(cell)) rep(FALSE, n) else
    flag(cell, "same raster cell", skip = dup)
  leaks <- if (length(leaks)) do.call(rbind, leaks) else
    data.frame(index = integer(0), set = integer(0), fold = integer(0),
               type = character(0), other_folds = character(0))

  overlap <- function(a, b) {
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    if (!length(a) || !length(b)) return(NA_real_)
    rg <- range(c(a, b))
    if (diff(rg) == 0) return(1)
    br <- seq(rg[1], rg[2], length.out = n_bins + 1)
    pa <- tabulate(findInterval(a, br, rightmost.closed = TRUE), n_bins)
    pb <- tabulate(findInterval(b, br, rightmost.closed = TRUE), n_bins)
    sum(pmin(pa / sum(pa), pb / sum(pb)))
  }

  rows <- list()
  for (s in unique(set[use])) {
    in_s <- use & set == s
    # typical spacing: distance from each record to its nearest other record
    Xs <- ds$X[in_s, , drop = FALSE]
    spacing <- if (nrow(Xs) > 1) {
      nn <- dbscan::kNN(Xs, k = 1)$dist[, 1]
      if (ds$lonlat) nn <- 2 * 6371.0088 * asin(pmin(1, nn / (2 * 6371.0088)))
      stats::median(nn)
    } else NA_real_
    for (f in sort(unique(folds[in_s]))) {
      te <- in_s & folds == f
      tr <- in_s & folds != f
      d <- nn_km(ds$X[te, , drop = FALSE], ds$X[tr, , drop = FALSE],
                 ds$lonlat)$dist
      row <- data.frame(set = s, fold = f, n_test = sum(te),
                        n_train = sum(tr))
      if (!is.null(occ)) {
        row$presences_test <- sum(occ[te] == 1)
        row$presences_train <- sum(occ[tr] == 1)
      }
      row$duplicates <- sum(dup[te])
      row$shared_cells <- if (is.null(cell)) NA_integer_ else sum(shared[te])
      row$min_dist_km <- if (length(d)) min(d) else NA_real_
      row$median_dist_km <- if (length(d)) stats::median(d) else NA_real_
      row$within_block <- if (is.null(block_size)) NA_real_ else
        mean(d < block_size)
      row$separation <- if (is.finite(spacing) && spacing > 0) {
        row$median_dist_km / spacing
      } else NA_real_
      if (!is.null(E)) {
        ov <- vapply(seq_len(ncol(E)), function(j) {
          overlap(E[tr, j], E[te, j])
        }, numeric(1))
        row$env_overlap <- mean(ov, na.rm = TRUE)
        row$env_overlap_min <- min(ov, na.rm = TRUE)
        row$env_overlap_var <- colnames(E)[which.min(ov)]
        lo <- apply(E[tr, , drop = FALSE], 2, min, na.rm = TRUE)
        hi <- apply(E[tr, , drop = FALSE], 2, max, na.rm = TRUE)
        Et <- E[te, , drop = FALSE]
        out <- sweep(Et, 2, lo, "<") | sweep(Et, 2, hi, ">")
        row$outside_range <- mean(apply(out, 1, any, na.rm = TRUE))
      }
      if (!is.null(occ)) {
        row$ssb <- ssb_value(ds, te & occ == 1, te & occ != 1,
                             tr & occ == 1)
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  tab <- do.call(rbind, rows)
  if (length(unique(set)) == 1) tab$set <- NULL
  structure(list(
    folds = tab, leaks = leaks, leakage = nrow(leaks) > 0,
    overall = list(n = sum(use), n_folds = length(unique(folds[use])),
                   duplicates = sum(dup), shared_cells = if (is.null(cell))
                     NA_integer_ else sum(shared),
                   cells_checked = !is.null(cell), block_size = block_size,
                   env_checked = !is.null(E)),
    points = data.frame(x = xy[, 1], y = xy[, 2], fold = folds,
                        leak = dup | shared, set = set),
    lonlat = ds$lonlat),
    class = "fold_report")
}

# Spatial sorting bias Dp / Da from logical selections (internal)
ssb_value <- function(ds, pres, bg, train) {
  if (!any(pres) || !any(bg) || !any(train)) return(NA_real_)
  dp <- nn_km(ds$X[pres, , drop = FALSE], ds$X[train, , drop = FALSE],
              ds$lonlat)$dist
  da <- nn_km(ds$X[bg, , drop = FALSE], ds$X[train, , drop = FALSE],
              ds$lonlat)$dist
  if (mean(da) == 0) return(NA_real_)
  mean(dp) / mean(da)
}

#' @export
print.fold_report <- function(x, digits = 2, ...) {
  t <- x$folds
  o <- x$overall
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  fmt <- function(v) formatC(v, digits = digits, format = "f")
  km <- function(v) paste0(signif(v, 3), " km")
  fl <- function(v) paste(sort(unique(v)), collapse = ", ")
  cat("Cross-validation audit: ", o$n, " records in ", o$n_folds,
      " folds", if ("set" %in% names(t)) paste0(", ",
      length(unique(t$set)), " sets"), "\n\nIN SHORT\n", sep = "")
  if (x$leakage) {
    lk <- x$leaks
    parts <- c(
      if (any(lk$type == "same coordinates")) paste0(
        sum(lk$type == "same coordinates"), " records have the same ",
        "coordinates as records in another fold (folds ",
        fl(lk$fold[lk$type == "same coordinates"]), ")"),
      if (any(lk$type == "same raster cell")) paste0(
        sum(lk$type == "same raster cell"), " records share a raster cell ",
        "with records in another fold (folds ",
        fl(lk$fold[lk$type == "same raster cell"]), ")"))
    if (length(parts) == 1) parts <- sub(" \\(folds [^)]*\\)$", "", parts)
    say(paste0("LEAKAGE FOUND in folds ", fl(lk$fold), ": ",
               paste(parts, collapse = "; "), ". A model tested on these ",
               "folds was trained on the same places and predictor values, ",
               "so their metrics are inflated. Keep one record per cell ",
               "(thin_by_cell()) before making folds, or assign folds by ",
               "block (spatial_folds())."))
  } else {
    say(paste0("No leakage found: no record has the same coordinates as",
               if (o$cells_checked) ", or shares a raster cell with,",
               " a record in another fold.",
               if (!o$cells_checked) " Raster cells were not checked; give `env`."))
  }
  sep <- t$separation
  if (any(!is.na(sep))) {
    ms <- stats::median(sep, na.rm = TRUE)
    grade <- if (ms < 1.5) {
      paste0("not separated: test records are about as close to training ",
             "records as neighbouring records are to each other, as in ",
             "random cross-validation, so spatial autocorrelation makes the ",
             "metrics optimistic")
    } else if (ms < 5) {
      "partly separated in space"
    } else "well separated in space"
    say(paste0("The folds are ", grade, ". Test records lie a median ",
               km(stats::median(t$median_dist_km, na.rm = TRUE)),
               " from the nearest training record (closest pair ",
               km(min(t$min_dist_km, na.rm = TRUE)), ")",
               if (!is.null(o$block_size)) paste0(
                 "; ", share_label(mean(t$within_block, na.rm = TRUE)),
                 " of test records have a training record closer than the ",
                 "block size (", o$block_size, " km): blocks separate whole ",
                 "areas, not every pair of records"),
               "."))
  }
  if (isTRUE(o$env_checked)) {
    mo <- mean(t$env_overlap, na.rm = TRUE)
    w <- which.min(t$env_overlap)
    say(paste0("Training and testing environments ",
               if (mo >= 0.8) "are similar" else if (mo >= 0.5)
                 "overlap partly" else "differ strongly",
               ": mean overlap of the predictor distributions ", fmt(mo),
               ", lowest in fold ", t$fold[w], " (", fmt(t$env_overlap[w]),
               "; least overlapping ", t$env_overlap_var[w], "). ",
               if (any(t$outside_range > 0, na.rm = TRUE)) paste0(
                 share_label(mean(t$outside_range, na.rm = TRUE)),
                 " of test records lie outside the training range of at ",
                 "least one predictor (most in fold ",
                 t$fold[which.max(t$outside_range)], "): there the folds ",
                 "test extrapolation.") else paste0(
                 "No test record lies outside the training range of any ",
                 "predictor.")))
  }
  if ("ssb" %in% names(t) && any(!is.na(t$ssb))) {
    m <- stats::median(t$ssb, na.rm = TRUE)
    say(paste0("Spatial sorting bias ", if (m >= 0.8) "is small" else
      if (m >= 0.5) "is moderate" else "is strong",
      ": test presences lie closer to training presences than test ",
      "background does (median SSB ", fmt(m), "). ",
      if (m < 0.8) "AUC is inflated; use spatial_sorting_bias() for a distance-matched background."))
  }
  cat("\nPer fold:\n")
  show <- t[intersect(c("set", "fold", "n_test", "n_train",
                        "presences_test", "duplicates", "shared_cells",
                        "median_dist_km", "within_block", "separation",
                        "env_overlap", "outside_range", "ssb"), names(t))]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "Leakage = the same coordinates or raster cell in training and ",
    "testing. separation = median distance from test to nearest training ",
    "record / median distance between neighbouring records; \"not ",
    "separated\" (< 1.5), \"partly\" (< 5) and \"well\" are rules of thumb. ",
    "env_overlap = overlap of the histograms of each predictor (1 = ",
    "identical), averaged over predictors; \"similar\" (>= 0.8) and ",
    "\"partly\" (>= 0.5) are rules of thumb. SSB = mean distance from test ",
    "presences to the nearest training presence / the same for test ",
    "background (Hijmans 2012: 1 = no bias, near 0 = extreme); the grades ",
    "0.8 and 0.5 are rules of thumb."), width = 78), sep = "\n")
  invisible(x)
}

#' Map of cross-validation folds
#'
#' The records coloured by fold, with leaking records (same coordinates or
#' raster cell in another fold) marked by black crosses. Legend entries give
#' the number and share of records in each fold.
#'
#' @param x A result of [fold_report()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Stock, A., Gregr, E.J., Chan, K.M.A. (2023) Data leakage
#'   jeopardizes ecological applications of machine learning. *Nature
#'   Ecology & Evolution* 7, 1743-1745. \doi{10.1038/s41559-023-02162-1}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
#' f <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)
#' plot_folds(fold_report(pts, f, occ = pts$occ))
plot_folds <- function(x, base_size = 11, file = NULL, width = 7,
                       height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "fold_report")) {
    stop("`x` must come from fold_report().", call. = FALSE)
  }
  d <- x$points[!is.na(x$points$fold), , drop = FALSE]
  if (length(unique(d$set)) > 1) d <- d[d$set == d$set[1], , drop = FALSE]
  lv <- sort(unique(d$fold))
  cnt <- table(factor(d$fold, lv))
  lab <- paste0("Fold ", lv, ": ", as.vector(cnt), " (",
                share_label(as.vector(cnt) / nrow(d)), ")")
  d$Fold <- factor(d$fold, lv, lab)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02", "#a6761d", "#666666")
  pal <- if (length(lv) <= 8) pal[seq_along(lv)] else
    grDevices::hcl.colors(length(lv), "Dark 3")
  nl <- sum(d$leak)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y)) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$Fold), size = 1.4,
                        alpha = 0.8) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, lab),
                                 name = NULL) +
    (if (isTRUE(x$lonlat)) ggplot2::coord_quickmap() else
      ggplot2::coord_equal()) +
    ggplot2::labs(x = NULL, y = NULL, title = "Cross-validation folds",
                  subtitle = paste0(nrow(d), " records in ", length(lv),
                                    " folds", if (nl) paste0(
                                      "; crosses: ", nl, " leaking records (",
                                      share_label(nl / nrow(d)), ")"), ".")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
  if (nl) {
    p <- p + ggplot2::geom_point(data = d[d$leak, , drop = FALSE],
                                 shape = 4, size = 2.6, colour = "black",
                                 stroke = 0.9)
  }
  p
}

#' Assign points to environmentally separated folds
#'
#' Spatial blocks keep nearby records apart, but folds can still share the
#' same environment; environmental blocking instead tests each model on
#' environmental conditions it was not trained on, which matters when the
#' model will be transferred to new conditions (Roberts et al. 2017; Valavi
#' et al. 2019; Koldasbayeva and Zaytsev 2025). The records are placed in
#' the same environmental space as [thin_environment()] (a PCA of the
#' standardised predictors, precipitation amounts log-transformed, from a
#' sample of cells) and clustered into `k` groups by k-means; each cluster is
#' a fold. Folds therefore differ in size, and a fold may lack presences
#' (a warning says so), in which case use fewer folds. Check the result with
#' [fold_report()].
#'
#' @inheritParams spatial_folds
#' @param env Predictor rasters ([terra::SpatRaster] or file paths); values
#'   are read only at the points and at a sample of `pca_n` cells.
#' @param log_vars,pca_n,pc_cumvar,pc_min Environmental space, as in
#'   [thin_environment()].
#' @param n_start Number of random starts of k-means.
#' @return An integer vector of fold numbers (1 to `k`; `NA` for points
#'   without predictor values), with attributes `centers` (cluster centres
#'   on the PCA axes), `size` and `env_space` (axes and variance explained).
#' @references Koldasbayeva, D., Zaytsev, A. (2025) Foundation for unbiased
#'   cross-validation of spatio-temporal models for species distribution
#'   modeling. *Ecological Informatics* 92, 103521.
#'   \doi{10.1016/j.ecoinf.2025.103521}
#'
#'   Roberts, D.R. et al. (2017) Cross-validation strategies for data with
#'   temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*
#'   40, 913-929.
#'
#'   Valavi, R., Elith, J., Lahoz-Monfort, J.J., Guillera-Arroita, G. (2019)
#'   blockCV: an R package for generating spatially or environmentally
#'   separated folds for k-fold cross-validation of species distribution
#'   models. *Methods in Ecology and Evolution* 10, 225-232.
#' @seealso [spatial_folds()], [fold_report()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
#' f <- environmental_folds(pts, sim$env, k = 4, occ = pts$occ, seed = 1)
#' table(f, pts$occ)
#' fold_report(pts, f, occ = pts$occ, env = sim$env)
environmental_folds <- function(points, env, k = 5, occ = NULL,
                                log_vars = NULL, pca_n = 2e5,
                                pc_cumvar = 0.95, pc_min = 2, n_start = 25,
                                seed = NULL, coords = NULL,
                                crs = "EPSG:4326") {
  pts <- as_points(points, coords, crs)
  n <- nrow(pts)
  if (!is.null(occ) && length(occ) != n) {
    stop("`occ` must have one value per point (", n, ").", call. = FALSE)
  }
  r <- as_grid(env)
  es <- env_pca(r, log_vars, pca_n, pc_cumvar, pc_min, "pca", seed)
  to_axes <- env_axes_fun(es$pca, names(r), es$log_vars, "pca", es$n_axes)
  pr <- if (terra::same.crs(pts, r)) pts else terra::project(pts, r)
  E <- env_at_cells(r, terra::cellFromXY(r, terra::crds(pr)))
  ok <- stats::complete.cases(E)
  if (sum(ok) < k) {
    stop("Fewer points with predictor values (", sum(ok), ") than folds.",
         call. = FALSE)
  }
  Z <- to_axes(E[ok, , drop = FALSE])
  km <- with_seed(seed, stats::kmeans(Z, centers = k, nstart = n_start,
                                      iter.max = 100))
  # number the folds along PC1, so the result does not depend on the seed's
  # labelling of the clusters
  relabel <- rank(km$centers[, 1], ties.method = "first")
  fold <- rep(NA_integer_, n)
  fold[ok] <- as.integer(relabel[km$cluster])
  if (any(!ok)) {
    warning(sum(!ok), " point(s) without values in all predictors get no ",
            "fold (NA).", call. = FALSE)
  }
  if (!is.null(occ)) {
    np <- tabulate(fold[occ == 1 & ok], k)
    if (any(np == 0)) {
      warning("Fold(s) ", paste(which(np == 0), collapse = ", "),
              " contain no presences: their environment has no records of ",
              "the species. Use fewer folds.", call. = FALSE)
    }
  }
  centers <- km$centers[order(relabel), , drop = FALSE]
  rownames(centers) <- seq_len(k)
  attr(fold, "centers") <- centers
  attr(fold, "size") <- tabulate(fold, k)
  attr(fold, "env_space") <- es$message
  fold
}

#' Spatial sorting bias and a distance-matched background
#'
#' Test presences usually lie closer to the training presences than test
#' background (or absence) points do, simply because records are clustered.
#' A model then scores well by predicting high values near known records,
#' and AUC is inflated: "spatial sorting bias" (Hijmans 2012). Its size is
#'
#' `SSB = Dp / Da`,
#'
#' where `Dp` is the mean distance from the test presences to the nearest
#' training presence and `Da` the same for the test background points. SSB
#' = 1 means no bias; values near 0 mean extreme bias (the background lies
#' much farther from the training presences than the test presences do).
#'
#' Pairwise distance sampling removes the bias: each test presence is paired
#' with the unused background point whose distance to the nearest training
#' presence is most similar; presences without a background point within
#' `tolerance` (33%, as in Hijmans 2012) of their own distance are dropped,
#' and each background point is used once. Evaluating on the returned
#' presences and matched background gives an AUC free of spatial sorting
#' bias (Hijmans 2012; `dismo::pwdSample()`).
#'
#' @param presence Test presences (data.frame with coordinates, `sf`,
#'   [terra::SpatVector] or CSV path).
#' @param background Test background or absence points.
#' @param train Training presences.
#' @param tolerance Largest relative difference between the distances of a
#'   presence and its matched background point.
#' @param n Number of background points matched to each presence.
#' @param seed Optional seed for the random order in which presences are
#'   matched.
#' @inheritParams extract_predictions
#' @return A list of class `spatial_sorting_bias`: `ssb`, `Dp`, `Da` (km),
#'   `presence_index` and `background_index` (rows of `presence` and
#'   `background` to evaluate on; `background_index[[i]]` pairs with
#'   `presence_index[i]` when `n = 1`), `ssb_matched` (SSB of the matched
#'   sample) and `n_dropped` (presences without a match).
#' @references Hijmans, R.J. (2012) Cross-validation of species distribution
#'   models: removing spatial sorting bias and calibration with a null
#'   model. *Ecology* 93, 679-688. \doi{10.1890/11-0826.1}
#' @seealso [fold_report()], which reports SSB per fold.
#' @export
#' @examples
#' set.seed(1)
#' train <- data.frame(lon = runif(100, 10, 12), lat = runif(100, 50, 52))
#' test <- data.frame(lon = runif(50, 10, 12), lat = runif(50, 50, 52))
#' bg <- data.frame(lon = runif(500, 8, 16), lat = runif(500, 48, 54))
#' s <- spatial_sorting_bias(test, bg, train, seed = 1)
#' s
#' # evaluate on test[s$presence_index, ] and bg[s$background_index, ]
spatial_sorting_bias <- function(presence, background, train,
                                 tolerance = 0.33, n = 1, seed = NULL,
                                 coords = NULL, crs = "EPSG:4326") {
  p <- as_points(presence, coords, crs)
  b <- as_points(background, coords, crs)
  t <- as_points(train, coords, crs)
  if (!terra::same.crs(b, p)) b <- terra::project(b, p)
  if (!terra::same.crs(t, p)) t <- terra::project(t, p)
  sp <- dist_space(p)
  X <- function(v) dist_space(v)$X
  dp <- nn_km(sp$X, X(t), sp$lonlat)$dist
  db <- nn_km(X(b), X(t), sp$lonlat)$dist
  if (!length(dp) || !length(db)) stop("No points.", call. = FALSE)
  # pairwise distance sampling
  free <- rep(TRUE, length(db))
  ord <- with_seed(seed, sample.int(length(dp)))
  keep <- integer(0)
  matched <- list()
  for (i in ord) {
    cand <- which(free)
    if (length(cand) < n) break
    dif <- abs(db[cand] - dp[i])
    pick <- cand[order(dif)[seq_len(n)]]
    tol <- tolerance * dp[i]
    if (all(abs(db[pick] - dp[i]) <= tol + 1e-12)) {
      keep <- c(keep, i)
      matched[[length(matched) + 1]] <- pick
      free[pick] <- FALSE
    }
  }
  o <- order(keep)
  keep <- keep[o]
  matched <- matched[o]
  bi <- if (n == 1) unlist(matched) else matched
  ub <- unlist(matched)
  structure(list(
    ssb = mean(dp) / mean(db), Dp = mean(dp), Da = mean(db),
    presence_index = keep, background_index = bi,
    ssb_matched = if (length(keep)) mean(dp[keep]) / mean(db[ub]) else
      NA_real_,
    n_presence = length(dp), n_background = length(db),
    n_dropped = length(dp) - length(keep), tolerance = tolerance, n = n),
    class = "spatial_sorting_bias")
}

#' @export
print.spatial_sorting_bias <- function(x, digits = 2, ...) {
  f <- function(v) formatC(v, digits = digits, format = "f")
  cat("Spatial sorting bias: SSB = ", f(x$ssb), " (Dp ", signif(x$Dp, 3),
      " km / Da ", signif(x$Da, 3), " km)\n\nIN SHORT\n", sep = "")
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  say(if (x$ssb >= 0.8) {
    "Spatial sorting bias is small: test presences and test background lie at similar distances from the training presences."
  } else if (x$ssb >= 0.5) {
    "Spatial sorting bias is moderate: test presences lie closer to the training presences than the test background does, so AUC is somewhat inflated."
  } else {
    "Spatial sorting bias is strong: test presences lie much closer to the training presences than the test background does, so AUC mostly measures distance to known records."
  })
  say(paste0("The distance-matched sample keeps ", length(x$presence_index),
             " of ", x$n_presence, " test presences and ",
             length(unlist(x$background_index)), " background points; its ",
             "SSB is ", f(x$ssb_matched), ". Evaluate on presence_index ",
             "and background_index for an AUC free of this bias."))
  if (x$n_dropped > x$n_presence / 2) {
    say(paste0("More than half of the presences found no background point at a similar distance: the background does not cover the distances of the presences, so the matched sample is small."))
  }
  cat("", strwrap(paste0(
    "SSB = Dp / Da: mean distance from test presences (Dp) and from test ",
    "background (Da) to the nearest training presence; 1 = no bias, near 0 ",
    "= extreme (Hijmans 2012). Matching: pairwise distance sampling, ",
    "tolerance ", round(100 * x$tolerance), "%. The grades 0.8 and 0.5 are ",
    "rules of thumb."), width = 78), sep = "\n")
  invisible(x)
}
