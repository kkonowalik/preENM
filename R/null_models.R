# Null-model significance of evaluation metrics -----------------------------------

#' Are the models better than expected from where records are collected?
#'
#' A model can reach a high AUC by predicting where people record, not
#' where the species lives. This test compares each metric with its
#' distribution under a null model in which the "presences" are random
#' locations drawn with the same sampling bias: the presence predictions
#' are replaced `n_null` times by the predictions at as many cells drawn
#' from `pool` (with probability proportional to `weight`, e.g. recording
#' effort), and the metrics are recomputed against the same absence or
#' background points. The p-value is the share of null values at least as
#' good as the observed one.
#'
#' Without `weight` the draws are uniform, which tests only whether the
#' model separates the records from random cells of the area (the null AUC
#' is then about 0.5). With recording effort as `weight` the question is
#' the stricter and more useful one: is the model better at finding the
#' species than at finding the places that are visited?
#'
#' This is a null model of the *locations*, applied to models that are
#' already fitted, so it works with predictions of any software. It differs
#' from the null models of Raes and ter Steege (2007) and Bohl et al.
#' (2019), which refit the model to each set of random points; those also
#' account for the flexibility of the algorithm and need the modelling
#' software.
#'
#' @param data A table of evaluation points as for [evaluate_models()]:
#'   `obs` (1 = presence, 0 = absence or background) and one column of
#'   predictions per model.
#' @param pool Predictions at the candidate cells for the random
#'   presences: a data.frame with the same model columns (e.g. from
#'   `terra::spatSample(models, 10000, na.rm = TRUE)` over the calibration
#'   area), or the model rasters themselves ([terra::SpatRaster]), from
#'   which up to `n_pool` cells are sampled.
#' @param weight Sampling weight of the pool cells: a column name of
#'   `pool`, one value per row, or (with rasters) a raster of recording
#'   effort, see [gbif_effort_rasters()]. `NULL` = uniform.
#' @param criteria Metrics to test and their direction, as in
#'   [evaluate_models()].
#' @param metrics Optional named list of your own metric functions
#'   `function(obs, pred)`.
#' @param n_null Number of null replicates.
#' @param n_pool Number of cells sampled when `pool` is a raster.
#' @param seed Seed.
#' @return A data.frame with one row per model and metric: `observed`,
#'   `null_mean`, `null_95` (the value reached or bettered by 5% of the
#'   null replicates) and `p_value`.
#' @references Raes, N., ter Steege, H. (2007) A null-model for
#'   significance testing of presence-only species distribution models.
#'   *Ecography* 30, 727-736. Bohl, C.L., Kass, J.M., Anderson, R.P. (2019)
#'   A new null model approach to quantify performance and significance for
#'   ecological niche models of species distributions. *Journal of
#'   Biogeography* 46, 1101-1111.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models, sim$presence, sim$absence)
#' # recording effort that follows the first environmental gradient
#' effort <- sim$env[[1]]
#' null_model_test(tab, sim$models, weight = effort, n_null = 99, seed = 1)
null_model_test <- function(data, pool, weight = NULL,
                            criteria = c(AUC = "max", MAEbal = "min"),
                            metrics = NULL, n_null = 499, n_pool = 20000,
                            seed = NULL) {
  validate_table(data)
  obs <- data$obs
  models <- setdiff(names(data), "obs")
  extra <- check_extra_metrics(metrics)
  funs <- metric_funs(extra)
  criteria <- check_criteria(criteria, names(funs))
  funs <- funs[names(criteria)]

  if (inherits(pool, "SpatRaster")) {
    miss <- setdiff(models, names(pool))
    if (length(miss)) {
      stop("`pool` lacks model layer(s): ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    r <- pool[[models]]
    if (inherits(weight, "SpatRaster")) {
      w <- weight[[1]]
      if (!isTRUE(terra::compareGeom(w, r, stopOnError = FALSE))) {
        w <- terra::resample(terra::project(w, terra::crs(r)), r)
      }
      names(w) <- ".weight"
      r <- c(r, w)
    }
    valid <- valid_cells(r)
    s <- with_seed(seed, sort(valid[sample.int(length(valid),
                                               min(length(valid), n_pool))]))
    pool <- as.data.frame(read_cells(r, s))
    names(pool) <- names(r)
    if (inherits(weight, "SpatRaster")) weight <- ".weight"
  }
  pool <- as.data.frame(pool)
  miss <- setdiff(models, names(pool))
  if (length(miss)) {
    stop("`pool` lacks model column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  w <- if (is.null(weight)) {
    rep(1, nrow(pool))
  } else if (is.character(weight)) {
    if (length(weight) != 1 || !weight %in% names(pool)) {
      stop("Weight column not found in `pool`.", call. = FALSE)
    }
    pool[[weight]]
  } else if (length(weight) == nrow(pool)) {
    as.numeric(weight)
  } else {
    stop("`weight` must be a column of `pool`, one value per row or a ",
         "raster.", call. = FALSE)
  }
  use <- stats::complete.cases(pool[models]) & !is.na(w) & w > 0
  if (any(w < 0, na.rm = TRUE)) stop("Negative weights.", call. = FALSE)
  P <- as.matrix(pool[use, models, drop = FALSE])
  w <- w[use]
  i1 <- which(obs == 1)
  if (nrow(P) < length(i1)) {
    stop("`pool` has fewer usable rows (", nrow(P), ") than there are ",
         "presences (", length(i1), ").", call. = FALSE)
  }

  observed <- vapply(models, function(m) compute_metrics(obs, data[[m]],
                                                         funs),
                     numeric(length(funs)))
  observed <- matrix(observed, nrow = length(funs),
                     dimnames = list(names(funs), models))
  # null[metric, model, replicate]
  null <- with_seed(seed, vapply(seq_len(n_null), function(b) {
    draw <- weighted_draw(w, length(i1))
    vapply(models, function(m) {
      pred <- data[[m]]
      pred[i1] <- P[draw, m]
      compute_metrics(obs, pred, funs)
    }, numeric(length(funs)))
  }, matrix(0, length(funs), length(models))))
  null <- array(null, c(length(funs), length(models), n_null))

  out <- do.call(rbind, lapply(seq_along(models), function(j) {
    do.call(rbind, lapply(seq_along(funs), function(i) {
      k <- names(funs)[i]
      o <- observed[i, j]
      v <- null[i, j, ]
      # larger = better, whatever the direction of the metric
      good <- function(z) switch(criteria[[k]], max = z, min = -z,
                                 zero = -abs(z))
      q <- stats::quantile(good(v), 0.95, na.rm = TRUE, names = FALSE)
      data.frame(model = models[j], metric = k, observed = o,
                 null_mean = mean(v, na.rm = TRUE),
                 null_95 = switch(criteria[[k]], max = q, min = -q,
                                  zero = -q),
                 p_value = (1 + sum(good(v) >= good(o), na.rm = TRUE)) /
                   (1 + sum(!is.na(v))))
    }))
  }))
  rownames(out) <- NULL
  attr(out, "n_null") <- n_null
  attr(out, "weighted") <- !is.null(weight)
  out
}
