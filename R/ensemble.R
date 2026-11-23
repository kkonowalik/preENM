# Ensembles of model predictions ------------------------------------------------------

#' Combine model predictions into an ensemble, with their disagreement
#'
#' Takes the prediction rasters of several models (from any software) and an
#' evaluation or selection, and returns the ensemble prediction together with
#' how much the models disagree in every cell. The disagreement is as much
#' the result as the ensemble: where the models disagree, the ensemble is
#' an average of different answers, not a confident prediction.
#'
#' **Which models.** By default only the models chosen by the evaluation:
#' the `selected` models of [select_models()], or the Pareto front
#' (`pareto_front == 1`) of [evaluate_models()]. Dominated or unselected
#' models are left out (`use = "all"` keeps every model of the
#' evaluation). Without an evaluation every supplied model is used.
#'
#' **How they are combined** (`method`):
#' * `"median"` (default): the median prediction per cell; robust to a
#'   single model with extreme values.
#' * `"mean"`: the unweighted mean.
#' * `"pcam"`: the median of the half of the models closest to their common
#'   trend: a principal component analysis of the model predictions (on a
#'   sample of cells) is computed and the models with the highest loadings
#'   on the first component, the consensus axis, are kept (Thuiller 2004;
#'   Araujo et al. 2005). Zhu et al. (2021) found it the most reliable
#'   consensus in transfers to new areas. Needs at least three models.
#' * `"weighted"`: the mean weighted by the compromise score of the
#'   evaluation, `1 - distance`, where `distance` is the normalised distance
#'   of the model to the ideal point on the evaluation criteria (as in
#'   [evaluate_models()] and [pareto_rank()]; all models of the evaluation
#'   define the scale), or by your own `weights`.
#'
#' **Why weighting is not the default.** Weights based on how well models
#' predict the data they were evaluated on reward interpolation, which is
#' not the same as predicting new areas or periods: in transfers across
#' continents, consensus models weighted by AUC on held-out data were worse
#' than the median of the models closest to their common trend (PCAm), and
#' the mean was pulled by models with extreme predictions (Zhu et al. 2021).
#' With spatially separated folds (see [fold_report()]) the evaluation is
#' closer to a transfer, but the evaluation criteria still do not measure
#' transferability.
#'
#' **Disagreement.** For the models used, per cell: the standard deviation
#' of their predictions (`sd`) and the coefficient of variation (`cv` =
#' sd / mean; `NA` where the mean is 0). The printout summarises where the
#' models disagree: overall, by level of the ensemble prediction, the area
#' above levels of disagreement, and which model departs most from the
#' ensemble.
#'
#' **Large rasters.** The rasters are read from disk and combined block by
#' block by terra, so stacks larger than memory work (above `big_cells`
#' cells every intermediate layer goes to temporary files); give `filename` to
#' write the result to a GeoTIFF (otherwise terra keeps it in memory when
#' it fits and in temporary files when it does not). The PCA of `"pcam"`
#' and the summaries use a random sample of `n_sample` cells.
#'
#' @param models Prediction rasters: a [terra::SpatRaster] with one layer
#'   per model (layer names = model names), or raster files and/or folders
#'   (model names = file names without extension; with `pattern` and
#'   `recursive` as in [crop_rasters()]). All must share the same grid.
#' @param evaluation Optional result of [select_models()] or
#'   [evaluate_models()] (or [evaluate_cv()] passed through
#'   [select_models()]); its `model` names must match the rasters.
#' @param method `"median"`, `"mean"`, `"pcam"` or `"weighted"`.
#' @param use `"selected"` (default; the Pareto front for an
#'   [evaluate_models()] result) or `"all"` models of the evaluation.
#' @param weights Optional named weights for `method = "weighted"` (names =
#'   models); default: the compromise score of the evaluation.
#' @param filename Optional GeoTIFF for the result (layers `ensemble`,
#'   `sd`, `cv`).
#' @param overwrite Overwrite `filename`.
#' @param n_sample Number of cells sampled for the PCA of `"pcam"` and for
#'   the summaries.
#' @param levels Levels of `sd` for the area summary (for predictions on
#'   the 0-1 scale; rules of thumb).
#' @param pattern,recursive File search in folders, as in [crop_rasters()].
#' @param seed Seed for the cell sample.
#' @param big_cells Rasters with more cells are combined with every
#'   intermediate layer in temporary files (terra option `todisk`), so
#'   memory use stays bounded.
#' @return A list of class `enm_ensemble`: `raster` (a
#'   [terra::SpatRaster] with layers `ensemble`, `sd` and `cv`), `models`
#'   (the models used), `left_out`, `method`, `weights` (for `"weighted"`),
#'   `loadings` (for `"pcam"`: loading of each model on the consensus axis
#'   and whether it was kept), `summary` (quantiles of `sd` and `cv`),
#'   `by_level` (mean `sd` per class of the ensemble prediction), `area`
#'   (area with `sd` above each level, km^2 and share), `departure` (mean
#'   absolute difference of each model from the ensemble) and `scale01`
#'   (whether the predictions are on the 0-1 scale).
#' @references Zhu, G., Fan, J., Peterson, A.T. (2021) Cautions in weighting
#'   individual ecological niche models in ensemble forecasting. *Ecological
#'   Modelling* 448, 109502. \doi{10.1016/j.ecolmodel.2021.109502}
#'
#'   Thuiller, W. (2004) Patterns and uncertainties of species' range shifts
#'   under climate change. *Global Change Biology* 10, 2020-2027.
#'
#'   Araujo, M.B., Whittaker, R.J., Ladle, R.J., Erhard, M. (2005) Reducing
#'   uncertainty in projections of extinction risk from climate change.
#'   *Global Ecology and Biogeography* 14, 529-538.
#' @seealso [plot_ensemble()], [select_models()], [evaluate_models()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' ev <- evaluate_models(sim$models, sim$presence, sim$absence)
#' sel <- select_models(ev, plot = FALSE)
#' ens <- ensemble_models(sim$models, sel)
#' ens
#' plot_ensemble(ens)
ensemble_models <- function(models, evaluation = NULL,
                            method = c("median", "mean", "pcam", "weighted"),
                            use = c("selected", "all"), weights = NULL,
                            filename = NULL, overwrite = FALSE,
                            n_sample = 1e5, levels = c(0.1, 0.2),
                            pattern = "\\.tif{1,2}$", recursive = FALSE,
                            seed = 1, big_cells = 1e7) {
  method <- match.arg(method)
  use <- match.arg(use)
  x <- ensemble_stack(models, pattern, recursive)
  available <- names(x)

  # models to use
  chosen <- available
  left_out <- character(0)
  score <- NULL
  if (!is.null(evaluation)) {
    if (!is.data.frame(evaluation) || !"model" %in% names(evaluation)) {
      stop("`evaluation` must come from select_models() or ",
           "evaluate_models().", call. = FALSE)
    }
    ev <- as.data.frame(evaluation)
    miss <- setdiff(ev$model, available)
    if (length(miss)) {
      stop("No raster for model(s) of the evaluation: ",
           paste(miss, collapse = ", "), ". Rasters: ",
           paste(available, collapse = ", "), call. = FALSE)
    }
    keep <- if (use == "all") rep(TRUE, nrow(ev)) else
      if ("selected" %in% names(ev)) ev$selected else
        if ("pareto_front" %in% names(ev)) ev$pareto_front %in% 1 else
          stop("`evaluation` has neither `selected` nor `pareto_front`.",
               call. = FALSE)
    chosen <- ev$model[keep]
    left_out <- c(ev$model[!keep], setdiff(available, ev$model))
    crit <- attr(evaluation, "criteria")
    if (!is.null(crit)) {
      crit <- crit[names(crit) %in% names(ev)]
      d <- pareto_core(to_minimise(ev, crit))$distance
      score <- stats::setNames(1 - d, ev$model)
    }
  }
  if (!length(chosen)) stop("No model to combine.", call. = FALSE)
  x <- x[[chosen]]
  # large stacks: keep intermediate layers in temporary files, so memory
  # stays bounded however many layers are combined
  if (terra::ncell(x) > big_cells) {
    old <- terra::terraOptions(print = FALSE)$todisk
    terra::terraOptions(todisk = TRUE)
    on.exit(terra::terraOptions(todisk = old), add = TRUE)
  }
  n <- length(chosen)
  if (n == 1) {
    warning("Only one model (", chosen, "): the ensemble is that model ",
            "and there is no disagreement to measure.", call. = FALSE)
  }

  # cell sample (PCA, summaries)
  smp <- with_seed(seed, sample_complete(x, n_sample))
  if (!nrow(smp)) stop("No cell has predictions of all models.",
                       call. = FALSE)

  w <- NULL
  loadings <- NULL
  if (method == "weighted") {
    w <- if (!is.null(weights)) weights else score
    if (is.null(w)) {
      stop("method = \"weighted\" needs `weights` or an evaluation with ",
           "criteria.", call. = FALSE)
    }
    if (is.null(names(w)) || !all(chosen %in% names(w))) {
      stop("`weights` must be named and cover the models: ",
           paste(chosen, collapse = ", "), call. = FALSE)
    }
    w <- pmax(w[chosen], 0)
    if (sum(w) == 0) {
      warning("All weights are 0; equal weights used.", call. = FALSE)
      w[] <- 1
    }
    w <- w / sum(w)
  }
  members <- chosen
  if (method == "pcam") {
    if (n < 3) stop("method = \"pcam\" needs at least three models.",
                    call. = FALSE)
    sds <- apply(smp, 2, stats::sd)
    if (any(sds == 0)) {
      stop("A model is constant over the sampled cells; use another method.",
           call. = FALSE)
    }
    pc <- stats::prcomp(smp, center = TRUE, scale. = TRUE)
    l1 <- pc$rotation[, 1]
    if (sum(l1) < 0) l1 <- -l1          # orient the consensus axis
    k <- ceiling(n / 2)
    members <- names(sort(l1, decreasing = TRUE))[seq_len(k)]
    loadings <- data.frame(model = names(l1), loading = unname(l1),
                           kept = names(l1) %in% members,
                           stringsAsFactors = FALSE)
    loadings <- loadings[order(-loadings$loading), , drop = FALSE]
    rownames(loadings) <- NULL
  }

  # combine (block by block, on disk when needed)
  ens <- switch(method,
                median = terra::median(x),
                mean = terra::mean(x),
                pcam = if (length(members) == 1) x[[members]] else
                  terra::median(x[[members]]),
                weighted = if (sum(w > 0) == 1) x[[names(w)[w > 0]]] else
                  terra::weighted.mean(x[[names(w)[w > 0]]],
                                       w = unname(w[w > 0])))
  if (n > 1) {
    # deviations from the mean, after shifting by the first model (the
    # standard deviation does not change): exactly 0 for identical models
    y <- x - x[[1]]
    sdv <- sqrt(sum((y - terra::mean(y))^2) / (n - 1))
    mn <- terra::mean(x)
    cvv <- terra::ifel(mn == 0, NA, sdv / mn)
  } else {
    sdv <- ens * 0
    cvv <- ens * 0
  }
  out <- c(ens, sdv, cvv)
  names(out) <- c("ensemble", "sd", "cv")
  if (!is.null(filename)) {
    out <- terra::writeRaster(out, filename, overwrite = overwrite,
                              gdal = c("COMPRESS=LZW", "TILED=YES",
                                       "BIGTIFF=IF_SAFER"))
  }

  # summaries on the sample
  e_s <- switch(method,
                median = apply(smp, 1, stats::median),
                mean = rowMeans(smp),
                pcam = apply(smp[, members, drop = FALSE], 1, stats::median),
                weighted = drop(smp %*% w))
  sd_s <- if (n > 1) apply(smp, 1, stats::sd) else rep(0, nrow(smp))
  mn_s <- rowMeans(smp)
  cv_s <- ifelse(mn_s == 0, NA_real_, sd_s / mn_s)
  scale01 <- all(smp >= 0 & smp <= 1)
  probs <- c(0.5, 0.9, 0.99, 1)
  summ <- data.frame(statistic = c("sd", "cv"),
                     rbind(stats::quantile(sd_s, probs, names = FALSE),
                           stats::quantile(cv_s, probs, na.rm = TRUE,
                                           names = FALSE)))
  names(summ)[-1] <- c("median", "q90", "q99", "max")
  rg <- range(e_s)
  br <- if (scale01) seq(0, 1, by = 0.2) else
    seq(rg[1], rg[2], length.out = 6)
  cls <- cut(e_s, br, include.lowest = TRUE,
             labels = paste(signif(br[-length(br)], 2), signif(br[-1], 2),
                            sep = "-"))
  by_level <- data.frame(
    ensemble = levels(cls),
    share_of_area = as.vector(table(cls)) / length(cls),
    mean_sd = as.vector(tapply(sd_s, cls, mean)),
    stringsAsFactors = FALSE)
  by_level <- by_level[by_level$share_of_area > 0, , drop = FALSE]
  departure <- data.frame(
    model = chosen,
    mean_abs_diff = colMeans(abs(smp - e_s)),
    member = chosen %in% members, stringsAsFactors = FALSE)
  departure <- departure[order(-departure$mean_abs_diff), , drop = FALSE]
  rownames(departure) <- NULL

  # area above each disagreement level (exact, on the raster)
  area <- NULL
  if (scale01 && n > 1) {
    cs <- terra::cellSize(sdv, mask = TRUE, unit = "km")
    tot <- terra::global(cs, "sum", na.rm = TRUE)[[1]]
    area <- do.call(rbind, lapply(levels, function(l) {
      a <- terra::global(terra::ifel(sdv > l, cs, 0), "sum",
                         na.rm = TRUE)[[1]]
      data.frame(sd_above = l, area_km2 = a, share = a / tot)
    }))
  }

  structure(list(raster = out, models = chosen, members = members,
                 left_out = left_out, method = method, weights = w,
                 loadings = loadings, summary = summ, by_level = by_level,
                 area = area, departure = departure, scale01 = scale01,
                 n_sample = nrow(smp)),
            class = "enm_ensemble")
}

# Values of all layers at up to `n` random cells with values in every layer
# (internal): random cell numbers read with terra::extract(), drawing more
# when many cells are empty (much faster than spatSample(na.rm = TRUE) on
# rasters read from files)
sample_complete <- function(x, n) {
  nc <- terra::ncell(x)
  if (nc <= n) {
    v <- as.matrix(terra::extract(x, seq_len(nc)))
    v <- v[stats::complete.cases(v), , drop = FALSE]
    colnames(v) <- names(x)
    return(v)
  }
  got <- matrix(numeric(0), 0, terra::nlyr(x))
  tried <- integer(0)
  for (k in 1:8) {
    left <- n - nrow(got)
    if (left <= 0 || length(tried) >= nc) break
    share <- if (length(tried)) max(nrow(got) / length(tried), 0.01) else 1
    m <- min(nc, ceiling(1.2 * left / share))
    cells <- unique(sample.int(nc, m))
    cells <- cells[!cells %in% tried]
    if (!length(cells)) next
    tried <- c(tried, cells)
    v <- as.matrix(terra::extract(x, sort(cells)))
    got <- rbind(got, v[stats::complete.cases(v), , drop = FALSE])
  }
  if (nrow(got) > n) got <- got[seq_len(n), , drop = FALSE]
  colnames(got) <- names(x)
  got
}

# Prediction rasters as one stack with model names (internal)
ensemble_stack <- function(models, pattern, recursive) {
  if (inherits(models, "SpatRaster")) {
    x <- models
  } else if (is.character(models)) {
    is_dir <- dir.exists(models)
    missing <- models[!is_dir & !file.exists(models)]
    if (length(missing)) {
      stop("Not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    files <- unique(c(models[!is_dir],
                      unlist(lapply(models[is_dir], list.files,
                                    pattern = pattern, full.names = TRUE,
                                    recursive = recursive,
                                    ignore.case = TRUE))))
    if (!length(files)) stop("No raster files found.", call. = FALSE)
    x <- tryCatch(terra::rast(files), error = function(e) {
      stop("The rasters cannot be stacked (do they share the same grid?): ",
           conditionMessage(e), call. = FALSE)
    })
    if (terra::nlyr(x) == length(files)) {
      names(x) <- tools::file_path_sans_ext(basename(files))
    }
  } else {
    stop("`models` must be a SpatRaster or raster files/folders.",
         call. = FALSE)
  }
  if (anyDuplicated(names(x))) {
    stop("Model names must be unique: ",
         paste(unique(names(x)[duplicated(names(x))]), collapse = ", "),
         call. = FALSE)
  }
  x
}

#' @export
print.enm_ensemble <- function(x, digits = 3, ...) {
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(v, digits = 2, format = "f")
  how <- switch(x$method, median = "median", mean = "unweighted mean",
                weighted = "weighted mean",
                pcam = "median of the models closest to the consensus axis (PCAm)")
  cat("Ensemble of ", length(x$models), " models: ", how, "\n", sep = "")
  cat("  Used: ", paste(x$models, collapse = ", "), "\n", sep = "")
  if (length(x$left_out)) {
    cat("  Left out (not selected or dominated): ",
        paste(x$left_out, collapse = ", "), "\n", sep = "")
  }
  if (x$method == "pcam") {
    cat("  Median of: ", paste(x$members, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(x$weights)) {
    cat("  Weights: ", paste(names(x$weights), f(x$weights), sep = " = ",
                             collapse = ", "), "\n", sep = "")
  }
  cat("\nIN SHORT\n")
  n <- length(x$models)
  if (n < 2) {
    say("One model only: there is no disagreement to measure.")
  } else {
    s <- x$summary
    msd <- s$median[s$statistic == "sd"]
    if (x$scale01) {
      grade <- if (msd < 0.05) "agree closely over most of the area" else
        if (msd < 0.1) "mostly agree" else "disagree over much of the area"
      say(paste0("The models ", grade, ": the median standard deviation ",
                 "between models is ", f(msd), " (90% of cells below ",
                 f(s$q90[s$statistic == "sd"]), ")."))
      if (!is.null(x$area)) {
        a <- x$area
        parts <- ifelse(a$area_km2 > 0, paste0(
          a$sd_above, " on ", share_label(a$share), " of the area (",
          formatC(round(a$area_km2), format = "d", big.mark = " "),
          " km2)"), paste0(a$sd_above, " nowhere"))
        say(paste0("The standard deviation exceeds ",
                   paste(parts, collapse = " and "), "."))
      }
    } else {
      say(paste0("The predictions are not on the 0-1 scale: compare the ",
                 "coefficient of variation (median ",
                 f(s$median[s$statistic == "cv"]), ") rather than the ",
                 "standard deviation."))
    }
    bl <- x$by_level
    w <- which.max(bl$mean_sd)
    say(paste0("Disagreement is highest where the ensemble predicts ",
               bl$ensemble[w], " (mean sd ", f(bl$mean_sd[w]), "; ",
               share_label(bl$share_of_area[w]), " of the area), ",
               if (w > 1 && w < nrow(bl)) "intermediate values" else
                 if (w == 1) "the lowest values" else "the highest values",
               "."))
    d <- x$departure
    say(paste0(d$model[1], " departs most from the ensemble (mean ",
               "absolute difference ", f(d$mean_abs_diff[1]), "; next ",
               d$model[min(2, nrow(d))], " ",
               f(d$mean_abs_diff[min(2, nrow(d))]), ")."))
  }
  if (x$method == "weighted") {
    say(paste0("Weights reward performance on the evaluation data, which ",
               "need not mean better predictions in new areas or periods ",
               "(Zhu et al. 2021)."))
  }
  cat("\nDisagreement by level of the ensemble prediction:\n")
  bl <- x$by_level
  bl$share_of_area <- share_label(bl$share_of_area)
  bl$mean_sd <- round(bl$mean_sd, digits)
  print(bl, row.names = FALSE)
  cat("", strwrap(paste0(
    "sd = standard deviation of the predictions of the models used, per ",
    "cell; cv = sd / mean. Shares and levels come from a sample of ",
    x$n_sample, " cells; areas from all cells. The grades (median sd ",
    "below 0.05 \"agree closely\", below 0.1 \"mostly agree\") and the ",
    "levels of sd are rules of thumb for predictions on the 0-1 scale."),
    width = 78), sep = "\n")
  invisible(x)
}

#' Plot an ensemble and the disagreement between its models
#'
#' The ensemble prediction and, next to it, the standard deviation (or
#' coefficient of variation) between the models used, from
#' [ensemble_models()]. Large rasters are drawn from a regular sample of
#' about `max_cells` cells.
#'
#' @param x A result of [ensemble_models()].
#' @param type `"both"` (side by side), `"ensemble"` or `"disagreement"`.
#' @param disagreement `"sd"` or `"cv"`.
#' @param max_cells Largest number of cells drawn.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches (default 13 x 5.5
#'   for both panels, 7 x 5.5 for one).
#' @return A ggplot, or for `type = "both"` invisibly a list of two.
#' @references Zhu, G., Fan, J., Peterson, A.T. (2021) Cautions in weighting
#'   individual ecological niche models in ensemble forecasting. *Ecological
#'   Modelling* 448, 109502. \doi{10.1016/j.ecolmodel.2021.109502}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' ens <- ensemble_models(sim$models)
#' plot_ensemble(ens)
plot_ensemble <- function(x, type = c("both", "ensemble", "disagreement"),
                          disagreement = c("sd", "cv"), max_cells = 5e5,
                          base_size = 11, file = NULL, width = NULL,
                          height = NULL) {
  if (!inherits(x, "enm_ensemble")) {
    stop("`x` must come from ensemble_models().", call. = FALSE)
  }
  type <- match.arg(type)
  disagreement <- match.arg(disagreement)
  if (is.null(width)) width <- if (type == "both") 13 else 7
  if (is.null(height)) height <- 5.5
  r <- x$raster[[c("ensemble", disagreement)]]
  if (terra::ncell(r) > max_cells) {
    r <- terra::spatSample(r, max_cells, method = "regular",
                           as.raster = TRUE)
  }
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
  lonlat <- terra::is.lonlat(r)
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  panel <- function(v, title, sub, pal, lab) {
    dd <- d[!is.na(d[[v]]), c("x", "y", v)]
    names(dd)[3] <- "value"
    ggplot2::ggplot(dd, ggplot2::aes(.data$x, .data$y, fill = .data$value)) +
      ggplot2::geom_raster() +
      ggplot2::scale_fill_distiller(palette = pal, direction = 1,
                                    name = lab) +
      (if (lonlat) ggplot2::coord_quickmap() else ggplot2::coord_equal()) +
      ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = sub) +
      theme_enm(base_size, ink, muted)
  }
  how <- switch(x$method, median = "Median", mean = "Mean",
                weighted = "Weighted mean", pcam = "PCAm median")
  n <- length(x$models)
  p_e <- panel("ensemble", "Ensemble prediction",
               paste0(how, " of ", n, " models",
                      if (x$method == "pcam") paste0(" (", length(x$members),
                                                     " closest to consensus)"),
                      "."), "YlGnBu", "Prediction")
  p_d <- panel(disagreement, "Disagreement between models",
               paste0(if (disagreement == "sd") "Standard deviation" else
                 "Coefficient of variation", " of the ", n,
                 " model predictions per cell."), "OrRd",
               if (disagreement == "sd") "SD" else "CV")
  p <- switch(type, ensemble = p_e, disagreement = p_d,
              list(ensemble = p_e, disagreement = p_d))
  if (!is.null(file)) save_figure(p, file, width, height)
  if (type == "both") {
    print_side_by_side(p)
    return(invisible(p))
  }
  p
}
