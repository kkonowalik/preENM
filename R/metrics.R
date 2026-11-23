# Evaluation metrics -----------------------------------------------------------
# Implemented in base R (no dependency on the 'Metrics' package) and numerically
# identical to Metrics::auc(), Metrics::mae() and Metrics::bias().

#' Area under the ROC curve
#'
#' Computed with the Mann-Whitney U statistic, so ties are handled correctly
#' and no threshold is needed.
#'
#' @param obs Numeric vector of observations (1 = presence, 0 = absence).
#' @param pred Numeric vector of predicted suitability.
#' @return A single number between 0 and 1.
#' @export
#' @examples
#' auc_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
auc_score <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(pred)
  (sum(r[obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Mean absolute error
#'
#' Mean absolute difference between the observations (0/1) and the predicted
#' suitability. Lower is better.
#'
#' @inheritParams auc_score
#' @export
#' @examples
#' mae_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
mae_score <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

#' Balanced mean absolute error (MAEbal)
#'
#' The mean absolute error computed separately for presences and for
#' absences/background points and then averaged:
#' `(mean(1 - pred | presence) + mean(pred | absence)) / 2`. The plain MAE
#' is dominated by whichever class has more points; with many background
#' points (e.g. 10 000) it mostly rewards predicting low values everywhere,
#' i.e. narrow models. MAEbal weighs both classes equally (as with a
#' prevalence of 0.5) and so does not depend on the number of background
#' points. Lower is better.
#'
#' @inheritParams auc_score
#' @export
#' @examples
#' obs <- c(1, 1, rep(0, 98))
#' narrow <- c(0.6, 0.6, rep(0, 98))
#' mae_score(obs, narrow)      # looks excellent: background dominates
#' mae_balanced(obs, narrow)   # presences count as much as background
mae_balanced <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (!any(obs == 1) || !any(obs == 0)) return(NA_real_)
  0.5 * mean(1 - pred[obs == 1]) + 0.5 * mean(pred[obs == 0])
}

#' Bias
#'
#' Mean of `obs - pred`. Positive values mean the model under-predicts
#' suitability on average and negative values mean it over-predicts. Values
#' closer to 0 are better.
#'
#' @inheritParams auc_score
#' @export
#' @examples
#' bias_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
bias_score <- function(obs, pred) mean(obs - pred, na.rm = TRUE)

#' True skill statistic (TSS)
#'
#' TSS = sensitivity + specificity - 1, ranging from -1 to 1 (0 = random).
#' By default the threshold that maximises TSS is used (the same as the
#' maximum of Youden's J). This is the usual way TSS is reported, but it is
#' slightly optimistic on hold-out data. Supply a `threshold` fixed beforehand
#' (e.g. from calibration data, see [threshold_p10()]) for a strict hold-out
#' value.
#'
#' @inheritParams auc_score
#' @param threshold Optional fixed threshold; predictions `>= threshold` count
#'   as presences. `NULL` (default) uses the TSS-maximising threshold.
#' @return A single number. The threshold used is stored in
#'   `attr(, "threshold")`.
#' @export
#' @examples
#' tss_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
#' tss_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1), threshold = 0.5)
tss_score <- function(obs, pred, threshold = NULL) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  if (!is.null(threshold)) {
    out <- mean(pred[obs == 1] >= threshold) +
      mean(pred[obs == 0] < threshold) - 1
    attr(out, "threshold") <- threshold
    return(out)
  }
  o <- order(pred, decreasing = TRUE)
  p <- pred[o]
  tp <- cumsum(obs[o] == 1)
  fp <- cumsum(obs[o] == 0)
  last <- c(p[-1] != p[-length(p)], TRUE)   # last index of each tied value
  tss <- tp[last] / n1 - fp[last] / n0
  best <- which.max(tss)
  out <- max(0, tss[best])
  attr(out, "threshold") <- if (tss[best] > 0) p[last][best] else Inf
  out
}

#' Continuous Boyce index
#'
#' Measures how well predictions match the *frequency* of presences, which
#' is especially suited to presence-only data (Hirzel et al. 2006). A moving
#' window slides along the range of predicted values; in each window the
#' ratio of the share of presences to the share of background points
#' (predicted/expected ratio) is computed. The index is the Spearman
#' correlation between this ratio and the window position: values close to 1
#' mean presences become steadily more frequent as suitability increases,
#' 0 means no better than random and negative values mean a wrong model. The
#' implementation reproduces `ecospat::ecospat.boyce()` (moving window,
#' Spearman correlation), which biomod2 also uses.
#'
#' Use background points that represent the whole study area (or all cells);
#' absences sampled away from presences make the index unreliable.
#'
#' @param presence Predicted values at presence points.
#' @param background Predicted values at background (or absence) points, or
#'   at all cells of the study area.
#' @param n_windows Number of window positions (default 100).
#' @param window_width Window width as a proportion of the range of the
#'   background predictions (default 0.1).
#' @return A single number between -1 and 1, or `NA` when it cannot be
#'   computed (e.g. too few distinct window ratios).
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A.
#'   (2006) Evaluating the ability of habitat suitability models to predict
#'   species presences. *Ecological Modelling* 199, 142-152.
#' @export
#' @examples
#' set.seed(1)
#' bg <- runif(1000)
#' pres <- rbeta(100, 4, 1)          # presences mostly at high suitability
#' boyce_index(pres, bg)
boyce_index <- function(presence, background, n_windows = 100,
                        window_width = 0.1) {
  presence <- presence[!is.na(presence)]
  background <- background[!is.na(background)]
  if (length(presence) < 2 || length(background) < 2) return(NA_real_)
  lo <- min(background, presence)
  hi <- max(background, presence)
  w <- (max(background) - min(background)) * window_width
  if (hi == lo || w <= 0 || hi - lo - w <= 0) return(NA_real_)
  start <- seq(from = lo, to = hi - w, by = (hi - lo - w) / n_windows)
  # reproduces ecospat exactly: the last window start is shifted by +1
  # (an ecospat "trick" for closed intervals that in effect drops it)
  start[length(start)] <- start[length(start)] + 1
  start_end <- start + w
  # share of values inside each closed window [start, start_end], counted
  # on the sorted values by binary search (all windows at once) instead of
  # comparing every value with every window
  share <- function(v) {
    v <- sort(v)
    (findInterval(start_end, v) - findInterval(start, v, left.open = TRUE)) /
      length(v)
  }
  # p > 0, e = 0 gives Inf (kept); 0/0 = NaN
  ratio <- round(share(presence) / share(background), 10)
  keep <- which(!is.nan(ratio))
  ratio <- ratio[keep]
  start <- start[keep]
  if (length(ratio) < 2) return(NA_real_)
  # drop successive duplicated ratios (flat stretches), exactly as ecospat
  r <- which(ratio != c(ratio[-1], TRUE))
  if (length(r) < 2) return(NA_real_)
  out <- suppressWarnings(stats::cor(ratio[r], start[r], method = "spearman"))
  if (is.na(out)) NA_real_ else out
}

#' Threshold at a given percentile of presence predictions
#'
#' The "10th percentile training presence" threshold: the suitability value
#' exceeded by 90% of presence points. Usually computed on the calibration
#' (training) presences and then applied to independent data, e.g. to compute
#' sensitivity of a hindcast in [temporal_validation()].
#'
#' @param pred_presence Predicted values at presence points.
#' @param p Proportion of presences allowed below the threshold (0.1 = 10%).
#' @export
#' @examples
#' threshold_p10(c(0.2, 0.5, 0.6, 0.8, 0.9))
threshold_p10 <- function(pred_presence, p = 0.1) {
  unname(stats::quantile(pred_presence, p, na.rm = TRUE))
}

# Compute all metrics for one model (internal)
metric_row <- function(obs, pred, extra = list()) {
  ok <- !is.na(pred)
  c(compute_metrics(obs, pred, metric_funs(extra)),
    n_presence = sum(obs[ok] == 1),
    n_absence = sum(obs[ok] == 0))
}

# Threshold rules ----------------------------------------------------------------

threshold_rule_names <- c("maxSSS", "maxTSS", "equal_ss", "maxMCC", "maxF",
                          "prevalence", "fixed_sens", "p10", "mutation_point")

#' Thresholds of a continuous prediction under several rules
#'
#' A continuous prediction is turned into a map of suitable and unsuitable
#' cells by a threshold, and the rules for choosing it disagree enough to
#' change the mapped range considerably. Instead of one rule, report several
#' together (Hellegers et al. 2025; Yu et al. 2024). For each rule the
#' function gives the threshold and the resulting sensitivity (share of
#' presences predicted suitable), specificity (share of absences or
#' background points predicted unsuitable), TSS and the predicted-suitable
#' share (share of absence/background points at or above the threshold;
#' with background points that represent the study area, the share of the
#' area mapped as suitable).
#'
#' Rules (`rules`):
#'
#' * `maxSSS`: maximum of sensitivity + specificity (Liu et al. 2013,
#'   recommended for presence-background data).
#' * `maxTSS`: maximum of TSS = sensitivity + specificity - 1. It is the
#'   same criterion as `maxSSS` and gives the same threshold; both are
#'   reported because both names are in use.
#' * `equal_ss`: sensitivity equal to specificity, i.e. the smallest
#'   difference between them (DSS in Hellegers et al. 2025; Liu et al.
#'   2005).
#' * `maxMCC`: maximum of Matthews' correlation coefficient,
#'   `(TP TN - FP FN) / sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN))`.
#' * `maxF`: maximum of the F-measure, `2 TP / (FN + 2 TP + FP)`.
#'
#'   MCC and F penalise overprediction (false positives relative to the
#'   predicted presences), so they give smaller ranges than TSS or DSS: in
#'   1671 European plants the median range under maxTSS was 1.3, 3.6 and
#'   9.2 times that under DSS, MCC and F (Hellegers et al. 2025). Both
#'   depend on the prevalence of the evaluation data, so with background
#'   points instead of absences they depend on how many background points
#'   there are; use them with true or pseudo-absences.
#' * `prevalence`: the predicted prevalence equals the prevalence. With
#'   `prevalence = NULL`, the observed prevalence of the evaluation data
#'   (presences / all points) and the predicted prevalence over all points
#'   (Liu et al. 2005); meaningful with true absences. With a given
#'   `prevalence` (e.g. an estimate of the share of the area the species
#'   occupies), the threshold at which that share of the absence/background
#'   points is predicted suitable.
#' * `fixed_sens`: the highest threshold that keeps the sensitivity at
#'   least at `sensitivity` (default 0.95): the presence prediction at that
#'   rank. The default differs from 0.90 because that is `p10`.
#' * `p10`: the 10th percentile of the presence predictions
#'   ([threshold_p10()], an interpolated quantile, so the sensitivity can
#'   be slightly below 0.90), i.e. about a fixed sensitivity of 0.90.
#'   Usually computed on the training presences; here on the presences
#'   given.
#' * `mutation_point`: the "mutation point" of Boyce's predicted-to-expected
#'   (P/E) curve, found by threshold regression (BTQR, Yu et al. 2024). The
#'   range of the predictions is cut into `n_classes` classes of equal width
#'   (100 in most studies Yu et al. review); for each class the P/E ratio is
#'   the share of the presences in the class divided by the share of the
#'   background (the area) in it, plotted against the mean prediction of
#'   the class. A threshold regression without intercept, `P/E = delta1 * s`
#'   up to the threshold and `delta2 * s` above it, is fitted for every
#'   candidate threshold among the class means between their 10th and 90th
#'   percentiles (10% trimming), and the one with the smallest residual sum
#'   of squares is the threshold. It needs no absences: the background
#'   stands for the study area. Yu et al. relate it, through source-sink
#'   theory, to the transition to high-quality habitat, so it is usually
#'   strict.
#'
#' Sensitivity and specificity are computed with `prediction >= threshold`
#' as suitable. Ties among candidate thresholds are resolved towards the
#' lowest threshold reaching the optimum.
#'
#' @param presence Predictions at presence points, or a table from
#'   [extract_predictions()] (a data.frame with `obs` and one column per
#'   model): then `background` is not needed and every model gets its rows.
#' @param background Predictions at absence or background points.
#' @param rules Rules to compute (default all, see Details).
#' @param sensitivity Sensitivity kept by the `fixed_sens` rule.
#' @param prevalence Optional prevalence for the `prevalence` rule.
#' @param n_classes Number of classes of the P/E curve for the
#'   `mutation_point` rule.
#' @return A data.frame of class `enm_thresholds` with `model` (if several),
#'   `rule`, `threshold`, `sensitivity`, `specificity`, `TSS`, `MCC`, `F`
#'   and `suitable_share`.
#' @references Hellegers, M., van Hinsberg, A., Lenoir, J., Dengler, J.,
#'   Huijbregts, M.A.J., Schipper, A.M. (2025) Multiple threshold-selection
#'   methods are needed to binarise species distribution model predictions.
#'   *Diversity and Distributions* 31, e70019. \doi{10.1111/ddi.70019}
#'
#'   Yu, B., Dai, W., Li, S., Wu, Z., Wang, J. (2024) A new threshold
#'   selection method for species distribution models with presence-only
#'   data: extracting the mutation point of the P/E curve by threshold
#'   regression. *Ecology and Evolution* 14, e11208.
#'   \doi{10.1002/ece3.11208}
#'
#'   Liu, C., Berry, P.M., Dawson, T.P., Pearson, R.G. (2005) Selecting
#'   thresholds of occurrence in the prediction of species distributions.
#'   *Ecography* 28, 385-393.
#'
#'   Liu, C., White, M., Newell, G. (2013) Selecting thresholds for the
#'   prediction of species occurrence with presence-only data. *Journal of
#'   Biogeography* 40, 778-789.
#' @seealso [threshold_agreement()], [plot_thresholds()], [threshold_p10()]
#' @export
#' @examples
#' set.seed(1)
#' pres <- rbeta(200, 5, 2)          # predictions at presences
#' bg <- rbeta(2000, 2, 4)           # predictions at background points
#' th <- thresholds(pres, bg)
#' th
#' plot_thresholds(th)
#'
#' # all models of an extract_predictions() table
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models, sim$presence, sim$absence)
#' thresholds(tab, rules = c("maxSSS", "p10", "mutation_point"))
thresholds <- function(presence, background = NULL,
                       rules = threshold_rule_names, sensitivity = 0.95,
                       prevalence = NULL, n_classes = 100) {
  bad <- setdiff(rules, threshold_rule_names)
  if (length(bad)) {
    stop("Unknown rule(s): ", paste(bad, collapse = ", "), ". Available: ",
         paste(threshold_rule_names, collapse = ", "), ".", call. = FALSE)
  }
  if (!is.numeric(sensitivity) || length(sensitivity) != 1 ||
      sensitivity <= 0 || sensitivity >= 1) {
    stop("`sensitivity` must be a number between 0 and 1.", call. = FALSE)
  }
  if (!is.null(prevalence) && (!is.numeric(prevalence) ||
                               length(prevalence) != 1 || prevalence <= 0 ||
                               prevalence >= 1)) {
    stop("`prevalence` must be a number between 0 and 1.", call. = FALSE)
  }
  if (is.data.frame(presence)) {
    validate_table(presence)
    obs <- presence$obs
    mods <- setdiff(names(presence), "obs")
    out <- do.call(rbind, lapply(mods, function(m) {
      d <- threshold_table(presence[[m]][obs == 1], presence[[m]][obs == 0],
                           rules, sensitivity, prevalence, n_classes)
      cbind(model = m, d, stringsAsFactors = FALSE)
    }))
  } else {
    if (is.null(background)) {
      stop("Give `background` (predictions at absence/background points).",
           call. = FALSE)
    }
    out <- threshold_table(presence, background, rules, sensitivity,
                           prevalence, n_classes)
  }
  rownames(out) <- NULL
  attr(out, "settings") <- list(sensitivity = sensitivity,
                                prevalence = prevalence)
  class(out) <- c("enm_thresholds", "data.frame")
  out
}

# One model (internal)
threshold_table <- function(pres, bg, rules, sensitivity, prevalence,
                            n_classes) {
  pres <- pres[!is.na(pres)]
  bg <- bg[!is.na(bg)]
  if (length(pres) < 2 || length(bg) < 2) {
    stop("At least two presence and two absence/background predictions ",
         "are needed.", call. = FALSE)
  }
  ps <- sort(pres)
  bs <- sort(bg)
  sens_at <- function(t) 1 - findInterval(t, ps, left.open = TRUE) /
    length(ps)
  spec_at <- function(t) findInterval(t, bs, left.open = TRUE) / length(bs)
  cand <- sort(unique(c(pres, bg)))
  se <- sens_at(cand)
  sp <- spec_at(cand)
  n1 <- length(ps)
  n0 <- length(bs)
  conf_metrics <- function(se, sp) {
    tp <- se * n1
    fn <- n1 - tp
    tn <- sp * n0
    fp <- n0 - tn
    den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    list(MCC = ifelse(den > 0, (tp * tn - fp * fn) / den, 0),
         F = ifelse(tp > 0, 2 * tp / (fn + 2 * tp + fp), 0))
  }
  cm <- conf_metrics(se, sp)
  first_max <- function(v) cand[which(v == max(v))[1]]
  rule_value <- function(r) {
    switch(r,
      maxSSS = first_max(round(se + sp, 12)),
      maxTSS = first_max(round(se + sp - 1, 12)),
      equal_ss = first_max(-round(abs(se - sp), 12)),
      maxMCC = first_max(round(cm$MCC, 12)),
      maxF = first_max(round(cm$F, 12)),
      prevalence = {
        if (is.null(prevalence)) {
          prev <- length(pres) / (length(pres) + length(bg))
          allv <- sort(c(pres, bg))
          pred_prev <- 1 - findInterval(cand, allv, left.open = TRUE) /
            length(allv)
          first_max(-round(abs(pred_prev - prev), 12))
        } else {
          first_max(-round(abs((1 - sp) - prevalence), 12))
        }
      },
      # the presence value at that rank, so that at least `sensitivity` of
      # the presences are at or above it
      fixed_sens = ps[floor((1 - sensitivity) * length(ps) + 1e-9) + 1],
      p10 = threshold_p10(pres, 0.1),
      mutation_point = mutation_point(pres, bg, n_classes))
  }
  thr <- vapply(rules, rule_value, numeric(1))
  se_t <- sens_at(thr)
  sp_t <- spec_at(thr)
  cm_t <- conf_metrics(se_t, sp_t)
  data.frame(rule = rules, threshold = unname(thr),
             sensitivity = se_t, specificity = sp_t, TSS = se_t + sp_t - 1,
             MCC = cm_t$MCC, F = cm_t$F, suitable_share = 1 - sp_t,
             stringsAsFactors = FALSE)
}

# Predicted-to-expected curve in classes of equal width (internal), as in
# Yu et al. (2024): P = share of the presences in the class, E = share of the
# background (the area); plotted against the mean prediction of the class
pe_curve <- function(presence, background, n_classes = 100) {
  rg <- range(c(presence, background))
  if (diff(rg) <= 0) return(data.frame(mean = numeric(0), ratio = numeric(0)))
  br <- seq(rg[1], rg[2], length.out = n_classes + 1)
  cp <- findInterval(presence, br, rightmost.closed = TRUE)
  cb <- findInterval(background, br, rightmost.closed = TRUE)
  P <- tabulate(cp, n_classes) / length(presence)
  E <- tabulate(cb, n_classes) / length(background)
  m <- as.vector(tapply(background, factor(cb, seq_len(n_classes)), mean))
  ok <- E > 0
  data.frame(mean = m[ok], ratio = P[ok] / E[ok])
}

# Mutation point of the P/E curve by threshold regression (BTQR, Yu et al.
# 2024) (internal): two regimes through the origin, 10% trimming
mutation_point <- function(presence, background, n_classes = 100,
                           trim = 0.1) {
  pe <- pe_curve(presence, background, n_classes)
  if (nrow(pe) < 5) return(NA_real_)
  s <- pe$mean
  y <- pe$ratio
  lim <- stats::quantile(s, c(trim, 1 - trim), names = FALSE)
  cand <- s[s >= lim[1] & s <= lim[2]]
  if (!length(cand)) return(NA_real_)
  ssr <- vapply(cand, function(g) {
    lo <- s <= g
    rss <- 0
    for (k in list(lo, !lo)) {
      if (!any(k)) next
      b <- sum(s[k] * y[k]) / sum(s[k]^2)
      rss <- rss + sum((y[k] - b * s[k])^2)
    }
    rss
  }, numeric(1))
  cand[which.min(ssr)]
}

#' @export
print.enm_thresholds <- function(x, digits = 3, ...) {
  d <- as.data.frame(unclass(x))
  pct <- function(z) paste0(round(100 * z), "%")
  by_model <- if ("model" %in% names(d)) split(d, factor(d$model,
                                                         unique(d$model))) else
    list(d)
  cat("Thresholds of the prediction under ", length(unique(d$rule)),
      " rules\n\n", sep = "")
  show <- d
  num <- vapply(show, is.numeric, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("\nIN SHORT\n")
  for (m in by_model) {
    m <- m[!is.na(m$threshold) & m$rule != "maxTSS", , drop = FALSE]
    if (nrow(m) < 2) next
    lo <- m[which.min(m$suitable_share), ]
    hi <- m[which.max(m$suitable_share), ]
    ratio <- if (lo$suitable_share > 0) hi$suitable_share /
      lo$suitable_share else Inf
    txt <- paste0(
      if (!is.null(m$model)) paste0(m$model[1], ": the") else "The",
      " area mapped as suitable ranges from ", pct(lo$suitable_share),
      " of the background (", lo$rule, ", threshold ",
      signif(lo$threshold, digits), ") to ", pct(hi$suitable_share), " (",
      hi$rule, ", threshold ", signif(hi$threshold, digits), "), ",
      if (is.finite(ratio)) sprintf("%.1f times as much", ratio) else
        "from none", ". ",
      if (is.finite(ratio) && ratio >= 1.5) {
        "The choice of rule changes the mapped range considerably: report several rules, or the cells suitable under all of them (threshold_agreement())."
      } else {
        "The rules largely agree on the mapped range."
      })
    cat(strwrap(txt, width = 78, exdent = 2, prefix = "  ", initial = "- "),
        sep = "\n")
  }
  cat(strwrap(paste0(
    "suitable_share = share of the absence/background points at or above the ",
    "threshold (with background that represents the study area, the share ",
    "of the area mapped as suitable). maxSSS and maxTSS are the same ",
    "criterion; p10 is about a fixed sensitivity of 0.90; MCC and F depend ",
    "on the number of absence/background points."), width = 78,
    prefix = "  "), sep = "\n")
  invisible(x)
}

#' How many threshold rules call each cell suitable?
#'
#' Applies every threshold of a [thresholds()] result to a prediction raster
#' and counts, per cell, the rules under which the cell is suitable
#' (`prediction >= threshold`). Cells suitable under all rules are the
#' robust core of the mapped range; cells suitable under at least one rule
#' its widest extent. The summary gives the area of each.
#'
#' Rules with the same threshold (e.g. `maxSSS` and `maxTSS`) are counted
#' once, so that the count is not inflated by duplicates.
#'
#' @param raster The prediction raster ([terra::SpatRaster] or file); one
#'   layer, or the layer named after the model.
#' @param thresholds A result of [thresholds()].
#' @param model The model (layer and rows of `thresholds`) when
#'   `thresholds` covers several models.
#' @return A list of class `threshold_agreement` with `count` (an integer
#'   raster: number of rules calling the cell suitable), `rules` (the
#'   distinct rules and thresholds used) and `summary` (cells and area,
#'   km2, suitable under all rules, under most, i.e. more than half, and
#'   under at least one).
#' @seealso [thresholds()], [plot_thresholds()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models, sim$presence, sim$absence)
#' th <- thresholds(tab)
#' ag <- threshold_agreement(sim$models[["GLM"]], th, model = "GLM")
#' ag
#' plot_thresholds(ag)
threshold_agreement <- function(raster, thresholds, model = NULL) {
  if (!inherits(thresholds, "enm_thresholds")) {
    stop("`thresholds` must come from thresholds().", call. = FALSE)
  }
  r <- as_grid(raster)
  d <- as.data.frame(unclass(thresholds))
  if ("model" %in% names(d)) {
    if (is.null(model)) {
      if (length(unique(d$model)) > 1) {
        stop("`thresholds` covers several models; choose one with `model`.",
             call. = FALSE)
      }
      model <- d$model[1]
    }
    d <- d[d$model == model, , drop = FALSE]
    if (!nrow(d)) stop("No thresholds for model \"", model, "\".",
                       call. = FALSE)
  }
  if (terra::nlyr(r) > 1) {
    if (is.null(model) || !model %in% names(r)) {
      stop("`raster` has several layers; give `model` naming one.",
           call. = FALSE)
    }
    r <- r[[model]]
  }
  d <- d[!is.na(d$threshold), , drop = FALSE]
  d <- d[!duplicated(signif(d$threshold, 10)), , drop = FALSE]
  if (!nrow(d)) stop("No threshold to apply.", call. = FALSE)
  count <- sum(terra::rast(lapply(d$threshold, function(t) r >= t)))
  count <- terra::mask(count, r)
  names(count) <- "rules_suitable"
  k <- nrow(d)
  area <- terra::cellSize(count, unit = "km", mask = TRUE)
  tab <- terra::zonal(area, count, fun = "sum", na.rm = TRUE)
  names(tab) <- c("rules", "area_km2")
  cells <- as.data.frame(terra::freq(count))
  n_cells <- stats::setNames(cells$count, cells$value)
  total <- sum(tab$area_km2)
  level <- function(sel) {
    a <- sum(tab$area_km2[sel(tab$rules)])
    c(cells = sum(n_cells[as.character(tab$rules[sel(tab$rules)])],
                  na.rm = TRUE), area_km2 = a, share = a / total)
  }
  summary <- data.frame(
    suitable_under = c("all rules", "most rules", "at least one rule"),
    rbind(level(function(v) v == k), level(function(v) v > k / 2),
          level(function(v) v >= 1)), row.names = NULL)
  structure(list(count = count, rules = d[c("rule", "threshold")],
                 summary = summary, n_rules = k),
            class = "threshold_agreement")
}

#' @export
print.threshold_agreement <- function(x, ...) {
  s <- x$summary
  pct <- function(z) paste0(round(100 * z), "%")
  cat("Agreement of ", x$n_rules, " threshold rules (",
      paste(x$rules$rule, collapse = ", "), ")\n\n", sep = "")
  show <- s
  show$area_km2 <- round(show$area_km2)
  show$share <- pct(show$share)
  print(show, row.names = FALSE)
  all <- s$area_km2[1]
  one <- s$area_km2[3]
  cat("\nIN SHORT\n")
  cat(strwrap(paste0(
    "Suitable under every rule: ", round(all), " km2 (", pct(s$share[1]),
    " of the area); under at least one rule: ", round(one), " km2 (",
    pct(s$share[3]), "). ",
    if (all > 0) sprintf(paste0("The widest mapped range is %.1f times the ",
                                "range on which all rules agree."),
                         one / all) else
      "No cell is suitable under all rules."), width = 78, exdent = 2,
    prefix = "  ", initial = "- "), sep = "\n")
  invisible(x)
}

#' Plot thresholds or their agreement
#'
#' For a [thresholds()] result: the predicted-suitable share under each rule
#' (one panel per model), labelled with the threshold, so the spread of the
#' mapped range between rules is visible at a glance. For a
#' [threshold_agreement()] result: the map of the number of rules that call
#' each cell suitable.
#'
#' @param x A result of [thresholds()] or [threshold_agreement()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @export
#' @examples
#' set.seed(1)
#' th <- thresholds(rbeta(200, 5, 2), rbeta(2000, 2, 4))
#' plot_thresholds(th)
plot_thresholds <- function(x, base_size = 11, file = NULL, width = 7,
                            height = 5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  if (inherits(x, "threshold_agreement")) {
    k <- x$n_rules
    df <- as.data.frame(x$count, xy = TRUE, na.rm = TRUE)
    names(df)[3] <- "n"
    s <- x$summary
    lab <- function(i) {
      sprintf("Suitable under %s: %s km\u00b2 (%s)",
              s$suitable_under[i],
              format(round(s$area_km2[i]), big.mark = ","),
              share_label(s$share[i]))
    }
    return(
      ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y,
                                       fill = factor(.data$n, 0:k))) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_manual(
          values = stats::setNames(c("#f0f0f0", grDevices::colorRampPalette(
            c("#c6dbef", "#08519c"))(k)), 0:k),
          name = "Rules calling\nthe cell suitable", drop = FALSE) +
        ggplot2::coord_equal() +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Agreement of threshold rules",
                      subtitle = paste(lab(1), lab(3), sep = "\n")) +
        theme_enm(base_size, ink, muted))
  }
  if (!inherits(x, "enm_thresholds")) {
    stop("`x` must come from thresholds() or threshold_agreement().",
         call. = FALSE)
  }
  d <- as.data.frame(unclass(x))
  d <- d[!is.na(d$threshold), , drop = FALSE]
  if (!"model" %in% names(d)) d$model <- "prediction"
  d$rule <- factor(d$rule, rev(threshold_rule_names[threshold_rule_names %in%
                                                      d$rule]))
  d$label <- sprintf("threshold %s", signif(d$threshold, 3))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$suitable_share, .data$rule)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$suitable_share,
                                       yend = .data$rule),
                          colour = "#c6dbef", linewidth = 2.5) +
    ggplot2::geom_point(colour = "#08519c", size = 2.6) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label), hjust = -0.15,
                       size = base_size / 3.6, colour = ink) +
    ggplot2::scale_x_continuous(
      labels = function(v) paste0(round(100 * v), "%"),
      expand = ggplot2::expansion(mult = c(0, 0.35))) +
    ggplot2::labs(x = "Share of background predicted suitable", y = NULL,
                  title = "Mapped range under different threshold rules",
                  subtitle = paste0("Bars: share of absence/background ",
                                    "points at or above the threshold.")) +
    theme_enm(base_size, ink, muted)
  if (length(unique(d$model)) > 1) p <- p + ggplot2::facet_wrap(~model)
  p
}
