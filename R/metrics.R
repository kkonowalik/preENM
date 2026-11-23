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
#' is especially suited to presence-only data (Hirzel et al. 2006). The
#' predicted/expected (P/E) ratio, the share of presences divided by the
#' share of background points, is estimated along the range of predicted
#' values; the index is the Spearman correlation between this ratio and the
#' prediction. Values close to 1 mean presences become steadily more frequent
#' as suitability increases, 0 means no better than random and negative
#' values mean a wrong model.
#'
#' **Two estimators of the P/E curve.**
#' * `smooth = "none"` (default): a moving window slides along the range of
#'   predicted values and the ratio is computed in each window. This
#'   reproduces `ecospat::ecospat.boyce()`, which biomod2 also uses. With few
#'   presences the ratios of neighbouring windows fluctuate, and every
#'   fluctuation lowers the correlation: in simulations the index was
#'   underestimated below about 50-200 presences (Liu et al. 2024).
#' * `smooth = "spline"`, `"gam"` or `"loess"`: the probability that a point
#'   is a presence rather than a background point is modelled as a smooth
#'   function of the prediction, from all points at once (Liu et al. 2024).
#'   Its odds are proportional to the P/E ratio, so the index is the
#'   Spearman correlation between the fitted curve and the prediction at
#'   `n_points` equally spaced values over the range of the predictions.
#'   `"spline"` is a penalised logistic regression spline (P-spline: cubic
#'   B-splines, 10 basis functions, second-order difference penalty chosen
#'   by UBRE, as mgcv does for binomial models), computed without extra
#'   packages; `"gam"` is the thin plate regression spline of Liu et al.
#'   (`mgcv::gam()`, `k = 10`; needs mgcv); both performed equally well in
#'   Liu et al. (2024). `"loess"` is a local regression of the 0/1 outcome
#'   (`span`), not evaluated by Liu et al.; it is less variable for flat
#'   curves but biased towards 0.
#'
#' The smoothed index is much less variable with few presences and agrees
#' with the moving window for models whose P/E ratio rises steadily. Two
#' cautions, which is why the moving window stays the default: for a model
#' no better than random the fitted curve is almost a straight line, and the
#' Spearman correlation of any straight line is +1 or -1, so the smoothed
#' index of a random model swings across the whole range (in 200 simulated
#' random models with 50 presences it was 0.9 or more in a third, the
#' moving window never; check the bootstrap interval and the test of a flat
#' curve reported with `details = TRUE`); good models often score exactly 1
#' and tie, which makes it a weaker Pareto criterion; and for curves that rise and then fall the two
#' estimators can differ, so compare models with the same estimator.
#'
#' Use background points that represent the whole study area (or all cells);
#' absences sampled away from presences make the index unreliable.
#'
#' @param presence Predicted values at presence points.
#' @param background Predicted values at background (or absence) points, or
#'   at all cells of the study area.
#' @param n_windows Number of window positions (default 100), for
#'   `smooth = "none"`.
#' @param window_width Window width as a proportion of the range of the
#'   background predictions (default 0.1), for `smooth = "none"`.
#' @param smooth Estimator of the P/E curve: `"none"` (moving window, the
#'   default), `"spline"`, `"gam"` or `"loess"`; see Details.
#' @param n_boot Number of bootstrap resamples for a confidence interval
#'   (presences and background resampled separately). With `n_boot > 0` the
#'   result has `details`.
#' @param conf Confidence level of the bootstrap interval.
#' @param details Return the P/E curve, the interval and diagnostics (a
#'   `boyce_index` object) instead of a single number.
#' @param n_points Number of equally spaced prediction values at which a
#'   smoothed curve is evaluated (501 as in Liu et al. 2024).
#' @param span Span of `smooth = "loess"`.
#' @param seed Optional random seed for the bootstrap.
#' @return A single number between -1 and 1, or `NA` when it cannot be
#'   computed (e.g. too few distinct window ratios). With `details = TRUE` or
#'   `n_boot > 0`, a list of class `boyce_index` with `index`, `lower`,
#'   `upper` (bootstrap interval, `NA` without bootstrap), `curve` (the P/E
#'   curve: `prediction`, `pe` and, with bootstrap, `pe_lower`, `pe_upper`),
#'   `boot` (the bootstrap values), `flat_p` (smoothed curves: p-value of a
#'   likelihood-ratio test against a flat P/E curve), `n_presence`,
#'   `n_background` and `smooth`; it has a print method and
#'   [plot_boyce()].
#' @references Hirzel, A.H., Le Lay, G., Helfer, V., Randin, C., Guisan, A.
#'   (2006) Evaluating the ability of habitat suitability models to predict
#'   species presences. *Ecological Modelling* 199, 142-152.
#'
#'   Liu, C., Newell, G., White, M., Machunter, J. (2024) Improving the
#'   estimation of the Boyce index using statistical smoothing methods for
#'   evaluating species distribution models with presence-only data.
#'   *Ecography* 2025, e07218. \doi{10.1111/ecog.07218}
#'
#'   Eilers, P.H.C., Marx, B.D. (1996) Flexible smoothing with B-splines and
#'   penalties. *Statistical Science* 11, 89-121.
#' @seealso [plot_boyce()]; the optional metric `Boyce_smooth` in
#'   [evaluate_models()].
#' @export
#' @examples
#' set.seed(1)
#' bg <- runif(1000)
#' pres <- rbeta(100, 4, 1)          # presences mostly at high suitability
#' boyce_index(pres, bg)
#' boyce_index(pres, bg, smooth = "spline")
#' # few presences: P/E curve, bootstrap interval and diagnostics
#' b <- boyce_index(pres[1:20], bg, smooth = "spline", n_boot = 99, seed = 1)
#' b
#' plot_boyce(b)
boyce_index <- function(presence, background, n_windows = 100,
                        window_width = 0.1,
                        smooth = c("none", "spline", "gam", "loess"),
                        n_boot = 0, conf = 0.95, details = n_boot > 0,
                        n_points = 501, span = 0.75, seed = NULL) {
  smooth <- match.arg(smooth)
  if (smooth == "gam" && !requireNamespace("mgcv", quietly = TRUE)) {
    stop("smooth = \"gam\" needs mgcv: install.packages(\"mgcv\"), or use ",
         "smooth = \"spline\".", call. = FALSE)
  }
  presence <- presence[!is.na(presence)]
  background <- background[!is.na(background)]
  rg <- range(c(presence, background))
  one <- function(p, b, curve = FALSE) {
    if (smooth == "none") {
      boyce_windows(p, b, n_windows, window_width, curve)
    } else {
      boyce_smooth(p, b, smooth, rg, n_points, span, curve)
    }
  }
  if (!details && n_boot == 0) return(one(presence, background)$index)

  est <- one(presence, background, curve = TRUE)
  out <- list(index = est$index, lower = NA_real_, upper = NA_real_,
              curve = est$curve, boot = NULL, flat_p = est$flat_p,
              n_presence = length(presence),
              n_background = length(background), smooth = smooth,
              conf = conf, n_boot = n_boot)
  if (n_boot > 0 && !is.na(est$index)) {
    bs <- with_seed(seed, lapply(seq_len(n_boot), function(i) {
      one(sample(presence, replace = TRUE),
          sample(background, replace = TRUE), curve = smooth != "none")
    }))
    out$boot <- vapply(bs, `[[`, numeric(1), "index")
    a <- (1 - conf) / 2
    q <- stats::quantile(out$boot, c(a, 1 - a), na.rm = TRUE, names = FALSE)
    out$lower <- q[1]
    out$upper <- q[2]
    if (smooth != "none") {
      pe <- vapply(bs, function(z) z$curve$pe, numeric(nrow(out$curve)))
      out$curve$pe_lower <- apply(pe, 1, stats::quantile, a, na.rm = TRUE,
                                  names = FALSE)
      out$curve$pe_upper <- apply(pe, 1, stats::quantile, 1 - a,
                                  na.rm = TRUE, names = FALSE)
    }
  }
  structure(out, class = "boyce_index")
}

# Moving-window estimator (internal); reproduces ecospat.boyce()
boyce_windows <- function(presence, background, n_windows, window_width,
                          curve = FALSE) {
  na <- list(index = NA_real_, curve = data.frame(prediction = numeric(0),
                                                  pe = numeric(0)),
             flat_p = NA_real_)
  if (length(presence) < 2 || length(background) < 2) return(na)
  lo <- min(background, presence)
  hi <- max(background, presence)
  w <- (max(background) - min(background)) * window_width
  if (hi == lo || w <= 0 || hi - lo - w <= 0) return(na)
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
  if (curve) {
    ok <- is.finite(ratio)
    na$curve <- data.frame(prediction = (start + w / 2)[ok], pe = ratio[ok])
  }
  keep <- which(!is.nan(ratio))
  ratio <- ratio[keep]
  start <- start[keep]
  if (length(ratio) < 2) return(na)
  # drop successive duplicated ratios (flat stretches), exactly as ecospat
  r <- which(ratio != c(ratio[-1], TRUE))
  if (length(r) < 2) return(na)
  out <- suppressWarnings(stats::cor(ratio[r], start[r], method = "spearman"))
  na$index <- if (is.na(out)) NA_real_ else out
  na
}

# Smoothed estimator (internal, Liu et al. 2024): P(presence | prediction)
# modelled from all points; its odds times n0/n1 is the P/E ratio. `rg` is
# the range of the original predictions (kept fixed in the bootstrap).
boyce_smooth <- function(presence, background, method, rg, n_points = 501,
                         span = 0.75, curve = FALSE) {
  grid <- seq(0, 1, length.out = n_points)
  na <- list(index = NA_real_, curve = data.frame(
    prediction = rg[1] + grid * diff(rg), pe = NA_real_), flat_p = NA_real_)
  n1 <- length(presence)
  n0 <- length(background)
  if (n1 < 2 || n0 < 2 || diff(rg) <= 0) return(na)
  x <- (c(presence, background) - rg[1]) / diff(rg)
  y <- rep(c(1, 0), c(n1, n0))
  if (length(unique(x)) < 4) return(na)
  dev_null <- -2 * (n1 * log(n1 / (n1 + n0)) + n0 * log(n0 / (n1 + n0)))
  if (method == "spline") {
    f <- pspline_logit(x, y)
    if (is.null(f)) return(na)
    eta <- drop(bspline_basis(grid) %*% f$beta)
    dev <- f$deviance
    edf <- f$edf
  } else if (method == "gam") {
    d <- data.frame(x = x, y = y)
    m <- tryCatch(suppressWarnings(mgcv::gam(
      y ~ s(x, k = 10, bs = "tp"), family = stats::binomial(), data = d)),
      error = function(e) NULL)
    if (is.null(m)) return(na)
    eta <- as.vector(stats::predict(m, data.frame(x = grid)))
    dev <- m$deviance
    edf <- sum(m$edf)
  } else {
    m <- tryCatch(suppressWarnings(stats::loess(
      y ~ x, span = span, degree = 2,
      control = stats::loess.control(surface = "direct"))),
      error = function(e) NULL)
    if (is.null(m)) return(na)
    p <- as.vector(stats::predict(m, data.frame(x = grid)))
    p <- pmin(pmax(p, 1e-9), 1 - 1e-9)
    eta <- stats::qlogis(p)
    pf <- pmin(pmax(stats::fitted(m), 1e-9), 1 - 1e-9)
    dev <- -2 * sum(y * log(pf) + (1 - y) * log(1 - pf))
    edf <- m$enp
  }
  idx <- suppressWarnings(stats::cor(eta, grid, method = "spearman"))
  na$index <- if (is.na(idx)) NA_real_ else idx
  if (curve) na$curve$pe <- exp(eta) * n0 / n1
  df <- max(edf - 1, 1)
  na$flat_p <- stats::pchisq(max(dev_null - dev, 0), df, lower.tail = FALSE)
  na
}

# Cubic B-spline basis on [0, 1] with k functions (internal; Eilers and
# Marx 1996, differences of truncated powers)
bspline_basis <- function(x, k = 10, deg = 3) {
  ndx <- k - deg
  dx <- 1 / ndx
  knots <- seq(-deg * dx, 1 + deg * dx, by = dx)
  tp <- outer(x, knots, function(x, t) (x - t)^deg * (x >= t))
  n <- length(knots)
  d <- diff(diag(n), differences = deg + 1) / (gamma(deg + 1) * dx^deg)
  b <- (-1)^(deg + 1) * tp %*% t(d)
  b[, seq_len(k), drop = FALSE]
}

# Penalised logistic regression spline (P-spline), smoothing parameter
# chosen by UBRE over a grid (internal)
pspline_logit <- function(x, y, k = 10) {
  b <- bspline_basis(x, k)
  dd <- diff(diag(k), differences = 2)
  pen <- crossprod(dd)
  fit <- function(lam) {
    eta <- rep(stats::qlogis(mean(y)), length(y))
    for (it in 1:100) {
      mu <- stats::plogis(eta)
      w <- pmax(mu * (1 - mu), 1e-10)
      z <- eta + (y - mu) / w
      beta <- solve(crossprod(b, b * w) + lam * pen, crossprod(b, w * z))
      new <- drop(b %*% beta)
      done <- max(abs(new - eta)) < 1e-8
      eta <- new
      if (done) break
    }
    mu <- pmin(pmax(stats::plogis(eta), 1e-12), 1 - 1e-12)
    w <- mu * (1 - mu)
    bwb <- crossprod(b, b * w)
    edf <- sum(diag(solve(bwb + lam * pen, bwb)))
    dev <- -2 * sum(y * log(mu) + (1 - y) * log(1 - mu))
    list(beta = beta, deviance = dev, edf = edf,
         ubre = dev / length(y) + 2 * edf / length(y) - 1)
  }
  best <- NULL
  for (l in 10^seq(-2, 8, by = 0.5)) {
    f <- tryCatch(fit(l), error = function(e) NULL)
    if (!is.null(f) && is.finite(f$ubre) &&
        (is.null(best) || f$ubre < best$ubre)) best <- f
  }
  best
}

#' @export
print.boyce_index <- function(x, digits = 3, ...) {
  est <- switch(x$smooth, none = "moving window",
                spline = "penalised regression spline",
                gam = "thin plate regression spline (mgcv)",
                loess = "loess")
  f <- function(v) formatC(v, digits = digits, format = "f")
  cat("Boyce index: ", f(x$index), sep = "")
  if (!is.na(x$lower)) {
    cat(" (", round(100 * x$conf), "% bootstrap interval ", f(x$lower),
        " to ", f(x$upper), ")", sep = "")
  }
  cat("\nP/E curve: ", est, "; ", x$n_presence, " presences, ",
      x$n_background, " background points\n", sep = "")
  cat("\nIN SHORT\n")
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  flat <- !is.na(x$flat_p) && x$flat_p > 0.05
  zero_in <- !is.na(x$lower) && x$lower <= 0 && x$upper >= 0
  if (is.na(x$index)) {
    say("The index cannot be computed (too few points or no range of predictions).")
  } else if (flat) {
    say(paste0("The P/E curve is not distinguishable from flat: presences ",
               "are about as frequent at every predicted value, so the ",
               "model is not shown to be better than random. A smoothed ",
               "index of a flat curve can be anywhere between -1 and 1; ",
               "do not read ", f(x$index), " as a result."))
  } else if (x$index < 0 && !zero_in) {
    say("Presences become rarer as the predicted suitability increases: the model runs against the data.")
  } else if (zero_in) {
    say("The interval includes 0: the model is not shown to be better than random.")
  } else {
    say(if (x$index >= 0.9) {
      "Presences become steadily more frequent as the predicted suitability increases."
    } else if (x$index >= 0.5) {
      "Presences become more frequent as the predicted suitability increases, with some reversals."
    } else {
      "Presences show only a weak tendency to become more frequent as the predicted suitability increases."
    })
  }
  if (!is.na(x$lower) && x$upper - x$lower > 0.5) {
    say("The bootstrap interval is wide: the index is uncertain with these data.")
  }
  if (x$smooth == "none" && x$n_presence < 50) {
    say(paste0("Fewer than 50 presences: the moving-window index tends to ",
               "be underestimated (Liu et al. 2024); a smoothed estimate ",
               "(smooth = \"spline\") is more stable."))
  }
  cat("", strwrap(paste0(
    "The P/E ratio is the share of presences divided by the share of ",
    "background points at a predicted value; the index is its Spearman ",
    "correlation with the prediction (1 = P/E rises steadily, 0 = random). ",
    "The wording \"steadily\" (index >= 0.9), \"with some reversals\" ",
    "(>= 0.5) and \"weak\" (below 0.5) is a rule of thumb.", if (!is.na(x$flat_p)) paste0(
      " Flat curve: likelihood-ratio test of the smooth against a constant, ",
      "p = ", formatC(x$flat_p, digits = 2, format = "g"), ".")),
    width = 78), sep = "\n")
  invisible(x)
}

#' Plot the P/E curve of the Boyce index
#'
#' The predicted-to-expected (P/E) ratio against the predicted value, as
#' estimated by [boyce_index()]: the moving-window ratios, or the smoothed
#' curve with its bootstrap band (on a log scale, as the smoothed ratio
#' spans orders of magnitude). A ratio above 1 (dotted line) means presences
#' are more frequent there than in the background.
#'
#' @param x A result of [boyce_index()] with `details = TRUE` or
#'   `n_boot > 0`.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Liu, C., Newell, G., White, M., Machunter, J. (2024)
#'   Improving the estimation of the Boyce index using statistical smoothing
#'   methods for evaluating species distribution models with presence-only
#'   data. *Ecography* 2025, e07218. \doi{10.1111/ecog.07218}
#' @export
#' @examples
#' set.seed(1)
#' b <- boyce_index(rbeta(30, 4, 1), runif(1000), smooth = "spline",
#'                  details = TRUE)
#' plot_boyce(b)
plot_boyce <- function(x, base_size = 11, file = NULL, width = 7,
                       height = 5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "boyce_index")) {
    stop("`x` must come from boyce_index(..., details = TRUE).",
         call. = FALSE)
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  d <- x$curve
  sub <- paste0(
    if (x$smooth == "none") "Moving window" else "Smoothed (" ,
    if (x$smooth != "none") paste0(x$smooth, ")"),
    "; ", x$n_presence, " presences, ", x$n_background,
    " background points.\nBoyce index ", sprintf("%.2f", x$index),
    if (!is.na(x$lower)) sprintf(" (%d%% interval %.2f to %.2f)",
                                 round(100 * x$conf), x$lower, x$upper),
    if ("pe_lower" %in% names(d)) "; band: bootstrap interval of the curve",
    ".")
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$prediction, .data$pe)) +
    ggplot2::geom_hline(yintercept = 1, colour = muted, linetype = "dotted")
  if ("pe_lower" %in% names(d)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$pe_lower,
                                               ymax = .data$pe_upper),
                                  fill = "#9ecae1", alpha = 0.6)
  }
  p <- p + ggplot2::geom_line(colour = "#08519c", linewidth = 0.7) +
    ggplot2::labs(x = "Predicted value", y = "P/E ratio",
                  title = "Predicted-to-expected ratio", subtitle = sub) +
    theme_enm(base_size, ink, muted)
  # a smoothed P/E ratio is positive and spans orders of magnitude
  if (x$smooth != "none") {
    top <- max(c(d$pe, d$pe_upper), na.rm = TRUE)
    p <- p + ggplot2::scale_y_log10(
      breaks = 10^(-2:4), labels = format(10^(-2:4), scientific = FALSE,
                                          drop0trailing = TRUE, trim = TRUE)) +
      ggplot2::coord_cartesian(ylim = c(0.01, max(top, 2))) +
      ggplot2::labs(y = "P/E ratio (log scale, cut at 0.01)")
  }
  p
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

# Calibration ----------------------------------------------------------------------

#' Brier score
#'
#' Mean squared difference between the predicted probability and the
#' observation (1 = presence, 0 = absence): `mean((pred - obs)^2)`, from 0
#' (perfect) to 1. This is the common binary form, half of Brier's original
#' definition summed over both classes (Silva Filho et al. 2023).
#'
#' @inheritParams auc_score
#' @return A single number.
#' @references Silva Filho, T., Song, H., Perello-Nieto, M.,
#'   Santos-Rodriguez, R., Kull, M., Flach, P. (2023) Classifier calibration:
#'   a survey on how to assess and improve predicted class probabilities.
#'   *Machine Learning* 112, 3211-3260. \doi{10.1007/s10994-023-06336-7}
#' @export
#' @examples
#' brier_score(c(1, 1, 0, 0), c(0.9, 0.6, 0.7, 0.1))
brier_score <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  mean((pred[ok] - obs[ok])^2)
}

# Logistic calibration (internal): slope of logit P(obs) on logit(pred), and
# the calibration intercept with the slope fixed at 1 (offset), each with its
# standard error. Predictions are clamped to [eps, 1 - eps].
calibration_fit <- function(obs, pred, eps = 1e-6) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  na <- c(slope = NA_real_, slope_se = NA_real_, intercept = NA_real_,
          intercept_se = NA_real_)
  if (length(obs) < 3 || length(unique(obs)) < 2) return(na)
  if (any(pred < 0 | pred > 1)) return(na)
  p <- pmin(pmax(pred, eps), 1 - eps)
  lp <- stats::qlogis(p)
  fit1 <- tryCatch(suppressWarnings(stats::glm.fit(
    cbind(1, lp), obs, family = stats::binomial())), error = function(e) NULL)
  fit0 <- tryCatch(suppressWarnings(stats::glm.fit(
    matrix(1, length(obs)), obs, offset = lp, family = stats::binomial())),
    error = function(e) NULL)
  se <- function(f, i) {
    if (is.null(f)) return(NA_real_)
    v <- tryCatch(chol2inv(f$qr$qr[seq_len(f$rank), seq_len(f$rank),
                                    drop = FALSE]), error = function(e) NULL)
    if (is.null(v)) NA_real_ else sqrt(v[i, i])
  }
  c(slope = if (is.null(fit1)) NA_real_ else unname(fit1$coefficients[2]),
    slope_se = se(fit1, 2),
    intercept = if (is.null(fit0)) NA_real_ else
      unname(fit0$coefficients[1]),
    intercept_se = se(fit0, 1))
}

calibration_slope_score <- function(obs, pred) {
  unname(calibration_fit(obs, pred)["slope"])
}
calibration_intercept_score <- function(obs, pred) {
  unname(calibration_fit(obs, pred)["intercept"])
}

# Bins of predictions (internal): equal width over the range, or equal count
calibration_bins <- function(pred, bins, binning) {
  if (binning == "width") {
    rg <- range(pred)
    if (diff(rg) == 0) return(rep(1L, length(pred)))
    br <- seq(rg[1], rg[2], length.out = bins + 1)
    findInterval(pred, br, rightmost.closed = TRUE, all.inside = TRUE)
  } else {
    r <- rank(pred, ties.method = "first")
    pmin(bins, as.integer(ceiling(r / length(pred) * bins)))
  }
}

#' Calibration curve (reliability diagram) of a prediction
#'
#' Calibration asks whether predicted probabilities can be read as
#' probabilities: among places predicted at 0.3, is the species found at
#' about 30%? Discrimination (AUC, TSS, Boyce) cannot tell, because it only
#' ranks the predictions; this is why MAE and balanced MAE are used next to
#' AUC (Konowalik and Nosol 2021). The machine-learning literature has
#' formalised the same idea into a standard toolkit (Silva Filho et al.
#' 2023), which this function implements:
#'
#' * **Reliability diagram.** The predictions are cut into `bins`, of equal
#'   width over their range or with equal numbers of points
#'   (`binning = "count"`). For each bin, the observed frequency of
#'   presences is plotted against the **mean prediction** of the bin (not
#'   the bin centre, which can make a well calibrated model look
#'   miscalibrated by up to half a bin), with a binomial (Clopper-Pearson)
#'   confidence interval. Points on the 1:1 line mean good calibration.
#' * **Calibration slope**: the slope of a logistic regression of the
#'   observations on the logit of the predictions. 1 = calibrated; below 1
#'   the predictions are too extreme (over-confident: high ones too high,
#'   low ones too low); above 1 too moderate (under-confident).
#' * **Calibration intercept** (calibration-in-the-large): the intercept of
#'   the same regression with the slope fixed at 1. 0 = calibrated; above 0
#'   the predictions are too low on average, below 0 too high.
#'
#' **Presence-background data.** With background points instead of true
#' absences, the "observed frequency" is the share of presences among the
#' presences and background points of a bin. It depends on how many
#' background points there are and on how well they represent the
#' environment of the study area (the environmental representativeness
#' effect; Jimenez-Valverde 2025, 2026). Then these quantities compare
#' models evaluated on the same presences and background; they do not
#' measure calibration in absolute terms. Predictions that are not
#' probabilities (e.g. raw Maxent output) cannot be calibrated.
#'
#' @param obs Observations (1 = presence, 0 = absence/background), or a
#'   table from [extract_predictions()] (a data.frame with `obs` and one
#'   column per model): then `pred` is not needed and every model gets its
#'   curve.
#' @param pred Predicted probabilities (between 0 and 1).
#' @param bins Number of bins.
#' @param binning `"width"` (equal-width bins over the range of the
#'   predictions) or `"count"` (equal numbers of points).
#' @param conf Confidence level of the intervals.
#' @return A list of class `calibration_curve` with `bins` (one row per
#'   model and bin: `n`, `mean_pred`, `observed`, `lower`, `upper`) and
#'   `summary` (one row per model: `n`, `slope`, `slope_lo`, `slope_hi`,
#'   `intercept`, `intercept_lo`, `intercept_hi`, `brier`, and the Brier
#'   decomposition of [brier_decompose()]).
#' @references Silva Filho, T., Song, H., Perello-Nieto, M.,
#'   Santos-Rodriguez, R., Kull, M., Flach, P. (2023) Classifier calibration:
#'   a survey on how to assess and improve predicted class probabilities.
#'   *Machine Learning* 112, 3211-3260. \doi{10.1007/s10994-023-06336-7}
#'
#'   Konowalik, K., Nosol, A. (2021) Evaluation metrics and validation of
#'   presence-only species distribution models based on distributional maps
#'   with varying coverage. *Scientific Reports* 11, 1482.
#'   \doi{10.1038/s41598-020-80062-1}
#'
#'   Jimenez-Valverde, A. (2025) Refining uniform discrimination metrics:
#'   towards a case-by-case weighting evaluation in species distribution
#'   models with presence-absence data. *Ecology and Evolution* 15, e72573.
#'   \doi{10.1002/ece3.72573}
#'
#'   Jimenez-Valverde, A. (2026) The environmental representativeness effect
#'   in species distribution model evaluation without absence data.
#'   *Biodiversity Informatics* 20, 74-80. \doi{10.17161/bi.v20i1.24860}
#' @seealso [plot_calibration()], [brier_decompose()], [evaluate_models()]
#' @export
#' @examples
#' set.seed(1)
#' p <- runif(2000)
#' y <- rbinom(2000, 1, p)                  # calibrated
#' calibration_curve(y, p)
#' # over-confident: the same ranking, more extreme probabilities
#' p_over <- plogis(2 * qlogis(p))
#' cc <- calibration_curve(y, p_over)
#' cc
#' plot_calibration(cc)
calibration_curve <- function(obs, pred = NULL, bins = 10,
                              binning = c("width", "count"), conf = 0.95) {
  binning <- match.arg(binning)
  if (!is.numeric(bins) || bins < 2) {
    stop("`bins` must be a number of at least 2.", call. = FALSE)
  }
  if (is.data.frame(obs)) {
    validate_table(obs)
    tab <- obs
    mods <- setdiff(names(tab), "obs")
    parts <- lapply(mods, function(m) {
      calibration_one(tab$obs, tab[[m]], bins, binning, conf, m)
    })
  } else {
    if (is.null(pred)) stop("Give `pred`.", call. = FALSE)
    parts <- list(calibration_one(obs, pred, bins, binning, conf, NULL))
  }
  structure(list(bins = do.call(rbind, lapply(parts, `[[`, "bins")),
                 summary = do.call(rbind, lapply(parts, `[[`, "summary")),
                 settings = list(bins = bins, binning = binning,
                                 conf = conf)),
            class = "calibration_curve")
}

# One model (internal)
calibration_one <- function(obs, pred, bins, binning, conf, model) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (!all(obs %in% c(0, 1))) stop("`obs` must be 0 or 1.", call. = FALSE)
  if (any(pred < 0 | pred > 1)) {
    stop("Predictions", if (!is.null(model)) paste0(" of ", model),
         " are not probabilities (outside 0-1): they cannot be calibrated.",
         call. = FALSE)
  }
  b <- calibration_bins(pred, bins, binning)
  ci <- function(k, n) {
    if (n == 0) return(c(NA_real_, NA_real_))
    stats::binom.test(k, n, conf.level = conf)$conf.int[1:2]
  }
  tb <- do.call(rbind, lapply(sort(unique(b)), function(i) {
    s <- b == i
    k <- sum(obs[s])
    n <- sum(s)
    iv <- ci(k, n)
    data.frame(bin = i, n = n, mean_pred = mean(pred[s]),
               observed = k / n, lower = iv[1], upper = iv[2])
  }))
  if (!is.null(model)) tb <- cbind(model = model, tb,
                                   stringsAsFactors = FALSE)
  cf <- calibration_fit(obs, pred)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  bd <- brier_decompose(obs, pred, bins = bins, binning = binning)
  sm <- data.frame(n = length(obs), slope = cf[["slope"]],
                   slope_lo = cf[["slope"]] - z * cf[["slope_se"]],
                   slope_hi = cf[["slope"]] + z * cf[["slope_se"]],
                   intercept = cf[["intercept"]],
                   intercept_lo = cf[["intercept"]] - z * cf[["intercept_se"]],
                   intercept_hi = cf[["intercept"]] + z * cf[["intercept_se"]],
                   as.list(bd[c("brier", "reliability", "resolution",
                                "uncertainty")]))
  if (!is.null(model)) sm <- cbind(model = model, sm,
                                   stringsAsFactors = FALSE)
  list(bins = tb, summary = sm)
}

# Plain words for one model (internal). Verdict first, numbers at the end.
calibration_words <- function(s, conf, name = NULL) {
  f <- function(v) sprintf("%.2f", v)
  pc <- paste0(round(100 * conf), "%")
  iv <- function(lo, hi) paste0(" (", pc, " interval ", f(lo), " to ",
                                f(hi), ")")
  numbers <- if (is.na(s$slope)) "" else paste0(
    " Slope ", f(s$slope), iv(s$slope_lo, s$slope_hi), ", intercept ",
    f(s$intercept), iv(s$intercept_lo, s$intercept_hi), ".")
  who <- if (is.null(name)) "The model" else name
  if (is.na(s$slope)) {
    return(paste0(who, ": calibration cannot be estimated (too few points ",
                  "or only one class)."))
  }
  if (s$slope_hi <= 0 || (s$slope_lo <= 0 && s$slope_hi >= 0)) {
    verdict <- if (s$slope_hi < 0) {
      "the predictions run against the observations (negative slope), so calibration cannot be judged"
    } else {
      "the predictions do not separate presences from absences (slope not different from 0), so calibration cannot be judged"
    }
    return(paste0(who, ": ", verdict, ".", numbers))
  }
  slope_ok <- s$slope_lo <= 1 && s$slope_hi >= 1
  int_ok <- s$intercept_lo <= 0 && s$intercept_hi >= 0
  slight <- s$slope >= 0.8 && s$slope <= 1.25
  spread <- if (slope_ok) {
    if (int_ok) "consistent with calibration: the spread" else "the spread"
  } else if (s$slope < 1) {
    paste0(if (slight) "slightly " else "", "over-confident: the ",
           "predictions are too extreme, high ones too high and low ones ",
           "too low")
  } else {
    paste0(if (slight) "slightly " else "", "under-confident: the ",
           "predictions are too moderate, too close to the middle")
  }
  if (slope_ok) spread <- paste(spread, "of the predictions fits the",
                                "observations")
  level <- if (int_ok) "on average they fit" else if (s$intercept > 0) {
    "on average they are too low"
  } else "on average they are too high"
  txt <- paste0(spread, if (slope_ok && int_ok) " and " else "; ", level,
                ".", numbers)
  if (!is.null(name)) return(paste0(name, ": ", txt))
  paste0(toupper(substr(txt, 1, 1)), substr(txt, 2, nchar(txt)))
}

#' @export
print.calibration_curve <- function(x, digits = 3, ...) {
  s <- x$summary
  multi <- "model" %in% names(s)
  cat("Calibration of the predictions (", x$settings$bins, " bins of equal ",
      if (x$settings$binning == "width") "width" else "count", ")\n\n",
      sep = "")
  show <- s[intersect(c("model", "n", "slope", "slope_lo", "slope_hi",
                        "intercept", "brier", "reliability", "resolution"),
                      names(s))]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("\nIN SHORT\n")
  for (i in seq_len(nrow(s))) {
    txt <- calibration_words(s[i, ], x$settings$conf,
                             if (multi) s$model[i] else NULL)
    cat(strwrap(txt, width = 76, exdent = 2, initial = "- "), sep = "\n")
  }
  if (nrow(s) && any(s$n < 100)) {
    cat(strwrap(paste0(
      "Fewer than 100 points", if (multi) paste0(" for ", paste(
        s$model[s$n < 100], collapse = ", ")), ": the bins are small ",
      "and the verdict is uncertain (wide intervals)."), width = 76,
      initial = "- ", exdent = 2), sep = "\n")
  }
  cat("", strwrap(paste0(
    "With background points instead of absences, compare only models ",
    "evaluated on the same presences and background: the values depend ",
    "on the number and environmental representativeness of the background ",
    "points (Jimenez-Valverde 2025, 2026)."), width = 78), sep = "\n")
  cat("", strwrap(paste0(
    "Definitions: the calibration slope is the slope of a logistic ",
    "regression of the observations on the logit of the predictions (1 = ",
    "calibrated, below 1 over-confident, above 1 under-confident); the ",
    "intercept is that of the same regression with the slope fixed at 1 ",
    "(0 = calibrated, above 0 predictions too low on average). ",
    "\"Slightly\" marks slopes between 0.8 and 1.25, a rule of thumb. ",
    "Brier = reliability - resolution + uncertainty (+ two small ",
    "within-bin terms, see brier_decompose())."), width = 78), sep = "\n")
  invisible(x)
}

#' Brier score decomposition
#'
#' Splits the Brier score into reliability (calibration loss: how far the
#' predictions are from the observed frequencies), resolution (how much the
#' observed frequencies of the bins differ from the overall frequency: the
#' information in the predictions) and uncertainty (the variance of the
#' observations, the same for every model):
#' `Brier = reliability - resolution + uncertainty` (Murphy 1973). With
#' predictions grouped into bins this holds only approximately, so two
#' within-bin terms are added (Stephenson et al. 2008), and
#' `Brier = reliability - resolution + uncertainty + within_bin_variance -
#' within_bin_covariance` exactly. Lower reliability and higher resolution
#' are better. See [calibration_curve()] for the presence-background caveat.
#'
#' @inheritParams calibration_curve
#' @return A named numeric vector: `brier`, `reliability`, `resolution`,
#'   `uncertainty`, `within_bin_variance`, `within_bin_covariance`.
#' @references Murphy, A.H. (1973) A new vector partition of the probability
#'   score. *Journal of Applied Meteorology* 12, 595-600.
#'
#'   Stephenson, D.B., Coelho, C.A.S., Jolliffe, I.T. (2008) Two extra
#'   components in the Brier score decomposition. *Weather and Forecasting*
#'   23, 752-757.
#'
#'   Silva Filho, T., Song, H., Perello-Nieto, M., Santos-Rodriguez, R.,
#'   Kull, M., Flach, P. (2023) Classifier calibration: a survey on how to
#'   assess and improve predicted class probabilities. *Machine Learning*
#'   112, 3211-3260. \doi{10.1007/s10994-023-06336-7}
#' @export
#' @examples
#' set.seed(1)
#' p <- runif(1000)
#' y <- rbinom(1000, 1, p)
#' bd <- brier_decompose(y, p)
#' bd
#' sum(bd[c("reliability", "uncertainty", "within_bin_variance")]) -
#'   sum(bd[c("resolution", "within_bin_covariance")])   # = Brier
brier_decompose <- function(obs, pred, bins = 10,
                            binning = c("width", "count")) {
  binning <- match.arg(binning)
  ok <- !is.na(obs) & !is.na(pred)
  o <- obs[ok]
  f <- pred[ok]
  n <- length(o)
  b <- calibration_bins(f, bins, binning)
  ob <- mean(o)
  fk <- stats::ave(f, b)
  ok_bin <- stats::ave(o, b)
  rel <- sum((fk - ok_bin)^2) / n
  res <- sum((ok_bin - ob)^2) / n
  unc <- ob * (1 - ob)
  wbv <- sum((f - fk)^2) / n
  wbc <- 2 * sum((o - ok_bin) * (f - fk)) / n
  c(brier = mean((f - o)^2), reliability = rel, resolution = res,
    uncertainty = unc, within_bin_variance = wbv,
    within_bin_covariance = wbc)
}

#' Reliability diagram
#'
#' The observed frequency of presences in each bin against the mean
#' prediction of the bin, with binomial confidence intervals (red bars) and
#' the 1:1 line of perfect calibration (dashed); point size shows the
#' number of points in the bin. One panel per model, labelled with the
#' calibration slope and intercept.
#'
#' @param x A result of [calibration_curve()].
#' @param base_size Base font size.
#' @param bar_colour Colour of the confidence intervals.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Silva Filho, T., Song, H., Perello-Nieto, M.,
#'   Santos-Rodriguez, R., Kull, M., Flach, P. (2023) Classifier calibration:
#'   a survey on how to assess and improve predicted class probabilities.
#'   *Machine Learning* 112, 3211-3260. \doi{10.1007/s10994-023-06336-7}
#' @export
#' @examples
#' set.seed(1)
#' p <- runif(2000)
#' plot_calibration(calibration_curve(rbinom(2000, 1, p), p))
plot_calibration <- function(x, base_size = 11, bar_colour = "#c8553d",
                             file = NULL, width = 7, height = 5.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
  if (!inherits(x, "calibration_curve")) {
    stop("`x` must come from calibration_curve().", call. = FALSE)
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  d <- x$bins
  s <- x$summary
  multi <- "model" %in% names(d)
  lab <- sprintf("slope %.2f, intercept %.2f", s$slope, s$intercept)
  if (multi) {
    d$panel <- factor(paste0(d$model, "\n",
                             lab[match(d$model, s$model)]),
                      paste0(s$model, "\n", lab))
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$mean_pred, .data$observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         colour = muted) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$lower,
                                        ymax = .data$upper),
                           width = 0, colour = bar_colour, linewidth = 0.5) +
    ggplot2::geom_line(colour = "#08519c", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = .data$n), colour = "#08519c") +
    ggplot2::scale_size_area(max_size = 4, name = "Points\nin the bin") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(x = "Mean prediction in the bin",
                  y = "Observed frequency of presences",
                  title = "Reliability diagram",
                  subtitle = paste0(
                    "Dashed: perfect calibration; bars: ",
                    round(100 * x$settings$conf), "% binomial intervals.",
                    if (!multi) paste0("\nCalibration ", lab[1], "."))) +
    theme_enm(base_size, ink, muted)
  if (multi) p <- p + ggplot2::facet_wrap(~panel)
  p
}
