test_that("AUC matches the Mann-Whitney definition, incl. ties", {
  obs <- c(1, 1, 1, 0, 0, 0, 0)
  pred <- c(0.9, 0.5, 0.5, 0.5, 0.2, 0.1, 0.6)
  # brute-force: P(pred_pres > pred_abs) + 0.5 * P(tie)
  pp <- pred[obs == 1]; pa <- pred[obs == 0]
  cmp <- outer(pp, pa, ">") + 0.5 * outer(pp, pa, "==")
  expect_equal(auc_score(obs, pred), mean(cmp))
  expect_equal(auc_score(c(1, 0), c(1, 0)), 1)
  expect_true(is.na(auc_score(c(1, 1), c(0.2, 0.3))))
})

test_that("MAE and bias match the original Metrics-package formulas", {
  obs <- c(1, 0, 1, 0)
  pred <- c(0.8, 0.3, 0.4, 0.1)
  expect_equal(mae_score(obs, pred), mean(abs(obs - pred)))
  expect_equal(bias_score(obs, pred), mean(obs - pred))
  if (requireNamespace("Metrics", quietly = TRUE)) {
    expect_equal(mae_score(obs, pred), Metrics::mae(obs, pred))
    expect_equal(bias_score(obs, pred), Metrics::bias(obs, pred))
    expect_equal(auc_score(obs, pred), Metrics::auc(obs, pred))
  }
})

test_that("TSS equals the brute-force maximum over thresholds", {
  set.seed(4)
  obs <- rbinom(300, 1, 0.3)
  pred <- round(plogis(obs * 1.5 + rnorm(300)), 2)   # rounded: many ties
  brute <- max(vapply(c(sort(unique(pred)), Inf), function(t) {
    mean(pred[obs == 1] >= t) + mean(pred[obs == 0] < t) - 1
  }, numeric(1)))
  expect_equal(as.numeric(tss_score(obs, pred)), brute)
  expect_equal(as.numeric(tss_score(c(1, 1, 0, 0), c(0.9, 0.8, 0.2, 0.1))), 1)
  expect_equal(as.numeric(tss_score(c(1, 0, 1, 0), c(0.9, 0.2, 0.3, 0.6),
                                    threshold = 0.5)), 0)
  expect_true(is.na(tss_score(c(1, 1), c(0.5, 0.6))))
})

test_that("Boyce index reproduces ecospat::ecospat.boyce()", {
  # reference values computed with ecospat 4.1.4 (rounded to 3 decimals)
  set.seed(10); bg <- runif(500); pr <- rbeta(60, 4, 1)
  expect_equal(round(boyce_index(pr, bg), 3), 0.961)
  set.seed(11); bg <- rbeta(300, 1, 3); pr <- runif(25)
  expect_equal(round(boyce_index(pr, bg), 3), 0.743)
  set.seed(12); bg <- round(runif(200), 2); pr <- round(rbeta(40, 2, 2), 2)
  expect_equal(round(boyce_index(pr, bg), 3), 0.268)
})

test_that("Boyce index handles degenerate input", {
  expect_true(is.na(boyce_index(0.5, runif(10))))           # 1 presence
  expect_true(is.na(boyce_index(rep(0.5, 5), rep(0.5, 50)))) # no range
  expect_true(is.na(boyce_index(c(NA, NA), runif(10))))
})

test_that("threshold_p10 is the 10th percentile", {
  expect_equal(threshold_p10(1:11 / 11), unname(quantile(1:11 / 11, 0.1)))
})

test_that("MAEbal weighs presences and background equally", {
  obs <- c(1, 1, rep(0, 98))
  narrow <- c(0.6, 0.6, rep(0, 98))
  expect_equal(mae_balanced(obs, narrow), 0.5 * 0.4 + 0)
  expect_lt(mae_score(obs, narrow), mae_balanced(obs, narrow))
  expect_true(is.na(mae_balanced(c(0, 0), c(0.1, 0.2))))
  # independent of the number of background points
  expect_equal(mae_balanced(c(1, 0), c(0.8, 0.3)),
               mae_balanced(c(1, rep(0, 50)), c(0.8, rep(0.3, 50))))
})

test_that("a user metric may replace a built-in one", {
  sim <- simulate_enm_data(seed = 1)
  my <- function(obs, pred) 42
  expect_message(r <- evaluate_models(sim$models, sim$presence, sim$absence,
                                      metrics = list(MAEbal = my),
                                      criteria = c(AUC = "max",
                                                   MAEbal = "min")),
                 "replace the built-in")
  expect_true(all(r$MAEbal == 42))
  expect_error(evaluate_models(sim$models, sim$presence, sim$absence,
                               metrics = list(obs = my)), "reserved")
})

test_that("boyce_index equals the window-by-window definition", {
  slow <- function(presence, background, n_windows = 100, window_width = 0.1) {
    lo <- min(background, presence)
    hi <- max(background, presence)
    w <- (max(background) - min(background)) * window_width
    start <- seq(from = lo, to = hi - w, by = (hi - lo - w) / n_windows)
    start[length(start)] <- start[length(start)] + 1
    end <- start + w
    ratio <- vapply(seq_along(start), function(i) {
      round(mean(presence >= start[i] & presence <= end[i]) /
              mean(background >= start[i] & background <= end[i]), 10)
    }, numeric(1))
    keep <- which(!is.nan(ratio))
    ratio <- ratio[keep]
    start <- start[keep]
    r <- which(ratio != c(ratio[-1], TRUE))
    suppressWarnings(stats::cor(ratio[r], start[r], method = "spearman"))
  }
  set.seed(7)
  for (i in 1:5) {
    bg <- round(runif(800), 2)            # many ties and window-edge values
    pres <- round(rbeta(60, 3, 1.5), 2)
    expect_identical(boyce_index(pres, bg), slow(pres, bg))
    expect_identical(boyce_index(pres, bg, 20, 0.25), slow(pres, bg, 20, 0.25))
  }
})

test_that("thresholds: every rule inside the range of the predictions", {
  set.seed(1)
  pres <- stats::rbeta(150, 5, 2)
  bg <- stats::rbeta(1500, 2, 4)
  th <- thresholds(pres, bg)
  expect_s3_class(th, "enm_thresholds")
  expect_equal(th$rule, c("maxSSS", "maxTSS", "equal_ss", "maxMCC", "maxF",
                          "prevalence", "fixed_sens", "p10",
                          "mutation_point"))
  rg <- range(c(pres, bg))
  expect_true(all(th$threshold >= rg[1] & th$threshold <= rg[2]))
  # definitions
  r <- function(x) th[th$rule == x, ]
  expect_equal(r("p10")$threshold, threshold_p10(pres))
  expect_gte(r("fixed_sens")$sensitivity, 0.95)
  expect_equal(r("maxSSS")$threshold, r("maxTSS")$threshold)
  expect_equal(r("maxTSS")$TSS, max(vapply(sort(unique(c(pres, bg))),
    function(t) mean(pres >= t) + mean(bg < t) - 1, numeric(1))))
  expect_lt(abs(r("equal_ss")$sensitivity - r("equal_ss")$specificity), 0.02)
  expect_equal(th$suitable_share, 1 - th$specificity)
  # a given prevalence: that share of the background is predicted suitable
  tp <- thresholds(pres, bg, rules = "prevalence", prevalence = 0.2)
  expect_equal(tp$suitable_share, 0.2, tolerance = 0.002)
  # MCC and F from the confusion matrix, checked by hand
  mcc <- function(t) {
    # as doubles: the product of four counts overflows integers
    tp <- as.numeric(sum(pres >= t)); fn <- as.numeric(sum(pres < t))
    tn <- as.numeric(sum(bg < t)); fp <- as.numeric(sum(bg >= t))
    (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  }
  f1 <- function(t) {
    tp <- sum(pres >= t)
    2 * tp / (sum(pres < t) + 2 * tp + sum(bg >= t))
  }
  cand <- sort(unique(c(pres, bg)))
  expect_equal(r("maxMCC")$MCC, max(vapply(cand, mcc, numeric(1)),
                                    na.rm = TRUE))
  expect_equal(r("maxF")$F, max(vapply(cand, f1, numeric(1))))
  expect_equal(r("p10")$MCC, mcc(r("p10")$threshold))
  # MCC and F penalise overprediction: smaller ranges than TSS (Hellegers
  # et al. 2025)
  expect_lt(r("maxMCC")$suitable_share, r("maxTSS")$suitable_share)
  expect_lt(r("maxF")$suitable_share, r("maxTSS")$suitable_share)
  expect_output(print(th), "times as much")
})

test_that("maxSSS and maxTSS agree on separable data", {
  pres <- c(0.7, 0.8, 0.9, 0.95)
  bg <- c(0.05, 0.1, 0.2, 0.3, 0.4)
  th <- thresholds(pres, bg, rules = c("maxSSS", "maxTSS", "equal_ss",
                                       "maxMCC", "maxF"))
  expect_equal(th$threshold, rep(0.7, 5))      # lowest perfect threshold
  expect_equal(th$TSS, rep(1, 5))
  expect_equal(th$MCC, rep(1, 5))
  expect_equal(th$F, rep(1, 5))
})

test_that("the mutation point finds the break of the P/E curve", {
  set.seed(2)
  bg <- stats::runif(20000)
  # presences: as frequent as the background below 0.6, ten times above
  p <- ifelse(bg < 0.6, 1, 10)
  pres <- sample(bg, 3000, prob = p, replace = TRUE)
  mp <- thresholds(pres, bg, rules = "mutation_point")$threshold
  expect_gt(mp, 0.5)
  expect_lt(mp, 0.7)
})

test_that("rules subset; tables of several models; errors", {
  sim <- simulate_enm_data(seed = 1)
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  th <- thresholds(tab, rules = c("p10", "maxSSS"))
  expect_equal(unique(th$rule), c("p10", "maxSSS"))
  expect_equal(sort(unique(th$model)), sort(names(sim$models)))
  expect_equal(nrow(th), 2 * terra::nlyr(sim$models))
  expect_error(thresholds(1:3, 1:3, rules = "nope"), "Unknown rule")
  expect_error(thresholds(1:3), "background")
  expect_error(thresholds(1:3, 1:3, sensitivity = 1.2), "sensitivity")
  expect_s3_class(plot_thresholds(th), "ggplot")
})

test_that("threshold_agreement counts the rules per cell", {
  r <- terra::rast(nrows = 2, ncols = 3, xmin = 0, xmax = 3, ymin = 0,
                   ymax = 2, crs = "EPSG:4326",
                   vals = c(0.1, 0.35, 0.5, 0.65, 0.9, NA))
  th <- structure(data.frame(rule = c("a", "b", "c"),
                             threshold = c(0.3, 0.6, 0.8),
                             sensitivity = NA, specificity = NA, TSS = NA,
                             suitable_share = NA),
                  class = c("enm_thresholds", "data.frame"))
  ag <- threshold_agreement(r, th)
  expect_equal(as.vector(terra::values(ag$count)), c(0, 1, 1, 2, 3, NA))
  s <- ag$summary
  expect_equal(s$cells, c(1, 2, 4))            # all, most (> 1.5), any
  expect_true(all(diff(s$area_km2) >= 0))
  expect_output(print(ag), "Suitable under every rule")
  # duplicated thresholds count once
  th2 <- th
  th2$threshold[2] <- 0.3
  expect_equal(threshold_agreement(r, th2)$n_rules, 2)
  expect_s3_class(plot_thresholds(ag), "ggplot")
  expect_error(threshold_agreement(r, data.frame()), "thresholds")
})

# Calibration ----------------------------------------------------------------

test_that("a calibrated predictor has slope ~1 and intercept ~0", {
  set.seed(42)
  p <- runif(20000)
  y <- rbinom(20000, 1, p)
  cc <- calibration_curve(y, p)
  s <- cc$summary
  expect_equal(s$slope, 1, tolerance = 0.05)
  expect_equal(s$intercept, 0, tolerance = 0.05)
  expect_true(s$slope_lo <= 1 && s$slope_hi >= 1)
  # observed frequencies follow the mean predictions of the bins
  expect_lt(max(abs(cc$bins$observed - cc$bins$mean_pred)), 0.04)
  expect_equal(sum(cc$bins$n), 20000)
  expect_output(print(cc), "Consistent with calibration")
})

test_that("over- and under-confident predictors are recognised", {
  set.seed(42)
  p <- runif(20000)
  y <- rbinom(20000, 1, p)
  # logit doubled: same ranking, too extreme -> slope 1/2
  over <- calibration_curve(y, plogis(2 * qlogis(p)))
  expect_equal(over$summary$slope, 0.5, tolerance = 0.05)
  expect_output(print(over), "Over-confident")
  # logit halved: too moderate -> slope 2
  under <- calibration_curve(y, plogis(0.5 * qlogis(p)))
  expect_equal(under$summary$slope, 2, tolerance = 0.05)
  expect_output(print(under), "Under-confident")
  # predictions too low on average -> positive intercept
  low <- calibration_curve(y, plogis(qlogis(p) - 1))
  expect_equal(low$summary$intercept, 1, tolerance = 0.05)
  expect_output(print(low), "too low")
  # no information -> no verdict
  rnd <- calibration_curve(y, runif(20000))
  expect_output(print(rnd), "cannot be judged")
})

test_that("binning by count gives equal bins; bins of known data", {
  set.seed(1)
  p <- runif(1000)
  y <- rbinom(1000, 1, p)
  cc <- calibration_curve(y, p, bins = 4, binning = "count")
  expect_equal(cc$bins$n, rep(250, 4))
  # two bins by hand
  obs <- c(0, 0, 1, 0, 1, 1)
  pred <- c(0.1, 0.2, 0.3, 0.7, 0.8, 0.9)
  b <- calibration_curve(obs, pred, bins = 2)$bins
  expect_equal(b$mean_pred, c(0.2, 0.8))
  expect_equal(b$observed, c(1 / 3, 2 / 3))
  ci <- stats::binom.test(1, 3)$conf.int
  expect_equal(c(b$lower[1], b$upper[1]), ci[1:2])
})

test_that("Brier decomposition sums to the Brier score", {
  set.seed(3)
  p <- runif(5000)
  y <- rbinom(5000, 1, plogis(1.7 * qlogis(p) + 0.3))
  for (bn in c("width", "count")) {
    bd <- brier_decompose(y, p, bins = 7, binning = bn)
    expect_equal(unname(bd["brier"]), mean((p - y)^2))
    expect_equal(unname(bd["reliability"] - bd["resolution"] +
                          bd["uncertainty"] + bd["within_bin_variance"] -
                          bd["within_bin_covariance"]),
                 unname(bd["brier"]))
    expect_equal(unname(bd["uncertainty"]), mean(y) * (1 - mean(y)))
  }
  # predictions constant within bins: Murphy's three terms are exact
  q <- rep(c(0.2, 0.5, 0.8), each = 100)
  o <- c(rep(1:0, c(20, 80)), rep(1:0, c(50, 50)), rep(1:0, c(80, 20)))
  bd <- brier_decompose(o, q, bins = 3)
  expect_equal(unname(bd["reliability"]), 0)
  expect_equal(unname(bd["within_bin_variance"]), 0)
  expect_equal(unname(bd["brier"]),
               unname(bd["uncertainty"] - bd["resolution"]))
  expect_equal(brier_score(c(1, 0), c(0.8, 0.4)), (0.04 + 0.16) / 2)
})

test_that("calibration works on tables, as criteria and in plots", {
  sim <- simulate_enm_data(seed = 1)
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  cc <- calibration_curve(tab)
  expect_equal(cc$summary$model, setdiff(names(tab), "obs"))
  expect_s3_class(plot_calibration(cc), "ggplot")
  expect_s3_class(plot_calibration(calibration_curve(tab$obs, tab[[2]])),
                  "ggplot")
  f <- tempfile(fileext = ".png")
  plot_calibration(cc, file = f)
  expect_true(file.exists(f))
  ev <- evaluate_models(data = tab, criteria = c(AUC = "max",
                                                 calibration_slope = "one",
                                                 brier = "min"))
  expect_true(all(c("calibration_slope", "brier") %in% names(ev)))
  expect_equal(ev$brier[match("GLM", ev$model)],
               brier_score(tab$obs, tab$GLM))
  expect_equal(ev$calibration_slope[match("GLM", ev$model)],
               cc$summary$slope[cc$summary$model == "GLM"])
  # metrics by name, without using them as criteria
  ev2 <- evaluate_models(data = tab, metrics = c("calibration_intercept"))
  expect_true("calibration_intercept" %in% names(ev2))
  expect_false("brier" %in% names(ev2))
  expect_error(evaluate_models(data = tab, metrics = "nope"), "Unknown")
  # "one": the model closest to 1 wins
  df <- data.frame(model = c("a", "b", "c"), s = c(0.5, 1.1, 2))
  expect_equal(pareto_rank(df, maximize = NULL, minimize = NULL,
                           one = "s")$model[1], "b")
  expect_error(calibration_curve(c(0, 1), c(0.5, 2)), "not probabilities")
})

# Smoothed Boyce index ---------------------------------------------------------

# presences drawn with probability f(prediction) from uniform predictions,
# so the true P/E curve is f / mean(f)
draw_presences <- function(n, f) {
  x <- runif(n * 60)
  x[runif(length(x)) < f(x)][seq_len(n)]
}
rising <- function(v) plogis(8 * (v - 0.6))

test_that("smoothed and moving-window Boyce agree on a large sample", {
  set.seed(5)
  bg <- runif(5000)
  pres <- draw_presences(2000, rising)
  b0 <- boyce_index(pres, bg)
  for (m in c("spline", "loess")) {
    expect_equal(boyce_index(pres, bg, smooth = m), b0, tolerance = 0.03)
  }
  # the smoothed P/E curve recovers the true one
  b <- boyce_index(pres, bg, smooth = "spline", details = TRUE)
  truth <- rising(b$curve$prediction) / mean(rising(runif(1e5)))
  mid <- b$curve$prediction > 0.3            # P/E well away from 0
  expect_lt(max(abs(b$curve$pe[mid] / truth[mid] - 1)), 0.2)
  expect_lt(b$flat_p, 1e-6)
})

test_that("the smoothed index varies less with few presences", {
  set.seed(6)
  bg <- runif(3000)
  r <- replicate(40, {
    p <- draw_presences(20, rising)
    c(none = boyce_index(p, bg), spline = boyce_index(p, bg, smooth = "spline"))
  })
  sds <- apply(r, 1, stats::sd, na.rm = TRUE)
  expect_lt(sds[["spline"]], sds[["none"]])
  # and the moving window underestimates it (Liu et al. 2024)
  expect_lt(mean(r["none", ], na.rm = TRUE), mean(r["spline", ]))
})

test_that("Boyce bootstrap interval covers the estimate; curve and plot", {
  set.seed(7)
  bg <- runif(2000)
  pres <- draw_presences(40, rising)
  b <- boyce_index(pres, bg, smooth = "spline", n_boot = 50, seed = 1)
  expect_s3_class(b, "boyce_index")
  expect_true(b$lower <= b$index && b$index <= b$upper)
  expect_length(b$boot, 50)
  expect_true(all(b$curve$pe_lower <= b$curve$pe_upper))
  expect_equal(nrow(b$curve), 501)
  expect_identical(boyce_index(pres, bg, smooth = "spline", n_boot = 50,
                               seed = 1)$boot, b$boot)
  expect_output(print(b), "IN SHORT")
  b0 <- boyce_index(pres, bg, n_boot = 50, seed = 1)
  expect_equal(b0$index, boyce_index(pres, bg))     # default unchanged
  expect_true(b0$lower <= b0$index && b0$index <= b0$upper)
  expect_s3_class(plot_boyce(b), "ggplot")
  expect_s3_class(plot_boyce(b0), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_boyce(b, file = f)
  expect_true(file.exists(f))
  expect_error(plot_boyce(0.5), "details")
})

test_that("Boyce: reversed and flat models are recognised", {
  set.seed(8)
  bg <- runif(3000)
  rev_p <- draw_presences(200, function(v) 1 - rising(v))
  expect_lt(boyce_index(rev_p, bg, smooth = "spline"), -0.9)
  flat <- boyce_index(runif(200), bg, smooth = "spline", details = TRUE)
  expect_gt(flat$flat_p, 0.01)
  expect_output(print(flat), "not distinguishable from flat")
  expect_true(is.na(boyce_index(0.5, runif(10), smooth = "spline")))
})

test_that("the internal spline matches mgcv; Boyce_smooth as a metric", {
  skip_if_not_installed("mgcv")
  set.seed(9)
  bg <- runif(2000)
  pres <- draw_presences(300, rising)
  expect_equal(boyce_index(pres, bg, smooth = "spline"),
               boyce_index(pres, bg, smooth = "gam"), tolerance = 0.02)
  sim <- simulate_enm_data(seed = 1)
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  ev <- evaluate_models(data = tab, criteria = c(AUC = "max",
                                                 Boyce_smooth = "max"))
  expect_equal(ev$Boyce_smooth[ev$model == "GLM"],
               boyce_index(tab$GLM[tab$obs == 1], tab$GLM[tab$obs == 0],
                           smooth = "spline"))
  expect_equal(ev$Boyce, evaluate_models(data = tab)$Boyce[
    match(ev$model, evaluate_models(data = tab)$model)])
})

test_that("thresholds() brings 0-100 predictions to 0-1", {
  set.seed(1)
  pres <- runif(100, 0.4, 1)
  bg <- runif(500)
  t1 <- thresholds(pres, bg, rules = c("p10", "maxSSS"))
  expect_message(t2 <- thresholds(pres * 100, bg * 100,
                                  rules = c("p10", "maxSSS")), "divided by 100")
  expect_equal(t2$threshold, t1$threshold)
  expect_warning(t3 <- thresholds(pres * 100, bg * 100, rules = "p10",
                                  rescale = FALSE), "not on the 0-1 scale")
  expect_equal(t3$threshold, 100 * t1$threshold[1])
})
