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
