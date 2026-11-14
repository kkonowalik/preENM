set.seed(5)
bg <- expand.grid(cell = 1:400, slice = 1:4,
                  model = c("good", "random", "reverse"),
                  stringsAsFactors = FALSE)
bg$pred <- runif(nrow(bg))
st <- expand.grid(record = 1:60, model = c("good", "random", "reverse"),
                  stringsAsFactors = FALSE)
st$slice <- rep(sample(1:4, 60, replace = TRUE), 3)
st$pred <- ifelse(st$model == "good", 0.999,
                  ifelse(st$model == "reverse", 0.0001, runif(nrow(st))))
st$period <- rep(ifelse(1:60 <= 40, "Holocene",
                        ifelse(1:60 <= 50, "Lateglacial", "Glacial")), 3)

test_that("temporal AUC: perfect = 1, random ~ 0.5, reversed = 0", {
  tv <- temporal_validation(st, bg)
  v <- setNames(tv$temporal_AUC, tv$model)
  expect_gt(v[["good"]], 0.99)
  expect_lt(v[["reverse"]], 0.01)
  expect_equal(v[["random"]], 0.5, tolerance = 0.1)
})

test_that("with one slice it equals the classical AUC", {
  b1 <- bg[bg$slice == 1 & bg$model == "random", ]
  s1 <- data.frame(record = 1:30, slice = 1, model = "random",
                   pred = round(runif(30), 1))
  b1$pred <- round(b1$pred, 1)                      # ties count half
  tv <- temporal_validation(s1, b1)
  expect_equal(tv$temporal_AUC,
               auc_score(c(rep(1, 30), rep(0, nrow(b1))),
                         c(s1$pred, b1$pred)))
})

test_that("records spanning several slices are averaged first", {
  b <- data.frame(model = "m", slice = rep(1:2, each = 4),
                  pred = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))
  s <- data.frame(model = "m", record = c("r1", "r1", "r2"),
                  slice = c(1, 2, 2), pred = c(0.35, 0.65, 0.9))
  # r1: slice 1 -> 3/4, slice 2 -> 2/4 => 0.625 ; r2: 1 ; mean = 0.8125
  expect_equal(temporal_validation(s, b)$temporal_AUC, 0.8125)
})

test_that("periods, combined periods, analog and sensitivity", {
  st$novel <- st$record %% 10 == 0
  tv <- temporal_validation(
    st, bg, period = "period", novel = "novel",
    combine = list(deep = c("Lateglacial", "Glacial")),
    thresholds = c(good = 0.5, random = 0.5, reverse = 0.5),
    prefix = "fossil")
  expect_true(all(c("fossil_AUC_Holocene", "fossil_AUC_Lateglacial",
                    "fossil_AUC_Glacial", "fossil_AUC_deep",
                    "fossil_AUC_analog", "fossil_sens",
                    "fossil_n_records") %in% names(tv)))
  r <- tv[tv$model == "random", ]
  runs <- attr(tv, "runs")
  # pooled: weighted by the number of records (10 Lateglacial, 10 Glacial)
  expect_equal(r$fossil_AUC_deep,
               mean(c(r$fossil_AUC_Lateglacial, r$fossil_AUC_Glacial)))
  expect_equal(tv$fossil_sens[tv$model == "good"], 1)
  expect_equal(tv$fossil_sens[tv$model == "reverse"], 0)
  expect_equal(r$fossil_n_records, 60)
})

test_that("runs are averaged with SD; thresholds as a data.frame", {
  st2 <- rbind(cbind(st, run = "a"), cbind(st, run = "b"))
  bg2 <- rbind(cbind(bg, run = "a"), cbind(bg, run = "b"))
  bg2$pred[bg2$run == "b"] <- runif(sum(bg2$run == "b"))
  thr <- data.frame(model = rep(c("good", "random", "reverse"), 2),
                    run = rep(c("a", "b"), each = 3), threshold = 0.5)
  tv <- temporal_validation(st2, bg2, run = "run", thresholds = thr)
  expect_true(all(c("temporal_AUC_sd", "temporal_n_runs") %in% names(tv)))
  expect_equal(nrow(attr(tv, "runs")), 6)
  expect_equal(unique(tv$temporal_n_runs), 2)
})

test_that("informative errors and warnings", {
  expect_error(temporal_validation(st, bg[bg$model != "good", ]),
               "No background")
  expect_warning(temporal_validation(st, bg[bg$slice != 4, ]),
                 "No background cells for slice")
  bad <- st
  bad$period[1] <- "Glacial"
  bad$period[bad$record == 1] <- c("Glacial", "Holocene", "Holocene")
  expect_error(temporal_validation(bad, bg, period = "period"), "constant")
  expect_error(temporal_validation(st, bg, record = "x"), "not found")
})

test_that("age_to_period uses inclusive upper limits", {
  expect_equal(age_to_period(c(0, 11700, 11701, 14700, 14701, 30000)),
               c("Holocene", "Holocene", "Lateglacial", "Lateglacial",
                 "Glacial", "Glacial"))
  expect_error(age_to_period(1, c(a = 10, b = 5)), "increasing")
})

test_that("anomaly_climate: data.frames and rasters agree", {
  cur <- data.frame(bio1 = c(5, 10), bio12 = c(600, 900))
  past <- data.frame(bio1 = c(-3, 2), bio12 = c(300, 5))
  pres <- data.frame(bio1 = c(4, 9), bio12 = c(650, 850))
  out <- anomaly_climate(cur, past, pres, ratio_vars = "bio12")
  expect_equal(out$bio1, c(5 - 7, 10 - 7))
  expect_equal(out$bio12[1], 600 * 301 / 651)
  expect_equal(out$bio12[2], 900 * 0.1)                # ratio clamped
  mk <- function(d) {
    r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 2)
    terra::values(r) <- as.matrix(d)
    names(r) <- names(d)
    r
  }
  ro <- anomaly_climate(mk(cur), mk(past), mk(pres), ratio_vars = "bio12")
  expect_equal(unname(terra::values(ro)[, "bio12"]), out$bio12)
  expect_equal(unname(terra::values(ro)[, "bio1"]), out$bio1)
  expect_error(anomaly_climate(cur, past["bio1"], pres), "same variable")
})

test_that("outside_range flags non-analog rows", {
  ref <- data.frame(t = c(0, 10), p = c(100, 900))
  x <- data.frame(p = c(500, 500, 950), t = c(5, -2, 5))
  expect_equal(outside_range(x, ref), c(FALSE, TRUE, TRUE))
})

test_that("temporal_validation gives predicted area and lift", {
  set.seed(2)
  bg <- expand.grid(cell = 1:400, slice = 1:2, model = c("good", "all"))
  bg$pred <- ifelse(bg$model == "all", 0.9, runif(nrow(bg)))
  st <- expand.grid(record = 1:30, model = c("good", "all"))
  st$slice <- rep(1:2, length.out = nrow(st))
  st$pred <- ifelse(st$model == "all", 0.9, runif(nrow(st), 0.6, 1))
  tv <- temporal_validation(st, bg, thresholds = c(good = 0.5, all = 0.5),
                            prefix = "f")
  good <- tv[tv$model == "good", ]
  all <- tv[tv$model == "all", ]
  expect_equal(all$f_sens, 1)
  expect_equal(all$f_area, 1)
  expect_equal(all$f_lift, 1)                  # predicts everything
  expect_equal(good$f_sens, 1)
  expect_equal(good$f_area, 0.5, tolerance = 0.06)
  expect_gt(good$f_lift, 1.8)
})

test_that("out-of-fold records are pooled across runs", {
  set.seed(3)
  bg <- expand.grid(cell = 1:300, slice = 1, run = 1:2, model = "m")
  bg$pred <- runif(nrow(bg))
  # run 1 scores 2 records (both poor), run 2 scores 18 records (good)
  st <- data.frame(model = "m", run = rep(1:2, c(2, 18)), record = 1:20,
                   slice = 1, pred = c(0.1, 0.1, rep(0.9, 18)))
  pooled <- temporal_validation(st, bg, run = "run")
  bg1 <- bg[bg$run == 1, ]
  manual <- function(x, v) (findInterval(x, sort(v), left.open = TRUE) +
                              findInterval(x, sort(v))) / (2 * length(v))
  sc <- c(manual(0.1, bg1$pred), manual(0.1, bg1$pred),
          rep(manual(0.9, bg$pred[bg$run == 2]), 18))
  expect_true(attr(pooled, "pooled"))
  expect_equal(pooled$temporal_AUC, mean(sc))
  expect_equal(pooled$temporal_n_records, 20)
  means <- temporal_validation(st, bg, run = "run", pool = FALSE)
  expect_false(attr(means, "pooled"))
  expect_equal(means$temporal_AUC, mean(c(mean(sc[1:2]), mean(sc[3:20]))))
  expect_lt(means$temporal_AUC, pooled$temporal_AUC)
  expect_equal(nrow(attr(pooled, "records")), 20)
  # records scored by every run: run values are averaged (auto)
  st2 <- rbind(transform(st, run = 1), transform(st, run = 2))
  both <- temporal_validation(st2, bg, run = "run")
  expect_false(attr(both, "pooled"))
})

test_that("record bootstrap is paired across models", {
  set.seed(5)
  n_rec <- 150
  bg <- expand.grid(cell = 1:500, slice = 1:3, model = c("A", "B", "C"))
  bg$pred <- runif(nrow(bg))
  # records differ a lot in difficulty; B is A minus a small constant
  base <- runif(n_rec, 0.5, 0.98)
  st <- data.frame(record = rep(1:n_rec, 3), slice = rep(sample(1:3, n_rec,
                   replace = TRUE), 3), model = rep(c("A", "B", "C"),
                   each = n_rec))
  st$pred <- c(base, base - 0.02, runif(n_rec))
  st$period <- rep(ifelse(1:n_rec <= 100, "Holocene", "Glacial"), 3)
  tv <- temporal_validation(st, bg, period = "period", n_boot = 300,
                            seed = 1, combine = list(deep = "Glacial"))
  expect_equal(attr(tv, "compare")[["metric"]], "temporal_AUC_deep")
  a <- tv[tv$model == "A", ]
  b <- tv[tv$model == "B", ]
  # separate intervals overlap, but the paired gap is clearly above 0
  expect_gt(a$temporal_AUC_deep_upper, b$temporal_AUC_deep_lower)
  expect_gt(b$temporal_AUC_deep_diff_lower, 0)
  expect_true(b$temporal_AUC_deep_worse)
  expect_gt(a$temporal_AUC_deep_p_best, 0.95)
  expect_equal(a$temporal_AUC_deep_diff, 0)
  expect_false(a$temporal_AUC_deep_worse)
  # a margin larger than the gap: not flagged
  tvm <- temporal_validation(st, bg, period = "period", n_boot = 300,
                             seed = 1, combine = list(deep = "Glacial"),
                             margin = 0.1)
  expect_false(tvm$temporal_AUC_deep_worse[tvm$model == "B"])
  expect_true(tvm$temporal_AUC_deep_worse[tvm$model == "C"])
  # reproducible; intervals around the estimate
  tv2 <- temporal_validation(st, bg, period = "period", n_boot = 300,
                             seed = 1, combine = list(deep = "Glacial"))
  expect_equal(tv$temporal_AUC_lower, tv2$temporal_AUC_lower)
  expect_true(all(tv$temporal_AUC_lower <= tv$temporal_AUC &
                    tv$temporal_AUC <= tv$temporal_AUC_upper))
  expect_equal(nrow(attr(tv, "bootstrap")), 300 * 3)
  expect_error(temporal_validation(st, bg, n_boot = 10, compare = "nope"),
               "must be one of")
})
