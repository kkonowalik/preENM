sim <- simulate_enm_data(seed = 3, n_presence = 120, n_absence = 400)
pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))

test_that("spatial_folds gives k balanced folds and keeps blocks intact", {
  f <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)
  expect_equal(sort(unique(as.integer(f))), 1:4)
  expect_true(all(tabulate(f[pts$occ == 1], 4) > 10))
  blk <- attr(f, "block")
  expect_true(all(tapply(f, blk, function(v) length(unique(v))) == 1))
  expect_equal(as.integer(f), as.integer(spatial_folds(
    pts, k = 4, block_size = 100, occ = pts$occ, seed = 1)))
})

test_that("spatial_folds does not change the user's random numbers", {
  set.seed(99); a <- runif(1)
  set.seed(99); spatial_folds(pts, 4, 100, occ = pts$occ, seed = 1)
  expect_equal(runif(1), a)
})

test_that("spatial_folds: errors and small-data warnings", {
  expect_error(spatial_folds(pts, k = 5), "block_size")
  expect_error(spatial_folds(pts, k = 5, block_size = 5000), "fewer than k")
  small <- rbind(cbind(sim$presence[1:6, ], occ = 1),
                 cbind(sim$absence[1:60, ], occ = 0))
  expect_warning(spatial_folds(small, k = 5, block_size = 60,
                               occ = small$occ, seed = 1), "presences")
})

test_that("spatial_folds accepts projected points", {
  v <- terra::project(terra::vect(pts, geom = c("lon", "lat"),
                                  crs = "EPSG:4326"), "EPSG:3035")
  f <- spatial_folds(v, k = 4, block_size = 100, occ = pts$occ, seed = 1)
  expect_length(f, nrow(pts))
})

# hold-out predictions of two GLMs on spatial folds, two background sets
make_holdout <- function() {
  pts$fold <- spatial_folds(pts, k = 4, block_size = 100, occ = pts$occ,
                            seed = 1)
  pts <- cbind(pts, terra::extract(sim$env, pts[c("lon", "lat")], ID = FALSE))
  pts$set <- ifelse(pts$occ == 1, 0, rep(1:2, length.out = nrow(pts)))
  out <- list()
  for (s in 1:2) for (f in 1:4) {
    d <- pts[pts$set %in% c(0, s), ]
    tr <- d[d$fold != f, ]
    te <- d[d$fold == f, ]
    good <- glm(occ ~ env1 + I(env1^2) + env2 + I(env2^2), binomial, tr)
    poor <- glm(occ ~ env1, binomial, tr)
    out[[length(out) + 1]] <- rbind(
      data.frame(model = "quadratic", set = s, fold = f, obs = te$occ,
                 pred = predict(good, te, type = "response")),
      data.frame(model = "linear", set = s, fold = f, obs = te$occ,
                 pred = predict(poor, te, type = "response")))
  }
  do.call(rbind, out)
}
ho <- suppressWarnings(make_holdout())  # glm: perfect separation is fine here

test_that("evaluate_cv averages per-run metrics correctly", {
  cv <- evaluate_cv(ho, set = "set")
  expect_s3_class(cv, "enm_cv")
  runs <- attr(cv, "runs")
  expect_equal(nrow(runs), 2 * 2 * 4)
  r1 <- ho[ho$model == "quadratic" & ho$set == 1 & ho$fold == 2, ]
  expect_equal(runs$AUC[runs$model == "quadratic" & runs$set == 1 &
                          runs$fold == 2], auc_score(r1$obs, r1$pred))
  q <- cv[cv$model == "quadratic", ]
  expect_equal(q$AUC, mean(runs$AUC[runs$model == "quadratic"]))
  expect_equal(q$AUC_sd, sd(runs$AUC[runs$model == "quadratic"]))
  expect_gt(q$AUC, cv$AUC[cv$model == "linear"])
  expect_true(all(c("AUC_sd_sets", "AUC_sd_folds") %in% names(cv)))
  expect_equal(ncol(attr(cv, "by_fold")), 1 + 4)
  expect_output(print(cv), "background sets")
})

test_that("evaluate_cv: fixed thresholds, factors, empty folds", {
  ho$thr <- 0.5
  a <- evaluate_cv(ho, set = "set")
  b <- evaluate_cv(ho, set = "set", threshold = "thr")
  expect_true(all(b$TSS <= a$TSS + 1e-12))
  ho$model <- factor(ho$model)
  expect_equal(sort(evaluate_cv(ho)$model), c("linear", "quadratic"))
  bad <- ho
  bad$obs[bad$fold == 1] <- 0
  expect_warning(evaluate_cv(bad), "no presences")
  expect_error(evaluate_cv(ho, fold = "nope"), "not found")
})

test_that("full chain: cv -> select -> plot", {
  cv <- evaluate_cv(ho, set = "set")
  s <- select_models(cv, min_metrics = c(AUC = 0.8))
  expect_equal(s$model[s$selected], "quadratic")
  expect_s3_class(plot_selection(s), "ggplot")
})

# Fold audit, environmental folds, spatial sorting bias --------------------------

grid_points <- function() {
  # 400 records on a regular 20 x 20 grid of 0.1 degree, half presences
  xy <- expand.grid(lon = 15 + (0:19) * 0.1 + 0.05,
                    lat = 50 + (0:19) * 0.1 + 0.05)
  xy$occ <- rep(c(1, 0), 200)
  xy
}

test_that("clean block folds pass the audit", {
  pts <- grid_points()
  # four quadrants: blocks of 1 degree
  f <- 1 + (pts$lon > 16) + 2 * (pts$lat > 51)
  r <- terra::rast(xmin = 15, xmax = 17, ymin = 50, ymax = 52,
                   resolution = 0.1, crs = "EPSG:4326")
  terra::values(r) <- seq_len(terra::ncell(r))
  names(r) <- "env1"
  a <- fold_report(pts, f, occ = pts$occ, env = r, block_size = 100)
  expect_false(a$leakage)
  expect_equal(nrow(a$leaks), 0)
  expect_equal(a$folds$n_test, rep(100, 4))
  expect_equal(a$folds$n_train, rep(300, 4))
  expect_equal(a$folds$presences_test, rep(50, 4))
  expect_true(all(a$folds$duplicates == 0 & a$folds$shared_cells == 0))
  # neighbouring grid records across the quadrant border are 0.1 degree
  # (about 7 to 11 km) apart
  expect_true(all(a$folds$min_dist_km > 6 & a$folds$min_dist_km < 12))
  expect_output(print(a), "No leakage found")
  expect_s3_class(plot_folds(a), "ggplot")
})

test_that("duplicate points and cells split across folds are leakage", {
  pts <- grid_points()
  f <- 1 + (pts$lon > 16) + 2 * (pts$lat > 51)
  # the same record twice, in folds 1 and 4
  dup <- rbind(pts, pts[1, ])
  fd <- c(f, 4)
  a <- fold_report(dup, fd, occ = dup$occ)
  expect_true(a$leakage)
  expect_equal(sort(a$leaks$index), c(1, 401))
  expect_equal(sort(a$leaks$fold), c(1, 4))
  expect_equal(a$folds$duplicates, c(1, 0, 0, 1))
  expect_output(print(a), "LEAKAGE FOUND in folds 1, 4")
  # a record moved inside the cell of another record, in another fold
  r <- terra::rast(xmin = 15, xmax = 17, ymin = 50, ymax = 52,
                   resolution = 0.1, crs = "EPSG:4326", vals = 1)
  same_cell <- rbind(pts, data.frame(lon = 15.07, lat = 50.07, occ = 0))
  a2 <- fold_report(same_cell, c(f, 2), occ = same_cell$occ, env = r)
  expect_true(a2$leakage)
  expect_equal(unique(a2$leaks$type), "same raster cell")
  expect_equal(sort(a2$leaks$index), c(1, 401))
  # within the same fold neither is leakage
  expect_false(fold_report(dup, c(f, 1))$leakage)
  expect_error(fold_report(pts, f[-1]), "one value per point")
})

test_that("random folds are reported as not separated", {
  set.seed(3)
  pts <- grid_points()
  a <- fold_report(pts, sample(rep(1:4, 100)))
  expect_lt(stats::median(a$folds$separation), 1.5)
  expect_output(print(a), "not separated")
  f <- 1 + (pts$lon > 16) + 2 * (pts$lat > 51)
  expect_gt(stats::median(fold_report(pts, f)$folds$separation),
            stats::median(a$folds$separation))
})

test_that("evaluate_cv audits the folds and warns on leakage", {
  pts <- grid_points()
  f <- 1 + (pts$lon > 16) + 2 * (pts$lat > 51)
  set.seed(4)
  d <- rbind(data.frame(model = "a", fold = f, obs = pts$occ,
                        pred = plogis(2 * pts$occ - 1 + rnorm(400)),
                        lon = pts$lon, lat = pts$lat),
             data.frame(model = "b", fold = f, obs = pts$occ,
                        pred = runif(400), lon = pts$lon, lat = pts$lat))
  cv <- evaluate_cv(d, coords = c("lon", "lat"))
  expect_s3_class(attr(cv, "audit"), "fold_report")
  expect_false(attr(cv, "audit")$leakage)
  expect_output(print(cv), "No leakage found")
  leaky <- rbind(d, transform(d[d$lon < 15.1 & d$lat < 50.1, ], fold = 4))
  expect_warning(cv2 <- evaluate_cv(leaky, coords = c("lon", "lat")),
                 "Leakage between cross-validation folds")
  expect_true(attr(cv2, "audit")$leakage)
  # metrics unchanged by the audit
  expect_equal(cv$AUC, evaluate_cv(d)$AUC)
})

test_that("environmental_folds separate distinct environments", {
  # two environments: west cold, east warm
  r <- terra::rast(xmin = 0, xmax = 4, ymin = 0, ymax = 2, resolution = 0.1,
                   crs = "EPSG:4326")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  t1 <- terra::rast(r, vals = ifelse(xy[, 1] < 2, 0, 10) +
                      (xy[, 2] - 1) * 0.1)
  t2 <- terra::rast(r, vals = ifelse(xy[, 1] < 2, 5, 20) + xy[, 2] * 0.05)
  env <- c(t1, t2)
  names(env) <- c("temp", "moist")
  set.seed(5)
  pts <- data.frame(lon = runif(200, 0.05, 3.95), lat = runif(200, 0.05, 1.95))
  f <- environmental_folds(pts, env, k = 2, seed = 1, pca_n = 500)
  # each fold holds exactly one side
  tab <- table(f, west = pts$lon < 2)
  expect_true(all(rowSums(tab > 0) == 1))
  expect_equal(nrow(attr(f, "centers")), 2)
  a <- fold_report(pts, f, env = env)
  expect_true(all(a$folds$env_overlap < 0.2))
  expect_warning(environmental_folds(pts, env, k = 2, seed = 1, pca_n = 500,
                                     occ = as.numeric(pts$lon < 2)),
                 "no presences")
})

test_that("spatial sorting bias: ~1 for a random background, < 1 far away", {
  set.seed(6)
  box <- function(n, x0, x1, y0, y1) {
    data.frame(lon = runif(n, x0, x1), lat = runif(n, y0, y1))
  }
  train <- box(300, 10, 14, 50, 53)
  test <- box(150, 10, 14, 50, 53)
  s <- spatial_sorting_bias(test, box(2000, 10, 14, 50, 53), train, seed = 1)
  expect_equal(s$ssb, 1, tolerance = 0.1)
  # background drawn far from the presences: Da >> Dp, SSB near 0
  far <- spatial_sorting_bias(test, box(2000, 0, 30, 40, 65), train, seed = 1)
  expect_lt(far$ssb, 0.2)
  expect_output(print(far), "strong")
  # pairwise distance sampling removes the bias
  expect_equal(far$ssb_matched, 1, tolerance = 0.15)
  expect_equal(length(far$presence_index), length(far$background_index))
  expect_false(anyDuplicated(far$background_index) > 0)
  # Dp by hand for one point: a test presence on top of a training record
  s0 <- spatial_sorting_bias(train[1:2, ], train[3:4, ], train[1:2, ])
  expect_equal(s0$Dp, 0)
  expect_equal(s0$ssb, 0)
})
