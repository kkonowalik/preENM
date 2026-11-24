# Niche overlap

overlap_env <- function() {
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 6, ymin = 45,
                   ymax = 51, crs = "EPSG:4326")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  set.seed(1)
  env <- c(terra::rast(r, vals = xy[, 1] * 3),
           terra::rast(r, vals = xy[, 2] - 45),
           terra::rast(r, vals = xy[, 1] + xy[, 2] + rnorm(nrow(xy), 0, 0.2)))
  names(env) <- c("temp", "rain", "mix")
  env
}
pts <- function(n, x, y, group, seed) {
  set.seed(seed)
  data.frame(lon = runif(n, x[1], x[2]), lat = runif(n, y[1], y[2]),
             group = group)
}

test_that("identical record sets give D = 1 and no dynamics", {
  env <- overlap_env()
  a <- pts(80, c(0.5, 3), c(46, 49), "a", 1)
  rec <- rbind(a, transform(a, group = "b"))
  no <- niche_overlap(rec, env, n_rep = 0)
  expect_equal(no$pairs$D, 1)
  expect_equal(no$pairs$I, 1)
  expect_equal(no$pairs$expansion, 0)
  expect_equal(no$pairs$unfilling, 0)
  # a named list gives the same as a grouping column
  nl <- niche_overlap(list(a = a, b = a), env, n_rep = 0)
  expect_equal(nl$pairs$D, 1)
})

test_that("disjoint environments give D near 0; tests reject", {
  env <- overlap_env()
  rec <- rbind(pts(60, c(0.2, 1.2), c(45.2, 46.5), "cold_dry", 2),
               pts(60, c(4.8, 5.8), c(49.5, 50.8), "warm_wet", 3))
  no <- niche_overlap(rec, env, n_rep = 49, seed = 1)
  expect_lt(no$pairs$D, 0.05)
  expect_lt(no$pairs$eq_p_lower, 0.05)
  expect_lt(no$pairs$eq_ses, -2)
  # both are less similar than random points of the other's range
  expect_lt(no$pairs$sim12_p_lower, 0.05)
  expect_output(print(no), "not equivalent")
  expect_output(print(no), "no or very limited overlap")
})

test_that("a random split of one set is not rejected", {
  env <- overlap_env()
  rec <- pts(160, c(0.5, 4), c(46, 50), "x", 4)
  set.seed(5)
  rec$group <- sample(c("x", "y"), nrow(rec), replace = TRUE)
  no <- niche_overlap(rec, env, n_rep = 99, seed = 2)
  expect_gt(no$pairs$D, 0.6)
  # the test is calibrated: random splits are rejected about 5% of the time
  p <- vapply(1:20, function(k) {
    set.seed(100 + k)
    rec$group <- sample(c("x", "y"), nrow(rec), replace = TRUE)
    niche_overlap(rec, env, n_rep = 39, seed = k,
                  tests = "equivalency")$pairs$eq_p_lower
  }, numeric(1))
  expect_lte(mean(p < 0.05), 0.15)
})

test_that("breadth, expansion and unfilling move as expected", {
  env <- overlap_env()
  narrow <- pts(80, c(2, 3), c(47, 48), "narrow", 6)
  wide <- pts(80, c(0.5, 5.5), c(45.5, 50.5), "wide", 7)
  no <- niche_overlap(rbind(narrow, wide), env, n_rep = 0)
  g <- no$groups
  expect_gt(g$breadth_B2[g$group == "wide"],
            g$breadth_B2[g$group == "narrow"])
  # relative to the narrow group, the wide one expands; nothing is unfilled
  p <- no$pairs
  expect_equal(c(p$group1, p$group2), c("narrow", "wide"))
  expect_gt(p$expansion, 0.3)
  expect_lt(p$unfilling, 0.05)
  expect_equal(p$expansion + p$stability, 1)
  # Levins B2 by hand: uniform occupancy over n cells = 1, one cell = 0
  expect_equal(levins_b2(matrix(1, 3, 3), matrix(TRUE, 3, 3)), 1)
  expect_equal(levins_b2(diag(c(1, 0, 0)), matrix(TRUE, 3, 3)), 0)
})

test_that("raw space, backgrounds per group and prediction overlap", {
  env <- overlap_env()
  rec <- rbind(pts(50, c(0.5, 2.5), c(46, 48), "a", 8),
               pts(50, c(2, 4), c(47, 49), "b", 9))
  raw <- niche_overlap(rec, env, space = "raw", vars = c("temp", "rain"),
                       n_rep = 0)
  expect_equal(raw$space$axes, c("temp", "rain"))
  areas <- list(a = terra::vect("POLYGON ((0 45, 3 45, 3 49, 0 49, 0 45))",
                                crs = "EPSG:4326"),
                b = terra::vect("POLYGON ((1.5 46, 4.5 46, 4.5 50, 1.5 50, 1.5 46))",
                                crs = "EPSG:4326"))
  bg <- niche_overlap(rec, env, background = areas, n_rep = 0)
  expect_equal(bg$settings$background, "given")
  expect_true(all(bg$groups$n_background > 100))
  # geographic overlap of predictions: identical layers -> 1, 1, 1
  s <- c(env[["temp"]], env[["temp"]])
  names(s) <- c("a", "b")
  po <- prediction_overlap(s / 20)
  expect_equal(c(po$D, po$I, po$rho), c(1, 1, 1))
  # two cells by hand: p1 = (0.25, 0.75), p2 = (0.75, 0.25): D = 0.5
  r2 <- terra::rast(nrows = 1, ncols = 2, vals = c(1, 3))
  r3 <- terra::rast(nrows = 1, ncols = 2, vals = c(3, 1))
  po2 <- prediction_overlap(c(r2, r3))
  expect_equal(po2$D, 0.5)
  expect_equal(po2$I, 1 - 0.5 * 2 * (sqrt(0.75) - sqrt(0.25))^2)
  no <- niche_overlap(rec, env, n_rep = 0, predictions = s / 20)
  expect_equal(no$pairs$geo_D, 1)
  expect_s3_class(plot_niche_overlap(no), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_niche_overlap(no, file = f)
  expect_true(file.exists(f))
  expect_error(niche_overlap(rec[rec$group == "a", ], env), "two groups")
})

test_that("overlap_scenarios reports robust and fragile verdicts", {
  env <- overlap_env()
  rec <- rbind(pts(60, c(0.2, 1.2), c(45.2, 46.5), "a", 10),
               pts(60, c(4.8, 5.8), c(49.5, 50.8), "b", 11))
  os <- overlap_scenarios(rec, env, thin_km = c(0, 20), grids = c(40, 80),
                          spaces = "pca", n_rep = 19, seed = 1)
  expect_equal(nrow(os$results), 4)
  expect_true(all(os$robustness$robust))
  expect_output(print(os), "robust to these choices")
  # a choice that decides the verdict is named
  near <- rbind(pts(60, c(1, 3), c(46, 48), "a", 12),
                pts(60, c(1.3, 3.3), c(46.3, 48.3), "b", 13))
  bgs <- list(whole = NULL,
              tiny = terra::vect("POLYGON ((1 46, 3.3 46, 3.3 48.3, 1 48.3, 1 46))",
                                 crs = "EPSG:4326"))
  os2 <- overlap_scenarios(near, env, backgrounds = bgs, spaces = "pca",
                           thin_km = 0, grids = 60, n_rep = 19, seed = 1)
  expect_equal(os2$settings$varying, "background")
  expect_error(overlap_scenarios(rec, env, backgrounds = list(NULL)),
               "named list")
})
