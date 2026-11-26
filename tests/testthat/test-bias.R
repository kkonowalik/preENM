bias_fixture <- function() {
  sim <- simulate_enm_data(seed = 1)
  r <- sim$env
  # clustered records: most in one corner, a few spread over M; one per
  # cell (so cell thinning removes nothing and random subsets of them are
  # cell-thinned sets)
  set.seed(4)
  cl <- data.frame(lon = stats::runif(500, 14.2, 16), lat = stats::runif(500,
                                                                         49.2, 50.5))
  sp <- data.frame(lon = stats::runif(40, 14, 24), lat = stats::runif(40, 49, 55))
  pts <- rbind(cl, sp)
  cells <- terra::cellFromXY(r, as.matrix(pts))
  pts <- pts[!duplicated(cells) & !is.na(cells), ]
  list(env = r, pts = pts)
}

test_that("a larger thinning radius keeps fewer records", {
  f <- bias_fixture()
  bs <- bias_sensitivity(f$pts, f$env, thin = c("none", "distance"),
                         thin_km = c(5, 10, 20, 40, 80), bg_n = 500,
                         n_null = 3, total_bins = 400, verbose = FALSE)
  t <- bs$thinning
  d <- t[t$thin_method == "distance", ]
  expect_equal(d$thin_level, paste0(c(5, 10, 20, 40, 80), "km"))
  expect_true(all(diff(d$n) <= 0))
  expect_lt(d$n[5], d$n[1])
  expect_equal(t$n[t$thinning == "none"], nrow(f$pts))
  # no two kept records closer than the radius
  k <- bs$records[["distance_40km"]]
  X <- km_space(k$lon, k$lat)
  expect_gte(min(dbscan::kNN(X, k = 1)$dist), 40)
  # spreading clustered records raises the Clark-Evans ratio
  expect_gt(d$clark_evans[5], t$clark_evans[t$thinning == "none"])
  expect_lt(t$clark_evans[t$thinning == "none"], 1)
})

test_that("environmental thinning keeps more environment than random removal", {
  f <- bias_fixture()
  bs <- bias_sensitivity(f$pts, f$env, thin = c("none", "cell",
                                                "environment"),
                         env_mult = 2, bg_n = 500, n_null = 10,
                         total_bins = 400, verbose = FALSE)
  t <- bs$thinning
  e <- t[t$thin_method == "environment", ]
  expect_equal(nrow(e), 1)
  expect_lt(e$n, nrow(f$pts))
  # a cell-thinned set of the same size (random subsets of the one-per-cell
  # records) covers less of the environment of M
  expect_gt(e$env_coverage, e$null_coverage)
  expect_gt(e$coverage_gain, 0)
  expect_true(e$niche_intact)
  expect_equal(t$n[t$thinning == "cell"], nrow(f$pts))
  expect_output(print(bs), "Removed no record: cell")
})

test_that("backgrounds, combinations, effects, print and plot", {
  f <- bias_fixture()
  ef <- terra::init(f$env[[1]], "x")
  ef <- ef - terra::global(ef, "min")[[1]]
  bs <- bias_sensitivity(f$pts, f$env, thin = c("none", "distance"),
                         thin_km = c(20, 60),
                         background = list(uniform = list(),
                                           effort = list(method = "effort",
                                                         effort = ef,
                                                         effort_floor = 0)),
                         bg_n = c(300, 1000), n_null = 3, total_bins = 400,
                         verbose = FALSE)
  expect_s3_class(bs, "bias_sensitivity")
  expect_equal(nrow(bs$backgrounds), 4)
  expect_equal(nrow(bs$designs), 3 * 4)
  b <- bs$backgrounds
  expect_equal(b$bg_shift[b$bg_method == "uniform"], c(0, 0))
  expect_true(all(b$bg_shift[b$bg_method == "effort"] > 0.1))
  expect_true(all(bs$designs$overlap >= 0 & bs$designs$overlap <= 1))
  expect_true(all(c("thinning method", "thinning radius / width",
                    "background strategy", "background size") %in%
                    bs$effects$choice))
  expect_output(print(bs), "IN SHORT")
  expect_s3_class(plot_bias_sensitivity(bs), "ggplot")
  expect_s3_class(plot_bias_sensitivity(bs, "bg_coverage", "overlap"),
                  "ggplot")
  expect_error(plot_bias_sensitivity(bs, "design"), "numeric")
  fl <- tempfile(fileext = ".png")
  plot_bias_sensitivity(bs, file = fl)
  expect_true(file.exists(fl))
  expect_error(bias_sensitivity(f$pts, f$env, background = list(list())),
               "named list")
})
