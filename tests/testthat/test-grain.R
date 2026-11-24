grain_env <- function() {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 20, ymin = 40,
                   ymax = 60, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  set.seed(1)
  a <- terra::rast(r, vals = xy[, 1] / 20 + stats::rnorm(nrow(xy), 0, 0.05))
  b <- terra::rast(r, vals = sin(xy[, 2] / 3) + stats::rnorm(nrow(xy), 0,
                                                             0.05))
  out <- c(a, b)
  names(out) <- c("a", "b")
  out
}
grain_records <- function(n = 300) {
  set.seed(2)
  data.frame(species = "sp", lon = stats::runif(n, 0.1, 19.9),
             lat = stats::runif(n, 40.1, 59.9))
}

test_that("aggregation on disk matches terra::aggregate()", {
  r <- grain_env()
  r[[1]][1:3] <- NA
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f)
  mixed <- c(terra::rast(f)[[2]], round(r[[1]] * 100))  # file + memory
  names(mixed) <- c("b", "a100")
  d <- tempfile()
  a <- aggregate_on_disk(mixed, 4, "average", d)
  ref <- terra::aggregate(mixed, 4, fun = "mean", na.rm = TRUE)
  expect_equal(names(a), c("b", "a100"))
  expect_equal(dim(a), dim(ref))
  expect_equal(as.vector(terra::ext(a)), as.vector(terra::ext(ref)))
  expect_equal(terra::values(a), terra::values(ref), tolerance = 1e-5,
               ignore_attr = TRUE)
  # integer averages are not rounded
  expect_true(any(terra::values(a[[2]]) %% 1 != 0, na.rm = TRUE))
})

test_that("one row per grain; coarser grains keep fewer records", {
  g <- grain_sensitivity(grain_records(), grain_env(), factors = c(1, 2, 4),
                         sizes = 2, n_sample = 1000, env_thin = FALSE,
                         uncertainty = NA, seed = 1, verbose = FALSE)
  expect_s3_class(g, "grain_sensitivity")
  expect_equal(nrow(g$table), 3)
  expect_equal(g$table$factor, c(1, 2, 4))
  expect_true(all(diff(g$table$records_cell) <= 0))
  expect_lt(g$table$records_cell[3], g$table$records_cell[1])
  # 300 records in 1600 cells keep most; in 100 cells at most 100 remain
  expect_lte(g$table$records_cell[3], 100)
  expect_true(all(is.na(g$table$unc_lost_share)))
  expect_output(print(g), "IN SHORT")
})

test_that("agreement is 1 when the same set is selected at every grain", {
  # two independent layers: the only passing set of two is "a b"
  g <- grain_sensitivity(grain_records(), grain_env(), factors = c(1, 2, 4),
                         sizes = 2, n_sample = 1000, env_thin = FALSE,
                         uncertainty = NA, seed = 1, verbose = FALSE)
  expect_equal(g$table$selected, rep("a b", 3))
  expect_equal(g$table$jaccard, c(1, 1, 1))
  expect_equal(g$set_agreement, 1)
  expect_equal(g$stable_factor, 4)
  expect_output(print(g), "same at every grain tested")
})

test_that("a selection that changes is detected", {
  expect_equal(jaccard_sets("a b c", "a b d"), 0.5)
  expect_equal(jaccard_sets("a b", "b a"), 1)
  expect_true(is.na(jaccard_sets(NA, "a")))
  # a third layer collinear with `a` only at coarse grains: at the finest
  # grain the set "a c" passes, at coarse grains it does not
  r <- grain_env()
  set.seed(5)
  noise <- terra::rast(r[[1]], vals = stats::rnorm(terra::ncell(r), 0, 1))
  c3 <- r[["a"]] + noise          # noise averages out when aggregated
  env <- c(r, c3)
  names(env) <- c("a", "b", "c")
  g <- grain_sensitivity(grain_records(), env, factors = c(1, 8),
                         sizes = 3, size = 3, thr_r = 0.7, n_sample = 1500,
                         env_thin = FALSE, uncertainty = NA, seed = 1,
                         verbose = FALSE)
  expect_equal(g$table$selected[1], "a b c")
  expect_true(is.na(g$table$selected[2]))
  expect_equal(g$stable_factor, 1)
  expect_gt(g$table$collinear_pairs[2], g$table$collinear_pairs[1])
})

test_that("records with no uncertainty are never lost; thinning runs", {
  rec <- grain_records(150)
  rec$coordinateUncertaintyInMeters <- 0
  g <- grain_sensitivity(rec, grain_env(), factors = c(1, 2), sizes = 2,
                         n_sample = 1000, seed = 1, verbose = FALSE)
  expect_equal(g$table$unc_lost, c(0L, 0L))
  expect_equal(g$table$unc_examined, c(150L, 150L))
  expect_true(all(g$table$records_env <= g$table$records_cell))
  expect_true(all(g$table$records_env > 0))
})

test_that("plot and write work", {
  g <- grain_sensitivity(grain_records(), grain_env(), factors = c(1, 2),
                         sizes = 2, n_sample = 1000, env_thin = FALSE,
                         uncertainty = NA, seed = 1, verbose = FALSE)
  expect_s3_class(plot_grain_sensitivity(g), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_grain_sensitivity(g, file = f)
  expect_true(file.exists(f))
  d <- tempfile()
  out <- write_grain_sensitivity(g, d)
  expect_true(all(file.exists(out)))
  expect_equal(nrow(utils::read.csv(out[1])), 2)
})
