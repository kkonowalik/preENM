test_that("check_sample_size applies the rules of thumb", {
  x <- check_sample_size(18, n_predictors = 6, k = 5)
  expect_equal(x$value, c(18, 3, 3.6))
  expect_equal(x$ok, c(TRUE, FALSE, FALSE))
  expect_match(x$note[1], "narrow-ranged")
  expect_match(x$note[2], "at most 1 predictors")
  expect_match(x$note[3], "at most 3 folds")
  expect_false(check_sample_size(10)$ok)
  expect_match(check_sample_size(40)$note, "widespread")
  # with rasters: occupied cells, not records
  sim <- simulate_enm_data(seed = 1, n_presence = 60)
  occ <- rbind(sim$presence, sim$presence)            # every record twice
  y <- check_sample_size(occ, sim$env)
  cells <- unique(terra::cellFromXY(sim$env, as.matrix(sim$presence)))
  expect_equal(y$value[1], length(cells))
  expect_equal(y$value[2], signif(length(cells) / terra::ncell(sim$env), 3))
  expect_equal(y$value[3], length(cells) / 2)
})

test_that("match_climate_period labels the records", {
  occ <- data.frame(lon = 1:6, lat = 50,
                    year = c("1950", "1981", "2010", "2012", "2021", NA))
  expect_message(m <- match_climate_period(occ, c(1981, 2010)),
                 "2 within, 1 before, 2 after, 1 unknown")
  expect_equal(m$period_match, c("before", "within", "within", "after",
                                 "after", "unknown"))
  expect_equal(m$years_outside, c(-31, 0, 0, 2, 11, NA))
  s <- attr(m, "summary")
  expect_equal(s$n, c(2L, 1L, 2L, 1L))
  expect_equal(s$median_year, c(1995.5, 1950, 2016.5, NA))
  m2 <- match_climate_period(occ, c(1981, 2010), tolerance = 5,
                             verbose = FALSE)
  expect_equal(m2$period_match[4], "within")
  expect_equal(m2$years_outside[5], 6)
  expect_error(match_climate_period(occ[1:2], c(1981, 2010)), "No year")
  expect_error(match_climate_period(occ, c(2010, 1981)), "first and the last")
})

test_that("suggest_block_size finds the range of a patchy layer", {
  # 200 x 200 km at 2 km; patches of 20 km and of 60 km, and a trend
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 2e5, ymin = 0,
                   ymax = 2e5, crs = "EPSG:3035")
  set.seed(1)
  patches <- function(km) {
    coarse <- terra::rast(terra::ext(r), resolution = km * 1000,
                          crs = terra::crs(r))
    terra::values(coarse) <- rnorm(terra::ncell(coarse))
    terra::resample(coarse, r, method = "near")
  }
  env <- c(patches(20), patches(50), terra::init(r, "x"))
  names(env) <- c("small", "large", "trend")
  b <- suggest_block_size(env, n = 1500)
  expect_s3_class(b, "block_size")
  rg <- stats::setNames(b$ranges$range_km, b$ranges$layer)
  # values are independent beyond one patch: between the side and the
  # diagonal of a patch
  expect_gt(rg[["small"]], 10)
  expect_lt(rg[["small"]], 20 * sqrt(2) + 5)
  expect_gt(rg[["large"]], 30)
  expect_lt(rg[["large"]], 50 * sqrt(2) + 5)
  expect_gt(rg[["large"]], rg[["small"]])
  expect_equal(b$block_size_km, signif(stats::median(rg), 3))
  expect_equal(unique(b$variogram$layer), names(env))
  expect_output(print(b), "Suggested block size")
  # usable directly in spatial_folds()
  pts <- as.data.frame(terra::spatSample(env, 200, xy = TRUE)[, c("x", "y")])
  f <- spatial_folds(pts, k = 3, block_size = b$block_size_km,
                     crs = "EPSG:3035", seed = 1)
  expect_equal(sort(unique(as.vector(f))), 1:3)
})

test_that("suggest_block_size works in longitude/latitude", {
  sim <- simulate_enm_data(seed = 1)
  b <- suggest_block_size(sim$env, n = 500)
  expect_true(is.finite(b$block_size_km))
  expect_true(all(b$variogram$distance_km <= b$settings$max_km))
})
