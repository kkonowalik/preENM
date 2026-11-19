cov_fixture <- function() {
  set.seed(1)
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 20, ymin = 40,
                   ymax = 60, crs = "EPSG:4326")
  temp <- 20 - (terra::init(r, "y") - 40)
  wet <- terra::init(r, "x") / 2 +
    terra::setValues(r, stats::rnorm(terra::ncell(r), 0, 0.3))
  env <- c(temp, wet)
  names(env) <- c("temp", "wet")
  cells <- sample(terra::ncell(r), 500)
  pts <- as.data.frame(terra::xyFromCell(r, cells))
  names(pts) <- c("lon", "lat")
  box <- function(y1, y2, x1 = 0, x2 = 20) {
    terra::as.polygons(terra::ext(x1, x2, y1, y2), crs = "EPSG:4326")
  }
  list(env = env, pts = pts, box = box)
}

test_that("a region with part of the range misses part of the niche", {
  f <- cov_fixture()
  small <- niche_coverage(f$pts, f$env, f$box(56, 60, 0, 6))
  half <- niche_coverage(f$pts, f$env, f$box(45, 55))
  all <- niche_coverage(f$pts, f$env, f$box(40.2, 59.8))
  expect_s3_class(half, "niche_coverage")
  expect_gt(small$summary$unfilling, 0.8)
  expect_equal(small$summary$advice, "nested model strongly advised")
  expect_gt(half$summary$unfilling, 0.2)
  expect_lt(half$summary$unfilling, 0.6)
  expect_lt(all$summary$unfilling, 0.05)
  expect_equal(all$summary$advice, "regional model sufficient")
  expect_lte(half$summary$not_available, half$summary$unfilling + 1e-9)
  v <- half$variables
  # the middle half covers about half of the temperature range, all of wet
  expect_equal(v$coverage[v$variable == "temp"], 0.5, tolerance = 0.1)
  expect_gt(v$coverage[v$variable == "wet"], 0.95)
  expect_equal(v$outside_records[v$variable == "temp"],
               mean(f$pts$lat < 45 | f$pts$lat > 55), tolerance = 0.05)
  expect_equal(half$summary$records_region,
               sum(f$pts$lat > 45 & f$pts$lat < 55))
  expect_output(print(half), "nested model worth considering")
  expect_s3_class(plot_niche_coverage(half), "ggplot")
  # stricter limits change the advice, not the numbers
  h2 <- niche_coverage(f$pts, f$env, f$box(45, 55), limits = c(0.6, 0.9))
  expect_equal(h2$summary$unfilling, half$summary$unfilling)
  expect_equal(h2$summary$advice, "regional model sufficient")
  expect_error(niche_coverage(f$pts, f$env, f$box(45, 55),
                              limits = c(0.5, 0.1)), "limits")
  expect_error(suppressWarnings(
    niche_coverage(f$pts, f$env, f$box(10, 11, 30, 31))), "region")
})
