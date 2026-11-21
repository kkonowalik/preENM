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

test_that("records outside the layers are reported; separate regional data", {
  f <- cov_fixture()
  region <- f$box(45, 55)
  env_m <- terra::mask(f$env, region)          # layers of the region only
  expect_warning(
    expect_message(bad <- niche_coverage(f$pts, env_m, region),
                   "outside the layers"),
    "trivially")
  expect_equal(bad$summary$advice, "cannot be judged")
  expect_gt(bad$summary$records_dropped, 0)
  expect_output(print(bad), "left out")
  # a national data set from the west of the region only
  natl <- f$pts[f$pts$lat > 45 & f$pts$lat < 55 & f$pts$lon < 10, ]
  d <- tempfile()
  nc <- niche_coverage(f$pts, f$env, region, region_points = natl,
                       file = file.path(d, "coverage"))
  expect_true(all(file.exists(file.path(d, c("coverage.png",
                                              "coverage.pdf")))))
  expect_equal(nc$summary$region_points, nrow(natl))
  base <- niche_coverage(f$pts, f$env, region)
  expect_gt(nc$summary$unfilling, base$summary$unfilling)
  expect_output(print(nc), "Regional records")
  p <- plot_niche_coverage(nc)
  expect_true(any(grepl("regional records", levels(p$data$group))))
})

test_that("every plotting function can save its plot", {
  f <- cov_fixture()
  nc <- niche_coverage(f$pts, f$env, f$box(45, 55))
  d <- tempfile()
  expect_message(plot_niche_coverage(nc, file = file.path(d, "a.pdf")),
                 "Saved")
  expect_true(file.exists(file.path(d, "a.pdf")))
  sim <- simulate_enm_data(seed = 1)
  ev <- evaluate_models(sim$models, sim$presence, sim$absence)
  expect_message(plot_pareto(ev, file = file.path(d, "pareto.png")), "Saved")
  sel <- select_models(as.data.frame(ev), plot = FALSE)
  expect_message(plot_selection(sel, file = file.path(d, "sel")), "Saved")
  expect_true(all(file.exists(file.path(d, c("sel.png", "sel.pdf")))))
  grDevices::pdf(NULL)
  expect_message(plot(sel, file = file.path(d, "both.pdf")), "Saved")
  grDevices::dev.off()
  expect_true(file.exists(file.path(d, "both.pdf")))
  expect_message(plot_models(sel, x = "MAE", y = "AUC",
                             file = file.path(d, "pm.jpeg")), "Saved")
  expect_error(plot_pareto(ev, file = file.path(d, "x.docx")), "Unsupported|Unknown")
})
