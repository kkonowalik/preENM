test_that("gbif_effort_sql builds a valid aggregation query", {
  q <- gbif_effort_sql(c(-32, 75, 34, 82), target = c("Silene acaulis",
                                                      "O'Brien"))
  expect_match(q, "FLOOR\\(decimalLatitude \\* 10\\) AS lat_bin")
  expect_match(q, "decimalLongitude >= -32 AND decimalLongitude < 75")
  expect_match(q, "'O''Brien'", fixed = TRUE)
  expect_match(q, "phylum = 'Tracheophyta'", fixed = TRUE)
  expect_match(gbif_effort_sql(c(0, 1, 0, 1), cell = 0.25), "\\* 4\\)")
  expect_match(gbif_effort_sql(c(0, 1, 0, 1)), "0 AS n_target")
  expect_error(gbif_effort_sql(c(1, 0, 0, 1)), "extent")
  expect_error(gbif_effort_sql(c(0, 1, 0, 1), cell = 0.3), "divide")
})

test_that("gbif_effort_rasters places counts in the right cells", {
  tab <- data.frame(lat_bin = c(450, 451, 999), lon_bin = c(100, 101, 0),
                    n_records = c(10, 3, 5), n_species = c(8, 2, 5))
  e <- gbif_effort_rasters(tab, c(9.5, 10.5, 44.5, 45.5))
  expect_equal(names(e), c("n_records", "n_species"))
  expect_equal(terra::extract(e, cbind(10.05, 45.05))$n_records, 10)
  expect_equal(terra::extract(e, cbind(10.15, 45.15))$n_species, 2)
  expect_equal(terra::global(e, "sum")[[1]], c(13, 10))  # outside dropped
  d <- tempfile("eff")
  f <- tempfile(fileext = ".csv")
  utils::write.csv(tab, f, row.names = FALSE)
  e2 <- gbif_effort_rasters(f, c(9.5, 10.5, 44.5, 45.5), path = d)
  expect_true(file.exists(file.path(d, "effort_n_species.tif")))
  expect_error(gbif_effort_rasters(data.frame(a = 1), c(0, 1, 0, 1)),
               "lat_bin")
})

test_that("download_gbif_effort asks for credentials and re-uses tables", {
  d <- tempfile("dl")
  dir.create(d)
  utils::write.csv(data.frame(lat_bin = 1, lon_bin = 1), file.path(d,
                   "effort_table.csv"), row.names = FALSE)
  expect_equal(download_gbif_effort("SELECT", d),
               file.path(d, "effort_table.csv"))
  skip_if_not_installed("rgbif", "3.8")
  withr_env <- Sys.getenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"),
                          unset = NA)
  on.exit(if (any(!is.na(withr_env))) {
    do.call(Sys.setenv, as.list(withr_env[!is.na(withr_env)]))
  })
  Sys.unsetenv(c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL"))
  expect_error(download_gbif_effort("SELECT", tempfile()), "Renviron")
})

test_that("ignorance_map is monotone decreasing in effort", {
  r <- terra::rast(nrows = 1, ncols = 6, xmin = 0, xmax = 6, ymin = 0,
                   ymax = 1, crs = "EPSG:4326")
  eff <- terra::rast(r, vals = c(0, 1, 5, 10, 50, 100))
  names(eff) <- "n_records"
  for (m in c("half", "log", "linear")) {
    v <- terra::values(ignorance_map(eff, method = m, half = 10)$raster)[, 1]
    expect_true(all(diff(v) < 0), info = m)
    expect_equal(v[1], 1, info = m)
  }
  # known values: O / (N + O); log and linear end at 0 for the maximum
  v <- terra::values(ignorance_map(eff, half = 10)$raster)[, 1]
  expect_equal(v, 10 / (c(0, 1, 5, 10, 50, 100) + 10))
  expect_equal(v[4], 0.5)
  vl <- terra::values(ignorance_map(eff, method = "log")$raster)[, 1]
  expect_equal(vl, 1 - log(c(0, 1, 5, 10, 50, 100) + 1) / log(101))
  # species observation index: records per species
  e2 <- c(eff, terra::rast(r, vals = c(0, 1, 5, 5, 5, 10)))
  names(e2) <- c("n_records", "n_species")
  vi <- terra::values(ignorance_map(e2, index = TRUE, half = 1)$raster)[, 1]
  expect_equal(vi, 1 / (c(0, 1, 1, 2, 10, 10) + 1))
})

test_that("ignorance_map annotates a suitability raster", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0,
                   ymax = 2, crs = "EPSG:4326")
  eff <- terra::rast(r, vals = c(0, 100, 0, 100))
  suit <- terra::rast(r, vals = c(0.9, 0.9, 0.1, 0.1))
  ig <- ignorance_map(eff, half = 1, suitability = suit)
  v <- terra::values(ig$raster)
  expect_equal(v[, "suitability_known"], c(NA, 0.9, NA, 0.1))
  expect_equal(as.vector(v[, "class"]), c(4, 3, 2, 1))
  expect_equal(ig$summary$classes$cells, c(1, 1, 1, 1))
  expect_output(print(ig), "50% lie where")
  expect_true(is.list(plot_ignorance(ig)))
  expect_s3_class(plot_ignorance(ig, type = "ignorance"), "ggplot")
  # effort within a radius: a zero cell next to recorded cells gains
  big <- terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 0.5, ymin = 45,
                     ymax = 45.5, crs = "EPSG:4326", vals = 0)
  big[13] <- 100
  s <- ignorance_map(big, half = 1, radius_km = 12)
  vs <- terra::values(s$raster)[, 1]
  expect_lt(vs[12], 1)
  expect_equal(vs[1], 1)
  expect_error(ignorance_map(eff, half = 0), "positive")
})
