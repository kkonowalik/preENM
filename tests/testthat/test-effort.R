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
