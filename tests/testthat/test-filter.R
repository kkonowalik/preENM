offline <- c("country_mismatch", "coordinatecleaner", "locality_is_country")
# scenario tests below use few points with round coordinates: switch off the
# decimals rule and the cluster rule (tested on their own)
scen <- c(offline, "coarse_decimals", "small_isolated_cluster")

run_fake <- function(settings, ...) {
  x <- make_fake_gbif()
  L <- make_fake_layers()
  filter_occurrences(
    x, settings,
    filter_layers(range = L$range, land = L$land, climate = L$climate,
                  built = L$built),
    exclude_ids = "13", keep_ids = "19", verbose = FALSE, ...)
}

test_that("settings: presets, overrides, test selection", {
  s <- filter_settings("strict")
  l <- filter_settings("lenient")
  expect_s3_class(s, "filter_settings")
  expect_lt(s$max_coord_unc_m, l$max_coord_unc_m)
  expect_true("small_isolated_cluster" %in% s$tests)
  expect_false("small_isolated_cluster" %in% l$tests)
  expect_true(is.na(s$climate_limit))              # species-specific: off
  s2 <- filter_settings("strict", max_coord_unc_m = 500,
                        drop_tests = "duplicate")
  expect_equal(s2$max_coord_unc_m, 500)
  expect_false("duplicate" %in% s2$tests)
  s3 <- filter_settings("lenient", add_tests = "small_isolated_cluster")
  expect_equal(s3$tests, filter_tests()[filter_tests() %in% s3$tests])
  expect_error(filter_settings(max_unc = 5), "Unknown setting")
  expect_error(filter_settings(drop_tests = "nope"), "Unknown test")
  expect_output(print(s), "preset: strict")
  # accented words in the (ASCII-escaped) patterns still match
  expect_true(all(grepl(s$garden_regex, c("Jard\u00edn bot\u00e1nico",
                                          "ogr\u00f3d", "Tr\u00e4dg\u00e5rd"),
                        ignore.case = TRUE)))
  expect_false(grepl(s$garden_regex, "rocky slope above the lake"))
})

test_that("every rule removes its record, with the first reason kept", {
  st <- filter_settings("strict", climate_limit = 22, min_year = 1900,
                        drop_tests = offline)
  res <- run_fake(st)
  fate <- record_fate(res, as.character(1:19))
  expected <- c("invalid_coordinates", "duplicate", "grid_cell_centre",
                "gbif_issue", "basis_of_record", "establishment_means",
                "coord_uncertainty", "coarse_decimals", "garden_text",
                "settlement_text", "KEPT", "too_old", "manual_exclude",
                "KEPT", "outside_range", "at_sea", "built_up",
                "climate_limit", "KEPT")
  expect_equal(fate$fate[order(as.integer(fate$id))], expected)
  # the gridded atlas dataset is recognised as a lattice
  atlas <- res$data$datasetKey == "ds-atlas"
  expect_true(all(res$data$reason[atlas] == "lattice_dataset"))
  expect_true(res$datasets$lattice[res$datasets$dataset == "ds-atlas"])
  expect_equal(nrow(kept(res)), 60 + 3)       # good records + 11, 14, 19
  expect_equal(nrow(kept(res)) + nrow(removed(res)), nrow(res$data))
  expect_output(print(res), "kept")
})

test_that("the original duplicate is kept, its copy removed", {
  res <- run_fake(filter_settings("strict", drop_tests = offline))
  d <- res$data
  first <- d$gbifID == "1001"
  expect_true(is.na(d$reason[first]))
  expect_equal(d$reason[d$gbifID == "2"], "duplicate")
})

test_that("keep_ids override any test and keep the old reason", {
  st <- filter_settings("strict", max_coord_unc_m = 10, drop_tests = offline)
  res <- run_fake(st)
  r19 <- res$data[res$data$gbifID == "19", ]
  expect_true(is.na(r19$reason))
  expect_equal(r19$reason_overridden, "coord_uncertainty")
})

test_that("tests without their data are skipped, not silently passed", {
  x <- make_fake_gbif()
  res <- filter_occurrences(x, filter_settings("strict", drop_tests = offline),
                            verbose = FALSE)
  sk <- res$log$test[res$log$status == "skipped"]
  expect_true(all(c("outside_range", "at_sea", "built_up", "climate_limit")
                  %in% sk))
  y <- x[c("decimalLongitude", "decimalLatitude")]
  res2 <- filter_occurrences(y, filter_settings(drop_tests = offline),
                             verbose = FALSE)
  sk2 <- res2$log$test[res2$log$status == "skipped"]
  expect_true(all(c("gbif_issue", "coord_uncertainty", "lattice_dataset")
                  %in% sk2))
})

test_that("custom column names and numeric coordinates work", {
  d <- data.frame(x = c(14.123456, 14.123456, 14.5), y = c(47.1, 47.1, 47.2),
                  unc = c(10, 10, 5000))
  res <- filter_occurrences(d, filter_settings(drop_tests = scen),
                            cols = c(lon = "x", lat = "y", uncertainty = "unc"),
                            verbose = FALSE)
  expect_equal(res$data$reason, c(NA, "duplicate", "coord_uncertainty"))
  expect_equal(res$data$lon_decimals, c(6L, 6L, 1L))
  expect_error(filter_occurrences(d["unc"], verbose = FALSE), "Longitude")
  expect_error(filter_occurrences(d, cols = c(longitude = "x"),
                                  verbose = FALSE), "Unknown column role")
})

test_that("lenient and strict differ as documented", {
  d <- data.frame(decimalLongitude = c(14.11111, 14.22222),
                  decimalLatitude = c(47.11111, 47.22222),
                  coordinateUncertaintyInMeters = c(1500, 30))
  s <- filter_occurrences(d, filter_settings("strict", drop_tests = scen),
                          verbose = FALSE)
  l <- filter_occurrences(d, filter_settings("lenient", drop_tests = scen),
                          verbose = FALSE)
  expect_equal(s$data$reason[1], "coord_uncertainty")
  expect_true(is.na(l$data$reason[1]))
})

test_that("grid block and stacked datasets are recognised", {
  n <- 120
  d <- data.frame(
    datasetKey = rep(c("block", "stack"), c(25, n)),
    decimalLongitude = sprintf("%.5f", c(runif(25, 14, 15),
                                         rep(c(14.11111, 14.22222), n / 2))),
    decimalLatitude = sprintf("%.5f", c(runif(25, 47, 48),
                                        rep(c(47.11111, 47.22222), n / 2))),
    coordinateUncertaintyInMeters = c(rep("1800", 25), rep(NA, n)),
    gbifID = as.character(seq_len(25 + n)))
  st <- filter_settings("lenient", dedup_digits = NA, drop_tests = scen)
  res <- filter_occurrences(d, st, verbose = FALSE)
  expect_true(all(res$data$reason[1:25] == "grid_uncertainty_block"))
  expect_true(all(res$data$reason[26:(25 + n)] == "stacked_dataset"))
})

test_that("isolation tests: small clusters and climate outliers", {
  set.seed(3)
  d <- data.frame(decimalLongitude = c(runif(30, 14, 14.5), 20, 20.01),
                  decimalLatitude = c(runif(30, 47, 47.5), 50, 50.01))
  st <- filter_settings("strict", drop_tests = c(offline, "coarse_decimals"))
  res <- filter_occurrences(d, st, verbose = FALSE)
  expect_equal(res$data$reason[31:32], rep("small_isolated_cluster", 2))
  expect_true(all(is.na(res$data$reason[1:30])))
  # protected when inside the range polygon
  rng <- terra::vect("POLYGON ((19 49, 21 49, 21 51, 19 51, 19 49))",
                     crs = "EPSG:4326")
  res2 <- filter_occurrences(d, st, filter_layers(range = rng), verbose = FALSE)
  expect_true(is.na(res2$data$reason[31]))
  # warm and isolated
  clim <- terra::rast(xmin = 13, xmax = 21, ymin = 46, ymax = 51, res = 0.05,
                      vals = 10, crs = "EPSG:4326")
  clim[terra::cellFromXY(clim, cbind(c(20, 20.01), c(50, 50.01)))] <- 20
  st3 <- filter_settings("strict", isolation_climate = 14, min_support = 5,
                         drop_tests = scen)
  res3 <- filter_occurrences(d, st3, filter_layers(climate = clim),
                             verbose = FALSE)
  expect_equal(res3$data$reason[31:32], rep("isolated_climate_outlier", 2))
})

test_that("soft climate limit keeps records with a cold refuge nearby", {
  clim <- terra::rast(xmin = 13, xmax = 16, ymin = 46, ymax = 48, res = 0.01,
                      vals = 20, crs = "EPSG:4326")
  clim[terra::cellFromXY(clim, cbind(14.05, 47))] <- 5     # cold peak
  d <- data.frame(decimalLongitude = c(14, 15.5), decimalLatitude = c(47, 47))
  st <- filter_settings("strict", climate_soft_limit = 16, refuge_limit = 11,
                        refuge_radius_km = 10, drop_tests = scen)
  res <- filter_occurrences(d, st, filter_layers(climate = clim),
                            verbose = FALSE)
  expect_equal(res$data$reason, c(NA, "climate_soft_limit"))
  expect_equal(res$data$climate_near, c(5, 20))
  # direction "below": too cold is bad, a warm refuge saves the record
  st2 <- filter_settings("strict", climate_direction = "below",
                         climate_limit = 0, drop_tests = scen)
  clim2 <- terra::init(clim, -5)
  res2 <- filter_occurrences(d, st2, filter_layers(climate = clim2),
                             verbose = FALSE)
  expect_equal(res2$data$reason, rep("climate_limit", 2))
})

test_that("water raster fallback flags records on water", {
  w <- terra::rast(xmin = 13, xmax = 16, ymin = 46, ymax = 48, res = 0.01,
                   vals = 0, crs = "EPSG:4326")
  w[terra::cellFromXY(w, cbind(15, 47))] <- 1      # open water
  w[terra::cellFromXY(w, cbind(14, 47))] <- 1
  w[terra::cellFromXY(w, cbind(14.01, 47))] <- 0.5  # lake shore nearby
  d <- data.frame(decimalLongitude = c(15, 14, 14.5),
                  decimalLatitude = c(47, 47, 47),
                  coordinateUncertaintyInMeters = c(10, 2000, 10))
  # lenient: uncertainty up to 2000 m, so the 2000 m record reaches the test;
  # its larger uncertainty widens the window to 3x3 cells, which reach land
  res <- filter_occurrences(d, filter_settings("lenient", drop_tests = scen),
                            filter_layers(water = w), verbose = FALSE)
  expect_equal(res$data$reason, c("on_water", NA, NA))
  expect_equal(res$data$land_near_pct, c(0, 100, 100))
})

test_that("GHS and cropland layers", {
  r <- terra::rast(xmin = 13, xmax = 16, ymin = 46, ymax = 48, res = 0.01,
                   vals = 0, crs = "EPSG:4326")
  gmean <- r
  gmax <- r
  crop <- r
  gmean[terra::cellFromXY(r, cbind(14, 47))] <- 300   # 3 % built
  gmax[terra::cellFromXY(r, cbind(14.5, 47))] <- 5000
  crop[terra::cellFromXY(r, cbind(15, 47))] <- 0.9
  d <- data.frame(decimalLongitude = c(14, 14.5, 15, 15.5), decimalLatitude = 47)
  res <- filter_occurrences(d, filter_settings("strict", drop_tests = scen),
                            filter_layers(ghs_mean = gmean, ghs_max = gmax,
                                          cropland = crop), verbose = FALSE)
  expect_equal(res$data$reason, c("ghs_dense", "ghs_building_nearby",
                                  "cropland", NA))
  expect_equal(res$data$ghs_mean_pct[1], 3)
})

test_that("distance to land and projection choice", {
  land <- terra::vect("POLYGON ((13.3 46, 18 46, 18 49, 13.3 49, 13.3 46))",
                      crs = "EPSG:4326")
  d <- distance_to_land(c(13.1, 14, 11), c(47.5, 47.5, 47.5), land)
  expect_equal(d[1], 0.2 * 111.32 * cos(47.5 * pi / 180) * 1000,
               tolerance = 0.01)
  expect_equal(d[2], 0)
  expect_equal(d[3], Inf)
  circ <- terra::vect("POLYGON ((-179 60, 179 60, 179 80, -179 80, -179 60))",
                      crs = "EPSG:4326")
  expect_match(buffer_crs(c(-170, 170), c(65, 70), circ, "auto"),
               "lat_0=90")
  small <- terra::vect("POLYGON ((10 45, 12 45, 12 47, 10 47, 10 45))",
                       crs = "EPSG:4326")
  expect_match(buffer_crs(c(10, 12), c(45, 47), small, "auto"), "lat_0=46")
})

test_that("reports, sensitivity, fate and files", {
  res <- run_fake(filter_settings("strict", climate_limit = 22,
                                  drop_tests = offline))
  rep <- filter_report(res)
  expect_equal(sum(rep$n), nrow(res$data))
  expect_equal(rep$outcome[1], "KEPT")
  byc <- filter_report(res, by = "countryCode")
  expect_true("KEPT" %in% names(byc))
  sens <- filter_sensitivity(res, "climate", c(15, 30))
  expect_equal(sens$removed, c(1, 0))     # the 25 degree record only
  expect_error(filter_sensitivity(res, "nope", 1), "not found")
  td <- tempfile()
  paths <- write_filter_results(res, td, "dryas")
  expect_true(all(file.exists(paths)))
  expect_equal(nrow(utils::read.csv(paths[["kept"]])), nrow(kept(res)))
  expect_message(record_fate(res, c("1", "999")), "Not found")
  expect_error(kept(data.frame()), "filter_occurrences")
})

test_that("CoordinateCleaner tests (when installed)", {
  skip_if_not_installed("CoordinateCleaner")
  skip_if_not_installed("countrycode")
  skip_if_not_installed("rnaturalearthdata")
  skip_on_cran()
  st <- filter_settings("strict", cc_tests = c("zeros", "equal", "gbif"))
  x <- make_fake_gbif()
  extra <- x[1, ]
  extra$gbifID <- "77"
  extra$decimalLongitude <- "15.00000"
  extra$decimalLatitude <- "15.00000"
  extra$countryCode <- NA
  res <- filter_occurrences(rbind(x, extra), st, keep_ids = "19",
                            verbose = FALSE)
  f <- record_fate(res, c("11", "14", "19", "77"))
  expect_equal(f$fate[f$id == "14"], "country_mismatch")    # stated IT
  expect_equal(f$fate[f$id == "11"], "locality_is_country")
  expect_equal(f$fate[f$id == "77"], "coordinatecleaner")   # lon == lat
  expect_equal(f$reason_overridden[f$id == "19"], "country_mismatch")
})

test_that("locality names in other languages are recognised", {
  skip_if_not_installed("countrycode")
  nm <- country_names()
  expect_true(all(normalise_text(c("Osterreich", "Norway", "Polska")) %in%
                    nm))
  expect_false("rocky slope" %in% nm)
})
