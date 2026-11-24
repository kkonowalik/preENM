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

test_that("record_ignorance: completeness, age, focal species, interpolation", {
  grid <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                      resolution = 1, crs = "EPSG:4326")
  # cell 1: 10 records, 5 species, 2 singletons -> completeness 1 - 2/10
  rec <- data.frame(lon = 0.5, lat = 0.5,
                    species = c("a", "a", "b", "b", "c", "c", "c", "d", "e",
                                "a"),
                    year = 2020)
  # cell 2: 10 old records of one species (completeness 1 - 0/10)
  rec <- rbind(rec, data.frame(lon = rep(1.5, 10), lat = 0.5,
                               species = "f", year = 1900))
  ig <- record_ignorance(rec, grid, min_records = 5, year = FALSE,
                         uncertainty = FALSE, radius_km = 10)
  r <- ig$records
  expect_equal(unique(r$completeness[r$cell == 1]), 0.8)
  expect_equal(unique(r$completeness[r$cell == 2]), 1)
  v <- terra::values(ig$raster)[, 1]
  # completeness only: cell 2 best (0), cell 1 worst (1), cell 3 empty and
  # farther than 10 km from any record: 1
  expect_equal(v, c(1, 0, 1))
  # with age the old cell loses
  ig2 <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                          radius_km = 10)
  v2 <- terra::values(ig2$raster)[, 1]
  expect_lt(v2[1], v2[2])
  expect_equal(ig2$settings$components, c("completeness", "age"))
  # focal species: a cell takes its best record of the species
  ig3 <- record_ignorance(rec, grid, species = "c", min_records = 5,
                          uncertainty = FALSE, radius_km = 10)
  r3 <- ig3$records
  expect_equal(unname(terra::values(ig3$raster)[1, 1]),
               min(r3$ignorance[r3$cell == 1 & r3$species == "c"]))
  # interpolation: within the radius the empty cell takes a weighted value
  ig4 <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                          radius_km = 500)
  v4 <- terra::values(ig4$raster)[, 1]
  expect_true(v4[3] < 1 && v4[3] >= min(v4[1:2]))
  # reliability by class of identifier
  rec$who <- rep(c("taxonomist", "amateur"), length.out = nrow(rec))
  ig5 <- record_ignorance(rec, grid, min_records = 5, year = FALSE,
                          uncertainty = FALSE, radius_km = 10,
                          reliability = c(taxonomist = 1, amateur = 0.75),
                          reliability_col = "who")
  expect_equal(sort(unique(ig5$records$reliability)), c(0.75, 1))
  expect_warning(record_ignorance(rec, grid, species = "zz",
                                  min_records = 5, radius_km = 10),
                 "No record of")
  expect_output(print(ig5), "from 20 records")
})

test_that("record_ignorance: spatio-environmental interpolation", {
  grid <- terra::rast(xmin = 0, xmax = 3, ymin = 0, ymax = 1,
                      resolution = 1, crs = "EPSG:4326")
  # well-known cell 1, poorly known cell 3, empty cell 2 in between
  rec <- rbind(data.frame(lon = 0.5, lat = 0.5,
                          species = rep(c("a", "b"), 5), year = 2020),
               data.frame(lon = 2.5, lat = 0.5,
                          species = paste0("s", 1:10), year = 2020))
  base <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                           year = FALSE, radius_km = 500)
  v0 <- terra::values(base$raster)[, 1]
  # geographically halfway: the mean of both cells
  expect_equal(v0[2], mean(v0[c(1, 3)]), tolerance = 1e-6)
  # environment: cell 2 like cell 1 -> its value moves towards cell 1
  env <- terra::rast(grid, vals = c(0, 0.1, 5))
  e1 <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                         year = FALSE, radius_km = 500, env = env)
  v1 <- terra::values(e1$raster)[, 1]
  expect_lt(abs(v1[2] - v1[1]), abs(v0[2] - v0[1]))
  # environment like cell 3 instead -> towards cell 3
  env3 <- terra::rast(grid, vals = c(0, 4.9, 5))
  e3 <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                         year = FALSE, radius_km = 500, env = env3)
  v3 <- terra::values(e3$raster)[, 1]
  expect_lt(abs(v3[2] - v3[3]), abs(v0[2] - v0[3]))
})

test_that("record_ignorance: focal-species cells carry half of the weight", {
  grid <- terra::rast(xmin = 0, xmax = 4, ymin = 0, ymax = 1,
                      resolution = 1, crs = "EPSG:4326")
  # cell 1: species "f" (focal); cell 2: other species; cells 3, 4 empty
  rec <- rbind(data.frame(lon = 0.5, lat = 0.5, species = "f",
                          year = 2000:2004),
               data.frame(lon = 1.5, lat = 0.5, species = paste0("o", 1:5),
                          year = 2020))
  ig <- record_ignorance(rec, grid, species = "f", min_records = 5,
                         uncertainty = FALSE, radius_km = 1000)
  v <- terra::values(ig$raster)[, 1]
  # cell 3 is nearer to cell 2, but the focal cell 1 has half the weight
  expect_equal(v[3], mean(v[1:2]), tolerance = 1e-6)
  ig0 <- record_ignorance(rec, grid, min_records = 5, uncertainty = FALSE,
                          radius_km = 1000)
  v0 <- terra::values(ig0$raster)[, 1]
  # without a focal species, inverse squared distance favours cell 2
  w <- c(1 / 2^2, 1 / 1^2)
  expect_equal(v0[3], sum(w * v0[1:2]) / sum(w), tolerance = 1e-3)
})
