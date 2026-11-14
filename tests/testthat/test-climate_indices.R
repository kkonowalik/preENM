make_months <- function(crs = "EPSG:4326") {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 30,
                   ymax = 38, crs = "EPSG:4326")
  season <- sin((1:12 - 4) / 12 * 2 * pi)          # warmest in July
  lat <- terra::init(r, "y")
  tmax <- terra::rast(lapply(1:12, function(m) 50 - lat + 10 * season[m]))
  tmin <- tmax - 12
  prec <- terra::rast(lapply(1:12, function(m) {
    20 + (38 - lat) * 10 * pmax(season[(m + 5) %% 12 + 1], 0)  # wet in Jan
  }))
  list(tmax = tmax, tmin = tmin, prec = prec)
}

test_that("radiation matches FAO-56", {
  expect_equal(round(extraterrestrial_radiation(-20 * pi / 180, j = 246), 1),
               32.2)
  expect_equal(extraterrestrial_radiation(75 * pi / 180, 12), 0)
})

test_that("climate_indices computes every index correctly", {
  d <- make_months()
  ci <- climate_indices(d$tmax, d$tmin, d$prec)
  expect_equal(names(ci), c("vpd", "pet", "cmd", "aridity", "heat_stress",
                            "sri", "synchronicity"))
  cell <- 1
  tx <- as.numeric(d$tmax[cell])
  tn <- as.numeric(d$tmin[cell])
  p <- as.numeric(d$prec[cell])
  tm <- (tx + tn) / 2
  lat <- terra::yFromCell(d$tmax, cell) * pi / 180
  es <- function(t) 0.6108 * exp(17.27 * t / (t + 237.3))
  expect_equal(ci$vpd[cell][[1]], mean((es(tx) + es(tn)) / 2 - es(tn)))
  days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  ra <- vapply(1:12, function(i) extraterrestrial_radiation(lat, i),
               numeric(1))
  pet <- pmax(0.0023 * 0.408 * ra * (tm + 17.8) * sqrt(tx - tn) * days, 0)
  expect_equal(ci$pet[cell][[1]], sum(pet), tolerance = 1e-6)
  expect_equal(ci$cmd[cell][[1]], sum(pmax(pet - p, 0)), tolerance = 1e-6)
  expect_equal(ci$aridity[cell][[1]], sum(p) / sum(pet), tolerance = 1e-6)
  expect_equal(ci$heat_stress[cell][[1]], sum(pmax(tx - 30, 0)))
  q <- function(m) ((m - 1 + 0:2) %% 12) + 1
  w <- which.max(vapply(1:12, function(m) mean(tm[q(m)]), 1))
  expect_equal(ci$sri[cell][[1]], sum(p[q(w)]) / sum(p))
  dd <- abs(which.max(tx) - which.max(p))
  expect_equal(ci$synchronicity[cell][[1]], (6 - min(dd, 12 - dd)) / 6)
  # warm season dry: rain six months from the heat
  expect_equal(ci$synchronicity[cell][[1]], 0)
  # subset of indices, written to files and re-read
  dir <- tempfile()
  ci2 <- climate_indices(d$tmax, d$tmin, d$prec, indices = c("pet", "sri"),
                         path = dir, prefix = "x_")
  expect_equal(names(ci2), c("pet", "sri"))
  expect_true(file.exists(file.path(dir, "x_pet.tif")))
  expect_equal(terra::values(ci2$pet), terra::values(ci$pet),
               tolerance = 1e-5)
  expect_error(climate_indices(d$tmax, d$tmin, d$prec, indices = "x"),
               "Unknown index")
  expect_error(climate_indices(d$tmax[[1:6]], d$tmin, d$prec), "12 monthly")
})

test_that("temperature units are checked", {
  d <- make_months()
  ref <- climate_indices(d$tmax, d$tmin, d$prec, indices = "heat_stress")
  expect_message(k <- climate_indices(d$tmax + 273.15, d$tmin + 273.15,
                                      d$prec, indices = "heat_stress"),
                 "Kelvin")
  expect_equal(terra::values(k), terra::values(ref), tolerance = 1e-6)
  k10 <- suppressMessages(climate_indices((d$tmax + 273.15) * 10,
                                          (d$tmin + 273.15) * 10, d$prec,
                                          indices = "heat_stress"))
  expect_equal(terra::values(k10), terra::values(ref), tolerance = 1e-5)
  expect_error(climate_indices(d$tmax / 10, d$tmin / 10, d$prec),
               "divided by 10")
  expect_s4_class(climate_indices(d$tmax / 10, d$tmin / 10, d$prec,
                                  temp_units = "C"), "SpatRaster")
  expect_error(climate_indices(d$tmax, d$tmin, d$prec - 1000), "Negative")
})

test_that("projected rasters give the same PET", {
  d <- make_months()
  ref <- climate_indices(d$tmax, d$tmin, d$prec, indices = "pet")
  pr <- lapply(d, function(x) terra::project(x, "EPSG:3035", res = 50000))
  ci <- climate_indices(pr$tmax, pr$tmin, pr$prec, indices = "pet")
  back <- terra::extract(ci, terra::project(terra::vect(cbind(4, 34),
                                                         crs = "EPSG:4326"),
                                             "EPSG:3035"))[[2]]
  expect_equal(back, terra::extract(ref, cbind(4, 34))[[1]], tolerance = 0.02)
})

test_that("monthly files are found and ordered", {
  dir <- tempfile()
  dir.create(dir)
  for (v in c("tas", "tasmax")) for (m in c(12, 1:11)) {
    file.create(file.path(dir, sprintf("CHELSA_%s_%02d_1981-2010_V.2.1.tif",
                                       v, m)))
  }
  f <- monthly_files(dir, "tas")
  expect_length(f, 12)
  expect_match(basename(f[1]), "_tas_01_")
  expect_match(basename(f[12]), "_tas_12_")
  expect_false(any(grepl("tasmax", f)))
  file.create(file.path(dir, "CHELSA_x_ssp126_pr_01_2071_2100_norm.tif"))
  expect_error(monthly_files(dir, "pr"), "Missing pr file")
  expect_error(monthly_files(dir, "tas", must_contain = "ssp585"), "Missing")
})

test_that("CHELSA file lookup", {
  expect_error(chelsa_monthly_urls("tas", "2071-2100"), "gcm")
  expect_equal(lonlat_window(c(1, 2, 3, 4)), c(1, 2, 3, 4))
  expect_error(lonlat_window(c(2, 1, 3, 4)), "xmin")
  v <- terra::vect(cbind(c(4e6, 4.1e6), c(3e6, 3.1e6)), crs = "EPSG:3035")
  w <- lonlat_window(v)
  expect_true(all(abs(w) <= 180))
  skip_on_cran()
  skip_on_ci()
  skip_if_offline("os.zhdk.cloud.switch.ch")
  u <- chelsa_monthly_urls("pr", "2071-2100", gcm = "UKESM1-0-LL",
                           ssp = "ssp370", months = 1:2)
  expect_equal(nrow(u), 2)
  expect_match(u$file[1], "ukesm1-0-ll.*ssp370_pr_01_2071_2100")
})
