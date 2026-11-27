# A local mirror of the SoilGrids folder layout: <var>/<var>_<depth>_mean.vrt
fake_soilgrids <- function(values = c(phh2o = 55, soc = 200, nitrogen = 100)) {
  root <- tempfile("sg")
  for (v in names(values)) {
    dir.create(file.path(root, v), recursive = TRUE)
    r <- terra::rast(nrows = 40, ncols = 40, xmin = 10, xmax = 11, ymin = 50,
                     ymax = 51, crs = "EPSG:4326")
    terra::values(r) <- values[[v]]
    r[1:40] <- NA                                   # a row of no data
    tif <- file.path(root, v, paste0(v, "_0-5cm_mean.tif"))
    terra::writeRaster(r, tif, NAflag = -32768, datatype = "INT2S")
    sf::gdal_utils("buildvrt", tif,
                   file.path(root, v, paste0(v, "_0-5cm_mean.vrt")),
                   quiet = TRUE)
  }
  b <- terra::rast(nrows = 40, ncols = 40, xmin = 10, xmax = 11, ymin = 50,
                   ymax = 51, crs = "EPSG:4326")
  terra::values(b) <- 150
  bf <- file.path(root, "BDTICM.tif")
  terra::writeRaster(b, bf)
  list(root = root, bedrock = bf)
}

test_that("layers are converted to conventional units and C:N is added", {
  fk <- fake_soilgrids()
  d <- tempfile("soil")
  s <- suppressMessages(download_soilgrids(
    d, extent = c(10.2, 10.8, 50.2, 50.8),
    vars = c("phh2o", "soc", "nitrogen"), resolution = 0.05,
    base_url = fk$root, bedrock_url = fk$bedrock))
  expect_equal(names(s), c("phh2o_0_5cm", "soc_0_5cm", "nitrogen_0_5cm",
                           "cn_0_5cm", "bedrock_depth"))
  v <- terra::values(s)
  v <- v[stats::complete.cases(v), , drop = FALSE]
  expect_gt(nrow(v), 0)
  expect_equal(unique(v[, "phh2o_0_5cm"]), 5.5, tolerance = 1e-6)   # pH x 10
  expect_equal(unique(v[, "soc_0_5cm"]), 20, tolerance = 1e-6)      # dg/kg
  expect_equal(unique(v[, "nitrogen_0_5cm"]), 1, tolerance = 1e-6)  # cg/kg
  expect_equal(unique(v[, "cn_0_5cm"]), 20, tolerance = 1e-6)
  expect_equal(unique(v[, "bedrock_depth"]), 150, tolerance = 1e-6)
  expect_equal(terra::res(s), c(0.05, 0.05))
  expect_true(all(file.exists(file.path(d, paste0(names(s), ".tif")))))
  # existing files are kept (no new download)
  expect_silent(download_soilgrids(d, extent = c(10.2, 10.8, 50.2, 50.8),
                                   vars = c("phh2o", "soc", "nitrogen"),
                                   resolution = 0.05, base_url = fk$root,
                                   bedrock_url = fk$bedrock))
})

test_that("a template puts the layers on the climate grid by averaging", {
  fk <- fake_soilgrids(c(phh2o = 60))
  tpl <- terra::rast(nrows = 4, ncols = 5, xmin = 10.25, xmax = 10.75,
                     ymin = 50.2, ymax = 50.6, crs = "EPSG:4326")
  s <- suppressMessages(download_soilgrids(
    tempfile("soil"), vars = "phh2o", template = tpl, bedrock = FALSE,
    base_url = fk$root))
  expect_true(terra::compareGeom(s, tpl, stopOnError = FALSE))
  expect_equal(unique(terra::values(s, mat = FALSE)), 6)
})

test_that("arguments are checked and an unreachable service stops politely", {
  expect_error(download_soilgrids(tempfile(), c(0, 1, 0, 1), vars = "ph"),
               "Unknown SoilGrids layer")
  expect_error(download_soilgrids(tempfile(), c(0, 1, 0, 1),
                                  depth = "0-25cm"), "0-5cm")
  expect_error(download_soilgrids(tempfile()), "extent")
  skip_on_cran()
  expect_error(download_soilgrids(tempfile(), c(0, 1, 0, 1),
                                  base_url = "https://soilgrids.invalid/data"),
               "cannot be reached")
})

test_that("soil layers are registered in the index table", {
  tab <- climate_index_names("soil", table = TRUE)
  expect_true(all(c("phh2o", "soc", "nitrogen", "cn", "cec", "clay",
                    "bedrock_depth") %in% tab$index))
  expect_true(all(tab$source == "soilgrids"))
  expect_equal(climate_index_names(source = "soilgrids"), tab$index)
  # a layer for every property download_soilgrids() offers
  expect_true(all(soilgrids_props()$var %in% tab$index))
  expect_false(any(tab$index %in% suppressMessages(
    climate_index_names(source = "monthly"))))
  m <- terra::rast(nrows = 2, ncols = 2, nlyrs = 12, vals = 10,
                   crs = "EPSG:4326")
  expect_error(climate_indices(m + 5, m - 5, m * 5, indices = "soil"),
               "download_soilgrids")
})

test_that("the real service returns a small region", {
  skip_on_cran()
  skip_on_ci()
  skip_if_offline("files.isric.org")
  d <- tempfile("soil")
  s <- suppressMessages(download_soilgrids(
    d, extent = c(19.9, 20.0, 49.25, 49.3), vars = "phh2o", bedrock = FALSE))
  ph <- terra::values(s, mat = FALSE)
  expect_true(all(ph[!is.na(ph)] > 3 & ph[!is.na(ph)] < 9))
})
