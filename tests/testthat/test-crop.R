# a CHELSA-like file: unsigned 16-bit integers with scale and offset
make_layers <- function(dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 20, ymin = 40,
                   ymax = 60, crs = "EPSG:4326")
  terra::values(r) <- 2700:3099
  raw <- file.path(dir, "raw.tmp.tif")
  terra::writeRaster(r, raw, datatype = "INT2U", overwrite = TRUE)
  sf::gdal_utils("translate", raw, file.path(dir, "CHELSA_TraCE21k_bio01_-100_V1.0.tif"),
                 options = c("-a_scale", "0.1", "-a_offset", "-273.15"))
  unlink(raw)
  f <- terra::init(r, function(n) seq_len(n) / 7)            # float, no NoData
  terra::writeRaster(f, file.path(dir, "CHELSA_TraCE21k_bio12_-100_V1.0.tif"),
                     overwrite = TRUE)
  terra::writeRaster(f, file.path(dir, "CHELSA_TraCE21k_bio010_-100_V1.0.tif"),
                     overwrite = TRUE)
  far <- terra::rast(nrows = 5, ncols = 5, xmin = 100, xmax = 105, ymin = 0,
                     ymax = 5, crs = "EPSG:4326", vals = 1)
  sub <- file.path(dir, "orog")
  dir.create(sub, showWarnings = FALSE)
  terra::writeRaster(far, file.path(sub, "CHELSA_TraCE21k_dem_-100.tif"),
                     overwrite = TRUE)
  dir
}
poly <- sf::st_sfc(sf::st_polygon(list(rbind(c(2.3, 42.2), c(12.6, 42.2),
  c(12.6, 55.4), c(2.3, 55.4), c(2.3, 42.2)))), crs = 4326)

test_that("values, data type, scale/offset and grid are kept exactly", {
  src <- make_layers(tempfile("src"))
  out <- tempfile("out")
  res <- crop_rasters(src, poly, out, variables = c("bio01", "bio12"),
                      verbose = FALSE)
  expect_equal(nrow(res), 2)                         # bio010 not selected
  expect_true(all(res$status == "written"))
  a <- terra::rast(file.path(src, "CHELSA_TraCE21k_bio01_-100_V1.0.tif"))
  b <- terra::rast(file.path(out, "CHELSA_TraCE21k_bio01_-100_V1.0.tif"))
  expect_equal(terra::datatype(b), "INT2U")
  expect_equal(as.vector(terra::scoff(b)), c(0.1, -273.15))
  ref <- terra::crop(a, terra::vect(poly), snap = "out")
  expect_true(terra::compareGeom(b, ref, stopOnError = FALSE))
  v <- terra::values(b)
  expect_equal(v[!is.na(v)], terra::values(ref)[!is.na(v)])
  expect_true(any(is.na(v)))                         # masked edge cells
  # the centre rule: the row centred at 55.5 N lies outside the polygon
  top <- terra::values(b)[seq_len(terra::ncol(b))]
  expect_true(all(is.na(top)))
  # float file without NoData: NaN outside, exact values inside
  fl <- terra::rast(file.path(out, "CHELSA_TraCE21k_bio12_-100_V1.0.tif"))
  fv <- terra::values(fl)
  expect_true(any(is.na(fv)))
  expect_equal(fv[!is.na(fv)],
               terra::values(terra::crop(terra::rast(file.path(
                 src, "CHELSA_TraCE21k_bio12_-100_V1.0.tif")),
                 terra::vect(poly), snap = "out"))[!is.na(fv)])
})

test_that("re-runs skip finished files; overwrite, all_touched, mask = FALSE", {
  src <- make_layers(tempfile("src"))
  out <- tempfile("out")
  r1 <- crop_rasters(src, poly, out, verbose = FALSE)
  expect_equal(sum(r1$status == "written"), 3)
  expect_equal(sum(r1$status == "no overlap"), 1)    # the far-away dem
  expect_true(is.na(r1$output[r1$status == "no overlap"]))
  r2 <- crop_rasters(src, poly, out, verbose = FALSE)
  expect_true(all(r2$status[r2$status != "no overlap"] == "skipped (exists)"))
  r3 <- crop_rasters(src, poly, out, overwrite = TRUE, all_touched = TRUE,
                     variables = "bio01", verbose = FALSE)
  expect_equal(r3$status, "written")
  at <- terra::values(terra::rast(r3$output))
  expect_false(any(is.na(at)))                       # every cell touched
  out2 <- tempfile("out")
  r4 <- crop_rasters(src, poly, out2, mask = FALSE, variables = "bio01",
                     verbose = FALSE)
  expect_false(any(is.na(terra::values(terra::rast(r4$output)))))
})

test_that("M objects, GeoPackages and other CRSs work", {
  src <- make_layers(tempfile("src"))
  set.seed(1)
  pts <- data.frame(lon = runif(200, 4, 10), lat = runif(200, 44, 52))
  m <- delimit_m(pts, link_km = 100, thin_km = 20)
  out <- tempfile("out")
  res <- crop_rasters(src, m, out, variables = "bio01", verbose = FALSE)
  expect_equal(res$status, "written")
  g <- write_m_area(m, tempfile("m"), "M")[["gpkg"]]
  res2 <- crop_rasters(src, g, tempfile("out"), variables = "bio01",
                       verbose = FALSE)
  expect_equal(res2$status, "written")
  # a projected raster: the polygon is projected to it
  pr <- terra::project(terra::rast(file.path(src,
                                             "CHELSA_TraCE21k_bio12_-100_V1.0.tif")),
                       "EPSG:3035")
  pdir <- tempfile("proj")
  dir.create(pdir)
  terra::writeRaster(pr, file.path(pdir, "bio12_3035.tif"))
  res3 <- crop_rasters(pdir, poly, tempfile("out"), verbose = FALSE)
  expect_equal(res3$status, "written")
  expect_true(terra::same.crs(terra::rast(res3$output), pr))
})

test_that("safety checks", {
  src <- make_layers(tempfile("src"))
  expect_error(crop_rasters(src, poly, src, verbose = FALSE),
               "must not be a folder")
  expect_error(crop_rasters(file.path(src, "nope.tif"), poly, tempfile(),
                            verbose = FALSE), "Not found")
  expect_error(crop_rasters(src, poly, tempfile(), variables = "zzz",
                            verbose = FALSE), "No raster files")
  bad <- file.path(src, "broken.tif")
  writeLines("not a raster", bad)
  expect_warning(r <- crop_rasters(bad, poly, tempfile(), verbose = FALSE),
                 "failed")
  expect_match(r$status, "^failed")
  expect_equal(default_nodata("INT2S"), "-32768")
  expect_equal(default_nodata("FLT4S"), "nan")
})
