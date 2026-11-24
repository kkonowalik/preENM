# Range metrics

laea <- paste("+proj=laea +lat_0=52 +lon_0=10 +x_0=0 +y_0=0 +datum=WGS84",
              "+units=m +no_defs")
# a 30 x 30 km landscape of 1 km cells, with `blocks` of suitable cells
landscape <- function(blocks) {
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30000, ymin = 0,
                   ymax = 30000, crs = laea, vals = 0)
  for (b in blocks) {
    rc <- expand.grid(row = b$rows, col = b$cols)
    r[terra::cellFromRowCol(r, rc$row, rc$col)] <- 1
  }
  r
}

test_that("EOO, AOO, area and landscape metrics of a known geometry", {
  # one 10 x 10 km block in rows/cols 3-12
  r <- landscape(list(list(rows = 3:12, cols = 3:12)))
  m <- range_metrics(r, equal_area = laea, threat_km = 10)$metrics
  expect_equal(m$area_km2, 100)
  expect_equal(m$EOO_km2, 100)
  # 2 km grid placed to give the fewest cells: 5 x 5 cells
  expect_equal(m$AOO_km2, 100)
  expect_equal(m$locations, 1)
  expect_equal(m$n_patches, 1)
  expect_equal(m$largest_patch_index, 100 * 100 / 900)
  # perimeter 40 km = 40000 m over 900 km2 = 90000 ha
  expect_equal(m$edge_density, 40000 / 90000)
  # two 2 x 2 km blocks, 6 km apart (cols 1-2 and 9-10, rows 1-2):
  # hull = rectangle 10 x 2 km
  r2 <- landscape(list(list(rows = 1:2, cols = 1:2),
                       list(rows = 1:2, cols = 9:10)))
  m2 <- range_metrics(r2, equal_area = laea, threat_km = 10)$metrics
  expect_equal(m2$area_km2, 8)
  expect_equal(m2$EOO_km2, 20)
  expect_equal(m2$AOO_km2, 8)
  expect_equal(m2$n_patches, 2)
  # edge to edge 6 km
  expect_equal(m2$mean_nn_km, 6)
  # EOO is never smaller than AOO (IUCN): a single cell
  r3 <- landscape(list(list(rows = 5, cols = 5)))
  m3 <- range_metrics(r3, equal_area = laea)$metrics
  expect_equal(m3$AOO_km2, 4)
  expect_equal(m3$EOO_km2, 4)
})

test_that("AOO subdivides coarse cells; projection and CRS handling", {
  # 10 km cells: one suitable cell covers 25 AOO cells
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 30000, ymin = 0,
                   ymax = 30000, crs = laea, vals = c(1, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(range_metrics(r, equal_area = laea)$metrics$AOO_km2, 100)
  # longitude/latitude: computed in an equal-area projection, with a message
  g <- terra::rast(nrows = 10, ncols = 10, xmin = 10, xmax = 11, ymin = 50,
                   ymax = 51, crs = "EPSG:4326", vals = 1)
  expect_message(mg <- range_metrics(g)$metrics, "equal-area")
  expect_equal(mg$area_km2, terra::expanse(g, unit = "km")$area,
               tolerance = 1e-6)
  expect_equal(mg$AOO_km2, mg$area_km2, tolerance = 0.1)
  # no CRS: a warning
  terra::crs(g) <- ""
  expect_warning(suppressMessages(range_metrics(g)), "no CRS")
})

test_that("strict thresholds give smaller ranges; change; occupied", {
  r <- terra::rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30000, ymin = 0,
                   ymax = 30000, crs = laea)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- exp(-((xy[, 1] - 15000)^2 + (xy[, 2] - 15000)^2) / 8e7)
  th <- c(lenient = 0.2, strict = 0.7)
  m <- range_metrics(r, th, equal_area = laea)$metrics
  s <- m[m$rule == "strict", ]
  l <- m[m$rule == "lenient", ]
  for (v in c("area_km2", "EOO_km2", "AOO_km2", "locations")) {
    expect_lt(s[[v]], l[[v]])
  }
  expect_gt(s$mean_patch_suitability, l$mean_patch_suitability)
  # a shrinking future: negative change under both rules
  fut <- r * 0.8
  x <- c(r, fut)
  names(x) <- c("now", "later")
  rm <- range_metrics(x, th, equal_area = laea)
  ch <- rm$change[rm$change$metric == "area_km2", ]
  expect_true(all(ch$change < 0))
  expect_output(print(rm), "uncertain in both directions")
  # occupied: of two patches, only the one with a record is kept
  r2 <- landscape(list(list(rows = 1:2, cols = 1:2),
                       list(rows = 1:2, cols = 9:10)))
  occ <- data.frame(x = 1500, y = 28500)            # in the first block
  mo <- range_metrics(r2, equal_area = laea, occurrences = occ,
                      occupied = TRUE, crs = laea)$metrics
  expect_equal(mo$area_km2, 4)
  expect_equal(mo$n_patches, 1)
})

test_that("Red List thresholds and their sensitivity to the rule", {
  expect_equal(iucn_class(c(5, 50, 1000, 3000, NA), c(10, 500, 2000)),
               c("CR", "EN", "VU", "none", NA))
  r <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 60000, ymin = 0,
                   ymax = 60000, crs = laea)
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- exp(-((xy[, 1] - 30000)^2 + (xy[, 2] - 30000)^2) /
                            3e8)
  # AOO about 2700 km2 (no VU) vs about 300 km2 (EN)
  rm <- range_metrics(r, c(lenient = 0.1, strict = 0.8), equal_area = laea)
  expect_equal(rm$categories$B2_AOO, c("none", "EN"))
  expect_output(print(rm), "The threshold rule changes the Red List")
  expect_error(range_metrics(r), "must be binary")
  f <- write_range_metrics(rm, tempfile())
  expect_true(all(file.exists(f)))
  expect_s3_class(plot_range_metrics(rm), "ggplot")
  # no suitable area is not a category
  z <- range_metrics(r, c(never = 2), equal_area = laea)
  expect_equal(z$categories$B2_AOO, "no suitable area")
})
