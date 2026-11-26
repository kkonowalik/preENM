# A plane falling eastwards: 10 m cells, gradient 0.1 (5.71 degrees)
plane <- function(ymin = 0) {
  r <- terra::rast(nrows = 7, ncols = 9, xmin = 500000, xmax = 500090,
                   ymin = ymin, ymax = ymin + 70, crs = "EPSG:32633")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- 100 - 0.1 * (xy[, 1] - 500000)
  r
}

test_that("slope, aspect terms and TRI match a plane computed by hand", {
  ti <- topo_indices(plane(), indices = c("slope", "northness", "eastness",
                                          "tri"))
  m <- function(n) terra::as.matrix(ti[[n]], wide = TRUE)[2:6, 2:8]
  expect_equal(unique(round(as.vector(m("slope")), 6)),
               round(atan(0.1) * 180 / pi, 6))
  expect_equal(max(abs(m("northness"))), 0, tolerance = 1e-12)
  expect_equal(unique(round(as.vector(m("eastness")), 8)),
               round(sin(atan(0.1)), 8))
  # 3 neighbours 1 m higher, 3 lower, 2 equal
  expect_equal(unique(round(as.vector(m("tri")), 8)), round(sqrt(6), 8))
})

test_that("flat ground has no aspect: northness and eastness are 0", {
  f <- plane()
  terra::values(f) <- 50
  ti <- topo_indices(f, indices = c("northness", "eastness"))
  v <- terra::values(ti)
  v <- v[stats::complete.cases(v), , drop = FALSE]
  expect_gt(nrow(v), 0)
  expect_equal(max(abs(v)), 0, tolerance = 1e-12)
})

test_that("TPI and DEV are 0 on a plane and known for a single peak", {
  ti <- topo_indices(plane(), indices = c("tpi", "dev"), radius = 2)
  expect_equal(names(ti), c("tpi_2", "dev_2"))
  inner <- terra::as.matrix(ti$tpi_2, wide = TRUE)[3:5, 3:7]
  expect_equal(max(abs(inner)), 0, tolerance = 1e-9)
  # a 10 m peak on flat ground; radius 1 = the 4 rook neighbours
  r <- terra::rast(nrows = 9, ncols = 9, xmin = 0, xmax = 90, ymin = 0,
                   ymax = 90, crs = "EPSG:32633")
  terra::values(r) <- 0
  r[41] <- 10
  tp <- topo_indices(r, indices = "tpi", radius = c(1, 1.5))
  expect_equal(tp$tpi_1[41][[1]], 10)
  expect_equal(tp$tpi_1[40][[1]], -10 / 4)
  expect_equal(tp$tpi_1[31][[1]], 0)          # diagonal: outside radius 1
  expect_equal(tp$`tpi_1.5`[31][[1]], -10 / 8) # inside radius 1.5
  # DEV of a cell beside the peak: (0 - 10/4) / sd(c(10, 0, 0, 0))
  dv <- topo_indices(r, indices = "dev", radius = 1)
  expect_equal(dv$dev_1[40][[1]], -2.5 / stats::sd(c(10, 0, 0, 0)))
  expect_true(is.na(dv$dev_1[41][[1]]))       # flat window
  # radius in km: 20 m = 2 cells
  km <- topo_indices(plane(), indices = "tpi", radius = 0.02,
                     radius_unit = "km")
  expect_equal(names(km), "tpi_0.02km")
  expect_equal(terra::values(km, mat = FALSE),
               terra::values(ti$tpi_2, mat = FALSE))
})

test_that("TWI on an inclined plane equals ln(n * cell / tan(slope))", {
  ti <- topo_indices(plane(), indices = "twi")
  tw <- terra::as.matrix(ti$twi, wide = TRUE)[2:6, 2:8]
  # column j receives j cells of 100 m2 across a 10 m contour
  expect_equal(tw[1, ], log(10 * (2:8) / 0.1), tolerance = 1e-9)
  expect_equal(tw[3, ], tw[1, ])
  # a V-shaped valley gathers more water at its floor
  r <- terra::rast(nrows = 15, ncols = 11, xmin = 0, xmax = 110, ymin = 0,
                   ymax = 150, crs = "EPSG:32633")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  terra::values(r) <- abs(xy[, 1] - 55) * 0.2 + xy[, 2] * 0.05
  v <- terra::as.matrix(topo_indices(r, indices = "twi")$twi, wide = TRUE)
  expect_true(all(v[3:13, 6] > v[3:13, 3]))
  expect_error(topo_indices(r, indices = "twi", max_cells = 10), "memory")
})

test_that("HLI follows McCune & Keon (2002) and folds by hemisphere", {
  north <- plane(5500000)
  ti <- topo_indices(north, indices = "hli")
  ll <- terra::project(terra::xyFromCell(north, 32), from = "EPSG:32633",
                       to = "EPSG:4326")
  L <- ll[, 2] * pi / 180
  S <- atan(0.1)
  A <- abs(180 - abs(90 - 225)) * pi / 180   # east-facing
  hand <- exp(-1.467 + 1.582 * cos(L) * cos(S) -
                1.5 * cos(A) * sin(S) * sin(L) - 0.262 * sin(L) * sin(S) +
                0.607 * sin(A) * sin(S))
  expect_equal(ti$hli[32][[1]], hand, tolerance = 1e-4)
  # a south-west slope is warmer than a north-east one in the north
  mk <- function(dx, dy, crs = "EPSG:4326", ymin = 45) {
    r <- terra::rast(nrows = 7, ncols = 7, xmin = 10, xmax = 10.07,
                     ymin = ymin, ymax = ymin + 0.07, crs = crs)
    xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
    terra::values(r) <- 1000 + 300 * (dx * (xy[, 1] - 10) +
                                       dy * (xy[, 2] - ymin))
    topo_indices(r, indices = "hli")$hli[25][[1]]
  }
  expect_gt(mk(1, 1), mk(-1, -1))          # facing SW vs NE at 45 N
  expect_gt(mk(1, -1, ymin = -45), mk(-1, 1, ymin = -45))  # NW vs SE at 45 S
})

test_that("terrain indices are registered and written to disk", {
  tab <- climate_index_table()
  expect_setequal(tab$index[tab$group == "terrain"],
                  c("slope", "northness", "eastness", "tri", "tpi", "dev",
                    "twi", "hli"))
  expect_true(all(tab$source[tab$group == "terrain"] == "dem"))
  expect_true("twi" %in% suppressMessages(climate_index_names("terrain")))
  d <- tempfile("topo")
  on.exit(unlink(d, recursive = TRUE))
  ti <- topo_indices(plane(), indices = c("slope", "tpi"), radius = 2,
                     path = d, prefix = "x_")
  expect_true(all(file.exists(file.path(d, c("x_slope.tif", "x_tpi_2.tif")))))
  expect_error(topo_indices(plane(), "slope", path = d, prefix = "x_"),
               "exists")
  expect_error(topo_indices(plane(), "wetness"), "Unknown")
})
