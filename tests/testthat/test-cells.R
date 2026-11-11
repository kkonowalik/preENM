r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 40,
                 ymax = 50, vals = 1, crs = "EPSG:4326")
r2 <- r
r2[terra::cellFromXY(r2, cbind(5.5, 45.5))] <- NA
stack <- c(r, r2)

test_that("thin_by_cell keeps one record per cell and centres it", {
  pts <- data.frame(lon = c(0.2, 0.7, 3.3, 5.5, 20),
                    lat = c(40.2, 40.6, 43.3, 45.5, 45),
                    source = c("GBIF", "own", "GBIF", "GBIF", "GBIF"))
  th <- suppressMessages(thin_by_cell(pts, stack))
  expect_equal(nrow(th), 2)               # 0.2/0.7 share a cell; 5.5 NA; 20 out
  expect_equal(th$source[1], "GBIF")      # first in input order wins
  expect_equal(th$lon, c(0.5, 3.5))
  expect_equal(th$lat_orig, c(40.2, 43.3))
  thp <- suppressMessages(thin_by_cell(pts, stack, priority = c("own",
                                                                "GBIF")))
  expect_equal(thp$source[thp$cell == th$cell[1]], "own")
  expect_equal(thp$lon_orig[thp$source == "own"], 0.7)
  nc <- suppressMessages(thin_by_cell(pts, stack, centre = FALSE))
  expect_equal(nc$lon, c(0.2, 3.3))
  expect_error(suppressMessages(thin_by_cell(pts, r, priority = "own",
                                             priority_col = "x")), "not found")
})

test_that("thin_by_cell works per species and with projected rasters", {
  pts <- data.frame(species = c("a", "b", "a"), lon = c(0.2, 0.3, 0.4),
                    lat = c(40.2, 40.3, 40.4))
  th <- suppressMessages(thin_by_cell(pts, r))
  expect_equal(sort(th$species), c("a", "b"))
  rp <- terra::project(r, "EPSG:3035")
  thp <- suppressMessages(thin_by_cell(pts, rp))
  expect_true(all(abs(thp$lon) <= 180))    # output in lon/lat
  expect_equal(nrow(thp), 2)
})

test_that("sample_background draws valid, distinct, reproducible cells", {
  bg <- suppressMessages(sample_background(stack, n = 30, n_sets = 3,
                                           seed = 1))
  expect_equal(as.vector(table(bg$set)), c(30, 30, 30))
  expect_false(any(duplicated(bg[c("set", "cell")])))
  bad <- terra::cellFromXY(stack, cbind(5.5, 45.5))
  expect_false(bad %in% bg$cell)          # NA in one layer: never drawn
  bg2 <- suppressMessages(sample_background(stack, n = 30, n_sets = 3,
                                            seed = 1))
  expect_equal(bg, bg2)
  expect_false(identical(bg$cell[bg$set == 1], bg$cell[bg$set == 2]))
  # excluding presence cells, by thin_by_cell output or coordinates
  th <- suppressMessages(thin_by_cell(data.frame(lon = 0.5, lat = 40.5), r))
  b3 <- suppressMessages(sample_background(r, n = 99, exclude = th, seed = 2))
  expect_false(th$cell %in% b3$cell)
  b4 <- suppressMessages(sample_background(r, n = 99, seed = 2,
                                           exclude = data.frame(lon = 0.5,
                                                                lat = 40.5)))
  expect_equal(sort(b3$cell), sort(b4$cell))
  expect_warning(suppressMessages(sample_background(r, n = 500)),
                 "Fewer valid cells")
})

test_that("thin_by_cell handles merged tables with extra coordinate columns", {
  own <- data.frame(lon = c(0.4, 3.4), lat = c(40.4, 43.4))
  gb <- data.frame(decimalLongitude = c(0.45, 7.5), decimalLatitude = c(40.45, 47.5))
  m <- suppressMessages(merge_occurrences(own = own, GBIF = gb))
  th <- suppressMessages(thin_by_cell(m, r, priority = c("own", "GBIF")))
  expect_equal(nrow(th), 3)
  expect_equal(th$source, c("own", "own", "GBIF"))
  bad <- rbind(m[1, ], m)
  bad$lon[1] <- NA
  expect_message(thin_by_cell(bad, r), "Without coordinates: 1")
})
