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

test_that("envstrat background evens the environmental bins", {
  # one gradient where 90% of the area has nearly the same climate
  g <- terra::rast(nrows = 100, ncols = 100, xmin = 0, xmax = 10, ymin = 40,
                   ymax = 50, crs = "EPSG:4326")
  x <- terra::crds(g, na.rm = FALSE)[, 1]
  v1 <- ifelse(x < 9, x / 100, x - 8)
  env <- c(terra::setValues(g, v1),
           terra::setValues(g, v1 + stats::rnorm(length(v1), 0, 0.001)))
  names(env) <- c("a", "b")
  uni <- suppressMessages(sample_background(env, n = 1000, seed = 1))
  es <- suppressMessages(sample_background(env, n = 1000, seed = 1,
                                           method = "envstrat",
                                           env_axes = 1, env_total_bins = 10,
                                           env_min_bin = 1))
  expect_true("bin" %in% names(es))
  expect_gt(mean(es$lon > 9), 2 * mean(uni$lon > 9))
  info <- attr(es, "background")
  expect_equal(info$method, "envstrat")
  expect_equal(info$n_axes, 1)
  es2 <- suppressMessages(sample_background(env, n = 1000, seed = 1,
                                            method = "envstrat",
                                            env_axes = 1, env_total_bins = 10,
                                            env_min_bin = 1))
  expect_equal(es$cell, es2$cell)
  expect_equal(anyDuplicated(es$cell), 0)
})

test_that("effort and enveffort backgrounds follow recording effort", {
  sim <- simulate_enm_data(seed = 1)
  eff <- terra::init(sim$env[[1]], "x")
  eff <- (eff - terra::global(eff, "min")[[1]])^3
  ux <- suppressMessages(sample_background(sim$env, n = 300, seed = 2))
  ef <- suppressMessages(sample_background(sim$env, n = 300, seed = 2,
                                           method = "effort", effort = eff,
                                           effort_floor = 0))
  expect_gt(mean(ef$lon), mean(ux$lon))
  ee <- suppressMessages(sample_background(sim$env, n = 300, seed = 2,
                                           method = "enveffort", effort = eff,
                                           env_axes = 2,
                                           env_total_bins = 100))
  es <- suppressMessages(sample_background(sim$env, n = 300, seed = 2,
                                           method = "envstrat", env_axes = 2,
                                           env_total_bins = 100))
  expect_true("bin" %in% names(ee))
  # within the same bins, cells with more effort are preferred
  set.seed(4)
  eff2 <- terra::init(sim$env[[1]], function(n) stats::rexp(n)^2)
  ee2 <- suppressMessages(sample_background(sim$env, n = 300, seed = 2,
                                            method = "enveffort",
                                            effort = eff2, env_axes = 2,
                                            env_total_bins = 100))
  e_of <- function(b) mean(terra::extract(eff2, cbind(b$lon, b$lat))[[1]])
  expect_gt(e_of(ee2), 1.5 * e_of(es))
  expect_equal(length(unique(ee2$bin)), length(unique(es$bin)),
               tolerance = 0.1)
  expect_error(suppressMessages(sample_background(sim$env, method = "effort")),
               "needs an `effort`")
})
