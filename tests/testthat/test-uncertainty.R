# a west-east gradient (flat from north to south) and a flat layer
make_gradient <- function() {
  r <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40000, ymin = 0,
                   ymax = 40000, crs = "EPSG:3035")
  g <- terra::init(r, "x") / 1000                  # 0.5 ... 39.5, 1 per km
  flat <- g * 0 + 7
  x <- c(g, flat)
  names(x) <- c("grad", "flat")
  x
}

test_that("the range within the uncertainty is measured cell by cell", {
  env <- make_gradient()
  occ <- data.frame(x = c(20500, 20500, 20500, 20500, 5500, 35500, 20500),
                    y = 20500,
                    unc = c(0, 400, 1200, 3000, 2000, NA, 60000))
  u <- uncertainty_environment(occ, env, uncertainty = "unc",
                               coords = c("x", "y"), points_crs = "EPSG:3035",
                               max_rel = 0.1, verbose = FALSE)
  expect_equal(nrow(u), 7)
  br <- attr(u, "breadth")
  expect_equal(unname(br["grad"]),
               diff(quantile(c(rep(20.5, 5), 5.5, 35.5), c(0.05, 0.95),
                             names = FALSE)))
  expect_equal(unname(br["flat"]), 0)
  # cells are 1 km: the point is in the cell centre, so a circle of 400 m
  # stays inside its cell, 1200 m reaches one cell on each side (range 2),
  # 3000 m three cells (range 6), 2000 m two cells (range 4)
  expect_equal(u$unc_cells[1:2], c(1L, 1L))
  expect_equal(u$unc_rel_grad[1:5] * br[["grad"]], c(0, 0, 2, 6, 4))
  expect_equal(u$unc_cells[3], 9L)
  # a flat layer never limits
  expect_true(all(u$unc_rel_flat[1:5] == 0))
  expect_equal(u$unc_limiting[3], "grad")
  expect_equal(u$unc_ok[1:5], u$unc_rel_grad[1:5] <= 0.1)
  expect_true(u$unc_ok[1] && !u$unc_ok[4])
  # no uncertainty: undecided; beyond max_radius_km: rejected, not examined
  expect_true(is.na(u$unc_ok[6]))
  expect_false(u$unc_ok[7])
  expect_true(is.na(u$unc_cells[7]))
  # an assumed uncertainty for records without one
  u2 <- uncertainty_environment(occ, env, uncertainty = "unc",
                                coords = c("x", "y"),
                                points_crs = "EPSG:3035",
                                na_uncertainty = 1200, verbose = FALSE)
  expect_equal(u2$unc_rel_grad[6] * br[["grad"]], 2)
})

test_that("a point off the cell centre reaches the nearer neighbour only", {
  env <- make_gradient()
  occ <- data.frame(x = 20900, y = 20500)          # 100 m from the east edge
  u <- uncertainty_environment(occ, env, uncertainty = 300,
                               coords = c("x", "y"), points_crs = "EPSG:3035",
                               verbose = FALSE)
  expect_equal(u$unc_cells, 2L)
  u <- uncertainty_environment(occ, env, uncertainty = 950,
                               coords = c("x", "y"), points_crs = "EPSG:3035",
                               verbose = FALSE)
  # the west side is reached (900 m), its corners are not (1030 m)
  expect_equal(u$unc_cells, 7L)
})

test_that("uncertainty_environment works in longitude/latitude", {
  sim <- simulate_enm_data(seed = 1, n_presence = 40)
  occ <- sim$presence
  occ$coordinateUncertaintyInMeters <- rep(c(10, 30000), 20)
  expect_message(u <- uncertainty_environment(occ, sim$env, max_rel = 0.2,
                                              max_radius_km = 40),
                 "usable")
  expect_true(all(u$unc_rel_max[occ$coordinateUncertaintyInMeters == 10] == 0))
  expect_true(all(u$unc_rel_max[occ$coordinateUncertaintyInMeters > 10] > 0))
  # a larger uncertainty never gives a smaller range
  occ$coordinateUncertaintyInMeters <- 15000
  u15 <- uncertainty_environment(occ, sim$env, verbose = FALSE,
                                 max_radius_km = 40)
  occ$coordinateUncertaintyInMeters <- 30000
  u30 <- uncertainty_environment(occ, sim$env, verbose = FALSE,
                                 max_radius_km = 40)
  expect_true(all(u30$unc_rel_max >= u15$unc_rel_max))
  expect_error(uncertainty_environment(sim$presence, sim$env),
               "No uncertainty column")
  expect_error(uncertainty_environment(occ, sim$env, uncertainty = -1),
               "Negative")
})

jitter_records <- function() {
  data.frame(lon = c(10, 10.5, 11, 11.5, 12, 12.5),
             lat = c(45, 45.5, 46, 46.5, 47, 47.5),
             unc = c(0, 5000, 20000, NA, 1000, 300000))
}
km_between <- function(lon1, lat1, lon2, lat2) {
  sqrt(((lon2 - lon1) * 111.32 * cos(lat1 * pi / 180))^2 +
         ((lat2 - lat1) * 111.32)^2)
}

test_that("jittered points lie within their radius; zero does not move", {
  j <- jitter_uncertainty(jitter_records(), n = 50, uncertainty = "unc",
                          seed = 1, verbose = FALSE)
  expect_s3_class(j, "jittered_records")
  expect_length(j$replicates, 50)
  all_d <- do.call(rbind, j$replicates)
  d <- km_between(all_d$lon_orig, all_d$lat_orig, all_d$lon, all_d$lat)
  rad <- jitter_records()$unc[all_d$record] / 1000
  moved <- all_d$record %in% c(2, 3, 5)
  expect_true(all(d[moved] <= rad[moved] * 1.001))
  # they do move, and fill the circle (uniform over the area: mean
  # distance 2/3 of the radius)
  expect_gt(mean(d[all_d$record == 3]), 0)
  expect_equal(mean(d[all_d$record == 3]) / 20, 2 / 3, tolerance = 0.1)
  # uncertainty 0, no uncertainty and radius above max_radius_km: fixed
  for (k in c(1, 4, 6)) {
    expect_true(all(all_d$lon[all_d$record == k] == jitter_records()$lon[k]))
    expect_true(all(all_d$lat[all_d$record == k] == jitter_records()$lat[k]))
  }
  expect_equal(j$records$status[c(1, 4, 6)],
               c("fixed: uncertainty 0", "fixed: no uncertainty",
                 "fixed: radius above max_radius_km"))
  # a default radius for records without one
  j2 <- jitter_uncertainty(jitter_records(), n = 5, uncertainty = "unc",
                           na_uncertainty = 2000, seed = 1, verbose = FALSE)
  expect_true(j2$records$moved[4])
  # the same seed gives the same replicates
  j3 <- jitter_uncertainty(jitter_records(), n = 50, uncertainty = "unc",
                           seed = 1, verbose = FALSE)
  expect_equal(j3$replicates[[7]]$lon, j$replicates[[7]]$lon)
})

test_that("jittered points stay inside M and in cells with values", {
  rec <- data.frame(lon = c(10, 10.02), lat = c(45, 45), unc = 30000)
  # M: a narrow strip; most of each circle lies outside it
  m <- sf::st_sfc(sf::st_polygon(list(rbind(c(9.9, 44.9), c(10.1, 44.9),
                                            c(10.1, 45.1), c(9.9, 45.1),
                                            c(9.9, 44.9)))), crs = 4326)
  j <- jitter_uncertainty(rec, n = 40, uncertainty = "unc", m = m,
                          seed = 2, verbose = FALSE)
  all_d <- do.call(rbind, j$replicates)
  pts <- sf::st_as_sf(all_d, coords = c("lon", "lat"), crs = 4326)
  expect_true(all(lengths(sf::st_intersects(pts, m)) > 0))
  expect_true(all(km_between(all_d$lon_orig, all_d$lat_orig, all_d$lon,
                             all_d$lat) <= 30.03))
  # cells without values are avoided
  r <- terra::rast(nrows = 40, ncols = 40, xmin = 9, xmax = 11, ymin = 44,
                   ymax = 46, vals = 1, crs = "EPSG:4326")
  r[terra::cellFromXY(r, cbind(10.05, 45.05))] <- NA
  r[terra::xFromCell(r, 1:terra::ncell(r)) > 10] <- NA   # east half empty
  j2 <- jitter_uncertainty(rec, r, n = 40, uncertainty = "unc", seed = 2,
                           verbose = FALSE)
  all2 <- do.call(rbind, j2$replicates)
  ok <- !all2$jitter_failed
  v <- terra::extract(r, cbind(all2$lon[ok], all2$lat[ok]))[, 1]
  expect_false(anyNA(v))
})

test_that("environmental spread: none on a constant layer, known on a ramp", {
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 9, xmax = 13,
                   ymin = 43, ymax = 47, crs = "EPSG:4326")
  terra::values(r) <- 1
  rec <- data.frame(lon = c(10, 11, 12), lat = 45, unc = 20000)
  j <- jitter_uncertainty(rec, r, n = 20, uncertainty = "unc", seed = 1,
                          verbose = FALSE)
  expect_true(all(j$records$rel_range_max == 0))
  expect_output(print(j), "almost the same dataset")
  # a steep east-west ramp: points 20 km apart differ, records vary
  ramp <- terra::init(r, "x")
  rec2 <- rbind(rec, data.frame(lon = c(9.5, 12.5), lat = 45, unc = 0))
  j2 <- jitter_uncertainty(rec2, ramp, n = 30, uncertainty = "unc",
                           seed = 1, verbose = FALSE)
  expect_true(all(j2$records$varies[1:3]))
  expect_equal(j2$records$rel_range_max[4:5], c(0, 0))
  expect_equal(j2$variables$share_above, 1)
  expect_output(print(j2), "worthwhile")
  # propagate = "variable": records in a uniform environment stay fixed
  both <- c(r, ramp)
  names(both) <- c("flat", "ramp")
  flat_rec <- data.frame(lon = 10, lat = 45, unc = 20000)
  jf <- jitter_uncertainty(flat_rec, r, n = 5, uncertainty = "unc",
                           propagate = "variable", seed = 1, verbose = FALSE)
  expect_false(jf$records$moved)
  expect_match(jf$records$status, "uniform")
})

test_that("uncertainty_environment() says keep, propagate or drop", {
  r <- terra::rast(nrows = 100, ncols = 100, xmin = 9, xmax = 13,
                   ymin = 43, ymax = 47, crs = "EPSG:4326")
  ramp <- terra::init(r, "x")
  rec <- data.frame(lon = c(10, 11, 12), lat = 45,
                    unc = c(100, 15000, 100000))
  u <- uncertainty_environment(rec, ramp, uncertainty = "unc",
                               verbose = FALSE)
  expect_equal(u$unc_action, c("keep", "propagate", "drop"))
})

test_that("coordinate_spread(): zero for identical replicates, matches sd", {
  sim <- simulate_enm_data(seed = 1)
  p <- sim$models[[1]]
  same <- coordinate_spread(c(p, p, p))
  expect_equal(max(terra::values(same$raster[["sd"]]), na.rm = TRUE), 0)
  reps <- c(p, p * 0.9, p * 0.8)
  cs <- coordinate_spread(reps, ensemble = ensemble_models(sim$models))
  v <- terra::values(c(p, p * 0.9, p * 0.8))
  expect_equal(as.vector(terra::values(cs$raster[["sd"]])),
               apply(v, 1, stats::sd), tolerance = 1e-6)
  expect_equal(as.vector(terra::values(cs$raster[["mean"]])), rowMeans(v),
               tolerance = 1e-6)
  expect_true(is.data.frame(cs$comparison))
  expect_output(print(cs), "IN SHORT")
  expect_error(coordinate_spread(p), "two")
  f <- tempfile(fileext = ".png")
  plot_coordinate_spread(cs, type = "spread", file = f)
  expect_true(file.exists(f))
})

test_that("plot_jitter() and write_jitter()", {
  sim <- simulate_enm_data(seed = 1, n_presence = 30)
  j <- jitter_uncertainty(sim$presence, sim$env, n = 3, uncertainty = 5000,
                          seed = 1, verbose = FALSE)
  expect_s3_class(plot_jitter(j), "ggplot")
  out <- write_jitter(j, tempfile())
  expect_length(out, 4)
  expect_equal(nrow(utils::read.csv(out[1])), 30)
})
