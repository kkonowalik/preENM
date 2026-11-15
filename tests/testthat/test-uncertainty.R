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
