make_m <- function() {
  r <- terra::rast(nrows = 50, ncols = 100, xmin = 0, xmax = 100, ymin = 0,
                   ymax = 50, crs = "EPSG:3035")
  a <- terra::init(r, "x")                         # 0.5 ... 99.5
  b <- terra::init(r, "y")                         # 0.5 ... 49.5
  x <- c(a, b)
  names(x) <- c("a", "b")
  x
}

test_that("truncation is found at the end where the records pile up", {
  env <- make_m()
  set.seed(1)
  # on `a` the species lives at the high end of M; on `b` in the middle
  occ <- data.frame(x = runif(300, 80, 100), y = rnorm(300, 25, 3))
  occ <- occ[occ$y > 0 & occ$y < 50, ]
  tr <- niche_truncation(occ, env, coords = c("x", "y"),
                         points_crs = "EPSG:3035")
  v <- tr$variables
  expect_equal(v$layer, c("a", "b"))
  expect_equal(v$truncated_high, c(TRUE, FALSE))
  expect_equal(v$truncated_low, c(FALSE, FALSE))
  # uniform between 80 and 100: a quarter of the records in the top 5% of
  # M, five times what M offers
  expect_equal(v$edge_high[1], 5, tolerance = 0.25)
  expect_equal(v$edge_low, c(0, 0))
  expect_equal(v$gap_high[1], 0, tolerance = 0.01)
  expect_gt(v$gap_low[1], 0.75)
  expect_equal(c(v$m_min[1], v$m_max[1]), c(0.5, 99.5))
  expect_null(tr$projections)
  expect_output(print(tr), "Truncated: a \\(high end\\)")
})

test_that("projections beyond a truncated end are flagged", {
  env <- make_m()
  set.seed(1)
  occ <- data.frame(x = runif(300, 80, 100), y = runif(300, 15, 35))
  up <- env
  up$a <- up$a + 20                 # a fifth of the cells above M on `a`
  down <- env
  down$a <- down$a - 20             # below M, where the niche is not cut
  down$b <- down$b + 10             # above M on `b`, which is not truncated
  tr <- niche_truncation(occ, env, list(up = up, down = down),
                         coords = c("x", "y"), points_crs = "EPSG:3035")
  p <- tr$projections
  expect_equal(nrow(p), 4)
  pa <- p[p$projection == "up" & p$layer == "a", ]
  expect_equal(pa$above_m, 0.2, tolerance = 0.01)
  expect_equal(pa$below_m, 0)
  expect_true(pa$risk_high)
  expect_false(pa$risk_low)
  da <- p[p$projection == "down" & p$layer == "a", ]
  expect_equal(da$below_m, 0.2, tolerance = 0.01)
  expect_false(da$risk_low || da$risk_high)
  db <- p[p$projection == "down" & p$layer == "b", ]
  expect_equal(db$above_m, 0.2, tolerance = 0.01)
  expect_false(db$risk_high)
  # outside the records is at least as common as outside M
  expect_true(all(p$above_records >= p$above_m))
  expect_output(print(tr), "a in up")
  g <- plot_niche_truncation(tr)
  expect_s3_class(g, "ggplot")
  expect_silent(ggplot2::ggplot_build(g))
  expect_equal(sort(unique(tr$ranges$data)),
               c("M", "down", "records", "up"))
  expect_error(niche_truncation(occ, env, list(env), coords = c("x", "y"),
                                points_crs = "EPSG:3035"), "named list")
  expect_error(niche_truncation(occ, env, list(p = env[[1]]),
                                coords = c("x", "y"),
                                points_crs = "EPSG:3035"), "lacks layer")
  expect_error(niche_truncation(occ[1:5, ], env, coords = c("x", "y"),
                                points_crs = "EPSG:3035"), "Fewer than 10")
})

# A temperature gradient west-east and a north-south gradient; the species
# occupies temp 6-14 (lon 6-14) and lat 46-50 (not truncated)
remedy_fixture <- function() {
  r <- terra::rast(nrows = 60, ncols = 200, xmin = 0, xmax = 20, ymin = 45,
                   ymax = 51, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  set.seed(1)
  env <- c(terra::rast(r, vals = xy[, 1] + stats::rnorm(nrow(xy), 0, 0.1)),
           terra::rast(r, vals = xy[, 2] + stats::rnorm(nrow(xy), 0, 0.1)))
  names(env) <- c("temp", "north")
  pres <- data.frame(lon = stats::runif(300, 6, 14),
                     lat = stats::runif(300, 46, 50))
  list(r = r, xy = xy, env = env, pres = pres)
}
box <- function(x0, x1) {
  sf::st_sfc(sf::st_polygon(list(rbind(c(x0, 45), c(x1, 45), c(x1, 51),
                                       c(x0, 51), c(x0, 45)))), crs = 4326)
}

test_that("a truncated niche is extended towards the missing conditions", {
  f <- remedy_fixture()
  m <- box(2, 10)                         # M cuts the range at lon 10
  p <- f$pres[f$pres$lon < 10, ]
  tr <- niche_truncation(p, crop_to_m(f$env, m))
  rem <- truncation_remedy(tr, p, f$env, m)
  e <- rem$ends
  expect_equal(e$layer, "temp")           # the right variable ...
  expect_equal(e$end, "high")             # ... and direction
  expect_equal(rem$remedy, "extend")
  expect_lte(rem$extension$increase, 3)
  # the proposed extent contains the missing range of temperature
  ext <- terra::vect(sf::st_transform(rem$extension$extent, 4326))
  v <- terra::values(terra::mask(terra::crop(f$env[["temp"]], ext), ext))
  # every quarter of the missing band is represented in it
  brk <- seq(e$band_lo, e$band_hi, length.out = 5)
  part <- cut(v[!is.na(v) & v > e$band_lo & v <= e$band_hi], brk)
  expect_true(all(table(part) >= 5))
  # and grows only eastwards
  bb <- sf::st_bbox(sf::st_transform(rem$extension$extent, 4326))
  expect_gte(bb[["xmin"]], 2 - 1e-6)
  expect_gt(bb[["xmax"]], 10.5)
  expect_output(print(rem), "Extend the calibration area")
  expect_s3_class(plot_truncation_remedy(rem), "ggplot")
  out <- write_truncation_remedy(rem, tempfile())
  expect_true(all(file.exists(out)))
  expect_equal(nrow(utils::read.csv(out[["child"]])), nrow(p))
})

test_that("distant conditions point to donor regions; absent ones to accept", {
  f <- remedy_fixture()
  m <- box(2, 10)
  p <- f$pres[f$pres$lon < 10, ]
  tr <- niche_truncation(p, crop_to_m(f$env, m))
  # warm conditions only far east, beyond a cold gap
  e2 <- f$env
  v <- terra::values(e2[["temp"]])[, 1]
  x <- f$xy[, 1]
  v[x > 10 & x < 16] <- 1
  v[x >= 16] <- x[x >= 16] - 6
  terra::values(e2[["temp"]]) <- v
  rem <- truncation_remedy(tr, p, e2, m)
  expect_equal(rem$remedy, "donors")
  expect_gt(rem$donors$distance_km[1], 300)
  expect_match(rem$donors$supplies[1], "temp high")
  expect_s3_class(rem$donor_cells, "sf")
  expect_output(print(rem), "Integrate donor regions")
  # nowhere warmer than M: accept
  e3 <- f$env
  v3 <- terra::values(e3[["temp"]])[, 1]
  v3[x > 10] <- 2
  terra::values(e3[["temp"]]) <- v3
  rem3 <- truncation_remedy(tr, p, e3, m)
  expect_equal(rem3$remedy, "accept")
  expect_output(print(rem3), "Accept the truncation")
})

test_that("an untruncated niche needs no remedy", {
  f <- remedy_fixture()
  m <- box(0, 20)
  tr <- niche_truncation(f$pres, crop_to_m(f$env, m))
  rem <- truncation_remedy(tr, f$pres, f$env, m)
  expect_equal(rem$remedy, "none")
  expect_equal(nrow(rem$ends), 0)
  expect_output(print(rem), "no remedy is needed")
  expect_error(plot_truncation_remedy(rem), "nothing")
})
