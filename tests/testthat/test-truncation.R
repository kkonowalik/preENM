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
