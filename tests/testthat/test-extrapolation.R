test_that("extrapolation_risk matches the definitions", {
  set.seed(1)
  ref <- cbind(a = runif(300), b = runif(300))
  r <- terra::rast(nrows = 6, ncols = 6, xmin = 0, xmax = 6, ymin = 0,
                   ymax = 6, nlyrs = 2)
  P <- cbind(a = runif(36, -0.5, 1.5), b = runif(36, -0.5, 1.5))
  P[1, ] <- NA
  terra::values(r) <- P
  names(r) <- c("a", "b")
  ex <- extrapolation_risk(ref, r, mop_share = 0.1)
  expect_equal(names(ex), c("nt1", "nt2", "mop", "n_outside", "most_novel"))
  v <- terra::values(ex)
  expect_true(all(is.na(v[1, ])))
  val <- function(i, k) unname(v[i, k])
  lo <- unname(apply(ref, 2, min))
  hi <- unname(apply(ref, 2, max))
  mu <- unname(colMeans(ref))
  sdv <- unname(apply(ref, 2, sd))
  S <- cov(ref)
  Rs <- scale(ref)
  for (i in 2:36) {
    p <- unname(P[i, ])
    under <- pmin((p - lo) / (hi - lo), 0)
    over <- pmin((hi - p) / (hi - lo), 0)
    expect_equal(val(i, "nt1"), sum(under + over))
    expect_equal(val(i, "n_outside"), sum(p < lo | p > hi))
    expect_equal(val(i, "nt2"), sqrt(mahalanobis(p, mu, S)) /
                   sqrt(max(mahalanobis(ref, mu, S))))
    d <- sqrt(colSums((t(Rs) - (p - mu) / sdv)^2))
    expect_equal(val(i, "mop"), mean(sort(d)[1:30]))
    worst <- unname(which.min(under + over))
    expect_equal(val(i, "most_novel"),
                 if (any(p < lo | p > hi)) worst else NA_real_)
  }
  expect_true(any(v[, "n_outside"] == 2, na.rm = TRUE))
  expect_true(any(v[, "n_outside"] == 0, na.rm = TRUE))
})

test_that("nt2 finds new combinations inside the ranges", {
  set.seed(2)
  a <- runif(500)
  ref <- cbind(a = a, b = a + rnorm(500, 0, 0.03))     # b follows a
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = 2)
  terra::values(r) <- rbind(c(0.5, 0.5), c(0.1, 0.9))  # usual, unusual
  names(r) <- c("a", "b")
  v <- terra::values(extrapolation_risk(ref, r))
  expect_equal(v[, "nt1"], c(0, 0))
  expect_lt(v[1, "nt2"], 1)
  expect_gt(v[2, "nt2"], 1)
  expect_gt(v[2, "mop"], v[1, "mop"])
})

test_that("extrapolation_risk takes rasters and checks its input", {
  sim <- simulate_enm_data(seed = 1)
  m <- terra::crop(sim$env, terra::ext(14, 21, 49, 55))
  f <- tempfile(fileext = ".tif")
  ex <- extrapolation_risk(m, sim$env, n_ref = 500, filename = f)
  expect_true(file.exists(f))
  inside <- terra::values(terra::crop(ex$n_outside, m))
  expect_true(all(inside == 0))                 # M itself is not novel
  expect_gt(max(terra::values(ex$n_outside)), 0)
  expect_error(extrapolation_risk(m, sim$env[[1]]), "lacks layer")
  expect_error(extrapolation_risk(m, sim$env, mop_share = 0), "mop_share")
})
