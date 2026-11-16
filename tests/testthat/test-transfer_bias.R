bias_fixture <- function() {
  sim <- simulate_enm_data(seed = 1)
  present <- sim$env
  future <- sim$env + 0.15
  names(future) <- names(present)
  # a climate envelope: 1 inside the range of the points in both variables
  fit <- function(p) {
    e <- terra::extract(present, as.matrix(p[, c("x", "y")]))
    lim <- apply(e, 2, range, na.rm = TRUE)
    inside <- function(r) {
      (r[[1]] >= lim[1, 1]) * (r[[1]] <= lim[2, 1]) *
        (r[[2]] >= lim[1, 2]) * (r[[2]] <= lim[2, 2])
    }
    out <- c(inside(present), inside(future))
    names(out) <- c("present", "future")
    out
  }
  list(sim = sim, present = present, future = future, fit = fit)
}

test_that("null occurrences are drawn uniformly, by effort and as a pattern", {
  f <- bias_fixture()
  u <- null_occurrences(f$present, f$sim$presence, n_sets = 4, seed = 2)
  expect_s3_class(u, "null_occurrences")
  expect_length(u, 4)
  cells <- terra::cellFromXY(f$present, as.matrix(f$sim$presence))
  expect_equal(nrow(u[[1]]), length(unique(cells)))
  expect_named(u[[1]], c("x", "y"))
  # distinct cells within a set, reproducible, different between sets
  cl <- terra::cellFromXY(f$present, as.matrix(u[[1]]))
  expect_false(anyDuplicated(cl) > 0)
  expect_identical(u, null_occurrences(f$present, f$sim$presence, n_sets = 4,
                                       seed = 2))
  expect_false(identical(u[[1]], u[[2]]))
  expect_output(print(u), "4 sets of")

  # effort that rises eastwards pulls the points eastwards
  b <- null_occurrences(f$present, n = 200, method = "bias",
                        weight = f$present[[1]]^3, n_sets = 3)
  un <- null_occurrences(f$present, n = 200, n_sets = 3)
  mx <- function(s) mean(vapply(s, function(d) mean(d$x), numeric(1)))
  expect_gt(mx(b), mx(un) + 1)
  # cells without effort are never drawn
  w <- f$present[[1]] > 0.5
  z <- null_occurrences(f$present, n = 100, method = "bias", weight = w,
                        n_sets = 2)
  expect_true(all(terra::extract(w, as.matrix(z[[1]]))[, 1] == 1))

  # the pattern keeps the distances between the points
  pts <- data.frame(lon = 18 + c(0, 0.3, 0.6, 0.2, 0.5),
                    lat = 52 + c(0, 0.1, 0.4, 0.5, 0.2))
  p <- null_occurrences(f$present, pts, method = "pattern", n_sets = 6,
                        seed = 3)
  expect_equal(attr(p, "method"), "pattern")
  km <- function(d) {
    sort(as.numeric(terra::distance(terra::vect(as.matrix(d),
                                                crs = "EPSG:4326"))))
  }
  d0 <- km(pts)
  whole <- which(vapply(p, function(d) {
    all(!is.na(terra::cellFromXY(f$present, as.matrix(d))))
  }, logical(1)) & attr(p, "redrawn") == 0)
  if (length(whole)) {
    expect_equal(km(p[[whole[1]]]), d0, tolerance = 0.03)
  }
  expect_false(isTRUE(all.equal(p[[1]]$x, pts$lon)))

  expect_error(null_occurrences(f$present), "points")
  expect_error(null_occurrences(f$present, n = 10, method = "bias"),
               "weight")
  expect_error(null_occurrences(f$present, n = 1e6), "fewer cells")
  expect_error(null_occurrences(f$present, n = 10, method = "pattern"),
               "pattern")
})

test_that("a change that every model predicts is called design-driven", {
  f <- bias_fixture()
  sets <- null_occurrences(f$present, f$sim$presence, n_sets = 29)
  # records spread over the whole area: the envelope of the real records is
  # that of random points
  spread <- null_occurrences(f$present, n = 100, n_sets = 1, seed = 99)[[1]]
  names(spread) <- c("lon", "lat")
  tb <- transfer_bias_test(fit = f$fit, points = spread, null_points = sets,
                           verbose = FALSE)
  expect_s3_class(tb, "transfer_bias")
  s <- tb$summary
  expect_equal(s$metric, c("mean_change", "share_declining"))
  m <- s[s$metric == "mean_change", ]
  expect_lt(m$null_mean, -0.1)            # decline for random points
  expect_lt(m$null_hi, 0)
  expect_equal(m$null_agreement, 1)
  expect_equal(m$verdict, "design-driven")
  expect_equal(m$signal, m$observed - m$null_mean)
  expect_output(print(tb), "does not depend on the")
  expect_output(print(tb), "IN SHORT")
  # maps: one stack per projection
  expect_named(tb$maps, "future")
  expect_equal(names(tb$maps$future),
               c("observed_change", "null_mean", "null_sd", "signal"))
  v <- terra::values(tb$maps$future)
  expect_equal(v[, "signal"], v[, "observed_change"] - v[, "null_mean"])
  expect_equal(mean(v[, "null_mean"]), m$null_mean)
  expect_equal(dim(tb$null), c(2L, 3L + 29L))
})

test_that("records with a narrow niche give a data-driven change", {
  f <- bias_fixture()
  sets <- null_occurrences(f$present, f$sim$presence, n_sets = 29)
  tb <- transfer_bias_test(fit = f$fit, points = f$sim$presence,
                           null_points = sets, threshold = 0.5,
                           verbose = FALSE)
  m <- tb$summary[tb$summary$metric == "mean_change", ]
  expect_equal(m$verdict, "data-driven")
  expect_lte(m$p_value, 2 / 30)
  expect_true("area_change" %in% tb$summary$metric)
  a <- tb$summary[tb$summary$metric == "area_change", ]
  expect_gte(a$observed, -1)
  expect_output(print(tb), "the records add")
  p <- plot_transfer_bias(tb)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot_transfer_bias(tb, "area_change"), "ggplot")
  expect_error(plot_transfer_bias(tb, scope = "region"), "No results")
})

test_that("predictions made elsewhere give the same result as `fit`", {
  f <- bias_fixture()
  sets <- null_occurrences(f$present, f$sim$presence, n_sets = 6)
  pts <- f$sim$presence
  names(pts) <- c("x", "y")
  obs <- f$fit(pts)
  nul <- lapply(sets, f$fit)
  region <- terra::as.polygons(terra::ext(14, 19, 49, 55), crs = "EPSG:4326")
  a <- transfer_bias_test(obs, nul, region = region)
  b <- transfer_bias_test(fit = f$fit, points = f$sim$presence,
                          null_points = sets, region = region,
                          verbose = FALSE)
  expect_equal(a$summary, b$summary)
  expect_setequal(a$summary$scope, c("area", "region"))
  expect_output(print(a), "\\[region\\]")
  # files work as well, and the baseline can be named
  d <- tempfile("null_models")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  files <- vapply(seq_along(nul), function(i) {
    fn <- file.path(d, paste0("null_", i, ".tif"))
    terra::writeRaster(nul[[i]], fn)
    fn
  }, character(1))
  c2 <- transfer_bias_test(obs, as.list(files), baseline = "present",
                           region = region, maps = FALSE)
  expect_equal(c2$summary, a$summary)
  expect_null(c2$maps)
  # p10 threshold from the points of each model
  e <- transfer_bias_test(obs, nul, points = f$sim$presence,
                          null_points = sets, threshold = "p10")
  expect_true("area_change" %in% e$summary$metric)

  expect_error(transfer_bias_test(obs), "observed")
  expect_error(transfer_bias_test(obs, nul[1]), "At least 2")
  expect_error(transfer_bias_test(obs[[1]], nul), "projection")
  expect_error(transfer_bias_test(obs, list(obs, obs[[1]])), "grid and the")
  expect_error(transfer_bias_test(obs, nul, baseline = "x"), "baseline")
  expect_error(transfer_bias_test(obs, nul, threshold = "p10"), "p10")
  expect_error(transfer_bias_test(fit = f$fit), "null_points")
})

test_that("environment_shift reports disappearing and novel conditions", {
  sim <- simulate_enm_data(seed = 1)
  warmer <- sim$env
  warmer[[1]] <- warmer[[1]] + 0.3
  es <- environment_shift(sim$env, list(warmer = warmer, same = sim$env),
                          trim = 0)
  expect_s3_class(es, "environment_shift")
  expect_equal(nrow(es), 6)
  w <- es[es$projection == "warmer", ]
  # env1 runs from 0 to 1: a shift of 0.3 removes and adds 30% of the cells
  expect_equal(w$disappearing[w$variable == "env1"], 0.3, tolerance = 0.02)
  expect_equal(w$novel[w$variable == "env1"], 0.3, tolerance = 0.02)
  expect_equal(w$disappearing[w$variable == "env2"], 0)
  expect_equal(w$disappearing[w$variable == "(any variable)"], 0.3,
               tolerance = 0.02)
  expect_gt(w$shift_sd[w$variable == "env1"], 0.9)
  s <- es[es$projection == "same", ]
  expect_true(all(s$disappearing == 0 & s$novel == 0))
  expect_output(print(es), "little change")
  expect_output(print(es), "predicted decline")
  # a single raster, and a missing layer
  one <- environment_shift(sim$env, warmer)
  expect_equal(unique(one$projection), "projection")
  expect_error(environment_shift(sim$env, warmer[[1]]), "lacks layer")
  expect_error(environment_shift(sim$env, warmer, trim = 0.5), "trim")
})
