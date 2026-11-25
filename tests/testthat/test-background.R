test_that("coverage increases with n and saturates; n stable across seeds", {
  sim <- simulate_enm_data(seed = 1)
  ns <- c(100, 250, 500, 1000, 2500, 4000)
  recs <- vapply(1:3, function(s) {
    bs <- background_sufficiency(sim$env, ns = ns, n_rep = 4,
                                 total_bins = 400, seed = s)
    cv <- bs$summary$coverage
    expect_true(all(diff(cv) >= -1e-9))          # monotonic
    expect_gt(cv[length(cv)], 0.99)              # saturates
    expect_lt(cv[1], 0.5)
    expect_true(bs$saturated)
    expect_true(all(diff(bs$summary$ks_max) < 0.02))
    bs$recommended_n
  }, numeric(1))
  expect_equal(length(unique(recs)), 1)
  expect_equal(recs[1], 2500)
})

test_that("a biased background covers less and does not converge", {
  sim <- simulate_enm_data(seed = 1)
  ef <- terra::init(sim$env[[1]], "x")
  ef <- (ef - terra::global(ef, "min")[[1]])^3
  b <- background_sufficiency(sim$env, ns = c(250, 1000, 2500), n_rep = 2,
                              total_bins = 400, seed = 1, method = "effort",
                              effort = ef, effort_floor = 0)
  u <- background_sufficiency(sim$env, ns = c(250, 1000, 2500), n_rep = 2,
                              total_bins = 400, seed = 1)
  expect_true(all(b$summary$coverage < u$summary$coverage))
  expect_gt(min(b$summary$ks_max), 0.2)
  expect_equal(b$settings$method, "effort")
  expect_output(print(b), "biased on purpose")
})

test_that("not saturated within the range is reported; print and plot", {
  sim <- simulate_enm_data(seed = 1)
  bs <- background_sufficiency(sim$env, ns = c(50, 100, 200), n_rep = 2,
                               total_bins = 400, seed = 1,
                               n_presences = 100)
  expect_false(bs$saturated)
  expect_output(print(bs), "has not saturated")
  expect_output(print(bs), "per presence")
  expect_s3_class(plot_background_sufficiency(bs), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_background_sufficiency(bs, file = f)
  expect_true(file.exists(f))
  expect_error(background_sufficiency(sim$env, ns = 100), "two")
})
