# Plausibility of response curves

grad <- seq(0, 30, length.out = 61)

test_that("monotonic, unimodal and noisy curves are classified", {
  curves <- rbind(
    data.frame(model = "m", variable = "up", value = grad,
               suitability = plogis((grad - 12) / 4)),
    data.frame(model = "m", variable = "down", value = grad,
               suitability = 1 - grad / 30),
    data.frame(model = "m", variable = "gauss", value = grad,
               suitability = exp(-(grad - 18)^2 / 30)))
  set.seed(1)
  noisy <- data.frame(model = "m", variable = "noisy", value = grad,
                      suitability = 0.5 + 0.2 * sin(grad / 5) +
                        rnorm(61, 0, 0.1))
  rc <- response_check(rbind(curves, noisy))
  s <- rc$shapes
  sh <- stats::setNames(s$shape, s$variable)
  expect_equal(unname(sh[c("up", "down", "gauss", "noisy")]),
               c("increasing", "decreasing", "unimodal", "irregular"))
  # optimum of the Gaussian at 18 (grid step 0.5)
  expect_equal(s$optimum[s$variable == "gauss"], 18)
  expect_equal(s$optimum_position[s$variable == "gauss"], 18 / 30)
  expect_equal(s$turning_points[s$variable == "gauss"], 1)
  expect_equal(s$wiggliness[s$variable == "down"], 1)
  expect_gt(s$wiggliness[s$variable == "noisy"], 3)
  expect_output(print(rc), "Implausible shapes: noisy")
})

test_that("flat, U-shaped and edges; tolerance ignores small noise", {
  curves <- rbind(
    data.frame(variable = "flat", value = grad, suitability = 0.5 +
                 0.01 * sin(grad)),
    data.frame(variable = "u", value = grad,
               suitability = 1 - exp(-(grad - 15)^2 / 40)),
    data.frame(variable = "edge", value = grad,
               suitability = exp(-(grad - 35)^2 / 200)),
    data.frame(variable = "small_noise", value = grad,
               suitability = exp(-(grad - 15)^2 / 40) +
                 0.01 * sin(grad * 3)))
  s <- response_check(curves)$shapes
  sh <- stats::setNames(s$shape, s$variable)
  expect_equal(unname(sh["flat"]), "flat")
  expect_equal(unname(sh["u"]), "U-shaped")
  expect_equal(unname(sh["small_noise"]), "unimodal")
  # optimum beyond the data: still rising at the upper edge only
  e <- s[s$variable == "edge", ]
  expect_true(e$rising_high)
  expect_false(e$rising_low)
  expect_equal(e$optimum_position, 1)
  # a unimodal curve inside the range rises at neither edge
  expect_false(any(unlist(s[s$variable == "small_noise",
                            c("rising_low", "rising_high")])))
})

test_that("records, calibrated range and projection", {
  curves <- data.frame(variable = "t", value = grad,
                       suitability = exp(-(grad - 35)^2 / 200))
  set.seed(2)
  occ <- data.frame(t = runif(300, 10, 30))
  bg <- data.frame(t = runif(1000, 0, 30))
  proj <- data.frame(t = c(rep(25, 70), rep(35, 30)))   # 30% beyond 30
  rc <- response_check(curves, occurrences = occ, background = bg,
                       projection = proj)
  s <- rc$shapes
  expect_equal(rc$limits$high, max(c(occ$t, bg$t)))
  expect_equal(s$beyond_high, 0.3)
  expect_equal(s$beyond_low, 0)
  expect_output(print(rc), "Part of the projection lies beyond")
  # a curve that runs against the records disagrees with them
  against <- data.frame(variable = "t", value = grad,
                        suitability = 1 - grad / 30)
  occ2 <- data.frame(t = rbeta(300, 5, 1) * 30)
  a <- response_check(against, occurrences = occ2, background = bg)
  expect_lt(a$shapes$agreement, 0)
  expect_output(print(a), "do not follow the records")
})

test_that("raw predictions are binned; models and runs; plot", {
  set.seed(3)
  d <- data.frame(temp = runif(20000, 0, 30), rain = runif(20000, 0, 30))
  d$pred <- exp(-(d$temp - 15)^2 / 40) * plogis((d$rain - 10) / 3)
  d$model <- "GLM"
  rc <- response_check(d, n_bins = 20)
  s <- rc$shapes
  expect_true(rc$settings$binned)
  expect_equal(s$shape[s$variable == "temp"], "unimodal")
  expect_equal(s$shape[s$variable == "rain"], "increasing")
  expect_equal(s$optimum[s$variable == "temp"], 15, tolerance = 0.1)
  curves <- data.frame(model = rep(c("A", "B"), each = 122),
                       run = rep(rep(1:2, each = 61), 2),
                       variable = "t", value = rep(grad, 4),
                       suitability = rep(exp(-(grad - 15)^2 / 40), 4))
  r2 <- response_check(curves)
  expect_equal(nrow(r2$shapes), 4)
  expect_s3_class(plot_response_check(r2), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_response_check(rc, file = f)
  expect_true(file.exists(f))
  expect_error(response_check(data.frame(a = 1)), "needs the columns")
})
