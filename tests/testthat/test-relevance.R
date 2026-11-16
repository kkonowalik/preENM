# three patchy layers and records of a species that lives where the first
# one is high
relevance_data <- function(seed = 1, n = 150) {
  sim <- simulate_enm_data(seed = 1)
  env <- random_variables(sim$env, n = 3, range_km = 100, seed = seed)
  names(env) <- c("temperature", "rainfall", "soil")
  set.seed(seed)
  cells <- sample(terra::ncell(env), n,
                  prob = terra::values(env$temperature)[, 1]^6)
  pts <- as.data.frame(terra::xyFromCell(env, cells))
  names(pts) <- c("lon", "lat")
  list(env = env, points = pts)
}

test_that("ks_distance equals the two-sample Kolmogorov-Smirnov statistic", {
  set.seed(1)
  x <- rnorm(40, 0.4)
  y <- rnorm(500)
  expect_equal(ks_distance(x, y),
               unname(suppressWarnings(stats::ks.test(x, y))$statistic))
  # with ties
  x <- round(x, 1)
  y <- round(y, 1)
  expect_equal(ks_distance(x, y),
               unname(suppressWarnings(stats::ks.test(x, y))$statistic))
  expect_equal(ks_distance(c(5, 6), 1:4), 1)
})

test_that("random_variables keeps the values and the empty cells", {
  sim <- simulate_enm_data(seed = 1)
  env <- sim$env
  env[1:50] <- NA
  for (m in c("field", "shuffle")) {
    rnd <- random_variables(env, n = 3, method = m, range_km = 150, seed = 2)
    expect_equal(names(rnd), paste0("random_", 1:3))
    expect_equal(unname(attr(rnd, "source")), c("env1", "env2", "env1"))
    for (i in 1:3) {
      a <- terra::values(rnd[[i]])[, 1]
      b <- terra::values(env[[attr(rnd, "source")[[i]]]])[, 1]
      expect_equal(is.na(a), is.na(b))
      expect_equal(sort(a), sort(b))
      expect_false(isTRUE(all.equal(a, b)))
    }
  }
  expect_equal(terra::values(random_variables(env, 2, seed = 5)),
               terra::values(random_variables(env, 2, seed = 5)))
  expect_error(random_variables(env, layers = "nope"), "not found")
  expect_error(random_variables(env, range_km = -1), "range_km")
})

test_that("field layers are spatially structured, shuffled ones are not", {
  sim <- simulate_enm_data(seed = 1)
  neighbour_cor <- function(r) {
    m <- terra::as.matrix(r, wide = TRUE)
    stats::cor(as.vector(m[, -1]), as.vector(m[, -ncol(m)]))
  }
  f <- random_variables(sim$env, n = 1, range_km = 100, seed = 1)
  s <- random_variables(sim$env, n = 1, method = "shuffle", seed = 1)
  expect_gt(neighbour_cor(f), 0.9)
  expect_lt(abs(neighbour_cor(s)), 0.1)
  # smaller patches, weaker correlation between neighbours
  f2 <- random_variables(sim$env, n = 1, range_km = 30, seed = 1)
  expect_lt(neighbour_cor(f2), neighbour_cor(f))
  tf <- tempfile(fileext = ".tif")
  w <- random_variables(sim$env, n = 2, range_km = 100, seed = 1,
                        filename = tf)
  expect_true(file.exists(tf))
  expect_equal(terra::nlyr(terra::rast(tf)), 2)
})

test_that("variable_relevance finds the driver and not the other layers", {
  d <- relevance_data()
  for (null in c("field", "shift")) {
    res <- variable_relevance(d$points, d$env, null = null, n_null = 99)
    expect_s3_class(res, "variable_relevance")
    expect_setequal(res$layer, names(d$env))
    expect_equal(res$layer[1], "temperature")
    expect_true(res$relevant[1])
    expect_equal(res$p_value[1], 0.01)
    expect_gt(res$shift_sd[1], 0.5)
    expect_lt(res$breadth_ratio[1], 0.8)
    expect_true(all(res$p_value > 0 & res$p_value <= 1))
    expect_equal(res$relevant, res$p_value <= 0.05)
    expect_equal(attr(res, "settings")$null, null)
    expect_output(print(res), "Not distinguishable from random")
  }
  # reproducible
  a <- variable_relevance(d$points, d$env, n_null = 19, seed = 3)
  b <- variable_relevance(d$points, d$env, n_null = 19, seed = 3)
  expect_equal(a, b)
})

test_that("variable_relevance does not flag layers unrelated to the records", {
  d <- relevance_data()
  set.seed(9)
  p <- vapply(1:12, function(s) {
    e <- random_variables(d$env, n = 2, range_km = 100, seed = 50 + s)
    variable_relevance(d$points, e, null = "shift", n_null = 39,
                       seed = s)$p_value
  }, numeric(2))
  # 24 unrelated layers: about 5% expected below 0.05
  expect_lte(sum(p <= 0.05), 4)
  expect_gt(stats::median(p), 0.2)
})

test_that("a smooth random layer beats white noise, hence the null", {
  d <- relevance_data()
  # D of unrelated smooth layers against D of unrelated white noise
  sm <- random_variables(d$env, n = 6, range_km = 100, seed = 11)
  wn <- random_variables(d$env, n = 6, method = "shuffle", seed = 11)
  r1 <- variable_relevance(d$points, sm, null = "shift", n_null = 19)
  r2 <- variable_relevance(d$points, wn, null = "shift", n_null = 19)
  expect_gt(mean(r1$D), 1.5 * mean(r2$D))
})

test_that("variable_relevance checks its input", {
  d <- relevance_data()
  expect_error(variable_relevance(d$points[1:5, ], d$env), "Fewer than 10")
  # an area that is a small part of the raster: shifts leave it
  e <- d$env
  cells <- terra::cellFromXY(e, as.matrix(d$points))
  keep <- unique(unlist(lapply(cells, function(i) i + c(-1, 0, 1))))
  keep <- keep[keep >= 1 & keep <= terra::ncell(e)]
  e[setdiff(seq_len(terra::ncell(e)), keep)] <- NA
  expect_error(suppressWarnings(
    variable_relevance(d$points, e, null = "shift", n_null = 19)), "field")
})

test_that("compare_to_random applies the majority rule", {
  mk <- function(alg, run, good, weak, bad, rnd) {
    data.frame(algorithm = alg, run = run,
               variable = c("good", "weak", "bad", "random_1", "random_2"),
               importance = c(good, weak, bad, rnd))
  }
  imp <- rbind(
    mk("glm", 1, 50, 12, 2, c(5, 8)), mk("glm", 2, 45, 6, 3, c(4, 7)),
    mk("rf", 1, 40, 10, 1, c(6, 3)),  mk("rf", 2, 42, 9, 2, c(5, 4)),
    mk("gbm", 1, 60, 3, 0, c(4, 4)),  mk("gbm", 2, 55, 2, 9, c(6, 3)))
  res <- compare_to_random(imp)
  expect_equal(res$variable, c("good", "weak", "bad"))
  expect_equal(res$glm, c(1, 0.5, 0))
  expect_equal(res$rf, c(1, 1, 0))
  expect_equal(res$gbm, c(1, 0, 0.5))
  expect_equal(res$algorithms, c(3, 2, 1))
  expect_equal(res$keep, c(TRUE, TRUE, FALSE))
  expect_equal(res$mean_importance[1], mean(c(50, 45, 40, 42, 60, 55)))
  lim <- attr(res, "random")
  expect_equal(nrow(lim), 6)
  expect_equal(lim$threshold[lim$algorithm == "glm" & lim$run == "1"], 8)

  # stricter share of runs
  res2 <- compare_to_random(imp, min_runs = 1)
  expect_equal(res2$keep[res2$variable == "weak"], FALSE)
  # the median random variable is easier to beat than the best one
  res3 <- compare_to_random(imp, threshold = 0.5)
  expect_true(all(res3[, c("glm", "rf", "gbm")] >=
                    res[match(res3$variable, res$variable),
                        c("glm", "rf", "gbm")]))
  # ties with random do not count as better
  tie <- data.frame(variable = c("a", "random_1"), importance = c(5, 5))
  expect_false(compare_to_random(tie)$keep)
  # other column names and one algorithm
  one <- imp[imp$algorithm == "glm", ]
  names(one) <- c("alg", "rep", "var", "gain")
  r1 <- compare_to_random(one, variable = "var", value = "gain", run = "rep")
  expect_equal(r1$all, c(1, 0.5, 0))
  expect_equal(r1$keep, c(TRUE, TRUE, FALSE))
})

test_that("compare_to_random checks its input", {
  imp <- data.frame(variable = c("a", "b"), importance = 1:2)
  expect_error(compare_to_random(imp), "No random variables")
  expect_error(compare_to_random(imp, variable = "x"), "not found")
  imp$variable <- c("random_1", "random_2")
  expect_error(compare_to_random(imp), "Only random")
  imp <- data.frame(variable = c("a", "random_1", "a"),
                    importance = c(1, 2, 3), run = c(1, 1, 2))
  expect_error(compare_to_random(imp), "no random variable")
  expect_error(compare_to_random(imp, threshold = 2), "threshold")
})
