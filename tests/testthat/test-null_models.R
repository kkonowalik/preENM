make_null_data <- function() {
  set.seed(1)
  n <- 4000
  effort <- runif(n)                     # where people record
  suit <- runif(n)                       # where the species lives
  pool <- data.frame(habitat = suit, bias = effort, noise = runif(n),
                     effort = effort)
  # presences: suitable places, recorded in proportion to effort
  pres <- sample(n, 150, prob = suit^3 * effort)
  abs <- sample(setdiff(seq_len(n), pres), 600)
  tab <- rbind(cbind(obs = 1, pool[pres, 1:3]), cbind(obs = 0, pool[abs, 1:3]))
  list(tab = tab, pool = pool)
}

test_that("a model of the sampling bias is not significant under effort", {
  d <- make_null_data()
  u <- null_model_test(d$tab, d$pool, n_null = 199, seed = 1)
  expect_equal(nrow(u), 6)
  expect_equal(unique(u$metric), c("AUC", "MAEbal"))
  auc <- function(x, m) x$p_value[x$model == m & x$metric == "AUC"]
  # uniform null: habitat and bias both look significant, noise does not
  expect_lt(auc(u, "habitat"), 0.05)
  expect_lt(auc(u, "bias"), 0.05)
  expect_gt(auc(u, "noise"), 0.05)
  expect_equal(u$null_mean[u$model == "noise" & u$metric == "AUC"], 0.5,
               tolerance = 0.05)
  # effort-weighted null: only the habitat model stays significant
  w <- null_model_test(d$tab, d$pool, weight = "effort", n_null = 199,
                       seed = 1)
  expect_lt(auc(w, "habitat"), 0.05)
  expect_gt(auc(w, "bias"), 0.05)
  expect_true(attr(w, "weighted"))
  # observed values are those of evaluate_models
  ev <- evaluate_models(data = d$tab)
  expect_equal(u$observed[u$metric == "AUC"],
               ev$AUC[match(c("habitat", "bias", "noise"), ev$model)])
  # reproducible, and p-values lie in (0, 1]
  expect_identical(u, null_model_test(d$tab, d$pool, n_null = 199, seed = 1))
  expect_true(all(u$p_value > 0 & u$p_value <= 1))
  expect_equal(min(u$p_value), 1 / 200)
})

test_that("null_model_test takes rasters and checks its input", {
  sim <- simulate_enm_data(seed = 1)
  tab <- extract_predictions(sim$models, sim$presence, sim$absence)
  x <- null_model_test(tab, sim$models, weight = sim$env[[1]], n_null = 19,
                       criteria = c(AUC = "max", bias = "zero"), seed = 1)
  expect_equal(nrow(x), 2 * terra::nlyr(sim$models))
  expect_true(all(c("observed", "null_mean", "null_95", "p_value") %in%
                    names(x)))
  expect_error(null_model_test(tab, data.frame(x = 1:3)), "lacks model")
  expect_error(null_model_test(tab, sim$models, weight = "w"),
               "Weight column")
  pool <- as.data.frame(terra::values(sim$models))[1:10, ]
  expect_error(null_model_test(tab, pool), "fewer usable rows")
})
