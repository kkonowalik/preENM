tab <- data.frame(
  model = c("DNN", "GLM", "GBM", "MARS", "GAM", "MAXNET", "ANN", "CTA", "RF",
            "SRE"),
  AUC = c(0.919, 0.918, 0.919, 0.916, 0.914, 0.918, 0.924, 0.901, 0.912,
          0.701),
  MAEbal = c(0.215, 0.218, 0.236, 0.215, 0.216, 0.332, 0.196, 0.203, 0.204,
             0.299),
  bias = c(-0.003, -0.041, -0.048, -0.039, -0.035, 0.101, -0.042, -0.041,
           -0.019, -0.183),
  TSS = c(0.72, 0.70, 0.71, 0.70, 0.69, 0.71, 0.74, 0.72, 0.57, 0),
  fossil_AUC_deep = c(0.773, 0.769, 0.757, 0.750, 0.762, 0.767, 0.709,
                      0.686, 0.675, 0.645),
  AUC_sd = 0.02)
crit <- c(AUC = "max", MAEbal = "min", fossil_AUC_deep = "max")

test_that("rule 'both' and the plot argument of select_models", {
  sel <- select_models(tab, min_metrics = c(AUC = 0.8, TSS = 0.6,
                                            fossil_AUC_deep = 0.75),
                       criteria = crit, tol = 0.01, plot = FALSE)
  expect_equal(sort(sel$model[sel$pareto]), c("ANN", "DNN", "GLM"))
  both <- select_models(tab, rule = "both", min_metrics = c(AUC = 0.8,
                        TSS = 0.6, fossil_AUC_deep = 0.75),
                        criteria = crit, tol = 0.01, plot = FALSE)
  expect_equal(sort(both$model[both$selected]), c("DNN", "GLM"))
  expect_warning(select_models(tab, rule = "both", min_metrics = c(AUC = 0.99),
                               criteria = crit, plot = FALSE), "No model")
})

test_that("zoom shows the region of the best models", {
  sel <- select_models(tab, criteria = crit, plot = FALSE)
  h <- zoom_hidden(sel, "MAEbal", "AUC")
  expect_setequal(h, c("MAXNET", "SRE"))
  p <- plot_selection(sel, "MAEbal", "AUC", zoom = TRUE)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$subtitle, "Outside this view: .*SRE")
  expect_match(p$labels$title, "zoom")
  pl <- plot_selection(sel, "MAEbal", "AUC",
                       zoom = list(x = c(0.19, 0.24), y = c(0.9, 0.93)))
  expect_match(pl$labels$subtitle, "MAXNET")
  # default axes: second and first criterion
  expect_match(plot_selection(sel)$labels$title, "AUC vs MAEbal")
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  out <- plot(sel)
  grDevices::dev.off()
  expect_named(out, c("all", "zoom"))
  # a selected model is drawn after (above) unselected ones
  d <- ggplot2::layer_data(plot_selection(sel, "MAEbal", "AUC"),
                           i = 3)
  expect_true(nrow(d) == nrow(tab))
})

test_that("plot_models: any metrics, plain tables, 2-D and 3-D", {
  sel <- select_models(tab, criteria = crit, plot = FALSE)
  p <- plot_models(sel, x = "bias", y = "AUC")
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$x, "closer to 0")
  pt <- plot_models(tab, x = "MAEbal", y = "AUC")       # plain data.frame
  expect_s3_class(pt, "ggplot")
  expect_error(plot_models(tab, x = "nope", y = "AUC"), "not found")
  tab$score <- runif(nrow(tab))
  expect_error(plot_models(tab[c("model", "score")], x = "score",
                           y = "score"), "directions")
  expect_s3_class(plot_models(tab, x = "score", y = "AUC",
                              directions = c(score = "max")), "ggplot")
  expect_equal(metric_direction("fossil_AUC_deep"), "max")
  expect_equal(metric_direction("yds_area"), "min")
  sim <- simulate_enm_data(seed = 1)
  ev <- evaluate_models(sim$models, sim$presence, sim$absence, n_boot = 20,
                        seed = 1)
  expect_s3_class(plot_models(ev, x = "MAEbal", y = "AUC", zoom = TRUE),
                  "ggplot")
  skip_if_not_installed("plotly")
  p3 <- plot_models(sel, x = "MAEbal", y = "AUC", z = "bias")
  expect_s3_class(p3, "plotly")
  b <- plotly::plotly_build(p3)
  expect_true(all(vapply(b$x$data, `[[`, "", "type") == "scatter3d"))
  expect_equal(b$x$layout$scene$xaxis$autorange, "reversed")   # MAEbal
  p2 <- plot_models(ev, interactive = TRUE)
  expect_s3_class(p2, "plotly")
  skip_if_not_installed("htmlwidgets")
  f <- tempfile(fileext = ".html")
  suppressMessages(plot_models(sel, x = "MAEbal", y = "AUC", z = "bias",
                               file = f, browse = FALSE))
  expect_true(file.exists(f))
})
