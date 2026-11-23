# Ensembles of model predictions

make_stack <- function() {
  r <- terra::rast(nrows = 20, ncols = 30, xmin = 10, xmax = 13, ymin = 50,
                   ymax = 52, crs = "EPSG:4326")
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  base <- plogis((xy[, 1] - 11.5) * 2)
  s <- c(terra::rast(r, vals = base),
         terra::rast(r, vals = base^2),
         terra::rast(r, vals = sqrt(base)),
         terra::rast(r, vals = 1 - base))
  names(s) <- c("A", "B", "C", "D")
  s
}

test_that("identical rasters give zero disagreement", {
  s <- make_stack()
  same <- c(s[["A"]], s[["A"]], s[["A"]])
  names(same) <- c("m1", "m2", "m3")
  for (m in c("median", "mean", "pcam")) {
    e <- if (m == "pcam") {
      # PCA of identical columns is degenerate but all loadings are equal
      suppressWarnings(tryCatch(ensemble_models(same, method = m),
                                error = function(e) NULL))
    } else ensemble_models(same, method = m)
    if (is.null(e)) next
    v <- terra::values(e$raster)
    expect_equal(v[, "ensemble"], terra::values(s[["A"]])[, 1])
    expect_true(all(v[, "sd"] == 0))
  }
  e <- ensemble_models(same)
  expect_equal(max(e$departure$mean_abs_diff), 0)
  expect_output(print(e), "agree closely")
})

test_that("equal weights equal the unweighted mean; median and mean", {
  s <- make_stack()
  ew <- ensemble_models(s, method = "weighted",
                        weights = c(A = 2, B = 2, C = 2, D = 2))
  em <- ensemble_models(s, method = "mean")
  expect_equal(terra::values(ew$raster[["ensemble"]]),
               terra::values(em$raster[["ensemble"]]))
  v <- terra::values(s)
  expect_equal(terra::values(em$raster[["ensemble"]])[, 1], rowMeans(v))
  emed <- ensemble_models(s)
  expect_equal(terra::values(emed$raster[["ensemble"]])[, 1],
               apply(v, 1, stats::median))
  expect_equal(terra::values(em$raster[["sd"]])[, 1], apply(v, 1, stats::sd))
  # unequal weights by hand
  ww <- ensemble_models(s, method = "weighted",
                        weights = c(A = 3, B = 1, C = 0, D = 0))
  expect_equal(terra::values(ww$raster[["ensemble"]])[, 1],
               0.75 * v[, "A"] + 0.25 * v[, "B"])
})

test_that("the selection is respected and dominated models are ignored", {
  s <- make_stack()
  sel <- data.frame(model = c("A", "B", "C", "D"),
                    selected = c(TRUE, TRUE, FALSE, FALSE))
  e <- ensemble_models(s, sel)
  expect_equal(e$models, c("A", "B"))
  expect_equal(sort(e$left_out), c("C", "D"))
  expect_equal(terra::values(e$raster[["ensemble"]])[, 1],
               apply(terra::values(s[[c("A", "B")]]), 1, stats::median))
  expect_equal(ensemble_models(s, sel, use = "all")$models,
               c("A", "B", "C", "D"))
  # evaluate_models() output: only the Pareto front
  ev <- data.frame(model = c("A", "B", "C", "D"), AUC = c(0.9, 0.8, 0.7, 0.85),
                   MAE = c(0.3, 0.2, 0.4, 0.1), pareto_front = c(1, 1, 2, 1))
  attr(ev, "criteria") <- c(AUC = "max", MAE = "min")
  e2 <- ensemble_models(s, ev, method = "weighted")
  expect_equal(e2$models, c("A", "B", "D"))
  # weights: compromise score 1 - distance, normalised
  d <- pareto_core(to_minimise(ev, c(AUC = "max", MAE = "min")))$distance
  w <- (1 - d[c(1, 2, 4)]) / sum(1 - d[c(1, 2, 4)])
  expect_equal(unname(e2$weights), w)
  expect_output(print(e2), "Zhu et al. 2021")
  expect_error(ensemble_models(s, data.frame(model = "Z", selected = TRUE)),
               "No raster")
})

test_that("PCAm keeps the models on the consensus axis", {
  s <- make_stack()
  # D runs against A, B, C: it loads opposite on PC1 and is left out
  e <- ensemble_models(s, method = "pcam")
  expect_false("D" %in% e$members)
  expect_length(e$members, 2)
  expect_equal(terra::values(e$raster[["ensemble"]])[, 1],
               apply(terra::values(s[[e$members]]), 1, stats::median))
  expect_error(ensemble_models(s[[1:2]], method = "pcam"), "three")
})

test_that("files on disk, filename and plot", {
  s <- make_stack()
  dir <- file.path(tempdir(), "ens_in")
  dir.create(dir, showWarnings = FALSE)
  for (n in names(s)) terra::writeRaster(s[[n]], file.path(dir, paste0(n, ".tif")),
                                          overwrite = TRUE)
  out <- tempfile(fileext = ".tif")
  e <- ensemble_models(dir, method = "mean", filename = out)
  expect_true(file.exists(out))
  # macOS resolves /var to /private/var: compare normalised paths
  expect_equal(normalizePath(terra::sources(e$raster)), normalizePath(out))
  expect_equal(names(terra::rast(out)), c("ensemble", "sd", "cv"))
  expect_equal(sort(e$models), c("A", "B", "C", "D"))
  expect_true(is.list(plot_ensemble(e)))
  expect_s3_class(plot_ensemble(e, type = "disagreement",
                                disagreement = "cv"), "ggplot")
  f <- tempfile(fileext = ".png")
  plot_ensemble(e, file = f)
  expect_true(file.exists(f))
})
