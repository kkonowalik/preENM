# Partitioning projection uncertainty

grid_r <- function() {
  terra::rast(nrows = 10, ncols = 12, xmin = 10, xmax = 13, ymin = 50,
              ymax = 52, crs = "EPSG:4326")
}

test_that("projections differing along one factor give it all variance", {
  r <- grid_r()
  set.seed(1)
  base <- terra::rast(r, vals = runif(terra::ncell(r), 0.2, 0.6))
  design <- expand.grid(gcm = c("g1", "g2", "g3"), scenario = c("s1", "s2"),
                        algorithm = c("a", "b"), stringsAsFactors = FALSE)
  # only the GCM changes the prediction
  shift <- c(g1 = 0, g2 = 0.1, g3 = 0.25)
  proj <- terra::rast(lapply(seq_len(nrow(design)), function(i) {
    base + shift[[design$gcm[i]]]
  }))
  pu <- partition_uncertainty(proj, design)
  v <- terra::values(pu$raster)
  expect_true(all(abs(v[, "share_gcm"] - 1) < 1e-6))
  expect_true(all(abs(v[, "share_scenario"]) < 1e-6))
  a <- pu$shares[pu$shares$scope == "area", ]
  expect_equal(a$share[a$source == "gcm"], 1)
  expect_equal(a$dominant_share[a$source == "gcm"], 1)
  expect_output(print(pu), "the choice of gcm explains")
  # total variance by hand: var(c(0, 0.1, 0.25) each 4 times)
  expect_equal(unname(v[1, "variance"]), stats::var(rep(shift, each = 4)),
               tolerance = 1e-6)
})

test_that("shares sum to 1 per cell and match an ANOVA by hand", {
  r <- grid_r()
  set.seed(2)
  design <- expand.grid(gcm = c("g1", "g2"), algorithm = c("a", "b", "c"),
                        run = 1:2, stringsAsFactors = FALSE)
  proj <- terra::rast(lapply(seq_len(nrow(design)), function(i) {
    terra::rast(r, vals = runif(terra::ncell(r)))
  }))
  pu <- partition_uncertainty(proj, design)
  v <- terra::values(pu$raster)
  sh <- v[, grep("^share_", colnames(v))]
  expect_true(all(abs(rowSums(sh) - 1) < 1e-6))   # stored as float32
  # cell 1 by hand with aov()
  y <- as.vector(terra::values(proj)[1, ])
  d <- data.frame(y = y, lapply(design, factor))
  ss <- summary(stats::aov(y ~ gcm + algorithm + run, data = d))[[1]][["Sum Sq"]]
  expect_equal(unname(v[1, c("share_gcm", "share_algorithm", "share_run",
                             "share_residual")]), ss / sum(ss),
               tolerance = 1e-6)
  expect_equal(terra::cats(pu$raster[["dominant"]])[[1]]$source,
               c("gcm", "algorithm", "run", "interactions/residual"))
})

test_that("region, scales, unbalanced designs and plot", {
  r <- grid_r()
  set.seed(3)
  design <- expand.grid(gcm = c("g1", "g2"), algorithm = c("a", "b"),
                        stringsAsFactors = FALSE)
  # west half: GCM matters; east half: algorithm matters
  xy <- terra::xyFromCell(r, seq_len(terra::ncell(r)))
  west <- xy[, 1] < 11.5
  proj <- terra::rast(lapply(seq_len(nrow(design)), function(i) {
    g <- design$gcm[i] == "g2"
    a <- design$algorithm[i] == "b"
    terra::rast(r, vals = 0.3 + ifelse(west, 0.3 * g + 0.01 * a,
                                       0.01 * g + 0.3 * a))
  }))
  reg <- terra::vect("POLYGON ((10 50, 11.5 50, 11.5 52, 10 52, 10 50))",
                     crs = "EPSG:4326")
  pu <- partition_uncertainty(proj, design, region = reg)
  rg <- pu$shares[pu$shares$scope == "region", ]
  expect_gt(rg$share[rg$source == "gcm"], 0.99)
  ar <- pu$shares[pu$shares$scope == "area", ]
  expect_equal(ar$share[ar$source == "gcm"], 0.5, tolerance = 0.02)
  # a layer written on 0-100 is rescaled
  p2 <- proj
  p2[[1]] <- p2[[1]] * 100
  expect_message(pu2 <- partition_uncertainty(p2, design), "divided by 100")
  expect_equal(terra::values(pu2$raster), terra::values(pu$raster))
  # unbalanced
  expect_warning(partition_uncertainty(proj[[c(1:4, 1)]],
                                       design[c(1:4, 1), ]),
                 "not balanced")
  expect_error(partition_uncertainty(proj, design[1:3, ]), "one row per layer")
  expect_true(is.list(plot_partition(pu)))
  f <- tempfile(fileext = ".png")
  plot_partition(pu, file = f)
  expect_true(file.exists(f))
})
