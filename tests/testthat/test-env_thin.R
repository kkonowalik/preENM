clim <- make_climate()$present[[c("bio01", "bio04", "bio12", "bio15")]]
clim$bio12 <- (clim$bio12 - terra::global(clim$bio12, "min")[[1]]) * 1000

# clustered records: many in one corner, few spread over the rest
make_pts <- function(seed = 1) {
  set.seed(seed)
  rbind(data.frame(lon = runif(250, 5, 12), lat = runif(250, 44, 48)),
        data.frame(lon = runif(60, -9, 34), lat = runif(60, 41, 69)))
}

test_that("thin_environment evens density without truncating ranges", {
  pts <- make_pts()
  th <- thin_environment(pts, clim, radius_km = 30, block_km = 200,
                         seed = 1, verbose = FALSE)
  expect_s3_class(th, "env_thinning")
  k <- kept(th)
  expect_lt(nrow(k), nrow(pts))
  expect_gt(nrow(k), 10)
  expect_equal(nrow(k) + nrow(removed(th)), nrow(th$data))
  rg <- th$report$ranges
  expect_equal(rg$min_all, rg$min_kept)
  expect_equal(rg$max_all, rg$max_kept)
  # the dense corner loses a larger share than the rest
  dense <- th$data$lon >= 5 & th$data$lon <= 12 & th$data$lat <= 48
  expect_lt(mean(th$data$kept[dense]), mean(th$data$kept[!dense]))
  # every removed record has a kept one within a bin on every axis
  expect_lte(unname(th$report$removed_distance["100%"]),
             sqrt(th$settings$n_axes) + 1e-9)
  expect_true(all(th$data$keep_reason[th$data$kept] %in%
                    c("protected", "extreme_minmax", "extreme_tail", "hull",
                      "env_bin", "block_quota")))
  expect_true("bio12" %in% th$settings$log_vars)
  expect_output(print(th), "kept")
})

test_that("one record per bin and the block quota", {
  pts <- make_pts(2)
  th <- thin_environment(pts, clim, radius_km = 30, block_min = 0,
                         protect_minmax = FALSE, protect_hull = FALSE,
                         verbose = FALSE)
  expect_equal(sum(th$data$kept), length(unique(th$data$bin)))
  expect_equal(max(table(th$data$bin[th$data$kept])), 1)
  thq <- thin_environment(pts, clim, radius_km = 30, block_km = 300,
                          block_min = 3, verbose = FALSE)
  expect_equal(thq$report$blocks_below_quota, 0)
  b <- thq$report$blocks
  expect_true(all(b$n_kept >= pmin(b$n_all, 3)))
})

test_that("protected records, raw axes, templates and errors", {
  pts <- make_pts(3)
  pts$source <- ifelse(seq_len(nrow(pts)) %in% 1:20, "sNMF", "GBIF")
  th <- thin_environment(pts, clim, radius_km = 30,
                         protect = list(source = "sNMF"), verbose = FALSE)
  expect_true(all(th$data$kept[th$data$source == "sNMF"]))
  expect_true(all(th$data$keep_reason[th$data$source == "sNMF"] ==
                    "protected"))
  raw <- thin_environment(pts, clim, axes = "raw", radius_km = 30,
                          verbose = FALSE)
  expect_equal(raw$settings$n_axes, 4)
  sub <- thin_environment(pts[1:100, ], clim, template = th, verbose = FALSE)
  expect_equal(sub$bin_widths, th$bin_widths)
  fixed <- thin_environment(pts, clim, radius_km = NA, window = 1,
                            verbose = FALSE)
  expect_s3_class(fixed, "env_thinning")
  # a 1 km radius on ~50 km cells covers only the record's own cell
  expect_error(thin_environment(pts, clim, radius_km = 1, verbose = FALSE),
               "bin width is 0")
  r1 <- thin_environment(pts, clim, radius_km = 30, pick = "random",
                         seed = 5, verbose = FALSE)
  r2 <- thin_environment(pts, clim, radius_km = 30, pick = "random",
                         seed = 5, verbose = FALSE)
  expect_equal(r1$data$kept, r2$data$kept)
  expect_error(thin_environment(pts[1, ], clim, verbose = FALSE), "Fewer")
  expect_error(thin_environment(pts, clim, max_per_bin = 0, verbose = FALSE),
               "max_per_bin")
  expect_error(kept(list()), "must come from")
})

test_that("plot and write the thinning", {
  th <- thin_environment(make_pts(), clim, radius_km = 30, verbose = FALSE)
  expect_s3_class(plot_env_thinning(th), "ggplot")
  d <- tempfile("envthin")
  f <- write_env_thinning(th, d, "dryas")
  expect_true(all(file.exists(f)))
  out <- utils::read.csv(file.path(d, "dryas.csv"))
  expect_equal(nrow(out), sum(th$data$kept))
  tpl <- readRDS(file.path(d, "dryas_pca.rds"))
  again <- thin_environment(make_pts(), clim, template = tpl, verbose = FALSE)
  expect_equal(again$data$kept, th$data$kept)
})
