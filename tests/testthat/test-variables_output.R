# 13 bio-like layers in five clusters, three periods, six tile regions
selection_fixture <- function() {
  sim <- simulate_enm_data(seed = 1)
  base <- sim$env[[1]]
  fld <- function(seed, km = 150) {
    (random_variables(base, n = 1, range_km = km, seed = seed) - 0.5) * 3
  }
  lat <- (terra::init(base, "y") - 52) / 3
  period <- function(shift, sd) {
    Tm <- -lat * 2 + fld(1) + shift
    S <- fld(2) + lat * sd
    P <- fld(3) - 0.3 * Tm
    D <- fld(4)
    Sp <- fld(5)
    e <- function(i, a = 0.25) a * fld(100 + i + round(shift * 10), 60)
    r <- c(Tm + e(1), S + e(2), Tm + S + e(3), Tm - S + e(4), 2 * S + e(5, .4),
           Tm + 0.8 * S + e(6), Tm - 0.9 * S + e(7), P + e(8),
           P + 0.5 * Sp + e(9), P - 0.6 * D + e(10), Sp + e(11, .5),
           P + 0.45 * Sp + e(12), P - 0.55 * D + e(13))
    names(r) <- sprintf("bio%02d", c(1, 4, 5, 6, 7, 10:17))
    r
  }
  list(present = period(0, 0.3), cold = period(-3, 1.2),
       cool = period(-1.5, -0.4))
}
fx <- selection_fixture()
tiles <- list(x = c(14, 17.5, 21, 24), y = c(49, 52, 55))
es_fx <- sample_environment(fx, n = 4000, regions = tiles, seed = 1,
                            verbose = FALSE)
vs_fx <- suppressMessages(select_variables(
  es_fx, sizes = 2:6, min_region_n = 400,
  require_one_of = list(
    temperature = c("bio01", "bio05", "bio06", "bio10", "bio11"),
    moisture = c("bio12", "bio13", "bio14", "bio16", "bio17"))))

test_that("best_sets explains why a size has no sets", {
  expect_true(all(2:4 %in% vs_fx$sets$n_vars))
  expect_false(any(5:6 %in% vs_fx$sets$n_vars))
  expect_message(b <- best_sets(vs_fx, n = 10, size = 6),
                 "No set of 6 variables passes")
  expect_message(best_sets(vs_fx, size = 6), "Sizes with passing sets: 2, 3, 4")
  expect_message(best_sets(vs_fx, size = 6), "clusters of collinear")
  expect_s3_class(b, "data.frame")
  expect_equal(nrow(b), 0)
  expect_true("set" %in% names(b))
  # a size that passes says nothing; a mix returns the passing size
  expect_silent(b4 <- best_sets(vs_fx, size = 4))
  expect_true(all(b4$n_vars == 4))
  expect_message(b46 <- best_sets(vs_fx, size = c(4, 6)), "No set of 6")
  expect_equal(b46, b4)
  expect_message(best_sets(vs_fx, size = 9), "not searched")
  # nothing passes at all
  none <- suppressWarnings(suppressMessages(
    select_variables(es_fx, sizes = 6, min_region_n = 400)))
  expect_message(expect_null(best_sets(none)), "No set passes")
  expect_match(recommend_sets(none, print = FALSE)[1], "No set passes")
})

test_that("the printed selection leads with the best sets and a summary", {
  out <- capture.output(print(vs_fx))
  i <- function(p) grep(p, out)[1]
  expect_lt(i("BEST SETS"), i("IN SHORT"))
  expect_lt(i("IN SHORT"), i("DETAILS"))
  expect_lt(i("DETAILS"), i("Sets surviving each filter"))
  expect_lt(i("Sets surviving each filter"), i("Regions, all variables"))
  expect_true(any(grepl("no compatible combination of 5 and 6", out)))
  expect_false(any(grepl("IN SHORT", capture.output(print(vs_fx,
                                                          advice = FALSE)))))
})

test_that("recommend_sets describes the selection", {
  txt <- recommend_sets(vs_fx, print = FALSE)
  best4 <- vs_fx$sets$set[vs_fx$sets$n_vars == 4][1]
  expect_match(txt[1], paste0("Suggested set: ", best4), fixed = TRUE)
  expect_match(txt[1], "4 variables, the largest size that passes")
  expect_match(txt[1], vs_fx$sets$worst_region[vs_fx$sets$set == best4],
               fixed = TRUE)
  expect_true(any(grepl("Close alternatives", txt)))
  expect_true(any(grepl("Best of the smaller sizes", txt)))
  expect_true(any(grepl("No set of 5 and 6 variables passes", txt)))
  expect_true(any(grepl("Most often in the Pareto sets", txt)))
  expect_true(any(grepl("excluded only by the regional screen", txt)))
  expect_match(txt[length(txt)], "ecological grounds")
  # the pair quoted as changing is one of those flagged
  ch <- grep("collinear in some data sets", txt, value = TRUE)
  expect_length(ch, 1)
  p <- vs_fx$pairs[vs_fx$pairs$collinear_only_in_some, ]
  expect_true(any(vapply(seq_len(nrow(p)), function(k) {
    grepl(paste(p$var1[k], "~", p$var2[k]), ch, fixed = TRUE)
  }, logical(1))))
  expect_match(ch, paste0("^", nrow(p), " pair"))
  expect_output(res <- recommend_sets(vs_fx), "^- Suggested set")
  expect_equal(res, txt)
  # regions too small are named
  small <- suppressMessages(select_variables(es_fx, sizes = 2,
                                             min_region_n = 700))
  expect_true(any(grepl("too few sampled cells", recommend_sets(
    small, print = FALSE))))
})

test_that("skip_regions leaves regions out of the screen", {
  full <- vs_fx$regions
  expect_true(all(full$screened))
  sk <- suppressMessages(select_variables(
    es_fx, sizes = 2:3, min_region_n = 400,
    skip_regions = c("cold | x3_y2", "x1_y1")))
  rg <- sk$regions
  off <- rg[!rg$screened, ]
  expect_setequal(paste(off$data_set, off$region, sep = " | "),
                  c("cold | x3_y2", "present | x1_y1", "cold | x1_y1",
                    "cool | x1_y1"))
  expect_true(all(off$why_not == "skipped"))
  expect_false(any(c("cold | x3_y2", "cool | x1_y1") %in%
                     names(sk$regional_rho)))
  expect_equal(length(sk$regional_rho), 18 - 4)
  # fewer regional limits, at least as many compatible pairs
  expect_gte(sum(sk$pairs$compatible), sum(vs_fx$pairs$compatible))
  expect_output(print(sk), "Not screened: .*skipped")
})

test_that("masks leave cells out of single data sets", {
  ice <- terra::as.polygons(terra::ext(14, 24, 52, 55), crs = "EPSG:4326")
  es <- sample_environment(fx, n = 2000, regions = tiles, seed = 1,
                           verbose = FALSE, exclude = list(cold = ice))
  expect_true(all(es$sets$cold$.y < 52))
  expect_true(any(es$sets$present$.y > 52))
  expect_true(any(es$sets$cool$.y > 52))
  # the common cells are unmasked in every data set
  expect_true(all(es$common$present$.y < 52))
  expect_equal(unname(es$masked[["cold"]]), 30 / 60, tolerance = 0.02)
  expect_true(is.na(es$masked[["present"]]))
  expect_output(print(es), "Masked")
  # the northern tiles have no cells in the masked period
  expect_false(any(grepl("y2", es$sets$cold$.region)))
  vs <- suppressMessages(select_variables(es, sizes = 2, min_region_n = 300))
  expect_false(any(vs$regions$data_set == "cold" &
                     grepl("y2", vs$regions$region)))
  # a raster mask on another grid: cells with a value other than NA and 0
  m <- terra::rast(nrows = 12, ncols = 20, xmin = 14, xmax = 24, ymin = 49,
                   ymax = 55, crs = "EPSG:4326")
  terra::values(m) <- 0
  m[terra::yFromCell(m, 1:terra::ncell(m)) > 52] <- 500
  es2 <- sample_environment(fx, n = 2000, seed = 1, verbose = FALSE,
                            exclude = list(cold = m))
  expect_true(all(es2$sets$cold$.y < 52))
  m[m == 0] <- NA
  es3 <- sample_environment(fx, n = 2000, seed = 1, verbose = FALSE,
                            exclude = list(cold = m))
  expect_equal(es3$sets$cold, es2$sets$cold)
  # keep is the complement; one mask applies to all data sets
  es4 <- sample_environment(fx, n = 2000, seed = 1, verbose = FALSE,
                            keep = ice)
  expect_true(all(vapply(es4$sets, function(d) all(d$.y > 52), logical(1))))
  expect_true(all(!is.na(es4$masked)))
  # a mask in another CRS, as sf
  ice_p <- sf::st_as_sf(terra::project(ice, "EPSG:3035"))
  es5 <- sample_environment(fx, n = 2000, seed = 1, verbose = FALSE,
                            exclude = list(cold = ice_p))
  expect_lt(max(es5$sets$cold$.y), 52.2)
  expect_error(sample_environment(fx, n = 100, exclude = list(zz = ice),
                                  verbose = FALSE), "named by data set")
  expect_error(sample_environment(fx, n = 100, exclude = list(cold = 1),
                                  verbose = FALSE), "polygons or a raster")
  expect_error(sample_environment(
    fx, n = 100, verbose = FALSE,
    exclude = terra::as.polygons(terra::ext(0, 40, 40, 60),
                                 crs = "EPSG:4326")), "after masking")
})

test_that("past periods without a mask give a warning", {
  old <- options(preENM.past_warning = TRUE)
  on.exit(options(old))
  d <- stats::setNames(fx, c("present", "LGM", "YD"))
  expect_warning(sample_environment(d, n = 300, regions = tiles, seed = 1,
                                    verbose = FALSE),
                 "LGM, YD look like past periods.*glaciated.*skip_regions")
  w <- tryCatch(sample_environment(d, n = 300, seed = 1, verbose = FALSE),
                warning = function(w) conditionMessage(w))
  expect_false(grepl("skip_regions", w))
  # a mask for the period, or neutral names, and no warning
  ice <- terra::as.polygons(terra::ext(14, 24, 54, 55), crs = "EPSG:4326")
  expect_warning(sample_environment(d, n = 300, seed = 1, verbose = FALSE,
                                    exclude = list(LGM = ice)),
                 "Data set\\(s\\) YD look")
  expect_no_warning(sample_environment(d, n = 300, seed = 1, verbose = FALSE,
                                       exclude = list(LGM = ice, YD = ice)))
  expect_no_warning(sample_environment(fx, n = 300, seed = 1,
                                       verbose = FALSE))
  expect_equal(looks_past(c("present", "ssp585_2070", "past", "LGM", "21ka",
                            "t_-190", "Younger_Dryas", "mid-Holocene", "YD"),
                          as.list(1:9)),
               c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE))
  # recognised from file names too
  expect_true(looks_past("a", list(c("x/CHELSA_TraCE21k_bio01_-190_V1.0.tif"))))
  options(preENM.past_warning = FALSE)
  expect_no_warning(sample_environment(d, n = 300, seed = 1, verbose = FALSE))
})

test_that("misplaced arguments are reported", {
  expect_warning(sample_environment(fx, n = 300, seed = 1, verbose = FALSE,
                                    bounds = list(fix_scale = TRUE)),
                 "not variables.*fix_scale = ... outside")
  expect_warning(sample_environment(fx, n = 300, seed = 1, verbose = FALSE,
                                    bounds = list(scd = c(0, 365))),
                 "not variables and are ignored: scd")
  expect_error(sample_environment(fx, n = 300, seed = 1, verbose = FALSE,
                                  bounds = c(0, 1)), "named list")
  expect_message(sample_environment(fx["present"], n = 300, seed = 1,
                                    ratio_vars = "bio12"),
                 "only applies to anomaly")
})

test_that("the figures of the selection can be drawn", {
  p <- plot_variable_sets(vs_fx)
  expect_s3_class(p, "ggplot")
  expect_true(all(table(p$data$row) > 0))
  # at most n rows per size, values as in the table
  expect_true(all(tapply(p$data$set, p$data$size,
                         function(z) length(unique(z))) <= 8))
  one <- p$data[p$data$set == vs_fx$sets$set[1] &
                  grepl("^worst_maxr", p$data$metric), ]
  expect_equal(one$value, vs_fx$sets$worst_maxr[1])
  expect_true(all(p$data$good >= 0 & p$data$good <= 1))
  expect_equal(nlevels(plot_variable_sets(vs_fx, n = 3)$data$row), 9)
  sc <- plot_variable_sets(vs_fx, type = "scatter")
  expect_s3_class(sc, "ggplot")
  expect_match(sc$labels$caption, "4 variables: 1 = ")
  # naming an axis gives the scatter
  expect_equal(plot_variable_sets(vs_fx, "novelty")$labels$x,
               "novelty (lower is better)")
  expect_error(plot_variable_sets(vs_fx, "nope"), "Unknown score")
  expect_silent(print(sc))

  t0 <- plot_variable_tree(vs_fx)
  expect_s3_class(t0, "ggplot")
  expect_match(t0$labels$subtitle, "best-ranked sets of each size")
  expect_silent(print(t0))
  expect_s3_class(plot_variable_tree(vs_fx, c("bio05", "bio15", "bio17")),
                  "ggplot")
  t2 <- plot_variable_tree(vs_fx, c(best = vs_fx$sets$set[1],
                                    mine = "bio05 bio15 bio17"))
  expect_match(t2$labels$subtitle, "the given sets")
  expect_silent(print(t2))
  expect_error(plot_variable_tree(vs_fx, "bio05 bio99"), "Unknown variable")
  d <- tempfile()
  f <- write_variable_selection(vs_fx, d, "sp")
  expect_true(all(file.exists(f)))
  expect_true(all(c("sp_best_sets.csv", "sp_summary.txt", "sp_sets.png",
                    "sp_sets_scatter.png", "sp_tree.png") %in%
                    list.files(d)))
  expect_equal(readLines(file.path(d, "sp_summary.txt")),
               recommend_sets(vs_fx, print = FALSE))
  expect_equal(utils::read.csv(file.path(d, "sp_best_sets.csv"))$set,
               best_sets(vs_fx, n = 10)$set)
})
