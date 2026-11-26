# a U-shaped range: two arms and a base, with a hollow between the arms
make_u <- function(seed = 2) {
  set.seed(seed)
  rbind(data.frame(lon = runif(150, 5, 8), lat = runif(150, 45, 55)),
        data.frame(lon = runif(150, 17, 20), lat = runif(150, 45, 55)),
        data.frame(lon = runif(100, 5, 20), lat = runif(100, 44, 46.5)))
}
inside_ll <- function(m, lon, lat) {
  p <- sf::st_transform(sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326),
                        m$settings$crs)
  lengths(sf::st_intersects(p, m$m)) > 0
}
hollow <- c(12.5, 51)

test_that("spanning tree equals single-linkage clustering", {
  set.seed(1)
  x <- runif(200, 0, 5e5)
  y <- runif(200, 0, 5e5)
  ms <- mst_edges(x, y)
  hc <- stats::hclust(stats::dist(cbind(x, y) / 1000), "single")
  expect_equal(sort(ms$mst$km), sort(hc$height))
  # collinear points and tiny inputs
  expect_equal(mst_edges(c(0, 1000, 3000), c(0, 0, 0))$mst$km, c(1, 2))
  expect_equal(nrow(mst_edges(0, 0)$mst), 0)
  # clusters at a linking distance
  cl <- link_components(200, ms$mst, 1e9)
  expect_equal(max(cl), 1)
})

test_that("automatic L is the last link before the largest jump", {
  h <- c(10, 20, 110, 120, 130, 400, 420)
  expect_equal(auto_link(h, 100, 600), 130)
  expect_true(is.na(auto_link(c(10, 20), 100, 600)))
  set.seed(3)
  main <- data.frame(lon = runif(400, 5, 15), lat = runif(400, 45, 50))
  outl <- data.frame(lon = c(16.5, 18, 19.6, 21.4, 23.5),
                     lat = c(50, 50, 50, 50, 50))           # growing gaps
  island <- data.frame(lon = runif(20, 30, 31), lat = runif(20, 55, 56))
  lk <- linking_distance(rbind(main, outl, island), link_search = c(50, 600))
  expect_s3_class(lk, "m_linking")
  expect_true(lk$auto)
  expect_equal(lk$L, auto_link(lk$links$km, 50, 600))
  expect_true(lk$L >= 50 && lk$L <= 600)
  expect_equal(lk$n_clusters,
               lk$table$n_clusters[lk$table$link_km == round(lk$L, 1)])
  expect_gte(lk$L2, lk$L)
  expect_output(print(lk), "Linking distance")
  expect_warning(linking_distance(main, link_search = c(2000, 3000)),
                 "Widen")
})

test_that("default Delaunay M leaves the hollow open; hulls and experts fill it", {
  u <- make_u()
  m <- delimit_m(u, link_km = 150)
  expect_s3_class(m, "m_area")
  expect_equal(m$settings$method, "delaunay")
  expect_equal(m$settings$closing_km, 75)               # L / 2
  expect_false(inside_ll(m, hollow[1], hollow[2]))
  expect_equal(m$report$records_inside_pct[3], 100)
  h <- delimit_m(u, link_km = 150, method = "cluster_hull")
  expect_true(inside_ll(h, hollow[1], hollow[2]))
  hc <- delimit_m(u, link_km = 150, method = "cluster_hull", hull_ratio = 0.3)
  expect_lt(geom_area(hc$m), geom_area(h$m))
  ch <- delimit_m(u, link_km = 150, method = "convex_hull")
  expect_true(inside_ll(ch, hollow[1], hollow[2]))
  b <- delimit_m(u, link_km = 150, method = "buffer", buffer_km = 30,
                 fill_holes = FALSE)
  expect_lt(geom_area(b$m), geom_area(m$m))
  # expert polygon over the hollow
  ex <- sf::st_sfc(sf::st_polygon(list(rbind(c(8, 49), c(17, 49), c(17, 53),
                                             c(8, 53), c(8, 49)))),
                   crs = 4326)
  me <- delimit_m(u, link_km = 150, add = ex,
                  add_sites = data.frame(lon = 26, lat = 58, buffer_km = 50))
  expect_true(inside_ll(me, hollow[1], hollow[2]))
  expect_true(inside_ll(me, 26, 58))
  expect_equal(me$settings$expert, c("add", "sites"))
  rep <- m_report(me)
  expect_gt(rep$area_km2[rep$part == "added_by_expert"], 1e5)
  expect_equal(rep$area_km2[3], rep$area_km2[1] + rep$area_km2[2],
               tolerance = 0.01)
  expect_identical(m$automatic, me$automatic)            # automatic unchanged
})

test_that("corridors join isolated parts within L2", {
  u <- make_u()
  pts <- rbind(cbind(u, time = "modern"),
               data.frame(lon = c(12, 13), lat = c(52, 53), time = "fossil"))
  m <- delimit_m(pts, link_km = 150, group = "time", corridors = TRUE)
  expect_gt(m$settings$L2, 150)
  expect_equal(length(m$corridors), 1)
  expect_true(inside_ll(m, hollow[1], hollow[2]))
  expect_true(all(c("records_inside_pct_modern", "records_inside_pct_fossil")
                  %in% names(m$report)))
})

test_that("holes are filled unless switched off", {
  a <- seq(0, 2 * pi, length.out = 200)[-200]
  ring <- data.frame(lon = 15 + 5 * cos(a), lat = 50 + 3 * sin(a))
  f <- delimit_m(ring, link_km = 100, method = "buffer", buffer_km = 40,
                 closing_km = 0)
  o <- delimit_m(ring, link_km = 100, method = "buffer", buffer_km = 40,
                 closing_km = 0, fill_holes = FALSE)
  expect_true(inside_ll(f, 15, 50))
  expect_false(inside_ll(o, 15, 50))
})

test_that("extent, exclude, land and thinning", {
  u <- make_u()
  m <- delimit_m(u, link_km = 150, extent = c(4, 21, 43, 50))
  bb <- sf::st_bbox(sf::st_transform(m$m, 4326))
  expect_lte(bb[["ymax"]], 50.01)
  expect_gt(m$settings$n_dropped, 0)
  exc <- sf::st_sfc(sf::st_polygon(list(rbind(c(16, 44), c(21, 44), c(21, 56),
                                              c(16, 56), c(16, 44)))),
                    crs = 4326)
  me <- delimit_m(u, link_km = 150, exclude = exc)
  expect_false(inside_ll(me, 18.5, 50))
  expect_true(all(me$all_points$lon < 16))
  land <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 40), c(12, 40), c(12, 60),
                                               c(0, 60), c(0, 40)))),
                     crs = 4326)
  ml <- delimit_m(u, link_km = 150, land = land)
  expect_false(inside_ll(ml, 18.5, 50))
  # thinning per group and source priority
  p <- data.frame(lon = c(10, 10.01, 10.02), lat = 50,
                  source = c("GBIF", "own", "GBIF"),
                  time = c("modern", "modern", "fossil"))
  pr <- m_prepare(p, NULL, NULL, 10, "time", c("own", "GBIF"), "source",
                  "auto", NULL, "EPSG:4326")
  expect_equal(sort(pr$thin$source), c("GBIF", "own"))
  expect_true("own" %in% pr$thin$source[pr$thin$time == "modern"])
})

test_that("projection choice and input checks", {
  expect_match(equal_area_crs(c(-170, 170), c(65, 70)), "lat_0=90")
  expect_match(equal_area_crs(c(10, 20), c(40, 50)), "lat_0=45 \\+lon_0=15")
  u <- make_u()
  expect_error(delimit_m(u, method = "alpha"), "should be one of")
  expect_error(delimit_m(u, link_km = 150,
                         add_sites = data.frame(lon = 1, lat = 2)),
               "buffer_km")
  expect_error(delimit_m(u, link_km = 150, priority = "a"), "not found")
  expect_error(delimit_m(u, extent = c(100, 110, 0, 1)), "No records")
  one <- delimit_m(data.frame(lon = 10, lat = 50), link_km = 100)
  expect_gt(geom_area(one$m), 0)
})

test_that("scenarios are compared with the reference", {
  u <- make_u()
  expect_warning(
    sc <- m_scenarios(u, link_km = 150,
                      scenarios = list(delaunay = list(),
                                       hull = list(method = "cluster_hull"),
                                       bad = list(method = "nope")),
                      reference = "delaunay"),
    "failed")
  expect_s3_class(sc, "m_scenarios")
  expect_equal(names(sc$areas), c("delaunay", "hull"))
  cmp <- sc$comparison
  expect_equal(cmp$jaccard_with_reference[cmp$scenario == "delaunay"], 1)
  expect_lt(cmp$jaccard_with_reference[cmp$scenario == "hull"], 1)
  expect_equal(cmp$share_in_reference_pct[cmp$scenario == "delaunay"], 100)
  expect_output(print(sc), "reference: delaunay")
  expect_error(m_scenarios(u, link_km = 150, reference = "zzz"), "reference")
  full <- m_scenarios(u, link_km = 150)
  expect_equal(nrow(full$comparison), 7)
})

test_that("reports, class shares, clipping, plots and files", {
  sim <- simulate_enm_data(seed = 1, n_presence = 150)
  m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
  expect_output(print(m), "Accessible area")
  rep <- m_report(m, raster = sim$env)
  expect_true("land_area_km2" %in% names(rep))
  expect_gt(rep$land_area_km2[3], 0)
  expect_lte(rep$land_area_km2[3], rep$area_km2[3] * 1.01)
  cls <- terra::classify(sim$env[["env1"]],
                         rbind(c(-Inf, 0.5, 1), c(0.5, Inf, 2)))
  sh <- m_class_shares(m, cls, labels = c("west", "east"))
  expect_equal(sum(sh$pct_M), 100, tolerance = 0.05)
  expect_equal(sh$label, c("west", "east"))
  expect_equal(length(koppen_geiger_labels()), 31)
  cp <- class_polygons(cls, "east", labels = c("west", "east"))
  expect_s3_class(cp, "sf")
  em <- crop_to_m(sim$env, m)
  expect_lt(sum(!is.na(terra::values(em[[1]]))),
            sum(!is.na(terra::values(sim$env[[1]]))) + 1)
  expect_s3_class(plot_m_area(m), "ggplot")
  expect_s3_class(plot_linking(m), "ggplot")
  expect_s3_class(plot_linking(m, "clusters"), "ggplot")
  expect_error(plot_linking(1), "linking_distance")
  td <- tempfile()
  f <- write_m_area(m, td, "test", raster = sim$env)
  expect_true(all(file.exists(f)))
  expect_true(all(c("M", "M_automatic", "records") %in%
                    sf::st_layers(f[["gpkg"]])$name))
  back <- crop_to_m(sim$env, f[["gpkg"]])      # clip with the written file
  expect_s4_class(back, "SpatRaster")
  sc <- m_scenarios(sim$presence, link_km = 120, thin_km = 5,
                    scenarios = list(a = list(), b = list(link_mult = 2)))
  expect_s3_class(plot_m_area(sc), "ggplot")
  f2 <- write_m_area(sc, td, "sc")
  expect_true(file.exists(f2[["comparison"]]))
  expect_error(m_report(1), "delimit_m")
})

test_that("adaptive M: clusters reach by their spread, singletons stay small", {
  set.seed(1)
  # a dense cluster of 60 records (about 1 x 1 degree) and one lone record
  # far away
  cl <- data.frame(lon = 20 + stats::runif(60), lat = 50 + stats::runif(60))
  pts <- rbind(cl, data.frame(lon = 26, lat = 54))
  m <- suppressMessages(delimit_m(pts, method = "adaptive", link_km = 150,
                                  thin_km = 1, close_km = 0,
                                  closing_km = 0, adaptive_min_km = 5))
  parts <- sf::st_cast(sf::st_geometry(m$m), "POLYGON")
  a <- sort(as.numeric(sf::st_area(parts)) / 1e6)
  expect_length(a, 2)
  # the lone record keeps a circle of 5 km
  expect_equal(a[1], pi * 5^2, tolerance = 0.05)
  # the cluster reaches about its mean pairwise distance (~55 km) beyond
  # its records
  expect_gt(a[2], 2 * 1e4)
  mw <- suppressMessages(delimit_m(pts, method = "adaptive", link_km = 150,
                                   thin_km = 1, close_km = 0,
                                   closing_km = 0, adaptive_mult = 0.5))
  expect_lt(as.numeric(sf::st_area(mw$m)), as.numeric(sf::st_area(m$m)))
})

sim_landscape <- function() {
  r <- terra::rast(nrows = 60, ncols = 90, xmin = 0, xmax = 9, ymin = 44,
                   ymax = 50, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  set.seed(1)
  a <- terra::rast(r, vals = 1 - abs(xy[, 2] - 47) / 3 +
                     stats::rnorm(nrow(xy), 0, 0.02))
  b <- terra::rast(r, vals = stats::rnorm(nrow(xy), 0, 0.05))
  env <- c(a, b)
  names(env) <- c("a", "b")
  env
}

test_that("m_simulate(): M holds all records and respects barriers", {
  env <- sim_landscape()
  set.seed(2)
  pts <- data.frame(lon = stats::runif(30, 1, 3), lat = stats::runif(30, 46,
                                                                      48))
  # a 1-degree wide barrier (10 cells) between x = 4 and 5, far wider than
  # the kernel (SD 1 cell)
  bar <- sf::st_sfc(sf::st_polygon(list(rbind(c(4, 43), c(5, 43), c(5, 51),
                                              c(4, 51), c(4, 43)))),
                    crs = 4326)
  ms <- m_simulate(pts, env, barrier = bar, n_rep = 4, n_events = 40,
                   kernel_sd = 1, seed = 1)
  expect_s3_class(ms, "m_simulation")
  inside <- lengths(sf::st_intersects(
    sf::st_as_sf(pts, coords = c("lon", "lat"), crs = 4326),
    sf::st_transform(ms$m, 4326))) > 0
  expect_true(all(inside))
  bb <- sf::st_bbox(sf::st_transform(ms$m, 4326))
  expect_lte(bb[["xmax"]], 4 + 1e-6)          # nothing beyond the barrier
  # without the barrier dispersal crosses
  m2 <- m_simulate(pts, env, n_rep = 4, n_events = 40, kernel_sd = 1,
                   seed = 1)
  expect_gt(sf::st_bbox(sf::st_transform(m2$m, 4326))[["xmax"]], 5)
  expect_output(print(ms), "IN SHORT")
  expect_s3_class(plot(ms), "ggplot")
  # as a method of delimit_m() and in m_scenarios()
  dm <- suppressMessages(delimit_m(pts, method = "simulation", env = env,
                                   link_km = 100, thin_km = 1,
                                   simulation = list(n_rep = 3,
                                                     n_events = 20)))
  expect_equal(dm$settings$method, "simulation")
  expect_s3_class(dm$simulation, "m_simulation")
  sc <- suppressMessages(m_scenarios(
    pts, link_km = 100, thin_km = 1,
    scenarios = list(delaunay = list(),
                     sim = list(method = "simulation", env = env,
                                simulation = list(n_rep = 3,
                                                  n_events = 20)))))
  expect_equal(sc$comparison$method, c("delaunay", "simulation"))
  expect_error(delimit_m(pts, method = "simulation", link_km = 100),
               "needs `env`")
})

test_that("m_minimal(): growing extents and a stabilisation point", {
  # a stationary landscape: short waves, the same mix of conditions
  # everywhere at a large enough scale
  r <- terra::rast(nrows = 300, ncols = 300, xmin = 0, xmax = 30,
                   ymin = 40, ymax = 70, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  set.seed(1)
  env <- c(terra::rast(r, vals = sin(xy[, 1] * 2.5) +
                         stats::rnorm(nrow(xy), 0, 0.1)),
           terra::rast(r, vals = cos(xy[, 2] * 2.5) +
                         stats::rnorm(nrow(xy), 0, 0.1)))
  names(env) <- c("a", "b")
  pts <- data.frame(lon = stats::rnorm(80, 15, 2),
                    lat = stats::rnorm(80, 55, 1.5))
  mm <- m_minimal(pts, env, steps = 8, beyond = c(1.5, 2, 3, 4), n = 5000)
  t <- mm$extents
  expect_true(all(diff(t$radius_km) > 0))         # monotone extents
  expect_true(all(diff(t$records_pct) >= 0))
  expect_false(is.null(mm$minimal))
  expect_lt(mm$minimal$radius_km, max(t$radius_km))
  expect_gt(t$ks_step[2], t$critical[2])          # small extents differ
  expect_s3_class(mm$m, "sfc")
  expect_output(print(mm), "stabilises")
  expect_s3_class(plot(mm), "ggplot")
  # a continental gradient keeps changing: no stabilisation
  g <- terra::rast(r, vals = xy[, 1] + xy[, 2])
  names(g) <- "g"
  m2 <- m_minimal(pts, g, steps = 6, beyond = c(1.5, 2), n = 5000)
  expect_null(m2$minimal)
  expect_output(print(m2), "did not stabilise")
})

test_that("m_availability() flags a truncated and a skewed predictor", {
  env <- sim_landscape()
  set.seed(3)
  pts <- data.frame(lon = stats::runif(60, 1, 8), lat = stats::runif(60, 44.5,
                                                                     49.5))
  # M that leaves out the north: the records' range of `a` reaches values
  # M does not have
  m <- sf::st_sfc(sf::st_polygon(list(rbind(c(0, 44), c(9, 44), c(9, 46),
                                            c(0, 46), c(0, 44)))),
                  crs = 4326)
  av <- m_availability(env, pts, m = m)
  va <- av$variables
  expect_true(va$distorted[va$layer == "a"])
  expect_lt(va$complete[va$layer == "a"], 1)
  # a skewed layer: most of M low, a thin tail where the records also are
  e2 <- env
  v <- terra::values(env[["a"]])
  terra::values(e2[["a"]]) <- exp(12 * v)
  av2 <- m_availability(e2, pts)
  w <- av2$variables
  expect_true(w$distorted[w$layer == "a"])
  expect_lt(w$well_represented[w$layer == "a"], 0.8)
  # an even layer is not flagged
  flat <- env
  terra::values(flat[["a"]]) <- stats::runif(terra::ncell(flat))
  av3 <- m_availability(flat, pts)
  expect_false(av3$variables$distorted[av3$variables$layer == "a"])
  expect_output(print(av), "distorted")
  expect_s3_class(plot_m_availability(av), "ggplot")
  mm <- suppressMessages(delimit_m(pts, link_km = 150, thin_km = 1))
  rp <- m_report(mm, env = env)
  expect_s3_class(attr(rp, "availability"), "m_availability")
})

test_that("m_simulate() with past climate: land bridges and extinction", {
  # two land masses separated by 2 degrees of sea today; in the past the
  # gap was land (an exposed shelf)
  r <- terra::rast(nrows = 40, ncols = 80, xmin = 0, xmax = 8, ymin = 44,
                   ymax = 48, crs = "EPSG:4326")
  set.seed(5)
  now <- c(terra::rast(r, vals = stats::rnorm(terra::ncell(r), 0, 0.1)),
           terra::rast(r, vals = stats::rnorm(terra::ncell(r), 0, 0.1)))
  names(now) <- c("a", "b")
  past <- now
  xs <- terra::xFromCell(r, seq_len(terra::ncell(r)))
  now[xs > 3 & xs < 5] <- NA
  pts <- data.frame(lon = stats::runif(30, 0.5, 2.5),
                    lat = stats::runif(30, 44.5, 47.5))
  stable <- m_simulate(pts, now, n_rep = 3, n_events = 30, kernel_sd = 3,
                       seed = 1)
  expect_lt(sf::st_bbox(sf::st_transform(stable$m, 4326))[["xmax"]], 3.1)
  cyc <- m_simulate(pts, now, past = past, n_rep = 3, n_events = 3,
                    kernel_sd = 3, seed = 1)
  expect_gt(sf::st_bbox(sf::st_transform(cyc$m, 4326))[["xmax"]], 5)
  expect_gt(cyc$report$accessed_not_today, 0)
  expect_equal(cyc$settings$n_scenarios, 30)
  # M is today's land only
  expect_false(any(!is.na(terra::extract(
    cyc$m_raster, cbind(4, 46))[[1]])))
  expect_output(print(cyc), "glacial cycle")
  # a glacial climate far outside the niche: extinction in every replicate
  cold <- now - 10
  names(cold) <- names(now)
  ext <- m_simulate(pts, now, past = cold, n_rep = 3, n_events = 2,
                    kernel_sd = 1, seed = 1)
  expect_equal(ext$report$extinct_replicates, 3)
  expect_error(m_simulate(pts, now, past = cold[[1]]), "lacks")
})

test_that("m_simulate() takes a raster barrier on another grid", {
  env <- sim_landscape()
  set.seed(2)
  pts <- data.frame(lon = stats::runif(30, 1, 3), lat = stats::runif(30, 46,
                                                                      48))
  br <- terra::rast(nrows = 12, ncols = 18, xmin = 0, xmax = 9, ymin = 44,
                    ymax = 50, crs = "EPSG:4326", vals = 0)
  xs <- terra::xFromCell(br, seq_len(terra::ncell(br)))
  br[xs > 4 & xs < 5.5] <- 1
  ms <- m_simulate(pts, env, barrier = br, n_rep = 3, n_events = 40,
                   kernel_sd = 1, seed = 1)
  expect_lte(sf::st_bbox(sf::st_transform(ms$m, 4326))[["xmax"]], 4 + 1e-6)
})
