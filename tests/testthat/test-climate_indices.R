make_months <- function(crs = "EPSG:4326") {
  r <- terra::rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8, ymin = 30,
                   ymax = 38, crs = "EPSG:4326")
  season <- sin((1:12 - 4) / 12 * 2 * pi)          # warmest in July
  lat <- terra::init(r, "y")
  tmax <- terra::rast(lapply(1:12, function(m) 50 - lat + 10 * season[m]))
  tmin <- tmax - 12
  prec <- terra::rast(lapply(1:12, function(m) {
    20 + (38 - lat) * 10 * pmax(season[(m + 5) %% 12 + 1], 0)  # wet in Jan
  }))
  list(tmax = tmax, tmin = tmin, prec = prec)
}

test_that("radiation matches FAO-56", {
  expect_equal(round(extraterrestrial_radiation(-20 * pi / 180, j = 246), 1),
               32.2)
  expect_equal(extraterrestrial_radiation(75 * pi / 180, 12), 0)
})

test_that("climate_indices computes every index correctly", {
  d <- make_months()
  ci <- climate_indices(d$tmax, d$tmin, d$prec)
  expect_equal(names(ci), climate_index_names())
  expect_equal(names(ci)[1:7], c("vpd", "pet", "cmd", "aridity",
                                 "heat_stress", "sri", "synchronicity"))
  cell <- 1
  tx <- as.numeric(d$tmax[cell])
  tn <- as.numeric(d$tmin[cell])
  p <- as.numeric(d$prec[cell])
  tm <- (tx + tn) / 2
  lat <- terra::yFromCell(d$tmax, cell) * pi / 180
  es <- function(t) 0.6108 * exp(17.27 * t / (t + 237.3))
  expect_equal(ci$vpd[cell][[1]], mean((es(tx) + es(tn)) / 2 - es(tn)))
  days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  ra <- vapply(1:12, function(i) extraterrestrial_radiation(lat, i),
               numeric(1))
  pet <- pmax(0.0023 * 0.408 * ra * (tm + 17.8) * sqrt(tx - tn) * days, 0)
  expect_equal(ci$pet[cell][[1]], sum(pet), tolerance = 1e-6)
  expect_equal(ci$cmd[cell][[1]], sum(pmax(pet - p, 0)), tolerance = 1e-6)
  expect_equal(ci$aridity[cell][[1]], sum(p) / sum(pet), tolerance = 1e-6)
  expect_equal(ci$heat_stress[cell][[1]], sum(pmax(tx - 30, 0)))
  q <- function(m) ((m - 1 + 0:2) %% 12) + 1
  w <- which.max(vapply(1:12, function(m) mean(tm[q(m)]), 1))
  expect_equal(ci$sri[cell][[1]], sum(p[q(w)]) / sum(p))
  dd <- abs(which.max(tx) - which.max(p))
  expect_equal(ci$synchronicity[cell][[1]], (6 - min(dd, 12 - dd)) / 6)
  # warm season dry: rain six months from the heat
  expect_equal(ci$synchronicity[cell][[1]], 0)
  # subset of indices, written to files and re-read
  dir <- tempfile()
  ci2 <- climate_indices(d$tmax, d$tmin, d$prec, indices = c("pet", "sri"),
                         path = dir, prefix = "x_")
  expect_equal(names(ci2), c("pet", "sri"))
  expect_true(file.exists(file.path(dir, "x_pet.tif")))
  expect_equal(terra::values(ci2$pet), terra::values(ci$pet),
               tolerance = 1e-5)
  expect_error(climate_indices(d$tmax, d$tmin, d$prec, indices = "x"),
               "Unknown index")
  expect_error(climate_indices(d$tmax[[1:6]], d$tmin, d$prec), "12 monthly")
})

test_that("temperature units are checked", {
  d <- make_months()
  ref <- climate_indices(d$tmax, d$tmin, d$prec, indices = "heat_stress")
  expect_message(k <- climate_indices(d$tmax + 273.15, d$tmin + 273.15,
                                      d$prec, indices = "heat_stress"),
                 "Kelvin")
  expect_equal(terra::values(k), terra::values(ref), tolerance = 1e-6)
  k10 <- suppressMessages(climate_indices((d$tmax + 273.15) * 10,
                                          (d$tmin + 273.15) * 10, d$prec,
                                          indices = "heat_stress"))
  expect_equal(terra::values(k10), terra::values(ref), tolerance = 1e-5)
  expect_error(climate_indices(d$tmax / 10, d$tmin / 10, d$prec),
               "divided by 10")
  expect_s4_class(climate_indices(d$tmax / 10, d$tmin / 10, d$prec,
                                  temp_units = "C"), "SpatRaster")
  expect_error(climate_indices(d$tmax, d$tmin, d$prec - 1000), "Negative")
})

test_that("projected rasters give the same PET", {
  d <- make_months()
  ref <- climate_indices(d$tmax, d$tmin, d$prec, indices = "pet")
  pr <- lapply(d, function(x) terra::project(x, "EPSG:3035", res = 50000))
  ci <- climate_indices(pr$tmax, pr$tmin, pr$prec, indices = "pet")
  back <- terra::extract(ci, terra::project(terra::vect(cbind(4, 34),
                                                         crs = "EPSG:4326"),
                                             "EPSG:3035"))[[2]]
  expect_equal(back, terra::extract(ref, cbind(4, 34))[[1]], tolerance = 0.02)
})

test_that("monthly files are found and ordered", {
  dir <- tempfile()
  dir.create(dir)
  for (v in c("tas", "tasmax")) for (m in c(12, 1:11)) {
    file.create(file.path(dir, sprintf("CHELSA_%s_%02d_1981-2010_V.2.1.tif",
                                       v, m)))
  }
  f <- monthly_files(dir, "tas")
  expect_length(f, 12)
  expect_match(basename(f[1]), "_tas_01_")
  expect_match(basename(f[12]), "_tas_12_")
  expect_false(any(grepl("tasmax", f)))
  file.create(file.path(dir, "CHELSA_x_ssp126_pr_01_2071_2100_norm.tif"))
  expect_error(monthly_files(dir, "pr"), "Missing pr file")
  expect_error(monthly_files(dir, "tas", must_contain = "ssp585"), "Missing")
})

test_that("CHELSA file lookup", {
  expect_error(chelsa_monthly_urls("tas", "2071-2100"), "gcm")
  expect_equal(lonlat_window(c(1, 2, 3, 4)), c(1, 2, 3, 4))
  expect_error(lonlat_window(c(2, 1, 3, 4)), "xmin")
  v <- terra::vect(cbind(c(4e6, 4.1e6), c(3e6, 3.1e6)), crs = "EPSG:3035")
  w <- lonlat_window(v)
  expect_true(all(abs(w) <= 180))
  skip_on_cran()
  skip_on_ci()
  skip_if_offline("os.zhdk.cloud.switch.ch")
  u <- chelsa_monthly_urls("pr", "2071-2100", gcm = "UKESM1-0-LL",
                           ssp = "ssp370", months = 1:2)
  expect_equal(nrow(u), 2)
  expect_match(u$file[1], "ukesm1-0-ll.*ssp370_pr_01_2071_2100")
})

test_that("the additional indices are computed correctly", {
  d <- make_months()
  ci <- climate_indices(d$tmax, d$tmin, d$prec, indices = "all",
                        gdd_base = 10, frost_threshold = 12)
  expect_equal(names(ci)[8:15],
               c("gdd", "cold_stress", "frost_months", "continentality",
                 "dry_months", "lgp", "emberger", "ellenberg"))
  days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  for (cell in c(1, 40, 64)) {
    tx <- as.numeric(d$tmax[cell])
    tn <- as.numeric(d$tmin[cell])
    p <- as.numeric(d$prec[cell])
    tm <- (tx + tn) / 2
    lat <- terra::yFromCell(d$tmax, cell) * pi / 180
    ra <- vapply(1:12, function(i) extraterrestrial_radiation(lat, i),
                 numeric(1))
    pet <- pmax(0.0023 * 0.408 * ra * (tm + 17.8) * sqrt(tx - tn) * days, 0)
    v <- function(k) ci[[k]][cell][[1]]
    expect_equal(v("gdd"), sum(pmax(tm - 10, 0) * days), tolerance = 1e-6)
    expect_equal(v("cold_stress"), sum(pmax(12 - tn, 0)), tolerance = 1e-6)
    expect_equal(v("frost_months"), sum(tn < 12))
    expect_equal(v("continentality"), max(tm) - min(tm), tolerance = 1e-6)
    expect_equal(v("dry_months"), sum(p < 2 * tm))
    expect_equal(v("lgp"), sum(tm >= 10 & p >= 0.5 * pet))
    expect_equal(v("emberger"), 2000 * sum(p) /
                   ((max(tx) + 273.15)^2 - (min(tn) + 273.15)^2),
                 tolerance = 1e-6)
    expect_equal(v("ellenberg"), 1000 * max(tm) / sum(p), tolerance = 1e-6)
  }
  # the indices vary across the grid (the checks above are not trivial)
  expect_gt(diff(range(terra::values(ci$frost_months))), 0)
  expect_gt(diff(range(terra::values(ci$dry_months))), 0)
})

test_that("Kira's indices, radiation and day length are correct", {
  d <- make_months()
  # colder in the north, so that some months are below the bases
  d$tmax <- d$tmax - 8
  d$tmin <- d$tmin - 8
  ci <- climate_indices(d$tmax, d$tmin, d$prec,
                        indices = c("kira_warmth", "kira_coldness", "light"),
                        gdd_base = 8)
  expect_equal(names(ci), c("kira_warmth", "kira_coldness", "ra", "ra_gs",
                            "daylength_start", "daylength_end",
                            "daylength_gs", "daylight_gs", "daylength_max",
                            "days_to_solstice"))
  days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  mid <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)
  for (cell in c(1, 40, 64)) {
    tm <- (as.numeric(d$tmax[cell]) + as.numeric(d$tmin[cell])) / 2
    lat <- unname(terra::yFromCell(d$tmax, cell)) * pi / 180
    ra <- vapply(1:12, function(i) extraterrestrial_radiation(lat, i),
                 numeric(1))
    v <- function(k) ci[[k]][cell][[1]]
    expect_equal(v("kira_warmth"), sum(pmax(tm - 5, 0)), tolerance = 1e-6)
    expect_equal(v("kira_coldness"), sum(pmin(tm - 5, 0)), tolerance = 1e-6)
    expect_equal(v("ra"), sum(ra * days) / sum(days), tolerance = 1e-6)
    expect_equal(v("ra_gs"), sum((ra * days)[tm >= 8]), tolerance = 1e-6)
    # first month at or above the base after the coldest month
    o <- (which.min(tm) - 1 + 0:11) %% 12 + 1
    s <- o[tm[o] >= 8][1]
    delta <- 0.409 * sin(2 * pi / 365 * mid[s] - 1.39)
    expect_false(is.na(s))
    expect_equal(v("daylength_start"),
                 24 / pi * acos(-tan(lat) * tan(delta)), tolerance = 1e-6)
    # the other day-length indices, from the same months
    dl <- day_length(lat, mid)
    g <- tm >= 8
    e <- rev(o)[tm[rev(o)] >= 8][1]
    expect_equal(v("daylength_end"), dl[e], tolerance = 1e-6)
    expect_equal(v("daylength_gs"), sum((dl * days)[g]) / sum(days[g]),
                 tolerance = 1e-6)
    expect_equal(v("daylight_gs"), sum((dl * days)[g]), tolerance = 1e-5)
    expect_equal(v("daylength_max"), day_length(lat, 172), tolerance = 1e-6)
    expect_gt(v("daylength_max"), v("daylength_gs"))
    # start between the middle of the last cold month and of the first
    # warm one
    before <- (s - 2) %% 12 + 1
    wrap <- function(d) ifelse(d %% 365 > 182.5, d %% 365 - 365, d %% 365)
    lo <- wrap(172 - mid[s])
    hi <- wrap(172 - mid[before])
    expect_gte(v("days_to_solstice"), lo - 1e-6)
    expect_lte(v("days_to_solstice"), hi + 1e-6)
  }
  expect_lt(min(terra::values(ci$kira_coldness)), 0)
  expect_gt(diff(range(terra::values(ci$daylength_start), na.rm = TRUE)), 0)
})

test_that("day length follows the season, the hemisphere and the base", {
  warm_july <- 10 * sin((1:12 - 4) / 12 * 2 * pi)
  f <- function(t, lat) season_start_daylength(rbind(c(t, lat * pi / 180)), 5)
  # mid-March, days before the equinox: a little under 12 hours
  expect_equal(f(c(rep(0, 2), rep(10, 10)), 50), 11.6, tolerance = 0.01)
  # spring start: longer days in April than in March
  expect_gt(f(c(0, 0, 0, rep(10, 9)), 50), f(c(0, 0, rep(10, 10)), 50))
  # southern hemisphere: season starts in the austral spring
  north <- f(4 + warm_july, 45)
  south <- f(4 + warm_july[c(7:12, 1:6)], -45)
  expect_equal(north, south, tolerance = 0.3)
  # always warm enough: the coldest month; never: NA
  expect_equal(f(20 + warm_july, 45), f(c(rep(0, 12)) + 20 + warm_july, 45))
  expect_lt(f(20 + warm_july, 45), 10)             # mid-winter day
  expect_true(is.na(f(rep(0, 12), 45)))
  expect_true(is.na(season_start_daylength(rbind(c(NA, rep(1, 12))), 5)))
})

test_that("day-length indices of the season behave as expected", {
  g <- function(t, lat, what) season_light(rbind(c(t, lat * pi / 180)), 5,
                                           what)
  expect_equal(day_length(0, 100), 12, tolerance = 1e-6)
  expect_equal(day_length(c(70, -70) * pi / 180, 172), c(24, 0))
  expect_equal(day_length(50 * pi / 180, 172), 16.2, tolerance = 0.02)
  # 0 degC until February, 10 degC from March: the base of 5 is reached
  # half-way between mid-February (day 46) and mid-March (day 74)
  t <- c(0, 0, rep(10, 10))
  expect_equal(g(t, 50, "solstice"), 172 - 60)
  # the same season six months later in the south
  expect_equal(g(t[c(7:12, 1:6)], -50, "solstice"), 355 - 242.5)
  # a season that starts after the solstice gives a negative value
  expect_lt(g(c(rep(0, 6), 10, 10, rep(0, 4)), 50, "solstice"), 0)
  # a later start is closer to the solstice
  expect_lt(g(c(0, 0, 0, 0, rep(10, 8)), 50, "solstice"),
            g(t, 50, "solstice"))
  # the season ends in the last warm month before the coldest one
  t2 <- c(-5, 0, 0, 10, 10, 10, 10, 10, 10, 0, 0, 0)
  expect_equal(g(t2, 50, "end"), day_length(50 * pi / 180, 258))
  expect_equal(g(t2, 50, "start"), day_length(50 * pi / 180, 105))
  expect_true(is.na(g(rep(0, 12), 50, "end")))
  expect_true(is.na(g(rep(0, 12), 50, "solstice")))
})

test_that("day-length indices handle cells without a growing season", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 60,
                   ymax = 62, nlyrs = 12)
  cold <- terra::rast(lapply(1:12, function(m) terra::init(r[[1]], 0) - 10))
  ci <- climate_indices(cold + 5, cold - 5, cold * 0 + 30,
                        indices = "light", temp_units = "C")
  expect_true(all(is.na(terra::values(ci$daylength_gs))))
  expect_true(all(terra::values(ci$daylight_gs) == 0))
  expect_true(all(is.na(terra::values(ci$days_to_solstice))))
  expect_true(all(terra::values(ci$daylength_max) > 18))
})

test_that("climate_index_names chooses indices and groups", {
  all <- climate_index_names()
  expect_length(all, 25)
  expect_equal(climate_index_names("cold"),
               c("cold_stress", "frost_months", "kira_coldness"))
  expect_equal(climate_index_names(drop = c("light", "emberger")),
               setdiff(all, c(climate_index_names("light"), "emberger")))
  expect_length(climate_index_names("light"), 8)
  tab <- climate_index_names(table = TRUE)
  expect_equal(names(tab), c("index", "group", "meaning", "unit",
                            "relevance", "reference"))
  expect_true(all(nzchar(unlist(tab))))
  expect_equal(tab$index, all)
  expect_error(climate_index_names("x"), "Unknown group")
  expect_error(climate_index_names(drop = "x"), "Unknown index")
  expect_equal(expand_indices(c("cold", "gdd")),
               c("gdd", "cold_stress", "frost_months", "kira_coldness"))
})
