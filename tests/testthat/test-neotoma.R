smp <- data.frame(
  siteid = c("1", "1", "1", "2", "3", "4", "5"),
  sitename = c("A", "A", "A", "B", "C", "D", "E"),
  lat = c(47, 47, 47, 60, 50, 46, 10), long = c(11, 11, 11, 10, 5, 8, 10),
  datasetid = "9",
  datasettype = c("pollen", "pollen", "pollen", "plant macrofossil",
                  "pollen", "pollen", "pollen"),
  variablename = c("Dryas", "Betula", "Dryas", "Dryas octopetala", "Dryas",
                   "Dryas", "Dryas"),
  value = c(3, 10, 2, 1, 0, 4, 2),
  age = c(12500, 12500, 12800, 25000, 9000, NA, 5000),
  ageolder = c(12700, NA, NA, NA, NA, NA, NA),
  ageyounger = c(12300, NA, NA, NA, NA, NA, NA),
  agetype = "Calendar years BP")

test_that("neotoma_records flags, windows and slices samples", {
  r <- neotoma_records(smp, "^Dryas", calibrate_c14 = FALSE)
  expect_equal(nrow(r), 6)                      # Betula dropped
  expect_equal(r$FLAG, c("ok", "ok", "older_than_max", "zero_value",
                         "no_age", "ok"))
  expect_equal(r$AGE_MIN[1], 12300)             # stated bounds
  expect_equal(r$AGE_MAX[2], 13300)             # +/- 500 default
  expect_equal(r$SLICE_MIN[1], trace21k_id(12700))
  expect_equal(r$RECORD_TYPE[3], "macrofossil")
  area <- sf::st_as_sfc("POLYGON ((0 40, 20 40, 20 70, 0 70, 0 40))",
                        crs = 4326)
  ra <- neotoma_records(smp, "^Dryas", area = area, calibrate_c14 = FALSE)
  expect_equal(ra$FLAG[6], "outside_area")
  expect_error(neotoma_records(smp["siteid"], "^Dryas"), "Required columns")
  expect_warning(neotoma_records(smp, "^Saxifraga", calibrate_c14 = FALSE),
                 "No samples")
})

test_that("radiocarbon ages are calibrated", {
  curve <- cbind(0:20000, 0:20000 * 0.9, 20)
  cal <- calibrate_radiocarbon(c(9000, NA, -30), c(50, 50, 50), curve = curve)
  expect_equal(cal$age[1], 10000, tolerance = 10)
  expect_lt(cal$younger[1], cal$age[1])
  expect_gt(cal$older[1], cal$age[1])
  expect_true(is.na(cal$age[2]))
  expect_equal(cal$age[3], -30)
  skip_if_not_installed("rintcal")
  s2 <- smp
  s2$agetype[1] <- "Radiocarbon years BP"
  s2$age[1] <- 10000
  r <- suppressMessages(neotoma_records(s2, "^Dryas"))
  expect_equal(r$AGE[1], 11480, tolerance = 100)   # IntCal20
  expect_equal(r$AGE_SOURCE[1], "IntCal20")
  expect_equal(r$FLAG[1], "ok")
})

test_that("site bins and their slices", {
  r <- neotoma_records(smp, "^Dryas", calibrate_c14 = FALSE)
  sb <- fossil_site_bins(r)
  expect_equal(nrow(sb), 2)
  expect_equal(sb$N_SAMPLES[sb$SITEID == "1"], 2)
  expect_equal(sb$PERIOD, c("Lateglacial", "Holocene"))
  sl <- site_bin_slices(sb)
  expect_equal(nrow(sl), sum(sb$SLICE_MAX - sb$SLICE_MIN + 1))
  expect_equal(trace21k_id(c(0, 2000, 12000, 30000, -70)),
               c(20L, 0L, -100L, -200L, 20L))
  expect_error(fossil_site_bins(r[r$FLAG != "ok", ]), "No usable")
})

test_that("download_neotoma re-uses its cache", {
  f <- tempfile(fileext = ".rds")
  saveRDS(smp, f)
  expect_equal(download_neotoma(1, cache = f), smp)
  if (!requireNamespace("neotoma2", quietly = TRUE)) {
    expect_error(download_neotoma(1), "Install neotoma2")
  }
})
