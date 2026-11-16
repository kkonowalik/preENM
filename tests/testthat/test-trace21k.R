test_that("trace21k_time converts between ka BP and time ids", {
  tm <- trace21k_time(c(0, 2, 21, 22))
  expect_equal(tm$time_id, c(20L, 0L, -190L, -200L))
  expect_equal(tm$year_start, c(1900, -100, -19100, -20100))
  expect_equal(tm$year_end, c(1990, -1, -19001, -20001))
  expect_equal(trace21k_time(time_id = c(20, -190))$kyr_bp, c(0, 21))
  expect_equal(trace21k_time(11.74)$time_id, -97L)
  expect_error(trace21k_time(23), "between 0 and 22")
  expect_error(trace21k_time(time_id = 21), "whole numbers")
  expect_error(trace21k_time(), "Give")
})

test_that("trace21k_urls builds the file names of the data set", {
  u <- trace21k_urls(kyr_bp = 21)
  expect_equal(nrow(u), 36)
  expect_equal(u$file[1], "CHELSA_TraCE21k_tasmax_1_-190_V1.0.tif")
  expect_equal(u$url[36], paste0(
    "https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_trace/pr/",
    "CHELSA_TraCE21k_pr_12_-190_V1.0.tif"))
  b <- trace21k_urls(c("bio1", "bio12"), kyr_bp = c(0, 21))
  expect_equal(b$file, c("CHELSA_TraCE21k_bio01_20_V1.0.tif",
                         "CHELSA_TraCE21k_bio12_20_V1.0.tif",
                         "CHELSA_TraCE21k_bio01_-190_V1.0.tif",
                         "CHELSA_TraCE21k_bio12_-190_V1.0.tif"))
  expect_true(all(grepl("/chelsa_trace/bio/", b$url)))
  expect_true(all(is.na(b$month)))
  expect_equal(nrow(trace21k_urls("bio", time_id = 0)), 19)
  expect_equal(nrow(trace21k_urls("pr", kyr_bp = 6, months = 6:8)), 3)
  expect_error(trace21k_urls("tas", kyr_bp = 1), "Unknown")
})

test_that("monthly_files finds months written without a leading zero", {
  d <- tempfile()
  dir.create(d)
  for (id in c(20, -190)) for (m in 1:12) {
    file.create(file.path(d, sprintf("CHELSA_TraCE21k_pr_%d_%d_V1.0.tif",
                                     m, id)))
  }
  f <- monthly_files(d, "pr", must_contain = "_-190_V1")
  expect_equal(basename(f), sprintf("CHELSA_TraCE21k_pr_%d_-190_V1.0.tif",
                                    1:12))
  expect_error(monthly_files(d, "pr"), "Several files")
  expect_error(monthly_files(d, "tasmax", "_20_V1"), "Missing")
})

test_that("download_chelsa_trace21k keeps existing files", {
  d <- tempfile()
  dir.create(d)
  u <- trace21k_urls("bio01", kyr_bp = 21)
  file.create(file.path(d, u$file))
  out <- download_chelsa_trace21k(d, "bio01", kyr_bp = 21)
  expect_equal(basename(out), u$file)
})
