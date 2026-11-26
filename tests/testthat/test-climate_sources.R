cs_data <- function() {
  sim <- simulate_enm_data(seed = 1, n_presence = 60)
  list(env = sim$env, occ = sim$presence)
}

test_that("identical products show zero disagreement", {
  d <- cs_data()
  cs <- climate_source_compare(list(A = d$env, B = d$env), d$occ, n = 400,
                               n_perm = 49)
  m <- cs$metrics
  expect_equal(m$offset, rep(0, nrow(m)))
  expect_equal(m$rmse, rep(0, nrow(m)))
  expect_equal(m$ccc, rep(1, nrow(m)))
  expect_equal(m$r, rep(1, nrow(m)))
  expect_equal(m$overlap, rep(1, nrow(m)))
  expect_equal(m$coverage, rep(1, nrow(m)))
  expect_equal(m$large, rep(0, nrow(m)))
  expect_true(all(cs$exposure$grade == "almost perfect"))
  expect_output(print(cs), "identical values")
})

test_that("a known offset and scale are recovered", {
  d <- cs_data()
  b <- d$env
  b[[1]] <- b[[1]] + 2
  b[[2]] <- b[[2]] * 2
  cs <- climate_source_compare(list(A = d$env, B = b), d$occ, n = 400,
                               n_perm = 0)
  m <- cs$metrics
  e1 <- m[m$variable == "env1", ]
  expect_equal(e1$offset, c(2, 2), tolerance = 1e-6)
  expect_equal(e1$rmse, c(2, 2), tolerance = 1e-6)
  expect_equal(e1$r, c(1, 1), tolerance = 1e-9)
  # Lin's CCC by hand at the records
  a <- cs$points$env1__A[cs$points$domain == "records"]
  bb <- a + 2
  ccc <- 2 * mean((a - mean(a)) * (bb - mean(bb))) /
    (mean((a - mean(a))^2) + mean((bb - mean(bb))^2) + 4)
  expect_equal(e1$ccc[e1$domain == "records"], ccc, tolerance = 1e-9)
  e2 <- m[m$variable == "env2", ]
  expect_equal(e2$r, c(1, 1), tolerance = 1e-9)
  expect_true(all(e2$ccc < 1))
  # an offset of several SD is reported as a unit mismatch
  expect_true(cs$exposure$unit_mismatch[cs$exposure$variable == "env1"])
  expect_output(print(cs), "units or scale")
})

test_that("a regional bias is located, and the helpers are correct", {
  d <- cs_data()
  b <- d$env
  x <- terra::init(b[[1]], "x")
  b[[1]] <- b[[1]] + 0.05 * (x - 14)   # grows eastwards
  cs <- climate_source_compare(list(A = d$env, B = b), n = 1000, n_perm = 0)
  g <- cs$geography[cs$geography$variable == "env1", ]
  expect_gt(g$rho_lon, 0.9)
  expect_equal(cs$settings$n_records, 0)
  expect_equal(unique(cs$metrics$domain), "M")
  # helpers
  expect_equal(harmonize_bio(c("bio_01", "wc2.1_30s_bio_12",
                               "CHELSA_bio19_1981-2010_V.2.1", "BIO5",
                               "bio25", "tmean")),
               c("bio1", "bio12", "bio19", "bio5", "bio25", "tmean"))
  expect_equal(ccc_grade(c(0.5, 0.92, 0.97, 0.995)),
               c("poor", "moderate", "substantial", "almost perfect"))
  set.seed(1)
  z <- stats::rnorm(500)
  expect_equal(overlap_index(z, z), 1)
  expect_lt(overlap_index(z, z + 10), 0.01)
  # Mantel: the same correlation matrix gives r = 1
  cm <- stats::cor(matrix(stats::rnorm(400), 100))
  expect_equal(mantel_cor(cm, cm, 0, 1)[1], 1)
})

test_that("names are harmonised, and inputs are checked", {
  d <- cs_data()
  a <- d$env
  names(a) <- c("bio_01", "bio_12")
  b <- d$env
  names(b) <- c("CHELSA_bio1_1981", "CHELSA_bio12_1981")
  cs <- climate_source_compare(list(WC = a, CH = b), d$occ, n = 200,
                               n_perm = 0)
  expect_equal(cs$settings$variables, c("bio1", "bio12"))
  expect_error(climate_source_compare(list(A = d$env)), "two or more")
  expect_error(climate_source_compare(list(A = d$env, B = b), n = 100),
               "not found in every product")
  expect_error(climate_source_compare(list(A = d$env, B = d$env),
                                      reference = "Z"), "reference")
})

test_that("plots and files are produced", {
  d <- cs_data()
  b <- d$env * 1.1
  cs <- climate_source_compare(list(A = d$env, B = b), d$occ, n = 300,
                               n_perm = 0)
  for (t in c("summary", "map", "scatter", "density")) {
    expect_s3_class(plot_climate_sources(cs, t), "ggplot")
  }
  f <- tempfile(fileext = ".png")
  plot_climate_sources(cs, "map", file = f, width = 6, height = 4)
  expect_true(file.exists(f))
  expect_error(plot_climate_sources(cs, "map", variable = "zz"), "variable")
  dir <- tempfile("cs")
  out <- write_climate_sources(cs, dir)
  expect_true(all(file.exists(out)))
  expect_equal(nrow(utils::read.csv(out[1])), nrow(cs$metrics))
})
