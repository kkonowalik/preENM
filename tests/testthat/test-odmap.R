test_that("odmap_table fills the elements from the workflow objects", {
  sim <- simulate_enm_data(seed = 1, n_presence = 150)
  occ <- data.frame(decimalLongitude = c(sim$presence$lon, 0, 15),
                    decimalLatitude = c(sim$presence$lat, 0, 51))
  occ <- rbind(occ, occ[1:5, ])
  fl <- filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
                           verbose = FALSE)
  m <- suppressMessages(delimit_m(sim$presence, link_km = 120, thin_km = 5))
  bg <- suppressMessages(sample_background(sim$env, 300, n_sets = 2,
                                           seed = 1))
  ev <- evaluate_models(sim$models, sim$presence, sim$absence, n_boot = 20,
                        seed = 1)
  sel <- select_models(ev, plot = FALSE)
  tr <- niche_truncation(sim$presence,
                         terra::crop(sim$env, terra::ext(14, 21, 49, 55)),
                         projections = list(whole = sim$env))
  pts <- rbind(cbind(sim$presence, occ = 1), cbind(sim$absence, occ = 0))
  folds <- suppressWarnings(spatial_folds(pts, k = 3, block_size = 100,
                                          occ = pts$occ, seed = 1))
  bs <- suggest_block_size(sim$env, n = 300)
  od <- odmap_table(filter = fl, m = m, records = sim$presence,
                    background = bg, predictors = names(sim$env),
                    truncation = tr, block_size = bs, folds = folds,
                    evaluation = ev, selection = sel,
                    fields = list("Focal Taxon" = "Testus exemplaris",
                                  "Transfer data: Data sources" = "CHELSA"))
  expect_equal(names(od), c("section", "subsection", "element", "value",
                            "filled_by"))
  expect_equal(unique(od$section), c("Overview", "Data", "Model",
                                     "Assessment", "Prediction"))
  expect_false(any(duplicated(od[c("subsection", "element")])))
  get <- function(sub, el) od$value[od$subsection == sub & od$element == el]
  expect_match(get("Biodiversity data", "Sample size"),
               "157 records before cleaning")
  expect_match(get("Biodiversity data", "Sample size"),
               "150 presence records used")
  expect_match(get("Biodiversity data", "Cleaning"), "duplicate \\(5\\)")
  expect_match(get("Biodiversity data", "Clipping"),
               "delaunay, linking distance 120 km")
  expect_match(get("Biodiversity data", "Background data"),
               "300 background points in each of 2 sets")
  expect_match(get("Biodiversity data", "Errors and biases"),
               "truncated on env1 \\(high end\\)")
  expect_match(get("Transfer data", "Quantification of Novelty"),
               "env1 in whole")
  expect_match(get("Data partitioning", "Validation data"),
               "3 folds, blocks of 100 km")
  expect_match(get("Analysis and Correction of non-independence",
                   "Spatial autocorrelation"), "km")
  expect_match(get("Performance statistics", "Performance on test data"),
               "20 stratified bootstrap")
  expect_match(get("Model selection - model averaging - ensembles",
                   "Model selection"), "selected: ")
  expect_match(get("Software", "Software"), "preENM")
  expect_equal(get("Focal Taxon", "Focal Taxon"), "Testus exemplaris")
  expect_equal(get("Transfer data", "Data sources"), "CHELSA")
  expect_equal(od$filled_by[od$element == "Focal Taxon"], "user")
  expect_equal(get("Authorship", "Study title"), "")
  # written to a file
  f <- tempfile(fileext = ".csv")
  expect_invisible(odmap_table(m = m, file = f))
  expect_equal(nrow(utils::read.csv(f)), nrow(od))
  # wrong input
  expect_error(odmap_table(fields = list("Data sources" = "x")),
               "several subsections")
  expect_error(odmap_table(fields = list("Nonsense" = "x")), "Unknown ODMAP")
  expect_error(odmap_table(m = 1), "delimit_m")
})
