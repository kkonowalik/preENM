# Taxonomic complexity

backbone <- function() {
  data.frame(
    taxonID = paste0("t", 1:9),
    scientificName = c("Aus bus", "Aus cus", "Aus bus var. dus",
                       "Aus bus subsp. eus", "Xus bus", "Yus bus", "Zus bus",
                       "Aus fus", "Aus cus var. gus"),
    taxonRank = c("SPECIES", "SPECIES", "VARIETY", "SUBSPECIES", "SPECIES",
                  "SPECIES", "SPECIES", "SPECIES", "VARIETY"),
    taxonomicStatus = c("Accepted", "Accepted", "Accepted", "Accepted",
                        "Synonym", "Synonym", "Synonym", "Accepted",
                        "Synonym"),
    acceptedNameUsageID = c("", "", "", "", "t1", "t1", "t1", "", "t2"),
    parentNameUsageID = c("", "", "t1", "t1", "", "", "", "", ""),
    stringsAsFactors = FALSE)
}

test_that("synonyms and infraspecific taxa are counted", {
  tc <- taxon_complexity(c("Aus bus", "Aus cus", "Aus fus", "Xus bus",
                           "Nobody"), backbone())
  # Aus bus: 3 synonyms (Xus, Yus, Zus), 2 accepted infraspecific taxa
  expect_equal(tc$n_synonyms, c(3, 1, 0, 3, NA))
  expect_equal(tc$n_infraspecific, c(2, 0, 0, 2, NA))
  # a synonym is followed to its accepted name
  expect_equal(tc$matched, c("Aus bus", "Aus cus", "Aus fus", "Aus bus", NA))
  expect_equal(tc$reliability[3], 1)
  expect_equal(tc$reliability[1], 1 / (1 + log(1 + 5)))
  # more complex names get lower weights
  expect_true(tc$reliability[1] < tc$reliability[2])
  only_syn <- taxon_complexity("Aus bus", backbone(), count = "synonyms")
  expect_equal(only_syn$reliability, 1 / (1 + log(1 + 3)))
  expect_error(taxon_complexity("Aus bus", backbone()[-1]), "taxonID")
})

test_that("fuzzy matching with WorldFlora agrees with exact matching", {
  skip_if_not_installed("WorldFlora")
  ex <- taxon_complexity(c("Aus bus", "Aus cus"), backbone())
  fz <- suppressWarnings(taxon_complexity(c("Aus bus", "Aus cus"),
                                          backbone(), fuzzy = TRUE))
  expect_equal(fz$n_synonyms, ex$n_synonyms)
})

test_that("record_ignorance takes taxon_complexity weights", {
  grid <- terra::rast(xmin = 0, xmax = 2, ymin = 0, ymax = 1,
                      resolution = 1, crs = "EPSG:4326")
  rec <- data.frame(lon = rep(c(0.5, 1.5), each = 6), lat = 0.5,
                    species = rep(c("Aus bus", "Aus fus"), each = 6))
  tc <- taxon_complexity(c("Aus bus", "Aus fus"), backbone())
  ig <- record_ignorance(rec, grid, reliability = tc, year = FALSE,
                         uncertainty = FALSE, min_records = 5,
                         radius_km = 10)
  r <- ig$records
  expect_equal(unique(r$reliability[r$species == "Aus bus"]),
               tc$reliability[1])
  # the taxonomically simple species' cell is the better known
  v <- terra::values(ig$raster)[, 1]
  expect_lt(v[2], v[1])
})
