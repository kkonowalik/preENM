test_that("read_gbif reads tab-delimited text and zip files as character", {
  td <- tempfile()
  dir.create(td)
  f <- file.path(td, "0001.csv")
  writeLines(c(paste("gbifID", "decimalLatitude", "decimalLongitude",
                     "locality", sep = "\t"),
               paste("1", "47.10", "14.5", "a \"quoted\" place", sep = "\t"),
               paste("2", "47.2", "14.6", "", sep = "\t")), f)
  x <- read_gbif(f)
  expect_true(all(vapply(x, is.character, logical(1))))
  expect_equal(x$decimalLatitude[1], "47.10")     # decimals kept as written
  expect_equal(x$locality[1], "a \"quoted\" place")
  expect_true(is.na(x$locality[2]))
  z <- file.path(td, "0001.zip")
  old <- setwd(td)
  on.exit(setwd(old))
  utils::zip(z, "0001.csv", flags = "-q")
  expect_equal(read_gbif(z), x)
  expect_error(read_gbif(file.path(td, "missing.csv")), "not found")
})

test_that("merge_occurrences standardises columns and drops bad rows", {
  a <- data.frame(Longitude = c(10.1, 10.2), Latitude = c(46.1, 46.2),
                  sample = c("s1", "s2"))
  b <- data.frame(decimalLongitude = c("10.1", "200"),
                  decimalLatitude = c("46.1", "46"), species = "Dryas",
                  gbifID = c("1", "2"))
  m <- suppressMessages(merge_occurrences(own = a, GBIF = b))
  expect_equal(names(m)[1:4], c("species", "lon", "lat", "source"))
  expect_equal(nrow(m), 3)                # lon 200 dropped
  expect_equal(m$source, c("own", "own", "GBIF"))
  expect_equal(m$species, c("unknown", "unknown", "Dryas"))
  expect_true(all(c("sample", "gbifID") %in% names(m)))
  m2 <- suppressMessages(merge_occurrences(own = a, species = "Dryas o."))
  expect_equal(unique(m2$species), "Dryas o.")
  expect_error(merge_occurrences(a), "name")
  f <- tempfile(fileext = ".csv")
  utils::write.csv2(a, f, row.names = FALSE)       # semicolon-separated
  m3 <- suppressMessages(merge_occurrences(file = f))
  expect_equal(m3$lon, a$Longitude)
})

test_that("invalid UTF-8 text is repaired", {
  bad <- "Krak\xf3w"                              # latin1 bytes for o-acute
  Encoding(bad) <- "unknown"
  expect_false(validUTF8(bad))
  a <- data.frame(lon = 19.9, lat = 50.06, place = bad)
  expect_message(m <- merge_occurrences(own = a), "repaired")
  expect_true(validUTF8(m$place))
  expect_equal(m$place, "Kraków")
})

test_that("parse_names() strips authors and markers, flags hybrids", {
  p <- parse_names(c("Aus bus L.", "Aus bus subsp. eus (Hook.) Sm.",
                     "Aus bus x cus", "Aus x hybridus", "Aus bus agg.",
                     "Aus sp.", "Aus bus L. ex Sm."))
  expect_equal(p$canonical, c("Aus bus", "Aus bus eus", "Aus bus",
                              "Aus hybridus", "Aus bus", "Aus", "Aus bus"))
  expect_equal(p$hybrid, c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(p$aggregate, c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE,
                              FALSE))
})

test_that("match_names(): synonym, fuzzy, genus-only and unmatched names", {
  occ <- data.frame(lon = 1:7, lat = 50,
                    scientificName = c("Aus bus", "Aus dus Sm.", "Aus buss",
                                       "Aus sp.", "Aus bus subsp. eus",
                                       "Zus bus", NA))
  m <- match_names(occ, backbone = name_backbone(), verbose = FALSE)
  expect_equal(nrow(m), 7)
  # a known synonym resolves to its accepted name
  expect_equal(m$name_accepted[2], "Aus bus")
  expect_true(m$name_synonym[2])
  expect_equal(m$name_status[2], "synonym")
  # a misspelled name is flagged as fuzzy
  expect_equal(m$name_match_type[3], "fuzzy")
  expect_true(m$name_fuzzy[3])
  expect_equal(m$name_accepted[3], "Aus bus")
  expect_lt(m$name_confidence[3], 100)
  # genus only
  expect_equal(m$name_rank[4], "genus")
  expect_true(is.na(m$name_accepted_species[4]))
  # infraspecific: its species
  expect_equal(m$name_rank[5], "subspecies")
  expect_equal(m$name_accepted_species[5], "Aus bus")
  # unknown genus and missing name
  expect_equal(m$name_match_type[6:7], c("none", "none"))
  expect_true(all(m$name_unmatched[6:7]))
  expect_equal(nrow(attr(m, "names_matched")), 6)
  # a species of a known genus that is not in the backbone: higher rank
  h <- match_names(data.frame(species = "Aus novus"),
                   backbone = name_backbone(), verbose = FALSE)
  expect_equal(h$name_match_type, "higherrank")
  expect_equal(h$name_rank, "genus")
})

test_that("online look-ups are cached", {
  calls <- 0
  local_mocked_bindings(match_gbif = function(u, kingdom, checklist) {
    calls <<- calls + length(u)
    match_rows(u, name_matched = u, name_accepted = u,
               name_accepted_species = u, name_rank = "SPECIES",
               name_status = "ACCEPTED", name_match_type = "EXACT",
               name_confidence = 99)
  })
  rm(list = ls(.name_cache), envir = .name_cache)
  f <- tempfile(fileext = ".rds")
  d <- data.frame(species = c("Qus aus", "Qus bus", "Qus aus"))
  a <- match_names(d, cache_file = f, verbose = FALSE)
  expect_equal(calls, 2)
  b <- match_names(d, cache_file = f, verbose = FALSE)
  expect_equal(calls, 2)
  expect_equal(a$name_accepted, b$name_accepted)
  expect_true(file.exists(f))
  # a new session reads the file
  rm(list = ls(.name_cache), envir = .name_cache)
  attr(.name_cache, "loaded_from") <- NULL
  match_names(d, cache_file = f, verbose = FALSE)
  expect_equal(calls, 2)
  rm(list = ls(.name_cache), envir = .name_cache)
  attr(.name_cache, "loaded_from") <- NULL
})

test_that("match_names() against GBIF", {
  skip_if_not_installed("rgbif")
  skip_on_cran()
  skip_on_ci()
  skip_if_offline()
  m <- match_names(data.frame(species = c("Potentilla tabernaemontani",
                                          "Leontodn hispidus", "Dryas")),
                   kingdom = "Plantae", cache = FALSE, verbose = FALSE)
  expect_equal(m$name_accepted[1], "Potentilla verna")
  expect_true(m$name_synonym[1])
  expect_true(m$name_fuzzy[2])
  expect_equal(m$name_rank[3], "genus")
})

test_that("the multiplication sign marks a hybrid", {
  p <- parse_names(c("Salix × rubens", "×Festulolium loliaceum"))
  expect_equal(p$canonical, c("Salix rubens", "Festulolium loliaceum"))
  expect_equal(p$hybrid, c(TRUE, TRUE))
})
