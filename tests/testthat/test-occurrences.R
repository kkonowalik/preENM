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
