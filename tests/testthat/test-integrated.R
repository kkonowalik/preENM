isdm_records <- function() {
  po <- data.frame(lon = c(10.1, 10.4, 11.2), lat = c(46.1, 46.3, 46.8),
                   species = "Aus bus")
  pa <- data.frame(lon = c(10.2, 10.6, 11.0, 11.3),
                   lat = c(46.2, 46.5, 46.6, 46.9), species = "Aus bus",
                   present = c(1, 0, 0, 1), visits = c(2, 1, 3, 1))
  ab <- data.frame(lon = c(10.2, 10.7), lat = c(46.2, 46.5),
                   species = "Aus bus", individualCount = c(4, 0))
  suppressMessages(merge_occurrences(
    GBIF = po, survey = pa, counts = ab,
    data_type = c(survey = "presence_absence", counts = "abundance")))
}

test_that("generic export: one dataset per type with coordinates and response", {
  ex <- export_integrated(isdm_records())
  expect_s3_class(ex, "integrated_data")
  expect_equal(names(ex$datasets), c("presence_only", "presence_absence",
                                     "abundance"))
  for (d in ex$datasets) {
    expect_true(all(c("species", "x", "y", "response", "data_type",
                      "source") %in% names(d)))
  }
  expect_equal(ex$datasets$presence_absence$response, c(1, 0, 0, 1))
  expect_equal(ex$datasets$abundance$response, c(4, 0))
  expect_equal(ex$summary$prevalence, c(NA, 0.5, 0.5))
  expect_output(print(ex), "IN SHORT")
})

test_that("intSDM export has the arguments of addStructured()", {
  ex <- export_integrated(isdm_records(), "intSDM", trials = "visits")
  expect_equal(unname(vapply(ex$export, `[[`, "", "datasetType")),
               c("PO", "PA", "Counts"))
  for (e in ex$export) {
    expect_true(all(c("dataStructured", "datasetType", "speciesName",
                      "datasetName") %in% names(e)))
    expect_s3_class(e$dataStructured, "sf")
    expect_true(e$speciesName %in% names(e$dataStructured))
  }
  expect_null(ex$export$presence_only$responseName)
  expect_equal(ex$export$presence_absence$responseName, "response")
  expect_equal(ex$export$presence_absence$trialsName, "trials")
  expect_equal(ex$export$presence_absence$dataStructured$trials,
               c(2, 1, 3, 1))
})

test_that("PointedSDMs export: response columns name the type", {
  ex <- export_integrated(isdm_records(), "PointedSDMs", crs = "EPSG:3035")
  d <- ex$export$data
  expect_true(all(vapply(d, inherits, logical(1), "sf")))
  # presence-only: no response column, so PointedSDMs sees it as PO
  expect_false(any(c("present", "counts", "response") %in%
                     names(d$presence_only)))
  expect_equal(d$presence_absence$present, c(1, 0, 0, 1))
  expect_equal(d$abundance$counts, c(4, 0))
  expect_equal(ex$export$args$responsePA, "present")
  expect_equal(ex$export$args$responseCounts, "counts")
  expect_equal(ex$export$args$Projection, "EPSG:3035")
  expect_equal(sf::st_crs(d$abundance), sf::st_crs("EPSG:3035"))
})

test_that("biomod2 export: one data set per type and species", {
  ex <- export_integrated(isdm_records(), "biomod2")
  e <- ex$export
  expect_equal(length(e), 3)
  expect_equal(unname(vapply(e, `[[`, "", "data.type")),
               c("binary", "binary", "count"))
  for (k in e) {
    expect_true(all(c("resp.name", "resp.var", "resp.xy", "data.type") %in%
                      names(k)))
    expect_equal(nrow(k$resp.xy), length(k$resp.var))
    expect_equal(k$resp.name, "Aus.bus")
  }
  ab <- isdm_records()
  ab$response[ab$data_type == "abundance"] <- c(2.5, 0)
  expect_equal(export_integrated(ab, "biomod2")$export[[3]]$data.type,
               "abundance")
})

test_that("absences in presence-only data are dropped; labels need a type", {
  x <- isdm_records()
  x$response[1] <- 0
  expect_message(ex <- export_integrated(x), "absence")
  expect_equal(ex$summary$records[1], 2)
  x$data_type[x$source == "GBIF"] <- "citizen"
  expect_error(export_integrated(x), "observation model")
  ex2 <- suppressMessages(export_integrated(
    x, types = c(citizen = "presence_only")))
  expect_equal(ex2$datasets[[1]]$data_type[1], "presence_only")
  expect_equal(names(export_integrated(isdm_records(),
                                       split = "source")$datasets),
               c("GBIF", "survey", "counts"))
  # a presence-absence dataset without absences is flagged
  y <- isdm_records()
  y$response[y$data_type == "presence_absence"] <- 1
  expect_output(print(export_integrated(y)), "no absences")
  f <- write_integrated(ex2, tempfile())
  expect_length(f, 4)
})
