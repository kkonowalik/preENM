# A small synthetic GBIF-like download in Austria with one problem of each
# kind, plus matching layers. Used by the filtering tests.
make_fake_gbif <- function() {
  set.seed(42)
  n_good <- 60
  good <- data.frame(
    gbifID = as.character(1000 + seq_len(n_good)),
    datasetKey = "ds-good", institutionCode = "HERB",
    occurrenceID = paste0("HERB:", seq_len(n_good)),
    species = "Dryas octopetala",
    decimalLongitude = sprintf("%.5f", runif(n_good, 13.5, 15.5)),
    decimalLatitude = sprintf("%.5f", runif(n_good, 47.2, 47.8)),
    coordinateUncertaintyInMeters = "30",
    basisOfRecord = "HUMAN_OBSERVATION", establishmentMeans = NA,
    issue = NA, countryCode = "AT", locality = "Rocky slope",
    habitat = NA, occurrenceRemarks = NA, year = "2015",
    stringsAsFactors = FALSE)
  bad <- function(id, lon, lat, ...) {
    r <- good[1, ]
    r$gbifID <- as.character(id)
    r$decimalLongitude <- lon
    r$decimalLatitude <- lat
    v <- list(...)
    r[names(v)] <- v
    r
  }
  probs <- rbind(
    bad(1, "abc", "47.5"),                                       # invalid
    bad(2, good$decimalLongitude[1], good$decimalLatitude[1]),   # duplicate
    bad(3, "14.01111", "47.31111",
        coordinateUncertaintyInMeters = "707.1"),                # grid centre
    bad(4, "14.02222", "47.32222", issue = "ZERO_COORDINATE"),
    bad(5, "14.03333", "47.33333", basisOfRecord = "LIVING_SPECIMEN"),
    bad(6, "14.04444", "47.34444", establishmentMeans = "INTRODUCED"),
    bad(7, "14.05555", "47.35555",
        coordinateUncertaintyInMeters = "5000"),                 # uncertainty
    bad(8, "14.1", "47.4"),                                      # coarse
    bad(9, "14.06666", "47.36666", locality = "Botanical garden"),
    bad(10, "14.07777", "47.37777", locality = "near the car park"),
    bad(11, "14.08888", "47.38888", locality = "Austria"),
    bad(12, "14.13579", "47.46802", year = "1850"),               # too old
    bad(13, "14.11111", "47.41111"),                              # manual
    bad(14, "14.95000", "47.45000", countryCode = "IT"),          # country
    bad(15, "17.50000", "47.50000", countryCode = NA),            # off range
    bad(16, "13.10000", "47.50000", countryCode = NA),            # at sea
    bad(17, "14.12345", "47.42345"),                              # built up
    bad(18, "14.23456", "47.43456"),                              # too warm
    bad(19, "14.34567", "47.44567", countryCode = "DE"))          # kept by id
  x <- rbind(good, probs)
  # a gridded atlas dataset: 36 records on a regular 0.1 degree lattice
  lat <- rep(seq(47.21, 47.71, by = 0.1), each = 6)
  lon <- rep(seq(13.61, 14.11, by = 0.1), 6)
  atlas <- good[rep(1, 36), ]
  atlas$gbifID <- as.character(5000 + 1:36)
  atlas$datasetKey <- "ds-atlas"
  atlas$institutionCode <- "ATLAS"
  atlas$decimalLongitude <- sprintf("%.5f", lon)
  atlas$decimalLatitude <- sprintf("%.5f", lat)
  atlas$coordinateUncertaintyInMeters <- NA
  rbind(x, atlas)
}

make_fake_layers <- function() {
  # land everywhere except a "sea" west of 13.3 E (inside the range)
  land <- terra::vect("POLYGON ((13.3 46, 18 46, 18 49, 13.3 49, 13.3 46))",
                      crs = "EPSG:4326")
  range <- terra::vect(
    "POLYGON ((13 47, 16 47, 16 48, 13 48, 13 47))", crs = "EPSG:4326")
  r <- terra::rast(xmin = 12, xmax = 18, ymin = 46, ymax = 49, res = 0.02,
                   crs = "EPSG:4326")
  clim <- terra::init(r, 10)
  clim[terra::cellFromXY(clim, cbind(14.23456, 47.43456))] <- 25  # too warm
  built <- terra::init(r, 0)
  built[terra::cellFromXY(built, cbind(14.12345, 47.42345))] <- 0.5  # 50 %
  list(land = land, range = range, climate = clim, built = built)
}
