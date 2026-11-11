# Settings and layers for occurrence filtering -----------------------------------

#' All filtering tests, in the order they are run
#'
#' Each record keeps the reason of the *first* test it fails, so the order
#' matters: cheap, unambiguous metadata tests come first and the climate tests
#' last. Returned by `filter_tests()`; see [filter_occurrences()] for what each
#' test does.
#'
#' @return Character vector of test names.
#' @export
#' @examples
#' filter_tests()
filter_tests <- function() {
  c("duplicate", "grid_cell_centre", "gbif_issue", "basis_of_record",
    "establishment_means", "coord_uncertainty", "coarse_decimals",
    "garden_text", "settlement_text", "locality_is_country", "too_old",
    "manual_exclude", "excluded_dataset", "grid_uncertainty_block",
    "lattice_dataset", "stacked_dataset", "country_mismatch",
    "coordinatecleaner", "outside_range", "at_sea", "no_climate_data",
    "built_up", "cropland", "ghs_dense", "ghs_building_nearby",
    "isolated_climate_outlier", "small_isolated_cluster", "climate_limit",
    "climate_soft_limit")
}

#' Settings for occurrence filtering
#'
#' Collects every threshold used by [filter_occurrences()] and the list of
#' tests to run. Two presets give sensible defaults:
#'
#' * **`"strict"`** - point quality for the presence records used in
#'   modelling: finer coordinate uncertainty, stricter settlement, sea and
#'   cropland limits, the stated-country check and the isolation tests.
#' * **`"lenient"`** - for delimiting the accessible area (the M of the BAM
#'   diagram) or other uses where losing good records is worse than keeping a
#'   few doubtful ones. No country check and no isolation tests.
#'
#' Any setting can be changed through `...`, e.g.
#' `filter_settings("strict", max_coord_unc_m = 500)`. Tests are switched off
#' with `drop_tests` and on with `add_tests`.
#'
#' **Species-specific thresholds have no default.** The climate tests
#' (`climate_limit`, `climate_soft_limit`, `refuge_limit`,
#' `isolation_climate`) depend on the species and are off until you set them;
#' they are in the units of the climate raster you supply (check whether it is
#' degrees C, Kelvin or scaled by 10). For the arctic-alpine *Dryas
#' octopetala* with CHELSA bio10 (mean temperature of the warmest quarter, deg
#' C), the strict values used were `climate_limit = 22`,
#' `climate_soft_limit = 16`, `refuge_limit = 11`, `refuge_radius_km = 10` and
#' `isolation_climate = 14`.
#'
#' @param preset `"strict"` (default) or `"lenient"`.
#' @param ... Settings to change; see the Settings section.
#' @param add_tests,drop_tests Test names (see [filter_tests()]) to add to or
#'   remove from the preset's list of tests.
#' @section Settings:
#' \describe{
#'   \item{`dedup_digits`}{Records of one species at the same location rounded
#'     to this many decimals are duplicates (default 4, about 10 m). `NA` = off.}
#'   \item{`grid_cell_sizes_m`}{Atlas grid squares (m). An uncertainty equal to
#'     half the diagonal of such a square (e.g. 707.1 m for 1 km) marks the
#'     centre of an atlas cell, not a locality.}
#'   \item{`bad_issues`, `bad_basis`, `bad_establishment`}{GBIF `issue`
#'     flags, `basisOfRecord` and `establishmentMeans` values to remove.}
#'   \item{`max_coord_unc_m`, `drop_na_uncertainty`}{Largest accepted
#'     coordinate uncertainty (m; strict 1000, lenient 2000), and whether
#'     records without a stated uncertainty are removed (default `FALSE`).}
#'   \item{`min_coord_decimals`}{Records with *both* coordinates written with
#'     fewer decimals are removed (default 3, i.e. 2 or fewer decimals, about
#'     1 km). 0 = off.}
#'   \item{`garden_regex`, `settlement_regex`}{Regular expressions searched in
#'     locality, habitat and remarks (many languages).}
#'   \item{`min_year`}{Records older than this year are removed (`NA` = off).}
#'   \item{`grid_min_unc_m`, `grid_block_n`}{Grid block rule: at least
#'     `grid_block_n` records of one dataset share the same uncertainty of at
#'     least `grid_min_unc_m` m.}
#'   \item{`lattice_min_coords`, `lattice_ratio`}{Lattice rule: a dataset with
#'     at least `lattice_min_coords` locations whose distinct latitudes and
#'     longitudes are each at most `lattice_ratio` of the locations lies on a
#'     regular mesh.}
#'   \item{`stack_min_records`, `stack_ratio`, `stack_coarse_share`,
#'     `stack_max_specimen_share`}{Stacked rule: a dataset of at least
#'     `stack_min_records` records with at least `stack_ratio` records per
#'     location, mostly coarse or unknown uncertainty and mostly not herbarium
#'     specimens re-uses a list of grid cells or sites.}
#'   \item{`country_buffer_m`, `country_check_skip`}{Tolerance (m) of the
#'     stated-country check and ISO2 codes never checked (default `"SJ"`,
#'     Svalbard, which has no polygon of its own).}
#'   \item{`cc_tests`, `centroids_rad`}{Tests passed to
#'     `CoordinateCleaner::clean_coordinates()` and its centroid radius (m).
#'     The `"urban"` test downloads Natural Earth urban areas (internet
#'     needed) unless you supply them in [filter_layers()].}
#'   \item{`range_buffer_km`, `range_crs`, `protect_in_range`}{Buffer around
#'     the range polygon (km), the projection used for buffering (`"auto"`
#'     chooses a polar projection for circumpolar data, otherwise one centred
#'     on the data), and whether records inside the polygon are exempt from
#'     the isolation tests.}
#'   \item{`sea_unc_factor`, `sea_tol_m`}{A record is at sea when it is
#'     farther from land than `max(sea_unc_factor * uncertainty, sea_tol_m)`.}
#'   \item{`sea_land_pct`, `sea_ring_unc_m`}{Fallback test with a water
#'     raster: on water when no cell within a window (own cell; 3x3 cells
#'     above `sea_ring_unc_m[1]`; 5x5 above `sea_ring_unc_m[2]`) has at least
#'     `sea_land_pct` % land.}
#'   \item{`max_built_pct`, `max_crop_pct`}{Largest accepted built-up and
#'     cropland percentage of the ~1 km cell.}
#'   \item{`ghs_mean_pct`, `ghs_max_m2`}{Largest accepted mean built-up %
#'     (settlement density) and built-up m2 of the most built 100 m cell
#'     within ~1 km (a building nearby), from GHS-BUILT-S.}
#'   \item{`climate_direction`}{`"above"` (default: too warm, too dry, ... is
#'     bad) or `"below"` (too cold is bad).}
#'   \item{`climate_limit`}{Hard limit: records beyond it are removed.}
#'   \item{`climate_soft_limit`, `refuge_radius_km`, `refuge_limit`}{Records
#'     beyond the soft limit survive only if a cell within `refuge_radius_km`
#'     is within `refuge_limit` (e.g. a valley floor below cold mountains).}
#'   \item{`drop_na_climate`}{Remove records without a climate value.}
#'   \item{`isolation_climate`, `support_radius_km`, `min_support`}{Records
#'     beyond `isolation_climate` with fewer than `min_support` other records
#'     within `support_radius_km` are removed.}
#'   \item{`cluster_eps_km`, `min_cluster_n`}{Records chained closer than
#'     `cluster_eps_km` form a cluster; clusters smaller than `min_cluster_n`
#'     are removed.}
#' }
#' @return A list of class `filter_settings`.
#' @export
#' @examples
#' s <- filter_settings("strict", max_coord_unc_m = 500,
#'                      drop_tests = "too_old")
#' s$max_coord_unc_m
#' s$tests
filter_settings <- function(preset = c("strict", "lenient"), ...,
                            add_tests = NULL, drop_tests = NULL) {
  preset <- match.arg(preset)
  common <- list(
    dedup_digits = 4,
    grid_cell_sizes_m = c(1000, 2000, 5000, 10000),
    bad_issues = c("COUNTRY_COORDINATE_MISMATCH", "COORDINATE_INVALID",
                   "COORDINATE_OUT_OF_RANGE", "ZERO_COORDINATE"),
    bad_basis = c("LIVING_SPECIMEN", "FOSSIL_SPECIMEN"),
    bad_establishment = c("MANAGED", "CULTIVATED", "INTRODUCED",
                          "NATURALISED", "INVASIVE"),
    drop_na_uncertainty = FALSE,
    min_coord_decimals = 3,
    garden_regex = paste0(
      "garden|botanic|cultivat|nurser|arboret|greenhouse|glasshouse|",
      "rockery|cemeter|churchyard|Garten|jardin|jard\u00edn|giardin|",
      "ogr\u00f3d|tr\u00e4dg\u00e5rd"),
    settlement_regex = paste0(
      "car[- ]?park|parking|parkplatz|parcheggio|aparcamiento|",
      "estacionamiento|stationnement|parkeringsplass|parkeringsplats|",
      "parkeerplaats|lay-?by|golf (course|club|links)|golfplatz|",
      "camp ?site|campsite|camping|airport|aerodrome|flugplatz|landfill|",
      "picnic (site|area)|service station"),
    min_year = NA,
    grid_min_unc_m = 1500, grid_block_n = 20,
    lattice_min_coords = 30, lattice_ratio = 0.3,
    stack_min_records = 100, stack_coarse_share = 0.5,
    stack_max_specimen_share = 0.5,
    country_check_skip = "SJ",
    cc_tests = c("centroids", "equal", "gbif", "institutions", "zeros",
                 "urban"),
    centroids_rad = 100,
    range_crs = "auto",
    protect_in_range = TRUE,
    sea_ring_unc_m = c(1000, 3000),
    climate_direction = "above",
    climate_limit = NA, climate_soft_limit = NA, refuge_limit = NA,
    isolation_climate = NA,
    drop_na_climate = TRUE,
    min_cluster_n = 3
  )
  presets <- list(
    strict = list(
      max_coord_unc_m = 1000, country_buffer_m = 20000, range_buffer_km = 25,
      sea_unc_factor = 0.25, sea_tol_m = 100, sea_land_pct = 20,
      max_built_pct = 3, max_crop_pct = 30, ghs_mean_pct = 0.5,
      ghs_max_m2 = 800, stack_ratio = 6, refuge_radius_km = 10,
      support_radius_km = 25, min_support = 5, cluster_eps_km = 75,
      tests = filter_tests()),
    lenient = list(
      max_coord_unc_m = 2000, country_buffer_m = 25000, range_buffer_km = 50,
      sea_unc_factor = 0.5, sea_tol_m = 500, sea_land_pct = 5,
      max_built_pct = 10, max_crop_pct = 60, ghs_mean_pct = 1,
      ghs_max_m2 = 2000, stack_ratio = 10, refuge_radius_km = 15,
      support_radius_km = 50, min_support = 2, cluster_eps_km = 100,
      tests = setdiff(filter_tests(), c("country_mismatch",
                                        "isolated_climate_outlier",
                                        "small_isolated_cluster")))
  )
  s <- c(common, presets[[preset]])
  user <- list(...)
  unknown <- setdiff(names(user), names(s))
  if (length(unknown)) {
    stop("Unknown setting(s): ", paste(unknown, collapse = ", "),
         ". See ?filter_settings.", call. = FALSE)
  }
  s[names(user)] <- user
  bad_tests <- setdiff(c(add_tests, drop_tests, s$tests), filter_tests())
  if (length(bad_tests)) {
    stop("Unknown test(s): ", paste(bad_tests, collapse = ", "),
         ". See filter_tests().", call. = FALSE)
  }
  s$tests <- filter_tests()[filter_tests() %in%
                              setdiff(union(s$tests, add_tests), drop_tests)]
  s$climate_direction <- match.arg(s$climate_direction, c("above", "below"))
  s$preset <- preset
  class(s) <- "filter_settings"
  s
}

#' @export
print.filter_settings <- function(x, ...) {
  cat("Occurrence filtering settings (preset: ", x$preset, ")\n", sep = "")
  cat("Tests (", length(x$tests), "): ", paste(x$tests, collapse = ", "),
      "\n\n", sep = "")
  y <- unclass(x)[setdiff(names(x), c("tests", "preset", "garden_regex",
                                      "settlement_regex"))]
  val <- vapply(y, function(v) paste(format(v), collapse = ", "),
                character(1))
  print(data.frame(setting = names(val), value = unname(val)),
        row.names = FALSE, right = FALSE)
  invisible(x)
}

#' External layers for occurrence filtering
#'
#' Collects the optional spatial data used by the layer-based tests of
#' [filter_occurrences()]. A test runs only when its layer is supplied. Each
#' argument accepts a file path or an object already in memory.
#'
#' @param range Species range polygon(s) (`sf`, [terra::SpatVector] or a path
#'   to e.g. a shapefile), any CRS. Used by `outside_range` and to exempt
#'   records inside the range from the isolation tests.
#' @param land Detailed land polygons, e.g. the GSHHG shoreline from
#'   [download_shoreline()] or OpenStreetMap land polygons. Used by `at_sea`.
#' @param water Raster of the water fraction or percentage of ~1 km cells
#'   (e.g. `download_landcover("water")`). Coarse fallback for `at_sea` when
#'   no `land` polygons are given, and used to pick candidate records for the
#'   shoreline test when both are given.
#' @param climate Climate raster for the climate tests (e.g. CHELSA bio10),
#'   in the units your thresholds use.
#' @param built,cropland Rasters of built-up and cropland fraction or
#'   percentage (e.g. from [download_landcover()]). Fractions (0-1) are
#'   converted to percent automatically.
#' @param ghs_mean,ghs_max GHS-BUILT-S derived rasters: mean and maximum
#'   built-up m2 per 100 m cell within ~1 km cells (see [aggregate_ghs()]).
#' @param urban Urban-area polygons for the CoordinateCleaner `urban` test.
#'   If `NULL`, Natural Earth urban areas are downloaded when that test runs.
#' @return A list of class `filter_layers`.
#' @export
#' @examples
#' \dontrun{
#' lay <- filter_layers(range = "range.shp",
#'                      land = download_shoreline("data/coast"),
#'                      climate = "CHELSA_bio10_1981-2010_V.2.1.tif",
#'                      built = download_landcover("built", "data"))
#' }
filter_layers <- function(range = NULL, land = NULL, water = NULL,
                          climate = NULL, built = NULL, cropland = NULL,
                          ghs_mean = NULL, ghs_max = NULL, urban = NULL) {
  out <- list(range = as_vector_layer(range, "range"),
              land = as_vector_layer(land, "land"),
              water = as_raster_layer(water, "water"),
              climate = as_raster_layer(climate, "climate"),
              built = as_raster_layer(built, "built"),
              cropland = as_raster_layer(cropland, "cropland"),
              ghs_mean = as_raster_layer(ghs_mean, "ghs_mean"),
              ghs_max = as_raster_layer(ghs_max, "ghs_max"),
              urban = urban)
  class(out) <- "filter_layers"
  out
}

#' @export
print.filter_layers <- function(x, ...) {
  have <- vapply(unclass(x), function(v) !is.null(v), logical(1))
  cat("Filtering layers supplied: ",
      if (any(have)) paste(names(x)[have], collapse = ", ") else "none",
      "\n", sep = "")
  invisible(x)
}

as_raster_layer <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) x <- terra::rast(x)
  if (inherits(x, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    x <- terra::rast(x)
  }
  if (!inherits(x, "SpatRaster")) {
    stop("`", what, "` must be a SpatRaster or a raster file path.",
         call. = FALSE)
  }
  if (terra::nlyr(x) > 1) {
    warning("`", what, "` has several layers; the first is used.",
            call. = FALSE)
    x <- x[[1]]
  }
  x
}

as_vector_layer <- function(x, what) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) x <- terra::vect(x)
  if (inherits(x, c("sf", "sfc"))) x <- terra::vect(x)
  if (!inherits(x, "SpatVector")) {
    stop("`", what, "` must be an sf object, a SpatVector or a file path.",
         call. = FALSE)
  }
  if (is.na(terra::crs(x)) || terra::crs(x) == "") {
    warning("`", what, "` has no CRS; assuming EPSG:4326.", call. = FALSE)
    terra::crs(x) <- "EPSG:4326"
  }
  x
}
