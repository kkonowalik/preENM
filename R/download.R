# Optional downloads of filtering layers ------------------------------------------

#' Download ~1 km land-cover layers (built-up, cropland, water)
#'
#' Downloads global ~1 km (30 arc-second) fractions of built-up area,
#' cropland or open water derived from ESA WorldCover, via
#' `geodata::landcover()` (needs the `geodata` package and internet; files are
#' cached in `path`). Use them in [filter_layers()] as `built`, `cropland`
#' and `water`.
#'
#' @param var One of `"built"`, `"cropland"`, `"water"` (or any other
#'   `geodata::landcover()` variable).
#' @param path Folder where the file is stored (and re-used next time).
#' @return A [terra::SpatRaster].
#' @export
#' @examples
#' \dontrun{
#' built <- download_landcover("built", "data")
#' water <- download_landcover("water", "data")
#' }
download_landcover <- function(var = c("built", "cropland", "water"), path) {
  if (!requireNamespace("geodata", quietly = TRUE)) {
    stop("Install the geodata package: install.packages(\"geodata\")",
         call. = FALSE)
  }
  var <- var[1]
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  geodata::landcover(var, path = path)
}

#' Download a detailed shoreline (GSHHG land polygons)
#'
#' Downloads the Global Self-consistent, Hierarchical, High-resolution
#' Geography (GSHHG) shoreline v2.3.7 (about 150 MB, once) and unzips the
#' land polygons (level 1) at the chosen resolution. The result is the
#' `land` layer of [filter_layers()] for the `at_sea` test. Any other land
#' polygon file (e.g. OpenStreetMap land polygons) works too.
#'
#' @param path Folder for the download (re-used next time).
#' @param resolution `"h"` high (~200 m, default), `"f"` full, `"i"`
#'   intermediate, `"l"` low or `"c"` crude.
#' @param urls Download addresses tried in turn.
#' @return The path of the land-polygon shapefile.
#' @references Wessel, P., Smith, W.H.F. (1996) A global, self-consistent,
#'   hierarchical, high-resolution shoreline database. *Journal of
#'   Geophysical Research* 101, 8741-8743.
#' @export
#' @examples
#' \dontrun{
#' land <- download_shoreline("data/coast")
#' lay <- filter_layers(land = land)
#' }
download_shoreline <- function(path, resolution = c("h", "f", "i", "l", "c"),
                               urls = c(
  "https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip",
  "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip")) {
  resolution <- match.arg(resolution)
  shp <- file.path(path, "GSHHS_shp", resolution,
                   paste0("GSHHS_", resolution, "_L1.shp"))
  if (file.exists(shp)) return(shp)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  zip <- file.path(path, "gshhg-shp-2.3.7.zip")
  ok <- file.exists(zip) && file.size(zip) > 1e6
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old), add = TRUE)
  for (u in urls) {
    if (ok) break
    message("Downloading the GSHHG shoreline (~150 MB, once) from ", u)
    ok <- tryCatch({
      utils::download.file(u, zip, mode = "wb", quiet = TRUE)
      file.exists(zip) && file.size(zip) > 1e6
    }, error = function(e) FALSE, warning = function(w) FALSE)
  }
  if (!ok) {
    stop("Download failed. Download gshhg-shp-2.3.7.zip by hand from one of ",
         paste(urls, collapse = " or "), ", put it in ", path,
         " and run this again.", call. = FALSE)
  }
  nm <- utils::unzip(zip, list = TRUE)$Name
  utils::unzip(zip, files = grep(paste0("GSHHS_shp/", resolution, "/"), nm,
                                 value = TRUE), exdir = path)
  if (!file.exists(shp)) stop("Unexpected archive content in ", zip,
                              call. = FALSE)
  shp
}

#' Aggregate GHS-BUILT-S to ~1 km settlement layers
#'
#' The JRC Global Human Settlement Layer GHS-BUILT-S (built-up surface, m2
#' per 100 m cell; e.g.
#' `GHS_BUILT_S_E2030_GLOBE_R2023A_54009_100_V1_0.tif`, ~2.8 GB, from
#' <https://ghsl.jrc.ec.europa.eu>) is too large to query directly. This
#' builds, once, two ~1 km (30 arc-second) layers in longitude/latitude:
#' the **mean** built-up m2 per 100 m cell (settlement density; divided by
#' 100 it is the built-up %) and the **maximum** (the most built 100 m cell,
#' i.e. a building nearby). Use them as `ghs_mean` and `ghs_max` in
#' [filter_layers()]. Needs the `sf` package (GDAL); takes several minutes.
#'
#' @param ghs_file Path to the GHS-BUILT-S GeoTIFF.
#' @param path Output folder.
#' @param extent Output extent `c(xmin, ymin, xmax, ymax)` in degrees;
#'   restrict it to your region to save time and space.
#' @param res Output resolution in degrees (default 30 arc-seconds).
#' @return A named character vector with the paths of the `mean` and `max`
#'   layers.
#' @export
#' @examples
#' \dontrun{
#' ghs <- aggregate_ghs("GHS_BUILT_S_E2030_GLOBE_R2023A_54009_100_V1_0.tif",
#'                      "data", extent = c(-25, 35, 45, 72))
#' lay <- filter_layers(ghs_mean = ghs[["mean"]], ghs_max = ghs[["max"]])
#' }
aggregate_ghs <- function(ghs_file, path, extent = c(-180, -60, 180, 85),
                          res = 1 / 120) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Install the sf package: install.packages(\"sf\")", call. = FALSE)
  }
  if (!file.exists(ghs_file)) stop("File not found: ", ghs_file, call. = FALSE)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  out <- c(mean = file.path(path, "ghs_built_mean_30as.tif"),
           max = file.path(path, "ghs_built_max_30as.tif"))
  for (k in names(out)) {
    if (file.exists(out[[k]])) next
    message("Building ", basename(out[[k]]), " (once; may take minutes) ...")
    tryCatch(
      sf::gdal_utils("warp", source = ghs_file, destination = out[[k]],
        options = c("-overwrite", "-t_srs", "EPSG:4326",
                    "-te", as.character(extent),
                    "-tr", as.character(res), as.character(res),
                    "-r", if (k == "mean") "average" else "max",
                    "-ovr", "NONE", "-ot", "Float32",
                    "-co", "COMPRESS=LZW", "-co", "TILED=YES",
                    "-multi", "-wo", "NUM_THREADS=ALL_CPUS")),
      error = function(e) {
        unlink(out[[k]])
        stop("GDAL warp failed: ", conditionMessage(e), call. = FALSE)
      })
  }
  out
}
