# Soil predictors from SoilGrids

# Layers of SoilGrids 2.0: conversion factors (mapped value / factor =
# conventional unit) from the SoilGrids FAQ (internal)
soilgrids_props <- function() {
  data.frame(
    var = c("phh2o", "soc", "nitrogen", "cec", "clay", "sand", "silt",
            "bdod", "cfvo", "ocd"),
    factor = c(10, 10, 100, 10, 10, 10, 10, 100, 10, 10),
    unit = c("pH", "g/kg", "g/kg", "cmol(c)/kg", "%", "%", "%", "kg/dm3",
             "vol%", "kg/m3"),
    stringsAsFactors = FALSE)
}

soilgrids_depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm",
                      "100-200cm")

#' Download SoilGrids soil layers for a region
#'
#' Fetches soil properties from SoilGrids 2.0 (Poggio et al. 2021; about
#' 250 m) for one or more depth intervals, plus the depth to bedrock from
#' SoilGrids 2017 (Shangguan et al. 2017), for a region only: as with
#' [download_chelsa_monthly()], only the cells of `extent` are read from
#' the global files over the internet (GDAL), and each layer is saved as a
#' GeoTIFF in `path`. Values are converted to conventional units (pH,
#' g/kg, %, cmol(c)/kg) with the factors of the SoilGrids documentation,
#' and the carbon to nitrogen ratio is added. Existing files are kept, so
#' an interrupted run continues where it stopped.
#'
#' Soil matters for plants and is missing from climate grids: in-situ soil
#' chemistry (B-horizon C:N, exchangeable calcium) predicted forest
#' understorey plants as well as topography did (Roe et al. 2022). Mapped
#' soil properties are a weaker stand-in for such measurements: SoilGrids is
#' itself a machine-learning prediction from soil profiles and covariates,
#' with large uncertainty where profiles are few (use `stat = "Q0.05"` and
#' `"Q0.95"`, or `"uncertainty"`, to see it), and exchangeable cations are
#' not mapped. pH is a poor indicator of fertility in acidic soils (Roe et
#' al. 2022); C:N, CEC and texture add to it.
#'
#' **Layers** (`vars`): `phh2o` pH in water; `soc` organic carbon (g/kg);
#' `nitrogen` total nitrogen (g/kg); `cec` cation exchange capacity at pH 7
#' (cmol(c)/kg); `clay`, `sand`, `silt` (%); `bdod` bulk density (kg/dm3);
#' `cfvo` coarse fragments (vol%); `ocd` organic carbon density (kg/m3).
#' `cn` is computed as `soc / nitrogen` when both are downloaded. `bedrock`
#' adds the absolute depth to bedrock (cm; SoilGrids 2017 `BDTICM`), a
#' single layer: in SoilGrids 2.0 it is not mapped.
#'
#' **Depths.** SoilGrids 2.0 maps six intervals: 0-5, 5-15, 15-30, 30-60,
#' 60-100 and 100-200 cm. The top one or two are the usual choice for
#' plants with shallow roots; deeper intervals follow trees and
#' deep-rooted shrubs.
#'
#' **Grid.** By default the layers are warped from the Homolosine
#' projection of SoilGrids to longitude/latitude at 7.5 arc-seconds (about
#' 250 m). With `template` (e.g. one layer of your climate stack) they are
#' warped directly onto its grid, averaging the 250 m cells, so they stack
#' with the climate layers without further resampling.
#'
#' **Availability.** The server (files.isric.org) is checked first; if it
#' cannot be reached the function stops with an explanation, and a layer
#' that fails half-way gives a warning and is fetched again on the next
#' run.
#'
#' @param path Folder for the files (created if needed).
#' @param extent Region: `c(xmin, xmax, ymin, ymax)` in degrees, a
#'   `SpatExtent`, or an `sf` / `SpatVector` / `SpatRaster` (its extent).
#'   Ignored when `template` is given.
#' @param vars SoilGrids 2.0 properties (see Details).
#' @param depth Depth interval(s), e.g. `"0-5cm"` (default) or
#'   `c("0-5cm", "5-15cm")`.
#' @param stat `"mean"` (default), `"Q0.5"` (median), `"Q0.05"`, `"Q0.95"`
#'   or `"uncertainty"` (the 90% prediction interval divided by the
#'   median, x 10).
#' @param bedrock Add the depth to bedrock (SoilGrids 2017)?
#' @param cn Add the carbon to nitrogen ratio (needs `soc` and
#'   `nitrogen`)?
#' @param template Optional raster whose grid (CRS, extent, resolution) the
#'   layers are warped onto, by averaging.
#' @param resolution Cell size in degrees without `template`.
#' @param prefix File name prefix.
#' @param overwrite Download again existing files?
#' @param base_url,bedrock_url Locations of the SoilGrids 2.0 folders and of
#'   the 2017 depth-to-bedrock file (change only for a mirror or a local
#'   copy).
#' @return A `SpatRaster` of the downloaded layers, named
#'   `<var>_<depth>` (e.g. `phh2o_0_5cm`, `cn_0_5cm`, `bedrock_depth`),
#'   invisibly.
#' @seealso `climate_index_names("soil", table = TRUE)` for the layers,
#'   their units and sources; [crop_rasters()] to cut other layers to the
#'   same region.
#' @references Poggio, L., de Sousa, L.M., Batjes, N.H., et al. (2021)
#'   SoilGrids 2.0: producing soil information for the globe with
#'   quantified spatial uncertainty. *SOIL* 7, 217-240.
#'   \doi{10.5194/soil-7-217-2021}
#'
#'   Roe, N.A., Ducey, M.J., Lee, T.D., et al. (2022) Soil chemical
#'   variables improve models of understorey plant species distributions.
#'   *Journal of Biogeography* 49, 753-766. \doi{10.1111/jbi.14344}
#'
#'   Shangguan, W., Hengl, T., Mendes de Jesus, J., Yuan, H., Dai, Y. (2017)
#'   Mapping the global depth to bedrock for land surface modeling.
#'   *Journal of Advances in Modeling Earth Systems* 9, 65-88.
#'   \doi{10.1002/2016MS000686}
#' @examples
#' \donttest{
#' # pH, organic carbon, nitrogen and C:N of the top 5 cm, and bedrock
#' if (interactive()) {
#'   soil <- download_soilgrids("soil", extent = c(19.8, 20.2, 49.2, 49.4),
#'                              vars = c("phh2o", "soc", "nitrogen"))
#'   soil
#' }
#' }
#' @export
download_soilgrids <- function(path, extent = NULL,
                               vars = c("phh2o", "soc", "nitrogen", "cec",
                                        "clay", "sand", "silt"),
                               depth = "0-5cm",
                               stat = c("mean", "Q0.5", "Q0.05", "Q0.95",
                                        "uncertainty"),
                               bedrock = TRUE, cn = TRUE, template = NULL,
                               resolution = 1 / 480, prefix = "",
                               overwrite = FALSE,
                               base_url = "https://files.isric.org/soilgrids/latest/data",
                               bedrock_url = "https://files.isric.org/soilgrids/former/2017-03-10/data/BDTICM_M_250m_ll.tif") {
  stat <- match.arg(stat)
  props <- soilgrids_props()
  bad <- setdiff(vars, props$var)
  if (length(bad)) {
    stop("Unknown SoilGrids layer: ", paste(bad, collapse = ", "),
         ". Choose from ", paste(props$var, collapse = ", "), ".",
         call. = FALSE)
  }
  bad <- setdiff(depth, soilgrids_depths)
  if (length(bad) || !length(depth)) {
    stop("Unknown depth: ", paste(bad, collapse = ", "), ". SoilGrids 2.0 ",
         "maps ", paste(soilgrids_depths, collapse = ", "), ".",
         call. = FALSE)
  }
  if (is.null(template) && is.null(extent)) {
    stop("Give `extent` (or `template`): the global layers are far too ",
         "large to download whole.", call. = FALSE)
  }
  if (!is.null(template)) {
    if (is.character(template)) template <- terra::rast(template)
    if (!inherits(template, "SpatRaster")) {
      stop("`template` must be a SpatRaster or a file path.", call. = FALSE)
    }
    template <- template[[1]]
    win <- NULL
  } else {
    win <- lonlat_window(extent)
    n <- prod(diff(win)[c(1, 3)] / resolution)
    if (n > 2e8) {
      warning("The region has about ", format(round(n), big.mark = ","),
              " cells per layer at this resolution: the download will be ",
              "slow and the files large. Consider a coarser `resolution` or ",
              "a `template`.", call. = FALSE)
    }
  }
  remote <- grepl("^https?://", base_url)
  if (remote) soilgrids_reachable(base_url)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  tag <- function(d) gsub("-", "_", d)
  stat_tag <- if (stat == "mean") "" else paste0("_", gsub("\\.", "", stat))
  files <- character(0)
  for (d in depth) {
    for (v in vars) {
      name <- paste0(v, "_", tag(d), stat_tag)
      out <- file.path(path, paste0(prefix, name, ".tif"))
      src <- paste0(base_url, "/", v, "/", v, "_", d, "_", stat, ".vrt")
      f <- soilgrids_get(src, out, win, template, resolution, remote,
                         factor = if (stat == "uncertainty") 1 else
                           props$factor[props$var == v],
                         name = name, overwrite = overwrite)
      files <- c(files, f)
    }
    if (cn && all(c("soc", "nitrogen") %in% vars) && stat != "uncertainty") {
      name <- paste0("cn_", tag(d), stat_tag)
      out <- file.path(path, paste0(prefix, name, ".tif"))
      fs <- file.path(path, paste0(prefix, c("soc", "nitrogen"), "_", tag(d),
                                   stat_tag, ".tif"))
      if (all(file.exists(fs)) && (overwrite || !file.exists(out))) {
        s <- terra::rast(fs[1])
        nn <- terra::rast(fs[2])
        r <- terra::ifel(nn > 0, s / nn, NA)
        names(r) <- name
        terra::writeRaster(r, out, overwrite = TRUE, names = name,
                           gdal = "COMPRESS=DEFLATE")
      }
      if (file.exists(out)) files <- c(files, out)
    }
  }
  if (bedrock) {
    name <- "bedrock_depth"
    out <- file.path(path, paste0(prefix, name, ".tif"))
    files <- c(files, soilgrids_get(bedrock_url, out, win, template,
                                    resolution,
                                    remote = grepl("^https?://", bedrock_url),
                                    factor = 1, name = name,
                                    overwrite = overwrite))
  }
  files <- files[file.exists(files)]
  if (!length(files)) {
    warning("No SoilGrids layer could be downloaded; run again to retry.",
            call. = FALSE)
    return(invisible(NULL))
  }
  invisible(terra::rast(files))
}

# Is the SoilGrids server reachable? (internal)
soilgrids_reachable <- function(base_url) {
  ok <- tryCatch({
    con <- url(paste0(base_url, "/"), open = "rb")
    close(con)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    stop("SoilGrids (", sub("^https?://([^/]+).*", "\\1", base_url),
         ") cannot be reached now: no internet connection, or the service ",
         "is down or under maintenance (see https://www.isric.org). Try ",
         "again later; files already downloaded are kept.", call. = FALSE)
  }
  invisible(TRUE)
}

# Warp one layer onto the target grid, convert units, save (internal)
soilgrids_get <- function(src, out, win, template, resolution, remote,
                          factor, name, overwrite) {
  if (file.exists(out) && !overwrite) return(out)
  message("Downloading ", basename(out))
  tmp <- paste0(out, ".part.tif")
  on.exit(unlink(tmp), add = TRUE)
  if (is.null(template)) {
    grid <- c("-t_srs", "EPSG:4326",
              "-te", num_arg(win[c(1, 3, 2, 4)]),
              "-tr", num_arg(c(resolution, resolution)),
              "-r", "bilinear")
  } else {
    e <- as.vector(terra::ext(template))
    tf <- tempfile(fileext = ".wkt")
    on.exit(unlink(tf), add = TRUE)
    writeLines(terra::crs(template), tf)
    grid <- c("-t_srs", tf,
              "-te", num_arg(e[c(1, 3, 2, 4)]),
              "-ts", ncol(template), nrow(template),
              "-r", "average")
  }
  ok <- tryCatch({
    sf::gdal_utils("warp", if (remote) paste0("/vsicurl/", src) else src,
                   tmp,
                   options = c(grid, "-ot", "Float32", "-co",
                               "COMPRESS=DEFLATE", "-co", "TILED=YES"),
                   config_options = c(GDAL_HTTP_MAX_RETRY = "5",
                                      GDAL_HTTP_RETRY_DELAY = "5"),
                   quiet = TRUE)
    TRUE
  }, error = function(e) {
    warning("Failed: ", basename(out), " (", conditionMessage(e),
            "); run again to retry.", call. = FALSE)
    FALSE
  })
  if (!ok || !file.exists(tmp)) return(character(0))
  r <- terra::rast(tmp) / factor
  names(r) <- name
  terra::writeRaster(r, out, overwrite = TRUE, names = name,
                     gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  out
}
