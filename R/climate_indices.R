# Derived climate variables from monthly climatologies ---------------------------

chelsa_base <- "https://os.zhdk.cloud.switch.ch/chelsav2"

#' URLs of CHELSA V2.1 monthly climatologies
#'
#' Finds the files of CHELSA V2.1 monthly climatologies (12 months per
#' variable) for the present (1981-2010) or a future period, climate model
#' (GCM) and scenario (SSP), by listing the CHELSA storage, so the exact file
#' names (which differ between present and future data) need not be known.
#'
#' @param vars Variables: `"tasmax"`, `"tasmin"`, `"tas"` (mean temperature)
#'   and `"pr"` (precipitation); any other CHELSA monthly variable, e.g.
#'   `"rsds"` or `"vpd"`, works too.
#' @param period `"1981-2010"` (present) or a future period, `"2011-2040"`,
#'   `"2041-2070"` or `"2071-2100"`.
#' @param gcm,ssp For future periods: the climate model, e.g.
#'   `"IPSL-CM6A-LR"`, `"MPI-ESM1-2-HR"`, `"UKESM1-0-LL"`, `"GFDL-ESM4"`,
#'   `"MRI-ESM2-0"`, and the scenario, `"ssp126"`, `"ssp370"` or `"ssp585"`.
#' @param months Months to return.
#' @return A data.frame with `var`, `month`, `url` and `file` (the file name).
#' @export
#' @examples
#' \dontrun{
#' chelsa_monthly_urls(c("tasmax", "pr"))
#' chelsa_monthly_urls("tas", "2071-2100", gcm = "MPI-ESM1-2-HR",
#'                     ssp = "ssp370")
#' }
chelsa_monthly_urls <- function(vars = c("tasmax", "tasmin", "tas", "pr"),
                                period = "1981-2010", gcm = NULL,
                                ssp = NULL, months = 1:12) {
  present <- period == "1981-2010"
  if (!present && (is.null(gcm) || is.null(ssp))) {
    stop("Future periods need `gcm` and `ssp`.", call. = FALSE)
  }
  rows <- lapply(vars, function(v) {
    dir <- if (present) {
      sprintf("GLOBAL/climatologies/%s/%s/", period, v)
    } else {
      sprintf("GLOBAL/climatologies/%s/%s/%s/%s/", period, toupper(gcm),
              tolower(ssp), v)
    }
    keys <- chelsa_list(dir)
    mm <- sprintf("%02d", months)
    hit <- vapply(mm, function(m) {
      k <- keys[grepl(paste0("_", v, "_", m, "_"), basename(keys))]
      if (length(k)) k[1] else NA_character_
    }, character(1))
    if (anyNA(hit)) {
      stop("Not found on the CHELSA server: ", v, " ", period,
           if (!present) paste0(" ", gcm, " ", ssp), ", month(s) ",
           paste(months[is.na(hit)], collapse = ", "),
           ". Check the names (e.g. gcm = \"MPI-ESM1-2-HR\", ssp = \"ssp370\").",
           call. = FALSE)
    }
    data.frame(var = v, month = months, url = paste0(chelsa_base, "/", hit),
               file = basename(hit), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# List the CHELSA storage under a prefix (S3 listing) (internal)
chelsa_list <- function(prefix) {
  u <- paste0(chelsa_base, "/?prefix=", prefix, "&max-keys=1000")
  txt <- tryCatch(paste(readLines(u, warn = FALSE), collapse = ""),
                  error = function(e) {
                    stop("Cannot reach the CHELSA server (", conditionMessage(e),
                         "). Check your internet connection or download the ",
                         "files by hand from chelsa-climate.org.",
                         call. = FALSE)
                  })
  keys <- regmatches(txt, gregexpr("<Key>[^<]+</Key>", txt))[[1]]
  gsub("</?Key>", "", keys)
}

#' Download CHELSA V2.1 monthly climatologies, optionally for a region only
#'
#' Downloads the 12 monthly files of each variable needed by
#' [climate_indices()]. With `extent`, only that region is downloaded:
#' GDAL reads the needed part of each global file over the internet and
#' writes it locally, keeping the stored values, scale and offset exactly
#' (no global 100+ MB files on disk). Files already present are skipped, so
#' an interrupted download can simply be run again.
#'
#' The global files are large (about 110 MB each, 48 files per period and
#' scenario for four variables); a region needs only its share, but at
#' least the full width of the rows it covers is read. Data:
#' Karger et al. (2017, 2023), CHELSA V2.1, <https://chelsa-climate.org>;
#' please cite them.
#'
#' @inheritParams chelsa_monthly_urls
#' @param path Folder for the files (created if needed).
#' @param extent Optional region: `c(xmin, xmax, ymin, ymax)` in degrees, a
#'   [terra::SpatExtent], or any object with an extent and CRS (`sf`,
#'   [terra::SpatVector], [terra::SpatRaster]), which is converted to
#'   longitude/latitude. `NULL` downloads the global files.
#' @param overwrite Download again files that exist.
#' @return The paths of the files, invisibly.
#' @export
#' @examples
#' \dontrun{
#' # present climate of the Alps, four variables x 12 months
#' f <- download_chelsa_monthly("chelsa_alps", extent = c(5, 17, 43, 49))
#' # a future scenario
#' download_chelsa_monthly("chelsa_alps", period = "2071-2100",
#'                         gcm = "MPI-ESM1-2-HR", ssp = "ssp370",
#'                         extent = c(5, 17, 43, 49))
#' }
download_chelsa_monthly <- function(path,
                                    vars = c("tasmax", "tasmin", "tas", "pr"),
                                    period = "1981-2010", gcm = NULL,
                                    ssp = NULL, extent = NULL,
                                    overwrite = FALSE) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  u <- chelsa_monthly_urls(vars, period, gcm, ssp)
  win <- if (is.null(extent)) NULL else lonlat_window(extent)
  out <- file.path(path, u$file)
  for (i in seq_len(nrow(u))) {
    if (file.exists(out[i]) && !overwrite) next
    message("Downloading ", u$file[i], " (", i, "/", nrow(u), ")")
    tmp <- paste0(out[i], ".part")
    ok <- tryCatch({
      if (is.null(win)) {
        utils::download.file(u$url[i], tmp, mode = "wb", quiet = TRUE)
      } else {
        sf::gdal_utils("translate", paste0("/vsicurl/", u$url[i]), tmp,
                       options = c("-of", "GTiff",
                                   "-projwin", num_arg(win[c(1, 4, 2, 3)]),
                                   "-co", "COMPRESS=DEFLATE",
                                   "-co", "TILED=YES"),
                       config_options = c(GDAL_HTTP_MAX_RETRY = "5",
                                          GDAL_HTTP_RETRY_DELAY = "5"),
                       quiet = TRUE)
      }
      TRUE
    }, error = function(e) {
      warning("Failed: ", u$file[i], " (", conditionMessage(e),
              "); run again to retry.", call. = FALSE)
      FALSE
    })
    if (ok) file.rename(tmp, out[i]) else unlink(tmp)
  }
  invisible(out[file.exists(out)])
}

# Extent as lon/lat c(xmin, xmax, ymin, ymax) (internal)
lonlat_window <- function(extent) {
  if (is.numeric(extent)) {
    e <- unname(extent)
    if (length(e) != 4 || e[1] >= e[2] || e[3] >= e[4]) {
      stop("`extent` must be c(xmin, xmax, ymin, ymax).", call. = FALSE)
    }
    return(e)
  }
  if (inherits(extent, "SpatExtent")) return(as.vector(extent))
  if (inherits(extent, c("sf", "sfc"))) extent <- terra::vect(extent)
  if (inherits(extent, c("SpatVector", "SpatRaster"))) {
    if (!terra::is.lonlat(extent)) {
      extent <- terra::project(terra::as.polygons(terra::ext(extent),
                                                  crs = terra::crs(extent)),
                               "EPSG:4326")
    }
    return(as.vector(terra::ext(extent)))
  }
  stop("Unsupported `extent`.", call. = FALSE)
}

#' Find the 12 monthly files of a variable
#'
#' Finds and orders the monthly files of one variable in a folder, e.g. the
#' files written by [download_chelsa_monthly()] or downloaded by hand. A file
#' belongs to variable `var` and month `MM` when its name contains
#' `_<var>_<MM>_` (so `tas` does not match `tasmax`).
#'
#' @param dir Folder with the files.
#' @param var Variable, e.g. `"tasmax"`.
#' @param must_contain Optional text(s) the file name must also contain
#'   (case-insensitive), e.g. the GCM and scenario: `c("mpi-esm1-2-hr",
#'   "ssp370")`, or `"1981-2010"`.
#' @return 12 file paths, January to December.
#' @export
#' @examples
#' d <- tempfile()
#' dir.create(d)
#' for (m in 1:12) file.create(file.path(d, sprintf("CHELSA_tas_%02d_1981-2010_V.2.1.tif", m)))
#' basename(monthly_files(d, "tas"))[1:2]
monthly_files <- function(dir, var, must_contain = NULL) {
  f <- list.files(dir, pattern = "\\.tif{1,2}$", full.names = TRUE,
                  ignore.case = TRUE)
  b <- tolower(basename(f))
  for (s in must_contain) {
    keep <- grepl(tolower(s), b, fixed = TRUE)
    f <- f[keep]
    b <- b[keep]
  }
  out <- vapply(sprintf("%02d", 1:12), function(m) {
    hit <- f[grepl(paste0("_", tolower(var), "_", m, "_"), b, fixed = TRUE)]
    if (length(hit) > 1) {
      stop("Several files for ", var, " month ", m, " in ", dir,
           ": use `must_contain` to choose, e.g. the GCM and scenario.",
           call. = FALSE)
    }
    if (length(hit)) hit else NA_character_
  }, character(1))
  if (anyNA(out)) {
    stop("Missing ", var, " file(s) for month(s) ",
         paste(which(is.na(out)), collapse = ", "), " in ", dir,
         if (length(must_contain)) paste0(" (containing ",
                                          paste(must_contain, collapse = ", "),
                                          ")"), ".", call. = FALSE)
  }
  unname(out)
}

#' Derived climate variables from monthly climatologies
#'
#' Computes ecologically meaningful variables that standard bioclimatic
#' variables lack, from the 12 monthly means of maximum, minimum and mean
#' temperature and monthly precipitation, e.g. CHELSA or WorldClim
#' climatologies (see [download_chelsa_monthly()]). They were chosen for
#' plants with different photosynthetic pathways (C3/C4: water deficit,
#' evaporative demand, heat, and whether rain falls in the warm season), but
#' are useful for any water- or heat-limited species.
#'
#' | Index | Meaning | Unit |
#' |---|---|---|
#' | `vpd` | mean vapour pressure deficit (FAO-56; dew point from Tmin) | kPa |
#' | `pet` | potential evapotranspiration, Hargreaves-Samani (FAO-56) | mm/year |
#' | `cmd` | climatic moisture deficit: sum of monthly PET - P where positive | mm/year |
#' | `aridity` | aridity index P / PET (UNEP; < 0.05 hyper-arid, > 0.65 humid) | - |
#' | `heat_stress` | sum over months of Tmax above `heat_threshold` | degC x months |
#' | `sri` | summer rain index: precipitation of the warmest quarter / annual | 0-1 |
#' | `synchronicity` | heat-rain synchrony: 1 = the warmest and wettest month coincide, 0 = six months apart | 0-1 |
#'
#' **Formulas.** Saturation vapour pressure
#' \eqn{e_s(T) = 0.6108 \exp(17.27 T / (T + 237.3))} (kPa); monthly VPD
#' \eqn{= (e_s(T_{max}) + e_s(T_{min}))/2 - e_s(T_{min} - d)} with
#' `tdew_offset` \eqn{d} (FAO-56 suggests 2-3 degC in arid climates), and the
#' annual mean of the 12 months. PET per month
#' \eqn{= 0.0023 \times 0.408 R_a (T_{mean} + 17.8) \sqrt{T_{max} - T_{min}}}
#' (mm/day, with extraterrestrial radiation \eqn{R_a} in MJ m-2 day-1 from
#' latitude and the mid-month day, FAO-56 eq. 21; 0.408 converts it to mm
#' of evaporation) times the days of the month, never below 0. The quarters
#' for `sri` are any three consecutive months, also across the year's end
#' (e.g. December-February), the warmest by mean temperature, as for
#' bioclim variables: `sri` equals bio18 / bio12, so do not use it together
#' with both. Ties in `synchronicity` take the first month.
#'
#' **Units.** Temperatures must be in degC and precipitation in mm per
#' month. CHELSA V2.1 stores temperature as K x 10 with a scale (0.1) and
#' offset (-273.15) in the files, which terra applies when reading, so the
#' values are already in degC: do **not** divide them by 10.
#' `temp_units = "auto"` checks the values and converts Kelvin (or K x 10
#' read without the scale) and stops if they look divided by 10 by mistake.
#'
#' Run it on files already cut to your region (e.g. with [crop_rasters()],
#' which keeps CHELSA's values exact), in their own longitude/latitude grid,
#' and project the results afterwards if needed: projecting 7 results is
#' faster than projecting 48 monthly layers, and the indices are computed
#' on unsmoothed values. Any CRS works; latitude is computed for each cell.
#'
#' @param tmax,tmin,prec Monthly maximum and minimum temperature and
#'   precipitation: 12-layer [terra::SpatRaster]s (January to December) or
#'   12 file paths (e.g. from [monthly_files()]).
#' @param tmean Monthly mean temperature; default `(tmax + tmin) / 2`.
#' @param indices Which indices to compute (default: all).
#' @param heat_threshold Temperature (degC) above which `heat_stress` sums.
#' @param tdew_offset Dew point below Tmin (degC) for `vpd`.
#' @param temp_units `"auto"`, `"C"`, `"K"` or `"K10"` (Kelvin x 10).
#' @param path Optional folder: each index is written to
#'   `<prefix><index>.tif` (recommended for large areas).
#' @param prefix File name prefix, e.g. `"CHELSA_MPI-ESM1-2-HR_ssp370_"`.
#' @param overwrite Overwrite existing files; otherwise existing files are
#'   kept and read (so interrupted runs continue).
#' @return A [terra::SpatRaster] with one layer per index.
#' @references Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998) Crop
#'   evapotranspiration. FAO Irrigation and Drainage Paper 56. Hargreaves,
#'   G.H., Samani, Z.A. (1985) Reference crop evapotranspiration from
#'   temperature. *Applied Engineering in Agriculture* 1, 96-99. Karger, D.N.
#'   et al. (2017) Climatologies at high resolution for the earth's land
#'   surface areas. *Scientific Data* 4, 170122.
#' @export
#' @examples
#' # a small synthetic grid: warm wet summers in the south, dry in the north
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 30,
#'                  ymax = 40, nlyrs = 12)
#' season <- sin((1:12 - 4) / 12 * 2 * pi)       # warmest in July
#' lat <- terra::init(r[[1]], "y")
#' tmax <- terra::rast(lapply(1:12, function(m) 45 - lat + 10 * season[m]))
#' tmin <- tmax - 12
#' prec <- terra::rast(lapply(1:12, function(m) 40 + (40 - lat) * 8 *
#'                                pmax(season[m], 0)))
#' ci <- climate_indices(tmax, tmin, prec)
#' ci
#' terra::plot(ci)
climate_indices <- function(tmax, tmin, prec, tmean = NULL,
                            indices = c("vpd", "pet", "cmd", "aridity",
                                        "heat_stress", "sri",
                                        "synchronicity"),
                            heat_threshold = 30, tdew_offset = 0,
                            temp_units = c("auto", "C", "K", "K10"),
                            path = NULL, prefix = "", overwrite = FALSE) {
  temp_units <- match.arg(temp_units)
  all_idx <- c("vpd", "pet", "cmd", "aridity", "heat_stress", "sri",
               "synchronicity")
  bad <- setdiff(indices, all_idx)
  if (length(bad)) {
    stop("Unknown index: ", paste(bad, collapse = ", "), ". Available: ",
         paste(all_idx, collapse = ", "), call. = FALSE)
  }
  monthly <- function(x, what) {
    if (is.character(x)) x <- terra::rast(x)
    if (!inherits(x, "SpatRaster") || terra::nlyr(x) != 12) {
      stop("`", what, "` must be 12 monthly layers (January to December).",
           call. = FALSE)
    }
    x
  }
  tmax <- monthly(tmax, "tmax")
  tmin <- monthly(tmin, "tmin")
  prec <- monthly(prec, "prec")
  if (!is.null(tmean)) tmean <- monthly(tmean, "tmean")
  if (!terra::compareGeom(tmax, tmin, stopOnError = FALSE) ||
      !terra::compareGeom(tmax, prec, stopOnError = FALSE)) {
    stop("tmax, tmin and prec must be on the same grid.", call. = FALSE)
  }
  conv <- temperature_converter(tmax, temp_units)
  tmax <- conv(tmax)
  tmin <- conv(tmin)
  tmean <- if (is.null(tmean)) (tmax + tmin) / 2 else conv(tmean)
  check_precipitation(prec)

  if (!is.null(path) && !dir.exists(path)) dir.create(path, recursive = TRUE)
  out <- list()
  # compute (or read an existing file of) one index
  make <- function(name, expr) {
    f <- if (is.null(path)) NULL else
      file.path(path, paste0(prefix, name, ".tif"))
    if (!is.null(f) && file.exists(f) && !overwrite) {
      r <- terra::rast(f)
    } else {
      r <- force(expr)
      names(r) <- name
      if (!is.null(f)) {
        r <- terra::writeRaster(r, f, overwrite = TRUE,
                                gdal = c("COMPRESS=DEFLATE"))
      }
    }
    names(r) <- name
    out[[name]] <<- r
    r
  }
  days <- c(31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

  if ("vpd" %in% indices) {
    es <- function(t) 0.6108 * exp(17.27 * t / (t + 237.3))
    make("vpd", terra::mean((es(tmax) + es(tmin)) / 2 -
                              es(tmin - tdew_offset)))
  }
  need_pet <- any(c("pet", "cmd", "aridity") %in% indices)
  if (need_pet) {
    phi <- latitude_radians(tmax[[1]])
    pet_m <- terra::rast(lapply(1:12, function(i) {
      ra <- extraterrestrial_radiation(phi, i)
      p <- 0.0023 * 0.408 * ra * (tmean[[i]] + 17.8) *
        sqrt(terra::clamp(tmax[[i]] - tmin[[i]], lower = 0, values = TRUE))
      terra::clamp(p * days[i], lower = 0, values = TRUE)
    }))
    pet_ann <- sum(pet_m)
    if ("pet" %in% indices) make("pet", pet_ann)
    if ("cmd" %in% indices) {
      make("cmd", sum(terra::clamp(pet_m - prec, lower = 0, values = TRUE)))
    }
    if ("aridity" %in% indices) {
      # PET below 1 mm/year (polar night) is set to 1 mm to stay finite
      make("aridity", sum(prec) / terra::clamp(pet_ann, lower = 1,
                                                values = TRUE))
    }
  }
  if ("heat_stress" %in% indices) {
    make("heat_stress", sum(terra::clamp(tmax - heat_threshold, lower = 0,
                                         values = TRUE)))
  }
  if ("sri" %in% indices) {
    q <- function(m) ((m - 1 + 0:2) %% 12) + 1
    tq <- terra::rast(lapply(1:12, function(m) terra::mean(tmean[[q(m)]])))
    pq <- terra::rast(lapply(1:12, function(m) sum(prec[[q(m)]])))
    warm <- terra::which.max(tq)
    p_warm <- terra::selectRange(pq, warm)
    tot <- sum(prec)
    make("sri", terra::ifel(tot > 0, p_warm / tot, NA))
  }
  if ("synchronicity" %in% indices) {
    d <- abs(terra::which.max(tmax) - terra::which.max(prec))
    d <- terra::ifel(d > 6, 12 - d, d)
    make("synchronicity", (6 - d) / 6)
  }
  terra::rast(out[intersect(indices, names(out))])
}

# Function converting temperatures to degC (internal)
temperature_converter <- function(x, units) {
  if (units == "auto") {
    v <- terra::spatSample(x[[c(1, 7)]], 2000, method = "regular",
                           na.rm = TRUE, as.df = FALSE)
    v <- as.vector(v)
    v <- v[is.finite(v)]
    if (!length(v)) stop("The temperature layers have no values.",
                         call. = FALSE)
    med <- stats::median(v)
    units <- if (med > 1500) "K10" else if (med > 150) "K" else "C"
    if (units == "C" && diff(stats::quantile(v, c(0.05, 0.95))) < 4 &&
        max(abs(v)) < 6) {
      stop("The temperatures look divided by 10 (they span ",
           paste(round(range(v), 2), collapse = " to "),
           " degC). CHELSA V2.1 values are already in degC when read with ",
           "terra; use them as they are, or set temp_units = \"C\" if ",
           "these values are right.", call. = FALSE)
    }
    if (units != "C") {
      message("Temperatures look like ", if (units == "K") "Kelvin" else
        "Kelvin x 10 (scale not applied)", "; converted to degC.")
    }
  }
  switch(units,
         C = function(r) r,
         K = function(r) r - 273.15,
         K10 = function(r) r * 0.1 - 273.15)
}

check_precipitation <- function(prec) {
  v <- as.vector(terra::spatSample(prec, 2000, method = "regular",
                                   na.rm = TRUE, as.df = FALSE))
  v <- v[is.finite(v)]
  if (length(v) && any(v < 0)) {
    stop("Negative precipitation values.", call. = FALSE)
  }
  if (length(v) && stats::quantile(v, 0.99) > 3000) {
    warning("Monthly precipitation above 3000 mm in 1% of cells: are the ",
            "values scaled (e.g. x 10 or x 100)? They must be in mm/month.",
            call. = FALSE)
  }
  invisible(TRUE)
}

# Latitude of every cell in radians, for any CRS (internal)
latitude_radians <- function(r) {
  if (terra::is.lonlat(r)) return(terra::init(r, "y") * pi / 180)
  xy <- c(terra::init(r, "x"), terra::init(r, "y"))
  wkt <- terra::crs(r)
  terra::app(xy, function(m) {
    out <- rep(NA_real_, nrow(m))
    ok <- stats::complete.cases(m)
    if (any(ok)) {
      out[ok] <- sf::sf_project(wkt, "EPSG:4326", m[ok, , drop = FALSE],
                                warn = FALSE)[, 2]
    }
    out * pi / 180
  })
}

# Extraterrestrial radiation (MJ m-2 day-1) for the middle of month i,
# FAO-56 equations 21-25 (internal)
extraterrestrial_radiation <- function(phi, i, j = NULL) {
  if (is.null(j)) {
    j <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)[i]
  }
  delta <- 0.409 * sin(2 * pi / 365 * j - 1.39)
  dr <- 1 + 0.033 * cos(2 * pi / 365 * j)
  arg <- -tan(phi) * tan(delta)
  arg <- if (inherits(arg, "SpatRaster")) {
    terra::clamp(arg, lower = -1, upper = 1, values = TRUE)
  } else pmin(pmax(arg, -1), 1)
  ws <- acos(arg)
  (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(phi) * sin(delta) + cos(phi) * cos(delta) * sin(ws))
}
