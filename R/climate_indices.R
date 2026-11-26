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
  download_rasters(u$url, file.path(path, u$file), extent, overwrite)
}

# Download raster files, whole or only a region of each (read over the
# internet by GDAL); existing files are kept (internal)
download_rasters <- function(url, out, extent = NULL, overwrite = FALSE) {
  win <- if (is.null(extent)) NULL else lonlat_window(extent)
  for (i in seq_along(url)) {
    if (file.exists(out[i]) && !overwrite) next
    message("Downloading ", basename(out[i]), " (", i, "/", length(url), ")")
    tmp <- paste0(out[i], ".part")
    ok <- tryCatch({
      if (is.null(win)) {
        utils::download.file(url[i], tmp, mode = "wb", quiet = TRUE)
      } else {
        sf::gdal_utils("translate", paste0("/vsicurl/", url[i]), tmp,
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
      warning("Failed: ", basename(out[i]), " (", conditionMessage(e),
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
#' `_<var>_<MM>_` (so `tas` does not match `tasmax`); the month may be
#' written `01` or `1`, as in CHELSA-TraCE21k.
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
    hit <- f[grepl(paste0("_", tolower(var), "_", m, "_"), b, fixed = TRUE) |
               grepl(paste0("_", tolower(var), "_", as.integer(m), "_"), b,
                     fixed = TRUE)]
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
#' | Index | Group | Meaning | Unit | Useful for | Source |
#' |---|---|---|---|---|---|
#' | `vpd` | `water` | mean vapour pressure deficit (dew point from Tmin) | kPa | dryness of the air: stomatal closure and transpiration of plants, desiccation of amphibians, snails and insects | Allen et al. 1998; Grossiord et al. 2020 |
#' | `pet` | `water` | potential evapotranspiration (Hargreaves-Samani) | mm/year | evaporative demand and available energy; plant productivity and species richness | Hargreaves & Samani 1985; Allen et al. 1998 |
#' | `cmd` | `water` | climatic moisture deficit: sum of monthly PET - P where positive | mm/year | drought not met by rain: forest-steppe limits, drought-sensitive trees, C3 versus C4 grasses | Stephenson 1990; Wang et al. 2012 |
#' | `aridity` | `water` | aridity index P / PET (< 0.05 hyper-arid, > 0.65 humid) | - | the humid-arid gradient: biome limits, dryland and desert taxa | UNEP 1992 |
#' | `heat_stress` | `heat` | sum over months of Tmax above `heat_threshold` | degC x months | exposure to hot months: upper thermal limits of alpine and boreal taxa and of ectotherms; C4 advantage | this package; Ehleringer et al. 1997 |
#' | `sri` | `season` | summer rain index: precipitation of the warmest quarter / annual (bio18 / bio12) | 0-1 | rain in the warm season: C4 versus C3 grasses, monsoon versus Mediterranean taxa | Teeri & Stowe 1976; Paruelo & Lauenroth 1996 |
#' | `synchronicity` | `season` | heat-rain synchrony: 1 = the warmest and the wettest month coincide, 0 = six months apart | 0-1 | whether water arrives when it is warm: summer-rain versus winter-rain floras, C3/C4 balance | this package; Paruelo & Lauenroth 1996 |
#' | `gdd` | `heat` | growing degree days: monthly Tmean above `gdd_base` x days | degC x days | heat sum for growth and development: cold range limits of plants, number of insect generations | Woodward 1987; Prentice et al. 1992 |
#' | `cold_stress` | `cold` | sum over months of Tmin below `frost_threshold` | degC x months | winter severity: frost-sensitive plants, survival of overwintering stages | Woodward 1987; Sakai & Larcher 1987 |
#' | `frost_months` | `cold` | number of months with Tmin below `frost_threshold` | months | duration of frost: evergreen versus deciduous plants, length of dormancy or hibernation | Woodward 1987; Sakai & Larcher 1987 |
#' | `continentality` | `heat` | Tmean of the warmest month - that of the coldest | degC | annual temperature amplitude: oceanic versus continental species | Rivas-Martinez et al. 2011 |
#' | `dry_months` | `water` | months with P < 2 x Tmean (the dry season of Walter-Lieth diagrams) | months | length of the dry season: Mediterranean and savanna taxa, drought-deciduous plants | Bagnouls & Gaussen 1953; Walter & Lieth 1960 |
#' | `lgp` | `season` | length of the growing period: months with Tmean >= `gdd_base` and P >= PET / 2 | months | time when both warmth and water allow growth: crops and their wild relatives, annual plants | FAO 1978 |
#' | `emberger` | `water` | Emberger's pluviothermic quotient 2000 P / (M^2 - m^2), M and m in Kelvin | - | Mediterranean bioclimatic stages (arid to humid); lower = drier | Emberger 1955; Daget 1977 |
#' | `ellenberg` | `water` | Ellenberg's quotient: Tmean of the warmest month / annual P x 1000 | degC/mm x 1000 | humidity gradient of temperate forests: beech dominates below about 30, oaks above | Ellenberg 1963; Fang & Lechowicz 2006 |
#' | `kira_warmth` | `heat` | Kira's warmth index: sum over months of Tmean - 5 where Tmean > 5 degC | degC x months | vegetation zones of East Asia (evergreen broad-leaved forest between 85 and 180) | Kira 1945; Kira 1991 |
#' | `kira_coldness` | `cold` | Kira's coldness index: sum over months of Tmean - 5 where Tmean < 5 degC (zero or negative) | degC x months | cold limit of evergreen broad-leaved trees in East Asia (about -10) | Kira 1945; Kira 1991 |
#' | `ra` | `light` | mean extraterrestrial radiation over the year (latitude only) | MJ m-2 day-1 | light at the top of the atmosphere; in a model a smooth stand-in for latitude | Allen et al. 1998 |
#' | `ra_gs` | `light` | extraterrestrial radiation received in the months with Tmean >= `gdd_base` | MJ m-2 | light available while it is warm enough to grow: productivity at high latitudes and elevations | this package; Allen et al. 1998 |
#' | `daylength_start` | `light` | day length in the middle of the month in which the growing season starts | hours | photoperiod at the onset of growth: bud burst, flowering of long-day plants, end of diapause | Koerner & Basler 2010; Way & Montgomery 2015 |
#' | `daylength_end` | `light` | day length in the middle of the month in which the growing season ends | hours | photoperiod at the end of the season: growth cessation and bud set of trees, start of diapause | Bradshaw & Holzapfel 2007; Way & Montgomery 2015 |
#' | `daylength_gs` | `light` | mean day length of the months with Tmean >= `gdd_base` | hours | photoperiod during growth: long-day versus short-day plants, daily activity time of diurnal animals | Saikkonen et al. 2012 |
#' | `daylight_gs` | `light` | hours of daylight in the months with Tmean >= `gdd_base` | hours | time for photosynthesis or foraging in the season: a short northern summer with long days against a long southern one | this package; Saikkonen et al. 2012 |
#' | `daylength_max` | `light` | day length on the summer solstice (latitude only; 24 = midnight sun) | hours | the longest day: species that need a critical day length to flower or to breed; like `ra` a stand-in for latitude | Allen et al. 1998; Bradshaw & Holzapfel 2007 |
#' | `days_to_solstice` | `light` | days from the start of the growing season to the summer solstice (negative = it starts later) | days | share of the season on lengthening days: phenological timing, species whose flowering is triggered before the solstice | this package; Koerner & Basler 2010 |
#'
#' All are computed by default; choose some by name or by group with
#' [climate_index_names()], e.g. `indices = climate_index_names("water")`
#' or `indices = c("cold", "gdd")`. "Useful for" gives examples of the
#' questions an index was made for or is commonly used in, not a limit on
#' its use; "Source" gives the origin of the index or of its ecological
#' use (see the references), and "this package" marks definitions that are
#' ours. Kira's indices always use 5 degC, whereas `gdd` uses `gdd_base`
#' and counts days; `continentality` uses mean, not extreme, temperatures,
#' unlike bio7. `ra` and `daylength_max` depend on latitude alone: in a
#' model they act as a smooth stand-in for latitude, so use them only with
#' that in mind. The other light indices combine latitude with the
#' temperature season, so they differ between a mountain and its foreland
#' and change with the climate of another period.
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
#' with both. Ties in `synchronicity` take the first month. In `emberger`,
#' M is the mean maximum temperature of the warmest month and m the mean
#' minimum of the coldest, P the annual precipitation. `gdd`, `lgp`,
#' `frost_months`, `cold_stress` and `dry_months` are computed from monthly
#' means, so they are approximations of the daily quantities: days above
#' the base within a month whose mean is below it (or frosts in a month
#' whose mean minimum is above the threshold) are not counted. The growing
#' season of `daylength_start` starts in the first month with Tmean at or
#' above `gdd_base`, counting from the coldest month (which is itself the
#' start where every month is warm enough); cells without such a month
#' get `NA`. Day length is \eqn{24 \omega_s / \pi} with the sunset hour
#' angle \eqn{\omega_s} of FAO-56 (eq. 25 and 34). The season ends in the
#' last month at or above `gdd_base` before the coldest month
#' (`daylength_end`). `daylength_gs` is the mean of the mid-month day
#' lengths of the months at or above `gdd_base`, weighted by their days
#' (`NA` without such a month), and `daylight_gs` their sum over those
#' days (0 without such a month). For `days_to_solstice` the start of the
#' season is placed between two mid-month days by linear interpolation of
#' the monthly mean temperatures (where every month is warm enough it is
#' the middle of the coldest month); the solstice is day 172 north and day
#' 355 south of the equator, and the result lies between -182 and 182.
#' With monthly means it is accurate to about one or two weeks. Several
#' indices are strongly correlated with each other and with bioclim
#' variables (e.g. `gdd` with bio1, `ellenberg` and `emberger` with
#' `aridity`); let [select_variables()] choose among them.
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
#' @param indices Which indices to compute: names of indices or of groups
#'   of them (see [climate_index_names()]), or `"all"` (default).
#' @param heat_threshold Temperature (degC) above which `heat_stress` sums.
#' @param gdd_base Base temperature (degC) of `gdd`, `lgp`, `ra_gs` and
#'   the day-length indices of the growing season: 5 suits
#'   temperate plants; use 0 for arctic-alpine plants or about 10 for
#'   thermophilous plants and many insects.
#' @param frost_threshold Temperature (degC) below which `cold_stress`
#'   sums and `frost_months` counts.
#' @param tdew_offset Dew point below Tmin (degC) for `vpd`.
#' @param temp_units `"auto"`, `"C"`, `"K"` or `"K10"` (Kelvin x 10).
#' @param path Optional folder: each index is written to
#'   `<prefix><index>.tif` (recommended for large areas).
#' @param prefix File name prefix, e.g. `"CHELSA_MPI-ESM1-2-HR_ssp370_"`.
#' @param overwrite Overwrite existing files; otherwise existing files are
#'   kept and read (so interrupted runs continue).
#' @return A [terra::SpatRaster] with one layer per index.
#' @references Allen, R.G., Pereira, L.S., Raes, D., Smith, M. (1998) Crop
#'   evapotranspiration. FAO Irrigation and Drainage Paper 56.
#'
#' Bagnouls, F., Gaussen, H. (1953) Saison seche et indice xerothermique.
#'   *Bulletin de la Societe d'Histoire Naturelle de Toulouse* 88, 193-239.
#'
#' Bradshaw, W.E., Holzapfel, C.M. (2007) Evolution of animal
#'   photoperiodism. *Annual Review of Ecology, Evolution, and Systematics*
#'   38, 1-25.
#'
#' Daget, P. (1977) Le bioclimat mediterraneen: analyse des formes
#'   climatiques par le systeme d'Emberger. *Vegetatio* 34, 87-103.
#'
#' Ehleringer, J.R., Cerling, T.E., Helliker, B.R. (1997) C4
#'   photosynthesis, atmospheric CO2, and climate. *Oecologia* 112, 285-299.
#'
#' Ellenberg, H. (1963) Vegetation Mitteleuropas mit den Alpen. Ulmer,
#'   Stuttgart.
#'
#' Emberger, L. (1955) A biogeographic classification of climates. *Recueil
#'   des Travaux des Laboratoires de Botanique, Geologie et Zoologie de la
#'   Faculte des Sciences de l'Universite de Montpellier, Serie Botanique* 7,
#'   3-43.
#'
#' Fang, J., Lechowicz, M.J. (2006) Climatic limits for the present
#'   distribution of beech (Fagus L.) species in the world. *Journal of
#'   Biogeography* 33, 1804-1819.
#'
#' FAO (1978) Report on the agro-ecological zones project. World Soil
#'   Resources Report 48.
#'
#' Grossiord, C. et al. (2020) Plant responses to rising vapor pressure
#'   deficit. *New Phytologist* 226, 1550-1566.
#'
#' Hargreaves, G.H., Samani, Z.A. (1985) Reference crop evapotranspiration
#'   from temperature. *Applied Engineering in Agriculture* 1, 96-99.
#'
#' Karger, D.N. et al. (2017) Climatologies at high resolution for the
#'   earth's land surface areas. *Scientific Data* 4, 170122.
#'
#' Kira, T. (1945) A new classification of climate in eastern Asia as the
#'   basis for agricultural geography. Horticultural Institute, Kyoto
#'   University.
#'
#' Kira, T. (1991) Forest ecosystems of East and Southeast Asia in a global
#'   perspective. *Ecological Research* 6, 185-200.
#'
#' Koerner, C., Basler, D. (2010) Phenology under global warming. *Science*
#'   327, 1461-1462.
#'
#' Paruelo, J.M., Lauenroth, W.K. (1996) Relative abundance of plant
#'   functional types in grasslands and shrublands of North America.
#'   *Ecological Applications* 6, 1212-1224.
#'
#' Prentice, I.C. et al. (1992) A global biome model based on plant
#'   physiology and dominance, soil properties and climate. *Journal of
#'   Biogeography* 19, 117-134.
#'
#' Rivas-Martinez, S., Rivas Saenz, S., Penas, A. (2011) Worldwide
#'   bioclimatic classification system. *Global Geobotany* 1, 1-634.
#'
#' Saikkonen, K. et al. (2012) Climate change-driven species' range shifts
#'   filtered by photoperiodism. *Nature Climate Change* 2, 239-242.
#'
#' Sakai, A., Larcher, W. (1987) Frost survival of plants. Springer,
#'   Berlin.
#'
#' Stephenson, N.L. (1990) Climatic control of vegetation distribution: the
#'   role of the water balance. *American Naturalist* 135, 649-670.
#'
#' Teeri, J.A., Stowe, L.G. (1976) Climatic patterns and the distribution
#'   of C4 grasses in North America. *Oecologia* 23, 1-12.
#'
#' UNEP (1992) World atlas of desertification. Edward Arnold, London.
#'
#' Walter, H., Lieth, H. (1960) Klimadiagramm-Weltatlas. Fischer, Jena.
#'
#' Wang, T., Hamann, A., Spittlehouse, D.L., Murdock, T.Q. (2012)
#'   ClimateWNA: high-resolution spatial climate data for western North
#'   America. *Journal of Applied Meteorology and Climatology* 51, 16-29.
#'
#' Way, D.A., Montgomery, R.A. (2015) Photoperiod constraints on tree
#'   phenology, performance and migration in a warming world. *Plant, Cell
#'   and Environment* 38, 1725-1736.
#'
#' Woodward, F.I. (1987) Climate and plant distribution. Cambridge
#'   University Press.
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
#' # only some: by group and by name, with another base temperature
#' names(climate_indices(tmax, tmin, prec, indices = c("cold", "gdd"),
#'                       gdd_base = 10))
climate_indices <- function(tmax, tmin, prec, tmean = NULL,
                            indices = "all",
                            heat_threshold = 30, gdd_base = 5,
                            frost_threshold = 0, tdew_offset = 0,
                            temp_units = c("auto", "C", "K", "K10"),
                            path = NULL, prefix = "", overwrite = FALSE) {
  temp_units <- match.arg(temp_units)
  indices <- expand_indices(indices)
  tab <- climate_index_table()
  dem <- intersect(indices, tab$index[tab$source == "dem"])
  if (length(dem)) {
    stop("Terrain indices (", paste(dem, collapse = ", "), ") come from an ",
         "elevation raster: use topo_indices().", call. = FALSE)
  }
  fl <- intersect(indices, tab$index[nzchar(tab$flag)])
  if (length(fl)) {
    message("Computing flagged variable(s) ", paste(fl, collapse = ", "),
            ": they combine temperature and precipitation over periods ",
            "that switch between neighbouring cells, which creates spatial ",
            "artefacts in some regions (Booth 2022). Check their maps ",
            "before use, or use t_wet, t_dry, p_warm and p_cold.")
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
  need_pet <- any(c("pet", "cmd", "aridity", "lgp") %in% indices)
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
  # smooth seasonal indices: monthly weights instead of the warmest or
  # wettest quarter, so that no period switches between neighbouring cells
  wmean <- function(v, w) {
    sw <- sum(w)
    terra::ifel(sw > 0, sum(v * w) / sw, terra::mean(v))
  }
  if ("t_wet" %in% indices) make("t_wet", wmean(tmean, prec))
  if ("t_dry" %in% indices) make("t_dry", wmean(tmean, max(prec) - prec))
  if ("p_warm" %in% indices) {
    make("p_warm", 3 * wmean(prec, tmean - min(tmean)))
  }
  if ("p_cold" %in% indices) {
    make("p_cold", 3 * wmean(prec, max(tmean) - tmean))
  }
  # bioclim variables (definitions of WorldClim; quarters are any three
  # consecutive months, also across the year's end)
  bio <- paste0("bio", 1:19)
  if (any(bio %in% indices)) {
    q <- function(m) ((m - 1 + 0:2) %% 12) + 1
    tq <- terra::rast(lapply(1:12, function(m) terra::mean(tmean[[q(m)]])))
    pq <- terra::rast(lapply(1:12, function(m) sum(prec[[q(m)]])))
    sd12 <- function(x) sqrt(sum((x - terra::mean(x))^2) / 11)
    b2 <- terra::mean(tmax - tmin)
    b7 <- max(tmax) - min(tmin)
    if ("bio1" %in% indices) make("bio1", terra::mean(tmean))
    if ("bio2" %in% indices) make("bio2", b2)
    if ("bio3" %in% indices) make("bio3", terra::ifel(b7 > 0, 100 * b2 / b7,
                                                      NA))
    if ("bio4" %in% indices) make("bio4", 100 * sd12(tmean))
    if ("bio5" %in% indices) make("bio5", max(tmax))
    if ("bio6" %in% indices) make("bio6", min(tmin))
    if ("bio7" %in% indices) make("bio7", b7)
    if ("bio8" %in% indices) {
      make("bio8", terra::selectRange(tq, terra::which.max(pq)))
    }
    if ("bio9" %in% indices) {
      make("bio9", terra::selectRange(tq, terra::which.min(pq)))
    }
    if ("bio10" %in% indices) make("bio10", max(tq))
    if ("bio11" %in% indices) make("bio11", min(tq))
    if ("bio12" %in% indices) make("bio12", sum(prec))
    if ("bio13" %in% indices) make("bio13", max(prec))
    if ("bio14" %in% indices) make("bio14", min(prec))
    if ("bio15" %in% indices) {
      make("bio15", 100 * sd12(prec + 1) / terra::mean(prec + 1))
    }
    if ("bio16" %in% indices) make("bio16", max(pq))
    if ("bio17" %in% indices) make("bio17", min(pq))
    if ("bio18" %in% indices) {
      make("bio18", terra::selectRange(pq, terra::which.max(tq)))
    }
    if ("bio19" %in% indices) {
      make("bio19", terra::selectRange(pq, terra::which.min(tq)))
    }
  }
  if ("gdd" %in% indices) {
    make("gdd", sum(terra::clamp(tmean - gdd_base, lower = 0,
                                 values = TRUE) * days))
  }
  if ("cold_stress" %in% indices) {
    make("cold_stress", sum(terra::clamp(frost_threshold - tmin, lower = 0,
                                         values = TRUE)))
  }
  if ("frost_months" %in% indices) {
    make("frost_months", sum(tmin < frost_threshold))
  }
  if ("continentality" %in% indices) {
    make("continentality", max(tmean) - min(tmean))
  }
  if ("dry_months" %in% indices) {
    make("dry_months", sum(prec < 2 * tmean))
  }
  if ("lgp" %in% indices) {
    make("lgp", sum(tmean >= gdd_base & prec >= 0.5 * pet_m))
  }
  if ("emberger" %in% indices) {
    hot <- max(tmax) + 273.15
    cold <- min(tmin) + 273.15
    # no annual temperature range: the quotient is undefined
    make("emberger", terra::ifel(hot > cold,
                                 2000 * sum(prec) / (hot^2 - cold^2), NA))
  }
  if ("ellenberg" %in% indices) {
    tot <- sum(prec)
    make("ellenberg", terra::ifel(tot > 0, 1000 * max(tmean) / tot, NA))
  }
  if ("kira_warmth" %in% indices) {
    make("kira_warmth", sum(terra::clamp(tmean - 5, lower = 0,
                                         values = TRUE)))
  }
  if ("kira_coldness" %in% indices) {
    make("kira_coldness", sum(terra::clamp(tmean - 5, upper = 0,
                                           values = TRUE)))
  }
  light <- c("ra", "ra_gs", "daylength_start", "daylength_end",
             "daylength_gs", "daylight_gs", "daylength_max",
             "days_to_solstice")
  if (any(light %in% indices)) {
    phi <- latitude_radians(tmax[[1]])
  }
  if (any(c("ra", "ra_gs") %in% indices)) {
    ra_m <- terra::rast(lapply(1:12, function(i) {
      extraterrestrial_radiation(phi, i)
    }))
    # cells without climate values get none either
    ra_m <- terra::mask(ra_m, tmean[[1]])
    if ("ra" %in% indices) make("ra", sum(ra_m * days) / sum(days))
    if ("ra_gs" %in% indices) {
      make("ra_gs", sum(ra_m * days * (tmean >= gdd_base)))
    }
  }
  if ("daylength_start" %in% indices) {
    make("daylength_start", terra::app(c(tmean, phi), season_start_daylength,
                                       base = gdd_base))
  }
  if ("daylength_end" %in% indices) {
    make("daylength_end", terra::app(c(tmean, phi), season_light,
                                     base = gdd_base, what = "end"))
  }
  if (any(c("daylength_gs", "daylight_gs") %in% indices)) {
    dl_m <- terra::rast(lapply(1:12, function(i) day_length(phi, mid_day[i])))
    grow <- tmean >= gdd_base
    hours <- sum(dl_m * days * grow)
    if ("daylength_gs" %in% indices) {
      n <- sum(grow * days)
      # no month warm enough: no growing season, no value
      make("daylength_gs", terra::ifel(n > 0, hours / n, NA))
    }
    if ("daylight_gs" %in% indices) make("daylight_gs", hours)
  }
  if ("daylength_max" %in% indices) {
    # the summer solstice of either hemisphere (declination 0.409 rad)
    make("daylength_max", terra::mask(day_length(abs(phi), 172), tmean[[1]]))
  }
  if ("days_to_solstice" %in% indices) {
    make("days_to_solstice", terra::app(c(tmean, phi), season_light,
                                        base = gdd_base, what = "solstice"))
  }
  terra::rast(out[intersect(indices, names(out))])
}

#' Names of the climate indices, by group
#'
#' The indices that [climate_indices()] can compute, to choose among them:
#' all of them, those of some groups, or all but some.
#'
#' | Group | Indices |
#' |---|---|
#' | `water` | `vpd`, `pet`, `cmd`, `aridity`, `dry_months`, `emberger`, `ellenberg` |
#' | `heat` | `heat_stress`, `gdd`, `kira_warmth`, `continentality` |
#' | `cold` | `cold_stress`, `frost_months`, `kira_coldness` |
#' | `season` | `sri`, `synchronicity`, `lgp` |
#' | `light` | `ra`, `ra_gs`, `daylength_start`, `daylength_end`, `daylength_gs`, `daylight_gs`, `daylength_max`, `days_to_solstice` |
#'
#' @param groups Groups to return (default: all).
#' @param drop Indices or groups to leave out.
#' @param table Return a data.frame with the group, meaning, unit,
#'   ecological relevance and literature source of each index instead of
#'   the names.
#' @return A character vector of index names, or a data.frame.
#' @export
#' @examples
#' climate_index_names()
#' climate_index_names("water")
#' climate_index_names(drop = c("light", "emberger"))
#' climate_index_names(table = TRUE)
climate_index_names <- function(groups = NULL, drop = NULL, table = FALSE,
                                flagged = FALSE, source = NULL) {
  tab <- climate_index_table()
  if (!is.null(source)) {
    source <- match.arg(source, c("monthly", "dem"))
    tab <- tab[tab$source == source, , drop = FALSE]
  }
  known <- unique(tab$group)
  bad <- setdiff(groups, known)
  if (length(bad)) {
    stop("Unknown group: ", paste(bad, collapse = ", "), ". Available: ",
         paste(known, collapse = ", "), call. = FALSE)
  }
  bad <- setdiff(drop, c(known, tab$index))
  if (length(bad)) {
    stop("Unknown index or group in `drop`: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }
  keep <- (is.null(groups) | tab$group %in% groups) &
    !(tab$index %in% drop | tab$group %in% drop)
  out_fl <- keep & nzchar(tab$flag) & !flagged
  if (any(out_fl) && !isTRUE(.preenm_env$flag_note)) {
    message("Left out: ", paste(tab$index[out_fl], collapse = ", "),
            ". They combine temperature and precipitation over the warmest, ",
            "coldest, wettest or driest quarter (or month), which can switch ",
            "between neighbouring cells and create spatial artefacts ",
            "(Booth 2022). flagged = TRUE includes them; the smooth ",
            "t_wet, t_dry, p_warm and p_cold measure the same contrasts ",
            "without switching. (Shown once per session.)")
    .preenm_env$flag_note <- TRUE
  }
  keep <- keep & !out_fl
  tab <- tab[keep, , drop = FALSE]
  rownames(tab) <- NULL
  if (table) tab else tab$index
}

# The indices in the order they are returned (internal)
climate_index_table <- function() {
  x <- c(
    "vpd",
    "water",
    "mean vapour pressure deficit (dew point from Tmin)",
    "kPa",
    "dryness of the air: stomatal closure and transpiration of plants, desiccation of amphibians, snails and insects",
    "Allen et al. 1998; Grossiord et al. 2020",
    "pet",
    "water",
    "potential evapotranspiration (Hargreaves-Samani)",
    "mm/year",
    "evaporative demand and available energy; plant productivity and species richness",
    "Hargreaves & Samani 1985; Allen et al. 1998",
    "cmd",
    "water",
    "climatic moisture deficit: sum of monthly PET - P where positive",
    "mm/year",
    "drought not met by rain: forest-steppe limits, drought-sensitive trees, C3 versus C4 grasses",
    "Stephenson 1990; Wang et al. 2012",
    "aridity",
    "water",
    "aridity index P / PET (< 0.05 hyper-arid, > 0.65 humid)",
    "-",
    "the humid-arid gradient: biome limits, dryland and desert taxa",
    "UNEP 1992",
    "heat_stress",
    "heat",
    "sum over months of Tmax above heat_threshold",
    "degC x months",
    "exposure to hot months: upper thermal limits of alpine and boreal taxa and of ectotherms; C4 advantage",
    "this package; Ehleringer et al. 1997",
    "sri",
    "season",
    "summer rain index: precipitation of the warmest quarter / annual (bio18 / bio12)",
    "0-1",
    "rain in the warm season: C4 versus C3 grasses, monsoon versus Mediterranean taxa",
    "Teeri & Stowe 1976; Paruelo & Lauenroth 1996",
    "synchronicity",
    "season",
    "heat-rain synchrony: 1 = the warmest and the wettest month coincide, 0 = six months apart",
    "0-1",
    "whether water arrives when it is warm: summer-rain versus winter-rain floras, C3/C4 balance",
    "this package; Paruelo & Lauenroth 1996",
    "gdd",
    "heat",
    "growing degree days: monthly Tmean above gdd_base x days",
    "degC x days",
    "heat sum for growth and development: cold range limits of plants, number of insect generations",
    "Woodward 1987; Prentice et al. 1992",
    "cold_stress",
    "cold",
    "sum over months of Tmin below frost_threshold",
    "degC x months",
    "winter severity: frost-sensitive plants, survival of overwintering stages",
    "Woodward 1987; Sakai & Larcher 1987",
    "frost_months",
    "cold",
    "number of months with Tmin below frost_threshold",
    "months",
    "duration of frost: evergreen versus deciduous plants, length of dormancy or hibernation",
    "Woodward 1987; Sakai & Larcher 1987",
    "continentality",
    "heat",
    "Tmean of the warmest month - that of the coldest",
    "degC",
    "annual temperature amplitude: oceanic versus continental species",
    "Rivas-Martinez et al. 2011",
    "dry_months",
    "water",
    "months with P < 2 x Tmean (the dry season of Walter-Lieth diagrams)",
    "months",
    "length of the dry season: Mediterranean and savanna taxa, drought-deciduous plants",
    "Bagnouls & Gaussen 1953; Walter & Lieth 1960",
    "lgp",
    "season",
    "length of the growing period: months with Tmean >= gdd_base and P >= PET / 2",
    "months",
    "time when both warmth and water allow growth: crops and their wild relatives, annual plants",
    "FAO 1978",
    "emberger",
    "water",
    "Emberger's pluviothermic quotient 2000 P / (M^2 - m^2), M and m in Kelvin",
    "-",
    "Mediterranean bioclimatic stages (arid to humid); lower = drier",
    "Emberger 1955; Daget 1977",
    "ellenberg",
    "water",
    "Ellenberg's quotient: Tmean of the warmest month / annual P x 1000",
    "degC/mm x 1000",
    "humidity gradient of temperate forests: beech dominates below about 30, oaks above",
    "Ellenberg 1963; Fang & Lechowicz 2006",
    "kira_warmth",
    "heat",
    "Kira's warmth index: sum over months of Tmean - 5 where Tmean > 5 degC",
    "degC x months",
    "vegetation zones of East Asia (evergreen broad-leaved forest between 85 and 180)",
    "Kira 1945; Kira 1991",
    "kira_coldness",
    "cold",
    "Kira's coldness index: sum over months of Tmean - 5 where Tmean < 5 degC (zero or negative)",
    "degC x months",
    "cold limit of evergreen broad-leaved trees in East Asia (about -10)",
    "Kira 1945; Kira 1991",
    "ra",
    "light",
    "mean extraterrestrial radiation over the year (latitude only)",
    "MJ m-2 day-1",
    "light at the top of the atmosphere; in a model a smooth stand-in for latitude",
    "Allen et al. 1998",
    "ra_gs",
    "light",
    "extraterrestrial radiation received in the months with Tmean >= gdd_base",
    "MJ m-2",
    "light available while it is warm enough to grow: productivity at high latitudes and elevations",
    "this package; Allen et al. 1998",
    "daylength_start",
    "light",
    "day length in the middle of the month in which the growing season starts",
    "hours",
    "photoperiod at the onset of growth: bud burst, flowering of long-day plants, end of diapause",
    "Koerner & Basler 2010; Way & Montgomery 2015",
    "daylength_end",
    "light",
    "day length in the middle of the month in which the growing season ends",
    "hours",
    "photoperiod at the end of the season: growth cessation and bud set of trees, start of diapause",
    "Bradshaw & Holzapfel 2007; Way & Montgomery 2015",
    "daylength_gs",
    "light",
    "mean day length of the months with Tmean >= gdd_base",
    "hours",
    "photoperiod during growth: long-day versus short-day plants, daily activity time of diurnal animals",
    "Saikkonen et al. 2012",
    "daylight_gs",
    "light",
    "hours of daylight in the months with Tmean >= gdd_base",
    "hours",
    "time for photosynthesis or foraging in the season: a short northern summer with long days against a long southern one",
    "this package; Saikkonen et al. 2012",
    "daylength_max",
    "light",
    "day length on the summer solstice (latitude only; 24 = midnight sun)",
    "hours",
    "the longest day: species that need a critical day length to flower or to breed; like ra a stand-in for latitude",
    "Allen et al. 1998; Bradshaw & Holzapfel 2007",
    "days_to_solstice",
    "light",
    "days from the start of the growing season to the summer solstice (negative = it starts later)",
    "days",
    "share of the season on lengthening days: phenological timing, species whose flowering is triggered before the solstice",
    "this package; Koerner & Basler 2010")
  m <- matrix(x, ncol = 6, byrow = TRUE)
  tab <- data.frame(index = m[, 1], group = m[, 2], meaning = m[, 3],
                    unit = m[, 4], relevance = m[, 5], reference = m[, 6],
                    stringsAsFactors = FALSE)
  tab <- rbind(tab, smooth_season_table(), bioclim_table(), terrain_table())
  comp <- paste0("combines temperature and precipitation over a period ",
                 "that switches between cells (artefacts; Booth 2022)")
  tab$flag <- ifelse(tab$index %in% c("bio8", "bio9", "bio18", "bio19",
                                      "sri", "synchronicity"), comp, "")
  tab$source <- ifelse(tab$group == "terrain", "dem", "monthly")
  rownames(tab) <- NULL
  tab
}

# Package-level state (internal)
.preenm_env <- new.env(parent = emptyenv())

# Smooth seasonal indices (internal)
smooth_season_table <- function() {
  data.frame(
    index = c("t_wet", "t_dry", "p_warm", "p_cold"),
    group = "season",
    meaning = c(
      "mean temperature weighted by monthly precipitation (smooth bio8)",
      "mean temperature weighted by monthly dryness, max(P) - P (smooth bio9)",
      "precipitation weighted by monthly warmth, T - min(T), x 3 (smooth bio18)",
      "precipitation weighted by monthly coldness, max(T) - T, x 3 (smooth bio19)"),
    unit = c("degC", "degC", "mm/quarter", "mm/quarter"),
    relevance = c(
      "temperature when water is available: the contrast of bio8 without switching quarters",
      "temperature of the dry season: the contrast of bio9 without switching quarters",
      "rain in the warm season: the contrast of bio18 (and sri) without switching quarters",
      "rain in the cold season: the contrast of bio19 without switching quarters"),
    reference = "this package; after Booth 2022",
    stringsAsFactors = FALSE)
}

# Bioclim variables (internal)
bioclim_table <- function() {
  m <- c(
    "annual mean temperature", "degC",
    "mean diurnal range: mean of monthly (Tmax - Tmin)", "degC",
    "isothermality: bio2 / bio7 x 100", "%",
    "temperature seasonality: SD of monthly means x 100", "degC x 100",
    "maximum temperature of the warmest month", "degC",
    "minimum temperature of the coldest month", "degC",
    "temperature annual range: bio5 - bio6", "degC",
    "mean temperature of the wettest quarter", "degC",
    "mean temperature of the driest quarter", "degC",
    "mean temperature of the warmest quarter", "degC",
    "mean temperature of the coldest quarter", "degC",
    "annual precipitation", "mm",
    "precipitation of the wettest month", "mm",
    "precipitation of the driest month", "mm",
    "precipitation seasonality: CV of monthly (P + 1) x 100", "%",
    "precipitation of the wettest quarter", "mm",
    "precipitation of the driest quarter", "mm",
    "precipitation of the warmest quarter", "mm",
    "precipitation of the coldest quarter", "mm")
  m <- matrix(m, ncol = 2, byrow = TRUE)
  data.frame(index = paste0("bio", 1:19), group = "bioclim",
             meaning = m[, 1], unit = m[, 2],
             relevance = "the standard bioclimatic variables, computed from the same monthly layers",
             reference = "Busby 1991; Hijmans et al. 2005; O'Donnell & Ignizio 2012",
             stringsAsFactors = FALSE)
}

# Terrain indices from an elevation raster, see topo_indices() (internal)
terrain_table <- function() {
  data.frame(
    index = c("slope", "northness", "eastness", "tri", "tpi", "dev",
              "twi", "hli"),
    group = "terrain",
    meaning = c(
      "slope",
      "northness: cos(aspect) x sin(slope), 0 on flat ground",
      "eastness: sin(aspect) x sin(slope), 0 on flat ground",
      "terrain ruggedness index (Riley): root of the summed squared differences to the 8 neighbours",
      "topographic position index: elevation minus the mean within a radius (one layer per radius)",
      "deviation from mean elevation: TPI / SD of elevation within the radius (scale-free TPI)",
      "topographic wetness index: ln(specific catchment area / tan(slope)), D8 flow",
      "heat load index: potential direct radiation from slope, aspect and latitude"),
    unit = c("degrees", "-1 to 1", "-1 to 1", "m", "m", "-", "-", "-"),
    relevance = c(
      "drainage, soil depth, erosion; steep rocky habitats",
      "cool, moist north-facing (south-facing in the south) slopes against warm, dry ones",
      "morning against afternoon sun",
      "habitat heterogeneity, rock outcrops, shelter",
      "ridges and convex slopes (positive) against valleys and hollows (negative); cold-air pools",
      "position relative to local relief, comparable between flat and mountainous terrain",
      "water accumulation: wetlands, moist valley floors, drought refugia",
      "warm and cold slopes: microclimate of plants and ectotherms"),
    reference = c(
      "Horn 1981", "this package; after Guisan et al. 1999",
      "this package; after Guisan et al. 1999", "Riley et al. 1999",
      "Weiss 2001; Guisan et al. 1999", "De Reu et al. 2013",
      "Beven & Kirkby 1979", "McCune & Keon 2002"),
    stringsAsFactors = FALSE)
}

# Index names from names, group names or "all" (internal)
expand_indices <- function(indices) {
  tab <- climate_index_table()
  # "all": every index computed from monthly climate, without the flagged
  if (identical(indices, "all")) {
    return(tab$index[tab$source == "monthly" & !nzchar(tab$flag)])
  }
  bad <- setdiff(indices, c(tab$index, tab$group))
  if (length(bad)) {
    stop("Unknown index: ", paste(bad, collapse = ", "), ". Available: ",
         paste(tab$index, collapse = ", "), "; or a group: ",
         paste(unique(tab$group), collapse = ", "), call. = FALSE)
  }
  # a group brings its unflagged members; flagged ones only by name
  tab$index[tab$index %in% indices |
              (tab$group %in% indices & !nzchar(tab$flag))]
}

# Middle day of each month (day of the year), as in FAO-56
mid_day <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)

# Day length (hours) at latitude `phi` (radians; numbers or a SpatRaster) on
# day `j` of the year: 24 / pi x the sunset hour angle, FAO-56 equations 24,
# 25 and 34 (internal)
day_length <- function(phi, j) {
  delta <- 0.409 * sin(2 * pi / 365 * j - 1.39)
  arg <- -tan(phi) * tan(delta)
  arg <- if (inherits(arg, "SpatRaster")) {
    terra::clamp(arg, lower = -1, upper = 1, values = TRUE)
  } else pmin(pmax(arg, -1), 1)
  24 / pi * acos(arg)
}

# Day length (hours) in the middle of the month in which the growing season
# starts: the first month at or above `base`, counting from the coldest
# month. `m` has one row per cell: 12 monthly mean temperatures and the
# latitude in radians (internal)
season_start_daylength <- function(m, base) season_light(m, base, "start")

# Light at the limits of the growing season, for cells in the rows of `m`
# (12 monthly mean temperatures and the latitude in radians). The season
# starts in the first month at or above `base` counting forwards from the
# coldest month, and ends in the first such month counting backwards from
# it. `what`: "start" / "end" = day length in the middle of that month;
# "solstice" = days from the start of the season to the summer solstice,
# with the start placed between two mid-month days by linear interpolation
# of the temperatures (internal)
season_light <- function(m, base, what = c("start", "end", "solstice")) {
  what <- match.arg(what)
  out <- rep(NA_real_, nrow(m))
  ok <- which(stats::complete.cases(m))
  if (!length(ok)) return(out)
  t <- m[ok, 1:12, drop = FALSE]
  phi <- m[ok, 13]
  row <- seq_along(ok)
  cold <- max.col(-t, ties.method = "first")
  hit_month <- rep(NA_integer_, length(ok))
  step <- if (what == "end") -1 else 1
  for (k in 0:11) {
    month <- (cold - 1 + step * k) %% 12 + 1
    hit <- is.na(hit_month) & t[cbind(row, month)] >= base
    hit_month[hit] <- month[hit]
  }
  if (what != "solstice") {
    out[ok] <- day_length(phi, mid_day[hit_month])
    return(out)
  }
  # the day on which the interpolated temperature reaches the base
  onset <- mid_day[hit_month]
  prev <- (hit_month - 2) %% 12 + 1
  rising <- !is.na(hit_month) & hit_month != cold
  i <- which(rising)
  t0 <- t[cbind(i, prev[i])]
  t1 <- t[cbind(i, hit_month[i])]
  gap <- (mid_day[hit_month[i]] - mid_day[prev[i]]) %% 365
  onset[i] <- mid_day[prev[i]] + (base - t0) / (t1 - t0) * gap
  solstice <- ifelse(phi >= 0, 172, 355)
  d <- (solstice - onset) %% 365
  out[ok] <- ifelse(d > 182.5, d - 365, d)
  out
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
