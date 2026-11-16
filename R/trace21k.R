# CHELSA-TraCE21k: climate of the last 22,000 years ---------------------------------

trace21k_base <- "https://os.zhdk.cloud.switch.ch/chelsav1/chelsa_trace"

#' Time steps of CHELSA-TraCE21k
#'
#' CHELSA-TraCE21k (Karger et al. 2023) has one layer per century, named
#' by a time id from 20 (1900-1990 CE, the present of the data set) down to
#' -200 (22,000 years before present). This function converts between
#' time ids and thousands of years before present (ka BP; present = 1950).
#'
#' @param kyr_bp Time in thousands of years before present, 0 to 22 (e.g.
#'   `21` for the Last Glacial Maximum, `12` for the Younger Dryas, `6` for
#'   the mid-Holocene); rounded to the century.
#' @param time_id Time ids, to convert them back: give either `kyr_bp` or
#'   `time_id`.
#' @return A data.frame with `time_id`, `kyr_bp` and the first and last
#'   year of the century (`year_start`, `year_end`; negative = BCE as
#'   written in the data set's documentation).
#' @references Karger, D.N., Nobis, M.P., Normand, S., Graham, C.H.,
#'   Zimmermann, N.E. (2023) CHELSA-TraCE21k: high-resolution (1 km)
#'   downscaled transient temperature and precipitation data since the Last
#'   Glacial Maximum. *Climate of the Past* 19, 439-456.
#'   \doi{10.5194/cp-19-439-2023}
#' @seealso [trace21k_id()] for the time id of ages in years.
#' @export
#' @examples
#' trace21k_time(c(0, 6, 12, 21))
#' trace21k_time(time_id = c(20, -190))
trace21k_time <- function(kyr_bp = NULL, time_id = NULL) {
  if (is.null(time_id)) {
    if (is.null(kyr_bp)) stop("Give `kyr_bp` or `time_id`.", call. = FALSE)
    if (!is.numeric(kyr_bp) || anyNA(kyr_bp) || any(kyr_bp < 0 | kyr_bp > 22)) {
      stop("`kyr_bp` must be between 0 and 22.", call. = FALSE)
    }
    time_id <- 20 - round(kyr_bp * 10)
  }
  if (!is.numeric(time_id) || anyNA(time_id) ||
      any(time_id != round(time_id)) || any(time_id < -200 | time_id > 20)) {
    stop("Time ids are whole numbers from -200 to 20.", call. = FALSE)
  }
  start <- 100 * (time_id - 1)
  data.frame(time_id = as.integer(time_id), kyr_bp = (20 - time_id) / 10,
             year_start = start,
             year_end = ifelse(time_id == 20, 1990, start + 99))
}

#' URLs of CHELSA-TraCE21k files
#'
#' The addresses of the files of CHELSA-TraCE21k V1.0 for one or more time
#' steps: the monthly layers needed by [climate_indices()] and the 19
#' bioclimatic variables.
#'
#' @param vars Monthly variables (`"tasmax"`, `"tasmin"`, `"pr"`) and/or
#'   bioclimatic variables (`"bio01"` to `"bio19"`, or `"bio"` for all 19).
#' @inheritParams trace21k_time
#' @param months Months of the monthly variables.
#' @return A data.frame with `var`, `month` (`NA` for bioclimatic
#'   variables), `time_id`, `kyr_bp`, `url` and `file`.
#' @inherit trace21k_time references
#' @export
#' @examples
#' head(trace21k_urls(kyr_bp = 21))
#' trace21k_urls(c("bio01", "bio12"), kyr_bp = c(0, 21))$file
trace21k_urls <- function(vars = c("tasmax", "tasmin", "pr"), kyr_bp = NULL,
                          time_id = NULL, months = 1:12) {
  tm <- trace21k_time(kyr_bp, time_id)
  if ("bio" %in% vars) {
    vars <- unique(c(setdiff(vars, "bio"), sprintf("bio%02d", 1:19)))
  }
  # bio1 -> bio01
  short <- grepl("^bio[0-9]$", vars)
  vars[short] <- sub("^bio", "bio0", vars[short])
  bad <- vars[!(vars %in% c("tasmax", "tasmin", "pr") |
                  vars %in% sprintf("bio%02d", 1:19))]
  if (length(bad)) {
    stop("Unknown CHELSA-TraCE21k variable: ", paste(bad, collapse = ", "),
         ". Available: tasmax, tasmin, pr, bio01 to bio19.", call. = FALSE)
  }
  if (!all(months %in% 1:12)) stop("`months` must be 1 to 12.", call. = FALSE)
  rows <- lapply(seq_len(nrow(tm)), function(i) {
    do.call(rbind, lapply(vars, function(v) {
      is_bio <- grepl("^bio", v)
      m <- if (is_bio) NA_integer_ else as.integer(months)
      file <- if (is_bio) {
        sprintf("CHELSA_TraCE21k_%s_%d_V1.0.tif", v, tm$time_id[i])
      } else {
        sprintf("CHELSA_TraCE21k_%s_%d_%d_V1.0.tif", v, m, tm$time_id[i])
      }
      data.frame(var = v, month = m, time_id = tm$time_id[i],
                 kyr_bp = tm$kyr_bp[i],
                 url = paste(trace21k_base, if (is_bio) "bio" else v, file,
                             sep = "/"),
                 file = file, stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows)
}

#' Download CHELSA-TraCE21k for past periods, optionally for a region only
#'
#' Downloads the files of CHELSA-TraCE21k V1.0 (Karger et al. 2023): the
#' climate of every century of the last 22,000 years at about 1 km,
#' downscaled from the TraCE-21ka simulation with changing ice sheets and
#' sea level. By default the monthly maximum and minimum temperature and
#' precipitation of the chosen time steps, from which [climate_indices()]
#' computes the same derived variables as for the present; with `vars =
#' "bio"` the 19 bioclimatic variables. As in [download_chelsa_monthly()],
#' with `extent` only that region is read from the global files over the
#' internet, existing files are kept, and an interrupted download
#' continues when run again.
#'
#' **Units.** Check them once after downloading (`terra::rast(file)`): the
#' data set stores temperatures in tenths of Kelvin, which
#' [climate_indices()] recognises and converts (`temp_units = "auto"`), and
#' precipitation in kg m-2 per month (= mm/month). Some bioclimatic layers
#' are stored multiplied by 10; [sample_environment()] with `fix_scale =
#' TRUE` detects layers whose scale differs between data sets. The data set
#' has no monthly mean temperature; [climate_indices()] then uses
#' `(tmax + tmin) / 2`. Compare past with present from the same source
#' (time step 0 ka BP of this data set, or an anomaly, see
#' [sample_environment()]) so that differences between sources are not
#' read as climate change.
#'
#' @inheritParams trace21k_urls
#' @param path Folder for the files (created if needed).
#' @param extent Optional region, as in [download_chelsa_monthly()].
#' @param overwrite Download again files that exist.
#' @return The paths of the files, invisibly.
#' @inherit trace21k_time references
#' @seealso [trace21k_time()], [monthly_files()], [climate_indices()]
#' @export
#' @examples
#' \dontrun{
#' # monthly climate of the Alps at the Last Glacial Maximum and today
#' download_chelsa_trace21k("trace_alps", kyr_bp = c(0, 21),
#'                          extent = c(5, 17, 43, 49))
#' # the files of one time step are told apart by its time id (-190)
#' f <- function(v) monthly_files("trace_alps", v, must_contain = "_-190_V1")
#' lgm <- climate_indices(f("tasmax"), f("tasmin"), f("pr"),
#'                        path = "indices", prefix = "TraCE_21ka_")
#' # bioclimatic variables of three periods
#' download_chelsa_trace21k("trace_alps", vars = "bio", kyr_bp = c(0, 12, 21),
#'                          extent = c(5, 17, 43, 49))
#' }
download_chelsa_trace21k <- function(path, vars = c("tasmax", "tasmin", "pr"),
                                     kyr_bp = NULL, time_id = NULL,
                                     months = 1:12, extent = NULL,
                                     overwrite = FALSE) {
  u <- trace21k_urls(vars, kyr_bp, time_id, months)
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  download_rasters(u$url, file.path(path, u$file), extent, overwrite)
}
