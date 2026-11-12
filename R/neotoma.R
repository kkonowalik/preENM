# Fossil records from Neotoma ----------------------------------------------------

#' Download Neotoma samples
#'
#' Fetches the full downloads of Neotoma datasets (package `neotoma2`) and
#' returns their samples, one row per sample and taxon, all columns as text.
#' Find the dataset IDs in the Neotoma Explorer (export, column
#' `datasetid`). Process the result with [neotoma_records()].
#'
#' @param dataset_ids Neotoma dataset IDs.
#' @param chunk Datasets per request.
#' @param cache Optional `.rds` file: re-used if it exists, written
#'   otherwise (delete it to download again).
#' @return A data.frame of samples.
#' @export
#' @examples
#' \dontrun{
#' smp <- download_neotoma(c(1435, 1436), cache = "neotoma_downloads.rds")
#' rec <- neotoma_records(smp, taxon_regex = "^Dryas")
#' }
download_neotoma <- function(dataset_ids, chunk = 25, cache = NULL) {
  if (!is.null(cache) && file.exists(cache)) {
    d <- readRDS(cache)
    if (nrow(d)) return(d)
  }
  if (!requireNamespace("neotoma2", quietly = TRUE)) {
    stop("Install neotoma2: install.packages(\"neotoma2\")", call. = FALSE)
  }
  ids <- sort(unique(as.integer(dataset_ids)))
  fetch <- function(p) {
    dl <- neotoma2::get_downloads(as.numeric(p))
    s <- as.data.frame(neotoma2::samples(dl))
    s[] <- lapply(s, as.character)
    s
  }
  # first dataset without error catching, so a real problem stops here
  first <- fetch(ids[1])
  parts <- split(ids[-1], ceiling(seq_along(ids[-1]) / chunk))
  errs <- character(0)
  rest <- lapply(seq_along(parts), function(k) {
    message("Downloading chunk ", k, "/", length(parts))
    tryCatch(fetch(parts[[k]]), error = function(e) {
      errs <<- c(errs, conditionMessage(e))
      NULL
    })
  })
  all <- c(list(first), Filter(Negate(is.null), rest))
  cols <- unique(unlist(lapply(all, names)))
  d <- do.call(rbind, lapply(all, function(x) {
    x[setdiff(cols, names(x))] <- NA_character_
    x[cols]
  }))
  if (length(errs)) {
    warning(length(errs), " chunk(s) failed (first error: ", errs[1],
            "); re-run to retry.", call. = FALSE)
  }
  if (!is.null(cache) && nrow(d)) saveRDS(d, cache)
  d
}

#' Prepare fossil samples from Neotoma for temporal validation
#'
#' Turns Neotoma samples (from [download_neotoma()] or any table with the
#' same columns) into dated records of your taxon: each sample with its
#' **own** age, an age window for matching it to climate time slices, and a
#' flag explaining why a sample is not usable. Uncalibrated radiocarbon ages
#' are calibrated with IntCal20 (package `rintcal`; see
#' [calibrate_radiocarbon()]).
#'
#' **Age window.** The stated older/younger bounds, or +/- `age_unc_default`
#' years, never narrower than +/- `age_unc_min`. The window is converted to
#' climate time slices by `slice_fun` (default: the century IDs of
#' CHELSA-TraCE21k, see [trace21k_id()]).
#'
#' **Flags:** `ok`, `outside_area` (outside `area`, e.g. M, which also
#' removes regions where other taxa of the genus occur), `zero_value`,
#' `no_age`, `not_calendar_age`, `older_than_max`.
#'
#' @param samples Sample table (columns found case-insensitively: siteid,
#'   sitename, lat, long/lon, datasetid, datasettype, sampleid, depth,
#'   variablename/taxonname, elementtype, value, units, age, ageolder,
#'   ageyounger, agetype).
#' @param taxon_regex Regular expression for the taxon names kept, e.g.
#'   `"^Dryas"`.
#' @param area Optional area: a raster (cells with values = inside) or
#'   polygons (`sf`, [terra::SpatVector] or path).
#' @param min_value Samples with `value <= min_value` are flagged.
#' @param age_unc_default,age_unc_min Age half-window (years) when none is
#'   stated, and its minimum.
#' @param max_age Oldest age covered by the climate data (cal BP; 22000 for
#'   CHELSA-TraCE21k).
#' @param calibrate_c14 Calibrate radiocarbon ages (needs `rintcal`).
#' @param c14_err_default Laboratory error (14C years) when none is given.
#' @param slice_fun Function converting ages (cal BP) to time-slice IDs.
#' @return A data.frame, one row per sample, with `SITEID`, `SITENAME`,
#'   `LAT`, `LON`, `DATASETID`, `DATASETTYPE`, `SAMPLEID`, `DEPTH`, `TAXON`,
#'   `ELEMENT`, `VALUE`, `UNITS`, `AGE`, `AGE_OLDER`, `AGE_YOUNGER`,
#'   `AGETYPE`, `AGE_SOURCE`, `AGE_MIN`, `AGE_MAX`, `SLICE`, `SLICE_MIN`,
#'   `SLICE_MAX`, `RECORD_TYPE` (pollen / macrofossil) and `FLAG`.
#' @export
#' @examples
#' smp <- data.frame(siteid = c("1", "1", "2"), sitename = c("A", "A", "B"),
#'                   lat = c(47, 47, 60), long = c(11, 11, 10),
#'                   datasetid = "9", datasettype = "pollen",
#'                   variablename = c("Dryas", "Betula", "Dryas octopetala"),
#'                   value = c(3, 10, 1), age = c(12500, 12500, 25000),
#'                   agetype = "Calendar years BP")
#' neotoma_records(smp, taxon_regex = "^Dryas", calibrate_c14 = FALSE)
neotoma_records <- function(samples, taxon_regex, area = NULL, min_value = 0,
                            age_unc_default = 500, age_unc_min = 100,
                            max_age = 22000, calibrate_c14 = TRUE,
                            c14_err_default = 100, slice_fun = trace21k_id) {
  s <- as.data.frame(samples, stringsAsFactors = FALSE)
  pick <- function(c) find_column(names(s), c)
  col <- list(siteid = pick("siteid"), sitename = pick("sitename"),
              lat = pick(c("lat", "latitude")),
              lon = pick(c("long", "lon", "longitude")),
              datasetid = pick("datasetid"),
              datasettype = pick("datasettype"),
              sampleid = pick("sampleid"), depth = pick("depth"),
              taxon = pick(c("variablename", "taxonname")),
              element = pick(c("elementtype", "element")),
              value = pick("value"), units = pick("units"),
              age = pick("age"), older = pick(c("ageolder", "age_older")),
              younger = pick(c("ageyounger", "age_younger")),
              agetype = pick("agetype"))
  need <- c("siteid", "lat", "lon", "taxon", "value", "age")
  miss <- need[vapply(col[need], is.na, logical(1))]
  if (length(miss)) {
    stop("Required columns not found: ", paste(miss, collapse = ", "),
         ". Columns are: ", paste(names(s), collapse = ", "), call. = FALSE)
  }
  g <- function(k) if (is.na(col[[k]])) rep(NA_character_, nrow(s)) else
    as.character(s[[col[[k]]]])
  num <- function(k) suppressWarnings(as.numeric(g(k)))
  d <- data.frame(SITEID = g("siteid"), SITENAME = g("sitename"),
                  LAT = num("lat"), LON = num("lon"),
                  DATASETID = g("datasetid"), DATASETTYPE = g("datasettype"),
                  SAMPLEID = g("sampleid"), DEPTH = num("depth"),
                  TAXON = g("taxon"), ELEMENT = g("element"),
                  VALUE = num("value"), UNITS = g("units"), AGE = num("age"),
                  AGE_OLDER = num("older"), AGE_YOUNGER = num("younger"),
                  AGETYPE = g("agetype"), stringsAsFactors = FALSE)
  d <- d[grepl(taxon_regex, d$TAXON, ignore.case = TRUE), , drop = FALSE]
  if (!nrow(d)) {
    warning("No samples of taxa matching '", taxon_regex, "'.", call. = FALSE)
    return(d)
  }
  d$AGE_SOURCE <- ifelse(is.na(d$AGE), NA_character_, "neotoma")

  # radiocarbon ages -> calendar ages (IntCal20)
  is_c14 <- !is.na(d$AGE) & !is.na(d$AGETYPE) &
    d$AGETYPE == "Radiocarbon years BP"
  if (calibrate_c14 && any(is_c14)) {
    if (!requireNamespace("rintcal", quietly = TRUE)) {
      message("Package rintcal not installed: radiocarbon ages are NOT ",
              "calibrated (install.packages(\"rintcal\")).")
    } else {
      idx <- which(is_c14)
      err <- ifelse(!is.na(d$AGE_OLDER[idx]) & !is.na(d$AGE_YOUNGER[idx]),
                    pmax((d$AGE_OLDER[idx] - d$AGE_YOUNGER[idx]) / 2, 1),
                    c14_err_default)
      cal <- calibrate_radiocarbon(d$AGE[idx], err)
      okc <- !is.na(cal$age)
      d$AGE[idx[okc]] <- cal$age[okc]
      d$AGE_OLDER[idx[okc]] <- cal$older[okc]
      d$AGE_YOUNGER[idx[okc]] <- cal$younger[okc]
      d$AGETYPE[idx[okc]] <- "Calendar years BP (IntCal20, calibrated by preENM)"
      d$AGE_SOURCE[idx[okc]] <- "IntCal20"
      message("Radiocarbon ages calibrated with IntCal20: ", sum(okc), " of ",
              length(idx))
    }
  }
  half <- pmax(ifelse(!is.na(d$AGE_OLDER) & !is.na(d$AGE_YOUNGER),
                      (d$AGE_OLDER - d$AGE_YOUNGER) / 2, age_unc_default),
               age_unc_min)
  now_bp <- 1950 - as.integer(format(Sys.Date(), "%Y"))
  d$AGE_MIN <- pmax(d$AGE - half, now_bp)
  d$AGE_MAX <- d$AGE + half
  d$SLICE <- slice_fun(d$AGE)
  d$SLICE_MIN <- slice_fun(d$AGE_MAX)          # older age = lower ID
  d$SLICE_MAX <- slice_fun(d$AGE_MIN)
  d$RECORD_TYPE <- ifelse(grepl("pollen|spore", paste(d$ELEMENT,
                                                      d$DATASETTYPE),
                                ignore.case = TRUE), "pollen", "macrofossil")

  inside <- rep(TRUE, nrow(d))
  if (!is.null(area)) {
    okxy <- !is.na(d$LON) & !is.na(d$LAT)
    inside[] <- FALSE
    if (any(okxy)) {
      if (inherits(area, "SpatRaster") ||
          (is.character(area) && grepl("\\.(tif|tiff|asc|grd|nc)$", area,
                                       ignore.case = TRUE))) {
        r <- as_grid(area)[[1]]
        inside[okxy] <- !is.na(extract_xy(r, cbind(d$LON[okxy],
                                                   d$LAT[okxy]), "EPSG:4326"))
      } else {
        g2 <- as_sfc(area, "area")
        p <- sf::st_as_sf(data.frame(x = d$LON[okxy], y = d$LAT[okxy]),
                          coords = c("x", "y"), crs = 4326)
        inside[okxy] <- lengths(sf::st_intersects(
          sf::st_transform(p, sf::st_crs(g2)), g2)) > 0
      }
    }
  }
  cal_age <- is.na(d$AGETYPE) | grepl("calendar|cal", d$AGETYPE,
                                      ignore.case = TRUE)
  d$FLAG <- ifelse(!inside, "outside_area",
            ifelse(is.na(d$VALUE) | d$VALUE <= min_value, "zero_value",
            ifelse(is.na(d$AGE), "no_age",
            ifelse(!cal_age, "not_calendar_age",
            ifelse(d$AGE_MIN > max_age, "older_than_max", "ok")))))
  rownames(d) <- NULL
  d
}

#' Calibrate radiocarbon ages with IntCal20
#'
#' The probability of every calendar age (5-year grid) given a radiocarbon
#' age and its error, using the IntCal20 curve (package `rintcal`, or a
#' curve you supply); returns the median and the 95% range.
#'
#' @param c14 Radiocarbon ages (14C years BP).
#' @param err Their 1-sigma laboratory errors.
#' @param curve Optional calibration curve: a matrix with calendar age BP,
#'   14C age and its 1-sigma error (default `rintcal::ccurve(1)`, IntCal20).
#' @return A data.frame with `age` (median, cal BP), `older` and `younger`
#'   (95% range). Ages of 0 or below (post-bomb) are returned unchanged.
#' @export
#' @examples
#' curve <- cbind(0:20000, 0:20000 * 0.9, 20)   # a toy curve
#' calibrate_radiocarbon(c(9000, 10800), c(50, 100), curve = curve)
calibrate_radiocarbon <- function(c14, err, curve = NULL) {
  if (is.null(curve)) {
    if (!requireNamespace("rintcal", quietly = TRUE)) {
      stop("Install rintcal or give a `curve`.", call. = FALSE)
    }
    curve <- rintcal::ccurve(1)
  }
  cc <- as.matrix(curve)[, 1:3]
  grid <- seq(0, max(cc[, 1]), by = 5)
  mu <- stats::approx(cc[, 1], cc[, 2], grid)$y
  sg <- stats::approx(cc[, 1], cc[, 3], grid)$y
  err <- rep_len(err, length(c14))
  res <- t(vapply(seq_along(c14), function(i) {
    if (is.na(c14[i])) return(c(NA_real_, NA_real_, NA_real_))
    if (c14[i] <= 0) return(c(c14[i], NA_real_, NA_real_))
    p <- stats::dnorm(c14[i], mu, sqrt(err[i]^2 + sg^2))
    p[!is.finite(p)] <- 0
    if (sum(p) == 0) return(c(NA_real_, NA_real_, NA_real_))
    cp <- cumsum(p) / sum(p)
    grid[c(which(cp >= 0.5)[1], which(cp >= 0.975)[1],
           which(cp >= 0.025)[1])]
  }, numeric(3)))
  data.frame(age = res[, 1], older = res[, 2], younger = res[, 3])
}

#' CHELSA-TraCE21k time-slice ID of an age
#'
#' The century IDs of CHELSA-TraCE21k: ID 20 = 1900-1990 CE (0 ka BP),
#' ID 0 = 100-1 BCE (2 ka BP), ID -200 = 22 ka BP; i.e.
#' `ID = round((2000 - age) / 100)`, limited to -200..20.
#'
#' @param age_bp Ages in calendar years BP.
#' @return Integer time-slice IDs.
#' @export
#' @examples
#' trace21k_id(c(0, 2000, 12000, 30000))
trace21k_id <- function(age_bp) {
  as.integer(pmax(-200, pmin(20, round((2000 - age_bp) / 100))))
}

#' One fossil record per site and time bin
#'
#' Summarises usable samples (`FLAG == "ok"`) into site x time-bin records
#' (e.g. per millennium), so that one long core does not dominate a
#' validation score. Each record keeps its age window and time-slice range;
#' expand them to one row per slice with [site_bin_slices()] for
#' [temporal_validation()].
#'
#' @param records Output of [neotoma_records()].
#' @param bin_years Width of the time bins (years).
#' @param periods Upper age limits of named periods (see [age_to_period()]).
#' @return A data.frame with `RECORD` (site_bin id), `SITEID`, `SITENAME`,
#'   `LAT`, `LON`, `BIN`, `N_SAMPLES`, `RECORD_TYPE` (macrofossil if any
#'   sample is), `TAXA`, `AGE_MEAN`, `AGE_MIN`, `AGE_MAX`, `SLICE_MIN`,
#'   `SLICE_MAX`, `MAX_VALUE` and `PERIOD`.
#' @export
#' @examples
#' smp <- data.frame(siteid = "1", sitename = "A", lat = 47, long = 11,
#'                   variablename = "Dryas", value = c(2, 5, 1),
#'                   age = c(12100, 12600, 15400),
#'                   agetype = "Calendar years BP")
#' fossil_site_bins(neotoma_records(smp, "^Dryas", calibrate_c14 = FALSE))
fossil_site_bins <- function(records, bin_years = 1000,
                             periods = c(Holocene = 11700,
                                         Lateglacial = 14700,
                                         Glacial = Inf)) {
  ok <- records[records$FLAG %in% "ok", , drop = FALSE]
  if (!nrow(ok)) stop("No usable records (FLAG == \"ok\").", call. = FALSE)
  ok$BIN <- floor(ok$AGE / bin_years) * bin_years
  key <- paste(ok$SITEID, ok$BIN, sep = "_")
  out <- do.call(rbind, lapply(split(seq_len(nrow(ok)), key), function(i) {
    x <- ok[i, , drop = FALSE]
    data.frame(RECORD = paste(x$SITEID[1], x$BIN[1], sep = "_"),
               SITEID = x$SITEID[1], SITENAME = x$SITENAME[1],
               LAT = x$LAT[1], LON = x$LON[1], BIN = x$BIN[1],
               N_SAMPLES = nrow(x),
               RECORD_TYPE = if (any(x$RECORD_TYPE == "macrofossil"))
                 "macrofossil" else "pollen",
               TAXA = paste(sort(unique(x$TAXON)), collapse = "; "),
               AGE_MEAN = mean(x$AGE), AGE_MIN = min(x$AGE_MIN),
               AGE_MAX = max(x$AGE_MAX), SLICE_MIN = min(x$SLICE_MIN),
               SLICE_MAX = max(x$SLICE_MAX), MAX_VALUE = max(x$VALUE),
               stringsAsFactors = FALSE)
  }))
  out$PERIOD <- age_to_period(out$AGE_MEAN, periods)
  out <- out[order(-out$BIN, out$SITEID), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Expand fossil site-bins to one row per climate time slice
#'
#' Each site-bin covers the time slices from `SLICE_MIN` to `SLICE_MAX`; the
#' result has one row per record and slice, ready to receive model
#' predictions under each slice's climate for [temporal_validation()].
#'
#' @param site_bins Output of [fossil_site_bins()].
#' @return A data.frame with `record`, `slice`, `lon`, `lat` and `period`.
#' @export
#' @examples
#' sb <- data.frame(RECORD = "s1_12000", LON = 11, LAT = 47, SLICE_MIN = -105,
#'                  SLICE_MAX = -102, PERIOD = "Lateglacial")
#' site_bin_slices(sb)
site_bin_slices <- function(site_bins) {
  rows <- lapply(seq_len(nrow(site_bins)), function(i) {
    sl <- seq(site_bins$SLICE_MIN[i], site_bins$SLICE_MAX[i])
    data.frame(record = site_bins$RECORD[i], slice = sl,
               lon = site_bins$LON[i], lat = site_bins$LAT[i],
               period = if ("PERIOD" %in% names(site_bins))
                 site_bins$PERIOD[i] else NA_character_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
