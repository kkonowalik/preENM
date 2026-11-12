# Working with filtering results ------------------------------------------------

#' @export
print.occ_filter <- function(x, ...) {
  d <- x$data
  n <- nrow(d)
  k <- sum(is.na(d$reason))
  cat("Occurrence filtering (preset: ", x$settings$preset, "): ", n,
      " records -> ", k, " kept (", round(100 * k / max(n, 1), 1), "%)\n\n",
      sep = "")
  lg <- x$log
  run <- lg[lg$status == "run" & lg$removed > 0, c("reason", "removed")]
  if (nrow(run)) {
    cat("Removed, by the first test failed (in test order):\n")
    print(run, row.names = FALSE)
  }
  sk <- lg[lg$status == "skipped", c("test", "note")]
  if (nrow(sk)) {
    cat("\nSkipped tests:\n")
    print(sk, row.names = FALSE, right = FALSE)
  }
  ov <- sum(!is.na(d$reason_overridden))
  if (ov) cat("\nKept by hand despite failing a test:", ov, "\n")
  invisible(x)
}

#' Kept or removed records of a filtering or thinning run
#'
#' @param x An `occ_filter` object from [filter_occurrences()] or an
#'   `env_thinning` from [thin_environment()].
#' @return A data.frame: the original columns plus the diagnostic columns
#'   added by the filtering (`reason`, `climate`, `built_pct`, ...) or the
#'   thinning (`env_*`, `ax_*`, `bin`, `block`, `keep_reason`).
#' @export
#' @examples
#' occ <- data.frame(decimalLongitude = c(10.123, 10.123, 0),
#'                   decimalLatitude = c(46.123, 46.123, 0))
#' res <- filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
#'                            verbose = FALSE)
#' kept(res)
#' removed(res)
kept <- function(x) UseMethod("kept")

#' @rdname kept
#' @export
removed <- function(x) UseMethod("removed")

#' @export
kept.default <- function(x) {
  stop("`x` must come from filter_occurrences() or thin_environment().",
       call. = FALSE)
}

#' @export
removed.default <- kept.default

#' @rdname kept
#' @export
kept.occ_filter <- function(x) {
  d <- x$data[is.na(x$data$reason), , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' @rdname kept
#' @export
removed.occ_filter <- function(x) {
  d <- x$data[!is.na(x$data$reason), , drop = FALSE]
  rownames(d) <- NULL
  d
}

check_filter <- function(x) {
  if (!inherits(x, "occ_filter")) {
    stop("`x` must be the output of filter_occurrences().", call. = FALSE)
  }
}

#' Summary tables of a filtering run
#'
#' Counts records per outcome (the reason of removal, or `KEPT`), optionally
#' split by a grouping column such as `countryCode`. Use it to see which tests
#' remove records where, e.g. whether one country loses far more than others.
#'
#' @param x An `occ_filter` object.
#' @param by Optional column name to split by (e.g. `"countryCode"`).
#' @return A data.frame: with `by = NULL`, columns `outcome` and `n`;
#'   otherwise one row per group and one column per outcome.
#' @export
#' @examples
#' occ <- data.frame(decimalLongitude = c(10.123, 10.123, 0, 11.456),
#'                   decimalLatitude = c(46.123, 46.123, 0, 46.789),
#'                   countryCode = c("AT", "AT", "IT", "IT"))
#' res <- filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
#'                            verbose = FALSE)
#' filter_report(res)
#' filter_report(res, by = "countryCode")
filter_report <- function(x, by = NULL) {
  check_filter(x)
  d <- x$data
  outcome <- ifelse(is.na(d$reason), "KEPT", d$reason)
  lev <- c("KEPT", intersect(x$log$reason, unique(outcome)))
  outcome <- factor(outcome, levels = unique(c(lev, outcome)))
  if (is.null(by)) {
    t <- table(outcome)
    return(data.frame(outcome = names(t), n = as.integer(t),
                      stringsAsFactors = FALSE))
  }
  if (!by %in% names(d)) stop("Column not found: ", by, call. = FALSE)
  t <- table(d[[by]], outcome, useNA = "ifany")
  grp <- rownames(t)
  grp[is.na(grp)] <- "NA"
  out <- as.data.frame(unclass(t), check.names = FALSE)
  out <- cbind(stats::setNames(data.frame(grp, stringsAsFactors = FALSE), by),
               out)
  rownames(out) <- NULL
  out
}

#' How many records would a threshold remove?
#'
#' A sensitivity table to choose thresholds: for each candidate threshold,
#' the number of records it would remove among those that passed every other
#' test (records kept, plus those removed only by the tests that use this
#' column).
#'
#' @param x An `occ_filter` object.
#' @param column Diagnostic column, e.g. `"climate"`, `"built_pct"`,
#'   `"crop_pct"`, `"ghs_mean_pct"`, `"ghs_max_m2"`, `"dist_land_m"`.
#' @param thresholds Numeric vector of thresholds to try.
#' @param direction `"above"`: values at or above the threshold would be
#'   removed (default); `"below"`: values at or below.
#' @param by Optional grouping column (e.g. `"countryCode"`): counts are
#'   then also given per group.
#' @param reasons Removal reasons whose records are also counted (default:
#'   the tests that use `column`).
#' @return A data.frame with `threshold`, `removed` and one column per group.
#' @export
#' @examples
#' occ <- data.frame(decimalLongitude = c(10, 11, 12), decimalLatitude = 46,
#'                   built = c(0.5, 4, 20))
#' res <- filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
#'                            verbose = FALSE)
#' res$data$built_pct <- res$data$built     # normally added by a built layer
#' filter_sensitivity(res, "built_pct", c(1, 5, 10))
filter_sensitivity <- function(x, column, thresholds, direction = c("above",
                                                                    "below"),
                               by = NULL, reasons = NULL) {
  check_filter(x)
  direction <- match.arg(direction)
  d <- x$data
  if (!column %in% names(d)) {
    stop("Column not found: ", column, ". Was its layer supplied?",
         call. = FALSE)
  }
  if (is.null(reasons)) {
    reasons <- switch(column,
      climate = c("climate_limit", "climate_soft_limit",
                  "isolated_climate_outlier", "no_climate_data"),
      built_pct = "built_up", crop_pct = "cropland",
      ghs_mean_pct = "ghs_dense", ghs_max_m2 = "ghs_building_nearby",
      dist_land_m = "at_sea", land_near_pct = "on_water", character(0))
  }
  pool <- is.na(d$reason) | d$reason %in% reasons
  v <- d[[column]]
  rows <- lapply(thresholds, function(th) {
    hit <- pool & !is.na(v) & (if (direction == "above") v >= th else v <= th)
    r <- data.frame(threshold = th, removed = sum(hit))
    if (!is.null(by)) {
      t <- table(factor(d[[by]][hit], levels = sort(unique(d[[by]][pool]))))
      r <- cbind(r, as.data.frame(as.list(t), check.names = FALSE))
    }
    r
  })
  do.call(rbind, rows)
}

#' Fate of selected records
#'
#' Shows whether chosen records were kept or which test removed them, with
#' their diagnostic values. Useful to check records you know well.
#'
#' @param x An `occ_filter` object.
#' @param ids Record IDs (values of the id column, GBIF `gbifID`).
#' @return A data.frame with the id, `fate` (`KEPT` or the reason) and the
#'   diagnostic columns; IDs not found are reported in a message (they may
#'   have been merged as duplicates).
#' @export
#' @examples
#' occ <- data.frame(gbifID = c("1", "2"), decimalLongitude = c(10.123, 0),
#'                   decimalLatitude = c(46.123, 0))
#' record_fate(filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
#'                            verbose = FALSE), c("1", "2", "9"))
record_fate <- function(x, ids) {
  check_filter(x)
  idc <- x$columns[["id"]]
  if (is.na(idc)) stop("No id column in the data.", call. = FALSE)
  d <- x$data
  i <- which(as.character(d[[idc]]) %in% as.character(ids))
  diag_cols <- intersect(c("reason_overridden", "climate", "climate_near",
                           "built_pct", "crop_pct", "ghs_mean_pct",
                           "ghs_max_m2", "land_near_pct", "dist_land_m",
                           "grid_centre", "in_range"), names(d))
  out <- data.frame(id = d[[idc]][i],
                    fate = ifelse(is.na(d$reason[i]), "KEPT", d$reason[i]),
                    d[i, diag_cols, drop = FALSE], stringsAsFactors = FALSE)
  miss <- setdiff(as.character(ids), as.character(d[[idc]]))
  if (length(miss)) message("Not found: ", paste(miss, collapse = ", "))
  rownames(out) <- NULL
  out
}

#' Write the results of a filtering run to CSV files
#'
#' Writes `<prefix>_kept.csv` (records that passed every test),
#' `<prefix>_removed.csv` (the rest, with the `reason`),
#' `<prefix>_log.csv` (records removed per test and skipped tests) and, for
#' GBIF data, `<prefix>_datasets.csv` (gridded-data diagnostics per dataset).
#'
#' @param x An `occ_filter` object.
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix, e.g. the species and preset.
#' @return The paths of the written files, invisibly.
#' @export
#' @examples
#' occ <- data.frame(decimalLongitude = c(10.123, 0), decimalLatitude = c(46.1, 0))
#' res <- filter_occurrences(occ, filter_settings(cc_tests = "zeros"),
#'                            verbose = FALSE)
#' write_filter_results(res, tempdir(), "example")
write_filter_results <- function(x, dir, prefix = "occurrences") {
  check_filter(x)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, "_", s, ".csv"))
  paths <- c(kept = f("kept"), removed = f("removed"), log = f("log"))
  utils::write.csv(kept(x), paths[["kept"]], row.names = FALSE)
  utils::write.csv(removed(x), paths[["removed"]], row.names = FALSE)
  utils::write.csv(x$log, paths[["log"]], row.names = FALSE)
  if (!is.null(x$datasets)) {
    paths[["datasets"]] <- f("datasets")
    utils::write.csv(x$datasets, paths[["datasets"]], row.names = FALSE)
  }
  invisible(paths)
}
