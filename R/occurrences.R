# Reading and merging occurrence records ----------------------------------------

#' Read a GBIF occurrence download
#'
#' Reads the tab-delimited table of a GBIF "simple" or Darwin Core occurrence
#' download, either the unzipped `.csv` file or the downloaded `.zip`. All
#' columns are read as text, so that coordinates keep the decimals exactly as
#' written (needed by the `coarse_decimals` test of [filter_occurrences()]).
#'
#' GBIF tables contain stray quotes and tabs in free-text fields, so quoting is
#' switched off, as GBIF recommends. The fast `readr` package is used when it
#' is installed.
#'
#' @param file Path to the GBIF `.csv` (tab-delimited) or `.zip` file.
#' @return A data.frame with all columns as character.
#' @export
#' @examples
#' \dontrun{
#' occ <- read_gbif("0020328-260928105237408.zip")
#' }
read_gbif <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  if (grepl("\\.zip$", file, ignore.case = TRUE)) {
    inside <- utils::unzip(file, list = TRUE)$Name
    tab <- inside[grepl("\\.(csv|txt)$", inside) &
                    !grepl("verbatim|citations|rights|meta", inside)]
    tab <- if ("occurrence.txt" %in% tab) "occurrence.txt" else tab[1]
    if (is.na(tab)) stop("No occurrence table found in ", file, call. = FALSE)
    td <- tempfile("gbif")
    dir.create(td)
    on.exit(unlink(td, recursive = TRUE), add = TRUE)
    utils::unzip(file, files = tab, exdir = td)
    file <- file.path(td, tab)
  }
  if (requireNamespace("readr", quietly = TRUE)) {
    x <- readr::read_tsv(file, col_types = readr::cols(.default = "c"),
                         quote = "", na = c("", "NA"),
                         show_col_types = FALSE, progress = FALSE)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  } else {
    x <- utils::read.delim(file, quote = "", colClasses = "character",
                           na.strings = c("", "NA"), comment.char = "",
                           check.names = FALSE, fill = TRUE)
  }
  x
}

#' Merge occurrence records from several sources
#'
#' Combines occurrence tables that use different column names (e.g. a GBIF
#' download and your own field or genetic samples) into one table with the
#' standard columns `species`, `lon`, `lat` and `source`, followed by all
#' original columns (as text). Coordinate and species columns are detected
#' automatically; records without valid coordinates are dropped.
#'
#' The order of the inputs matters for [thin_by_cell()]: when two records fall
#' in the same raster cell, the record from the source listed first is kept
#' (use `priority` there to change this).
#'
#' @param ... Named data.frames or paths to CSV files (comma-separated, or
#'   semicolon-separated with decimal commas), e.g. `merge_occurrences(sNMF = "samples.csv",
#'   GBIF = cleaned)`. Names become the `source` column.
#' @param species Optional single species name written to every record. If
#'   `NULL`, the species column of each table is used when found.
#' @param encoding Encoding assumed for text that is not valid UTF-8 (common
#'   in CSV files exported from Excel or some databases); such text is
#'   converted to UTF-8. Try `"CP1250"` for Central European files if the
#'   repaired characters look wrong.
#' @return A data.frame with `species`, `lon`, `lat`, `source` and the other
#'   columns of the inputs.
#' @export
#' @examples
#' a <- data.frame(Longitude = c(10.1, 10.2), Latitude = c(46.1, 46.2))
#' b <- data.frame(decimalLongitude = c(10.1, 11), decimalLatitude = c(46.1, 47),
#'                 species = "Dryas octopetala", gbifID = c("1", "2"))
#' merge_occurrences(own = a, GBIF = b, species = "Dryas octopetala")
merge_occurrences <- function(..., species = NULL, encoding = "latin1") {
  inputs <- list(...)
  if (!length(inputs) || is.null(names(inputs)) || any(names(inputs) == "")) {
    stop("Give every input a name, e.g. merge_occurrences(own = a, GBIF = b).",
         call. = FALSE)
  }
  parts <- lapply(names(inputs), function(src) {
    x <- inputs[[src]]
    if (is.character(x) && length(x) == 1) x <- read_any_csv(x)
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    xy <- guess_coords(names(x))
    sp <- find_column(names(x), c("species", "scientificname", "taxon"))
    x[] <- lapply(x, as.character)
    names(x) <- fix_utf8(names(x), encoding)
    n_fixed <- 0L
    for (j in seq_along(x)) {
      v <- fix_utf8(x[[j]], encoding)
      n_fixed <- n_fixed + attr(v, "n_fixed")
      attr(v, "n_fixed") <- NULL
      x[[j]] <- v
    }
    if (n_fixed) {
      message(src, ": ", n_fixed, " text value(s) with invalid UTF-8 ",
              "repaired (read as ", encoding, ")")
    }
    out <- data.frame(
      species = if (!is.null(species)) species else if (!is.na(sp)) x[[sp]]
      else NA_character_,
      lon = as_number(x[[xy[1]]]),
      lat = as_number(x[[xy[2]]]),
      source = src, stringsAsFactors = FALSE)
    rest <- x[setdiff(names(x), c("species", "lon", "lat", "source"))]
    cbind(out, rest)
  })
  all_cols <- unique(unlist(lapply(parts, names)))
  parts <- lapply(parts, function(p) {
    p[setdiff(all_cols, names(p))] <- NA_character_
    p[all_cols]
  })
  out <- do.call(rbind, parts)
  out$species[is.na(out$species)] <- "unknown"
  bad <- is.na(out$lon) | is.na(out$lat) | abs(out$lon) > 180 |
    abs(out$lat) > 90
  if (any(bad)) {
    message("Records without valid coordinates (dropped): ", sum(bad))
  }
  out <- out[!bad, , drop = FALSE]
  rownames(out) <- NULL
  message("Merged records: ", paste(names(table(out$source)), table(out$source),
                                    sep = " = ", collapse = ", "))
  out
}

# Read a CSV that may be comma- or semicolon-separated (internal)
read_any_csv <- function(f) {
  if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
  first <- readLines(f, n = 1, warn = FALSE)
  semi <- lengths(regmatches(first, gregexpr(";", first))) >
    lengths(regmatches(first, gregexpr(",", first)))
  if (semi) {
    utils::read.csv2(f, check.names = FALSE, colClasses = "character",
                     na.strings = c("", "NA"))
  } else {
    utils::read.csv(f, check.names = FALSE, colClasses = "character",
                    na.strings = c("", "NA"))
  }
}

# Numbers written with a decimal point or a decimal comma (internal)
as_number <- function(v) suppressWarnings(as.numeric(sub(",", ".", v,
                                                          fixed = TRUE)))

# Make text valid UTF-8, re-reading invalid strings in `from` (internal)
fix_utf8 <- function(x, from = "latin1") {
  if (!is.character(x)) return(x)
  Encoding(x) <- "unknown"
  bad <- !is.na(x) & !validUTF8(x)
  if (any(bad)) {
    y <- iconv(x[bad], from = from, to = "UTF-8")
    y[is.na(y)] <- iconv(x[bad][is.na(y)], "UTF-8", "UTF-8", sub = "?")
    x[bad] <- y
  }
  Encoding(x) <- "UTF-8"
  attr(x, "n_fixed") <- sum(bad)
  x
}

# First matching column name, case-insensitive (internal)
find_column <- function(nms, candidates) {
  i <- match(tolower(candidates), tolower(nms))
  i <- i[!is.na(i)]
  if (length(i)) nms[i[1]] else NA_character_
}
