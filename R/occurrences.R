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

# Taxonomic name matching ------------------------------------------------------

#' Match the scientific names of occurrence records to a taxonomic backbone
#'
#' [filter_occurrences()] tests many things about a record, but not whether
#' its name is right, and different cleaning pipelines lead to different
#' downstream conclusions (Fuehrding-Potschkat et al. 2022; Van Eupen et al.
#' 2021). Records of a species in GBIF arrive under synonyms, outdated
#' names, misspellings, genus-only or aggregate names, hybrid formulas and
#' infraspecific names; checking them against a backbone and keeping only
#' records that can be assigned confidently to the accepted taxon is part of
#' a transparent cleaning workflow (Lipinski 2026). `match_names()` resolves
#' every distinct name once and adds the result to each record; the name
#' tests of [filter_occurrences()] (`name_unmatched`,
#' `name_synonym_conflict`, `rank_conflict`) then use it.
#'
#' **Backbones**
#' * `"gbif"` (default): the GBIF name matching service through the rgbif
#'   package (`rgbif::name_backbone_checklist()`); with rgbif 3.9 or later
#'   another checklist can be chosen with `checklist` (a GBIF checklist
#'   key).
#' * `"checklistbank"`: a dataset of ChecklistBank (default `"3LR"`, the
#'   latest Catalogue of Life release) through its name-match API; needs
#'   the jsonlite package.
#' * A Darwin Core taxon table (data.frame or path): e.g. the World Flora
#'   Online backbone (`classification.csv`, from
#'   `WorldFlora::WFO.download()`), a GBIF backbone `Taxon.tsv`, or your
#'   own checklist. Matched locally, without internet.
#'
#' **What is returned for each record**
#' * `name_query`: the name looked up (from `name_col`).
#' * `name_matched`: the backbone name it matched, without authorship.
#' * `name_accepted`, `name_accepted_species`: the accepted name (a synonym
#'   is followed to it) and its species (for infraspecific taxa).
#' * `name_rank`, `name_status`: rank and status of the matched name.
#' * `name_match_type`: `"exact"`, `"fuzzy"` (the matched name is spelled
#'   differently from the query), `"higherrank"` (only a higher taxon, e.g.
#'   the genus, was found) or `"none"`.
#' * `name_confidence`: the backbone's confidence (0-100; GBIF), or, for a
#'   local table, 100 minus the share of edited characters (`NA` for
#'   ChecklistBank, which does not report one).
#' * Flags: `name_synonym`, `name_hybrid` (a hybrid formula or nothotaxon,
#'   from the name), `name_aggregate` (an aggregate, *agg.* or *s.l.*),
#'   `name_higherrank`, `name_fuzzy` and `name_unmatched`.
#'
#' **Cache.** Names are looked up once per distinct name; results are kept
#' for the R session and, with `cache_file`, in an RDS file, so that
#' repeated runs over a large download do not query the service again.
#'
#' @param x Occurrence records: a data.frame (e.g. from [read_gbif()] or
#'   [merge_occurrences()]) or a path to a GBIF download.
#' @param name_col Column with the names. Default: the first of
#'   `verbatimScientificName` (the name as recorded, in GBIF downloads),
#'   `scientificName`, `species` and `taxon` that exists. In GBIF downloads
#'   `scientificName` and `species` were already matched by GBIF; matching
#'   the verbatim name shows what that interpretation did.
#' @param backbone `"gbif"`, `"checklistbank"` or a Darwin Core taxon table
#'   (data.frame or file) with `taxonID`, `scientificName`, `taxonRank`,
#'   `taxonomicStatus` and `acceptedNameUsageID` (and optionally
#'   `parentNameUsageID`, `genus`, `specificEpithet`,
#'   `infraspecificEpithet`).
#' @param kingdom Optional kingdom, which avoids homonyms between plants and
#'   animals (e.g. `"Plantae"`); GBIF only.
#' @param checklist A GBIF checklist key (rgbif 3.9 or later) or a
#'   ChecklistBank dataset key; `NULL` for the default of each service.
#' @param max_dist Largest number of edited characters for a fuzzy match to
#'   a local table (at most a tenth of the name is also required).
#' @param cache Keep the results of online look-ups for the session.
#' @param cache_file Optional RDS file to keep them across sessions.
#' @param verbose Print a summary.
#' @return `x` with the `name_*` columns added. `attr(, "names_matched")`
#'   holds one row per distinct name.
#' @references Fuehrding-Potschkat, P., Kreft, H., Ickert-Bond, S.M. (2022)
#'   Influence of different data cleaning solutions of point-occurrence
#'   records on downstream macroecological diversity models. *Ecology and
#'   Evolution* 12, e9168. \doi{10.1002/ece3.9168}
#'
#'   Van Eupen, C., Maes, D., Herremans, M., Swinnen, K.R.R., Somers, B.,
#'   Luca, S. (2021) The impact of data quality filtering of opportunistic
#'   citizen science data on species distribution model performance.
#'   *Ecological Modelling* 444, 109453.
#'   \doi{10.1016/j.ecolmodel.2021.109453}
#'
#'   Lipinski, S. (2026) Spatial bias in open biodiversity data: how GBIF
#'   record quality shapes conservation analyses. *Conservation* 6, 40.
#'   \doi{10.3390/conservation6020040}
#' @seealso [filter_occurrences()], [filter_settings()], [taxon_complexity()]
#' @export
#' @examples
#' # a small local backbone in Darwin Core form
#' bb <- data.frame(
#'   taxonID = c("1", "2", "3", "4", "5"),
#'   scientificName = c("Aus bus L.", "Aus cus Mill.", "Aus dus Sm.",
#'                      "Aus bus subsp. eus Hook.", "Aus"),
#'   taxonRank = c("species", "species", "species", "subspecies", "genus"),
#'   taxonomicStatus = c("accepted", "accepted", "synonym", "accepted",
#'                       "accepted"),
#'   acceptedNameUsageID = c("", "", "1", "", ""),
#'   parentNameUsageID = c("5", "5", "", "1", ""))
#' occ <- data.frame(lon = 1:5, lat = 50,
#'                   scientificName = c("Aus bus", "Aus dus", "Aus buss",
#'                                      "Aus sp.", "Aus bus subsp. eus"))
#' m <- match_names(occ, backbone = bb)
#' m[, c("name_query", "name_accepted", "name_match_type", "name_rank")]
#' \donttest{
#' # GBIF (needs rgbif and internet)
#' if (requireNamespace("rgbif", quietly = TRUE)) {
#'   match_names(data.frame(species = "Potentilla tabernaemontani"),
#'               kingdom = "Plantae")
#' }
#' }
match_names <- function(x, name_col = NULL, backbone = "gbif",
                        kingdom = NULL, checklist = NULL, max_dist = 2,
                        cache = TRUE, cache_file = NULL, verbose = TRUE) {
  if (is.character(x) && length(x) == 1) x <- read_gbif(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (is.null(name_col)) {
    name_col <- find_column(names(x), c("verbatimScientificName",
                                        "scientificName", "species",
                                        "taxon"))
    if (is.na(name_col)) {
      stop("No name column found: give it in `name_col`.", call. = FALSE)
    }
  } else if (!name_col %in% names(x)) {
    stop("Column not found: ", name_col, call. = FALSE)
  }
  q <- trimws(gsub("\\s+", " ", as.character(x[[name_col]])))
  q[q == ""] <- NA
  u <- unique(q[!is.na(q)])
  local <- is.data.frame(backbone) ||
    (is.character(backbone) && length(backbone) == 1 &&
       !backbone %in% c("gbif", "checklistbank"))
  service <- if (local) "local" else backbone
  res <- if (!length(u)) {
    match_rows(character(0))
  } else if (local) {
    match_local(u, read_backbone(backbone), max_dist)
  } else {
    key <- paste(service, checklist, kingdom, u, sep = "|")
    store <- name_cache(cache_file)
    have <- if (cache) key %in% names(store) else rep(FALSE, length(u))
    out <- vector("list", length(u))
    if (any(have)) out[have] <- mget(key[have], envir = store)
    if (any(!have)) {
      fresh <- if (service == "gbif") {
        match_gbif(u[!have], kingdom, checklist)
      } else {
        match_clb(u[!have], if (is.null(checklist)) "3LR" else checklist)
      }
      rows <- split(fresh, seq_len(nrow(fresh)))
      out[!have] <- rows
      if (cache) {
        for (i in seq_along(rows)) assign(key[!have][i], rows[[i]],
                                          envir = store)
        if (!is.null(cache_file)) {
          saveRDS(as.list(store), cache_file)
        }
      }
    }
    do.call(rbind, out)
  }
  res$name_query <- u
  res <- finish_matches(res)
  rownames(res) <- NULL
  i <- match(q, res$name_query)
  add <- res[i, setdiff(names(res), "name_query"), drop = FALSE]
  # records without a name: not matched
  if (any(is.na(i))) {
    add$name_match_type[is.na(i)] <- "none"
    for (k in c("name_synonym", "name_hybrid", "name_aggregate",
                "name_higherrank", "name_fuzzy")) add[[k]][is.na(i)] <- FALSE
    add$name_unmatched[is.na(i)] <- TRUE
  }
  out <- cbind(x, name_query = q, add, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "names_matched") <- res
  attr(out, "backbone") <- service
  if (verbose) say_matches(out, res, service)
  out
}

# Package-level cache of online name look-ups (internal)
.name_cache <- new.env(parent = emptyenv())
name_cache <- function(file) {
  if (!is.null(file) && file.exists(file) &&
      !isTRUE(attr(.name_cache, "loaded_from") == file)) {
    old <- readRDS(file)
    for (k in names(old)) assign(k, old[[k]], envir = .name_cache)
    attr(.name_cache, "loaded_from") <- file
  }
  .name_cache
}

# Empty or filled result rows (internal)
match_rows <- function(n_or_names, ...) {
  n <- if (is.numeric(n_or_names)) n_or_names else length(n_or_names)
  d <- data.frame(name_matched = rep(NA_character_, n),
                  name_accepted = NA_character_,
                  name_accepted_species = NA_character_,
                  name_rank = NA_character_, name_status = NA_character_,
                  name_match_type = "none", name_confidence = NA_real_,
                  stringsAsFactors = FALSE)
  v <- list(...)
  for (k in names(v)) d[[k]] <- v[[k]]
  d
}

# Canonical form of scientific names: genus, epithet and infraspecific
# epithet without authorship and rank markers, with hybrid and aggregate
# flags (internal)
parse_names <- function(x) {
  markers <- c("subsp.", "subsp", "ssp.", "ssp", "var.", "var", "f.", "fo.",
               "forma", "subvar.", "subf.", "nothosubsp.", "nothovar.",
               "nothossp.", "cv.", "convar.", "morph", "race")
  agg <- c("agg.", "agg", "aggr.", "aggr", "s.l.", "s.lat.", "sensu",
           "lato", "s.str.", "s.s.", "coll.", "group", "complex")
  stop_at <- c("sp.", "sp", "spp.", "spp", "spec.", "indet.", "cf.", "aff.",
               "?")
  particles <- c("ex", "de", "von", "van", "in", "et", "al.", "la", "le", "da",
                 "du", "der", "den", "d'", "del", "della", "fil.", "f", "and",
                 "non", "emend.", "nom.", "illeg.", "nud.")
  one <- function(s) {
    if (is.na(s) || !nzchar(s)) return(c(NA_character_, "FALSE", "FALSE"))
    s <- gsub("\u00d7", " x ", s, fixed = TRUE, useBytes = TRUE)
    tk <- strsplit(trimws(gsub("\\s+", " ", s)), " ", fixed = TRUE)[[1]]
    hyb <- FALSE
    ag <- FALSE
    if (length(tk) && tolower(tk[1]) == "x") {
      hyb <- TRUE
      tk <- tk[-1]
    }
    if (!length(tk)) return(c(NA_character_, as.character(hyb), as.character(ag)))
    out <- sub("^([A-Za-z-]+).*$", "\\1", tk[1])
    out <- paste0(toupper(substr(out, 1, 1)), tolower(substring(out, 2)))
    for (t in tk[-1]) {
      tl <- tolower(t)
      if (tl == "x") {
        hyb <- TRUE
        if (length(out) >= 2) break
        next
      }
      if (tl %in% agg) {
        ag <- TRUE
        next
      }
      if (tl %in% stop_at) break
      if (tl %in% markers) next
      if (grepl("^[a-z][a-z-]+$", t) && !tl %in% particles &&
          length(out) < 3) {
        out <- c(out, t)
      }
    }
    c(paste(out, collapse = " "), as.character(hyb), as.character(ag))
  }
  m <- vapply(x, one, character(3), USE.NAMES = FALSE)
  data.frame(canonical = m[1, ], hybrid = m[2, ] == "TRUE",
             aggregate = m[3, ] == "TRUE", stringsAsFactors = FALSE)
}

# Harmonised rank names (internal)
norm_rank <- function(r) {
  r <- tolower(as.character(r))
  r[r %in% c("ssp", "subsp", "subsp.")] <- "subspecies"
  r[r %in% c("var", "var.")] <- "variety"
  r[r %in% c("f", "f.", "forma")] <- "form"
  r
}

# Order of ranks, coarse to fine (internal)
rank_level <- function(r) {
  lv <- c(domain = 1, kingdom = 2, subkingdom = 2.5, phylum = 3,
          subphylum = 3.5, class = 4, subclass = 4.5, order = 5,
          suborder = 5.5, family = 6, subfamily = 6.3, tribe = 6.5,
          subtribe = 6.7, genus = 7, subgenus = 7.3, section = 7.5,
          series = 7.7, species = 8, subspecies = 9, variety = 9.5,
          subvariety = 9.7, form = 10, subform = 10.5, infraspecies = 9,
          infraspecific_name = 9, unranked = NA)
  unname(lv[norm_rank(r)])
}

# Flags and harmonised fields shared by all backbones (internal)
finish_matches <- function(d) {
  p <- parse_names(d$name_query)
  d$name_rank <- norm_rank(d$name_rank)
  d$name_status <- tolower(d$name_status)
  d$name_hybrid <- p$hybrid
  d$name_aggregate <- p$aggregate
  mt <- tolower(d$name_match_type)
  mt[mt %in% c("variant", "fuzzy", "canonical")] <- "fuzzy"
  mt[is.na(mt) | mt %in% c("ambiguous", "none", "")] <- "none"
  # fuzzy only when the matched name is spelled differently from the
  # query (services also report authorship variants as variants)
  same <- !is.na(d$name_matched) & !is.na(p$canonical) &
    tolower(parse_names(d$name_matched)$canonical) == tolower(p$canonical)
  mt[mt == "fuzzy" & same] <- "exact"
  mt[mt == "exact" & !same & !is.na(d$name_matched) &
       rank_level(d$name_rank) >= 8 & !is.na(p$canonical)] <- "fuzzy"
  d$name_match_type <- mt
  d$name_synonym <- !is.na(d$name_status) & d$name_status != "accepted" &
    !is.na(d$name_accepted) & d$name_accepted != d$name_matched
  d$name_higherrank <- mt == "higherrank"
  d$name_fuzzy <- mt == "fuzzy"
  d$name_unmatched <- mt == "none"
  for (k in c("name_matched", "name_accepted", "name_accepted_species")) {
    v <- d[[k]]
    ok <- !is.na(v)
    v[ok] <- parse_names(v[ok])$canonical
    d[[k]] <- v
  }
  d[c("name_query", "name_matched", "name_accepted", "name_accepted_species",
      "name_rank", "name_status", "name_match_type", "name_confidence",
      "name_synonym", "name_hybrid", "name_aggregate", "name_higherrank",
      "name_fuzzy", "name_unmatched")]
}

# GBIF name matching via rgbif (internal)
match_gbif <- function(u, kingdom, checklist) {
  if (!requireNamespace("rgbif", quietly = TRUE)) {
    stop("Matching against GBIF needs rgbif: install.packages(\"rgbif\")",
         call. = FALSE)
  }
  p <- parse_names(u)
  send <- ifelse(p$aggregate & !is.na(p$canonical), p$canonical, u)
  nd <- data.frame(scientificName = send, stringsAsFactors = FALSE)
  if (!is.null(kingdom)) nd$kingdom <- kingdom
  args <- list(name_data = nd)
  if (!is.null(checklist)) args$checklistKey <- checklist
  m <- as.data.frame(do.call(rgbif::name_backbone_checklist, args))
  if ("is_alternative" %in% names(m)) {
    m <- m[is.na(m$is_alternative) | !as.logical(m$is_alternative), ,
           drop = FALSE]
  }
  if (nrow(m) != length(u)) {
    m <- m[match(send, m$verbatim_name), , drop = FALSE]
  }
  col <- function(k) if (k %in% names(m)) as.character(m[[k]]) else
    rep(NA_character_, nrow(m))
  status <- col("status")
  rank <- col("rank")
  matched <- ifelse(is.na(col("canonicalName")), col("scientificName"),
                    col("canonicalName"))
  syn <- !is.na(status) & status != "ACCEPTED"
  acc <- ifelse(syn, col("acceptedScientificName"), matched)
  # older rgbif: the accepted name of a species-level synonym is its
  # species; otherwise look the accepted usage up
  need <- which(syn & is.na(acc))
  if (length(need)) {
    sp_rank <- toupper(rank[need]) == "SPECIES"
    acc[need[sp_rank]] <- col("species")[need[sp_rank]]
    for (i in need[!sp_rank]) {
      acc[i] <- tryCatch(rgbif::name_usage(
        key = col("acceptedUsageKey")[i])$data$canonicalName[1],
        error = function(e) NA_character_)
    }
  }
  fine <- rank_level(rank) >= 8
  match_rows(u, name_matched = matched, name_accepted = acc,
             name_accepted_species = ifelse(!is.na(fine) & fine,
                                            col("species"), NA),
             name_rank = rank, name_status = status,
             name_match_type = col("matchType"),
             name_confidence = suppressWarnings(as.numeric(
               col("confidence"))))
}

# ChecklistBank name matching (internal)
match_clb <- function(u, dataset) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Matching against ChecklistBank needs jsonlite: ",
         "install.packages(\"jsonlite\")", call. = FALSE)
  }
  p <- parse_names(u)
  send <- ifelse(p$aggregate & !is.na(p$canonical), p$canonical, u)
  rows <- lapply(send, function(s) {
    url <- paste0("https://api.checklistbank.org/dataset/",
                  utils::URLencode(dataset, reserved = TRUE),
                  "/match/nameusage?q=", utils::URLencode(s, reserved = TRUE))
    d <- tryCatch(jsonlite::fromJSON(url, simplifyVector = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || !isTRUE(d$match) || is.null(d$usage)) {
      return(match_rows(1))
    }
    us <- d$usage
    cl <- us$classification
    status <- us$status
    acc <- us$name
    if (!identical(status, "accepted") && length(cl)) acc <- cl[[1]]$name
    sp <- NA_character_
    for (c1 in cl) if (identical(c1$rank, "species")) {
      sp <- c1$name
      break
    }
    if (identical(us$rank, "species") && identical(status, "accepted")) {
      sp <- us$name
    }
    match_rows(1, name_matched = us$name, name_accepted = acc,
               name_accepted_species = sp, name_rank = us$rank,
               name_status = status, name_match_type = d$type)
  })
  do.call(rbind, rows)
}

# Read a Darwin Core backbone (internal)
read_backbone <- function(b) {
  if (is.character(b)) {
    if (!file.exists(b)) stop("Backbone file not found: ", b, call. = FALSE)
    b <- utils::read.delim(b, quote = "", encoding = "UTF-8",
                           colClasses = "character")
  }
  b <- as.data.frame(b, stringsAsFactors = FALSE)
  need <- c("taxonID", "scientificName", "taxonRank", "taxonomicStatus",
            "acceptedNameUsageID")
  miss <- setdiff(need, names(b))
  if (length(miss)) {
    stop("The backbone lacks the Darwin Core column(s) ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  }
  for (k in intersect(c(need, "parentNameUsageID", "genus",
                        "specificEpithet", "infraspecificEpithet"),
                      names(b))) {
    v <- as.character(b[[k]])
    v[is.na(v)] <- ""
    b[[k]] <- v
  }
  if (!"parentNameUsageID" %in% names(b)) b$parentNameUsageID <- ""
  # canonical names: from the epithet columns where present, else parsed
  canon <- if (all(c("genus", "specificEpithet") %in% names(b))) {
    trimws(gsub("\\s+", " ", paste(b$genus, b$specificEpithet,
                                   if ("infraspecificEpithet" %in% names(b))
                                     b$infraspecificEpithet else "")))
  } else rep("", nrow(b))
  bad <- !nzchar(canon) | !grepl(" ", canon) &
    rank_level(b$taxonRank) %in% c(8, 9, 9.5, 10)
  bad[is.na(bad)] <- TRUE
  if (any(bad)) canon[bad] <- parse_names(b$scientificName[bad])$canonical
  b$.canon <- tolower(canon)
  b$.status <- tolower(b$taxonomicStatus)
  b$.rank <- norm_rank(b$taxonRank)
  b
}

# Match names to a local Darwin Core table (internal)
match_local <- function(u, b, max_dist) {
  p <- parse_names(u)
  qc <- tolower(p$canonical)
  pref <- ifelse(b$.status == "accepted", 0, 1)
  ord <- order(pref)
  first_hit <- function(cands) {
    i <- which(b$.canon %in% cands)
    if (!length(i)) return(NA_integer_)
    i[order(pref[i])][1]
  }
  genus_rows <- which(b$.rank == "genus")
  n <- length(u)
  row <- rep(NA_integer_, n)
  type <- rep("none", n)
  conf <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    if (is.na(qc[k])) next
    h <- first_hit(qc[k])
    if (!is.na(h)) {
      row[k] <- h
      type[k] <- "exact"
      conf[k] <- 100
      next
    }
    g <- strsplit(qc[k], " ", fixed = TRUE)[[1]][1]
    if (grepl(" ", qc[k])) {
      # fuzzy: names of the same genus, then names starting alike
      pool <- which(startsWith(b$.canon, paste0(g, " ")))
      if (!length(pool)) {
        pool <- which(startsWith(b$.canon, substr(qc[k], 1, 1)) &
                        grepl(" ", b$.canon))
      }
      if (length(pool)) {
        d <- utils::adist(qc[k], b$.canon[pool])[1, ]
        lim <- min(max_dist, floor(0.1 * nchar(qc[k])) + 1)
        best <- which(d == min(d))
        if (min(d) <= lim) {
          cand <- pool[best]
          row[k] <- cand[order(pref[cand])][1]
          type[k] <- "fuzzy"
          conf[k] <- round(100 * (1 - min(d) / nchar(qc[k])))
          next
        }
      }
    }
    gh <- genus_rows[b$.canon[genus_rows] == g]
    if (length(gh)) {
      row[k] <- gh[order(pref[gh])][1]
      type[k] <- if (grepl(" ", qc[k])) "higherrank" else "exact"
      conf[k] <- if (grepl(" ", qc[k])) NA_real_ else 100
    }
  }
  ok <- !is.na(row)
  matched <- rep(NA_character_, n)
  rank <- status <- acc <- acc_sp <- rep(NA_character_, n)
  matched[ok] <- b$scientificName[row[ok]]
  rank[ok] <- b$.rank[row[ok]]
  status[ok] <- b$.status[row[ok]]
  # follow synonyms to the accepted taxon
  arow <- row
  syn <- ok & status != "accepted" & nzchar(b$acceptedNameUsageID[row])
  arow[which(syn)] <- match(b$acceptedNameUsageID[row[which(syn)]],
                            b$taxonID)
  aok <- !is.na(arow)
  acc[aok] <- b$scientificName[arow[aok]]
  # species of the accepted taxon: itself, or its parent for infraspecific
  for (k in which(aok)) {
    r <- arow[k]
    lv <- rank_level(b$.rank[r])
    if (!is.na(lv) && lv == 8) acc_sp[k] <- b$scientificName[r]
    if (!is.na(lv) && lv > 8) {
      pr <- match(b$parentNameUsageID[r], b$taxonID)
      acc_sp[k] <- if (!is.na(pr)) b$scientificName[pr] else {
        w <- strsplit(parse_names(b$scientificName[r])$canonical, " ")[[1]]
        paste(w[1:2], collapse = " ")
      }
    }
  }
  match_rows(u, name_matched = matched, name_accepted = acc,
             name_accepted_species = acc_sp, name_rank = rank,
             name_status = status, name_match_type = type,
             name_confidence = conf)
}

# Summary of a name matching in plain words (internal)
say_matches <- function(out, res, service) {
  n <- nrow(out)
  cnt <- function(f) sum(out[[f]] %in% TRUE)
  tab <- sort(table(out$name_accepted_species[!out$name_unmatched]),
              decreasing = TRUE)
  message(
    "Names matched (", service, "): ", n, " records, ", nrow(res),
    " distinct names. Exact: ", sum(out$name_match_type %in% "exact"),
    " records; fuzzy (spelled differently, check them): ", cnt("name_fuzzy"),
    "; synonyms resolved to an accepted name: ", cnt("name_synonym"),
    "; only a higher rank found: ", cnt("name_higherrank"),
    "; not found: ", cnt("name_unmatched"),
    "; hybrid names: ", cnt("name_hybrid"), "; aggregates: ",
    cnt("name_aggregate"), ".",
    if (length(tab)) paste0(" Accepted species: ", paste0(
      utils::head(names(tab), 5), " (", utils::head(as.integer(tab), 5),
      ")", collapse = ", "), if (length(tab) > 5) ", ..." else "", "."))
}
