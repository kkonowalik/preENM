# Input for integrated species distribution models ------------------------------

#' Shape occurrence data for integrated species distribution models
#'
#' Integrated species distribution models (ISDMs) combine presence-only,
#' presence-absence and count data, each with its own observation model,
#' in one joint likelihood (Morera-Pujol et al. 2023; Zulian et al. 2021).
#' Model-based integration improved models of data-deficient, narrow-ranged
#' species compared with presence-absence data alone, though only slightly,
#' and improved them more clearly compared with presence-only data
#' (Makinen et al. 2026). preENM fits no model; `export_integrated()`
#' shapes the cleaned records into the input that the modelling packages
#' expect, one dataset per data type (or per source), with its
#' coordinates, response and the metadata they need. It does not depend on
#' those packages.
#'
#' **Formats**
#' * `"generic"`: a list of data.frames with `species`, `x`, `y`,
#'   `response`, `data_type`, `source` (and `trials`).
#' * `"intSDM"` (Mostert et al. 2025): one list per dataset with the
#'   arguments of the workflow's `$addStructured()`: `dataStructured` (an
#'   `sf` object, so intSDM transforms it to the workflow's projection),
#'   `datasetType` (`"PO"`, `"PA"` or `"Counts"`), `responseName`,
#'   `trialsName`, `speciesName` and `datasetName`. Call
#'   `do.call(workflow$addStructured, d)` for each.
#' * `"PointedSDMs"`: `data`, a named list of `sf` objects for
#'   `startISDM()` / `startSpecies()`, in which presence-absence datasets
#'   have the response column `present`, count datasets `counts` and
#'   presence-only datasets no response column (that is how PointedSDMs
#'   recognises them); and `args`, the matching arguments (`responsePA`,
#'   `responseCounts`, `trialsPA`, `Projection`, `speciesName`).
#' * `"biomod2"` (Gueguen et al. 2026): one list per dataset and species
#'   with the arguments of `BIOMOD_FormatingData()`: `resp.name` (with
#'   spaces and underscores replaced by dots, as biomod2 does), `resp.var`,
#'   `resp.xy` and `data.type` (`"binary"` for presence-only and
#'   presence-absence data, `"count"` for whole-number counts,
#'   `"abundance"` otherwise). biomod2 models each data type separately
#'   and does not integrate them; presence-only datasets need its
#'   pseudo-absence arguments (`PA.nb.rep`, ...).
#'
#' **Checks.** Records without a response are dropped; presence-only
#' records with response 0 (absences inside presence-only data) are
#' dropped, and reported. The printout warns about presence-absence
#' datasets without absences or without presences (they cannot inform a
#' presence-absence likelihood), about very low or high prevalence and
#' about small datasets.
#'
#' @param x Records with `data_type` and `response` (from
#'   [merge_occurrences()], after [filter_occurrences()]: its result or
#'   [kept()] of it).
#' @param format `"generic"`, `"intSDM"`, `"PointedSDMs"` or `"biomod2"`.
#' @param split One dataset per `"data_type"` (default) or per `"source"`.
#'   Integrated models usually give each dataset its own intercept, so
#'   sources collected with different protocols are better kept apart.
#' @param types For data types with your own labels, the observation model
#'   of each, e.g. `c(atlas = "presence_absence", iNaturalist =
#'   "presence_only")`.
#' @param trials Optional column with the number of trials (visits) of
#'   presence-absence records.
#' @param crs CRS of the exported coordinates (default WGS84
#'   longitude/latitude); give the projection of your model.
#' @param coords,points_crs Coordinate columns and their CRS (default the
#'   `lon`, `lat` of [merge_occurrences()], WGS84).
#' @return An object of class `integrated_data`: `datasets` (the generic
#'   data.frames), `export` (the structure for `format`), `summary` (one
#'   row per dataset: type, records, presences, absences, prevalence,
#'   species), `dropped` and `settings`. See [write_integrated()].
#' @references Morera-Pujol, V., Mostert, P.S., Murphy, K.J., Burkitt, T.,
#'   Coad, B., McMahon, B.J., Nieuwenhuis, M., Morelle, K., Ward, A.I.,
#'   Ciuti, S. (2023) Bayesian species distribution models integrate
#'   presence-only and presence-absence data to predict deer distribution
#'   and relative abundance. *Ecography* 2023, e06451.
#'   \doi{10.1111/ecog.06451}
#'
#'   Zulian, V., Miller, D.A.W., Ferraz, G. (2021) Integrating
#'   citizen-science and planned-survey data improves species distribution
#'   estimates. *Diversity and Distributions* 27, 2498-2509.
#'   \doi{10.1111/ddi.13416}
#'
#'   Mostert, P.S., O'Hara, R.B., Finstad, A.G., Isaac, N.J.B. et al.
#'   (2025) intSDM: an R package for building a reproducible workflow for
#'   the field of integrated species distribution models. *Ecology and
#'   Evolution* 15, e71029. \doi{10.1002/ece3.71029}
#'
#'   Gueguen, M. et al. (2026) biomod2: extending presence-absence species
#'   distribution models to multiple data types. *Ecography* 2026,
#'   e08588. \doi{10.1002/ecog.08588}
#'
#'   Makinen, J., Cohen, J., Jetz, W. (2026) Model-based data integration
#'   improves species distribution models for data-deficient and
#'   narrow-ranged hummingbird species. *Diversity and Distributions*,
#'   e70157. \doi{10.1111/ddi.70157}
#' @seealso [merge_occurrences()], [filter_occurrences()]
#' @export
#' @examples
#' po <- data.frame(lon = c(10.1, 10.4, 11.2), lat = c(46.1, 46.3, 46.8),
#'                  species = "Aus bus")
#' pa <- data.frame(lon = c(10.2, 10.6, 11.0, 11.3),
#'                  lat = c(46.2, 46.5, 46.6, 46.9),
#'                  species = "Aus bus", present = c(1, 0, 0, 1))
#' m <- merge_occurrences(GBIF = po, survey = pa,
#'                        data_type = c(survey = "presence_absence"))
#' ex <- export_integrated(m, "PointedSDMs")
#' ex
#' names(ex$export$data)
#' ex$export$args
export_integrated <- function(x, format = c("generic", "intSDM",
                                            "PointedSDMs", "biomod2"),
                              split = c("data_type", "source"),
                              types = NULL, trials = NULL,
                              crs = "EPSG:4326", coords = NULL,
                              points_crs = "EPSG:4326") {
  format <- match.arg(format)
  split <- match.arg(split)
  if (inherits(x, "occ_filter")) x <- kept(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!nrow(x)) stop("No records.", call. = FALSE)
  if (!"data_type" %in% names(x)) x$data_type <- "presence_only"
  if (!"response" %in% names(x)) {
    if (any(x$data_type != "presence_only")) {
      stop("A `response` column is needed for presence-absence and ",
           "abundance data (see merge_occurrences()).", call. = FALSE)
    }
    x$response <- 1
  }
  if (!"source" %in% names(x)) x$source <- x$data_type
  if (!"species" %in% names(x)) x$species <- "species"
  xy_cols <- if (is.null(coords)) {
    if (all(c("lon", "lat") %in% names(x))) c("lon", "lat") else
      guess_coords(names(x))
  } else coords
  base_types <- c("presence_only", "presence_absence", "abundance")
  dt <- ifelse(is.na(x$data_type), "presence_only", x$data_type)
  obs <- dt
  if (!is.null(types)) {
    bad <- setdiff(types, base_types)
    if (length(bad)) stop("`types` must map to presence_only, ",
                          "presence_absence or abundance.", call. = FALSE)
    hit <- dt %in% names(types)
    obs[hit] <- unname(types[dt[hit]])
  }
  unknown <- setdiff(unique(obs), base_types)
  if (length(unknown)) {
    stop("Data type(s) ", paste(unknown, collapse = ", "), " need an ",
         "observation model: give `types`, e.g. types = c(",
         unknown[1], " = \"presence_only\").", call. = FALSE)
  }
  resp <- as_number(x$response)
  if (!is.null(trials) && !trials %in% names(x)) {
    stop("Column not found: ", trials, call. = FALSE)
  }
  # records that cannot be used
  drop_na <- is.na(resp)
  drop_po0 <- obs == "presence_only" & !is.na(resp) & resp == 0
  dropped <- data.frame(
    reason = c("no response", "absence in presence-only data"),
    n = c(sum(drop_na), sum(drop_po0)), stringsAsFactors = FALSE)
  if (any(drop_na | drop_po0)) {
    message("Dropped: ", sum(drop_na), " record(s) without a response, ",
            sum(drop_po0), " absence(s) in presence-only data")
  }
  ok <- !drop_na & !drop_po0
  x <- x[ok, , drop = FALSE]
  obs <- obs[ok]
  dt <- dt[ok]
  resp <- resp[ok]
  if (!nrow(x)) stop("No usable records.", call. = FALSE)

  # coordinates in the target CRS
  p <- terra::vect(cbind(as_number(x[[xy_cols[1]]]),
                         as_number(x[[xy_cols[2]]])), crs = points_crs)
  if (!terra::same.crs(p, terra::vect(cbind(0, 0), crs = crs))) {
    p <- terra::project(p, crs)
  }
  xy <- terra::crds(p)
  key <- if (split == "data_type") dt else x$source
  groups <- unique(key)
  ds <- lapply(groups, function(g) {
    i <- which(key == g)
    d <- data.frame(species = x$species[i], x = xy[i, 1], y = xy[i, 2],
                    response = resp[i], data_type = obs[i],
                    source = x$source[i], stringsAsFactors = FALSE)
    if (!is.null(trials)) d$trials <- as_number(x[[trials]][i])
    if (length(unique(d$data_type)) > 1) {
      stop("Dataset '", g, "' mixes observation models (",
           paste(unique(d$data_type), collapse = ", "), "); use split = ",
           "\"data_type\" or give each source one type.", call. = FALSE)
    }
    rownames(d) <- NULL
    d
  })
  names(ds) <- make.names(groups, unique = TRUE)

  summ <- do.call(rbind, lapply(names(ds), function(k) {
    d <- ds[[k]]
    t <- d$data_type[1]
    pres <- if (t == "abundance") sum(d$response > 0) else
      sum(d$response == 1)
    data.frame(dataset = k, data_type = t, records = nrow(d),
               presences = pres, absences = nrow(d) - pres,
               prevalence = if (t == "presence_only") NA_real_ else
                 pres / nrow(d),
               species = length(unique(d$species)),
               stringsAsFactors = FALSE)
  }))

  export <- switch(format,
    generic = ds,
    intSDM = lapply(names(ds), function(k) {
      d <- ds[[k]]
      t <- d$data_type[1]
      sfd <- to_sf(d, crs, keep = c("species", if (t != "presence_only")
        "response", if ("trials" %in% names(d) && t == "presence_absence")
          "trials"))
      list(dataStructured = sfd,
           datasetType = c(presence_only = "PO", presence_absence = "PA",
                           abundance = "Counts")[[t]],
           responseName = if (t == "presence_only") NULL else "response",
           trialsName = if ("trials" %in% names(d) &&
                            t == "presence_absence") "trials" else NULL,
           speciesName = "species", datasetName = k)
    }),
    PointedSDMs = {
      dat <- lapply(names(ds), function(k) {
        d <- ds[[k]]
        t <- d$data_type[1]
        if (t == "presence_absence") names(d)[names(d) == "response"] <-
            "present"
        if (t == "abundance") names(d)[names(d) == "response"] <- "counts"
        keep <- c("species", switch(t, presence_absence = "present",
                                    abundance = "counts", NULL),
                  if (t == "presence_absence" && "trials" %in% names(d))
                    "trials")
        to_sf(d, crs, keep)
      })
      names(dat) <- names(ds)
      list(data = dat,
           args = list(responsePA = "present", responseCounts = "counts",
                       trialsPA = if (!is.null(trials)) "trials" else NULL,
                       Projection = crs, speciesName = "species"))
    },
    biomod2 = {
      out <- list()
      for (k in names(ds)) {
        d <- ds[[k]]
        t <- d$data_type[1]
        for (sp in unique(d$species)) {
          e <- d[d$species == sp, , drop = FALSE]
          nm <- gsub("[ _]+", ".", sp)
          out[[paste(k, nm, sep = ":")]] <- list(
            resp.name = nm, resp.var = e$response,
            resp.xy = data.frame(x = e$x, y = e$y),
            data.type = if (t != "abundance") "binary" else
              if (all(e$response == round(e$response))) "count" else
                "abundance")
        }
      }
      out
    })
  if (format == "intSDM") names(export) <- names(ds)
  structure(list(datasets = ds, export = export, summary = summ,
                 dropped = dropped,
                 settings = list(format = format, split = split, crs = crs,
                                 types = types)),
            class = "integrated_data")
}

# A data.frame with x, y as an sf object keeping some columns (internal)
to_sf <- function(d, crs, keep) {
  sf::st_as_sf(d[c(intersect(keep, names(d)), "x", "y")],
               coords = c("x", "y"), crs = crs)
}

#' @export
print.integrated_data <- function(x, ...) {
  s <- x$summary
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  lab <- c(presence_only = "presence-only", presence_absence =
             "presence-absence", abundance = "count/abundance")
  cat("Input for integrated SDMs (", x$settings$format, "): ", nrow(s),
      " dataset(s), ", sum(s$records), " records\n\nIN SHORT\n", sep = "")
  tt <- table(factor(s$data_type, names(lab)))
  say(paste0(paste0(as.integer(tt[tt > 0]), " ", lab[names(tt)[tt > 0]],
                    collapse = ", "), " dataset(s).",
             if (length(unique(s$data_type)) < 2)
               " Only one kind of data: nothing to integrate; a model of that type alone does the same." else ""))
  for (i in seq_len(nrow(s))) {
    r <- s[i, ]
    say(paste0(r$dataset, ": ", r$records, " ", lab[[r$data_type]],
               " records", if (r$data_type != "presence_only") paste0(
                 " (", r$presences, " with presence, ", r$absences,
                 if (r$data_type == "abundance") " zero counts" else
                   " absences", "; prevalence ",
                 share_label(r$prevalence), ")") else "",
               if (r$species > 1) paste0(", ", r$species, " species") else "",
               "."))
    if (r$data_type == "presence_absence" &&
        (r$absences == 0 || r$presences == 0)) {
      say(paste0("Warning: ", r$dataset, " has no ",
                 if (r$absences == 0) "absences" else "presences",
                 ", so it cannot inform a presence-absence likelihood; ",
                 "treat it as presence-only or add the missing sites."))
    } else if (r$data_type != "presence_only" &&
               (r$prevalence < 0.05 || r$prevalence > 0.95)) {
      say(paste0("Prevalence in ", r$dataset, " is extreme (",
                 share_label(r$prevalence), "): the dataset carries little ",
                 "information on where the species is absent or present."))
    }
    if (r$records < 30) {
      say(paste0(r$dataset, " has fewer than 30 records: the parameters ",
                 "that belong to it alone (its intercept and observation ",
                 "process, or a separate model) will be poorly estimated."))
    }
  }
  if (any(x$dropped$n > 0)) {
    say(paste0("Dropped before export: ", paste0(
      x$dropped$reason[x$dropped$n > 0], " (", x$dropped$n[x$dropped$n > 0],
      ")", collapse = "; "), "."))
  }
  cat("\n")
  print(s, row.names = FALSE)
  cat("", strwrap(paste0(
    "prevalence = share of records with presence (or a count above 0). ",
    "The limits (fewer than 30 records, prevalence below 5% or above 95%) ",
    "are rules of thumb. Coordinates are in ", x$settings$crs, "."),
    width = 76), sep = "\n")
  invisible(x)
}

#' Write the datasets for integrated SDMs to files
#'
#' One CSV per dataset (`species`, `x`, `y`, `response`, `data_type`,
#' `source`) and `<prefix>_datasets.csv` with the summary, for software
#' outside R or for archiving the exact input of a model.
#'
#' @param x An `integrated_data` from [export_integrated()].
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @return Invisibly, the paths written.
#' @references Mostert, P.S. et al. (2025) intSDM: an R package for
#'   building a reproducible workflow for the field of integrated species
#'   distribution models. *Ecology and Evolution* 15, e71029.
#'   \doi{10.1002/ece3.71029}
#' @export
#' @examples
#' m <- merge_occurrences(GBIF = data.frame(lon = 10, lat = 46),
#'                        survey = data.frame(lon = c(10.2, 10.3),
#'                                            lat = 46, present = c(1, 0)),
#'                        data_type = c(survey = "presence_absence"))
#' write_integrated(export_integrated(m), tempdir())
write_integrated <- function(x, dir, prefix = "isdm") {
  if (!inherits(x, "integrated_data")) {
    stop("`x` must come from export_integrated().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- file.path(dir, paste0(prefix, "_", names(x$datasets), ".csv"))
  for (i in seq_along(f)) {
    utils::write.csv(x$datasets[[i]], f[i], row.names = FALSE)
  }
  fs <- file.path(dir, paste0(prefix, "_datasets.csv"))
  utils::write.csv(x$summary, fs, row.names = FALSE)
  invisible(c(f, fs))
}
