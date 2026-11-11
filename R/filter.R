# Rule-based filtering of occurrence records ------------------------------------

#' Filter occurrence records with transparent, ordered rules
#'
#' Removes records that are unsuitable for niche modelling: cultivated
#' plants, rounded or default coordinates, atlas-grid cell centres, points in
#' towns, car parks or the sea, records far outside the known range and
#' climatic or geographic outliers. It is designed especially for **GBIF
#' downloads** (it uses GBIF columns such as `issue`, `basisOfRecord`,
#' `datasetKey` and `coordinateUncertaintyInMeters`), but works on any table
#' with coordinates: tests whose columns are missing are skipped.
#'
#' Every test can be switched on or off and every threshold changed (see
#' [filter_settings()]). Nothing is removed silently: each record keeps the
#' **first** test it failed in the `reason` column, and is kept only if it
#' passes every active test. Tests that need a layer (range polygon,
#' shoreline, climate, land cover) run only when the layer is supplied through
#' [filter_layers()].
#'
#' @section Tests, in the order they run:
#' **Records and metadata**
#' * `invalid_coordinates` (always): missing or impossible coordinates.
#' * `duplicate`: same species at the same location (rounded to
#'   `dedup_digits` decimals).
#' * `grid_cell_centre`: coordinate uncertainty equal to half the diagonal of
#'   an atlas grid square (e.g. 707.1 m for 1 km squares, as in BSBI monads):
#'   the point is a cell centre, not a locality.
#' * `gbif_issue`, `basis_of_record`, `establishment_means`: GBIF flags,
#'   living or fossil specimens, cultivated or introduced plants.
#' * `coord_uncertainty`: uncertainty above `max_coord_unc_m` (and, with
#'   `drop_na_uncertainty`, records without one: `no_coord_uncertainty`).
#' * `coarse_decimals`: both coordinates written with few decimals.
#' * `garden_text`, `settlement_text`: gardens, car parks, campsites, etc.
#'   mentioned in locality, habitat or remarks.
#' * `locality_is_country`: the locality is only a country name (a default
#'   position). Needs the `countrycode` package.
#' * `too_old`: before `min_year`.
#' * `manual_exclude`, `excluded_dataset`: your own lists.
#'
#' **Gridded datasets** (computed per GBIF `datasetKey`, before duplicates
#' are removed, because repeated coordinates are the signal)
#' * `grid_uncertainty_block`, `lattice_dataset`, `stacked_dataset`: atlas or
#'   grid data recognised from a block of identical large uncertainties, a
#'   regular mesh of coordinates, or many records stacked on few locations.
#'
#' **CoordinateCleaner** (needs the `CoordinateCleaner` package; please cite
#' Zizka et al. 2019)
#' * `country_mismatch`: outside the stated country, with a tolerance buffer
#'   (`CoordinateCleaner::cc_coun()` against Natural Earth `iso_a3_eh`).
#' * `coordinatecleaner`: `CoordinateCleaner::clean_coordinates()` with the
#'   tests in `cc_tests` (country centroids, identical lon/lat, GBIF
#'   headquarters, institutions, zeros, urban areas).
#'
#' **Layers**
#' * `outside_range`: farther than `range_buffer_km` from the range polygon.
#' * `at_sea`: farther from land than allowed by the coordinate uncertainty,
#'   measured to detailed land polygons; with only a water raster, the coarse
#'   fallback reason is `on_water`.
#' * `no_climate_data`: no value in the climate raster.
#' * `built_up`, `cropland`, `ghs_dense`, `ghs_building_nearby`: settlements
#'   and fields.
#'
#' **Isolation and climate** (run last, on records that passed everything
#' else; records inside the range polygon are exempt from the isolation tests
#' when `protect_in_range = TRUE`)
#' * `isolated_climate_outlier`: beyond `isolation_climate` and with fewer
#'   than `min_support` other records within `support_radius_km`.
#' * `small_isolated_cluster`: clusters of fewer than `min_cluster_n`
#'   records (chained closer than `cluster_eps_km`).
#' * `climate_limit`: beyond the hard `climate_limit`.
#' * `climate_soft_limit`: beyond `climate_soft_limit` with no cell within
#'   `refuge_limit` in `refuge_radius_km` (e.g. a warm valley floor is kept
#'   when cold mountains are nearby; a lowland park is not).
#'
#' Finally, records in `keep_ids` are kept whatever they failed (the reason is
#' moved to `reason_overridden`).
#'
#' @param x A data.frame of occurrences (e.g. from [read_gbif()]) or a path to
#'   a GBIF download.
#' @param settings A [filter_settings()] object (default: the strict preset).
#' @param layers A [filter_layers()] object with the optional spatial layers.
#' @param exclude_ids,keep_ids Record IDs (values of the `id` column, GBIF
#'   `gbifID`) to remove or to keep regardless of any test.
#' @param exclude_datasets GBIF `datasetKey` values to remove.
#' @param exclude_label_regex Regular expression on the dataset label
#'   (institution code or `occurrenceID` prefix), e.g. `"SIVIM|INFOFLORA"`.
#' @param cols Named character vector to override the column names, e.g.
#'   `c(lon = "x", lat = "y", uncertainty = "unc_m")`. Roles and their GBIF
#'   defaults: `id` (gbifID), `dataset` (datasetKey), `occurrence_id`,
#'   `institution` (institutionCode), `issue`, `basis` (basisOfRecord),
#'   `establishment` (establishmentMeans), `country` (countryCode, ISO2),
#'   `uncertainty` (coordinateUncertaintyInMeters), `year`, `locality`,
#'   `habitat`, `remarks` (occurrenceRemarks), `species`, `lon`
#'   (decimalLongitude), `lat` (decimalLatitude). Matching is
#'   case-insensitive.
#' @param verbose Print progress and the number of records each test removes.
#' @return An object of class `occ_filter`: a list with `data` (all records,
#'   with `reason` = `NA` for kept records and diagnostic columns such as
#'   `climate`, `built_pct`, `dist_land_m`), `log` (records removed per test,
#'   and skipped tests with the reason), `datasets` (one row per dataset with
#'   the gridded-data diagnostics) and `settings`. Use [kept()], [removed()],
#'   [filter_report()], [filter_sensitivity()] and [write_filter_results()].
#' @references Zizka, A. et al. (2019) CoordinateCleaner: standardized
#'   cleaning of occurrence records from biological collection databases.
#'   *Methods in Ecology and Evolution* 10, 744-751.
#'   \doi{10.1111/2041-210X.13152}
#' @export
#' @examples
#' occ <- data.frame(
#'   gbifID = as.character(1:6), species = "Dryas octopetala",
#'   decimalLongitude = c("10.12345", "10.12345", "11.2", "12.54321",
#'                        "13.11111", "14.22222"),
#'   decimalLatitude = c("46.12345", "46.12345", "47.1", "46.54321",
#'                       "47.11111", "46.33333"),
#'   coordinateUncertaintyInMeters = c("30", "30", "50", "707.1", "5000", NA),
#'   basisOfRecord = c(rep("HUMAN_OBSERVATION", 5), "LIVING_SPECIMEN"))
#' # (the default CoordinateCleaner tests include "urban", which downloads
#' #  Natural Earth urban areas; here only offline tests are used)
#' st <- filter_settings("strict", cc_tests = c("centroids", "equal", "zeros"))
#' res <- filter_occurrences(occ, st)
#' res
#' kept(res)
filter_occurrences <- function(x, settings = filter_settings(),
                               layers = filter_layers(),
                               exclude_ids = NULL, keep_ids = NULL,
                               exclude_datasets = NULL,
                               exclude_label_regex = NULL,
                               cols = NULL, verbose = TRUE) {
  if (!inherits(settings, "filter_settings")) {
    stop("`settings` must come from filter_settings().", call. = FALSE)
  }
  if (!inherits(layers, "filter_layers")) {
    stop("`layers` must come from filter_layers().", call. = FALSE)
  }
  if (is.character(x) && length(x) == 1) x <- read_gbif(x)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  st <- settings
  cn <- resolve_columns(x, cols)
  if (is.na(cn[["lon"]]) || is.na(cn[["lat"]])) {
    stop("Longitude/latitude columns not found. Give them with ",
         "`cols = c(lon = \"...\", lat = \"...\")`.", call. = FALSE)
  }
  n <- nrow(x)
  say <- function(...) if (verbose) message(...)
  get_chr <- function(role) {
    if (is.na(cn[[role]])) NULL else as.character(x[[cn[[role]]]])
  }

  # --- bookkeeping ----------------------------------------------------------
  reason <- rep(NA_character_, n)
  log <- list()
  add_log <- function(test, rsn, k, status, note = "") {
    log[[length(log) + 1]] <<- data.frame(test = test, reason = rsn,
                                          removed = k, status = status,
                                          note = note,
                                          stringsAsFactors = FALSE)
  }
  flag <- function(bad, rsn, test = rsn) {
    bad <- !is.na(bad) & bad
    new <- is.na(reason) & bad
    reason[new] <<- rsn
    add_log(test, rsn, sum(new), "run")
    say(sprintf("  %-26s removed: %7d", rsn, sum(new)))
  }
  skip <- function(test, note) {
    if (test %in% st$tests) {
      add_log(test, test, NA_integer_, "skipped", note)
      say(sprintf("  %-26s skipped: %s", test, note))
    }
  }
  on <- function(test) test %in% st$tests
  alive <- function() is.na(reason)

  # --- coordinates, uncertainty, decimals ------------------------------------
  lon_txt <- trimws(get_chr("lon"))
  lat_txt <- trimws(get_chr("lat"))
  lon <- suppressWarnings(as.numeric(lon_txt))
  lat <- suppressWarnings(as.numeric(lat_txt))
  valid <- !is.na(lon) & !is.na(lat) & abs(lon) <= 180 & abs(lat) <= 90
  unc_chr <- get_chr("uncertainty")
  unc <- if (is.null(unc_chr)) rep(NA_real_, n) else
    suppressWarnings(as.numeric(unc_chr))
  n_dec <- function(t) {
    t[is.na(t)] <- ""
    ifelse(grepl("\\.", t), nchar(sub("^[^.]*\\.", "", t)), 0L)
  }
  out <- x
  out$lat_decimals <- n_dec(lat_txt)
  out$lon_decimals <- n_dec(lon_txt)

  say("Filtering ", n, " records (preset: ", st$preset, ")")
  flag(!valid, "invalid_coordinates")

  # --- dataset diagnostics (before duplicates are removed) -------------------
  digits <- if (is.na(st$dedup_digits)) 4 else st$dedup_digits
  lat_r <- round(lat, digits)
  lon_r <- round(lon, digits)
  ds <- get_chr("dataset")
  ds_stats <- NULL
  if (!is.null(ds)) {
    inst <- get_chr("institution")
    occid <- get_chr("occurrence_id")
    lab <- if (!is.null(inst)) inst else rep(NA_character_, n)
    if (!is.null(occid)) {
      lab <- ifelse(is.na(lab), sub("[-:_/ ]*[0-9].*$", "", occid), lab)
    }
    out$dataset_label <- substr(lab, 1, 40)
    unc_key <- ifelse(is.na(unc), "NA", format(unc, scientific = FALSE))
    ds_key <- ifelse(is.na(ds), "NA", ds)
    out$block_n <- as.integer(stats::ave(seq_len(n), ds_key, unc_key,
                                         FUN = length))
    basis <- get_chr("basis")
    ds_stats <- dataset_diagnostics(ds_key[valid], out$dataset_label[valid],
                                    lat_r[valid], lon_r[valid], unc[valid],
                                    if (is.null(basis)) NULL else basis[valid],
                                    st)
    out$lattice_dataset <- ds_stats$lattice[match(ds_key, ds_stats$dataset)]
    out$stacked_dataset <- ds_stats$stacked[match(ds_key, ds_stats$dataset)]
  }

  # --- 1. duplicates ----------------------------------------------------------
  if (on("duplicate")) {
    if (is.na(st$dedup_digits)) {
      skip("duplicate", "dedup_digits is NA")
    } else {
      sp <- get_chr("species")
      key <- paste(if (is.null(sp)) "" else sp, lat_r, lon_r)
      dup <- rep(FALSE, n)
      dup[valid] <- duplicated(key[valid])
      flag(dup, "duplicate")
    }
  }

  # --- 2. metadata --------------------------------------------------------------
  has_unc <- !is.null(unc_chr)
  if (on("grid_cell_centre")) {
    if (!has_unc) skip("grid_cell_centre", "no uncertainty column") else {
      h <- st$grid_cell_sizes_m * sqrt(2) / 2
      gc <- Reduce(`|`, lapply(h, function(v) {
        !is.na(unc) & abs(unc - v) <= pmax(0.6, 0.003 * v)
      }))
      out$grid_centre <- gc
      flag(gc, "grid_cell_centre")
    }
  }
  test_values <- function(test, role, values, rsn = test, upper = FALSE) {
    if (!on(test)) return(invisible())
    v <- get_chr(role)
    if (is.null(v)) return(skip(test, paste0("no ", role, " column")))
    if (upper) v <- toupper(v)
    flag(v %in% values, rsn, test)
  }
  if (on("gbif_issue")) {
    iss <- get_chr("issue")
    if (is.null(iss)) skip("gbif_issue", "no issue column") else
      flag(grepl(paste(st$bad_issues, collapse = "|"), iss), "gbif_issue")
  }
  test_values("basis_of_record", "basis", st$bad_basis, upper = TRUE)
  test_values("establishment_means", "establishment", st$bad_establishment,
              upper = TRUE)
  if (on("coord_uncertainty")) {
    if (!has_unc) skip("coord_uncertainty", "no uncertainty column") else {
      flag(!is.na(unc) & unc > st$max_coord_unc_m, "coord_uncertainty")
      if (isTRUE(st$drop_na_uncertainty)) {
        flag(is.na(unc), "no_coord_uncertainty", "coord_uncertainty")
      }
    }
  }
  if (on("coarse_decimals")) {
    if (st$min_coord_decimals <= 0) skip("coarse_decimals", "set to 0") else
      flag(out$lat_decimals < st$min_coord_decimals &
             out$lon_decimals < st$min_coord_decimals, "coarse_decimals")
  }
  text <- NULL
  txt_cols <- Filter(Negate(is.null), lapply(c("locality", "habitat",
                                               "remarks"), get_chr))
  if (length(txt_cols)) {
    txt_cols <- lapply(txt_cols, function(v) ifelse(is.na(v), "", v))
    text <- do.call(paste, txt_cols)
  }
  for (tt in c("garden_text", "settlement_text")) {
    if (!on(tt)) next
    if (is.null(text)) {
      skip(tt, "no locality, habitat or remarks column")
    } else {
      rx <- if (tt == "garden_text") st$garden_regex else st$settlement_regex
      flag(grepl(rx, text, ignore.case = TRUE), tt)
    }
  }
  if (on("locality_is_country")) {
    loc <- get_chr("locality")
    if (is.null(loc)) {
      skip("locality_is_country", "no locality column")
    } else if (!requireNamespace("countrycode", quietly = TRUE)) {
      skip("locality_is_country", "install the countrycode package")
    } else {
      flag(normalise_text(loc) %in% country_names(), "locality_is_country")
    }
  }
  if (on("too_old")) {
    yr <- get_chr("year")
    if (is.na(st$min_year)) skip("too_old", "min_year is NA") else
      if (is.null(yr)) skip("too_old", "no year column") else
        flag(suppressWarnings(as.numeric(yr)) < st$min_year, "too_old")
  }
  ids <- get_chr("id")
  if (on("manual_exclude") && length(exclude_ids)) {
    if (is.null(ids)) skip("manual_exclude", "no id column") else
      flag(ids %in% as.character(exclude_ids), "manual_exclude")
  }
  if (on("excluded_dataset") &&
      (length(exclude_datasets) || !is.null(exclude_label_regex))) {
    if (is.null(ds)) skip("excluded_dataset", "no dataset column") else {
      flag(ds %in% exclude_datasets, "excluded_dataset")
      if (!is.null(exclude_label_regex)) {
        flag(grepl(exclude_label_regex, out$dataset_label,
                   ignore.case = TRUE), "excluded_label", "excluded_dataset")
      }
    }
  }

  # --- 3. gridded datasets --------------------------------------------------------
  for (tt in c("grid_uncertainty_block", "lattice_dataset",
               "stacked_dataset")) {
    if (!on(tt)) next
    if (is.null(ds)) {
      skip(tt, "no dataset column")
      next
    }
    bad <- switch(tt,
      grid_uncertainty_block = !is.na(unc) & unc >= st$grid_min_unc_m &
        out$block_n >= st$grid_block_n,
      lattice_dataset = out$lattice_dataset,
      stacked_dataset = out$stacked_dataset)
    flag(bad, tt)
  }

  # --- 4. CoordinateCleaner -------------------------------------------------------
  if (on("country_mismatch")) {
    run_country_check(get_chr("country"), lon, lat, alive() & valid, st,
                      flag, skip)
  }
  if (on("coordinatecleaner")) {
    run_cc(lon, lat, alive() & valid, st, layers$urban, flag, skip, say)
  }

  # --- 5. range polygon -------------------------------------------------------------
  in_core <- rep(FALSE, n)
  if (!is.null(layers$range)) {
    rr <- range_membership(lon, lat, valid, layers$range, st)
    in_core <- rr$core
    out$in_range <- rr$core
    out$in_range_zone <- rr$zone
    say("  range polygon: ", sum(rr$core, na.rm = TRUE), " records inside, ",
        sum(rr$zone, na.rm = TRUE), " inside + ", st$range_buffer_km, " km")
    if (on("outside_range")) flag(valid & !rr$zone, "outside_range")
  } else {
    skip("outside_range", "no range polygon in layers")
  }

  # --- 6. sea, climate, settlement --------------------------------------------------
  if (on("at_sea")) {
    if (is.null(layers$land) && is.null(layers$water)) {
      skip("at_sea", "no land polygons or water raster in layers")
    } else {
      sea <- sea_test(lon, lat, unc, valid & alive(), layers, st, say)
      out$water_near_max <- sea$water_max
      out$land_near_pct <- sea$land_max
      out$dist_land_m <- sea$dist
      if (!is.null(layers$land)) {
        flag(!is.na(sea$dist) &
               sea$dist > pmax(st$sea_unc_factor * ifelse(is.na(unc), 0, unc),
                               st$sea_tol_m), "at_sea")
      }
      if (!is.null(layers$water)) {
        flag(is.na(sea$dist) & !is.na(sea$land_max) &
               sea$land_max < st$sea_land_pct, "on_water", "at_sea")
      }
    }
  }
  clim <- NULL
  if (!is.null(layers$climate)) {
    clim <- rep(NA_real_, n)
    clim[valid] <- extract_at(layers$climate, lon[valid], lat[valid])
    out$climate <- clim
  }
  if (on("no_climate_data")) {
    if (is.null(clim)) skip("no_climate_data", "no climate layer") else
      if (isTRUE(st$drop_na_climate)) flag(valid & is.na(clim),
                                           "no_climate_data")
  }
  ras_pct <- function(r, what) {
    v <- rep(NA_real_, n)
    v[valid] <- extract_at(r, lon[valid], lat[valid])
    as_percent(v, what, say)
  }
  if (!is.null(layers$built)) out$built_pct <- ras_pct(layers$built, "built")
  if (!is.null(layers$cropland)) {
    out$crop_pct <- ras_pct(layers$cropland, "cropland")
  }
  if (!is.null(layers$ghs_mean)) {
    v <- rep(NA_real_, n)
    v[valid] <- extract_at(layers$ghs_mean, lon[valid], lat[valid])
    out$ghs_mean_pct <- v / 100       # m2 per 10,000 m2 (100 m cell) -> %
  }
  if (!is.null(layers$ghs_max)) {
    v <- rep(NA_real_, n)
    v[valid] <- extract_at(layers$ghs_max, lon[valid], lat[valid])
    out$ghs_max_m2 <- v
  }
  layer_test <- function(test, col, limit, lay) {
    if (!on(test)) return(invisible())
    if (is.null(out[[col]])) return(skip(test, paste("no", lay, "layer")))
    flag(out[[col]] >= limit, test)
  }
  layer_test("built_up", "built_pct", st$max_built_pct, "built")
  layer_test("cropland", "crop_pct", st$max_crop_pct, "cropland")
  layer_test("ghs_dense", "ghs_mean_pct", st$ghs_mean_pct, "ghs_mean")
  layer_test("ghs_building_nearby", "ghs_max_m2", st$ghs_max_m2, "ghs_max")

  # --- 7. isolation ------------------------------------------------------------------
  sgn <- if (st$climate_direction == "above") 1 else -1
  exempt <- if (isTRUE(st$protect_in_range)) in_core & !is.na(in_core) else
    rep(FALSE, n)
  if (on("isolated_climate_outlier")) {
    if (is.null(clim)) {
      skip("isolated_climate_outlier", "no climate layer")
    } else if (is.na(st$isolation_climate)) {
      skip("isolated_climate_outlier", "isolation_climate is NA")
    } else {
      a <- which(alive() & valid)
      nb <- neighbour_counts(lon[a], lat[a], st$support_radius_km)
      bad <- rep(FALSE, n)
      bad[a] <- !is.na(clim[a]) & sgn * clim[a] > sgn * st$isolation_climate &
        nb < st$min_support & !exempt[a]
      flag(bad, "isolated_climate_outlier")
    }
  }
  if (on("small_isolated_cluster")) {
    a <- which(alive() & valid)
    if (length(a)) {
      sz <- cluster_sizes(lon[a], lat[a], st$cluster_eps_km)
      bad <- rep(FALSE, n)
      bad[a] <- sz < st$min_cluster_n & !exempt[a]
      flag(bad, "small_isolated_cluster")
    }
  }

  # --- 8. climate limits (last) ------------------------------------------------------
  if (on("climate_limit")) {
    if (is.null(clim)) skip("climate_limit", "no climate layer") else
      if (is.na(st$climate_limit)) skip("climate_limit", "climate_limit is NA")
    else flag(sgn * clim > sgn * st$climate_limit, "climate_limit")
  }
  if (on("climate_soft_limit")) {
    if (is.null(clim)) {
      skip("climate_soft_limit", "no climate layer")
    } else if (is.na(st$climate_soft_limit) || is.na(st$refuge_limit)) {
      skip("climate_soft_limit", "climate_soft_limit or refuge_limit is NA")
    } else {
      cand <- which(alive() & !is.na(clim) &
                      sgn * clim > sgn * st$climate_soft_limit)
      near <- rep(NA_real_, n)
      say("  records beyond the soft limit, searched for a refuge within ",
          st$refuge_radius_km, " km: ", length(cand))
      if (length(cand)) {
        near[cand] <- vapply(cand, function(i) {
          window_extreme(layers$climate, lon[i], lat[i], st$refuge_radius_km,
                         if (sgn > 0) min else max)
        }, numeric(1))
      }
      out$climate_near <- near
      bad <- rep(FALSE, n)
      bad[cand] <- is.na(near[cand]) |
        sgn * near[cand] > sgn * st$refuge_limit
      flag(bad, "climate_soft_limit")
    }
  }

  # --- 9. manual keep list -------------------------------------------------------------
  out$reason_overridden <- NA_character_
  if (length(keep_ids)) {
    if (is.null(ids)) {
      warning("`keep_ids` given but no id column found.", call. = FALSE)
    } else {
      k <- ids %in% as.character(keep_ids) & valid
      out$reason_overridden[k] <- reason[k]
      say("  kept by hand: ", sum(k), " (", sum(k & !is.na(reason)),
          " would otherwise have been removed)")
      reason[k] <- NA_character_
    }
  }
  out$reason <- reason

  res <- list(data = out, log = do.call(rbind, log), datasets = ds_stats,
              settings = st, columns = cn)
  class(res) <- "occ_filter"
  say("Kept ", sum(is.na(reason)), " of ", n, " records")
  res
}

# GBIF default column names and user overrides -> actual column names
resolve_columns <- function(x, cols = NULL) {
  roles <- c(id = "gbifID", dataset = "datasetKey",
             occurrence_id = "occurrenceID", institution = "institutionCode",
             issue = "issue", basis = "basisOfRecord",
             establishment = "establishmentMeans", country = "countryCode",
             uncertainty = "coordinateUncertaintyInMeters", year = "year",
             locality = "locality", habitat = "habitat",
             remarks = "occurrenceRemarks", species = "species",
             lon = "decimalLongitude", lat = "decimalLatitude")
  if (!is.null(cols)) {
    unknown <- setdiff(names(cols), names(roles))
    if (length(unknown)) {
      stop("Unknown column role(s) in `cols`: ",
           paste(unknown, collapse = ", "), ". Roles: ",
           paste(names(roles), collapse = ", "), call. = FALSE)
    }
    roles[names(cols)] <- cols
  }
  out <- vapply(roles, function(nm) find_column(names(x), nm), character(1))
  if (is.na(out[["lon"]]) || is.na(out[["lat"]])) {
    xy <- tryCatch(guess_coords(names(x)), error = function(e) NULL)
    if (!is.null(xy)) out[c("lon", "lat")] <- xy
  }
  out
}

# Per-dataset statistics for the gridded-data tests (internal)
dataset_diagnostics <- function(ds, label, lat_r, lon_r, unc, basis, st) {
  idx <- split(seq_along(ds), ds)
  rows <- lapply(names(idx), function(k) {
    i <- idx[[k]]
    u <- unc[i]
    tu <- table(u)
    data.frame(
      dataset = k,
      label = label[i][!is.na(label[i])][1],
      n_records = length(i),
      n_coords = nrow(unique(cbind(lat_r[i], lon_r[i]))),
      n_lat = length(unique(lat_r[i])),
      n_lon = length(unique(lon_r[i])),
      coarse_share = mean(is.na(u) | u >= st$grid_min_unc_m),
      specimen_share = if (is.null(basis)) NA_real_ else
        mean(basis[i] %in% "PRESERVED_SPECIMEN"),
      modal_unc = if (length(tu)) as.numeric(names(tu)[which.max(tu)]) else
        NA_real_,
      modal_unc_share = if (length(tu)) max(tu) / length(i) else NA_real_,
      stringsAsFactors = FALSE)
  })
  d <- do.call(rbind, rows)
  d$lattice <- d$n_coords >= st$lattice_min_coords &
    d$n_lat / d$n_coords <= st$lattice_ratio &
    d$n_lon / d$n_coords <= st$lattice_ratio
  spec_ok <- is.na(d$specimen_share) |
    d$specimen_share < st$stack_max_specimen_share
  d$stacked <- d$n_records >= st$stack_min_records &
    d$n_records / d$n_coords >= st$stack_ratio &
    d$coarse_share >= st$stack_coarse_share & spec_ok
  d[order(-d$n_records), , drop = FALSE]
}

# Lower case, no accents or punctuation (internal)
normalise_text <- function(x) {
  x <- if (requireNamespace("stringi", quietly = TRUE)) {
    stringi::stri_trans_general(x, "Latin-ASCII")
  } else {
    iconv(x, "UTF-8", "ASCII//TRANSLIT", sub = "")
  }
  x <- tolower(x)
  trimws(gsub(" +", " ", gsub("[^a-z ]", " ", x)))
}

# Country names in many languages, normalised (internal; needs countrycode)
country_names <- function() {
  cl <- countrycode::codelist
  cols <- grep("^(country\\.name|cldr\\.(name|short))\\.", names(cl),
               value = TRUE)
  cols <- cols[!grepl("regex", cols)]
  nm <- normalise_text(unlist(cl[cols], use.names = FALSE))
  unique(nm[!is.na(nm) & nzchar(nm)])
}

# Convert 0-1 fractions to percent (internal)
as_percent <- function(v, label, say = message) {
  mx <- suppressWarnings(max(v, na.rm = TRUE))
  if (!is.finite(mx)) {
    warning(label, ": every extracted value is NA; that test cannot work.",
            call. = FALSE)
  } else if (mx <= 1) {
    say("  ", label, ": values look like 0-1 fractions, converted to percent")
    v <- v * 100
  }
  v
}
