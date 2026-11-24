# Range metrics for conservation assessment -------------------------------------------

#' Range metrics of a predicted distribution
#'
#' The usual end of a modelling workflow is a suitability map; what an
#' assessment needs are the measures of the range derived from it: its
#' extent and occupancy as used by the IUCN Red List (criterion B), its area
#' and its landscape structure (Curd et al. 2023). Today these are obtained
#' with ConR (Dauby et al. 2017) or by hand in a GIS. `range_metrics()`
#' computes them from binary maps, under every threshold rule given, and
#' reports how much they depend on the rule.
#'
#' **Input.** Binary rasters (1 = suitable, 0 = unsuitable; one layer per
#' projection, e.g. present and future), or continuous predictions with
#' thresholds: a [thresholds()] result (all its rules are used) or a named
#' vector, e.g. `c(maxTSS = 0.42, p10 = 0.31)`. Each layer is a projection;
#' each projection is binarised under each rule.
#'
#' **Metrics** (per projection and rule):
#' * `area_km2`: total suitable area.
#' * `EOO_km2`, extent of occurrence: the area of the minimum convex
#'   polygon around the suitable cells, in an equal-area projection; as in
#'   the IUCN guidelines, it is set to the AOO when smaller.
#' * `AOO_km2`, area of occupancy: the number of cells of an `aoo_km` grid
#'   (2 x 2 km, the IUCN standard) that contain suitable area, times their
#'   area, in an equal-area projection. Cells coarser than the grid are
#'   subdivided first. The grid starts at the lower-left edge of the
#'   suitable cells and is shifted (`shifts` positions per axis); the
#'   smallest count is kept, as ConR does for locations.
#' * `locations`: the number of cells of a `threat_km` grid (the extent of a
#'   single threatening event; e.g. 10 km, Dauby et al. 2017) containing
#'   suitable area, also minimised over grid origins.
#' * Landscape structure (Curd et al. 2023): `n_patches` (8-neighbour
#'   connected suitable cells), `mean_patch_km2`, `median_patch_km2`,
#'   `largest_patch_km2`, `largest_patch_index` (largest patch as a
#'   percentage of the landscape, all cells with values), `edge_density`
#'   (m of edge between suitable and unsuitable cells per ha of
#'   landscape, as in FRAGSTATS), `mean_nn_km` (mean distance from each
#'   patch to the nearest other patch, edge to edge) and, for continuous
#'   input, `mean_patch_suitability`.
#'
#' **Threshold sensitivity.** AOO and EOO under different binarisation
#' rules can cross Red List thresholds, so the category a map supports may
#' depend on the rule rather than on the species. The printout reports the
#' range of each metric over the rules and whether the criterion-B
#' thresholds they meet differ.
#'
#' **Occupied patches.** The predicted suitable area is not the occupied
#' area: suitable but unoccupied areas inflate EOO and AOO. With
#' `occupied = TRUE` only patches containing at least one record
#' (`occurrences`) are kept, closer to the IUCN practice of estimating
#' occupancy from known, inferred or projected sites.
#'
#' **Red List thresholds** (IUCN 2012, criterion B): EOO below 100, 5000
#' and 20000 km^2 (CR, EN, VU; B1), AOO below 10, 500 and 2000 km^2 (B2),
#' and 1, at most 5 and at most 10 locations (condition a). Meeting a
#' threshold is not a listing: criterion B also needs at least two of the
#' conditions (a) to (c), and an assessment needs evidence on threats and
#' trends that a model cannot give.
#'
#' **Projected changes.** With several projections, the change of every
#' metric from the `baseline` is reported. Species distribution models
#' foresaw fewer than a third of the extirpations observed over a century
#' of climate and land-use change, missed many colonisations, and predicted
#' extirpations that had not happened where change was largest (Auffret et
#' al. 2024): projected range changes are uncertain in both directions.
#'
#' Memory: the suitable cells are read (with their coordinates); for
#' large rasters at fine resolution aggregate first, e.g. to 1 km.
#'
#' @param x Binary or continuous prediction rasters ([terra::SpatRaster]
#'   or file paths): one layer per projection.
#' @param thresholds For continuous input: a [thresholds()] result or a
#'   named numeric vector (names = rules).
#' @param model For a [thresholds()] result of several models, the model to
#'   use.
#' @param occurrences Optional records (data.frame with coordinates, `sf`,
#'   SpatVector or CSV path), for `occupied`.
#' @param occupied Keep only patches containing at least one record.
#' @param aoo_km Side of the AOO grid (km).
#' @param threat_km Side of the grid of locations (km).
#' @param shifts Number of grid origins tried per axis.
#' @param baseline The projection (layer name or number) changes are
#'   computed from.
#' @param equal_area Optional equal-area CRS for EOO, AOO and locations
#'   (default: a Lambert azimuthal equal-area projection centred on the
#'   raster).
#' @param rescale For continuous input, bring the predictions to the 0-1
#'   scale of the thresholds (0-100 or 0-1000 are divided; see
#'   [extract_predictions()]).
#' @inheritParams extract_predictions
#' @return A list of class `range_metrics`: `metrics` (one row per
#'   projection and rule), `change` (with several projections: absolute and
#'   relative change of each metric from the baseline, per rule), `patches`
#'   (area of every patch), `categories` (the criterion-B thresholds met by
#'   EOO, AOO and locations) and `settings`.
#' @references Curd, A., Chevalier, M., Vasquez, M., Boye, A., Firth, L.B.,
#'   Marzloff, M.P., Bricheno, L.M., Burrows, M.T., Bush, L.E., Cordier, C.,
#'   Davies, A.J., Green, J.A.M., Hawkins, S.J., Lima, F.P., Meneghesso, C.,
#'   Mieszkowska, N., Seabra, R., Dubois, S.F. (2023) Applying landscape
#'   metrics to species distribution model predictions to characterize
#'   internal range structure and associated changes. *Global Change
#'   Biology* 29, 631-647. \doi{10.1111/gcb.16496}
#'
#'   Dauby, G., Stevart, T., Droissart, V., Cosiaux, A., Deblauwe, V.,
#'   Simo-Droissart, M., Sosef, M.S.M., Lowry, P.P., Schatz, G.E.,
#'   Gereau, R.E., Couvreur, T.L.P. (2017) ConR: an R package to assist
#'   large-scale multispecies preliminary conservation assessments using
#'   distribution data. *Ecology and Evolution* 7, 11292-11303.
#'   \doi{10.1002/ece3.3704}
#'
#'   Omar, K., Elgamal, I. (2021) IUCN Red List and species distribution
#'   models as tools for the conservation of poorly known species: a case
#'   study of endemic plants *Micromeria serbaliana* and *Veronica kaiseri*
#'   in South Sinai, Egypt. *Kew Bulletin* 76, 477-496.
#'   \doi{10.1007/s12225-021-09953-4}
#'
#'   Auffret, A.G., Nenzen, H.K., Polaina, E. (2024) Underprediction of
#'   extirpation and colonisation following climate and land-use change
#'   using species distribution models. *Diversity and Distributions* 30,
#'   e13834. \doi{10.1111/ddi.13834}
#'
#'   IUCN (2012) *IUCN Red List Categories and Criteria: Version 3.1*.
#'   Second edition. IUCN, Gland and Cambridge.
#' @seealso [thresholds()], [write_range_metrics()],
#'   [plot_range_metrics()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models[["GLM"]], sim$presence,
#'                            sim$absence)
#' th <- thresholds(tab, rules = c("maxSSS", "p10", "maxMCC"))
#' rm <- range_metrics(sim$models[["GLM"]], th, threat_km = 50)
#' rm
range_metrics <- function(x, thresholds = NULL, model = NULL,
                          occurrences = NULL, occupied = FALSE,
                          aoo_km = 2, threat_km = 10, shifts = 4,
                          baseline = 1, equal_area = NULL, rescale = TRUE,
                          coords = NULL, crs = "EPSG:4326") {
  r <- as_grid(x)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") {
    e <- as.vector(terra::ext(r))
    if (all(abs(e[1:2]) <= 180) && all(abs(e[3:4]) <= 90)) {
      warning("The raster has no CRS; its coordinates look like degrees, ",
              "so EPSG:4326 is assumed. Set the CRS (terra::crs(x) <- ...) ",
              "to be sure.", call. = FALSE)
      terra::crs(r) <- "EPSG:4326"
    } else {
      warning("The raster has no CRS; coordinates are taken as metres in ",
              "an equal-area projection. Set the CRS to be sure.",
              call. = FALSE)
    }
  }
  lonlat <- isTRUE(terra::is.lonlat(r))
  if (lonlat && is.null(equal_area)) {
    message("The raster is in longitude/latitude: EOO, AOO and locations ",
            "are computed in a Lambert azimuthal equal-area projection ",
            "centred on it.")
  }
  if (is.null(equal_area)) {
    e <- as.vector(terra::ext(if (lonlat) r else terra::project(
      terra::as.polygons(terra::ext(r), crs = terra::crs(r)), "EPSG:4326")))
    equal_area <- sprintf(paste0("+proj=laea +lat_0=%.4f +lon_0=%.4f ",
                                 "+x_0=0 +y_0=0 +datum=WGS84 +units=m ",
                                 "+no_defs"), mean(e[3:4]), mean(e[1:2]))
  }
  if (!lonlat && !(is.na(terra::crs(r)) || terra::crs(r) == "") &&
      !is.null(equal_area)) {
    # projected input: used as it is for areas; the equal-area CRS only
    # for EOO/AOO/locations
  }
  proj <- names(r)
  if (anyDuplicated(proj)) proj <- make.unique(proj)
  names(r) <- proj
  b0 <- if (is.character(baseline)) match(baseline, proj) else baseline
  if (is.na(b0) || b0 < 1 || b0 > length(proj)) {
    stop("`baseline` is not a layer of `x`.", call. = FALSE)
  }

  # rules
  cont <- !is.null(thresholds)
  if (cont) {
    if (is.data.frame(thresholds)) {
      th <- as.data.frame(thresholds)
      if ("model" %in% names(th)) {
        mods <- unique(th$model)
        if (length(mods) > 1) {
          if (is.null(model)) {
            stop("`thresholds` has several models (",
                 paste(mods, collapse = ", "), "); choose one with `model`.",
                 call. = FALSE)
          }
          th <- th[th$model == model, , drop = FALSE]
        }
      }
      rules <- stats::setNames(th$threshold, th$rule)
    } else if (is.numeric(thresholds) && !is.null(names(thresholds))) {
      rules <- thresholds
    } else {
      stop("`thresholds` must come from thresholds() or be a named vector.",
           call. = FALSE)
    }
    rules <- rules[!is.na(rules)]
    if (!length(rules)) stop("No threshold.", call. = FALSE)
    smp <- sample_complete(r, 20000)
    h <- harmonise_scales(as.data.frame(smp), rescale, "predictions")
    if (any(h$divisor != 1)) {
      r <- r / unname(h$divisor[proj])
      names(r) <- proj
    }
  } else {
    smp <- sample_complete(r, 20000)
    if (!all(smp %in% c(0, 1))) {
      stop("Without `thresholds`, `x` must be binary (0/1). Give ",
           "thresholds for continuous predictions.", call. = FALSE)
    }
    rules <- c(binary = 0.5)
  }

  occ_cells <- NULL
  if (occupied) {
    if (is.null(occurrences)) stop("`occupied = TRUE` needs `occurrences`.",
                                   call. = FALSE)
    p <- as_points(occurrences, coords, crs)
    if (!terra::same.crs(p, r)) p <- terra::project(p, r)
    occ_cells <- unique(stats::na.omit(terra::cellFromXY(r,
                                                         terra::crds(p))))
  }

  land_km2 <- terra::global(terra::cellSize(r[[1]], mask = TRUE, unit = "km"),
                            "sum", na.rm = TRUE)[1, 1]
  rows <- list()
  patch_rows <- list()
  for (i in seq_along(proj)) {
    for (k in seq_along(rules)) {
      b <- r[[i]] >= rules[[k]]
      m <- binary_range_metrics(b, if (cont) r[[i]] else NULL, equal_area,
                                aoo_km, threat_km, shifts, occ_cells,
                                land_km2, lonlat)
      rows[[length(rows) + 1]] <- cbind(
        data.frame(projection = proj[i], rule = names(rules)[k],
                   threshold = if (cont) unname(rules[[k]]) else NA_real_,
                   stringsAsFactors = FALSE), m$metrics)
      if (nrow(m$patches)) {
        patch_rows[[length(patch_rows) + 1]] <- cbind(
          projection = proj[i], rule = names(rules)[k], m$patches,
          stringsAsFactors = FALSE)
      }
    }
  }
  met <- do.call(rbind, rows)
  rownames(met) <- NULL
  cats <- data.frame(projection = met$projection, rule = met$rule,
                     B1_EOO = iucn_class(met$EOO_km2, c(100, 5000, 20000)),
                     B2_AOO = iucn_class(met$AOO_km2, c(10, 500, 2000)),
                     a_locations = iucn_class(met$locations, c(1.5, 5.5,
                                                               10.5)),
                     stringsAsFactors = FALSE)
  # no suitable area is not a category
  cats[met$area_km2 == 0, c("B1_EOO", "B2_AOO", "a_locations")] <-
    "no suitable area"
  change <- NULL
  if (length(proj) > 1) {
    vars <- c("area_km2", "EOO_km2", "AOO_km2", "locations", "n_patches",
              "mean_patch_km2", "largest_patch_index", "edge_density")
    change <- do.call(rbind, lapply(setdiff(proj, proj[b0]), function(pj) {
      do.call(rbind, lapply(names(rules), function(ru) {
        a <- met[met$projection == proj[b0] & met$rule == ru, vars]
        z <- met[met$projection == pj & met$rule == ru, vars]
        d <- unlist(z) - unlist(a)
        rel <- ifelse(unlist(a) == 0, NA_real_, d / unlist(a))
        data.frame(projection = pj, rule = ru, metric = vars,
                   baseline = unlist(a), value = unlist(z), change = d,
                   relative = rel, stringsAsFactors = FALSE,
                   row.names = NULL)
      }))
    }))
  }
  structure(list(metrics = met, change = change,
                 patches = if (length(patch_rows)) do.call(rbind, patch_rows)
                   else NULL,
                 categories = cats,
                 settings = list(aoo_km = aoo_km, threat_km = threat_km,
                                 occupied = occupied, baseline = proj[b0],
                                 continuous = cont, rules = rules,
                                 landscape_km2 = land_km2,
                                 equal_area = equal_area)),
            class = "range_metrics")
}

# Criterion-B class of values (internal): below the 1st limit CR, 2nd EN,
# 3rd VU, else "none"
iucn_class <- function(v, lim) {
  out <- ifelse(v < lim[1], "CR", ifelse(v < lim[2], "EN",
                                         ifelse(v < lim[3], "VU", "none")))
  out[is.na(v)] <- NA_character_
  out
}

# Metrics of one binary map (internal)
binary_range_metrics <- function(b, cont, ea, aoo_km, threat_km, shifts,
                                 occ_cells, land_km2, lonlat) {
  empty <- data.frame(area_km2 = 0, EOO_km2 = 0, AOO_km2 = 0, locations = 0,
                      n_patches = 0, mean_patch_km2 = NA_real_,
                      median_patch_km2 = NA_real_,
                      largest_patch_km2 = NA_real_,
                      largest_patch_index = 0, edge_density = 0,
                      mean_nn_km = NA_real_,
                      mean_patch_suitability = NA_real_)
  pt <- terra::patches(b, directions = 8, zeroAsNA = TRUE)
  cells <- as.vector(terra::cells(pt))
  if (!length(cells)) {
    return(list(metrics = empty, patches = data.frame()))
  }
  pid <- terra::extract(pt, cells)[[1]]
  if (!is.null(occ_cells)) {
    keep <- unique(stats::na.omit(terra::extract(pt, occ_cells)[[1]]))
    sel <- pid %in% keep
    cells <- cells[sel]
    pid <- pid[sel]
    if (!length(cells)) return(list(metrics = empty, patches = data.frame()))
    # the map of the kept patches, for edges
    b <- terra::ifel(is.na(b), NA, 0)
    v <- terra::values(b, mat = FALSE)
    v[cells] <- 1
    b <- terra::setValues(b, v)
  }
  cs <- terra::extract(terra::cellSize(b, unit = "km"), cells)[[1]]
  area <- sum(cs)
  pa <- tapply(cs, pid, sum)

  # coordinates in the equal-area projection (km)
  xy <- terra::xyFromCell(b, cells)
  pxy <- terra::crds(terra::project(terra::vect(xy, crs = terra::crs(b)),
                                    ea)) / 1000
  side <- sqrt(cs)
  # EOO: hull of the cells (their corners)
  eoo <- if (length(cells) >= 1) {
    hv <- if (nrow(pxy) >= 3) grDevices::chull(pxy) else seq_len(nrow(pxy))
    cx <- rep(pxy[hv, 1], each = 4) + c(-1, 1, 1, -1) * rep(side[hv], each = 4) / 2
    cy <- rep(pxy[hv, 2], each = 4) + c(-1, -1, 1, 1) * rep(side[hv], each = 4) / 2
    h2 <- grDevices::chull(cx, cy)
    px <- cx[h2]
    py <- cy[h2]
    abs(sum(px * c(py[-1], py[1]) - c(px[-1], px[1]) * py)) / 2
  } else 0
  # points inside every cell, spaced at most half the AOO grid
  k <- pmax(1L, ceiling(side / (aoo_km / 2)))
  if (sum(k^2) > 2e7) {
    stop("Too many sub-cells for the AOO grid; aggregate the raster ",
         "first.", call. = FALSE)
  }
  sub <- lapply(seq_along(k), function(j) {
    o <- ((seq_len(k[j]) - 0.5) / k[j] - 0.5) * side[j]
    cbind(pxy[j, 1] + rep(o, k[j]), pxy[j, 2] + rep(o, each = k[j]))
  })
  sp <- do.call(rbind, sub)
  # grid origins: the lower-left edge of the suitable cells, then shifted
  # by fractions of the grid (the first try is aligned with the range)
  x0 <- min(pxy[, 1] - side / 2)
  y0 <- min(pxy[, 2] - side / 2)
  grid_count <- function(g) {
    off <- seq(0, g, length.out = shifts + 1)[seq_len(shifts)]
    best <- Inf
    for (ox in x0 + off) for (oy in y0 + off) {
      ix <- floor((sp[, 1] - ox) / g)
      iy <- floor((sp[, 2] - oy) / g)
      key <- sort.int((ix - min(ix)) * (max(iy) - min(iy) + 1) +
                        (iy - min(iy)), method = "radix")
      best <- min(best, 1 + sum(key[-1] != key[-length(key)]))
    }
    best
  }
  aoo <- grid_count(aoo_km) * aoo_km^2
  loc <- grid_count(threat_km)
  if (eoo < aoo) eoo <- aoo                       # IUCN guidelines

  # edges between suitable and unsuitable cells (inside the landscape)
  res <- terra::res(b)
  wh <- matrix(c(0, 0, 0, 0, -1, 1, 0, 0, 0), 3, byrow = TRUE)
  wv <- matrix(c(0, 0, 0, 0, -1, 0, 0, 1, 0), 3, byrow = TRUE)
  eh <- abs(terra::focal(b, wh, fun = "sum", na.rm = FALSE)) == 1
  ev <- abs(terra::focal(b, wv, fun = "sum", na.rm = FALSE)) == 1
  if (lonlat) {
    h_km <- res[2] * 110.574
    w_km <- terra::init(b, "y")
    w_km <- res[1] * 111.320 * cos((w_km - res[2] / 2) * pi / 180)
    edge_km <- terra::global(eh, "sum", na.rm = TRUE)[1, 1] * h_km +
      terra::global(ev * w_km, "sum", na.rm = TRUE)[1, 1]
  } else {
    u <- terra::linearUnits(b)
    if (!is.finite(u) || u <= 0) u <- 1
    edge_km <- (terra::global(eh, "sum", na.rm = TRUE)[1, 1] * res[2] +
                  terra::global(ev, "sum", na.rm = TRUE)[1, 1] * res[1]) *
      u / 1000
  }
  ed <- if (land_km2 > 0) edge_km * 1000 / (land_km2 * 100) else NA_real_

  # nearest other patch, edge to edge (km): only boundary cells of the
  # patches can be nearest to another patch
  nn <- NA_real_
  ids <- sort(unique(pid))
  if (length(ids) > 1) {
    # a suitable cell is on a boundary when one of its 4 neighbours is not
    # suitable (or outside the raster)
    rc <- terra::rowColFromCell(pt, cells)
    nr <- terra::nrow(pt)
    nc <- terra::ncol(pt)
    nb <- function(dr, dc) {
      rr <- rc[, 1] + dr
      cc <- rc[, 2] + dc
      inside <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
      out <- rep(FALSE, length(cells))
      out[inside] <- ((rr[inside] - 1) * nc + cc[inside]) %in% cells
      out
    }
    bi <- which(!(nb(-1, 0) & nb(1, 0) & nb(0, -1) & nb(0, 1)))
    # every patch keeps at least one cell
    bi <- sort(union(bi, which(!duplicated(pid))))
    bxy <- pxy[bi, , drop = FALSE]
    bp <- pid[bi]
    kk <- min(30, nrow(bxy) - 1)
    q <- dbscan::kNN(bxy, k = kk)
    by_patch <- split(seq_along(bp), bp)
    d_min <- vapply(as.character(ids), function(id) {
      i <- by_patch[[id]]
      if (!length(i)) return(NA_real_)
      other <- bp[q$id[i, , drop = FALSE]] != as.numeric(id)
      dd <- q$dist[i, , drop = FALSE][other]
      if (length(dd)) return(min(dd))
      # no other patch among the neighbours: search the other patches
      j <- which(bp != as.numeric(id))
      if (length(j) < 2) {
        return(min(sqrt(outer(bxy[i, 1], bxy[j, 1], "-")^2 +
                          outer(bxy[i, 2], bxy[j, 2], "-")^2)))
      }
      min(dbscan::kNN(bxy[j, , drop = FALSE], k = 1,
                      query = bxy[i, , drop = FALSE])$dist)
    }, numeric(1))
    # centre to centre minus one cell side ~ edge to edge
    nn <- mean(pmax(0, d_min - stats::median(side)), na.rm = TRUE)
  }
  ms <- NA_real_
  if (!is.null(cont)) {
    sv <- terra::extract(cont, cells)[[1]]
    ms <- mean(tapply(sv, pid, mean))
  }
  met <- data.frame(area_km2 = area, EOO_km2 = eoo, AOO_km2 = aoo,
                    locations = loc, n_patches = length(pa),
                    mean_patch_km2 = mean(pa),
                    median_patch_km2 = stats::median(pa),
                    largest_patch_km2 = max(pa),
                    largest_patch_index = 100 * max(pa) / land_km2,
                    edge_density = ed, mean_nn_km = nn,
                    mean_patch_suitability = ms)
  list(metrics = met,
       patches = data.frame(patch = as.integer(names(pa)),
                            area_km2 = as.vector(pa)))
}

#' @export
print.range_metrics <- function(x, digits = 1, ...) {
  m <- x$metrics
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(round(v), format = "d", big.mark = " ")
  nr <- length(s$rules)
  cat("Range metrics: ", length(unique(m$projection)), " projection(s) x ",
      nr, " threshold rule(s)", if (s$occupied) "; occupied patches only",
      "\n\nIN SHORT\n", sep = "")
  for (pj in unique(m$projection)) {
    t <- m[m$projection == pj, , drop = FALSE]
    cg <- x$categories[x$categories$projection == pj, , drop = FALSE]
    rng <- function(v, unit = " km2") {
      if (nr == 1 || diff(range(v)) == 0) paste0(f(v[1]), unit) else
        paste0(f(min(v)), unit, " (", t$rule[which.min(v)], ") to ",
               f(max(v)), unit, " (", t$rule[which.max(v)], ")")
    }
    say(paste0(pj, ": suitable area ", rng(t$area_km2), "; EOO ",
               rng(t$EOO_km2), "; AOO ", rng(t$AOO_km2), "; ",
               rng(t$locations, " locations"), "; ",
               if (length(unique(t$n_patches)) == 1)
                 paste0(t$n_patches[1], " patch", if (t$n_patches[1] != 1)
                   "es") else paste0(min(t$n_patches), " to ",
                                     max(t$n_patches), " patches"),
               "."))
    for (crit in c("B1_EOO", "B2_AOO", "a_locations")) {
      v <- cg[[crit]]
      lab <- c(B1_EOO = "EOO (B1)", B2_AOO = "AOO (B2)",
               a_locations = "the number of locations (condition a)")[crit]
      if (length(unique(v)) > 1) {
        say(paste0("The threshold rule changes the Red List threshold met ",
                   "by ", lab, ": ", paste0(t$rule, " ", v, collapse = ", "),
                   ". The category would depend on the rule, not on the ",
                   "species; report the range."))
      }
    }
  }
  if (!is.null(x$change)) {
    ch <- x$change[x$change$metric %in% c("area_km2", "AOO_km2"), ]
    for (pj in unique(ch$projection)) {
      a <- ch[ch$projection == pj & ch$metric == "area_km2", ]
      say(paste0("From ", s$baseline, " to ", pj, ", the suitable area ",
                 "changes by ", if (nrow(a) == 1 || diff(range(a$relative)) == 0)
                   share_signed(a$relative[1]) else paste0(
                     share_signed(min(a$relative)), " to ",
                     share_signed(max(a$relative)), " over the rules"), "."))
    }
    say(paste0("Projected changes are uncertain in both directions: over a ",
               "century of climate and land-use change, models foresaw ",
               "fewer than a third of the observed extirpations, missed ",
               "many colonisations, and predicted extirpations that had not ",
               "happened where change was largest (Auffret et al. 2024)."))
  }
  if (!s$occupied) {
    say("EOO and AOO are computed from the whole predicted suitable area, which is larger than the occupied area; occupied = TRUE keeps only patches with records, closer to IUCN practice.")
  }
  cat("\n")
  show <- m[c("projection", "rule", "threshold", "area_km2", "EOO_km2",
              "AOO_km2", "locations", "n_patches", "largest_patch_index",
              "edge_density", "mean_nn_km")]
  if (!s$continuous) show$threshold <- NULL
  if (length(unique(show$projection)) == 1) show$projection <- NULL
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "EOO = minimum convex polygon of the suitable cells; AOO = ",
    s$aoo_km, " x ", s$aoo_km, " km grid cells with suitable area; ",
    "locations = ", s$threat_km, " km grid cells (the extent of one ",
    "threat); both grids are placed to give the fewest cells, in an ",
    "equal-area projection. largest_patch_index = largest patch, % of the ",
    "landscape; edge_density = m of edge per ha of landscape; mean_nn_km = ",
    "mean distance to the nearest other patch. IUCN (2012) criterion B: EOO ",
    "< 100 / 5000 / 20000 km2 and AOO < 10 / 500 / 2000 km2 (CR / EN / ",
    "VU), locations 1 / <= 5 / <= 10; a threshold met is not a listing, ",
    "which also needs conditions (a)-(c) and evidence on threats."),
    width = 78), sep = "\n")
  invisible(x)
}

# Signed percentage (internal)
share_signed <- function(v) {
  if (is.na(v)) return("NA")
  paste0(if (v > 0) "+" else "", round(100 * v), "%")
}

#' Plot range metrics under each threshold rule
#'
#' One metric (by default AOO) under each threshold rule, per projection,
#' with the IUCN criterion-B thresholds as dashed lines, so you can see
#' whether the rule moves the range across a Red List threshold.
#'
#' @param x A result of [range_metrics()].
#' @param metric `"AOO_km2"`, `"EOO_km2"`, `"locations"`, `"area_km2"` or
#'   another column of `x$metrics`.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references IUCN (2012) *IUCN Red List Categories and Criteria: Version
#'   3.1*. Second edition. IUCN, Gland and Cambridge.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models[["GLM"]], sim$presence,
#'                            sim$absence)
#' th <- thresholds(tab, rules = c("maxSSS", "p10", "maxMCC"))
#' plot_range_metrics(range_metrics(sim$models[["GLM"]], th))
plot_range_metrics <- function(x, metric = "AOO_km2", base_size = 11,
                               file = NULL, width = 7, height = 4.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "range_metrics")) {
    stop("`x` must come from range_metrics().", call. = FALSE)
  }
  m <- x$metrics
  if (!metric %in% names(m)) stop("Unknown metric: ", metric, call. = FALSE)
  m$rule <- factor(m$rule, unique(m$rule))
  m$projection <- factor(m$projection, unique(m$projection))
  lim <- switch(metric, AOO_km2 = c(CR = 10, EN = 500, VU = 2000),
                EOO_km2 = c(CR = 100, EN = 5000, VU = 20000),
                locations = c(CR = 1, EN = 5, VU = 10), NULL)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02")
  np <- length(unique(m$projection))
  pal <- if (np <= 6) pal[seq_len(np)] else
    grDevices::hcl.colors(np, "Dark 3")
  lab <- c(AOO_km2 = "AOO (km2)", EOO_km2 = "EOO (km2)",
           locations = "Number of locations", area_km2 = "Suitable area (km2)")
  p <- ggplot2::ggplot(m, ggplot2::aes(.data$rule, .data[[metric]],
                                       fill = .data$projection)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                      width = 0.7) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::labs(x = "Threshold rule",
                  y = if (metric %in% names(lab)) lab[[metric]] else metric,
                  title = paste0(if (metric %in% names(lab)) lab[[metric]]
                                 else metric, " under each threshold rule"),
                  subtitle = if (!is.null(lim))
                    "Dashed: IUCN criterion B thresholds within the plotted range." else
                      NULL) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
  p <- p + ggplot2::scale_y_continuous(labels = function(v) {
    formatC(v, format = "d", big.mark = " ")
  })
  if (!is.null(lim)) {
    top <- max(m[[metric]], na.rm = TRUE)
    # thresholds far below the values would pile up at the bottom
    lim <- lim[lim <= top * 1.5 & lim >= top / 40]
    if (!length(lim)) {
      p <- p + ggplot2::labs(subtitle = paste0(
        "All IUCN criterion B thresholds (up to ",
        formatC(max(switch(metric, AOO_km2 = 2000, EOO_km2 = 20000,
                           locations = 10)), format = "d", big.mark = " "),
        ") lie far below these values."))
    }
    if (length(lim)) {
      p <- p + ggplot2::geom_hline(yintercept = lim, linetype = "dashed",
                                   colour = "#8a8a86") +
        ggplot2::annotate("text", x = 0.5, y = lim, label = names(lim),
                          hjust = 0, vjust = -0.4, colour = "#5c5c58",
                          size = base_size / 3.4)
    }
  }
  if (np == 1) p <- p + ggplot2::guides(fill = "none")
  p
}

#' Write range metrics to files
#'
#' Writes `<prefix>.csv` (metrics per projection and rule),
#' `<prefix>_categories.csv` (criterion-B thresholds met),
#' `<prefix>_change.csv` (with several projections), `<prefix>_patches.csv`
#' (area of every patch) and `<prefix>_AOO.png`.
#'
#' @param x A `range_metrics` object.
#' @param dir Output folder.
#' @param prefix File name prefix.
#' @return The paths of the written files, invisibly.
#' @references IUCN (2012) *IUCN Red List Categories and Criteria: Version
#'   3.1*. Second edition. IUCN, Gland and Cambridge.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' tab <- extract_predictions(sim$models[["GLM"]], sim$presence,
#'                            sim$absence)
#' rm <- range_metrics(sim$models[["GLM"]], thresholds(tab, rules = "p10"))
#' write_range_metrics(rm, tempdir())
write_range_metrics <- function(x, dir, prefix = "range_metrics") {
  if (!inherits(x, "range_metrics")) {
    stop("`x` must come from range_metrics().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, s))
  out <- c(metrics = f(".csv"), categories = f("_categories.csv"))
  utils::write.csv(x$metrics, out[["metrics"]], row.names = FALSE)
  utils::write.csv(x$categories, out[["categories"]], row.names = FALSE)
  if (!is.null(x$change)) {
    out <- c(out, change = f("_change.csv"))
    utils::write.csv(x$change, out[["change"]], row.names = FALSE)
  }
  if (!is.null(x$patches)) {
    out <- c(out, patches = f("_patches.csv"))
    utils::write.csv(x$patches, out[["patches"]], row.names = FALSE)
  }
  out <- c(out, plot = f("_AOO.png"))
  ggplot2::ggsave(out[["plot"]], plot_range_metrics(x), width = 7,
                  height = 4.5, dpi = 150)
  invisible(out)
}
