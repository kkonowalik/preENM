# One point per raster cell and background points ------------------------------

#' Keep one record per raster cell, centred in the cell
#'
#' Reduces presence records to one per cell of the rasters you will model with
#' (per species), and moves each kept record to the centre of its cell. This
#' removes the sampling bias of many records in the same cell and makes every
#' record match exactly one cell of the predictors. Records outside the
#' rasters, or in cells where *any* layer has no value, are dropped.
#'
#' When several records fall in one cell, the first one is kept. With
#' `priority`, records from preferred sources win (e.g. your own verified
#' samples over GBIF records); otherwise the input order decides.
#'
#' This is geographic thinning to the raster resolution. It is not the same
#' as environmental thinning (thinning in climate space).
#'
#' @param x Occurrence records: a data.frame with coordinates (e.g. from
#'   [merge_occurrences()] or [kept()]), a CSV path, an `sf` object or a
#'   [terra::SpatVector].
#' @param raster The predictor layers: a [terra::SpatRaster] (several layers
#'   allowed) or file paths. They define the grid; all layers must share
#'   extent, resolution and CRS.
#' @param by Column with species names: one record per cell *and species*.
#'   Ignored if absent.
#' @param priority Optional character vector of values of `priority_col`,
#'   from most to least preferred, e.g. `c("sNMF", "GBIF")`.
#' @param priority_col Column holding the source (default `"source"`, as
#'   created by [merge_occurrences()]).
#' @param centre Move kept records to the centre of their cell (default
#'   `TRUE`); the original coordinates are kept in `lon_orig`, `lat_orig`.
#' @param crs CRS of the output coordinates (default WGS84 longitude/latitude).
#' @param coords Optional names of the x/longitude and y/latitude columns of
#'   `x`; detected automatically if `NULL`.
#' @param points_crs CRS of the coordinates in `x` (default WGS84
#'   longitude/latitude).
#' @return A data.frame of the kept records with `lon`, `lat` (cell centres),
#'   `cell` (cell number), `lon_orig`, `lat_orig` and the other columns of
#'   `x`.
#' @export
#' @examples
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 40,
#'                  ymax = 50, vals = 1)
#' pts <- data.frame(lon = c(0.2, 0.7, 5.5, 20), lat = c(40.2, 40.6, 45.5, 45),
#'                   source = c("GBIF", "own", "GBIF", "GBIF"))
#' thin_by_cell(pts, r, priority = c("own", "GBIF"))
thin_by_cell <- function(x, raster, by = "species", priority = NULL,
                         priority_col = "source", centre = TRUE,
                         crs = "EPSG:4326", coords = NULL,
                         points_crs = "EPSG:4326") {
  r <- as_grid(raster)
  if (is.character(x) && length(x) == 1) x <- read_any_csv(x)
  if (is.data.frame(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
    lon_in <- as_number(df[[xy_cols[1]]])
    lat_in <- as_number(df[[xy_cols[2]]])
    has_xy <- !is.na(lon_in) & !is.na(lat_in)
    p <- terra::vect(cbind(lon_in[has_xy], lat_in[has_xy]), crs = points_crs)
  } else {
    p <- as_points(x, coords, points_crs)
    df <- as.data.frame(p)
    has_xy <- rep(TRUE, nrow(df))
  }
  # every vector below has one element per row of df
  n <- nrow(df)
  xy0 <- matrix(NA_real_, n, 2)
  cell <- rep(NA_real_, n)
  if (any(has_xy)) {
    xy0[has_xy, ] <- terra::crds(terra::project(p, "EPSG:4326"))
    pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
    cell[has_xy] <- terra::cellFromXY(r, terra::crds(pr))
  }
  inside <- !is.na(cell)
  ok <- rep(FALSE, n)
  if (any(inside)) {
    v <- terra::extract(r, cell[inside])
    ok[inside] <- stats::complete.cases(v)
  }
  message(if (any(!has_xy)) paste0("Without coordinates: ", sum(!has_xy),
                                   "; "),
          "Outside the raster: ", sum(has_xy & !inside),
          "; inside but without values in all layers: ", sum(inside & !ok),
          "; remaining: ", sum(ok))
  d <- df[ok, , drop = FALSE]
  d$cell <- cell[ok]
  d$lon_orig <- xy0[ok, 1]
  d$lat_orig <- xy0[ok, 2]
  ord <- seq_len(nrow(d))
  if (!is.null(priority)) {
    if (!priority_col %in% names(d)) {
      stop("`priority` given but column `", priority_col, "` not found.",
           call. = FALSE)
    }
    pr_rank <- match(d[[priority_col]], priority)
    pr_rank[is.na(pr_rank)] <- length(priority) + 1
    ord <- order(pr_rank, ord)
  }
  d <- d[ord, , drop = FALSE]
  grp <- if (!is.null(by) && by %in% names(d)) d[[by]] else ""
  keep <- !duplicated(paste(grp, d$cell))
  message("Records in a cell already taken (dropped): ", sum(!keep))
  d <- d[keep, , drop = FALSE]
  if (centre) {
    ctr <- terra::xyFromCell(r, d$cell)
    ctr <- cell_coords(ctr, r, crs)
    d$lon <- ctr[, 1]
    d$lat <- ctr[, 2]
    if (nrow(d)) {
      shift <- sqrt(((d$lon - d$lon_orig) * cos(d$lat * pi / 180))^2 +
                      (d$lat - d$lat_orig)^2) * 111.2
      if (terra::is.lonlat(terra::vect(cbind(0, 0), crs = crs))) {
        message("Largest shift to a cell centre: ~", round(max(shift), 2),
                " km")
      }
    }
  } else {
    d$lon <- d$lon_orig
    d$lat <- d$lat_orig
  }
  first <- intersect(c(by, "lon", "lat"), names(d))
  d <- d[c(first, setdiff(names(d), first))]
  rownames(d) <- NULL
  message("Thinned records: ", nrow(d))
  d
}

#' Background (pseudo-absence) points: uniform, environmental, effort-based
#'
#' Draws background points from cells of the predictor rasters, one point per
#' cell (without replacement), only from cells where *all* layers have
#' values. Several independent sets can be drawn (e.g. for the `set` argument
#' of [evaluate_cv()]); each uses its own seed (`seed + set number`).
#' Restrict the rasters to the accessible area (M) beforehand.
#'
#' **Methods**
#' * `"uniform"` (default): random cells; every cell equally likely.
#' * `"envstrat"`: uniform in environmental space. A PCA of the layers
#'   (precipitation amounts log-transformed) is cut into bins: all PCs up to
#'   `env_cumvar` of the variance, about `env_total_bins` bins in total
#'   (`env_axes = 2, env_total_bins = 1600` gives the earlier PC1 x PC2,
#'   40 x 40 version). Every occupied bin is equally likely, so rare climates
#'   are represented as well as common ones (bins with fewer than
#'   `env_min_bin` cells are down-weighted). It mirrors environmentally
#'   thinned presences ([thin_environment()]).
#' * `"effort"`: probability proportional to recording effort, e.g. the
#'   number of plant species recorded in GBIF per cell (see
#'   [gbif_effort_sql()]), or to the records of a target group of companion
#'   species (a "target-group background"). Corrects for the fact that "not
#'   recorded" means "absent" only where people recorded.
#' * `"enveffort"`: the `envstrat` bins with equal total probability per
#'   bin, but *within* a bin cells with more recording effort are preferred
#'   (their missing record of the species is a more reliable absence).
#' * `"target_group"`: probability proportional to the density of the
#'   records of a target group: species of the same genus, family or
#'   another group recorded by the same observers and methods, supplied as
#'   records in `target` (e.g. a GBIF download built with
#'   [gbif_target_group_sql()]). Target-group backgrounds corrected
#'   sampling bias best among four common methods, though they sometimes
#'   overcompensated (Barber et al. 2021). With `target_bandwidth_km = 0`
#'   the density is the number of target records per cell (the classic
#'   target-group background); with a bandwidth it is a Gaussian kernel
#'   density estimate, which spreads the effort to neighbouring cells
#'   (TGOB+ of Balej et al. 2025, their recommendation when no independent
#'   data are available; no single bandwidth was best, so try a few, e.g.
#'   1, 2 and 4 cell widths). With `target_observer` (the observer column
#'   of `target` and `presences`) each observer's density is normalised and
#'   weighted by that observer's share of the records of your species, as
#'   in the presence-weighted observer-oriented approach (PW-OOA) of Balej
#'   et al. (2025), which was best when checked against independent data.
#'   `target_stratify = TRUE` combines the density with the environmental
#'   bins of `"enveffort"`: every occupied bin keeps the same total
#'   probability and the target density decides within it (bins without
#'   target records are not sampled unless `target_floor > 0`).
#'
#' No method is best for every species and data set; compare the
#' backgrounds on the same hold-out folds (and independent records, if you
#' have them), with `"uniform"` as the reference.
#'
#' For `"envstrat"`, `"effort"` and `"enveffort"` the points are drawn from
#' a random pool of `env_pool` valid cells (all cells if fewer), by weighted
#' sampling without replacement.
#'
#' @param raster The predictor layers ([terra::SpatRaster] or file paths),
#'   clipped to the accessible area. For the environmental methods use the
#'   variables of your model.
#' @param n Number of points per set (default 10000). If fewer cells are
#'   available, all are used, with a warning.
#' @param n_sets Number of independent sets.
#' @param exclude Optional presence records whose cells are excluded: a
#'   data.frame with a `cell` column (as returned by [thin_by_cell()]) or
#'   with coordinates, or a vector of cell numbers.
#' @param seed Optional random seed.
#' @param crs CRS of the output coordinates (default WGS84 longitude/latitude).
#' @param method `"uniform"`, `"envstrat"`, `"effort"` or `"enveffort"`.
#' @param effort Recording-effort raster for `"effort"` and `"enveffort"`
#'   (any resolution and CRS; `NA` counts as 0), e.g. `effort_n_species`
#'   from [gbif_effort_rasters()].
#' @param effort_floor Added to every cell's effort, so cells without records
#'   keep a small chance.
#' @param effort_power Exponent of the effort (1 = proportional; 0.5 = a
#'   softer correction).
#' @param target For `"target_group"`: records of the target group (a
#'   data.frame with coordinates, a CSV path, `sf` or SpatVector; a column
#'   `n`, `count` or `n_target` gives the number of records at each point,
#'   e.g. from an aggregated GBIF table), or a raster of record counts.
#'   Leave your own species out of it.
#' @param target_bandwidth_km Standard deviation of the Gaussian kernel
#'   (km); 0 = records per cell.
#' @param target_floor Added to the target density of every cell (0: cells
#'   without target records are never drawn when the bandwidth is 0).
#' @param target_observer Name of the observer column (e.g. `recordedBy`)
#'   in `target` and `presences`, for observer weighting.
#' @param presences Records of your species, with the observer column, for
#'   the observer weights.
#' @param target_stratify Combine the target density with the environmental
#'   bins (as `"enveffort"` does with `effort`).
#' @param env_pool Size of the random pool of cells for the weighted methods.
#' @param env_cumvar,env_axes PCs used for the bins: those explaining
#'   `env_cumvar` of the variance, or exactly `env_axes` PCs.
#' @param env_total_bins Approximate total number of bins.
#' @param env_min_bin Bins with fewer pool cells are down-weighted.
#' @param log_vars Layers log-transformed for the PCA (default: bioclim
#'   precipitation amounts recognised by name).
#' @return A data.frame with `set`, `cell`, `lon`, `lat` (cell centres) and,
#'   for the environmental methods, `bin`. Attribute `"background"` holds the
#'   method and, for the environmental methods, the PCA and bin settings.
#' @references Barber, R.A., Ball, S.G., Morris, R.K.A., Gilbert, F. (2021)
#'   Target-group backgrounds prove effective at correcting sampling bias in
#'   Maxent models. *Diversity and Distributions* 28, 128-141.
#'   \doi{10.1111/ddi.13442}
#'
#'   Balej, P., Moudry, V., Prajzlerova, D., Gabor, L., Sillero, N.,
#'   Rocchini, D., Simova, P. (2025) Species-observer link and kernel
#'   density estimation of background points allow for sampling bias
#'   correction in bird species distribution models. *Ecography* 2025,
#'   e08202. \doi{10.1002/ecog.08202}
#' @export
#' @examples
#' r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10, ymin = 40,
#'                  ymax = 50, vals = 1)
#' bg <- sample_background(r, n = 50, n_sets = 3, seed = 1)
#' table(bg$set)
#'
#' # uniform in environmental space, with recording effort within bins
#' sim <- simulate_enm_data(seed = 1)
#' effort <- terra::init(sim$env[[1]], function(n) rpois(n, 5))
#' bg2 <- sample_background(sim$env, n = 500, method = "enveffort",
#'                          effort = effort, seed = 1)
#' head(bg2)
sample_background <- function(raster, n = 10000, n_sets = 1, exclude = NULL,
                              seed = NULL, crs = "EPSG:4326",
                              method = c("uniform", "envstrat", "effort",
                                         "enveffort", "target_group"),
                              effort = NULL, effort_floor = 1,
                              effort_power = 1, target = NULL,
                              target_bandwidth_km = 0, target_floor = 0,
                              target_observer = NULL, presences = NULL,
                              target_stratify = FALSE, env_pool = 1e6,
                              env_cumvar = 0.95, env_axes = NULL,
                              env_total_bins = 3000, env_min_bin = 5,
                              log_vars = NULL) {
  method <- match.arg(method)
  r <- as_grid(raster)
  if (method %in% c("effort", "enveffort") && is.null(effort)) {
    stop("method = \"", method, "\" needs an `effort` raster.", call. = FALSE)
  }
  if (method == "target_group" && is.null(target)) {
    stop("method = \"target_group\" needs `target`: records of the target ",
         "group (e.g. from gbif_target_group_sql()).", call. = FALSE)
  }
  ok <- if (terra::nlyr(r) > 1) sum(r) else r
  valid <- terra::cells(ok)
  if (is.list(valid)) valid <- valid[[1]]
  if (!is.null(exclude)) {
    ex <- if (is.numeric(exclude) && is.null(dim(exclude))) {
      exclude
    } else if (is.data.frame(exclude) && "cell" %in% names(exclude)) {
      exclude$cell
    } else {
      p <- as_points(exclude)
      if (!terra::same.crs(p, r)) p <- terra::project(p, r)
      terra::cellFromXY(r, terra::crds(p))
    }
    valid <- setdiff(valid, ex)
  }
  message("Cells available for background: ", length(valid))
  n_take <- min(n, length(valid))
  if (n_take < n) {
    warning("Fewer valid cells (", length(valid), ") than n; using all.",
            call. = FALSE)
  }
  info <- list(method = method)
  if (method == "uniform") {
    pool <- valid
    w <- NULL
    bin <- NULL
  } else {
    pool <- with_seed(seed, sort(valid[sample.int(length(valid),
                                                  min(length(valid),
                                                      env_pool))]))
    bin <- NULL
    if (method %in% c("envstrat", "enveffort") ||
        (method == "target_group" && target_stratify)) {
      eb <- env_bins(r, pool, env_cumvar, env_axes, env_total_bins, log_vars)
      pool <- pool[eb$ok]
      bin <- eb$bin
      cnt <- as.numeric(table(bin)[bin])
      info <- c(info, eb[c("n_axes", "bins_per_axis", "n_bins", "log_vars",
                           "variance")])
      message("Environmental bins: PC1-PC", eb$n_axes, ", ", eb$bins_per_axis,
              " bins per axis, ", eb$n_bins, " occupied")
    }
    if (method %in% c("effort", "enveffort")) {
      ef <- as_grid(effort)[[1]]
      e <- extract_xy(ef, terra::xyFromCell(r, pool), terra::crs(r))
      e[is.na(e)] <- 0
      e <- (pmax(e, 0) + effort_floor)^effort_power
      info$effort <- c(floor = effort_floor, power = effort_power)
    }
    if (method == "target_group") {
      td <- target_density(r, pool, target, target_bandwidth_km,
                           target_observer, presences)
      e <- (td$density + target_floor)^effort_power
      info$target <- td$info
      info$target$floor <- target_floor
      info$target$stratify <- target_stratify
      if (!target_stratify) bin <- NULL
    }
    stratified <- function(e) {
      tot <- stats::ave(e, bin, FUN = sum)
      ifelse(tot > 0, e / tot, 0) * pmin(1, cnt / env_min_bin)
    }
    w <- switch(method,
      envstrat = 1 / pmax(cnt, env_min_bin),
      effort = e,
      enveffort = stratified(e),
      target_group = if (target_stratify) stratified(e) else e)
    keep <- w > 0
    pool <- pool[keep]
    w <- w[keep]
    if (!is.null(bin)) bin <- bin[keep]
    n_take <- min(n_take, length(pool))
  }
  sets <- lapply(seq_len(n_sets), function(k) {
    s <- with_seed(if (is.null(seed)) NULL else seed + k, {
      if (is.null(w)) sample.int(length(pool), n_take) else
        weighted_draw(w, n_take)
    })
    cl <- pool[s]
    xy <- cell_coords(terra::xyFromCell(r, cl), r, crs)
    d <- data.frame(set = k, cell = cl, lon = xy[, 1], lat = xy[, 2])
    if (!is.null(bin)) d$bin <- bin[s]
    d
  })
  out <- do.call(rbind, sets)
  rownames(out) <- NULL
  attr(out, "background") <- info
  out
}

# Density of target-group records at pool cells (internal): records per
# cell (bandwidth 0) or a Gaussian kernel density in records per cell,
# optionally per observer, normalised and weighted by the observer's share
# of the focal species' records (PW-OOA of Balej et al. 2025)
target_density <- function(r, pool, target, bw_km, observer, presences) {
  tp <- target_points(target, r)
  xy <- tp$xy
  cnt <- tp$n
  obs <- NULL
  w_obs <- NULL
  if (!is.null(observer)) {
    obs_t <- if (is.null(tp$data)) NA else find_column(names(tp$data),
                                                       observer)
    if (is.na(obs_t)) {
      stop("Observer column `", observer, "` not found in `target`.",
           call. = FALSE)
    }
    if (is.null(presences)) {
      stop("Observer weighting needs `presences` (records of your ",
           "species with the observer column).", call. = FALSE)
    }
    pr <- as.data.frame(presences, stringsAsFactors = FALSE)
    obs_p <- find_column(names(pr), observer)
    if (is.na(obs_p)) {
      stop("Observer column `", observer, "` not found in `presences`.",
           call. = FALSE)
    }
    pp <- as_points(pr)
    if (!terra::same.crs(pp, r)) pp <- terra::project(pp, r)
    pc <- terra::cellFromXY(r, terra::crds(pp))
    # an observer's records in one cell count once (Balej et al. 2025)
    u <- !duplicated(paste(pr[[obs_p]], pc)) & !is.na(pc) &
      !is.na(pr[[obs_p]])
    tab <- table(pr[[obs_p]][u])
    w_obs <- tab / sum(tab)
    obs <- as.character(tp$data[[obs_t]])
  }
  geo <- terra::is.lonlat(r)
  pxy <- terra::xyFromCell(r, pool)
  res <- terra::res(r)
  cell_km2 <- if (geo) {
    (res[1] * 111.32 * cos(pxy[, 2] * pi / 180)) * (res[2] * 111.32)
  } else rep(prod(res) / 1e6 * (terra::linearUnits(r)^2), length(pool))
  group <- if (is.null(obs)) rep("all", nrow(xy)) else obs
  keep_g <- if (is.null(w_obs)) "all" else names(w_obs)
  dens <- numeric(length(pool))
  for (g in keep_g) {
    j <- which(group == g)
    if (!length(j)) next
    if (bw_km <= 0) {
      tc <- terra::cellFromXY(r, xy[j, , drop = FALSE])
      ok <- !is.na(tc)
      d <- tapply(cnt[j][ok], tc[ok], sum)
      v <- as.numeric(d[match(pool, as.numeric(names(d)))])
      v[is.na(v)] <- 0
    } else {
      v <- kernel_at(pxy, xy[j, , drop = FALSE], cnt[j], bw_km, geo) *
        cell_km2
    }
    if (!is.null(w_obs)) {
      s <- sum(v)
      v <- if (s > 0) v / s * as.numeric(w_obs[[g]]) else v
    }
    dens <- dens + v
  }
  if (!any(dens > 0)) {
    warning("No target-group record falls in the background cells: check ",
            "the extent and CRS of `target`.", call. = FALSE)
  }
  list(density = dens,
       info = list(n_records = sum(cnt), bandwidth_km = bw_km,
                   observers = if (is.null(w_obs)) NULL else length(w_obs)))
}

# Target-group records as coordinates in the raster's CRS with counts
# (internal)
target_points <- function(target, r) {
  if (inherits(target, "SpatRaster")) {
    t1 <- target[[1]]
    if (terra::ncell(t1) > 5e7) {
      stop("The target raster is very large; give the records as points.",
           call. = FALSE)
    }
    v <- terra::values(t1, mat = FALSE)
    i <- which(!is.na(v) & v > 0)
    xy <- terra::xyFromCell(t1, i)
    p <- terra::vect(xy, crs = terra::crs(t1))
    if (!terra::same.crs(p, r)) p <- terra::project(p, r)
    return(list(xy = terra::crds(p), n = v[i], data = NULL))
  }
  if (is.character(target) && length(target) == 1) {
    target <- read_any_csv(target)
  }
  d <- if (is.data.frame(target) && !inherits(target, "sf")) {
    as.data.frame(target, stringsAsFactors = FALSE)
  } else NULL
  if (!is.null(d) && all(c("lat_bin", "lon_bin") %in% tolower(names(d)))) {
    stop("An aggregated GBIF table: turn it into a raster with ",
         "gbif_effort_rasters() and give its count layer as `target`.",
         call. = FALSE)
  }
  if (!is.null(d)) {
    xy_cols <- guess_coords(names(d))
    d <- d[!is.na(as_number(d[[xy_cols[1]]])) &
             !is.na(as_number(d[[xy_cols[2]]])), , drop = FALSE]
  }
  p <- as_points(if (is.null(d)) target else d)
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  ncol_ <- if (is.null(d)) NA else find_column(names(d), c("n", "count",
                                                          "n_target"))
  n <- if (is.na(ncol_)) rep(1, nrow(terra::crds(p))) else
    pmax(as_number(d[[ncol_]]), 0)
  n[is.na(n)] <- 0
  list(xy = terra::crds(p), n = n, data = d)
}

# Gaussian kernel density (records per km2) at query points Q from points
# X with weights w and standard deviation h km, summed over neighbours
# within 4 h (internal); processed in chunks
kernel_at <- function(Q, X, w, h, geo) {
  to_space <- function(xy) {
    if (geo) {
      lon <- xy[, 1] * pi / 180
      lat <- xy[, 2] * pi / 180
      cbind(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat)) * 6371.0088
    } else xy / 1000
  }
  Xs <- to_space(X)
  Qs <- to_space(Q)
  out <- numeric(nrow(Qs))
  if (nrow(Xs) == 1) Xs <- rbind(Xs, Xs + 1e9)    # dbscan needs two
  w <- c(w, 0)[seq_len(nrow(Xs))]
  step <- 2e4
  for (s0 in seq(1, nrow(Qs), by = step)) {
    i <- s0:min(nrow(Qs), s0 + step - 1)
    nn <- dbscan::frNN(Xs, eps = 4 * h, query = Qs[i, , drop = FALSE],
                       sort = FALSE)
    out[i] <- vapply(seq_along(i), function(k) {
      id <- nn$id[[k]]
      if (!length(id)) return(0)
      sum(w[id] * exp(-nn$dist[[k]]^2 / (2 * h^2)))
    }, numeric(1))
  }
  out / (2 * pi * h^2)
}

# Weighted sampling without replacement (Efraimidis-Spirakis; fast for
# millions of cells) (internal)
weighted_draw <- function(w, n) order(stats::rexp(length(w)) / w)[seq_len(n)]

# Environmental bins of pool cells (internal)
env_bins <- function(r, cells, cumvar, n_axes, total_bins, log_vars) {
  M <- read_cells(r, cells)
  colnames(M) <- names(r)
  ok <- stats::complete.cases(M)
  M <- M[ok, , drop = FALSE]
  lv <- if (is.null(log_vars)) default_log_vars(names(r)) else
    intersect(log_vars, names(r))
  M <- log_tf(pmax_cols(M, lv), lv)
  pc <- stats::prcomp(M, scale. = TRUE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  k <- if (!is.null(n_axes)) min(n_axes, ncol(M)) else
    which(cumsum(ve) >= cumvar - 1e-9)[1]
  nb <- max(4, round(total_bins^(1 / k)))
  bin <- do.call(paste, lapply(seq_len(k), function(j) {
    cut(pc$x[, j], nb, labels = FALSE)
  }))
  list(ok = ok, bin = bin, n_axes = k, bins_per_axis = nb,
       n_bins = length(unique(bin)), log_vars = lv,
       variance = sum(ve[seq_len(k)]))
}

pmax_cols <- function(M, cols) {
  if (length(cols)) M[, cols] <- pmax(M[, cols], 0)
  M
}

# Raster values at x/y given in `crs` (internal)
extract_xy <- function(r, xy, crs) {
  p <- terra::vect(xy, crs = crs)
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  terra::extract(r, p, ID = FALSE)[[1]]
}

# Read predictor layers as one SpatRaster (internal)
as_grid <- function(raster) {
  if (is.character(raster)) raster <- terra::rast(raster)
  if (inherits(raster, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    raster <- terra::rast(raster)
  }
  if (!inherits(raster, "SpatRaster")) {
    stop("`raster` must be a SpatRaster or raster file path(s).",
         call. = FALSE)
  }
  raster
}

# Cell-centre coordinates in the output CRS (internal)
cell_coords <- function(xy, r, crs) {
  if (!nrow(xy)) return(xy)
  if (terra::same.crs(r, crs)) return(xy)
  terra::crds(terra::project(terra::vect(xy, crs = terra::crs(r)), crs))
}
