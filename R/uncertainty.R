# Coordinate uncertainty in environmental terms ----------------------------------

#' Coordinate uncertainty in environmental terms
#'
#' Asks, for every record, how much the predictors vary within its coordinate
#' uncertainty, instead of judging the uncertainty in metres. A record with
#' an uncertainty of 5 km on a plain has almost the same climate wherever
#' it truly lies and is safe to use; one with 1 km on a steep slope may not
#' be. Records are usually dropped above a fixed uncertainty (see
#' `max_coord_unc_m` in [filter_settings()]), which discards many usable
#' records, often the older ones and those from poorly sampled regions, and
#' so narrows the sampled niche (Smith et al. 2023; Marcer et al. 2022).
#'
#' For each record, the cells that its uncertainty circle touches are found
#' (the record's own cell always counts) and the range of each layer among
#' them is divided by the species' niche breadth on that layer: the
#' distance between the `breadth_probs` quantiles (default 5% and 95%) of
#' the values at all records. A relative range of 0.1 means that the true
#' value may differ by up to a tenth of the breadth of the species' niche
#' on that variable. `unc_ok` is `TRUE` when no layer exceeds `max_rel`.
#'
#' Use the layers you will model with (or their main gradients, e.g.
#' temperature and precipitation), at the resolution you will model at.
#' Typical use: run [filter_occurrences()] with a generous
#' `max_coord_unc_m` (or none), then keep the records with `unc_ok`.
#'
#' **Decide or propagate.** Instead of dropping the records whose
#' environment varies, they can be kept and their uncertainty propagated
#' into the models with [jitter_uncertainty()], which samples the same
#' quantity (the range of the predictors within the radius, relative to the
#' niche breadth) by moving the records in replicate datasets.
#' `unc_action` says which: `"keep"` (the environment is uniform within the
#' radius; the record can be used as it is), `"propagate"` (it varies;
#' drop the record, or keep it and propagate) and `"drop"` (radius above
#' `max_radius_km`, not examined).
#' Records without an uncertainty get `NA` (decide yourself, or set
#' `na_uncertainty` to an assumed value in metres).
#'
#' @param points Records: a data.frame, CSV path, `sf` or
#'   [terra::SpatVector].
#' @param env Predictor layers ([terra::SpatRaster] or file paths).
#' @param uncertainty The coordinate uncertainty in metres: a column name
#'   (default: `coordinateUncertaintyInMeters` or `uncertainty`, whichever
#'   exists), one value per record, or a single value for all.
#' @param max_rel Largest acceptable range within the uncertainty, as a
#'   share of the species' niche breadth, on any layer.
#' @param breadth_probs Quantiles of the values at the records that define
#'   the niche breadth.
#' @param max_radius_km Records with a larger uncertainty are not examined
#'   and get `unc_ok = FALSE` (their circles cover very many cells).
#' @param na_uncertainty Uncertainty (m) assumed for records without one;
#'   `NA` leaves them undecided.
#' @param coords Names of the coordinate columns, if not recognised.
#' @param points_crs CRS of the coordinates.
#' @param verbose Print a summary.
#' @return The records (data.frame) with added columns: `unc_m` (the
#'   uncertainty used), `unc_cells` (cells touched by the circle),
#'   `unc_rel_<layer>` (range within the circle / niche breadth, one column
#'   per layer), `unc_rel_max` (the largest of them), `unc_limiting` (its
#'   layer), `unc_ok` and `unc_action` (`"keep"`, `"propagate"` or
#'   `"drop"`). `attr(, "breadth")` holds the niche breadth of
#'   each layer.
#' @references Smith, A.B., Murphy, S.J., Henderson, D., Erickson, K.D.
#'   (2023) Including imprecisely georeferenced specimens improves accuracy
#'   of species distribution models and estimates of niche breadth. *Global
#'   Ecology and Biogeography* 32, 342-355. Marcer, A., Chapman, A.D.,
#'   Wieczorek, J.R., Pico, F.X., Uribe, F., Waller, J., Arino, A.H. (2022)
#'   Uncertainty matters: ascertaining where specimens in natural history
#'   collections come from and its implications for predicting species
#'   distributions. *Ecography* 2022, e06025.
#' @seealso [jitter_uncertainty()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 60)
#' occ <- sim$presence
#' set.seed(2)
#' occ$coordinateUncertaintyInMeters <- sample(c(100, 5000, 30000, NA), 60,
#'                                             replace = TRUE)
#' u <- uncertainty_environment(occ, sim$env, max_radius_km = 40)
#' table(uncertainty_m = u$unc_m, usable = u$unc_ok, useNA = "ifany")
uncertainty_environment <- function(points, env, uncertainty = NULL,
                                    max_rel = 0.1,
                                    breadth_probs = c(0.05, 0.95),
                                    max_radius_km = 20, na_uncertainty = NA,
                                    coords = NULL, points_crs = "EPSG:4326",
                                    verbose = TRUE) {
  say <- function(...) if (verbose) message(...)
  r <- as_grid(env)
  vars <- names(r)
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  # every vector below has one element per row of df
  if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
    x_in <- as_number(df[[xy_cols[1]]])
    y_in <- as_number(df[[xy_cols[2]]])
  } else {
    p <- as_points(points, coords, points_crs)
    df <- as.data.frame(p)
    points_crs <- terra::crs(p)
    x_in <- terra::crds(p)[, 1]
    y_in <- terra::crds(p)[, 2]
  }
  n <- nrow(df)
  has_xy <- !is.na(x_in) & !is.na(y_in)

  # the uncertainty of each record, in metres
  if (is.null(uncertainty)) {
    uncertainty <- find_column(names(df), c("coordinateUncertaintyInMeters",
                                            "uncertainty", "unc_m"))
    if (is.na(uncertainty)) {
      stop("No uncertainty column found: name it in `uncertainty`.",
           call. = FALSE)
    }
  }
  u <- if (is.character(uncertainty)) {
    if (length(uncertainty) != 1 || !uncertainty %in% names(df)) {
      stop("Column not found: ", paste(uncertainty, collapse = ", "),
           call. = FALSE)
    }
    as_number(df[[uncertainty]])
  } else if (length(uncertainty) %in% c(1, n)) {
    rep_len(as.numeric(uncertainty), n)
  } else {
    stop("`uncertainty` must be a column name, one value per record or a ",
         "single value.", call. = FALSE)
  }
  u[is.na(u)] <- na_uncertainty
  if (any(u < 0, na.rm = TRUE)) stop("Negative uncertainty.", call. = FALSE)

  # position in the raster's CRS, own cell and its values
  px <- py <- lat <- cell <- rep(NA_real_, n)
  if (any(has_xy)) {
    p <- terra::vect(cbind(x_in[has_xy], y_in[has_xy]), crs = points_crs)
    xy <- terra::crds(if (terra::same.crs(p, r)) p else terra::project(p, r))
    px[has_xy] <- xy[, 1]
    py[has_xy] <- xy[, 2]
    lat[has_xy] <- terra::crds(terra::project(p, "EPSG:4326"))[, 2]
    cell[has_xy] <- terra::cellFromXY(r, xy)
  }
  E0 <- env_at_cells(r, cell)
  breadth <- apply(E0, 2, function(v) {
    if (all(is.na(v))) NA_real_ else
      diff(stats::quantile(v, breadth_probs, na.rm = TRUE, names = FALSE))
  })

  # size of a cell in km at each record
  rs <- terra::res(r)
  if (terra::is.lonlat(r)) {
    ns <- rs[2] * 111.32
    ew <- rs[1] * 111.32 * pmax(cos(lat * pi / 180), 0.01)
  } else {
    ns <- rs[2] / 1000
    ew <- rep(rs[1] / 1000, n)
  }
  rad <- u / 1000
  too_far <- !is.na(rad) & rad > max_radius_km
  scan <- which(!is.na(cell) & !is.na(rad) & !too_far)
  rmin <- rmax <- matrix(NA_real_, n, length(vars),
                         dimnames = list(NULL, vars))
  n_cells <- rep(NA_integer_, n)
  if (length(scan)) {
    rc <- terra::rowColFromCell(r, cell[scan])
    ctr <- terra::xyFromCell(r, cell[scan])
    # the point's offset from its cell centre, in cells
    ox <- (px[scan] - ctr[, 1]) / rs[1]
    oy <- (py[scan] - ctr[, 2]) / rs[2]
    Dr <- ceiling(max(rad[scan]) / ns) + 1
    Dc <- ceiling(max(rad[scan] / ew[scan])) + 1
    if ((2 * Dr + 1) * (2 * Dc + 1) > 4e5) {
      stop("The largest uncertainty covers too many cells (window of ",
           2 * Dr + 1, " x ", 2 * Dc + 1, " cells): lower `max_radius_km`.",
           call. = FALSE)
    }
    lo <- matrix(Inf, length(scan), length(vars))
    hi <- matrix(-Inf, length(scan), length(vars))
    cnt <- integer(length(scan))
    for (dr in -Dr:Dr) for (dc in -Dc:Dc) {
      # distance (km) from the point to the rectangle of the neighbouring
      # cell: the circle touches the cell when it is within the radius
      gx <- pmax(0, abs(dc - ox) - 0.5) * ew[scan]
      gy <- pmax(0, abs(-dr - oy) - 0.5) * ns
      inc <- which(sqrt(gx^2 + gy^2) <= rad[scan] | (dr == 0 & dc == 0))
      if (!length(inc)) next
      X <- env_at_cells(r, terra::cellFromRowCol(r, rc[inc, 1] + dr,
                                                 rc[inc, 2] + dc))
      lo[inc, ] <- pmin(lo[inc, , drop = FALSE], X, na.rm = TRUE)
      hi[inc, ] <- pmax(hi[inc, , drop = FALSE], X, na.rm = TRUE)
      cnt[inc] <- cnt[inc] + stats::complete.cases(X)
    }
    lo[!is.finite(lo)] <- NA
    hi[!is.finite(hi)] <- NA
    rmin[scan, ] <- lo
    rmax[scan, ] <- hi
    n_cells[scan] <- cnt
  }
  rel <- sweep(rmax - rmin, 2, breadth, "/")
  rel[is.nan(rel)] <- 0                       # no range and no breadth
  colnames(rel) <- paste0("unc_rel_", vars)
  known <- rowSums(!is.na(rel)) > 0
  worst <- rep(NA_real_, n)
  limiting <- rep(NA_character_, n)
  if (any(known)) {
    z <- rel[known, , drop = FALSE]
    z[is.na(z)] <- -Inf
    worst[known] <- apply(z, 1, max)
    limiting[known] <- vars[max.col(z, ties.method = "first")]
  }
  ok <- worst <= max_rel
  ok[too_far] <- FALSE

  out <- df
  out$unc_m <- u
  out$unc_cells <- n_cells
  out <- cbind(out, as.data.frame(rel))
  out$unc_rel_max <- worst
  out$unc_limiting <- limiting
  out$unc_ok <- ok
  # keep (uniform environment), propagate (varies, but within reach of
  # jitter_uncertainty()) or drop (radius above max_radius_km)
  out$unc_action <- ifelse(too_far, "drop", ifelse(is.na(ok), NA_character_,
                                                   ifelse(ok, "keep",
                                                          "propagate")))
  attr(out, "breadth") <- breadth
  say("Records: ", n, "; usable (range within the uncertainty <= ",
      max_rel, " of the niche breadth on every layer): ",
      sum(ok, na.rm = TRUE), "; too uncertain: ", sum(!ok, na.rm = TRUE),
      if (any(too_far)) paste0(" (", sum(too_far), " above ", max_radius_km,
                               " km, not examined)"),
      "; undecided (no uncertainty, coordinates or values): ",
      sum(is.na(ok)))
  if (any(!ok & !too_far, na.rm = TRUE)) {
    t <- sort(table(limiting[!is.na(ok) & !ok & !too_far]),
              decreasing = TRUE)
    say("Limiting layer of the records that are too uncertain: ",
        paste0(names(t), " (", as.integer(t), ")", collapse = ", "))
  }
  out
}

# Propagating coordinate uncertainty --------------------------------------------

#' Replicate datasets with each record moved within its coordinate uncertainty
#'
#' [uncertainty_environment()] decides whether an imprecise record is usable
#' (keep or drop). The stronger alternative is to *propagate* the
#' uncertainty: fit the model once per replicate dataset in which every
#' record is placed at a random position within its own uncertainty radius,
#' and look at how much the predictions vary. Models built from positionally
#' uncertain records can remain ecologically interpretable: in simulations,
#' positional error lowered predictive performance more than it changed
#' variable importance and response shapes, especially for spatially
#' homogeneous predictors (Gabor et al. 2023). Placement within the radius
#' can also be optimised against model performance (Wu et al. 2026); that
#' needs the model and is left to the modelling software.
#'
#' **Placement.** Each record is moved to a point drawn uniformly over the
#' area of its uncertainty circle (distance `R * sqrt(u)`, random
#' direction). With `m`, `land` and `require_values`, points outside M, off
#' land or in cells without values in every layer of `env` are drawn again
#' (up to `max_tries` times, as Gabor et al. 2023 re-drew points shifted
#' outside the study area); a record with no valid position stays where it
#' was (`jitter_failed`). Records with an uncertainty of 0 do not move;
#' records without one move by `na_uncertainty` metres, or stay fixed when
#' it is `NA`.
#'
#' **Is it worth it?** For each record, the range of each layer over the
#' replicates is divided by the species' niche breadth on that layer (the
#' distance between the `breadth_probs` quantiles of the values at the
#' records), the same scale as [uncertainty_environment()]. Records whose
#' environment changes by more than `max_rel` of the breadth are the ones
#' that make propagating worth the extra model runs; when there are almost
#' none, the replicates are environmentally the same dataset.
#'
#' **Link with [uncertainty_environment()].** Both measure how much the
#' predictors vary within a record's uncertainty, one analytically (every
#' cell the circle touches), one by sampling. `propagate = "variable"` uses
#' the analytical check to move only the records whose environment varies
#' (`unc_ok = FALSE`), keeping environmentally uniform records fixed, so
#' that the replicates differ only where it matters; with
#' [uncertainty_environment()] the workflow becomes keep (uniform), propagate
#' (variable) or drop (radius above `max_radius_km`).
#'
#' @param points Records: a data.frame, CSV path, `sf` or
#'   [terra::SpatVector], or the result of [uncertainty_environment()] (its
#'   `unc_m` is then used).
#' @param env Optional predictor layers ([terra::SpatRaster] or files), for
#'   the `require_values` constraint and the summary.
#' @param n Number of replicate datasets.
#' @param uncertainty Coordinate uncertainty in metres: a column name
#'   (default: `unc_m`, `coordinateUncertaintyInMeters` or `uncertainty`),
#'   one value per record, or a single value.
#' @param na_uncertainty Radius (m) for records without an uncertainty;
#'   `NA` (default) keeps them fixed.
#' @param max_radius_km Records with a larger radius are not moved (they are
#'   reported; drop them or decide yourself).
#' @param m Optional accessible area: polygons (`sf`, SpatVector, path) or
#'   an `m_area` from [delimit_m()]; jittered points must fall inside.
#' @param land Optional land polygons or a raster (cells other than `NA`
#'   and 0 are land); jittered points must fall on land.
#' @param require_values Jittered points must fall in cells with values in
#'   every layer of `env`.
#' @param propagate `"all"` (default) moves every record with an
#'   uncertainty; `"variable"` only those whose environment varies within
#'   the radius by [uncertainty_environment()] (needs `env`).
#' @param max_rel,breadth_probs As in [uncertainty_environment()].
#' @param max_tries Draws per record and replicate before giving up.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Optional random seed.
#' @param verbose Print progress.
#' @return An object of class `jittered_records`: `replicates` (a list of
#'   `n` data.frames, each the records with jittered coordinates in the
#'   original coordinate columns, `lon_orig`, `lat_orig`, `record` (row of
#'   the input), `replicate`, `jitter_m`, `moved` and `jitter_failed`),
#'   `records` (one row per record: `unc_m`, `moved`, failures, and with
#'   `env` the relative range of each layer over the replicates,
#'   `rel_range_max`, its `limiting` layer and `varies`), `variables` (per
#'   layer: median and 90th percentile of the relative range, share of
#'   records above `max_rel`), `breadth` and `settings`. Loop over
#'   `x$replicates` to fit one model per replicate, then summarise the
#'   predictions with [coordinate_spread()].
#' @references Gabor, L., Jetz, W., Zarzo-Arias, A., Winner, K., Yanco, S.,
#'   Pinkert, S., Marsh, C.J., Rogan, M.S., Makinen, J., Rocchini, D.,
#'   Bartak, V., Malavasi, M., Balej, P., Moudry, V. (2023) Species
#'   distribution models affected by positional uncertainty in species
#'   occurrences can still be ecologically interpretable. *Ecography* 2023,
#'   e06358. \doi{10.1111/ecog.06358}
#'
#'   Wu, H., Liu, Y., Wu, Z., Abulike, K., Wang, J. (2026) A simulated
#'   annealing optimization framework to reduce positional uncertainty and
#'   improve performance of species distribution models. *Environmental
#'   Modelling & Software* 204, 107096. \doi{10.1016/j.envsoft.2026.107096}
#' @seealso [uncertainty_environment()], [coordinate_spread()],
#'   [plot_jitter()], [write_jitter()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 60)
#' occ <- sim$presence
#' occ$coordinateUncertaintyInMeters <- rep(c(0, 2000, 20000), 20)
#' j <- jitter_uncertainty(occ, sim$env, n = 5, seed = 1)
#' j
#' head(j$replicates[[1]])
#' plot_jitter(j)
#' \donttest{
#' # fit one model per replicate in your software, e.g.
#' # for (d in j$replicates) fit_and_predict(d)
#' }
jitter_uncertainty <- function(points, env = NULL, n = 10, uncertainty = NULL,
                               na_uncertainty = NA, max_radius_km = 50,
                               m = NULL, land = NULL, require_values = TRUE,
                               propagate = c("all", "variable"),
                               max_rel = 0.1, breadth_probs = c(0.05, 0.95),
                               max_tries = 100, coords = NULL,
                               points_crs = "EPSG:4326", seed = NULL,
                               verbose = TRUE) {
  propagate <- match.arg(propagate)
  say <- function(...) if (verbose) message(...)
  if (n < 1) stop("`n` must be at least 1.", call. = FALSE)
  r <- if (is.null(env)) NULL else as_grid(env)
  if (propagate == "variable" && is.null(r)) {
    stop("propagate = \"variable\" needs `env`.", call. = FALSE)
  }
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
  } else {
    p <- as_points(points, coords, points_crs)
    points_crs <- terra::crs(p)
    df <- as.data.frame(p)
    xy_cols <- c("lon", "lat")
    df$lon <- terra::crds(p)[, 1]
    df$lat <- terra::crds(p)[, 2]
  }
  x0 <- as_number(df[[xy_cols[1]]])
  y0 <- as_number(df[[xy_cols[2]]])
  N <- nrow(df)
  if (any(is.na(x0) | is.na(y0))) {
    stop("Records without coordinates: remove them first.", call. = FALSE)
  }
  geo <- terra::is.lonlat(terra::vect(cbind(0, 0), crs = points_crs))

  # uncertainty (m) of each record
  if (is.null(uncertainty)) {
    uncertainty <- find_column(names(df), c("unc_m",
                                            "coordinateUncertaintyInMeters",
                                            "uncertainty"))
    if (is.na(uncertainty)) {
      stop("No uncertainty column found: name it in `uncertainty`.",
           call. = FALSE)
    }
  }
  u <- if (is.character(uncertainty)) {
    if (length(uncertainty) != 1 || !uncertainty %in% names(df)) {
      stop("Column not found: ", paste(uncertainty, collapse = ", "),
           call. = FALSE)
    }
    as_number(df[[uncertainty]])
  } else if (length(uncertainty) %in% c(1, N)) {
    rep_len(as.numeric(uncertainty), N)
  } else {
    stop("`uncertainty` must be a column name, one value per record or a ",
         "single value.", call. = FALSE)
  }
  no_unc <- is.na(u)
  u[no_unc] <- na_uncertainty
  if (any(u < 0, na.rm = TRUE)) stop("Negative uncertainty.", call. = FALSE)
  too_far <- !is.na(u) & u / 1000 > max_radius_km
  move <- !is.na(u) & u > 0 & !too_far
  status <- ifelse(is.na(u), "fixed: no uncertainty",
                   ifelse(too_far, "fixed: radius above max_radius_km",
                          ifelse(u == 0, "fixed: uncertainty 0", "moved")))

  # only records whose environment varies (analytical check)
  unc_ok <- rep(NA, N)
  if (propagate == "variable") {
    ue <- uncertainty_environment(
      data.frame(.x = x0, .y = y0), r, uncertainty = ifelse(is.na(u), 0, u),
      max_rel = max_rel, breadth_probs = breadth_probs,
      max_radius_km = max_radius_km, coords = c(".x", ".y"),
      points_crs = points_crs, verbose = FALSE)
    unc_ok <- ue$unc_ok
    uniform <- move & unc_ok %in% TRUE
    move[uniform] <- FALSE
    status[uniform] <- "fixed: environment uniform within the radius"
  }

  # constraints
  m_sfc <- NULL
  if (!is.null(m)) {
    if (is.character(m) && length(m) == 1 && grepl("\\.gpkg$", m) &&
        "M" %in% sf::st_layers(m)$name) {
      m <- sf::st_read(m, layer = "M", quiet = TRUE)
    }
    m_sfc <- sf::st_union(if (inherits(m, "m_area")) m$m else as_sfc(m, "m"))
  }
  land_r <- NULL
  land_sfc <- NULL
  if (!is.null(land)) {
    if (inherits(land, "SpatRaster") ||
        (is.character(land) && grepl("\\.(tif|tiff|nc|grd|asc)$", land[1],
                                     ignore.case = TRUE))) {
      land_r <- as_grid(land)[[1]]
    } else {
      land_sfc <- sf::st_union(as_sfc(land, "land"))
    }
  }
  inside <- function(sfc, X, Y) {
    pts <- sf::st_as_sf(data.frame(x = X, y = Y), coords = c("x", "y"),
                        crs = points_crs)
    pts <- sf::st_transform(pts, sf::st_crs(sfc))
    lengths(sf::st_intersects(pts, sfc)) > 0
  }
  valid <- function(X, Y) {
    ok <- rep(TRUE, length(X))
    if (!is.null(m_sfc)) ok <- ok & inside(m_sfc, X, Y)
    if (!is.null(land_sfc) && any(ok)) {
      ok[ok] <- inside(land_sfc, X[ok], Y[ok])
    }
    if (!is.null(land_r) && any(ok)) {
      v <- values_xy(land_r, cbind(X[ok], Y[ok]), points_crs)[, 1]
      ok[ok] <- !is.na(v) & v != 0
    }
    if (require_values && !is.null(r) && any(ok)) {
      v <- values_xy(r, cbind(X[ok], Y[ok]), points_crs)
      ok[ok] <- stats::complete.cases(v)
    }
    ok
  }
  offset <- function(i, k) {
    rad <- u[i] * sqrt(stats::runif(k))
    ang <- stats::runif(k, 0, 2 * pi)
    dx <- rad * cos(ang)
    dy <- rad * sin(ang)
    if (geo) {
      cbind(x0[i] + dx / (111320 * pmax(cos(y0[i] * pi / 180), 0.01)),
            y0[i] + dy / 111320)
    } else cbind(x0[i] + dx, y0[i] + dy)
  }

  # draw all replicates at once; re-draw the invalid ones
  say("Jittering ", sum(move), " of ", N, " records in ", n,
      " replicate(s)")
  idx <- rep(seq_len(N), times = n)
  X <- rep(x0, n)
  Y <- rep(y0, n)
  failed <- rep(FALSE, N * n)
  todo <- which(move[idx])
  with_seed(seed, {
    tries <- 0
    while (length(todo) && tries < max_tries) {
      tries <- tries + 1
      xy <- offset(idx[todo], length(todo))
      ok <- valid(xy[, 1], xy[, 2])
      X[todo[ok]] <- xy[ok, 1]
      Y[todo[ok]] <- xy[ok, 2]
      todo <- todo[!ok]
    }
  })
  failed[todo] <- TRUE
  dist_m <- if (geo) {
    sqrt(((X - rep(x0, n)) * 111320 * cos(rep(y0, n) * pi / 180))^2 +
           ((Y - rep(y0, n)) * 111320)^2)
  } else sqrt((X - rep(x0, n))^2 + (Y - rep(y0, n))^2)

  reps <- lapply(seq_len(n), function(k) {
    j <- (k - 1) * N + seq_len(N)
    d <- df
    d[[xy_cols[1]]] <- X[j]
    d[[xy_cols[2]]] <- Y[j]
    d$lon_orig <- x0
    d$lat_orig <- y0
    d$record <- seq_len(N)
    d$replicate <- k
    d$jitter_m <- dist_m[j]
    d$moved <- move & !failed[j]
    d$jitter_failed <- failed[j]
    d
  })
  rec <- data.frame(record = seq_len(N), unc_m = u, status = status,
                    moved = move,
                    failed_replicates = as.integer(tabulate(idx[failed],
                                                            nbins = N)),
                    stringsAsFactors = FALSE)
  if (propagate == "variable") rec$unc_ok <- unc_ok

  # how much the environment at each record varies over the replicates
  vars_tab <- NULL
  breadth <- NULL
  if (!is.null(r)) {
    E0 <- values_xy(r, cbind(x0, y0), points_crs)
    breadth <- apply(E0, 2, function(v) {
      if (all(is.na(v))) NA_real_ else
        diff(stats::quantile(v, breadth_probs, na.rm = TRUE, names = FALSE))
    })
    E <- values_xy(r, cbind(X, Y), points_crs)
    vars <- colnames(E)
    rel <- vapply(vars, function(v) {
      val <- matrix(E[, v], nrow = N)          # records x replicates
      val <- cbind(E0[, v], val)               # the original position too
      rg <- apply(val, 1, function(z) {
        if (all(is.na(z))) NA_real_ else diff(range(z, na.rm = TRUE))
      })
      rg / breadth[[v]]
    }, numeric(N))
    rel <- matrix(rel, nrow = N, dimnames = list(NULL, vars))
    rel[is.nan(rel)] <- 0
    for (v in vars) rec[[paste0("rel_range_", v)]] <- rel[, v]
    z <- rel
    z[is.na(z)] <- -Inf
    rec$rel_range_max <- apply(z, 1, max)
    rec$rel_range_max[!is.finite(rec$rel_range_max)] <- NA
    rec$limiting <- ifelse(is.na(rec$rel_range_max), NA_character_,
                           vars[max.col(z, ties.method = "first")])
    rec$varies <- !is.na(rec$rel_range_max) & rec$rel_range_max > max_rel
    mv <- rec$moved
    vars_tab <- data.frame(
      layer = vars, breadth = unname(breadth[vars]),
      median_rel_range = vapply(vars, function(v)
        stats::median(rel[mv, v], na.rm = TRUE), numeric(1)),
      q90_rel_range = vapply(vars, function(v)
        stats::quantile(rel[mv, v], 0.9, na.rm = TRUE, names = FALSE),
        numeric(1)),
      share_above = vapply(vars, function(v)
        mean(rel[mv, v] > max_rel, na.rm = TRUE), numeric(1)),
      row.names = NULL, stringsAsFactors = FALSE)
    if (!any(mv)) vars_tab[, -(1:2)] <- NA_real_
  }
  out <- list(replicates = reps, records = rec, variables = vars_tab,
              breadth = breadth,
              settings = list(n = n, coords = xy_cols, crs = points_crs,
                              max_rel = max_rel, max_radius_km = max_radius_km,
                              na_uncertainty = na_uncertainty,
                              propagate = propagate,
                              constraints = c(m = !is.null(m),
                                              land = !is.null(land),
                                              values = require_values &&
                                                !is.null(r))))
  class(out) <- "jittered_records"
  out
}

#' @export
print.jittered_records <- function(x, digits = 2, ...) {
  s <- x$settings
  rec <- x$records
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  N <- nrow(rec)
  nm <- sum(rec$moved)
  cat("Jittered records: ", N, " records x ", s$n, " replicate(s); ",
      nm, " moved within their uncertainty\n\nIN SHORT\n", sep = "")
  if (nm == 0) {
    say("No record was moved (no uncertainty above 0 within max_radius_km): the replicates are identical and there is nothing to propagate.")
  } else if (!is.null(x$variables)) {
    nv <- sum(rec$varies[rec$moved])
    sh <- nv / nm
    say(paste0(nv, " of the ", nm, " moved records (", share_label(sh),
               ") change their environment by more than ", s$max_rel,
               " of the niche breadth on at least one layer across the ",
               "replicates. ",
               if (sh < 0.05) {
                 paste0("The replicates are environmentally almost the same ",
                        "dataset: propagating is unlikely to be worth the ",
                        "extra model runs.")
               } else if (sh < 0.25) {
                 paste0("Propagating matters for a minority of records; ",
                        "it is worth a few replicates (e.g. 10) to see ",
                        "whether the predictions move.")
               } else {
                 paste0("Propagating is worthwhile: fit one model per ",
                        "replicate and summarise the predictions with ",
                        "coordinate_spread().")
               }))
    w <- which.max(x$variables$share_above)
    if (length(w) && x$variables$share_above[w] > 0) {
      say(paste0("The most sensitive layer is ", x$variables$layer[w],
                 " (", share_label(x$variables$share_above[w]),
                 " of the moved records above the limit)."))
    }
    if (s$n < 10) {
      say(paste0("Only ", s$n, " replicates: the ranges are underestimated;",
                 " use 10 or more for the models."))
    }
  } else {
    say(paste0(nm, " records were moved. Give `env` to see how much their ",
               "environment changes, the quantity that tells whether ",
               "propagating is worth the extra model runs."))
  }
  fx <- table(rec$status[!rec$moved])
  if (length(fx)) {
    say(paste0("Not moved: ", paste0(names(fx), " (", as.integer(fx), ")",
                                     collapse = "; "), "."))
  }
  fl <- sum(rec$failed_replicates > 0)
  if (fl) {
    say(paste0(fl, " record(s) found no valid position within their radius ",
               "in some replicates (constraints: ",
               paste(names(s$constraints)[s$constraints], collapse = ", "),
               ") and stayed at their original position there."))
  }
  if (!is.null(x$variables)) {
    cat("\n")
    v <- x$variables
    num <- vapply(v, is.double, logical(1))
    v[num] <- lapply(v[num], signif, 3)
    print(v, row.names = FALSE)
  }
  cat("", strwrap(paste0(
    "Points are drawn uniformly over the area of each record's uncertainty ",
    "circle", if (any(s$constraints)) paste0(", re-drawn until they ",
                                              "satisfy the constraints (",
                                              paste(names(s$constraints)[
                                                s$constraints],
                                                collapse = ", "), ")"),
    ". rel_range = range of a layer at a record over the replicates (and ",
    "its original position) / the niche breadth on that layer (5-95% ",
    "quantiles at the records), as in uncertainty_environment(); ",
    "share_above = share of moved records above ", s$max_rel, ". The ",
    "limits of 5% and 25% of records in the advice are rules of thumb, not ",
    "values from the literature. Gabor et al. (2023) found that positional ",
    "error lowered model performance more than it changed variable ",
    "importance and response shapes."), width = 76), sep = "\n")
  invisible(x)
}

#' Plot how much the environment of each record varies when jittered
#'
#' One point per record: its coordinate uncertainty against the largest
#' relative range of the predictors over the replicates, coloured by
#' whether it exceeds `max_rel`.
#'
#' @param x A `jittered_records` from [jitter_uncertainty()] (with `env`).
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Gabor, L. et al. (2023) Species distribution models affected
#'   by positional uncertainty in species occurrences can still be
#'   ecologically interpretable. *Ecography* 2023, e06358.
#'   \doi{10.1111/ecog.06358}
#' @export
#' @examples
#' # see ?jitter_uncertainty
plot_jitter <- function(x, base_size = 11, file = NULL, width = 7,
                        height = 4.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "jittered_records")) {
    stop("`x` must come from jitter_uncertainty().", call. = FALSE)
  }
  rec <- x$records
  if (is.null(rec$rel_range_max)) {
    stop("No environmental summary: run jitter_uncertainty() with `env`.",
         call. = FALSE)
  }
  d <- rec[rec$moved & !is.na(rec$rel_range_max), , drop = FALSE]
  if (!nrow(d)) stop("No moved record to plot.", call. = FALSE)
  lim <- x$settings$max_rel
  lv <- c("above", "below")
  cnt <- c(sum(d$varies), sum(!d$varies))
  lab <- paste0(c(paste0("Varies more than ", lim, " of the breadth"),
                  paste0("Varies up to ", lim, " of the breadth")), ": ",
                cnt, " (", share_label(cnt / nrow(d)), ")")
  d$cls <- factor(ifelse(d$varies, lv[1], lv[2]), lv, lab)
  ggplot2::ggplot(d, ggplot2::aes(.data$unc_m / 1000, .data$rel_range_max,
                                  colour = .data$cls)) +
    ggplot2::geom_hline(yintercept = lim, linetype = "dashed",
                        colour = "#8a8a86") +
    ggplot2::geom_point(size = 1.8, alpha = 0.8) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c("#d95f02", "#1b9e77"), lab), name = NULL, drop = FALSE) +
    ggplot2::labs(x = "Coordinate uncertainty (km)",
                  y = "Largest relative range of the predictors",
                  title = "Environment of the records across replicates",
                  subtitle = paste0(nrow(d), " moved records, ",
                                    x$settings$n, " replicates; dashed: ",
                                    lim, " of the niche breadth")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom", legend.direction = "vertical")
}

#' Write jittered replicate datasets to CSV files
#'
#' @param x A `jittered_records` from [jitter_uncertainty()].
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @return Invisibly, the paths: one file per replicate
#'   (`<prefix>_rep01.csv`, ...) and `<prefix>_records.csv` with the
#'   per-record summary.
#' @references Gabor, L. et al. (2023) Species distribution models affected
#'   by positional uncertainty in species occurrences can still be
#'   ecologically interpretable. *Ecography* 2023, e06358.
#'   \doi{10.1111/ecog.06358}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 20)
#' j <- jitter_uncertainty(sim$presence, n = 2, uncertainty = 1000,
#'                         seed = 1, verbose = FALSE)
#' write_jitter(j, tempdir())
write_jitter <- function(x, dir, prefix = "jittered") {
  if (!inherits(x, "jittered_records")) {
    stop("`x` must come from jitter_uncertainty().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  w <- nchar(length(x$replicates))
  f <- file.path(dir, sprintf(paste0("%s_rep%0", max(2, w), "d.csv"),
                              prefix, seq_along(x$replicates)))
  for (k in seq_along(f)) {
    utils::write.csv(x$replicates[[k]], f[k], row.names = FALSE)
  }
  fr <- file.path(dir, paste0(prefix, "_records.csv"))
  utils::write.csv(x$records, fr, row.names = FALSE)
  invisible(c(f, fr))
}

#' Spread of predictions attributable to coordinate error
#'
#' Closes the loop of [jitter_uncertainty()]: given one prediction raster
#' per replicate dataset (the same model, settings and predictors, fitted
#' to each replicate), returns the per-cell mean, standard deviation and
#' coefficient of variation. Because only the positions of the records
#' differ between the replicates, the spread is the part of the prediction
#' uncertainty due to coordinate error alone. It is computed in the same
#' way as the disagreement layer of [ensemble_models()] (the replicates
#' play the role of the models, combined by their mean), so the two can be
#' compared cell by cell: with `ensemble`, the summary says where
#' coordinate error moves the predictions more than the choice of model.
#'
#' @param predictions One prediction raster per replicate: a
#'   [terra::SpatRaster], raster files or a folder.
#' @param ensemble Optional result of [ensemble_models()] on the same grid,
#'   for the comparison with the disagreement between models.
#' @param levels Levels of `sd` for the area summary (for predictions on
#'   the 0-1 scale).
#' @param n_sample Cells sampled for the summaries.
#' @param rescale Bring all predictions to the 0-1 scale first (see
#'   [ensemble_models()]).
#' @param filename Optional file for the raster (layers `mean`, `sd`,
#'   `cv`).
#' @param overwrite Overwrite `filename`.
#' @param pattern,recursive File search in folders.
#' @param seed Seed of the cell sample.
#' @return An object of class `coordinate_spread`: `raster` (layers
#'   `mean`, `sd`, `cv`), `summary`, `area`, `comparison` (with
#'   `ensemble`), `n_replicates` and `n_sample`. See
#'   [plot_coordinate_spread()].
#' @references Gabor, L. et al. (2023) Species distribution models affected
#'   by positional uncertainty in species occurrences can still be
#'   ecologically interpretable. *Ecography* 2023, e06358.
#'   \doi{10.1111/ecog.06358}
#' @seealso [jitter_uncertainty()], [ensemble_models()],
#'   [partition_uncertainty()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' # stand-ins for predictions from models fitted to 5 jittered replicates
#' reps <- c(sim$models[[1]], sim$models[[1]] * 0.95,
#'           sim$models[[1]] * 0.9 + 0.05, sim$models[[1]] * 0.98,
#'           sim$models[[1]] * 0.92 + 0.02)
#' cs <- coordinate_spread(reps, ensemble = ensemble_models(sim$models))
#' cs
#' plot_coordinate_spread(cs)
coordinate_spread <- function(predictions, ensemble = NULL,
                              levels = c(0.05, 0.1), n_sample = 1e5,
                              rescale = TRUE, filename = NULL,
                              overwrite = FALSE, pattern = "\\.tif{1,2}$",
                              recursive = FALSE, seed = 1) {
  x <- ensemble_stack(predictions, pattern, recursive, unique_names = FALSE)
  if (anyDuplicated(names(x))) {
    names(x) <- sprintf("replicate_%02d", seq_len(terra::nlyr(x)))
  }
  if (terra::nlyr(x) < 2) {
    stop("At least two replicate predictions are needed.", call. = FALSE)
  }
  e <- suppressWarnings(suppressMessages(ensemble_models(
    x, method = "mean", n_sample = n_sample, levels = levels, seed = seed,
    rescale = rescale)))
  out <- e$raster
  names(out) <- c("mean", "sd", "cv")
  if (!is.null(filename)) {
    out <- terra::writeRaster(out, filename, overwrite = overwrite,
                              gdal = c("COMPRESS=LZW", "TILED=YES",
                                       "BIGTIFF=IF_SAFER"))
  }
  comp <- NULL
  if (!is.null(ensemble)) {
    if (!inherits(ensemble, "enm_ensemble")) {
      stop("`ensemble` must come from ensemble_models().", call. = FALSE)
    }
    both <- c(out[["sd"]], ensemble$raster[["sd"]])
    names(both) <- c("coordinates", "models")
    smp <- with_seed(seed, sample_complete(both, n_sample))
    comp <- data.frame(
      median_sd_coordinates = stats::median(smp[, 1]),
      median_sd_models = stats::median(smp[, 2]),
      share_coordinates_larger = mean(smp[, 1] > smp[, 2]),
      median_ratio = stats::median(ifelse(smp[, 2] > 0, smp[, 1] / smp[, 2],
                                          NA), na.rm = TRUE))
  }
  structure(list(raster = out, summary = e$summary, area = e$area,
                 comparison = comp, scales = e$scales,
                 n_replicates = terra::nlyr(x), n_sample = e$n_sample,
                 scale01 = e$scale01),
            class = "coordinate_spread")
}

#' @export
print.coordinate_spread <- function(x, digits = 3, ...) {
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(v, digits = digits, format = "f")
  s <- x$summary
  msd <- s$median[s$statistic == "sd"]
  cat("Prediction spread from coordinate error: ", x$n_replicates,
      " replicate predictions\n\nIN SHORT\n", sep = "")
  say(paste0("Moving the records within their coordinate uncertainty ",
             "changes the predictions by a median standard deviation of ",
             f(msd), " per cell (90% of cells below ",
             f(s$q90[s$statistic == "sd"]), ")",
             if (x$scale01) paste0(
               ": ", if (msd < 0.02) "negligible" else if (msd < 0.05)
                 "small" else if (msd < 0.1) "moderate" else "large",
               " on the 0-1 scale") else "", "."))
  if (!is.null(x$area)) {
    for (i in seq_len(nrow(x$area))) {
      say(paste0("Cells with sd above ", x$area$sd_above[i], ": ",
                 share_label(x$area$share[i]), " of the area (",
                 formatC(round(x$area$area_km2[i]), format = "d",
                         big.mark = " "), " km2)."))
    }
  }
  cmp <- x$comparison
  if (!is.null(cmp)) {
    say(paste0("Compared with the disagreement between models (median sd ",
               f(cmp$median_sd_models), "), coordinate error moves the ",
               "predictions more in ", share_label(cmp$share_coordinates_larger),
               " of the cells (median ratio ", f(cmp$median_ratio), "). ",
               if (cmp$share_coordinates_larger > 0.5)
                 "Coordinate error is the larger source of uncertainty here."
               else "The choice of model is the larger source here."))
  }
  if (x$n_replicates < 10) {
    say(paste0("Only ", x$n_replicates, " replicates: the spread is ",
               "uncertain; 10 or more are better."))
  }
  cat("\n")
  print(s, row.names = FALSE, digits = digits)
  cat("", strwrap(paste0(
    "sd and cv are computed per cell over the replicate predictions as in ",
    "ensemble_models(); the replicates differ only in the positions of the ",
    "records, so the spread is due to coordinate error alone. The grades ",
    "(median sd below 0.02 negligible, 0.05 small, 0.1 moderate) are rules ",
    "of thumb for predictions on the 0-1 scale."), width = 76), sep = "\n")
  invisible(x)
}

#' Plot the prediction spread due to coordinate error
#'
#' Maps the mean of the replicate predictions and their standard deviation
#' (or coefficient of variation), from [coordinate_spread()].
#'
#' @param x A result of [coordinate_spread()].
#' @param type `"both"`, `"mean"` or `"spread"`.
#' @param spread `"sd"` or `"cv"`.
#' @param max_cells Largest number of cells drawn.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to.
#' @param width,height Size of the saved plot in inches (default 13 x 5.5
#'   for both panels, 7 x 5.5 for one).
#' @return A ggplot, or for `type = "both"` invisibly a list of two.
#' @references Gabor, L. et al. (2023) Species distribution models affected
#'   by positional uncertainty in species occurrences can still be
#'   ecologically interpretable. *Ecography* 2023, e06358.
#'   \doi{10.1111/ecog.06358}
#' @export
#' @examples
#' # see ?coordinate_spread
plot_coordinate_spread <- function(x, type = c("both", "mean", "spread"),
                                   spread = c("sd", "cv"), max_cells = 5e5,
                                   base_size = 11, file = NULL, width = NULL,
                                   height = NULL) {
  if (!inherits(x, "coordinate_spread")) {
    stop("`x` must come from coordinate_spread().", call. = FALSE)
  }
  type <- match.arg(type)
  spread <- match.arg(spread)
  e <- structure(list(raster = c(x$raster[["mean"]], x$raster[[spread]]),
                      method = "mean", models = seq_len(x$n_replicates)),
                 class = "enm_ensemble")
  names(e$raster) <- c("ensemble", spread)
  one <- function(tp) {
    plot_ensemble(e, type = tp, disagreement = spread, max_cells = max_cells,
                  base_size = base_size)
  }
  relabel <- function(q, title, sub) {
    q + ggplot2::labs(title = title, subtitle = sub)
  }
  n <- x$n_replicates
  pm <- function(q) relabel(q, "Mean prediction",
                            paste0("Mean of ", n, " replicate predictions."))
  ps <- function(q) relabel(q, "Spread due to coordinate error",
                            paste0(if (spread == "sd") "Standard deviation"
                                   else "Coefficient of variation",
                                   " of the ", n,
                                   " replicate predictions per cell."))
  if (type == "both") {
    p <- list(mean = pm(one("ensemble")), spread = ps(one("disagreement")))
    if (!is.null(file)) {
      save_figure(p, file, if (is.null(width)) 13 else width,
                  if (is.null(height)) 5.5 else height)
    }
    print_side_by_side(p)
    return(invisible(p))
  }
  p <- if (type == "mean") pm(one("ensemble")) else ps(one("disagreement"))
  if (!is.null(file)) {
    save_figure(p, file, if (is.null(width)) 7 else width,
                if (is.null(height)) 5.5 else height)
  }
  p
}

# Values of all layers at coordinates in `crs` (internal): a matrix with
# one column per layer, NA outside the raster
values_xy <- function(r, xy, crs) {
  p <- terra::vect(xy, crs = crs)
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  v <- as.matrix(terra::extract(r, p, ID = FALSE))
  colnames(v) <- names(r)
  v
}
