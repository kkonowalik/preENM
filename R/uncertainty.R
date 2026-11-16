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
#'   layer) and `unc_ok`. `attr(, "breadth")` holds the niche breadth of
#'   each layer.
#' @references Smith, A.B., Murphy, S.J., Henderson, D., Erickson, K.D.
#'   (2023) Including imprecisely georeferenced specimens improves accuracy
#'   of species distribution models and estimates of niche breadth. *Global
#'   Ecology and Biogeography* 32, 342-355. Marcer, A., Chapman, A.D.,
#'   Wieczorek, J.R., Pico, F.X., Uribe, F., Waller, J., Arino, A.H. (2022)
#'   Uncertainty matters: ascertaining where specimens in natural history
#'   collections come from and its implications for predicting species
#'   distributions. *Ecography* 2022, e06025.
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
