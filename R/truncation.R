# Niche truncation report ---------------------------------------------------------

#' Is the niche truncated by the calibration area?
#'
#' A check to run before modelling. A model can only learn the part of a
#' species' response that the calibration area (M) contains. If the records
#' are still frequent where the conditions of M end, the limit of the niche
#' on that variable lies outside M: the response is *truncated*, and
#' projecting the model to conditions beyond M (another period, scenario or
#' region) extrapolates a curve whose end was never seen (Thuiller et al.
#' 2004; Chevalier et al. 2021). The function reports, for every layer and
#' for both ends of its range:
#'
#' * `gap_low`, `gap_high`: the distance between the most extreme record and
#'   the most extreme conditions of M, as a share of M's range (0 = the
#'   records reach the limit of M).
#' * `edge_low`, `edge_high`: how often records fall in the outer `edge`
#'   (default 5%) of M's conditions, relative to how common those conditions
#'   are in M. About 0 means that the species avoids the margin, so its
#'   limit lies inside M; 1 means that it uses the margin as often as M
#'   offers it; above 1, that it prefers it.
#' * `truncated_low`, `truncated_high`: `edge_*` at or above `edge_ratio`.
#'
#' With `projections` it also gives the share of cells of each projection
#' that lie outside the conditions of M and outside those of the records,
#' and marks as `risk_low` / `risk_high` the layers that are truncated at an
#' end beyond which the projection has more than `min_share` of its cells.
#' For those, consider a larger M (so that the limit is sampled), another
#' variable, or clamping and a cautious reading of the projection.
#'
#' The flags are a screening aid, not a test: a species may truly live up
#' to the edge of M because M ends at a coast or a mountain crest. Records
#' concentrated near an edge because of sampling bias inflate `edge_*`;
#' run it on thinned records ([thin_by_cell()], [thin_environment()]).
#'
#' **What next.** [truncation_remedy()] takes this result and layers
#' extending beyond M, and proposes what to do: extend M towards the
#' missing conditions, look for records in donor regions, or accept the
#' truncation; it also prepares the datasets of a spatially nested model.
#' With records of the whole range, [niche_coverage()] measures how much of
#' the niche a region misses.
#'
#' @param points Records: a data.frame, CSV path, `sf` or
#'   [terra::SpatVector].
#' @param env Predictor layers of the calibration area, e.g. from
#'   [crop_to_m()] ([terra::SpatRaster] or file paths).
#' @param projections Optional named list of rasters to project to (other
#'   periods, scenarios or regions) with the same layer names.
#' @param edge Share of M's conditions that counts as its margin.
#' @param edge_ratio `edge_*` at or above which an end is called truncated.
#' @param min_share Share of a projection's cells beyond M above which a
#'   truncated end is a risk.
#' @param n Number of cells sampled from M and from each projection.
#' @param coords Names of the coordinate columns, if not recognised.
#' @param points_crs CRS of the coordinates.
#' @param seed Seed for the cell samples.
#' @return A `niche_truncation` object: a list with `variables` (one row
#'   per layer), `projections` (one row per projection and layer, or
#'   `NULL`), `ranges` (the quantiles drawn by [plot_niche_truncation()])
#'   and `settings`.
#' @references Thuiller, W., Brotons, L., Araujo, M.B., Lavorel, S. (2004)
#'   Effects of restricting environmental range of data to project current
#'   and future species distributions. *Ecography* 27, 165-172. Chevalier,
#'   M., Broennimann, O., Cornuault, J., Guisan, A. (2021) Data integration
#'   methods to account for spatial niche truncation effects in regional
#'   projections of species distribution. *Ecological Applications* 31,
#'   e02427.
#' @seealso [plot_niche_truncation()], [truncation_remedy()],
#'   [niche_coverage()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' # a calibration area that is too small: only the western part, where
#' # env1 stops at 0.7 although the species goes on to the east
#' m <- terra::crop(sim$env, terra::ext(14, 21, 49, 55))
#' tr <- niche_truncation(sim$presence, m,
#'                        projections = list(whole_area = sim$env))
#' tr
#' plot_niche_truncation(tr)
niche_truncation <- function(points, env, projections = NULL, edge = 0.05,
                             edge_ratio = 1, min_share = 0.01, n = 1e5,
                             coords = NULL, points_crs = "EPSG:4326",
                             seed = 1) {
  r <- as_grid(env)
  vars <- names(r)
  if (edge <= 0 || edge >= 0.5) stop("`edge` must be between 0 and 0.5.",
                                     call. = FALSE)
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  if (is.data.frame(points) && !inherits(points, "sf")) {
    xy_cols <- if (is.null(coords)) guess_coords(names(points)) else coords
    x <- as_number(points[[xy_cols[1]]])
    y <- as_number(points[[xy_cols[2]]])
    ok <- !is.na(x) & !is.na(y)
    p <- terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else {
    p <- as_points(points, coords, points_crs)
  }
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  E <- env_at_cells(r, terra::cellFromXY(r, terra::crds(p)))
  E <- E[stats::complete.cases(E), , drop = FALSE]
  if (nrow(E) < 10) {
    stop("Fewer than 10 records with values in all layers.", call. = FALSE)
  }
  cell_sample <- function(x) {
    valid <- valid_cells(x)
    if (!length(valid)) stop("A raster has no cells with values in all ",
                             "layers.", call. = FALSE)
    s <- with_seed(seed, sort(valid[sample.int(length(valid),
                                               min(length(valid), n))]))
    M <- read_cells(x, s)
    colnames(M) <- names(x)
    M
  }
  M <- cell_sample(r)
  probs <- c(0, 0.05, 0.5, 0.95, 1)
  qtab <- function(X, what) {
    do.call(rbind, lapply(vars, function(v) {
      q <- stats::quantile(X[, v], probs, names = FALSE)
      data.frame(layer = v, data = what, min = q[1], q05 = q[2],
                 median = q[3], q95 = q[4], max = q[5])
    }))
  }
  ranges <- rbind(qtab(M, "M"), qtab(E, "records"))

  variables <- do.call(rbind, lapply(vars, function(v) {
    m <- M[, v]
    e <- E[, v]
    span <- diff(range(m))
    q <- stats::quantile(m, c(edge, 1 - edge), names = FALSE)
    # records in the margin of M, relative to how common the margin is in M
    # (its share is `edge`, or more when many cells share the same value)
    low <- mean(e <= q[1]) / mean(m <= q[1])
    high <- mean(e >= q[2]) / mean(m >= q[2])
    data.frame(layer = v, m_min = min(m), m_max = max(m), rec_min = min(e),
               rec_max = max(e),
               gap_low = if (span > 0) max(0, min(e) - min(m)) / span else 0,
               gap_high = if (span > 0) max(0, max(m) - max(e)) / span else 0,
               edge_low = low, edge_high = high,
               truncated_low = low >= edge_ratio,
               truncated_high = high >= edge_ratio)
  }))

  proj <- NULL
  if (length(projections)) {
    if (is.null(names(projections)) || any(names(projections) == "")) {
      stop("`projections` must be a named list.", call. = FALSE)
    }
    proj <- do.call(rbind, lapply(names(projections), function(nm) {
      x <- as_grid(projections[[nm]])
      miss <- setdiff(vars, names(x))
      if (length(miss)) {
        stop("Projection '", nm, "' lacks layer(s): ",
             paste(miss, collapse = ", "), call. = FALSE)
      }
      P <- cell_sample(x[[vars]])
      ranges <<- rbind(ranges, qtab(P, nm))
      do.call(rbind, lapply(seq_along(vars), function(i) {
        v <- vars[i]
        d <- data.frame(projection = nm, layer = v,
                        below_m = mean(P[, v] < variables$m_min[i]),
                        above_m = mean(P[, v] > variables$m_max[i]),
                        below_records = mean(P[, v] < variables$rec_min[i]),
                        above_records = mean(P[, v] > variables$rec_max[i]))
        d$risk_low <- variables$truncated_low[i] & d$below_m > min_share
        d$risk_high <- variables$truncated_high[i] & d$above_m > min_share
        d
      }))
    }))
  }
  rownames(ranges) <- NULL
  out <- list(variables = variables, projections = proj, ranges = ranges,
              settings = list(edge = edge, edge_ratio = edge_ratio,
                              min_share = min_share, n_records = nrow(E),
                              n_cells = nrow(M)))
  class(out) <- "niche_truncation"
  out
}

#' @export
print.niche_truncation <- function(x, digits = 3, ...) {
  st <- x$settings
  cat("Niche truncation: ", st$n_records, " records, ", st$n_cells,
      " cells of M; margin = outer ", 100 * st$edge, "% of M's conditions\n",
      sep = "")
  v <- x$variables
  ends <- c(if (any(v$truncated_low))
    paste(v$layer[v$truncated_low], "(low end)"),
    if (any(v$truncated_high))
      paste(v$layer[v$truncated_high], "(high end)"))
  cat("  Truncated: ", if (length(ends)) paste(ends, collapse = ", ") else
    "none", "\n\n", sep = "")
  num <- vapply(v, is.double, logical(1))
  v[num] <- lapply(v[num], signif, digits = digits)
  print(v, row.names = FALSE, ...)
  if (!is.null(x$projections)) {
    p <- x$projections
    risky <- p[p$risk_low | p$risk_high, , drop = FALSE]
    cat("\nProjections (share of cells outside M and outside the records)\n")
    cat("  Truncated ends that a projection goes beyond: ",
        if (nrow(risky)) paste0(risky$layer, " in ", risky$projection,
                                collapse = ", ") else "none", "\n\n",
        sep = "")
    num <- vapply(p, is.double, logical(1))
    p[num] <- lapply(p[num], signif, digits = digits)
    print(p, row.names = FALSE, ...)
  }
  invisible(x)
}

#' Plot of a niche truncation report
#'
#' For each layer, the conditions of the calibration area (M), of the
#' records and of each projection on one scale (0 to 1 = the range of M):
#' thin lines span the minimum to the maximum, thick lines the 5th to the
#' 95th percentile, points mark the medians. Records whose thick line
#' reaches an end of M are truncated there; a projection that extends
#' beyond 0 or 1 asks the model for conditions M does not contain.
#'
#' @param x A `niche_truncation` object from [niche_truncation()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- terra::crop(sim$env, terra::ext(14, 21, 49, 55))
#' tr <- niche_truncation(sim$presence, m,
#'                        projections = list(whole_area = sim$env))
#' plot_niche_truncation(tr)
plot_niche_truncation <- function(x, base_size = 12,
                                  file = NULL, width = 9, height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
  if (!inherits(x, "niche_truncation")) {
    stop("`x` must come from niche_truncation().", call. = FALSE)
  }
  d <- x$ranges
  v <- x$variables
  i <- match(d$layer, v$layer)
  span <- v$m_max[i] - v$m_min[i]
  span[span == 0] <- 1
  for (k in c("min", "q05", "median", "q95", "max")) {
    d[[k]] <- (d[[k]] - v$m_min[i]) / span
  }
  lev <- unique(d$data)
  d$data <- factor(d$data, levels = rev(lev))
  d$layer <- factor(d$layer, levels = v$layer)
  flag <- ifelse(v$truncated_low & v$truncated_high, "truncated at both ends",
                 ifelse(v$truncated_low, "truncated at the low end",
                        ifelse(v$truncated_high, "truncated at the high end",
                               "")))
  lab <- stats::setNames(ifelse(flag == "", as.character(v$layer),
                                paste0(v$layer, ": ", flag)), v$layer)
  cols <- c(M = "#8a8a86", records = "#1f1f1c")
  extra <- setdiff(lev, names(cols))
  pal <- c("#c8553d", "#2e6f9e", "#3f8f5f", "#b0812a", "#7a5ea8", "#4f9aa3")
  cols <- c(cols, stats::setNames(rep_len(pal, length(extra)), extra))
  ggplot2::ggplot(d, ggplot2::aes(y = .data$data, colour = .data$data)) +
    ggplot2::geom_vline(xintercept = c(0, 1), colour = "#b9b8b3",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$min, xend = .data$max,
                                       yend = .data$data), linewidth = 0.5) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$q05, xend = .data$q95,
                                       yend = .data$data), linewidth = 2.2) +
    ggplot2::geom_point(ggplot2::aes(x = .data$median), colour = "white",
                        size = 1.1) +
    ggplot2::facet_wrap(~layer, ncol = 1, labeller = ggplot2::as_labeller(lab)) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::labs(x = "Position within the range of M (0 = minimum, 1 = maximum)",
                  y = NULL,
                  subtitle = "Thin: minimum to maximum; thick: 5th to 95th percentile; dot: median") +
    theme_enm(base_size, "#1f1f1c", "#6b6b66") +
    ggplot2::theme(strip.text = ggplot2::element_text(hjust = 0,
                                                      face = "bold"))
}

# Remedies for a truncated niche ---------------------------------------------------

#' What to do about a niche truncated by the calibration area
#'
#' [niche_truncation()] tells whether the records crowd against the limits
#' of M; `truncation_remedy()` says what to do about it, using predictor
#' layers that extend beyond M. A model can only learn the part of a
#' response that its calibration data contain, and models transferred to
#' other regions or times lose accuracy with distance and extrapolation
#' (Rousseau & Betts 2022), with niche change (Liu et al. 2022), and for
#' invasive species in their new ranges (Nguyen & Leung 2022). Two remedies
#' use data from beyond the region: integrating data from other regions
#' (Chevalier et al. 2021) and spatially nested hierarchical models, in
#' which a model of the whole range informs the regional one (Goicolea et
#' al. 2025). This function proposes and prepares; it fits nothing.
#'
#' 1. **The missing tail.** For each truncated end (from `trunc`), the
#'    conditions just beyond M's limit: a band of `tail` times M's range
#'    on that variable, with the other variables within the range of the
#'    records (so that the conditions are a plausible extension of the
#'    niche, not any cell beyond the limit).
#' 2. **Extension.** Cells in that band are sampled from `env` (`n` random
#'    cells). The extension reaches the smallest distance from M at which
#'    each of `tail_bins` parts of every band has at least `min_cells`
#'    sampled cells; the proposed extent is the convex hull of M and the
#'    band cells within that distance (so M grows only towards the missing
#'    conditions), within the area of `env`. The report gives the distance
#'    and how many times larger the extent is than M.
#' 3. **Donor regions.** Cells of the bands outside the proposed extent,
#'    or all of them when no extension covers the bands, are grouped into
#'    regions (cells within `donor_km` of each other), with their size,
#'    distance from M and the ends they supply: places to look for records
#'    of the species.
#' 4. **Advice, per truncated end.** Extend the extent when the missing
#'    conditions are adjacent: the band begins within `max_gap_km` of M
#'    (default a tenth of the radius of a circle of M's area) and the
#'    extension covering it is at most `max_increase` times M's area; otherwise integrate the donor regions;
#'    when the conditions occur nowhere in `env` (e.g. only in a future
#'    climate), accept the truncation and state it. The main remedy is
#'    "donors" if any end needs them, else "extend", else "accept".
#' 5. **Nested design.** `nested` holds a *parent* (whole range: the
#'    records of `parent_records`, or of `records` when that is all there
#'    is, within the extended extent and the donor regions, to be modelled
#'    at a coarser resolution, `parent_fact` cells aggregated, with
#'    climatic variables) and a *child* (the records in M at full
#'    resolution, with all predictors). Fit the parent model first, project
#'    it onto the child's grid, and either use the projection as a
#'    covariate of the child model ("covariate" strategy) or combine the
#'    two predictions by their geometric mean ("multiply"); the covariate
#'    strategy was the more accurate of the two (Goicolea et al. 2025).
#'    [write_truncation_remedy()] writes both datasets and the extents.
#'
#' With records of the whole range, [niche_coverage()] measures how much of
#' the niche the region misses (unfilling) and whether a nested model is
#' worth it; this function answers where the missing conditions are and
#' prepares the data.
#'
#' @param trunc A `niche_truncation` object from [niche_truncation()].
#' @param records The records (as given to [niche_truncation()]).
#' @param env Predictor layers covering a region larger than M (same layer
#'   names), e.g. the continent; [terra::SpatRaster] or files.
#' @param m The calibration area: an `m_area` from [delimit_m()] or polygons.
#' @param parent_records Optional records of the whole range (e.g. GBIF),
#'   for the parent dataset.
#' @param tail Width of the missing band, as a share of M's range.
#' @param tail_bins,min_cells Parts of the band that must each have
#'   `min_cells` sampled cells in the proposed extent.
#' @param max_increase Largest area of the proposed extent, as a multiple of
#'   M's, for the extension to be advised.
#' @param max_gap_km Largest distance between M and the nearest cell of the
#'   missing conditions for them to count as adjacent (default: a tenth of
#'   the radius of a circle with M's area).
#' @param donor_km Distance linking cells into one donor region.
#' @param parent_fact Aggregation factor suggested for the parent's layers.
#' @param n Random cells sampled from `env`.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Random seed.
#' @return An object of class `truncation_remedy`: `remedy` (the main
#'   remedy: `"none"`, `"extend"`, `"donors"` or `"accept"`), `ends` (one
#'   row per truncated end: the band, its cells in `env`, whether the
#'   extension covers it, and the remedy for that end),
#'   `extension` (distance, area, increase, the extent as `sfc`), `donors`
#'   (one row per region) and `donor_cells` (`sf` points), `nested` (parent
#'   and child datasets) and `settings`. See [plot_truncation_remedy()].
#' @references Chevalier, M., Broennimann, O., Cornuault, J., Guisan, A.
#'   (2021) Data integration methods to account for spatial niche
#'   truncation effects in regional projections of species distribution.
#'   *Ecological Applications* 31, e02427. \doi{10.1002/eap.2427}
#'
#'   Goicolea, T., Adde, A., Broennimann, O., Garcia-Vinas, J.I., Gaston, A.,
#'   Aroca-Fernandez, M.J., Guisan, A., Mateo, R.G. (2025) Spatially-nested
#'   hierarchical species distribution models to overcome niche truncation
#'   in national-scale studies. *Ecography* 2025, e07328.
#'   \doi{10.1111/ecog.07328}
#'
#'   Rousseau, J.S., Betts, M.G. (2022) Factors influencing transferability
#'   in species distribution models. *Ecography* 2022, e06060.
#'   \doi{10.1111/ecog.06060}
#'
#'   Liu, C., Wolter, C., Courchamp, F., Roura-Pascual, N., Jeschke, J.M.
#'   (2022) Biological invasions reveal how niche change affects the
#'   transferability of species distribution models. *Ecology* 103, e3719.
#'   \doi{10.1002/ecy.3719}
#'
#'   Nguyen, D., Leung, B. (2022) How well do species distribution models
#'   predict occurrences in exotic ranges? *Global Ecology and
#'   Biogeography* 31, 1051-1065. \doi{10.1111/geb.13482}
#' @seealso [niche_truncation()], [niche_coverage()], [delimit_m()]
#' @export
#' @examples
#' \donttest{
#' # a west-east temperature gradient; the species lives at 6-14 degrees
#' r <- terra::rast(nrows = 30, ncols = 100, xmin = 0, xmax = 20,
#'                  ymin = 45, ymax = 51, crs = "EPSG:4326")
#' xy <- terra::crds(r, na.rm = FALSE)
#' set.seed(1)
#' env <- c(terra::rast(r, vals = xy[, 1] + rnorm(nrow(xy), 0, 0.1)),
#'          terra::rast(r, vals = xy[, 2]))
#' names(env) <- c("temp", "north")
#' pres <- data.frame(lon = runif(150, 6, 10), lat = runif(150, 46, 50))
#' # a calibration area that ends at 10 degrees east
#' m <- sf::st_sfc(sf::st_polygon(list(rbind(c(2, 45), c(10, 45), c(10, 51),
#'                                           c(2, 51), c(2, 45)))),
#'                 crs = 4326)
#' tr <- niche_truncation(pres, crop_to_m(env, m))
#' rem <- truncation_remedy(tr, pres, env, m)
#' rem
#' plot_truncation_remedy(rem)
#' }
truncation_remedy <- function(trunc, records, env, m, parent_records = NULL,
                              tail = 0.25, tail_bins = 4, min_cells = 5,
                              max_increase = 3, max_gap_km = NULL,
                              donor_km = 100,
                              parent_fact = 4, n = 2e5, coords = NULL,
                              points_crs = "EPSG:4326", seed = 1) {
  if (!inherits(trunc, "niche_truncation")) {
    stop("`trunc` must come from niche_truncation().", call. = FALSE)
  }
  r <- as_grid(env)
  v <- trunc$variables
  miss <- setdiff(v$layer, names(r))
  if (length(miss)) {
    stop("`env` lacks layer(s) of `trunc`: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  r <- r[[v$layer]]
  g <- if (inherits(m, "m_area")) m$m else as_sfc(m, "m")
  g <- union_geom(sf::st_transform(g, 4326))
  # records (lon/lat)
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  to_ll <- function(x) {
    p <- if (is.data.frame(x) && !inherits(x, "sf")) {
      df <- as.data.frame(x, stringsAsFactors = FALSE)
      xy <- if (is.null(coords)) guess_coords(names(df)) else coords
      lo <- as_number(df[[xy[1]]])
      la <- as_number(df[[xy[2]]])
      ok <- !is.na(lo) & !is.na(la)
      terra::vect(cbind(lo[ok], la[ok]), crs = points_crs)
    } else as_points(x, coords, points_crs)
    xy <- terra::crds(terra::project(p, "EPSG:4326"))
    data.frame(lon = xy[, 1], lat = xy[, 2])
  }
  rec <- to_ll(records)
  ea <- equal_area_crs(rec$lon, rec$lat)
  g_ea <- sf::st_transform(g, ea)
  area_m <- geom_area(g_ea) / 1e6

  ends <- rbind(
    data.frame(layer = v$layer[v$truncated_low], end = rep("low",
                                                           sum(v$truncated_low))),
    data.frame(layer = v$layer[v$truncated_high], end = rep("high",
                                                            sum(v$truncated_high))))
  settings <- list(tail = tail, tail_bins = tail_bins, min_cells = min_cells,
                   max_increase = max_increase, donor_km = donor_km,
                   parent_fact = parent_fact, m_area_km2 = area_m)
  if (!nrow(ends)) {
    return(structure(list(remedy = "none", ends = ends, extension = NULL,
                          donors = NULL, donor_cells = NULL, nested = NULL,
                          settings = settings),
                     class = "truncation_remedy"))
  }

  # sample of env with values and positions
  smp <- suppressWarnings(suppressMessages(
    sample_background(r, n = n, seed = seed)))
  X <- env_at_cells(r, smp$cell)
  ok <- stats::complete.cases(X)
  smp <- smp[ok, , drop = FALSE]
  X <- X[ok, , drop = FALSE]
  pts <- sf::st_as_sf(smp[c("lon", "lat")], coords = c("lon", "lat"),
                      crs = 4326)
  pts_ea <- sf::st_transform(pts, ea)
  inside <- lengths(sf::st_intersects(pts_ea, g_ea)) > 0
  # records' range of every variable (for "plausible" combinations)
  rr <- v[c("layer", "rec_min", "rec_max")]
  in_rec_range <- function(skip) {
    keep <- rep(TRUE, nrow(X))
    for (k in setdiff(v$layer, skip)) {
      lo <- rr$rec_min[rr$layer == k]
      hi <- rr$rec_max[rr$layer == k]
      sp <- (hi - lo) * 0.05
      keep <- keep & X[, k] >= lo - sp & X[, k] <= hi + sp
    }
    keep
  }
  dist_km <- rep(0, nrow(X))
  out_i <- which(!inside)
  ends$band_lo <- NA_real_
  ends$band_hi <- NA_real_
  ends$cells <- 0L
  tail_rows <- list()
  for (i in seq_len(nrow(ends))) {
    k <- ends$layer[i]
    mm <- v[v$layer == k, ]
    span <- mm$m_max - mm$m_min
    band <- if (ends$end[i] == "high") c(mm$m_max, mm$m_max + tail * span)
            else c(mm$m_min - tail * span, mm$m_min)
    ends$band_lo[i] <- band[1]
    ends$band_hi[i] <- band[2]
    hit <- which(!inside & X[, k] > band[1] - 1e-12 & X[, k] <= band[2] &
                   in_rec_range(k))
    if (ends$end[i] == "low") {
      hit <- which(!inside & X[, k] >= band[1] & X[, k] < band[2] + 1e-12 &
                     in_rec_range(k))
    }
    ends$cells[i] <- length(hit)
    if (length(hit)) {
      part <- pmin(tail_bins, 1 + floor((X[hit, k] - band[1]) /
                                          (diff(band) / tail_bins)))
      tail_rows[[i]] <- data.frame(row = hit, end = i, part = part)
    }
  }
  tails <- do.call(rbind, tail_rows)
  ext <- NULL
  donors <- NULL
  donor_cells <- NULL
  if (!is.null(tails) && nrow(tails)) {
    u <- unique(tails$row)
    d <- as.numeric(sf::st_distance(pts_ea[u, ], g_ea)) / 1000
    dist_km[u] <- d
    tails$dist <- dist_km[tails$row]
    # smallest buffer at which every part of every band has min_cells cells
    need <- expand.grid(end = which(ends$cells > 0), part = seq_len(tail_bins))
    dneed <- vapply(seq_len(nrow(need)), function(j) {
      dd <- sort(tails$dist[tails$end == need$end[j] &
                              tails$part == need$part[j]])
      if (length(dd) >= min_cells) dd[min_cells] else Inf
    }, numeric(1))
    covered_end <- tapply(is.finite(dneed), need$end, all)
    ends$covered <- FALSE
    ends$covered[as.integer(names(covered_end))] <- as.vector(covered_end)
    d_star <- if (all(is.finite(dneed))) max(dneed) else
      if (any(is.finite(dneed))) max(dneed[is.finite(dneed)]) else NA
    if (is.finite(d_star)) {
      # directional: the hull of M and the band cells within d_star, so
      # that M grows only towards the missing conditions
      near <- unique(tails$row[tails$dist <= d_star + 1e-9])
      buf <- union_geom(sf::st_convex_hull(sf::st_union(c(
        g_ea, sf::st_geometry(pts_ea[near, ])))))
      # within the area the layers cover
      # the layers' extent, with vertices along its edges in the layers'
      # own coordinates (so that parallels stay parallels when projected)
      x0 <- terra::xmin(r)
      x1 <- terra::xmax(r)
      y0 <- terra::ymin(r)
      y1 <- terra::ymax(r)
      sx <- seq(x0, x1, length.out = 100)
      sy <- seq(y0, y1, length.out = 100)
      ring <- rbind(cbind(sx, y0), cbind(x1, sy), cbind(rev(sx), y1),
                    cbind(x0, rev(sy)))
      ring <- rbind(ring, ring[1, ])
      bb <- sf::st_sfc(sf::st_polygon(list(unname(ring))),
                       crs = sf::st_crs(terra::crs(r)))
      bb <- sf::st_make_valid(sf::st_transform(bb, ea))
      buf <- union_geom(c(safe_intersection(buf, bb), g_ea))
      a_ext <- geom_area(buf) / 1e6
      ext <- list(distance_km = d_star, area_km2 = a_ext,
                  increase = a_ext / area_m,
                  complete = all(is.finite(dneed)),
                  extent = sf::st_transform(buf, 4326))
    }
    # ends the extension resolves: covered within max_increase
    # adjacent: the band begins next to M (no gap of other conditions)
    if (is.null(max_gap_km)) max_gap_km <- 0.1 * sqrt(area_m / pi)
    gap <- tapply(tails$dist, tails$end, min)
    ends$gap_km <- NA_real_
    ends$gap_km[as.integer(names(gap))] <- as.vector(gap)
    ext_ok <- !is.null(ext) && ext$increase <= max_increase
    ends$remedy <- ifelse(ends$cells == 0, "accept",
                          ifelse(ends$covered & ext_ok &
                                   !is.na(ends$gap_km) &
                                   ends$gap_km <= max_gap_km,
                                 "extend", "donors"))
    # donor regions: band cells of the ends the extension does not resolve
    far <- tails[ends$remedy[tails$end] == "donors", , drop = FALSE]
    if (nrow(far)) {
      fu <- unique(far$row)
      xy <- km_space(smp$lon[fu], smp$lat[fu])
      cl <- if (length(fu) > 1) dbscan::dbscan(xy, eps = donor_km,
                                               minPts = 1)$cluster else 1L
      sz <- tabulate(cl)
      keep_cl <- which(sz >= min_cells)
      if (length(keep_cl)) {
        near <- vapply(keep_cl, function(cc) min(dist_km[fu[cl == cc]]),
                       numeric(1))
        keep_cl <- keep_cl[order(near)]
        donors <- do.call(rbind, lapply(seq_along(keep_cl), function(j) {
          rows <- fu[cl == keep_cl[j]]
          e <- unique(far$end[far$row %in% rows])
          data.frame(donor = j, cells = length(rows),
                     lon = stats::median(smp$lon[rows]),
                     lat = stats::median(smp$lat[rows]),
                     distance_km = min(dist_km[rows]),
                     supplies = paste(paste(ends$layer[e], ends$end[e]),
                                      collapse = ", "),
                     stringsAsFactors = FALSE)
        }))
        lab <- match(cl, keep_cl)
        dc <- data.frame(lon = smp$lon[fu], lat = smp$lat[fu], donor = lab)
        dc <- dc[!is.na(dc$donor), , drop = FALSE]
        donor_cells <- sf::st_as_sf(dc, coords = c("lon", "lat"),
                                    crs = 4326)
      }
    }
    # a band with cells but no donor region large enough: accept
    if (is.null(donors)) ends$remedy[ends$remedy == "donors"] <- "accept"
  } else {
    ends$covered <- FALSE
    ends$gap_km <- NA_real_
    ends$remedy <- "accept"
  }
  settings$max_gap_km <- if (is.null(max_gap_km)) 0.1 * sqrt(area_m / pi)
                         else max_gap_km
  remedy <- if (any(ends$remedy == "donors")) "donors" else
    if (any(ends$remedy == "extend")) "extend" else "accept"

  # nested design: parent (whole range, coarse) and child (M, fine)
  pr <- if (is.null(parent_records)) rec else to_ll(parent_records)
  parts <- list(g)
  if (!is.null(ext) && any(ends$remedy == "extend")) {
    parts <- c(parts, list(ext$extent))
  }
  if (!is.null(donor_cells)) {
    db <- sf::st_buffer(sf::st_transform(sf::st_union(donor_cells), ea),
                        donor_km * 500)
    parts <- c(parts, list(sf::st_transform(db, 4326)))
  }
  parent_extent <- union_geom(do.call(c, parts))
  in_parent <- lengths(sf::st_intersects(
    sf::st_as_sf(pr, coords = c("lon", "lat"), crs = 4326),
    parent_extent)) > 0
  in_child <- lengths(sf::st_intersects(
    sf::st_as_sf(rec, coords = c("lon", "lat"), crs = 4326), g)) > 0
  nested <- list(
    parent = list(records = pr[in_parent, , drop = FALSE],
                  extent = parent_extent, aggregation = parent_fact,
                  from_parent_records = !is.null(parent_records)),
    child = list(records = rec[in_child, , drop = FALSE], extent = g))
  structure(list(remedy = remedy, ends = ends, extension = ext,
                 donors = donors, donor_cells = donor_cells,
                 nested = nested, settings = settings),
            class = "truncation_remedy")
}

#' @export
print.truncation_remedy <- function(x, digits = 3, ...) {
  s <- x$settings
  e <- x$ends
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(round(v), format = "d", big.mark = " ")
  cat("Remedy for a truncated niche: ", nrow(e), " truncated end(s)\n\n",
      "IN SHORT\n", sep = "")
  if (x$remedy == "none") {
    say("No variable is truncated: no remedy is needed.")
    return(invisible(x))
  }
  say(paste0("Truncated: ", paste0(e$layer, " (", e$end, " end)",
                                   collapse = ", "), "."))
  ext <- x$extension
  lab <- function(i) paste0(e$layer[i], " (", e$end[i], ")",
                            collapse = ", ")
  ie <- which(e$remedy == "extend")
  idn <- which(e$remedy == "donors")
  ia <- which(e$remedy == "accept")
  if (length(ie)) {
    say(paste0("Extend the calibration area for ", lab(ie), ": the missing ",
               "conditions lie next to M. Extending M by up to ",
               f(ext$distance_km), " km towards them covers the band of ",
               s$tail * 100,
               "% of M's range beyond ", if (length(ie) > 1) "these ends"
               else "this end", ", at ",
               formatC(ext$increase, format = "f", digits = 1),
               " times the area of M (", f(ext$area_km2), " km2)."))
  }
  if (length(idn)) {
    d <- x$donors
    k <- seq_len(min(3, nrow(d)))
    say(paste0("Integrate donor regions for ", lab(idn), ": ",
               "the missing conditions begin ",
               f(min(e$gap_km[idn], na.rm = TRUE)), " km from M (more than ",
               f(s$max_gap_km), " km), or no small extension covers them, ",
               "so they are not adjacent",
               ". The conditions occur in ", nrow(d), " region(s) elsewhere: ",
               "look for records of the species there and fit a nested ",
               "model (x$nested). Nearest: ", paste0(
                 "region ", d$donor[k], " at ", f(d$distance_km[k]),
                 " km (", d$supplies[k], ")", collapse = "; "), "."))
  }
  if (length(ia)) {
    say(paste0("Accept the truncation of ", lab(ia), " and state it: the ",
               "conditions beyond ", if (length(ia) > 1) "these ends" else
                 "this end", " occur nowhere (or too rarely) in the layers ",
               "given, e.g. only in a future climate or beyond the layers. ",
               "The response there is unknown; projections into those ",
               "conditions extrapolate, so clamp or mark them ",
               "(extrapolation_risk())."))
  }
  if (x$remedy != "extend" && !x$nested$parent$from_parent_records) {
    say(paste0("The parent dataset holds only the records given; add ",
               "records of the whole range (parent_records, e.g. GBIF) ",
               "before fitting a nested model."))
  }
  say(paste0("With records of the whole range, niche_coverage() measures ",
             "how much of the niche the region misses."))
  cat("\n")
  show <- e[c("layer", "end", "band_lo", "band_hi", "cells", "gap_km",
              "covered", "remedy")]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], signif, digits)
  print(show, row.names = FALSE)
  if (!is.null(x$donors)) {
    cat("\nDonor regions:\n")
    d <- x$donors
    d$lon <- round(d$lon, 2)
    d$lat <- round(d$lat, 2)
    d$distance_km <- round(d$distance_km)
    print(d, row.names = FALSE)
  }
  cat("", strwrap(paste0(
    "band = conditions within ", s$tail * 100, "% of M's range beyond the ",
    "truncated end, with the other variables in the records' range; cells ",
    "= sampled cells of the layers in the band outside M; covered = every ",
    "one of ", s$tail_bins, " parts of the band has ", s$min_cells,
    " cells within the extension; gap_km = distance from M to the nearest ",
    "cell of the band. The limits (", s$max_increase, " times M's area, ",
    "a gap of ", round(s$max_gap_km), " km) are rules of thumb. Nested design (x$nested): fit the ",
    "parent on the whole-range records at a coarser resolution, project ",
    "it to the child's grid, and use it as a covariate of the child model ",
    "(or combine the two predictions by their geometric mean)."),
    width = 76), sep = "\n")
  invisible(x)
}

#' Map the remedies for a truncated niche
#'
#' M, the proposed extension and the sampled cells of the missing
#' conditions, coloured by donor region, from [truncation_remedy()].
#'
#' @param x A `truncation_remedy` object.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Goicolea, T. et al. (2025) Spatially-nested hierarchical
#'   species distribution models to overcome niche truncation in
#'   national-scale studies. *Ecography* 2025, e07328.
#'   \doi{10.1111/ecog.07328}
#' @export
#' @examples
#' # see ?truncation_remedy
plot_truncation_remedy <- function(x, base_size = 11, file = NULL,
                                   width = 7, height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "truncation_remedy")) {
    stop("`x` must come from truncation_remedy().", call. = FALSE)
  }
  if (x$remedy == "none") stop("No truncation: nothing to map.",
                               call. = FALSE)
  p <- ggplot2::ggplot()
  ext_used <- !is.null(x$extension) && any(x$ends$remedy == "extend")
  if (ext_used) {
    p <- p + ggplot2::geom_sf(data = sf::st_sf(geometry = x$extension$extent),
                              fill = "#7570b3", alpha = 0.15,
                              colour = "#7570b3", linetype = "dashed")
  }
  p <- p + ggplot2::geom_sf(data = sf::st_sf(geometry = x$nested$child$extent),
                            fill = "#1b9e77", alpha = 0.25,
                            colour = "#1b9e77")
  if (!is.null(x$donor_cells)) {
    dc <- x$donor_cells
    cnt <- table(dc$donor)
    lv <- paste0("Region ", names(cnt), ": ", as.integer(cnt), " cells")
    dc$region <- factor(lv[match(dc$donor, as.integer(names(cnt)))], lv)
    pal <- grDevices::hcl.colors(max(3, length(lv)), "Dark 3")[seq_along(lv)]
    p <- p + ggplot2::geom_sf(data = dc, ggplot2::aes(colour = .data$region),
                              size = 0.8) +
      ggplot2::scale_colour_manual(values = pal, name = "Donor regions")
  }
  rec <- x$nested$child$records
  p + ggplot2::geom_point(data = rec, ggplot2::aes(.data$lon, .data$lat),
                          size = 0.6, colour = "#3a3a38") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Where the missing conditions are",
                  subtitle = paste0(
                    "Green: M with its records",
                    if (ext_used) "; dashed: proposed extension" else "",
                    if (!is.null(x$donor_cells))
                      ";\ncoloured points: cells of the missing conditions"
                    else "")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
}

#' Write the parent and child datasets of a nested design
#'
#' Writes the records of the parent (whole range) and the child (M) as CSV
#' files, their extents and the donor regions as a GeoPackage, and a short
#' README with what to do with them.
#'
#' @param x A `truncation_remedy` object.
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @return Invisibly, the paths written.
#' @references Goicolea, T. et al. (2025) Spatially-nested hierarchical
#'   species distribution models to overcome niche truncation in
#'   national-scale studies. *Ecography* 2025, e07328.
#'   \doi{10.1111/ecog.07328}
#' @export
#' @examples
#' # see ?truncation_remedy
write_truncation_remedy <- function(x, dir, prefix = "nested") {
  if (!inherits(x, "truncation_remedy")) {
    stop("`x` must come from truncation_remedy().", call. = FALSE)
  }
  if (x$remedy == "none") stop("No truncation: nothing to write.",
                               call. = FALSE)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, "_", s))
  out <- c(parent = f("parent_records.csv"), child = f("child_records.csv"),
           gpkg = f("extents.gpkg"), readme = f("README.txt"))
  utils::write.csv(x$nested$parent$records, out[["parent"]],
                   row.names = FALSE)
  utils::write.csv(x$nested$child$records, out[["child"]], row.names = FALSE)
  if (file.exists(out[["gpkg"]])) unlink(out[["gpkg"]])
  sf::st_write(sf::st_sf(part = "parent",
                         geometry = x$nested$parent$extent),
               out[["gpkg"]], layer = "parent_extent", quiet = TRUE)
  sf::st_write(sf::st_sf(part = "child", geometry = x$nested$child$extent),
               out[["gpkg"]], layer = "child_extent", quiet = TRUE,
               append = TRUE)
  if (!is.null(x$donor_cells)) {
    sf::st_write(x$donor_cells, out[["gpkg"]], layer = "donor_cells",
                 quiet = TRUE, append = TRUE)
  }
  writeLines(c(
    "Spatially nested design (Goicolea et al. 2025, Ecography e07328)",
    "",
    paste0("Remedy advised: ", x$remedy),
    "1. Parent: fit a model on parent_records.csv within parent_extent,",
    paste0("   with climatic predictors aggregated by about ",
           x$nested$parent$aggregation, " cells (coarser resolution)."),
    "2. Project the parent model onto the child's grid (child_extent).",
    "3. Child: fit a model on child_records.csv within child_extent at full",
    "   resolution, with all predictors and the parent projection as a",
    "   covariate ('covariate' strategy), or fit it without and combine the",
    "   two predictions by their geometric mean ('multiply').",
    "4. Report both strategies; the covariate one was the more accurate in",
    "   Goicolea et al. (2025)."), out[["readme"]])
  invisible(out)
}
