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
#' @seealso [plot_niche_truncation()]
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
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- terra::crop(sim$env, terra::ext(14, 21, 49, 55))
#' tr <- niche_truncation(sim$presence, m,
#'                        projections = list(whole_area = sim$env))
#' plot_niche_truncation(tr)
plot_niche_truncation <- function(x, base_size = 12) {
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
