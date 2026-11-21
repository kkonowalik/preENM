# Regional study of a widespread species: is a nested model needed? --------------

#' How much of a species' niche does a region contain?
#'
#' Many studies model a widespread species only in a region of interest, e.g.
#' a country. A model calibrated there learns only the conditions the
#' species meets in that region. If the species also lives in conditions
#' that the region lacks or where it has no records, the response curves
#' are truncated, and projections (to the future, or to warmer or drier
#' parts of the region) extrapolate. Spatially nested (hierarchical) models
#' address this by combining a model of the whole range, often at a coarser
#' resolution, with the regional model (Pearson et al. 2004; Mateo et al.
#' 2019; Chevalier et al. 2021; Goicolea et al. 2025). This function tells
#' you how much a regional model would miss, from the records of the whole
#' range and layers covering it:
#'
#' * `unfilling`: the share of the species' niche that the regional records
#'   do not occupy, the measure used by Goicolea et al. (2025). The space of
#'   the first two principal components of the whole-range records is cut
#'   into `n_bins` x `n_bins` cells; the niche is the cells that hold at
#'   least `min_records` whole-range records, and `unfilling` is the share
#'   of these cells without a regional record. **Every cell counts once**,
#'   however many records it holds.
#' * `not_available`: the share of the niche cells whose conditions do not
#'   occur in the region at all. No regional model can learn these; they
#'   matter when the region moves into them (e.g. under warming).
#' * `unfilling_records`, `not_available_records`: the same counted by
#'   records, i.e. the share of whole-range records in such cells. These
#'   are biased when most records come from the region (often the case:
#'   regions of interest are usually the best sampled): they can never
#'   exceed the share of records outside the region, so a large part of the
#'   niche known from few records elsewhere hardly counts. They are reported
#'   for comparison; the advice uses the niche cells. Sparse records still
#'   occupy fewer cells than dense ones would, so even `unfilling` is
#'   conservative where the rest of the range is poorly sampled.
#' * Per variable: the share of the whole-range records outside the range
#'   of the regional records (`outside_records`) and outside the conditions
#'   of the region (`outside_region`), and the part of the range of the
#'   whole-range records covered by the regional ones (`coverage`).
#'
#' The printed advice uses `unfilling`: below the first of `limits`
#' (default 10%) a regional model sees nearly the whole niche; between the
#' limits a nested model is worth considering, above the second (50%, "very
#' poorly filled" in Goicolea et al. 2025) it is strongly advised. A low
#' `unfilling` with a high `not_available` share cannot happen; a high
#' `unfilling` with a low `not_available` share means that the conditions
#' exist in the region but the species has no records there (absent,
#' undersampled, or excluded by other factors), which a nested model does
#' not solve by itself.
#'
#' **The records and layers must cover the whole range.** Records outside
#' the layers have no values and are left out (their number is reported);
#' if only the regional records remain, the region trivially contains
#' their niche, and the function warns. Use thinned records of the whole
#' range (e.g. GBIF, filtered with [filter_occurrences()] and
#' [thin_by_cell()]) and the same climate variables as for the regional
#' model, for the whole range; the layers may be coarser.
#' [niche_truncation()] then checks, for the chosen calibration area,
#' along which predictors and at which end the response is cut.
#'
#' @param points Records of the whole range (data.frame, CSV path, `sf` or
#'   [terra::SpatVector]).
#' @param env Layers covering the whole range ([terra::SpatRaster] or
#'   files).
#' @param region The region of interest: polygons (`sf`, `SpatVector`,
#'   file) or a raster (cells with values), e.g. a country or the M of a
#'   regional model.
#' @param region_points Optional regional records from another source, e.g.
#'   a national survey (the usual case for nested models, where the
#'   regional model uses more precise regional data). Then `unfilling` and
#'   the per-variable values compare these with `points`; by default the
#'   records of `points` that fall in `region` are the regional records.
#' @param probs Quantiles of the whole-range records that define their range
#'   per variable (default the central 98%, so that a few outliers do not
#'   decide it).
#' @param n_bins Number of bins per principal component for `unfilling`.
#' @param min_records Records a cell of the environmental space needs to
#'   count as part of the niche; e.g. 2 ignores cells with a single,
#'   possibly erroneous record.
#' @param limits Two `unfilling` values separating "regional model
#'   sufficient", "nested model worth considering" and "nested model
#'   strongly advised".
#' @param n Number of cells sampled from the region.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @param seed Seed for the cell sample.
#' @param file Optional file(s) to save the plot to (see
#'   [plot_niche_coverage()]), e.g. `c("coverage.png", "coverage.pdf")`.
#' @param width,height Size of the saved plot in inches.
#' @return A `niche_coverage` object: a list with `summary` (records in
#'   total and in the region, `niche_cells`, `unfilling`, `not_available`,
#'   `unfilling_records`, `not_available_records`, `advice`), `bins` (one
#'   row per niche cell: centre, records, whether it has regional records
#'   and occurs in the region),
#'   `variables` (one row per layer), `pca` (the PCA of the whole-range
#'   records) and the scores used by [plot_niche_coverage()].
#' @references Goicolea, T., Adde, A., Broennimann, O., Garcia-Vinas, J.I.,
#'   Gaston, A., Aroca-Fernandez, M.J., Guisan, A., Mateo, R.G. (2025)
#'   Spatially-nested hierarchical species distribution models to overcome
#'   niche truncation in national-scale studies. *Ecography* 2025, e07328.
#'   \doi{10.1111/ecog.07328}
#'
#'   Chevalier, M., Broennimann, O., Cornuault, J., Guisan, A. (2021) Data
#'   integration methods to account for spatial niche truncation effects in
#'   regional projections of species distribution. *Ecological
#'   Applications* 31, e02427.
#'
#'   Mateo, R.G., Gaston, A., Aroca-Fernandez, M.J., Broennimann, O.,
#'   Guisan, A., Saura, S., Garcia-Vinas, J.I. (2019) Hierarchical species
#'   distribution models in support of vegetation conservation at the
#'   landscape scale. *Journal of Vegetation Science* 30, 386-396.
#'
#'   Pearson, R.G., Dawson, T.P., Liu, C. (2004) Modelling species
#'   distributions in Britain: a hierarchical integration of climate and
#'   land-cover data. *Ecography* 27, 285-298.
#' @seealso [niche_truncation()], [plot_niche_coverage()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 200)
#' # the region of interest: the western third of the range
#' west <- terra::as.polygons(terra::ext(14, 17, 49, 55), crs = "EPSG:4326")
#' nc <- niche_coverage(sim$presence, sim$env, west)
#' nc
#' plot_niche_coverage(nc)
#' # plot_niche_coverage(nc, file = c("coverage.png", "coverage.pdf"))
niche_coverage <- function(points, env, region, region_points = NULL,
                           probs = c(0.01, 0.99), n_bins = 20,
                           min_records = 1,
                           limits = c(0.1, 0.5), n = 1e5, coords = NULL,
                           points_crs = "EPSG:4326", seed = 1, file = NULL,
                           width = 7.5, height = 5.5) {
  r <- as_grid(env)
  vars <- names(r)
  if (length(limits) != 2 || any(limits < 0 | limits > 1) ||
      limits[1] >= limits[2]) {
    stop("`limits` must be two increasing shares between 0 and 1.",
         call. = FALSE)
  }
  # records and their values; records outside the layers are counted
  values_at <- function(p) {
    xy <- points_xy(p, r, coords, points_crs)
    cell <- terra::cellFromXY(r, xy)
    n_in <- nrow(xy)
    cell <- cell[!is.na(cell)]
    E <- read_cells(r, cell)
    colnames(E) <- vars
    keep <- stats::complete.cases(E)
    list(E = E[keep, , drop = FALSE], cell = cell[keep],
         dropped = n_in - sum(keep), total = n_in)
  }
  w <- values_at(points)
  E <- w$E
  cell <- w$cell
  if (nrow(E) < 10) {
    stop("Fewer than 10 records with values in all layers (", w$dropped,
         " of ", w$total, " records lie outside the layers or have missing ",
         "values). The layers must cover the whole range of the species.",
         call. = FALSE)
  }
  if (w$dropped > 0) {
    message(w$dropped, " of ", w$total, " records lie outside the layers or ",
            "have missing values and were left out.")
  }
  inside <- in_mask(r, cell, region)
  if (is.null(region_points)) {
    if (sum(inside) < 3) {
      stop("Fewer than 3 records in `region`.", call. = FALSE)
    }
    if (all(inside)) {
      warning(
        "All ", nrow(E), " records lie in `region`, so the region trivially ",
        "contains their niche. Give records of the whole range of the ",
        "species (e.g. from GBIF) and layers that cover it",
        if (w$dropped > 0) paste0("; ", w$dropped, " records were left out ",
                                  "because the layers do not cover them"),
        ".", call. = FALSE)
    }
    Er <- E[inside, , drop = FALSE]
  } else {
    wr <- values_at(region_points)
    Er <- wr$E
    if (nrow(Er) < 3) {
      stop("Fewer than 3 `region_points` with values in all layers.",
           call. = FALSE)
    }
    if (wr$dropped > 0) {
      message(wr$dropped, " of ", wr$total, " regional records lie outside ",
              "the layers or have missing values and were left out.")
    }
  }
  # conditions of the region
  valid <- valid_cells(r)
  in_reg <- valid[in_mask(r, valid, region)]
  if (!length(in_reg)) stop("`region` covers no cell with values.",
                            call. = FALSE)
  smp <- with_seed(seed, if (length(in_reg) > n) {
    sort(in_reg[sample.int(length(in_reg), n)])
  } else in_reg)
  R <- read_cells(r, smp)
  colnames(R) <- vars

  # per variable
  per_var <- do.call(rbind, lapply(vars, function(v) {
    g <- stats::quantile(E[, v], probs, names = FALSE)
    rr <- range(Er[, v])
    re <- range(R[, v])
    span <- diff(g)
    cover <- if (span > 0) {
      max(0, min(g[2], rr[2]) - max(g[1], rr[1])) / span
    } else 1
    data.frame(variable = v, range_low = g[1], range_high = g[2],
               region_records_low = rr[1], region_records_high = rr[2],
               coverage = cover,
               outside_records = mean(E[, v] < rr[1] | E[, v] > rr[2]),
               outside_region = mean(E[, v] < re[1] | E[, v] > re[2]),
               stringsAsFactors = FALSE)
  }))

  # joint niche: first two principal components of the whole-range records
  pca <- stats::prcomp(E, scale. = TRUE)
  k <- min(2, ncol(E))
  S <- stats::predict(pca, E)[, seq_len(k), drop = FALSE]
  SR <- stats::predict(pca, R)[, seq_len(k), drop = FALSE]
  lim <- apply(S, 2, range)
  bin_of <- function(M) {
    b <- vapply(seq_len(k), function(j) {
      br <- seq(lim[1, j], lim[2, j], length.out = n_bins + 1)
      out <- findInterval(M[, j], br, rightmost.closed = TRUE,
                          all.inside = TRUE)
      out[M[, j] < lim[1, j] | M[, j] > lim[2, j]] <- NA
      out
    }, numeric(nrow(M)))
    b <- matrix(b, ncol = k)
    id <- do.call(paste, as.data.frame(b))
    id[!stats::complete.cases(b)] <- NA_character_
    id
  }
  b_all <- bin_of(S)
  SRec <- stats::predict(pca, Er)[, seq_len(k), drop = FALSE]
  b_reg <- unique(stats::na.omit(bin_of(SRec)))
  b_env <- unique(stats::na.omit(bin_of(SR)))
  # by niche space: every cell of the environmental space with at least
  # `min_records` whole-range records counts once, however many records it
  # holds, so densely sampled parts of the range do not dominate
  tab <- table(b_all)
  occ <- names(tab)[tab >= min_records]
  if (!length(occ)) stop("No cell of the environmental space has ",
                         min_records, " records; lower `min_records`.",
                         call. = FALSE)
  unfilling <- mean(!occ %in% b_reg)
  not_available <- mean(!occ %in% b_env)
  # by records: the share of whole-range records in such cells
  unfilling_records <- mean(!b_all %in% b_reg)
  not_available_records <- mean(!b_all %in% b_env)
  cell_xy <- do.call(rbind, strsplit(occ, " "))
  width_pc <- (lim[2, ] - lim[1, ]) / n_bins
  bins <- data.frame(
    pc1 = lim[1, 1] + (as.numeric(cell_xy[, 1]) - 0.5) * width_pc[1],
    pc2 = if (k > 1) lim[1, 2] + (as.numeric(cell_xy[, 2]) - 0.5) *
      width_pc[2] else 0,
    records = as.integer(tab[occ]), regional = occ %in% b_reg,
    available = occ %in% b_env)
  advice <- if (is.null(region_points) && all(inside)) {
    "cannot be judged"
  } else if (unfilling < limits[1]) {
    "regional model sufficient"
  } else if (unfilling < limits[2]) {
    "nested model worth considering"
  } else "nested model strongly advised"

  summary <- data.frame(records = nrow(E), records_dropped = w$dropped,
                        records_region = sum(inside),
                        share_region = mean(inside),
                        region_points = if (is.null(region_points)) NA else
                          nrow(Er),
                        niche_cells = length(occ), unfilling = unfilling,
                        not_available = not_available,
                        unfilling_records = unfilling_records,
                        not_available_records = not_available_records,
                        advice = advice, stringsAsFactors = FALSE)
  out <- structure(
    list(summary = summary, variables = per_var, pca = pca,
         scores = data.frame(pc1 = S[, 1], pc2 = if (k > 1) S[, 2] else 0,
                             region = inside),
         region_scores = data.frame(pc1 = SR[, 1],
                                    pc2 = if (k > 1) SR[, 2] else 0),
         region_point_scores = if (is.null(region_points)) NULL else
           data.frame(pc1 = SRec[, 1], pc2 = if (k > 1) SRec[, 2] else 0),
         bins = bins, bin_width = width_pc,
         settings = list(probs = probs, n_bins = n_bins, limits = limits,
                         min_records = min_records,
                         separate_region_points = !is.null(region_points))),
    class = "niche_coverage")
  if (!is.null(file)) plot_niche_coverage(out, file = file, width = width,
                                          height = height)
  out
}

#' @export
print.niche_coverage <- function(x, digits = 3, ...) {
  s <- x$summary
  sep <- isTRUE(x$settings$separate_region_points)
  pct <- function(z) paste0(round(100 * z), "%")
  cat("How much of the species' niche does the region contain?\n\n")
  cat("Records of the whole range: ", s$records,
      if (!is.null(s$records_dropped) && s$records_dropped > 0) {
        paste0(" (", s$records_dropped, " more left out: outside the layers ",
               "or missing values)")
      }, "\n", sep = "")
  cat("  of these inside the region: ", s$records_region, " (",
      pct(s$share_region), ")\n", sep = "")
  if (sep) cat("Regional records (`region_points`): ", s$region_points, "\n",
               sep = "")
  reg <- if (sep) "the regional records" else "the records inside the region"
  cat("\nPer variable:\n")
  v <- x$variables[c("variable", "coverage", "outside_records",
                     "outside_region")]
  v[-1] <- lapply(v[-1], round, digits)
  print(v, row.names = FALSE)
  cat(strwrap(paste0(
    "coverage = part of the range of the whole-range records that ", reg,
    " cover; outside_records = share of whole-range records beyond the range ",
    "of ", reg, "; outside_region = share beyond the conditions of the ",
    "region. Each variable is judged alone, so these shares are smaller ",
    "than those of the niche below, which counts combinations of ",
    "conditions."), width = 78, prefix = "  "), sep = "\n")
  cat("\nIN SHORT\n")
  n_un <- round(s$unfilling * s$niche_cells)
  out_share <- if (sep) NA else 1 - s$share_region
  txt <- c(
    paste0(pct(s$unfilling), " of the species' niche has no record in the ",
           "region: ", n_un, " of the ", s$niche_cells, " cells of ",
           "environmental space that the species occupies (unfilling)."),
    paste0(pct(s$not_available), " of the niche lies in conditions that do ",
           "not occur anywhere in the region (not available)."),
    paste0("Counted by records instead, ", pct(s$unfilling_records),
           " of the whole-range records lie in the unfilled part",
           if (!is.na(out_share)) paste0(
             "; this share cannot exceed the share of records outside the ",
             "region (", pct(out_share), ")"),
           ". It understates the gap when most records come from the ",
           "region, so the advice uses the niche space."),
    paste0("Advice: ", s$advice, ". ", switch(
      s$advice,
      "regional model sufficient" = paste0(
        "A model calibrated in the region sees nearly the whole niche."),
      "nested model worth considering" = paste0(
        "A regional model misses part of the niche; combine it with a ",
        "model of the whole range (nested, or hierarchical, model), above ",
        "all for projections."),
      "nested model strongly advised" = paste0(
        "A regional model sees less than half of the niche; its response ",
        "curves are truncated and projections extrapolate. Calibrate a ",
        "model of the whole range and nest the regional model in it."),
      "cannot be judged" = paste0(
        "All records lie in the region, so the comparison says nothing. ",
        "Give records of the whole range of the species and layers that ",
        "cover it."))),
    if (s$unfilling - s$not_available > 0.1) {
      paste0(pct(s$unfilling - s$not_available), " of the niche exists in ",
             "the region but has no records there: check for undersampling ",
             "before blaming the extent.")
    },
    NULL)
  for (t in txt) {
    cat(strwrap(t, width = 78, exdent = 2, prefix = "  ", initial = "- "),
        sep = "\n")
  }
  cat(strwrap(paste0(
    "Niche: the cells of the space of the first two principal components ",
    "of the variables (cut into ", x$settings$n_bins, " x ",
    x$settings$n_bins, ") that hold at least ", x$settings$min_records,
    " whole-range record(s); each counts once."), width = 78,
    prefix = "  "), sep = "\n")
  invisible(x)
}

#' Plot the niche coverage of a region
#'
#' The environmental space of the first two principal components of the
#' records of the whole range: the conditions of the region (beige), the
#' niche cells without regional records (red frames; their share is
#' `unfilling`), the
#' whole-range records outside the region (grey) and inside it (blue), and
#' separate regional records if given (`region_points`, red crosses); the
#' legend gives their numbers. Grey points outside the shaded area are
#' conditions the region does not have; grey points inside it without blue
#' points nearby are conditions the region has but where the species has no
#' regional records.
#'
#' @param x Result of [niche_coverage()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @export
#' @examples
#' # see niche_coverage()
plot_niche_coverage <- function(x, base_size = 11,
                                file = NULL, width = 7.5, height = 5.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
  if (!inherits(x, "niche_coverage")) {
    stop("`x` must come from niche_coverage().", call. = FALSE)
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  s <- x$scores
  n_in <- sum(s$region)
  n_out <- sum(!s$region)
  lab_in <- sprintf("whole-range records inside the region (%d)", n_in)
  lab_out <- sprintf("whole-range records outside the region (%d)", n_out)
  lab_reg <- sprintf("regional records (%d)",
                     if (is.null(x$region_point_scores)) 0L else
                       nrow(x$region_point_scores))
  s$group <- factor(ifelse(s$region, lab_in, lab_out),
                    levels = c(lab_out, lab_in, lab_reg))
  s <- s[order(s$region), , drop = FALSE]
  if (!is.null(x$region_point_scores)) {
    rp <- x$region_point_scores
    rp$group <- factor(lab_reg, levels = levels(s$group))
    rp$region <- NA
    s <- rbind(s, rp[names(s)])
  }
  cols <- stats::setNames(c("#9a9a95", "#1c5cab", "#c8553d"),
                          c(lab_out, lab_in, lab_reg))
  shapes <- stats::setNames(c(16, 16, 4), c(lab_out, lab_in, lab_reg))
  # conditions of the region as tiles of a 40 x 40 grid
  rs <- x$region_scores
  rg <- rbind(range(c(rs$pc1, x$scores$pc1)), range(c(rs$pc2, x$scores$pc2)))
  rt_w <- pmax((rg[, 2] - rg[, 1]) / 40, 1e-9)
  ix <- floor((rs$pc1 - rg[1, 1]) / rt_w[1])
  iy <- floor((rs$pc2 - rg[2, 1]) / rt_w[2])
  rt <- unique(data.frame(ix, iy))
  rt <- data.frame(pc1 = rg[1, 1] + (rt$ix + 0.5) * rt_w[1],
                   pc2 = rg[2, 1] + (rt$iy + 0.5) * rt_w[2])
  b <- x$bins
  if (!is.null(b)) {
    b <- b[!b$regional, , drop = FALSE]
    lab_gap <- sprintf("niche cells without regional records (%d of %d)",
                       nrow(b), nrow(x$bins))
    b$what <- factor(lab_gap, levels = lab_gap)
  }
  ve <- summary(x$pca)$importance[2, ]
  ax <- function(j) sprintf("PC%d (%.0f%% of the variance)", j, 100 * ve[j])
  ggplot2::ggplot(s, ggplot2::aes(.data$pc1, .data$pc2)) +
    ggplot2::geom_tile(data = rt, inherit.aes = FALSE,
                       ggplot2::aes(.data$pc1, .data$pc2), width = rt_w[1],
                       height = rt_w[2], fill = "#e9e2cf", colour = NA) +
    (if (!is.null(b) && nrow(b)) {
      ggplot2::geom_tile(data = b, inherit.aes = FALSE,
                         ggplot2::aes(.data$pc1, .data$pc2, fill = .data$what),
                         width = x$bin_width[1],
                         height = if (length(x$bin_width) > 1)
                           x$bin_width[2] else 1,
                         colour = "#c8553d", linewidth = 0.25, alpha = 0.35)
    }) +
    (if (!is.null(b) && nrow(b)) {
      ggplot2::scale_fill_manual(values = stats::setNames("#f2b8a8",
                                                          levels(b$what)),
                                 name = NULL)
    }) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group,
                                     shape = .data$group,
                                     size = .data$group), alpha = 0.8) +
    ggplot2::scale_size_manual(values = stats::setNames(c(1.4, 1.4, 0.9),
                                                        c(lab_out, lab_in,
                                                          lab_reg)),
                               name = NULL, drop = TRUE) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::scale_shape_manual(values = shapes, name = NULL, drop = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1, order = 1),
                    shape = ggplot2::guide_legend(ncol = 1, order = 1),
                    size = ggplot2::guide_legend(ncol = 1, order = 1),
                    fill = ggplot2::guide_legend(order = 2)) +
    ggplot2::labs(
      x = ax(1), y = ax(min(2, length(ve))),
      title = "How much of the niche does the region contain?",
      subtitle = sprintf(paste0("Beige: conditions that occur in the ",
                                "region; red frames: niche cells without ",
                                "regional records.\n%.0f%% of the niche has ",
                                "no record in the region, %.0f%% is not ",
                                "available there.\nAdvice: %s."),
                         100 * x$summary$unfilling,
                         100 * x$summary$not_available,
                         x$summary$advice)) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom")
}
