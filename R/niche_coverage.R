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
#' * `unfilling`: the share of the species' niche (records of the whole
#'   range in the space of the first two principal components) that the
#'   regional records do not occupy. This is the measure of Goicolea et al.
#'   (2025), there computed with `ecospat.niche.similarity.test()`; here the
#'   environmental space is cut into `n_bins` x `n_bins` cells and the share
#'   of whole-range records in cells without a regional record is counted.
#' * `not_available`: the share of whole-range records in conditions that do
#'   not occur in the region at all. No regional model can learn these;
#'   they matter when the region moves into them (e.g. under warming).
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
#' Use thinned records of the whole range (e.g. GBIF, filtered with
#' [filter_occurrences()] and [thin_by_cell()]) and the same climate
#' variables as for the regional model; the layers may be coarser.
#' [niche_truncation()] then checks, for the chosen calibration area,
#' along which predictors and at which end the response is cut.
#'
#' @param points Records of the whole range (data.frame, CSV path, `sf` or
#'   [terra::SpatVector]).
#' @param env Layers covering the whole range ([terra::SpatRaster] or
#'   files).
#' @param region The region of interest: polygons (`sf`, `SpatVector`,
#'   file) or a raster (cells with values).
#' @param probs Quantiles of the whole-range records that define their range
#'   per variable (default the central 98%, so that a few outliers do not
#'   decide it).
#' @param n_bins Number of bins per principal component for `unfilling`.
#' @param limits Two `unfilling` values separating "regional model
#'   sufficient", "nested model worth considering" and "nested model
#'   strongly advised".
#' @param n Number of cells sampled from the region.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @param seed Seed for the cell sample.
#' @return A `niche_coverage` object: a list with `summary` (records in
#'   total and in the region, `unfilling`, `not_available`, `advice`),
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
niche_coverage <- function(points, env, region, probs = c(0.01, 0.99),
                           n_bins = 20, limits = c(0.1, 0.5), n = 1e5,
                           coords = NULL, points_crs = "EPSG:4326",
                           seed = 1) {
  r <- as_grid(env)
  vars <- names(r)
  if (length(limits) != 2 || any(limits < 0 | limits > 1) ||
      limits[1] >= limits[2]) {
    stop("`limits` must be two increasing shares between 0 and 1.",
         call. = FALSE)
  }
  # records of the whole range and their values
  xy <- points_xy(points, r, coords, points_crs)
  cell <- terra::cellFromXY(r, xy)
  ok <- !is.na(cell)
  cell <- cell[ok]
  E <- read_cells(r, cell)
  colnames(E) <- vars
  keep <- stats::complete.cases(E)
  E <- E[keep, , drop = FALSE]
  cell <- cell[keep]
  if (nrow(E) < 10) {
    stop("Fewer than 10 records with values in all layers.", call. = FALSE)
  }
  inside <- in_mask(r, cell, region)
  if (sum(inside) < 3) {
    stop("Fewer than 3 records in `region`.", call. = FALSE)
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
  Er <- E[inside, , drop = FALSE]
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
  b_reg <- unique(b_all[inside])
  b_env <- unique(stats::na.omit(bin_of(SR)))
  unfilling <- mean(!b_all %in% b_reg)
  not_available <- mean(!b_all %in% b_env)
  advice <- if (unfilling < limits[1]) {
    "regional model sufficient"
  } else if (unfilling < limits[2]) {
    "nested model worth considering"
  } else "nested model strongly advised"

  summary <- data.frame(records = nrow(E), records_region = sum(inside),
                        share_region = mean(inside), unfilling = unfilling,
                        not_available = not_available, advice = advice,
                        stringsAsFactors = FALSE)
  structure(list(summary = summary, variables = per_var, pca = pca,
                 scores = data.frame(pc1 = S[, 1],
                                     pc2 = if (k > 1) S[, 2] else 0,
                                     region = inside),
                 region_scores = data.frame(pc1 = SR[, 1],
                                            pc2 = if (k > 1) SR[, 2] else 0),
                 settings = list(probs = probs, n_bins = n_bins,
                                 limits = limits)),
            class = "niche_coverage")
}

#' @export
print.niche_coverage <- function(x, digits = 3, ...) {
  s <- x$summary
  cat("Niche coverage of the region: ", s$records_region, " of ", s$records,
      " records (", round(100 * s$share_region), "%) in the region\n\n",
      sep = "")
  v <- x$variables[c("variable", "coverage", "outside_records",
                     "outside_region")]
  v[-1] <- lapply(v[-1], round, digits)
  print(v, row.names = FALSE)
  pct <- function(z) paste0(round(100 * z), "%")
  cat("\nIN SHORT\n")
  txt <- c(
    paste0("Unfilling: ", pct(s$unfilling), " of the species' niche (records ",
           "of the whole range, first two principal components) has no ",
           "regional record; ", pct(s$not_available), " lies in conditions ",
           "that do not occur in the region at all."),
    paste0("Advice: ", s$advice, ". ", switch(
      s$advice,
      "regional model sufficient" = paste0(
        "A model calibrated in the region sees nearly the whole niche."),
      "nested model worth considering" = paste0(
        "A regional model misses part of the niche; combine it with a ",
        "model of the whole range (nested / hierarchical model), above all ",
        "for projections."),
      "nested model strongly advised" = paste0(
        "A regional model sees less than half of the niche; its response ",
        "curves are truncated and projections extrapolate. Calibrate a ",
        "model of the whole range and nest the regional model in it."))),
    if (s$unfilling - s$not_available > 0.1) {
      paste0("Part of the unfilled niche (", pct(s$unfilling -
                                                   s$not_available),
             ") exists in the region but has no records there: check for ",
             "undersampling before blaming the extent.")
    })
  for (t in txt) {
    cat(strwrap(t, width = 78, exdent = 2, prefix = "  ", initial = "- "),
        sep = "\n")
  }
  invisible(x)
}

#' Plot the niche coverage of a region
#'
#' The environmental space of the first two principal components of the
#' records of the whole range: the conditions of the region (shaded), the
#' records of the whole range (grey) and those in the region (blue). Grey
#' points outside the shaded area are conditions the region does not have;
#' grey points inside it without blue points nearby are conditions the
#' region has but where the species has no regional records.
#'
#' @param x Result of [niche_coverage()].
#' @param base_size Base font size.
#' @return A ggplot object.
#' @export
#' @examples
#' # see niche_coverage()
plot_niche_coverage <- function(x, base_size = 11) {
  if (!inherits(x, "niche_coverage")) {
    stop("`x` must come from niche_coverage().", call. = FALSE)
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  s <- x$scores
  s$group <- factor(ifelse(s$region, "records in the region",
                           "records elsewhere"),
                    levels = c("records elsewhere", "records in the region"))
  s <- s[order(s$region), , drop = FALSE]
  ve <- summary(x$pca)$importance[2, ]
  ax <- function(j) sprintf("PC%d (%.0f%% of the variance)", j, 100 * ve[j])
  ggplot2::ggplot(s, ggplot2::aes(.data$pc1, .data$pc2)) +
    ggplot2::geom_bin_2d(data = x$region_scores, bins = 40,
                         fill = "#e9e2cf", colour = NA) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$group), size = 1.4,
                        alpha = 0.8) +
    ggplot2::scale_colour_manual(values = c(
      "records elsewhere" = "#9a9a95", "records in the region" = "#1c5cab"),
      name = NULL) +
    ggplot2::labs(
      x = ax(1), y = ax(min(2, length(ve))),
      title = "How much of the niche does the region contain?",
      subtitle = sprintf(paste0("Shaded: conditions of the region.\n",
                                "Unfilling %.0f%%, not available in the ",
                                "region %.0f%%: %s."),
                         100 * x$summary$unfilling,
                         100 * x$summary$not_available,
                         x$summary$advice)) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom")
}
