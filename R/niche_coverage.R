# Regional study of a widespread species: is a nested model needed? --------------

# Categories of cells of the environmental space (internal)
coverage_levels <- c("filled", "elsewhere, available", "elsewhere, absent",
                     "region only", "available, unoccupied", "unoccupied",
                     "outside")

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
#' range and layers covering it.
#'
#' **The niche.** The records of the whole range are placed in the space of
#' the first two principal components of the variables (after `log_vars`
#' are log-transformed), which is cut into `n_bins` x `n_bins` cells. The
#' niche is the cells holding at least `min_records` whole-range records;
#' **every cell counts once**, however many records it holds, so that the
#' well-sampled region does not dominate.
#'
#' * `unfilling`: the share of the niche cells without a regional record,
#'   the measure used by Goicolea et al. (2025).
#' * `not_available`: the share of the niche cells whose conditions do not
#'   occur in the region at all. No regional model can learn these; they
#'   matter when the region moves into them (e.g. under warming).
#' * `unfilling_records`, `not_available_records`: the same counted by
#'   records. They can never exceed the share of records outside the
#'   region, so with a well-sampled region they understate the gap; given
#'   for comparison only.
#' * Per variable: the part of the range of the whole-range records covered
#'   by the regional ones (`coverage`), and the shares of whole-range
#'   records beyond the range of the regional records (`outside_records`)
#'   and beyond the conditions of the region (`outside_region`).
#'
#' **How serious is the gap?** Two statistics, both in the printout:
#'
#' * A bootstrap interval of `unfilling` (`n_boot` resamples of the
#'   records): how much it depends on the particular records.
#' * A null comparison: the unfilling left by random samples of the
#'   whole-range records of the same size as the regional records
#'   (`n_null` samples). If the region leaves much more unfilled than a
#'   random sample of the same size (large `unfilling_ses`, small
#'   `unfilling_p`), the gap comes from the conditions of the region, not
#'   from its number of records.
#'
#' **Robustness.** Single records in remote cells (outliers, other taxa,
#' errors) can drive `unfilling`. The printout gives it for `min_records`
#' 1, 2, 3 and 5 and says when the advice changes; check the records in
#' unfilled cells with [unfilled_records()] (e.g. their taxonomy).
#'
#' **Advice** (from `unfilling` and `limits`), graded: below 10% a small
#' gap (a regional model sees nearly the whole niche); 10-25% a moderate gap
#' (a nested model worth considering for projections); 25-50% a large gap
#' (nested model advised); above 50%, "very poorly filled" in Goicolea et
#' al. (2025), a very large gap (strongly advised). Only the last limit is
#' from the paper. A nested model matters most for projections: future (or
#' past) conditions of the region that lie outside the regional niche but
#' inside the conditions the species occupies elsewhere; a regional model
#' predicts losses there that the whole-range model does not.
#' [classify_points()] shows how much of the region is in that situation.
#' This assumes no strong local adaptation of the regional populations.
#'
#' **Memory.** Values are read only at the records and at `n` random
#' locations in the region, so large (e.g. global 30 arc-second) layers can
#' be used.
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
#' @param log_vars Variables log-transformed (`log1p`) before the PCA, e.g.
#'   skewed precipitation amounts: names, `"auto"` for bioclim
#'   precipitation amounts recognised by name (bio12-14, bio16-19), or
#'   `NULL` (default) for none.
#' @param probs Quantiles of the whole-range records that define their range
#'   per variable (default the central 98%, so that a few outliers do not
#'   decide it).
#' @param n_bins Number of bins per principal component, or `"auto"`
#'   (default): 20, coarser with fewer than about 1,200 records (at least
#'   8), so that occupied cells hold several records on average. The
#'   printout also gives the unfilling for a coarser and a finer grid.
#' @param min_records Records a cell needs to count as part of the niche;
#'   e.g. 2 ignores cells with a single, possibly erroneous record.
#' @param limits Three `unfilling` values separating a small gap
#'   ("regional model sufficient"), a moderate one ("nested model worth
#'   considering"), a large one ("nested model advised") and a very large
#'   one ("nested model strongly advised"). Rules of thumb, except 50%,
#'   "very poorly filled" in Goicolea et al. (2025).
#' @param n_boot,n_null Bootstrap resamples and random samples for the
#'   statistics (0 to skip).
#' @param n Number of cells sampled from the region.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @param seed Seed for the cell sample and the statistics.
#' @param file Optional file(s) to save the plot to (see
#'   [plot_niche_coverage()]), e.g. `c("coverage.png", "coverage.pdf")`.
#' @param width,height Size of the saved plot in inches.
#' @return A `niche_coverage` object: a list with `summary` (records,
#'   `niche_cells`, `unfilling` with `unfilling_lo`/`unfilling_hi`,
#'   `unfilling_null`, `unfilling_ses`, `unfilling_p`, `not_available`,
#'   the record-based shares and `advice`), `robustness` (unfilling for
#'   several `min_records`), `variables` (one row per layer), `bins` (one
#'   row per cell with records: centre, records, regional records, regional
#'   density percentile, category), `records` (the whole-range records with
#'   coordinates and category, see [unfilled_records()]) and what
#'   [classify_points()] and [plot_niche_coverage()] need.
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
#' @seealso [plot_niche_coverage()], [classify_points()],
#'   [unfilled_records()], [niche_truncation()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 200)
#' # the region of interest: the western third of the range
#' west <- terra::as.polygons(terra::ext(14, 17, 49, 55), crs = "EPSG:4326")
#' nc <- niche_coverage(sim$presence, sim$env, west, n_boot = 49,
#'                      n_null = 49)
#' nc
#' plot_niche_coverage(nc)
#' head(unfilled_records(nc))
#' # plot_niche_coverage(nc, file = c("coverage.png", "coverage.pdf"))
niche_coverage <- function(points, env, region, region_points = NULL,
                           log_vars = NULL, probs = c(0.01, 0.99),
                           n_bins = "auto", min_records = 1,
                           limits = c(0.1, 0.25, 0.5), n_boot = 199, n_null = 199,
                           n = 1e5, coords = NULL, points_crs = "EPSG:4326",
                           seed = 1, file = NULL, width = 14, height = 6) {
  r <- as_grid(env)
  vars <- names(r)
  limits <- check_limits(limits)
  lv <- if (identical(log_vars, "auto")) default_log_vars(vars) else
    intersect(log_vars, vars)
  miss_lv <- setdiff(log_vars, c(vars, "auto"))
  if (length(miss_lv)) {
    stop("`log_vars` not among the layers: ", paste(miss_lv, collapse = ", "),
         call. = FALSE)
  }
  # records and their values; records outside the layers are counted
  values_at <- function(p) {
    xy <- points_xy(p, r, coords, points_crs)
    cell <- terra::cellFromXY(r, xy)
    n_in <- nrow(xy)
    ok <- !is.na(cell)
    xy <- xy[ok, , drop = FALSE]
    cell <- cell[ok]
    # values at the points only (not block by block over the raster)
    E <- if (length(cell)) as.matrix(terra::extract(r, xy)) else
      matrix(numeric(0), 0, length(vars))
    colnames(E) <- vars
    keep <- stats::complete.cases(E)
    list(E = E[keep, , drop = FALSE], cell = cell[keep],
         xy = xy[keep, , drop = FALSE], dropped = n_in - sum(keep),
         total = n_in)
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
  reg_obj <- coverage_region(region, r)
  inside <- inside_region(w$xy, reg_obj, r)
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
  # conditions of the region: random locations inside it (memory grows
  # with `n`, not with the size of the rasters)
  R <- with_seed(seed, sample_region(r, reg_obj, n))
  if (!nrow(R)) stop("`region` covers no cell with values.", call. = FALSE)
  colnames(R) <- vars

  # per variable (original units)
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

  # the space: PCA of the (transformed) whole-range records, cut into cells
  tf <- function(M) {
    if (length(lv)) {
      if (any(M[, lv] < -1, na.rm = TRUE)) {
        stop("Values below -1 in `log_vars`: ", paste(lv, collapse = ", "),
             call. = FALSE)
      }
      M[, lv] <- log1p(M[, lv])
    }
    M
  }
  pca <- stats::prcomp(tf(E), scale. = TRUE)
  k <- min(2, ncol(E))
  score <- function(M) {
    s <- stats::predict(pca, tf(M))[, seq_len(k), drop = FALSE]
    if (k == 1) s <- cbind(s, 0)
    s
  }
  S <- score(E)
  if (identical(n_bins, "auto")) {
    # 20 x 20, coarser with few records, so that occupied cells hold
    # several records on average
    n_bins <- as.integer(min(20, max(8, round(sqrt(nrow(E) / 3)))))
  }
  if (!is.numeric(n_bins) || n_bins < 3) {
    stop("`n_bins` must be \"auto\" or a number of at least 3.",
         call. = FALSE)
  }
  lim <- apply(S, 2, range)
  grid <- list(lim = lim, n_bins = n_bins,
               width = pmax((lim[2, ] - lim[1, ]) / n_bins, 1e-12))
  cell_all <- grid_cells(S, grid)
  cell_reg <- grid_cells(score(Er), grid)
  cell_env <- grid_cells(score(R), grid)
  nc2 <- n_bins^2
  n_all <- tabulate(cell_all, nc2)
  n_reg <- tabulate(stats::na.omit(cell_reg), nc2)
  avail <- tabulate(stats::na.omit(cell_env), nc2) > 0
  measure <- function(n_all, n_reg, mr = min_records) {
    niche <- n_all >= mr
    if (!any(niche)) return(c(unfilling = NA_real_, not_available = NA_real_))
    c(unfilling = mean(n_reg[niche] == 0),
      not_available = mean(!avail[niche]))
  }
  m <- measure(n_all, n_reg)
  if (is.na(m[["unfilling"]])) {
    stop("No cell of the environmental space has ", min_records,
         " records; lower `min_records`.", call. = FALSE)
  }
  in_gap <- n_reg[cell_all] == 0 & n_all[cell_all] >= min_records
  unfilling_records <- mean(in_gap)
  not_available_records <- mean(!avail[cell_all] &
                                  n_all[cell_all] >= min_records)
  robustness <- data.frame(min_records = sort(unique(c(1, 2, 3, 5,
                                                       min_records))))
  robustness$niche_cells <- vapply(robustness$min_records, function(mr) {
    sum(n_all >= mr)
  }, numeric(1))
  robustness$unfilling <- vapply(robustness$min_records, function(mr) {
    measure(n_all, n_reg, mr)[["unfilling"]]
  }, numeric(1))
  robustness$advice <- vapply(robustness$unfilling, coverage_advice,
                              character(1), limits = limits)
  # the same at coarser and finer grids
  SRg <- score(Er)
  SEg <- score(R)
  resolution <- data.frame(n_bins = sort(unique(c(
    max(3, round(n_bins * 0.66)), n_bins, round(n_bins * 1.5)))))
  resolution$niche_cells <- NA_real_
  resolution$unfilling <- NA_real_
  for (i in seq_len(nrow(resolution))) {
    g <- list(lim = lim, n_bins = resolution$n_bins[i],
              width = pmax((lim[2, ] - lim[1, ]) / resolution$n_bins[i],
                           1e-12))
    na <- tabulate(grid_cells(S, g), g$n_bins^2)
    nr <- tabulate(stats::na.omit(grid_cells(SRg, g)), g$n_bins^2)
    niche <- na >= min_records
    resolution$niche_cells[i] <- sum(niche)
    resolution$unfilling[i] <- mean(nr[niche] == 0)
  }

  # how serious: bootstrap of the records, and random samples of the same
  # size as the regional records
  sep <- !is.null(region_points)
  boot <- null <- numeric(0)
  stats_out <- with_seed(seed, {
    if (n_boot > 0) {
      cr <- if (sep) stats::na.omit(cell_reg) else NULL
      boot <- vapply(seq_len(n_boot), function(b) {
        i <- sample.int(length(cell_all), replace = TRUE)
        na <- tabulate(cell_all[i], nc2)
        nr <- if (sep) {
          tabulate(cr[sample.int(length(cr), replace = TRUE)], nc2)
        } else tabulate(cell_all[i][inside[i]], nc2)
        measure(na, nr)[["unfilling"]]
      }, numeric(1))
    }
    if (n_null > 0) {
      n_r <- min(sum(!is.na(cell_reg)), length(cell_all))
      null <- vapply(seq_len(n_null), function(b) {
        nr <- tabulate(cell_all[sample.int(length(cell_all), n_r)], nc2)
        measure(n_all, nr)[["unfilling"]]
      }, numeric(1))
    }
    list(boot = boot, null = null)
  })
  boot <- stats_out$boot
  null <- stats_out$null
  u <- m[["unfilling"]]
  ses <- if (length(null) > 1 && stats::sd(null) > 0) {
    (u - mean(null)) / stats::sd(null)
  } else NA_real_
  advice <- if (!sep && all(inside)) "cannot be judged" else
    coverage_advice(u, limits)

  # categories of the cells and of the records
  cat_cell <- coverage_category(seq_len(nc2), n_all, n_reg, avail,
                                min_records)
  dens <- coverage_density(n_reg)
  used <- which(n_all > 0 | n_reg > 0)
  centre <- function(id) {
    ix <- (id - 1) %/% n_bins + 1
    iy <- (id - 1) %% n_bins + 1
    cbind(lim[1, 1] + (ix - 0.5) * grid$width[1],
          lim[1, 2] + (iy - 0.5) * grid$width[2])
  }
  cc <- centre(used)
  bins <- data.frame(cell = used, pc1 = cc[, 1], pc2 = cc[, 2],
                     records = n_all[used], regional_records = n_reg[used],
                     regional_density = dens[used], available = avail[used],
                     category = factor(cat_cell[used], coverage_levels))
  ll <- terra::crds(terra::project(terra::vect(w$xy, crs = terra::crs(r)),
                                   "EPSG:4326"))
  records <- data.frame(id = seq_len(nrow(E)), lon = ll[, 1], lat = ll[, 2],
                        pc1 = S[, 1], pc2 = S[, 2], in_region = inside,
                        niche_cell = cell_all,
                        category = factor(cat_cell[cell_all],
                                          coverage_levels))

  summary <- data.frame(
    records = nrow(E), records_dropped = w$dropped,
    records_region = sum(inside), share_region = mean(inside),
    region_points = if (sep) nrow(Er) else NA,
    niche_cells = sum(n_all >= min_records), unfilling = u,
    unfilling_lo = if (length(boot)) stats::quantile(boot, 0.025,
                                                     names = FALSE) else NA,
    unfilling_hi = if (length(boot)) stats::quantile(boot, 0.975,
                                                     names = FALSE) else NA,
    unfilling_null = if (length(null)) mean(null) else NA,
    unfilling_ses = ses,
    unfilling_p = if (length(null)) (1 + sum(null >= u)) /
      (length(null) + 1) else NA,
    not_available = m[["not_available"]],
    unfilling_records = unfilling_records,
    not_available_records = not_available_records,
    advice = advice, stringsAsFactors = FALSE)
  # the region as polygons in longitude/latitude, for the map
  region_ll <- tryCatch(coverage_region_outline(reg_obj),
                        error = function(e) NULL)
  out <- structure(
    list(summary = summary, robustness = robustness,
         resolution = resolution, variables = per_var,
         bins = bins, records = records, pca = pca,
         scores = data.frame(pc1 = S[, 1], pc2 = S[, 2], region = inside),
         region_scores = data.frame(pc1 = score(R)[, 1],
                                    pc2 = score(R)[, 2]),
         region_point_scores = if (sep) {
           sr <- score(Er)
           data.frame(pc1 = sr[, 1], pc2 = sr[, 2])
         },
         region_outline = region_ll,
         fit = list(vars = vars, log_vars = lv, k = k, grid = grid,
                    n_all = n_all, n_reg = n_reg, avail = avail,
                    min_records = min_records, tf = tf, pca = pca,
                    medians = apply(E, 2, stats::median),
                    q = apply(E, 2, stats::quantile, probs = c(0.01, 0.99)),
                    crs = terra::crs(r)),
         boot = boot, null = null,
         settings = list(probs = probs, n_bins = n_bins, limits = limits,
                         min_records = min_records, log_vars = lv,
                         separate_region_points = sep)),
    class = "niche_coverage")
  if (!is.null(file)) plot_niche_coverage(out, file = file, width = width,
                                          height = height)
  out
}

# Cell numbers of scores in a grid (internal); NA outside it
grid_cells <- function(S, grid) {
  lim <- grid$lim
  nb <- grid$n_bins
  ij <- vapply(1:2, function(j) {
    i <- floor((S[, j] - lim[1, j]) / grid$width[j]) + 1
    i[S[, j] == lim[2, j]] <- nb
    i[S[, j] < lim[1, j] | S[, j] > lim[2, j]] <- NA
    i
  }, numeric(nrow(S)))
  ij <- matrix(ij, ncol = 2)
  as.integer((ij[, 1] - 1) * nb + ij[, 2])
}

# Category of grid cells (internal)
coverage_category <- function(id, n_all, n_reg, avail, min_records) {
  out <- rep("outside", length(id))
  ok <- !is.na(id)
  niche <- n_all[id[ok]] >= min_records
  reg <- n_reg[id[ok]] > 0
  av <- avail[id[ok]]
  out[ok] <- ifelse(niche & reg, "filled",
             ifelse(niche & av, "elsewhere, available",
             ifelse(niche, "elsewhere, absent",
             ifelse(reg, "region only",
             ifelse(av, "available, unoccupied", "unoccupied")))))
  out
}

# Percentile of the number of regional records of a cell among the cells
# with regional records (internal); 0 for cells without
coverage_density <- function(n_reg) {
  out <- numeric(length(n_reg))
  f <- n_reg > 0
  if (any(f)) out[f] <- rank(n_reg[f], ties.method = "max") / sum(f)
  out
}

advice_levels <- c("regional model sufficient",
                   "nested model worth considering",
                   "nested model advised", "nested model strongly advised")
gap_levels <- c("small gap", "moderate gap", "large gap", "very large gap")

coverage_advice <- function(u, limits) {
  if (is.na(u)) return(NA_character_)
  advice_levels[findInterval(u, limits) + 1]
}

# Three limits; two (as before version 2.22.0) get a middle one (internal)
check_limits <- function(limits) {
  if (length(limits) == 2) limits <- c(limits[1], mean(limits), limits[2])
  if (length(limits) != 3 || any(limits < 0 | limits > 1) ||
      any(diff(limits) <= 0)) {
    stop("`limits` must be three increasing shares between 0 and 1.",
         call. = FALSE)
  }
  limits
}

# The region as polygons in the CRS of `r`, or as a raster (internal)
coverage_region <- function(region, r) {
  if (is.character(region)) {
    region <- if (grepl("\\.(tif|tiff|asc|nc|vrt|img)$", region[1],
                        ignore.case = TRUE)) terra::rast(region) else
      terra::vect(region)
  }
  if (inherits(region, c("sf", "sfc"))) region <- terra::vect(region)
  if (inherits(region, "SpatVector")) {
    if (!terra::same.crs(region, r)) region <- terra::project(region,
                                                              terra::crs(r))
    return(region)
  }
  if (inherits(region, "SpatRaster")) return(region[[1]])
  stop("`region` must be polygons or a raster (object or file).",
       call. = FALSE)
}

# Which points (x/y in the CRS of `r`) lie in the region (internal)
inside_region <- function(xy, reg, r) {
  if (!nrow(xy)) return(logical(0))
  if (inherits(reg, "SpatVector")) {
    p <- terra::vect(xy, crs = terra::crs(r))
    return(as.vector(terra::is.related(p, reg, "intersects")))
  }
  if (!terra::same.crs(reg, r)) {
    xy <- terra::crds(terra::project(terra::vect(xy, crs = terra::crs(r)),
                                     terra::crs(reg)))
  }
  v <- terra::extract(reg, xy)[[1]]
  !is.na(v) & v != 0
}

# Values of `r` at up to `n` random cells inside the region (internal)
sample_region <- function(r, reg, n) {
  e <- terra::intersect(terra::ext(r), terra::ext(
    if (inherits(reg, "SpatVector")) reg else
      terra::project(terra::as.polygons(terra::ext(reg),
                                        crs = terra::crs(reg)),
                     terra::crs(r))))
  if (is.null(e)) return(matrix(numeric(0), 0, terra::nlyr(r)))
  rows <- range(terra::rowFromY(r, c(e$ymin, e$ymax) +
                                  c(0.5, -0.5) * terra::res(r)[2]),
                na.rm = TRUE)
  cols <- range(terra::colFromX(r, c(e$xmin, e$xmax) +
                                  c(0.5, -0.5) * terra::res(r)[1]),
                na.rm = TRUE)
  n_win <- (diff(rows) + 1) * (diff(cols) + 1)
  out <- list()
  got <- 0
  for (try in 1:6) {
    k <- min(n_win, ceiling(3 * (n - got)) + 100)
    cell <- if (n_win <= 3 * n) {
      terra::cellFromRowColCombine(r, rows[1]:rows[2], cols[1]:cols[2])
    } else {
      unique(terra::cellFromRowCol(
        r, sample.int(diff(rows) + 1, k, replace = TRUE) + rows[1] - 1,
        sample.int(diff(cols) + 1, k, replace = TRUE) + cols[1] - 1))
    }
    xy <- terra::xyFromCell(r, cell)
    xy <- xy[inside_region(xy, reg, r), , drop = FALSE]
    if (nrow(xy)) {
      v <- as.matrix(terra::extract(r, xy))
      v <- v[stats::complete.cases(v), , drop = FALSE]
      out[[length(out) + 1]] <- v
      got <- got + nrow(v)
    }
    if (got >= n || n_win <= 3 * n) break
  }
  v <- do.call(rbind, out)
  if (is.null(v)) return(matrix(numeric(0), 0, terra::nlyr(r)))
  v <- unique(v)
  if (nrow(v) > n) v <- v[sample.int(nrow(v), n), , drop = FALSE]
  v
}

# The region as an sf polygon in longitude/latitude (internal)
coverage_region_outline <- function(reg) {
  if (inherits(reg, "SpatRaster")) {
    a <- if (terra::ncell(reg) > 4e6) {
      terra::aggregate(reg, ceiling(sqrt(terra::ncell(reg) / 4e6)),
                       fun = "max", na.rm = TRUE)
    } else reg
    reg <- terra::as.polygons(!is.na(a) & a != 0)
    reg <- reg[reg[[1]][, 1] == 1]
  }
  v <- terra::project(reg, "EPSG:4326")
  sf::st_as_sf(terra::aggregate(v))
}

#' @export
print.niche_coverage <- function(x, digits = 3, ...) {
  s <- x$summary
  sep <- isTRUE(x$settings$separate_region_points)
  pct <- function(z) paste0(round(100 * z), "%")
  wrap <- function(t, ...) cat(strwrap(t, width = 78, ...), sep = "\n")
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
  wrap(paste0(
    "coverage = part of the range of the whole-range records that ", reg,
    " cover; outside_records = share of whole-range records beyond the range ",
    "of ", reg, "; outside_region = share beyond the conditions of the ",
    "region. Each variable is judged alone, so these shares are smaller ",
    "than those of the niche below, which counts combinations of ",
    "conditions."), prefix = "  ")
  rb <- x$robustness
  rs <- x$resolution
  if (!is.null(rb)) {
    cat("\nUnfilling when a cell needs more records to count as niche:\n")
    show <- rb
    show$unfilling <- pct(show$unfilling)
    print(show, row.names = FALSE)
  }
  if (!is.null(rs)) {
    cat("\nUnfilling with a coarser or finer grid (cells per axis):\n")
    show <- rs
    show$unfilling <- pct(show$unfilling)
    print(show, row.names = FALSE)
  }
  lim <- x$settings$limits
  if (length(lim) == 2) lim <- c(lim[1], mean(lim), lim[2])
  n_un <- round(s$unfilling * s$niche_cells)
  out_share <- if (sep) NA else 1 - s$share_region
  ci <- if (!is.null(s$unfilling_lo) && !is.na(s$unfilling_lo)) {
    paste0(" (95% interval ", pct(s$unfilling_lo), "-", pct(s$unfilling_hi),
           ")")
  } else ""
  null_txt <- if (!is.null(s$unfilling_null) && !is.na(s$unfilling_null)) {
    big <- !is.na(s$unfilling_p) && s$unfilling_p <= 0.05
    paste0(
      "A random sample of the whole-range records of the same size as the ",
      "regional records would leave ", pct(s$unfilling_null), " unfilled",
      if (!is.na(s$unfilling_ses)) sprintf(" (SES = %.1f, p = %s)",
                                           s$unfilling_ses,
                                           signif(s$unfilling_p, 2)), ". ",
      if (big) {
        "The gap comes from the conditions of the region, not from its number of records."
      } else {
        "The gap is about what the number of regional records alone would leave; more regional records may close it."
      })
  }
  # how robust: the range over min_records and grid sizes
  range_txt <- NULL
  us <- c(s$unfilling, rb$unfilling, rs$unfilling)
  us <- us[!is.na(us)]
  if (length(us) > 1 && diff(range(us)) >= 0.02) {
    g <- gap_levels[findInterval(range(us), lim) + 1]
    range_txt <- paste0(
      "Depending on how many records a cell needs (",
      paste(range(rb$min_records), collapse = " to "), ") and on the grid (",
      paste(range(rs$n_bins), collapse = " to "), " cells per axis), the ",
      "unfilling is ", pct(min(us)), " to ", pct(max(us)), ": ",
      if (g[1] == g[2]) paste0("a ", g[1], " in every case.") else
        paste0("from a ", g[1], " to a ", g[2], ". Sparse records in remote ",
               "cells decide the upper value; check them with ",
               "unfilled_records() before deciding."))
  }
  advice_txt <- switch(
    s$advice,
    "regional model sufficient" = paste0(
      "Small gap: a regional model sees nearly the whole niche; a nested ",
      "model is not needed."),
    "nested model worth considering" = paste0(
      "Moderate gap: a regional model misses some of the niche. It is ",
      "usually acceptable for present-day maps; for projections a nested ",
      "model is worth considering."),
    "nested model advised" = paste0(
      "Large gap: a regional model misses a large part of the niche. A ",
      "nested model is advised, above all for projections."),
    "nested model strongly advised" = paste0(
      "Very large gap: a regional model sees less than half of the niche ",
      "(\"very poorly filled\" in Goicolea et al. 2025); its response ",
      "curves are truncated and projections extrapolate. A nested model is ",
      "strongly advised."),
    "cannot be judged" = paste0(
      "All records lie in the region, so the comparison says nothing. Give ",
      "records of the whole range of the species and layers that cover it."),
    "")
  txt <- c(
    paste0(pct(s$unfilling), ci, " of the species' niche has no record in ",
           "the region: ", n_un, " of the ", s$niche_cells, " cells of ",
           "environmental space that the species occupies (unfilling)."),
    null_txt,
    paste0(pct(s$not_available), " of the niche lies in conditions that do ",
           "not occur anywhere in the region (not available)."),
    paste0("Counted by records instead, ", pct(s$unfilling_records),
           " of the whole-range records lie in the unfilled part",
           if (!is.na(out_share)) paste0(
             "; this share cannot exceed the share of records outside the ",
             "region (", pct(out_share), ")"),
           ". It understates the gap when most records come from the ",
           "region, so the advice uses the niche cells."),
    paste0("Advice: ", s$advice, ". ", advice_txt),
    range_txt,
    if (s$unfilling - s$not_available > 0.1) {
      paste0(pct(s$unfilling - s$not_available), " of the niche exists in ",
             "the region but has no records there: check for undersampling ",
             "before blaming the extent.")
    },
    paste0("What a nested model adds: a model of the whole range knows ",
           "conditions that the regional records lack. When the region ",
           "becomes warmer or drier, its future conditions may lie outside ",
           "the regional niche but inside the conditions the species ",
           "occupies elsewhere; a regional model then predicts losses that ",
           "the whole-range model does not. To see how much of the region ",
           "is in that situation, place its projected conditions with ",
           "classify_points(). This assumes that the regional populations ",
           "can live in the conditions of the populations elsewhere, or be ",
           "replaced by them; local adaptation would make the nested model ",
           "too optimistic."))
  cat("\nIN SHORT\n")
  for (t in txt) wrap(t, exdent = 2, prefix = "  ", initial = "- ")
  wrap(paste0(
    "Niche: the cells of the space of the first two principal components ",
    "of the variables",
    if (length(x$settings$log_vars)) paste0(" (",
      paste(x$settings$log_vars, collapse = ", "), " log-transformed)"),
    ", cut into ", x$settings$n_bins, " x ", x$settings$n_bins, ", that ",
    "hold at least ", x$settings$min_records, " whole-range record(s); each ",
    "counts once. Gaps: small below ", pct(lim[1]), ", moderate to ",
    pct(lim[2]), ", large to ", pct(lim[3]), ", very large above (rules of ",
    "thumb; only the last is from Goicolea et al. 2025). Records in ",
    "unfilled cells: unfilled_records()."), prefix = "  ")
  invisible(x)
}

#' Records of the whole range in the unfilled part of the niche
#'
#' The whole-range records that lie in cells of the environmental space
#' without regional records, with their coordinates, so that they can be
#' mapped and checked (e.g. taxonomy, coordinates, other subspecies).
#'
#' @param x Result of [niche_coverage()].
#' @param all Return all records with their category instead.
#' @return A data.frame with `id` (row of the records used), `lon`, `lat`,
#'   `pc1`, `pc2`, `in_region`, `niche_cell` and `category`
#'   (`"elsewhere, available"`: conditions that occur in the region;
#'   `"elsewhere, absent"`: conditions the region lacks).
#' @export
#' @examples
#' # see niche_coverage()
unfilled_records <- function(x, all = FALSE) {
  if (!inherits(x, "niche_coverage")) {
    stop("`x` must come from niche_coverage().", call. = FALSE)
  }
  r <- x$records
  if (all) return(r)
  out <- r[r$category %in% c("elsewhere, available", "elsewhere, absent"), ,
           drop = FALSE]
  out$category <- droplevels(out$category)
  rownames(out) <- NULL
  out
}

#' Place new points in the niche space of a niche_coverage result
#'
#' Classifies points (e.g. fossil or evaluation sites, sites of a future
#' projection, or cells of the region under another climate) in the
#' environmental space and grid fitted by [niche_coverage()], and says, for
#' each, whether its conditions are known from regional records. This is the
#' decisive test of whether a nested model can help: a regional model
#' misses conditions that are known only from elsewhere; where its errors
#' fall into cells already filled by regional records, extending the model
#' to the whole range will not fix them.
#'
#' Categories of the cell a point falls into:
#'
#' * `filled`: niche cell with regional records.
#' * `elsewhere, available`: niche cell known only from records outside
#'   the region; the conditions occur in the region.
#' * `elsewhere, absent`: niche cell known only from elsewhere; the region
#'   lacks these conditions today.
#' * `region only`: cell with regional records but fewer than `min_records`
#'   whole-range records.
#' * `available, unoccupied`: no records, but conditions of the region.
#' * `unoccupied`: inside the space of the whole-range records, but no
#'   records and not conditions of the region.
#' * `outside`: beyond the range of the scores of the whole-range records,
#'   i.e. conditions the species is not known from at all.
#'
#' **Units.** The new values must be in the units of the layers used in
#' [niche_coverage()]. Points whose scores lie far beyond the scores of the
#' whole-range records (more than `far` times their range from its centre)
#' and variables whose median lies far outside the range of the
#' whole-range records give a warning; a common cause is a variable in
#' Kelvin in one data set and in degrees Celsius in the other.
#'
#' @param x Result of [niche_coverage()].
#' @param values A data.frame with one column per variable of `x` (the
#'   layer names), e.g. values at fossil sites under the climate of their
#'   time. Or give `points` and `env`.
#' @param points,env Points (data.frame, CSV path, `sf`, `SpatVector`) and
#'   layers (same variable names) to extract the values from.
#' @param predicted Optional logical vector: whether the model predicted each
#'   point as suitable. Gives a table of misses per category.
#' @param far Points farther than `far` times the half range of the
#'   whole-range scores from their centre give a units warning.
#' @param coords,points_crs Coordinate columns and CRS of `points`.
#' @return A data.frame of class `niche_classes`: `pc1`, `pc2`, `cell`,
#'   `category`, `regional_records` (in the point's cell),
#'   `regional_density` (percentile of that number among the cells with
#'   regional records) and, if given, `predicted`. With `predicted`, the
#'   table of misses is in `attr(x, "misses")`.
#' @seealso [niche_coverage()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 200)
#' west <- terra::as.polygons(terra::ext(14, 17, 49, 55), crs = "EPSG:4326")
#' nc <- niche_coverage(sim$presence, sim$env, west, n_boot = 0, n_null = 0)
#' # the region under a warmer climate
#' warm <- sim$env + 0.2
#' names(warm) <- names(sim$env)
#' cells <- terra::spatSample(terra::crop(warm, west), 200, na.rm = TRUE)
#' cl <- classify_points(nc, values = cells,
#'                       predicted = runif(nrow(cells)) > 0.3)
#' cl
classify_points <- function(x, values = NULL, points = NULL, env = NULL,
                            predicted = NULL, far = 3, coords = NULL,
                            points_crs = "EPSG:4326") {
  if (!inherits(x, "niche_coverage") || is.null(x$fit)) {
    stop("`x` must come from niche_coverage() (version 2.21.0 or later).",
         call. = FALSE)
  }
  f <- x$fit
  vars <- f$vars
  if (is.null(values)) {
    if (is.null(points) || is.null(env)) {
      stop("Give `values`, or `points` and `env`.", call. = FALSE)
    }
    r <- as_grid(env)
    miss <- setdiff(vars, names(r))
    if (length(miss)) stop("`env` lacks layer(s): ",
                           paste(miss, collapse = ", "), call. = FALSE)
    xy <- points_xy(points, r, coords, points_crs)
    values <- terra::extract(r[[vars]], xy)
  }
  values <- as.data.frame(values)
  miss <- setdiff(vars, names(values))
  if (length(miss)) {
    stop("`values` lacks column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  V <- as.matrix(values[vars])
  ok <- stats::complete.cases(V)
  # units: medians far outside the range of the whole-range records
  med <- apply(V[ok, , drop = FALSE], 2, stats::median)
  span <- f$q[2, ] - f$q[1, ]
  off <- (med < f$q[1, ] - 3 * span) | (med > f$q[2, ] + 3 * span)
  if (any(off, na.rm = TRUE)) {
    warning("The median of ", paste(vars[which(off)], collapse = ", "),
            " (", paste(signif(med[which(off)], 4), collapse = ", "),
            ") is far outside the values of the whole-range records (",
            paste(signif(f$medians[which(off)], 4), collapse = ", "),
            "): different units (e.g. Kelvin and degrees C)?", call. = FALSE)
  }
  S <- matrix(NA_real_, nrow(V), 2)
  if (any(ok)) {
    s <- stats::predict(f$pca, f$tf(V[ok, , drop = FALSE]))[, seq_len(f$k),
                                                             drop = FALSE]
    if (f$k == 1) s <- cbind(s, 0)
    S[ok, ] <- s
  }
  lim <- f$grid$lim
  ctr <- colMeans(lim)
  half <- (lim[2, ] - lim[1, ]) / 2
  far_pts <- ok & (abs(S[, 1] - ctr[1]) > far * half[1] |
                     abs(S[, 2] - ctr[2]) > far * max(half[2], 1e-12))
  if (any(far_pts)) {
    warning(sum(far_pts), " of ", sum(ok), " points lie more than ", far,
            " times the range of the whole-range records away in the ",
            "niche space: check that their values are in the same units as ",
            "the layers of niche_coverage().", call. = FALSE)
  }
  id <- grid_cells(S, f$grid)
  cat_pt <- coverage_category(id, f$n_all, f$n_reg, f$avail, f$min_records)
  cat_pt[!ok] <- NA
  dens <- coverage_density(f$n_reg)
  out <- data.frame(pc1 = S[, 1], pc2 = S[, 2], cell = id,
                    category = factor(cat_pt, coverage_levels),
                    regional_records = ifelse(is.na(id), 0L, f$n_reg[id]),
                    regional_density = ifelse(is.na(id), 0, dens[id]))
  out$regional_records[!ok] <- NA
  out$regional_density[!ok] <- NA
  if (!is.null(predicted)) {
    if (length(predicted) != nrow(out)) {
      stop("`predicted` must have one value per point.", call. = FALSE)
    }
    out$predicted <- as.logical(predicted)
    keep <- !is.na(out$category) & !is.na(out$predicted)
    tab <- table(category = droplevels(out$category[keep]),
                 predicted = factor(ifelse(out$predicted[keep], "predicted",
                                           "missed"),
                                    c("missed", "predicted")))
    misses <- data.frame(category = rownames(tab),
                         missed = as.vector(tab[, "missed"]),
                         predicted = as.vector(tab[, "predicted"]),
                         stringsAsFactors = FALSE)
    misses$share_of_misses <- misses$missed / max(sum(misses$missed), 1)
    misses$miss_rate <- misses$missed / pmax(misses$missed +
                                                misses$predicted, 1)
    attr(out, "misses") <- misses
  }
  class(out) <- c("niche_classes", "data.frame")
  out
}

#' @rdname classify_points
#' @param object Result of [niche_coverage()].
#' @param newdata Values (a data.frame with the variables) or points with
#'   `env`.
#' @param ... Passed to [classify_points()].
#' @export
predict.niche_coverage <- function(object, newdata, env = NULL, ...) {
  if (is.null(env)) classify_points(object, values = newdata, ...) else
    classify_points(object, points = newdata, env = env, ...)
}

#' @export
print.niche_classes <- function(x, ...) {
  pct <- function(z) paste0(round(100 * z), "%")
  d <- as.data.frame(unclass(x))
  cat("Points in the niche space:", sum(!is.na(d$category)), "\n\n")
  t <- table(d$category)
  t <- t[t > 0]
  print(data.frame(category = names(t), points = as.vector(t),
                   share = pct(as.vector(t) / sum(t))), row.names = FALSE)
  m <- attr(x, "misses")
  if (!is.null(m) && sum(m$missed) > 0) {
    cat("\nMisses (points the model did not predict) per category:\n")
    show <- m
    show$share_of_misses <- pct(show$share_of_misses)
    show$miss_rate <- pct(show$miss_rate)
    print(show, row.names = FALSE)
    filled <- sum(m$missed[m$category == "filled"]) / sum(m$missed)
    elsewhere <- sum(m$missed[grepl("^elsewhere", m$category)]) /
      sum(m$missed)
    cat("\nIN SHORT\n")
    txt <- c(
      paste0(pct(filled), " of the misses lie in cells filled by regional ",
             "records: the regional model knows these conditions, so a ",
             "nested model would not help there (look at the model, the ",
             "variables or the records instead)."),
      paste0(pct(elsewhere), " lie in cells known only from records ",
             "outside the region: a nested model could help these."))
    for (s in txt) {
      cat(strwrap(s, width = 78, exdent = 2, prefix = "  ", initial = "- "),
          sep = "\n")
    }
  }
  invisible(x)
}

#' Plot the niche coverage of a region
#'
#' Two panels. **Niche space**: the space of the first two principal
#' components of the records of the whole range, with the conditions of the
#' region (green), the niche cells without regional records (red frames;
#' their share is `unfilling`), the whole-range records outside the region
#' (grey) and inside it (blue), and separate regional records if given
#' (`region_points`, red crosses); the legends give numbers and shares.
#' **Map**:
#' the whole-range records coloured by the category of their niche cell
#' (filled by regional records; known only from elsewhere, with conditions
#' available in the region or absent from it), with the region outlined.
#'
#' @param x Result of [niche_coverage()].
#' @param type `"both"` (default; the two panels side by side), `"pca"` or
#'   `"map"`.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches (default: 14 x 6
#'   for both panels, 7.5 x 5.5 for one).
#' @return For one panel a ggplot object; for both, the two ggplots in a
#'   list (invisibly), after drawing them.
#' @export
#' @examples
#' # see niche_coverage()
plot_niche_coverage <- function(x, type = c("both", "pca", "map"),
                                base_size = 11, file = NULL, width = NULL,
                                height = NULL) {
  if (!inherits(x, "niche_coverage")) {
    stop("`x` must come from niche_coverage().", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.null(width)) width <- if (type == "both") 14 else 7.5
  if (is.null(height)) height <- if (type == "both") 6 else 5.5
  if (type == "both" && is.null(x$records)) type <- "pca"
  p_pca <- if (type != "map") coverage_pca_plot(x, base_size)
  p_map <- if (type != "pca") coverage_map_plot(x, base_size)
  p <- switch(type, pca = p_pca, map = p_map, list(pca = p_pca, map = p_map))
  if (!is.null(file)) save_figure(p, file, width, height)
  if (type == "both") {
    print_side_by_side(p)
    return(invisible(p))
  }
  p
}

coverage_pca_plot <- function(x, base_size) {
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  s <- x$scores
  n_in <- sum(s$region)
  n_out <- sum(!s$region)
  n_tot <- n_in + n_out
  lab_in <- sprintf("whole-range records inside the region (%d, %s)",
                    n_in, share_label(n_in / n_tot))
  lab_out <- sprintf("whole-range records outside the region (%d, %s)",
                     n_out, share_label(n_out / n_tot))
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
  bw <- if (!is.null(x$fit)) x$fit$grid$width else x$bin_width
  if (!is.null(b)) {
    if (!"category" %in% names(b)) {
      b$category <- ifelse(b$regional, "filled", "elsewhere, available")
    }
    b <- b[b$category %in% c("elsewhere, available", "elsewhere, absent"), ,
           drop = FALSE]
    lab_gap <- sprintf(paste0("niche cells without regional records ",
                              "(%d of %d, %s)"), nrow(b),
                       x$summary$niche_cells,
                       share_label(nrow(b) / x$summary$niche_cells))
    b$what <- factor(rep(lab_gap, nrow(b)), levels = lab_gap)
  }
  ve <- summary(x$pca)$importance[2, ]
  ax <- function(j) sprintf("PC%d (%.0f%% of the variance)", j, 100 * ve[j])
  ggplot2::ggplot(s, ggplot2::aes(.data$pc1, .data$pc2)) +
    ggplot2::geom_tile(data = rt, inherit.aes = FALSE,
                       ggplot2::aes(.data$pc1, .data$pc2), width = rt_w[1],
                       height = rt_w[2], fill = "#c7e9c0", colour = NA) +
    (if (!is.null(b) && nrow(b)) {
      ggplot2::geom_tile(data = b, inherit.aes = FALSE,
                         ggplot2::aes(.data$pc1, .data$pc2, fill = .data$what),
                         width = bw[1], height = bw[2],
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
      title = "The niche of the whole range and the region",
      subtitle = sprintf(paste0("Green: conditions that occur in the ",
                                "region; red frames: niche cells without ",
                                "regional records.\n%.0f%% of the niche has ",
                                "no record in the region; %.0f%% lies in ",
                                "conditions absent from the region."),
                         100 * x$summary$unfilling,
                         100 * x$summary$not_available)) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical")
}

coverage_map_plot <- function(x, base_size) {
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  d <- x$records
  lev <- c("filled", "elsewhere, available", "elsewhere, absent",
           "not in the niche")
  cat4 <- as.character(d$category)
  cat4[!cat4 %in% lev[1:3]] <- lev[4]
  n <- table(factor(cat4, lev))
  labs <- c(
    filled = "niche cell filled by regional records",
    `elsewhere, available` = "known only from elsewhere; conditions in the region",
    `elsewhere, absent` = "known only from elsewhere; conditions not in the region",
    `not in the niche` = "cell with too few records")
  labs <- sprintf("%s (%d, %s)", labs, as.vector(n),
                  share_label(as.vector(n) / max(sum(n), 1)))
  names(labs) <- lev
  d$cat <- factor(cat4, lev, labels = labs)
  d <- d[order(d$cat), , drop = FALSE]
  cols <- stats::setNames(c("#1c5cab", "#e0912f", "#c0392b", "#b4b3ae"),
                          labs)
  p <- ggplot2::ggplot()
  land <- NULL
  if (requireNamespace("rnaturalearth", quietly = TRUE)) {
    land <- tryCatch(rnaturalearth::ne_countries(scale = 110,
                                                 returnclass = "sf"),
                     error = function(e) NULL)
  }
  if (!is.null(land)) {
    p <- p + ggplot2::geom_sf(data = land, fill = "#e3e3e3", colour = "white",
                              linewidth = 0.2)
  }
  if (!is.null(x$region_outline)) {
    p <- p + ggplot2::geom_sf(data = x$region_outline, fill = "#c7e9c0",
                              colour = "#41ab5d", linewidth = 0.4,
                              alpha = 0.8)
  }
  xr <- range(d$lon)
  yr <- range(d$lat)
  pad <- function(r) r + c(-1, 1) * max(diff(r) * 0.05, 1)
  p + ggplot2::geom_point(data = d, ggplot2::aes(.data$lon, .data$lat,
                                                 colour = .data$cat),
                          size = 0.9, alpha = 0.85) +
    (if (!is.null(x$region_outline)) {
      ggplot2::geom_sf(data = x$region_outline, fill = NA,
                       colour = "#238b45", linewidth = 0.5)
    }) +
    ggplot2::scale_colour_manual(values = cols, name = NULL, drop = TRUE) +
    ggplot2::guides(colour = ggplot2::guide_legend(
      ncol = 1, override.aes = list(size = 2.5))) +
    ggplot2::coord_sf(xlim = pad(xr), ylim = pad(yr), expand = FALSE,
                      crs = "EPSG:4326", default_crs = "EPSG:4326") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Whole-range records by niche cell",
                  subtitle = paste0("Green: the region.\nOrange, red: ",
                                    "records in niche cells without ",
                                    "regional records.")) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom")
}

# A share as a percentage label: one decimal below 1% (so that a few
# records do not show as 0%), none above (internal)
share_label <- function(p) {
  ifelse(p > 0 & p < 0.01, sprintf("%.1f%%", 100 * p),
         sprintf("%.0f%%", 100 * p))
}
