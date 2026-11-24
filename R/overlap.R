# Niche overlap in environmental space ----------------------------------------------

#' Niche overlap between groups of records
#'
#' Compares the environmental niches of two or more groups of records:
#' species, lineages, cytotypes, or native and introduced ranges (a
#' grouping column handles within-species comparisons, e.g. Zhang and
#' Kubota 2021). Follows the framework of Broennimann et al. (2012), as in
#' ecospat and ENMTools (Warren et al. 2021):
#'
#' 1. **Environmental space.** The predictor values of the backgrounds of
#'    all groups are pooled and summarised by a principal component
#'    analysis (the PCA of [thin_environment()], precipitation amounts
#'    log-transformed); the first two axes form the space (`space =
#'    "pca"`). With `space = "raw"`, two standardised predictors are used
#'    instead (`vars`).
#' 2. **Densities.** The space is cut into a `grid` x `grid` lattice; the
#'    densities of records and of available environments (the background)
#'    of each group are estimated with a Gaussian kernel (reference
#'    bandwidth, as in ecospat). The occupancy is the density of records
#'    divided by the density of the environment (`correct = TRUE`,
#'    Broennimann et al. 2012), so that common environments do not look
#'    preferred.
#' 3. **Overlap.** Schoener's `D = 1 - 0.5 * sum(|p1 - p2|)` and Warren's
#'    `I = 1 - 0.5 * sum((sqrt(p1) - sqrt(p2))^2)`, from the two
#'    occupancies scaled to sum to 1: 0 = no overlap, 1 = identical. This is
#'    the corrected Hellinger index `I_cor` of Rodder and Engler (2011),
#'    used by ecospat and ENMTools; the original form `1 - H / 2` (Warren et
#'    al. 2008) cannot fall below 0.29.
#' 4. **Breadth.** Levins' B2 (standardised inverse Simpson index of the
#'    occupancy over the cells with available environment): 0 = one cell,
#'    1 = all cells equally.
#' 5. **Niche dynamics** of the second group relative to the first (Guisan
#'    et al. 2014): `expansion` (share of the second group's occupancy in
#'    environments the first does not occupy), `stability` (the rest) and
#'    `unfilling` (share of the first group's occupancy in environments the
#'    second does not occupy), within the environments available to both.
#'    `marginal_env` removes rare environments (e.g. 0.25 for the 75th
#'    percentile used by Lo Parrino et al. 2023), `marginal_records` the
#'    most marginal records (e.g. 0.05 for 95% of the records).
#'
#' **Tests** (`n_rep` replicates each), reported with a null distribution,
#' p-values and a standardised effect size (`(observed - null mean) / null
#' SD`), as in [null_model_test()] and [transfer_bias_test()]:
#' * **Equivalency** (identity) test: the records of both groups are pooled
#'   and randomly re-split into groups of the same sizes (Warren et al.
#'   2008). `p_lower` small = the niches are less alike than random splits:
#'   not equivalent.
#' * **Similarity** (background) test, in both directions: the occupancy of
#'   one group is compared with that of random points from the background
#'   of the other (Warren et al. 2008, ENMTools). `p_higher` small = more
#'   similar than expected from the environments available (niche
#'   conservatism); `p_lower` small = less similar (divergence).
#'
#' **Geographic overlap.** ENMTools first measured overlap between model
#' predictions in geographic space. With `predictions` (one suitability
#' raster per group, layer names = groups), D, I and the Spearman
#' correlation of the suitabilities are added ([prediction_overlap()]).
#' They describe where the groups could co-occur in the area mapped, not
#' their niches: they depend on the extent and on which environments it
#' contains, so niche comparisons are better made in environmental space
#' (Broennimann et al. 2012).
#'
#' **Sensitivity.** Whether a niche shift is detected can depend on the
#' background, the variables, the space and the thinning more than on the
#' species (Lo Parrino et al. 2023); [overlap_scenarios()] repeats the
#' analysis across such choices.
#'
#' @param records Records: a data.frame (or CSV path, `sf`, SpatVector)
#'   with coordinates and a `group` column, or a named list of record
#'   tables, one per group.
#' @param env Predictor rasters shared by all groups.
#' @param group Column with the group of each record.
#' @param background Available environment of each group: `NULL` (the
#'   whole raster, the same for all groups), one area for all groups, a
#'   named list of areas per group (polygons, `sf`, SpatVector, an
#'   `m_area` from [delimit_m()], or a raster whose non-`NA`, non-0 cells
#'   are the area), or a data.frame of background points with coordinates
#'   and the `group` column.
#' @param space `"pca"` or `"raw"`.
#' @param vars For `space = "raw"`, the two predictors (default: the first
#'   two layers).
#' @param log_vars Predictors to log-transform before the PCA (default:
#'   Bioclim precipitation amounts, as in [thin_environment()]).
#' @param grid Number of grid cells per axis of the density estimate.
#' @param correct Divide the density of records by the density of the
#'   environment (Broennimann et al. 2012).
#' @param marginal_records Share of records at the lowest densities
#'   left out of the occupied envelope used for expansion, stability,
#'   unfilling and the plot: cells below the density at that quantile of
#'   the records are not occupied. Default 0.05: the niche holds 95% of the
#'   records (Lo Parrino et al. 2023); with 0 the kernel tails make nearly
#'   the whole space "occupied". D, I and breadth use the full densities,
#'   as ecospat does by default.
#' @param truncate Also compute D, I and breadth from densities set to 0
#'   outside the occupied envelope (`marginal_records`), as ecospat does
#'   with `th.sp > 0`. Default `FALSE`: in simulations with known overlap
#'   (uniform niches in a uniform environment; 100 replicates of 20, 50
#'   and 200 records per group), truncation lowered D in every case, e.g.
#'   for a true D of 0.5 a mean of 0.40 instead of 0.48 with 200 records
#'   and 0.32 instead of 0.46 with 20, and raised its error; it only
#'   brought D of fully separate niches from 0.007 to exactly 0. The
#'   equivalency and similarity tests use the same choice.
#' @param marginal_env Share of the rarest available environments ignored
#'   in the niche dynamics (0 = none; 0.25 = above the 25th percentile).
#' @param n_rep Number of replicates of each test (0 = no tests).
#' @param tests `"equivalency"` and/or `"similarity"`.
#' @param n_background Number of background points sampled per group.
#' @param predictions Optional suitability rasters, one layer per group
#'   (named as the groups), for overlap in geographic space.
#' @param seed Optional random seed.
#' @inheritParams extract_predictions
#' @return A list of class `niche_overlap`: `pairs` (one row per pair of
#'   groups: `D`, `I`, `expansion`, `stability`, `unfilling` and, with
#'   tests, `p`, `ses` and verdicts, and with `predictions`, `geo_D`,
#'   `geo_I`, `geo_rho`), `groups` (records, background points and breadth
#'   per group), `null` (the null distributions), `densities` (the grids),
#'   `space` (axes, variance explained, loadings) and `settings`.
#' @references Broennimann, O., Fitzpatrick, M.C., Pearman, P.B.,
#'   Petitpierre, B., Pellissier, L., Yoccoz, N.G., Thuiller, W., Fortin,
#'   M.-J., Randin, C., Zimmermann, N.E., Graham, C.H., Guisan, A. (2012)
#'   Measuring ecological niche overlap from occurrence and spatial
#'   environmental data. *Global Ecology and Biogeography* 21, 481-497.
#'   \doi{10.1111/j.1466-8238.2011.00698.x}
#'
#'   Warren, D.L., Glor, R.E., Turelli, M. (2008) Environmental niche
#'   equivalency versus conservatism: quantitative approaches to niche
#'   evolution. *Evolution* 62, 2868-2883.
#'   \doi{10.1111/j.1558-5646.2008.00482.x}
#'
#'   Warren, D.L. et al. (2021) ENMTools 1.0: an R package for comparative
#'   ecological biogeography. *Ecography* 44, 504-511.
#'   \doi{10.1111/ecog.05485}
#'
#'   Guisan, A., Petitpierre, B., Broennimann, O., Daehler, C., Kueffer, C.
#'   (2014) Unifying niche shift studies: insights from biological
#'   invasions. *Trends in Ecology & Evolution* 29, 260-269.
#'   \doi{10.1016/j.tree.2014.02.009}
#'
#'   Lo Parrino, E., Falaschi, M., Manenti, R., Ficetola, G.F. (2023) All
#'   that changes is not shift: methodological choices influence niche shift
#'   detection in freshwater invasive species. *Ecography* 2023, e06432.
#'   \doi{10.1111/ecog.06432}
#'
#'   Zhang, Z., Kubota, S. (2021) Integrating intraspecific variation in
#'   species distribution models by quantifying niche differentiation.
#'   *Biological Journal of the Linnean Society* 133, 187-201.
#'   \doi{10.1093/biolinnean/blab021}
#'
#'   Rodder, D., Engler, J.O. (2011) Quantitative metrics of overlaps in
#'   Grinnellian niches: advances and possible drawbacks. *Global Ecology
#'   and Biogeography* 20, 915-927. \doi{10.1111/j.1466-8238.2011.00659.x}
#' @seealso [overlap_scenarios()], [plot_niche_overlap()],
#'   [prediction_overlap()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pres <- sim$presence
#' # two "lineages": west and east of the area
#' pres$group <- ifelse(pres$lon < 19, "west", "east")
#' no <- niche_overlap(pres, sim$env, n_rep = 19, seed = 1)
#' no
#' plot_niche_overlap(no)
niche_overlap <- function(records, env, group = "group", background = NULL,
                          space = c("pca", "raw"), vars = NULL,
                          log_vars = NULL, grid = 100, correct = TRUE,
                          marginal_records = 0.05, truncate = FALSE,
                          marginal_env = 0, n_rep = 100, tests = c("equivalency", "similarity"),
                          n_background = 10000, predictions = NULL,
                          seed = NULL, coords = NULL, crs = "EPSG:4326") {
  space <- match.arg(space)
  tests <- match.arg(tests, c("equivalency", "similarity"),
                     several.ok = TRUE)
  r <- as_grid(env)
  occ <- overlap_records(records, group, r, coords, crs)
  groups <- unique(occ$group)
  if (length(groups) < 2) stop("At least two groups are needed.",
                               call. = FALSE)
  bg <- with_seed(seed, overlap_background(background, groups, r,
                                           n_background, coords, crs,
                                           group))
  # environmental space from the pooled backgrounds
  B <- do.call(rbind, bg)
  B <- B[stats::complete.cases(B), , drop = FALSE]
  if (space == "pca") {
    es <- env_pca(r, log_vars, axes = "pca", values = B, pc_min = 2)
    to_axes <- env_axes_fun(es$pca, colnames(B), es$log_vars, "pca", 2)
    ve <- es$variance
    sp_info <- list(axes = c("PC1", "PC2"), variance = ve[1:2],
                    loadings = es$pca$rotation[, 1:2, drop = FALSE])
  } else {
    if (is.null(vars)) vars <- colnames(B)[1:2]
    if (length(vars) != 2 || !all(vars %in% colnames(B))) {
      stop("`vars` must name two predictors.", call. = FALSE)
    }
    mu <- colMeans(B[, vars])
    sdv <- apply(B[, vars], 2, stats::sd)
    to_axes <- function(X) {
      A <- sweep(sweep(as.matrix(X[, vars, drop = FALSE]), 2, mu, "-"), 2,
                 sdv, "/")
      colnames(A) <- vars
      A
    }
    sp_info <- list(axes = vars, variance = NULL, loadings = NULL)
  }
  S_occ <- lapply(groups, function(g) {
    to_axes(occ$values[occ$group == g, , drop = FALSE])
  })
  S_bg <- lapply(bg, function(b) {
    to_axes(b[stats::complete.cases(b), , drop = FALSE])
  })
  names(S_occ) <- names(S_bg) <- groups
  ext <- range_ext(c(S_occ, S_bg))
  dens <- lapply(groups, function(g) {
    group_density(S_occ[[g]], S_bg[[g]], ext, grid, marginal_records,
                  truncate)
  })
  names(dens) <- groups

  pairs <- utils::combn(groups, 2, simplify = FALSE)
  null <- list()
  rows <- lapply(pairs, function(pr) {
    a <- dens[[pr[1]]]
    b <- dens[[pr[2]]]
    ov <- overlap_DI(a, b, correct)
    dy <- niche_dynamics(a, b, marginal_env)
    row <- data.frame(group1 = pr[1], group2 = pr[2], n1 = nrow(S_occ[[pr[1]]]),
                      n2 = nrow(S_occ[[pr[2]]]), D = ov[["D"]], I = ov[["I"]],
                      expansion = dy[["expansion"]],
                      stability = dy[["stability"]],
                      unfilling = dy[["unfilling"]],
                      stringsAsFactors = FALSE)
    key <- paste(pr, collapse = " vs ")
    if (n_rep > 0) {
      if ("equivalency" %in% tests) {
        nd <- with_seed(if (is.null(seed)) NULL else seed + 1,
                        equivalency_null(S_occ[[pr[1]]], S_occ[[pr[2]]], a, b,
                                         ext, grid, marginal_records, correct,
                                         n_rep, truncate))
        null[[key]]$equivalency <<- nd
        row <- cbind(row, test_columns(ov[["D"]], nd[, "D"], "eq"))
      }
      if ("similarity" %in% tests) {
        n12 <- with_seed(if (is.null(seed)) NULL else seed + 2,
                         similarity_null(a, S_bg[[pr[2]]], nrow(S_occ[[pr[2]]]),
                                         b, ext, grid, marginal_records,
                                         correct, n_rep, truncate = truncate))
        n21 <- with_seed(if (is.null(seed)) NULL else seed + 3,
                         similarity_null(b, S_bg[[pr[1]]], nrow(S_occ[[pr[1]]]),
                                         a, ext, grid, marginal_records,
                                         correct, n_rep, swap = TRUE,
                                         truncate = truncate))
        null[[key]]$similarity_12 <<- n12
        null[[key]]$similarity_21 <<- n21
        row <- cbind(row, test_columns(ov[["D"]], n12[, "D"], "sim12"),
                     test_columns(ov[["D"]], n21[, "D"], "sim21"))
      }
    }
    if (!is.null(predictions)) {
      po <- prediction_overlap(predictions[[pr]])
      row$geo_D <- po$D[1]
      row$geo_I <- po$I[1]
      row$geo_rho <- po$rho[1]
    }
    row
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL
  gr <- data.frame(group = groups,
                   n_records = vapply(S_occ, nrow, integer(1)),
                   n_background = vapply(S_bg, nrow, integer(1)),
                   breadth_B2 = vapply(dens, function(d) {
                     levins_b2(if (correct) d$z_cor else d$z_uncor, d$Z > 0)
                   }, numeric(1)), stringsAsFactors = FALSE)
  rownames(gr) <- NULL
  structure(list(pairs = tab, groups = gr, null = null, densities = dens,
                 scores = list(records = S_occ, background = S_bg),
                 space = c(sp_info, list(type = space, extent = ext)),
                 settings = list(grid = grid, correct = correct,
                                 marginal_records = marginal_records,
                                 truncate = truncate,
                                 marginal_env = marginal_env, n_rep = n_rep,
                                 tests = tests,
                                 background = if (is.null(background))
                                   "whole area" else "given")),
            class = "niche_overlap")
}

# Records as group + predictor values (internal)
overlap_records <- function(records, group, r, coords, crs) {
  if (is.list(records) && !is.data.frame(records) &&
      !inherits(records, c("sf", "SpatVector"))) {
    if (is.null(names(records)) || any(names(records) == "")) {
      stop("A list of record tables must be named by group.", call. = FALSE)
    }
    pts <- lapply(names(records), function(g) {
      p <- as_points(records[[g]], coords, crs)
      list(p = p, g = rep(g, nrow(p)))
    })
    p <- do.call(rbind, lapply(pts, `[[`, "p"))
    g <- unlist(lapply(pts, `[[`, "g"))
  } else {
    if (is.character(records) && length(records) == 1) {
      records <- read_any_csv(records)
    }
    df <- as.data.frame(records)
    if (!group %in% names(df)) {
      stop("No `", group, "` column in the records.", call. = FALSE)
    }
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
    x <- as_number(df[[xy_cols[1]]])
    y <- as_number(df[[xy_cols[2]]])
    ok <- !is.na(x) & !is.na(y) & !is.na(df[[group]])
    p <- terra::vect(cbind(x[ok], y[ok]), crs = crs)
    g <- as.character(df[[group]][ok])
  }
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  v <- as.matrix(terra::extract(r, p, ID = FALSE))
  keep <- stats::complete.cases(v)
  if (any(!keep)) {
    message(sum(!keep), " record(s) without predictor values were dropped.")
  }
  list(values = v[keep, , drop = FALSE], group = g[keep])
}

# Background predictor values per group (internal)
overlap_background <- function(background, groups, r, n, coords, crs,
                               group) {
  whole <- function() sample_complete(r, n)
  in_area <- function(a) {
    if (inherits(a, "m_area")) a <- a$m
    if (inherits(a, c("sf", "sfc"))) a <- terra::vect(a)
    if (is.character(a)) a <- terra::vect(a)
    if (inherits(a, "SpatVector")) {
      if (!terra::same.crs(a, r)) a <- terra::project(a, r)
      rr <- terra::mask(terra::crop(r, a, snap = "out"), a)
    } else if (inherits(a, "SpatRaster")) {
      m <- region_mask(a, r)
      rr <- terra::mask(r, m)
    } else {
      stop("Unknown background area.", call. = FALSE)
    }
    v <- sample_complete(rr, n)
    if (!nrow(v)) stop("A background area holds no cell with values.",
                       call. = FALSE)
    v
  }
  if (is.null(background)) {
    w <- whole()
    out <- rep(list(w), length(groups))
  } else if (is.data.frame(background) && group %in% names(background)) {
    xy_cols <- if (is.null(coords)) guess_coords(names(background)) else
      coords
    out <- lapply(groups, function(g) {
      b <- background[background[[group]] == g, , drop = FALSE]
      p <- terra::vect(cbind(as_number(b[[xy_cols[1]]]),
                             as_number(b[[xy_cols[2]]])), crs = crs)
      if (!terra::same.crs(p, r)) p <- terra::project(p, r)
      as.matrix(terra::extract(r, p, ID = FALSE))
    })
  } else if (is.list(background) && !inherits(background, c("sf", "sfc",
                                                             "m_area")) &&
             !is.null(names(background))) {
    miss <- setdiff(groups, names(background))
    if (length(miss)) {
      stop("No background for group(s): ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    out <- lapply(groups, function(g) in_area(background[[g]]))
  } else {
    a <- in_area(background)
    out <- rep(list(a), length(groups))
  }
  names(out) <- groups
  out
}

# Common extent of the space (internal)
range_ext <- function(lst) {
  A <- do.call(rbind, lst)
  c(range(A[, 1]), range(A[, 2]))
}

# Gaussian kernel density of points on a grid x grid lattice over `ext`,
# reference bandwidth on axes rescaled to 0-1, as adehabitatHR's "href"
# used by ecospat (internal). Returns a grid x grid matrix (rows = axis 1).
kernel_grid <- function(S, ext, grid) {
  x <- (S[, 1] - ext[1]) / max(ext[2] - ext[1], 1e-12)
  y <- (S[, 2] - ext[3]) / max(ext[4] - ext[3], 1e-12)
  n <- length(x)
  if (n < 2) {
    h <- 0.05
  } else {
    h <- sqrt(0.5 * (stats::var(x) + stats::var(y))) * n^(-1 / 6)
    if (!is.finite(h) || h <= 0) h <- 0.05
  }
  g <- seq(0, 1, length.out = grid)
  kx <- stats::dnorm(outer(x, g, "-"), sd = h)
  ky <- stats::dnorm(outer(y, g, "-"), sd = h)
  crossprod(kx, ky) / n
}

# Density of records at the records (nearest grid cell) (internal)
at_points <- function(m, S, ext, grid) {
  i <- 1 + round((S[, 1] - ext[1]) / max(ext[2] - ext[1], 1e-12) * (grid - 1))
  j <- 1 + round((S[, 2] - ext[3]) / max(ext[4] - ext[3], 1e-12) * (grid - 1))
  m[cbind(pmin(pmax(i, 1), grid), pmin(pmax(j, 1), grid))]
}

# Densities of one group (internal): Z environment, z records (masked by
# the environment), occupancy uncorrected and corrected, occupied cells
group_density <- function(So, Sb, ext, grid, marginal_records,
                          truncate = FALSE) {
  Zd <- kernel_grid(Sb, ext, grid)
  Zd[Zd < max(Zd) * 1e-6] <- 0
  Z <- Zd * nrow(Sb) / sum(Zd)
  out <- occupancy(kernel_grid(So, ext, grid), So, Z, ext, grid,
                   marginal_records, truncate)
  out$Z <- Z
  out
}

# Occupancy of records from their kernel density `zd` within the
# environment Z (internal): the occupied envelope `w` (cells above the
# density at the most marginal records) and, with `truncate`, the density
# set to 0 outside it before scaling (ecospat's th.sp)
occupancy <- function(zd, S, Z, ext, grid, marginal_records, truncate) {
  zd[Z <= 0] <- 0
  th <- if (marginal_records > 0) {
    stats::quantile(at_points(zd, S, ext, grid), marginal_records,
                    names = FALSE)
  } else max(zd) * 1e-3
  w <- zd >= th & Z > 0
  if (truncate && marginal_records > 0) zd[!w] <- 0
  z <- zd * nrow(S) / max(sum(zd), 1e-300)
  z_uncor <- z / max(z, 1e-300)
  zc <- ifelse(Z > 0, z / Z, 0)
  list(z = z, z_uncor = z_uncor, z_cor = zc / max(zc, 1e-300), w = w)
}

# Schoener's D and Warren's I (corrected form) (internal)
overlap_DI <- function(a, b, correct) {
  p1 <- if (correct) a$z_cor else a$z_uncor
  p2 <- if (correct) b$z_cor else b$z_uncor
  p1 <- p1 / sum(p1)
  p2 <- p2 / sum(p2)
  c(D = 1 - 0.5 * sum(abs(p1 - p2)),
    I = 1 - 0.5 * sum((sqrt(p1) - sqrt(p2))^2))
}

# Levins' B2 over the cells with available environment (internal)
levins_b2 <- function(z, avail) {
  p <- z[avail]
  if (length(p) < 2 || sum(p) <= 0) return(NA_real_)
  p <- p / sum(p)
  (1 / sum(p^2) - 1) / (length(p) - 1)
}

# Expansion, stability, unfilling of b relative to a, in environments
# available to both (ecospat.niche.dyn.index logic) (internal)
niche_dynamics <- function(a, b, marginal_env) {
  q1 <- if (marginal_env > 0) stats::quantile(a$Z[a$Z > 0], marginal_env,
                                              names = FALSE) else 0
  q2 <- if (marginal_env > 0) stats::quantile(b$Z[b$Z > 0], marginal_env,
                                              names = FALSE) else 0
  both <- (a$Z > q1) & (b$Z > q2)
  w1 <- a$w & both
  w2 <- b$w & both
  exp_c <- w2 & !w1
  stab <- w1 & w2
  unf <- w1 & !w2
  e <- sum(b$z_uncor[exp_c])
  s2 <- sum(b$z_uncor[stab])
  u <- sum(a$z_uncor[unf])
  s1 <- sum(a$z_uncor[stab])
  c(expansion = if (e + s2 > 0) e / (e + s2) else 0,
    stability = if (e + s2 > 0) s2 / (e + s2) else 0,
    unfilling = if (u + s1 > 0) u / (u + s1) else 0)
}

# Null distribution of the equivalency test (internal)
equivalency_null <- function(S1, S2, a, b, ext, grid, mr, correct, n_rep,
                             truncate = FALSE) {
  pool <- rbind(S1, S2)
  n1 <- nrow(S1)
  out <- t(vapply(seq_len(n_rep), function(i) {
    k <- sample.int(nrow(pool))
    da <- redensity(pool[k[seq_len(n1)], , drop = FALSE], a, ext, grid, mr,
                    truncate)
    db <- redensity(pool[k[-seq_len(n1)], , drop = FALSE], b, ext, grid, mr,
                    truncate)
    overlap_DI(da, db, correct)
  }, numeric(2)))
  colnames(out) <- c("D", "I")
  out
}

# Null distribution of the similarity test: `fixed` against random points
# from the background of the other group (internal)
similarity_null <- function(fixed, Sb_other, n_other, other, ext, grid, mr,
                            correct, n_rep, swap = FALSE, truncate = FALSE) {
  out <- t(vapply(seq_len(n_rep), function(i) {
    pts <- Sb_other[sample.int(nrow(Sb_other), n_other,
                               replace = nrow(Sb_other) < n_other), ,
                    drop = FALSE]
    dr <- redensity(pts, other, ext, grid, mr, truncate)
    if (swap) overlap_DI(dr, fixed, correct) else
      overlap_DI(fixed, dr, correct)
  }, numeric(2)))
  colnames(out) <- c("D", "I")
  out
}

# Occupancy of new points with the environment of an existing group
# (internal)
redensity <- function(S, g, ext, grid, mr, truncate = FALSE) {
  out <- occupancy(kernel_grid(S, ext, grid), S, g$Z, ext, grid, mr,
                   truncate)
  out$Z <- g$Z
  out
}

# Test summary columns (internal)
test_columns <- function(obs, null, prefix) {
  null <- null[!is.na(null)]
  n <- length(null)
  sdn <- stats::sd(null)
  out <- data.frame(
    p_lower = (1 + sum(null <= obs)) / (n + 1),
    p_higher = (1 + sum(null >= obs)) / (n + 1),
    ses = if (is.finite(sdn) && sdn > 0) (obs - mean(null)) / sdn else
      NA_real_,
    null_mean = mean(null))
  names(out) <- paste0(prefix, "_", names(out))
  out
}

#' Overlap of model predictions in geographic space
#'
#' Schoener's D, Warren's I and the Spearman rank correlation between
#' suitability rasters (ENMTools: Warren et al. 2008, 2021): each raster is
#' scaled to sum to 1 over the cells where all have values. D and I measure
#' how similar the amounts of suitability are cell by cell, i.e. where the
#' groups could meet in the mapped area; the correlation measures whether
#' they respond in the same direction. They depend on the extent mapped;
#' to compare niches, use [niche_overlap()] in environmental space.
#'
#' @param x Suitability rasters ([terra::SpatRaster] or files), two or more
#'   layers.
#' @param n_sample Cells sampled for the correlation.
#' @return A data.frame with one row per pair of layers: `layer1`,
#'   `layer2`, `D`, `I`, `rho`.
#' @references Warren, D.L., Glor, R.E., Turelli, M. (2008) Environmental
#'   niche equivalency versus conservatism: quantitative approaches to niche
#'   evolution. *Evolution* 62, 2868-2883.
#'   \doi{10.1111/j.1558-5646.2008.00482.x}
#'
#'   Warren, D.L. et al. (2021) ENMTools 1.0: an R package for comparative
#'   ecological biogeography. *Ecography* 44, 504-511.
#'   \doi{10.1111/ecog.05485}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' prediction_overlap(sim$models[[c("GLM", "GBM", "Random")]])
prediction_overlap <- function(x, n_sample = 1e5) {
  r <- as_grid(x)
  if (terra::nlyr(r) < 2) stop("Give at least two layers.", call. = FALSE)
  # cells with values in all layers
  m <- sum(r)
  r <- terra::mask(r, m)
  s <- terra::global(r, "sum", na.rm = TRUE)[, 1]
  if (any(s <= 0)) stop("A layer sums to 0.", call. = FALSE)
  smp <- sample_complete(r, n_sample)
  out <- do.call(rbind, utils::combn(seq_len(terra::nlyr(r)), 2,
                                     simplify = FALSE, FUN = function(k) {
    a <- r[[k[1]]] / s[k[1]]
    b <- r[[k[2]]] / s[k[2]]
    data.frame(layer1 = names(r)[k[1]], layer2 = names(r)[k[2]],
               D = 1 - 0.5 * terra::global(abs(a - b), "sum",
                                           na.rm = TRUE)[1, 1],
               I = 1 - 0.5 * terra::global((sqrt(a) - sqrt(b))^2, "sum",
                                           na.rm = TRUE)[1, 1],
               rho = suppressWarnings(stats::cor(smp[, k[1]], smp[, k[2]],
                                                 method = "spearman")),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

# Plain grade of an overlap (internal; the classes of Rodder and Engler
# 2011, p. 918: 0-0.2 no or very limited, 0.2-0.4 low, 0.4-0.6 moderate,
# 0.6-0.8 high, 0.8-1 very high)
overlap_grade <- function(d) {
  as.character(cut(d, c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
                   c("no or very limited", "low", "moderate", "high",
                     "very high")))
}

#' @export
print.niche_overlap <- function(x, digits = 2, alpha = 0.05, ...) {
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(v, digits = digits, format = "f")
  s <- x$settings
  cat("Niche overlap: ", nrow(x$groups), " groups in ",
      if (x$space$type == "pca") "the PCA environmental space" else
        paste0("the space of ", paste(x$space$axes, collapse = " and ")),
      " (", s$grid, " x ", s$grid, " grid; background: ", s$background,
      ")\n\nIN SHORT\n", sep = "")
  for (i in seq_len(nrow(x$pairs))) {
    p <- x$pairs[i, ]
    nm <- paste0(p$group1, " vs ", p$group2)
    say(paste0(nm, ": ", overlap_grade(p$D), " overlap (D = ", f(p$D),
               ", I = ", f(p$I), "). ", share_label(p$expansion), " of ",
               p$group2, "'s niche lies in environments that ", p$group1,
               " does not occupy, and ", share_label(p$unfilling), " of ",
               p$group1, "'s niche in environments that ", p$group2,
               " does not occupy."))
    if (!is.null(p$eq_p_lower)) {
      say(if (p$eq_p_lower < alpha) {
        paste0("Equivalency test: the niches are not equivalent ",
               "(p = ", f(p$eq_p_lower), "): the two groups are less alike ",
               "than random splits of their pooled records.")
      } else {
        paste0("Equivalency test: the niches cannot be told apart from ",
               "random splits of the pooled records (p = ",
               f(p$eq_p_lower), "); with few records the test has little ",
               "power.")
      })
    }
    if (!is.null(p$sim12_p_higher)) {
      sim_txt <- function(ph, pl, from, to) {
        if (ph < alpha) {
          paste0(from, " is more similar to ", to, " than expected from the ",
                 "environments available to ", to, " (p = ", f(ph), ")")
        } else if (pl < alpha) {
          paste0(from, " is less similar to ", to, " than expected from the ",
                 "environments available to ", to, " (p = ", f(pl), ")")
        } else {
          paste0(from, " is as similar to ", to, " as random points of its ",
                 "range (p = ", f(min(ph, pl)), ")")
        }
      }
      say(paste0("Similarity test: ", sim_txt(p$sim12_p_higher,
                                               p$sim12_p_lower, p$group1,
                                               p$group2), "; ",
                 sim_txt(p$sim21_p_higher, p$sim21_p_lower, p$group2,
                         p$group1), "."))
    }
    if (!is.null(p$geo_D)) {
      say(paste0("In geographic space (model predictions), D = ",
                 f(p$geo_D), ", I = ", f(p$geo_I), ", rho = ", f(p$geo_rho),
                 "; this depends on the area mapped, not only on the ",
                 "niches."))
    }
  }
  small <- x$groups$group[x$groups$n_records < 30]
  if (length(small)) {
    say(paste0("Fewer than 30 records in ", paste(small, collapse = ", "),
               ": densities and tests are uncertain."))
  }
  if (isTRUE(s$truncate)) {
    say(paste0("D, I and breadth are computed from densities cut at the ",
               "occupied envelope (truncate = TRUE). In simulations with ",
               "known overlap this lowered D, by up to a third with few ",
               "records; compare with truncate = FALSE before reading the ",
               "grade."))
  }
  if (s$n_rep > 0 && s$n_rep < 99) {
    say(paste0("Only ", s$n_rep, " replicates: p-values cannot be smaller ",
               "than ", f(1 / (s$n_rep + 1)), "; use 100 or more."))
  }
  cat("\n")
  show <- x$pairs[intersect(c("group1", "group2", "D", "I", "expansion",
                              "stability", "unfilling", "eq_p_lower",
                              "eq_ses", "sim12_p_higher", "sim21_p_higher",
                              "geo_D", "geo_I", "geo_rho"),
                            names(x$pairs))]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits + 1)
  print(show, row.names = FALSE)
  g <- x$groups
  g$breadth_B2 <- round(g$breadth_B2, digits + 1)
  cat("\nGroups:\n")
  print(g, row.names = FALSE)
  cat("", strwrap(paste0(
    "D = Schoener's overlap, I = Warren's (corrected) overlap of the ",
    if (s$correct) "occupancies (record density / environment density)"
    else "record densities", if (isTRUE(s$truncate))
      ", set to 0 outside the occupied envelope" else "", "; 0 = none, 1 = identical. Grades (below 0.2 ",
    "no or very limited, 0.2-0.4 low, 0.4-0.6 moderate, 0.6-0.8 high, above 0.8 ",
    "very high) are those suggested by Rodder and Engler (2011). breadth_B2 = Levins' ",
    "standardised breadth (0-1). Tests: equivalency p_lower small = not ",
    "equivalent; similarity p_higher small = more similar than expected ",
    "(conservatism), p_lower small = less similar. ses = (observed - null ",
    "mean) / null SD."), width = 78), sep = "\n")
  invisible(x)
}

#' Sensitivity of niche overlap to methodological choices
#'
#' Whether a niche shift is detected can depend more on methodological
#' choices than on the organisms: in freshwater invaders, background
#' selection, the variables and the environmental space changed overlap
#' and test results more than biological factors did (Lo Parrino et al.
#' 2023). `overlap_scenarios()` runs [niche_overlap()] over a grid of
#' choices, as [m_scenarios()] does for the accessible area, and says
#' whether the conclusion holds across them:
#'
#' * `backgrounds`: named alternatives for the available environment (e.g.
#'   the accessible areas of each group against the whole raster);
#' * `spaces`: the PCA environmental space or two raw predictors;
#' * `thin_km`: records unthinned (0) or thinned to one per cell of a grid
#'   of that size (km), which reduces the weight of intensively sampled
#'   places;
#' * `grids`: the resolution of the density estimate.
#'
#' For each pair of groups and each test, the printout reports in how many
#' scenarios the verdict is the same and, when it changes, which choice
#' changes it: the share of pairs of scenarios that differ only in that
#' choice and give different verdicts.
#'
#' @inheritParams niche_overlap
#' @param backgrounds Named list of `background` specifications for
#'   [niche_overlap()]; use `NULL` as an element for the whole raster, e.g.
#'   `list(M = list(A = m_A, B = m_B), whole = NULL)`.
#' @param spaces `"pca"` and/or `"raw"`.
#' @param thin_km Thinning distances (km); 0 = no thinning.
#' @param grids Grid sizes of the density estimate.
#' @param alpha Significance level for the verdicts.
#' @param ... Passed to [niche_overlap()] (e.g. `n_rep`, `vars`,
#'   `marginal_env`).
#' @return A list of class `overlap_scenarios`: `results` (one row per
#'   scenario and pair: the choices, `D`, `I`, `expansion`, `unfilling`,
#'   p-values and verdicts), `robustness` (per pair and test: number of
#'   scenarios per verdict and the choice that changes it most) and
#'   `settings`.
#' @references Lo Parrino, E., Falaschi, M., Manenti, R., Ficetola, G.F.
#'   (2023) All that changes is not shift: methodological choices influence
#'   niche shift detection in freshwater invasive species. *Ecography* 2023,
#'   e06432. \doi{10.1111/ecog.06432}
#'
#'   Broennimann, O. et al. (2012) Measuring ecological niche overlap from
#'   occurrence and spatial environmental data. *Global Ecology and
#'   Biogeography* 21, 481-497. \doi{10.1111/j.1466-8238.2011.00698.x}
#' @seealso [niche_overlap()], [m_scenarios()]
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_enm_data(seed = 1)
#' pres <- sim$presence
#' pres$group <- ifelse(pres$lon < 19, "west", "east")
#' os <- overlap_scenarios(pres, sim$env, thin_km = c(0, 30),
#'                         grids = c(50, 100), n_rep = 49, seed = 1)
#' os
#' }
overlap_scenarios <- function(records, env, group = "group",
                              backgrounds = list(whole = NULL),
                              spaces = c("pca", "raw"), thin_km = c(0, 10),
                              grids = c(50, 100), alpha = 0.05, seed = NULL,
                              coords = NULL, crs = "EPSG:4326", ...) {
  if (!is.list(backgrounds) || is.null(names(backgrounds)) ||
      any(names(backgrounds) == "")) {
    stop("`backgrounds` must be a named list (NULL = the whole raster).",
         call. = FALSE)
  }
  r <- as_grid(env)
  rec <- overlap_record_table(records, group, coords, crs)
  sc <- expand.grid(background = names(backgrounds), space = spaces,
                    thin_km = thin_km, grid = grids,
                    stringsAsFactors = FALSE)
  res <- lapply(seq_len(nrow(sc)), function(i) {
    s <- sc[i, ]
    d <- if (s$thin_km > 0) thin_grid_km(rec, s$thin_km) else rec
    no <- niche_overlap(d, r, group = "group",
                        background = backgrounds[[s$background]],
                        space = s$space, grid = s$grid, seed = seed,
                        coords = c("x", "y"), crs = crs, ...)
    p <- no$pairs
    out <- cbind(s[rep(1, nrow(p)), , drop = FALSE],
                 p[c("group1", "group2", "n1", "n2", "D", "I", "expansion",
                     "unfilling")])
    if ("eq_p_lower" %in% names(p)) {
      out$eq_p <- p$eq_p_lower
      out$equivalency <- ifelse(p$eq_p_lower < alpha, "not equivalent",
                                "equivalent")
    }
    if ("sim12_p_higher" %in% names(p)) {
      v <- function(ph, pl) ifelse(ph < alpha, "more similar",
                                   ifelse(pl < alpha, "less similar",
                                          "as expected"))
      out$similarity_12 <- v(p$sim12_p_higher, p$sim12_p_lower)
      out$similarity_21 <- v(p$sim21_p_higher, p$sim21_p_lower)
    }
    out
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  factors <- c("background", "space", "thin_km", "grid")
  varying <- factors[vapply(factors, function(k) {
    length(unique(res[[k]])) > 1
  }, logical(1))]
  tests <- intersect(c("equivalency", "similarity_12", "similarity_21"),
                     names(res))
  rob <- do.call(rbind, lapply(split(res, paste(res$group1, res$group2)),
                               function(d) {
    do.call(rbind, lapply(tests, function(tt) {
      tb <- table(d[[tt]])
      driver <- NA_character_
      flip <- NA_real_
      if (length(tb) > 1 && length(varying)) {
        fl <- vapply(varying, function(k) {
          others <- setdiff(varying, k)
          key <- if (length(others)) {
            do.call(paste, c(d[others], sep = "|"))
          } else rep("all", nrow(d))
          ch <- tapply(d[[tt]], key, function(v) length(unique(v)) > 1)
          mean(ch)
        }, numeric(1))
        driver <- names(fl)[which.max(fl)]
        flip <- max(fl)
      }
      data.frame(group1 = d$group1[1], group2 = d$group2[1], test = tt,
                 scenarios = nrow(d),
                 verdicts = paste(paste0(names(tb), ": ", as.vector(tb)),
                                  collapse = "; "),
                 robust = length(tb) == 1, driver = driver,
                 driver_flips = flip, stringsAsFactors = FALSE)
    }))
  }))
  rownames(rob) <- NULL
  structure(list(results = res, robustness = rob,
                 settings = list(varying = varying, alpha = alpha)),
            class = "overlap_scenarios")
}

# Records as one table with x, y, group (internal)
overlap_record_table <- function(records, group, coords, crs) {
  if (is.list(records) && !is.data.frame(records) &&
      !inherits(records, c("sf", "SpatVector"))) {
    if (is.null(names(records))) {
      stop("A list of record tables must be named by group.", call. = FALSE)
    }
    return(do.call(rbind, lapply(names(records), function(g) {
      p <- as_points(records[[g]], coords, crs)
      xy <- terra::crds(p)
      data.frame(x = xy[, 1], y = xy[, 2], group = g,
                 stringsAsFactors = FALSE)
    })))
  }
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  df <- as.data.frame(records)
  xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
  data.frame(x = as_number(df[[xy_cols[1]]]), y = as_number(df[[xy_cols[2]]]),
             group = as.character(df[[group]]), stringsAsFactors = FALSE)
}

# One record per cell of a grid of `km` per group, in an equal-area
# projection (internal)
thin_grid_km <- function(rec, km) {
  p <- terra::vect(as.matrix(rec[c("x", "y")]), crs = "EPSG:4326")
  laea <- sprintf(paste0("+proj=laea +lat_0=%.4f +lon_0=%.4f +x_0=0 +y_0=0 ",
                         "+datum=WGS84 +units=m +no_defs"),
                  mean(rec$y), mean(rec$x))
  xy <- terra::crds(terra::project(p, laea)) / 1000
  key <- paste(rec$group, floor(xy[, 1] / km), floor(xy[, 2] / km))
  rec[!duplicated(key), , drop = FALSE]
}

#' @export
print.overlap_scenarios <- function(x, digits = 2, ...) {
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(v, digits = digits, format = "f")
  r <- x$results
  ns <- nrow(unique(r[c("background", "space", "thin_km", "grid")]))
  cat("Niche overlap across ", ns, " scenarios (varying: ",
      paste(x$settings$varying, collapse = ", "), ")\n\nIN SHORT\n",
      sep = "")
  lab <- c(equivalency = "Equivalency test", similarity_12 = "Similarity",
           similarity_21 = "Similarity")
  for (key in unique(paste(r$group1, r$group2))) {
    d <- r[paste(r$group1, r$group2) == key, , drop = FALSE]
    g1 <- d$group1[1]
    g2 <- d$group2[1]
    rg <- range(d$D)
    say(paste0(g1, " vs ", g2, ": D ranges from ", f(rg[1]), " to ", f(rg[2]),
               " across the scenarios (", if (overlap_grade(rg[1]) ==
                                               overlap_grade(rg[2]))
                 paste0(overlap_grade(rg[1]), " overlap throughout") else
                   paste0(overlap_grade(rg[1]), " to ", overlap_grade(rg[2]),
                          " overlap"), ")."))
    rb <- x$robustness[x$robustness$group1 == g1 &
                         x$robustness$group2 == g2, , drop = FALSE]
    for (i in seq_len(nrow(rb))) {
      t <- rb[i, ]
      what <- switch(t$test, equivalency = "Equivalency test",
                     similarity_12 = paste0("Similarity of ", g1, " to ", g2),
                     similarity_21 = paste0("Similarity of ", g2, " to ", g1))
      if (t$robust) {
        say(paste0(what, ": the same verdict in all ", t$scenarios,
                   " scenarios (", sub(": .*", "", t$verdicts), "): robust ",
                   "to these choices."))
      } else {
        say(paste0(what, ": NOT robust (", t$verdicts, "). The verdict ",
                   "depends on the analysis, mostly on ", t$driver,
                   " (it changes in ", share_label(t$driver_flips),
                   " of the comparisons that differ only in ", t$driver,
                   "); report it as uncertain."))
      }
    }
  }
  cat("\n")
  show <- r
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits + 1)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "Verdicts at alpha = ", x$settings$alpha, ". A conclusion that holds ",
    "in every scenario does not depend on these choices; one that changes ",
    "is a result of the choices as much as of the organisms (Lo Parrino et ",
    "al. 2023)."), width = 78), sep = "\n")
  invisible(x)
}

#' Plot niche overlap in environmental space
#'
#' The environments available to the two groups (grey; darker where both
#' have them) and their niches in the shared space: environments occupied
#' by both (stability), only by the first group (unfilling) and only by the
#' second (expansion), with the 50% contour of each group's occupancy.
#' Legend entries give the number and share of occupied cells.
#'
#' @param x A result of [niche_overlap()].
#' @param pair Which pair of groups (row of `x$pairs`).
#' @param limits The part of the space shown: `"available"` (default, all
#'   environments available to either group), `"niches"` (zoom to the
#'   occupied envelopes of the pair, for niches much narrower than the
#'   available environment), or `c(xmin, xmax, ymin, ymax)` in axis units
#'   (e.g. to show several pairs or analyses on the same scale).
#' @param margin Space added around the `"available"` or `"niches"`
#'   window, as a share of its width and height (not beyond the grid of
#'   the densities, which ends at the most extreme record or background
#'   point).
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Broennimann, O. et al. (2012) Measuring ecological niche
#'   overlap from occurrence and spatial environmental data. *Global Ecology
#'   and Biogeography* 21, 481-497. \doi{10.1111/j.1466-8238.2011.00698.x}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' pres <- sim$presence
#' pres$group <- ifelse(pres$lon < 19, "west", "east")
#' plot_niche_overlap(niche_overlap(pres, sim$env, n_rep = 0))
plot_niche_overlap <- function(x, pair = 1, limits = "available",
                               margin = 0.1, base_size = 11, file = NULL,
                               width = 7.5, height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "niche_overlap")) {
    stop("`x` must come from niche_overlap().", call. = FALSE)
  }
  if (!(is.numeric(limits) && length(limits) == 4) &&
      !(is.character(limits) && length(limits) == 1 &&
        limits %in% c("available", "niches"))) {
    stop("`limits` must be \"available\", \"niches\" or c(xmin, xmax, ",
         "ymin, ymax).", call. = FALSE)
  }
  p <- x$pairs[pair, ]
  a <- x$densities[[p$group1]]
  b <- x$densities[[p$group2]]
  e <- x$space$extent
  g <- nrow(a$Z)
  gx <- seq(e[1], e[2], length.out = g)
  gy <- seq(e[3], e[4], length.out = g)
  d <- expand.grid(i = seq_len(g), j = seq_len(g))
  d$x <- gx[d$i]
  d$y <- gy[d$j]
  avail <- (a$Z > 0) + (b$Z > 0)
  w1 <- a$w
  w2 <- b$w
  cat_v <- ifelse(w1 & w2, "Stability", ifelse(w1, paste0("Unfilling (only ",
                                                          p$group1, ")"),
                                               ifelse(w2, paste0(
                                                 "Expansion (only ", p$group2,
                                                 ")"), NA)))
  d$avail <- as.vector(avail)
  d$cat <- as.vector(cat_v)
  occ <- d[!is.na(d$cat), , drop = FALSE]
  lv <- c(paste0("Unfilling (only ", p$group1, ")"), "Stability",
          paste0("Expansion (only ", p$group2, ")"))
  cnt <- table(factor(occ$cat, lv))
  lab <- paste0(lv, ": ", as.vector(cnt), " cells (",
                share_label(as.vector(cnt) / max(1, sum(cnt))), ")")
  occ$cat <- factor(occ$cat, lv, lab)
  # ColorBrewer RdYlBu (diverging): the first group's niche blue, the
  # second's red, the shared niche in the middle (yellow)
  pal <- stats::setNames(c("#4575b4", "#fee090", "#d73027"), lab)
  dc <- d
  dc$z1 <- as.vector(a$z_uncor)
  dc$z2 <- as.vector(b$z_uncor)
  # plotted window: the whole space, the niches of the pair, or given
  lim <- if (is.numeric(limits)) limits else if (limits == "available") {
    av <- d[d$avail > 0, , drop = FALSE]
    c(range(av$x), range(av$y))
  } else {
    oc <- d[as.vector(w1 | w2), , drop = FALSE]
    if (!nrow(oc)) oc <- d
    c(range(oc$x), range(oc$y))
  }
  if (!is.numeric(limits)) {
    dx <- diff(lim[1:2]) * margin
    dy <- diff(lim[3:4]) * margin
    lim <- c(max(e[1], lim[1] - dx), min(e[2], lim[2] + dx),
             max(e[3], lim[3] - dy), min(e[4], lim[4] + dy))
  }
  axl <- if (!is.null(x$space$variance)) {
    paste0(x$space$axes, " (", round(100 * x$space$variance, 1), "%)")
  } else x$space$axes
  ggplot2::ggplot() +
    ggplot2::geom_raster(data = d[d$avail > 0, ],
                         ggplot2::aes(.data$x, .data$y,
                                      alpha = factor(.data$avail)),
                         fill = "#969696") +
    ggplot2::scale_alpha_manual(values = c("1" = 0.25, "2" = 0.5),
                                labels = c("1" = "Available to one group",
                                           "2" = "Available to both"),
                                name = NULL) +
    ggplot2::geom_raster(data = occ, ggplot2::aes(.data$x, .data$y,
                                                  fill = .data$cat)) +
    ggplot2::scale_fill_manual(values = pal, name = NULL) +
    ggplot2::geom_contour(data = dc, ggplot2::aes(.data$x, .data$y,
                                                  z = .data$z1),
                          breaks = 0.5, colour = "#313695",
                          linewidth = 0.6) +
    ggplot2::geom_contour(data = dc, ggplot2::aes(.data$x, .data$y,
                                                  z = .data$z2),
                          breaks = 0.5, colour = "#a50026",
                          linewidth = 0.6, linetype = "dashed") +
    ggplot2::labs(x = axl[1], y = axl[2],
                  title = paste0("Niches of ", p$group1, " and ", p$group2),
                  subtitle = paste0(
                    "D = ", sprintf("%.2f", p$D), ", I = ",
                    sprintf("%.2f", p$I), "; lines: 50% of the occupancy ",
                    "of ", p$group1, " (solid) and ", p$group2,
                    " (dashed).")) +
    ggplot2::coord_cartesian(xlim = lim[1:2], ylim = lim[3:4],
                             expand = FALSE) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
}
