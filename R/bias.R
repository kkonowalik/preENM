# Sensitivity of the data design to sampling-bias corrections -------------------

#' How do sampling-bias corrections change the data you will model?
#'
#' preENM offers several corrections of sampling bias (thinning in
#' geographic or environmental space, one record per cell, and backgrounds
#' weighted by environment, effort or a target group), and the literature
#' gives no rule for choosing among them: there is no optimal spatial
#' filtering distance (Lamboley & Fourcade 2024), thinning often does no
#' better than removing the same number of records at random (Ten Caten &
#' Dallas 2023), the best correction depends on the structure of the bias
#' and can even make predictions worse for narrow niches (Inman et al.
#' 2021), and without independent test data the effect of a correction on
#' model performance cannot be judged reliably; results with and without
#' correction should be presented (Baker et al. 2024).
#'
#' `bias_sensitivity()` therefore makes a statement about the *data*, not
#' about a model: it applies a grid of corrections and reports how each
#' changes the design. It fits nothing and reports no predictive
#' performance.
#'
#' **Records** (each thinning design):
#' * `n`, `n_share`: records kept, and their share of all records.
#' * `env_coverage`: the share of the environmental bins of M (a PCA of
#'   `ref_n` cells of M, cut into about `total_bins` bins) that hold a
#'   record; `niche_kept`: the share of the bins occupied by all records
#'   that are still occupied.
#' * `null_coverage`: the same coverage for random subsets of the same size
#'   (`n_null` draws), and `coverage_gain` = `env_coverage` -
#'   `null_coverage`: does the thinning keep more of the environment than
#'   losing the same number of records at random?
#' * `nn_km`, `clark_evans`: mean nearest-neighbour distance and its ratio
#'   to the expectation for random points in M (Clark-Evans ratio; below 1
#'   clustered, about 1 random, above 1 regular).
#' * `range_kept`: the smallest share, over the layers, of the range of the
#'   records' values that the kept records still span; `gap_max`,
#'   `n_truncated`: the largest gap between the records and the limits of M
#'   and the number of truncated ends, as [niche_truncation()] reports them.
#' * `shift`: 1 - Schoener's D between the environmental bins of the kept
#'   and of all records: 0 = the same distribution in environmental space,
#'   1 = no overlap. How far the design moved from the uncorrected one.
#' * `niche_intact`: `range_kept` and `niche_kept` both at least
#'   `intact_limit` (default 0.95; a rule of thumb).
#'
#' **Backgrounds** (each strategy and size, drawn with
#' [sample_background()]): `bg_coverage` (bins of M covered), `bg_shift`
#' (1 - D against a uniform background of the same size: how far the
#' strategy moves the background) and, for each combination with a
#' thinning design, `overlap` (Schoener's D between the kept records and
#' the background in environmental space: a corrected background that
#' carries the records' sampling bias overlaps them more).
#'
#' @param records Presence records: a data.frame with coordinates, a CSV
#'   path, `sf` or [terra::SpatVector].
#' @param env Predictor layers clipped to M ([terra::SpatRaster] or files).
#' @param thin Thinning methods: `"none"` (the uncorrected records),
#'   `"cell"` (one per raster cell, [thin_by_cell()]), `"distance"`
#'   (spatial filtering: no two records closer than each of `thin_km`) and
#'   `"environment"` ([thin_environment()] on the cell-thinned records, at
#'   each of `env_mult`).
#' @param thin_km Filtering distances (km) for `"distance"`.
#' @param env_mult Bin-width multipliers for `"environment"` (its
#'   `width_mult`).
#' @param background Named list of background strategies: each a list of
#'   arguments for [sample_background()], e.g. `list(uniform = list(),
#'   target = list(method = "target_group", target = tg))`.
#' @param bg_n Background sizes.
#' @param ref_n Cells sampled from M for its environmental space.
#' @param total_bins Approximate number of environmental bins.
#' @param n_null Random subsets per thinning design for the null.
#' @param intact_limit Limit of `range_kept` and `niche_kept` for an intact
#'   niche.
#' @param env_thin_args Further arguments for [thin_environment()].
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Random seed.
#' @param verbose Print progress.
#' @return An object of class `bias_sensitivity`: `designs` (one row per
#'   thinning design and background: the measures above), `thinning` (one
#'   row per thinning design), `backgrounds` (one row per background),
#'   `effects` (how much each choice changes each measure), `records` (the
#'   kept records of each thinning design) and `settings`. See
#'   [plot_bias_sensitivity()].
#' @references Baker, D.J., Maclean, I.M.D., Gaston, K.J. (2024) Effective
#'   strategies for correcting spatial sampling bias in species
#'   distribution models without independent test data. *Diversity and
#'   Distributions* 30, e13802. \doi{10.1111/ddi.13802}
#'
#'   Inman, R., Franklin, J., Esque, T., Nussear, K. (2021) Comparing sample
#'   bias correction methods for species distribution modeling using
#'   virtual species. *Ecosphere* 12, e03422. \doi{10.1002/ecs2.3422}
#'
#'   Ten Caten, C., Dallas, T. (2023) Thinning occurrence points does not
#'   improve species distribution model performance. *Ecosphere* 14,
#'   e4703. \doi{10.1002/ecs2.4703}
#'
#'   Lamboley, Q., Fourcade, Y. (2024) No optimal spatial filtering distance
#'   for mitigating sampling bias in ecological niche models. *Journal of
#'   Biogeography* 51, 1783-1794. \doi{10.1111/jbi.14854}
#' @seealso [thin_by_cell()], [thin_environment()], [sample_background()],
#'   [niche_truncation()], [m_scenarios()]
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_enm_data(seed = 1, n_presence = 300)
#' bs <- bias_sensitivity(sim$presence, sim$env, thin_km = c(10, 30),
#'                        env_mult = 1, bg_n = 1000, n_null = 3,
#'                        total_bins = 400, seed = 1, verbose = FALSE)
#' bs
#' plot_bias_sensitivity(bs)
#' }
bias_sensitivity <- function(records, env,
                             thin = c("none", "cell", "distance",
                                      "environment"),
                             thin_km = c(5, 10, 25, 50), env_mult = c(1, 2),
                             background = list(uniform = list()),
                             bg_n = c(1000, 10000), ref_n = 1e5,
                             total_bins = 1000, n_null = 10,
                             intact_limit = 0.95, env_thin_args = list(),
                             coords = NULL, points_crs = "EPSG:4326",
                             seed = 1, verbose = TRUE) {
  thin <- match.arg(thin, c("none", "cell", "distance", "environment"),
                    several.ok = TRUE)
  if (!"none" %in% thin) thin <- c("none", thin)
  say <- function(...) if (verbose) message(...)
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  r <- as_grid(env)
  if (is.null(names(background)) || any(!nzchar(names(background)))) {
    stop("`background` must be a named list of sample_background() ",
         "argument lists.", call. = FALSE)
  }
  # records with values, in the raster's CRS
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  p <- if (is.data.frame(records) && !inherits(records, "sf")) {
    df <- as.data.frame(records, stringsAsFactors = FALSE)
    xy <- if (is.null(coords)) guess_coords(names(df)) else coords
    x <- as_number(df[[xy[1]]])
    y <- as_number(df[[xy[2]]])
    ok <- !is.na(x) & !is.na(y)
    terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else as_points(records, coords, points_crs)
  ll <- terra::crds(terra::project(p, "EPSG:4326"))
  pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
  cell <- terra::cellFromXY(r, terra::crds(pr))
  E_all <- env_at_cells(r, cell)
  ok <- stats::complete.cases(E_all)
  if (sum(ok) < 10) stop("Fewer than 10 records with values in all layers.",
                         call. = FALSE)
  base <- data.frame(lon = ll[ok, 1], lat = ll[ok, 2], cell = cell[ok])
  E_all <- E_all[ok, , drop = FALSE]
  vars <- colnames(E_all)
  N <- nrow(base)

  # environment of M
  say("Sampling the environment of M")
  ref <- quiet(sample_background(r, n = ref_n, seed = seed))
  R <- env_at_cells(r, ref$cell)
  R <- R[stats::complete.cases(R), , drop = FALSE]
  sp <- env_space(R, total_bins)
  ref_bins <- unique(sp$ref_bin)
  m_area_km2 <- terra::expanse(r[[1]], unit = "km")$area
  rec_bin <- sp$bin_of(E_all)
  all_bins <- unique(rec_bin)
  schoener <- function(b1, b2) {
    u <- union(b1, b2)
    p1 <- tabulate(match(b1, u), length(u)) / length(b1)
    p2 <- tabulate(match(b2, u), length(u)) / length(b2)
    1 - 0.5 * sum(abs(p1 - p2))
  }
  rng_all <- apply(E_all, 2, function(v) diff(range(v)))
  geo <- terra::is.lonlat(r)

  # thinning designs: indices of the kept records
  say("Thinning")
  designs <- list(none = seq_len(N))
  if ("cell" %in% thin || "environment" %in% thin) {
    cl <- which(!duplicated(base$cell))
    if ("cell" %in% thin) designs$cell <- cl
  }
  if ("distance" %in% thin) {
    X3 <- km_space(base$lon, base$lat)
    for (d in sort(thin_km)) {
      designs[[paste0("distance_", d, "km")]] <-
        with_seed(seed, thin_distance(X3, d))
    }
  }
  if ("environment" %in% thin) {
    ea <- env_thin_args
    if (is.null(ea$radius_km) && is.null(ea$window)) {
      # records at cell resolution: accurate to about one cell
      rs <- terra::res(r)
      ea$radius_km <- if (geo) {
        sqrt(rs[1] * 111.32 * cos(stats::median(base$lat) * pi / 180) *
               rs[2] * 111.32)
      } else sqrt(prod(rs)) / 1000
    }
    for (k in env_mult) {
      th <- tryCatch(quiet(do.call(thin_environment, c(list(
        points = base[cl, c("lon", "lat")], env = r, width_mult = k,
        seed = seed, verbose = FALSE), ea))),
        error = function(e) {
          warning("thin_environment() failed at width_mult = ", k, ": ",
                  conditionMessage(e), call. = FALSE)
          NULL
        })
      if (!is.null(th)) {
        designs[[paste0("environment_x", k)]] <- cl[th$data$kept]
      }
    }
  }

  rec_measures <- function(idx) {
    E <- E_all[idx, , drop = FALSE]
    b <- rec_bin[idx]
    n <- length(idx)
    nn <- if (n > 1) {
      X <- km_space(base$lon[idx], base$lat[idx])
      mean(dbscan::kNN(X, k = 1)$dist[, 1])
    } else NA_real_
    expct <- 0.5 * sqrt(m_area_km2 / n)
    tr <- tryCatch(quiet(niche_truncation(base[idx, c("lon", "lat")], r,
                                          n = min(ref_n, 2e4),
                                          seed = seed)$variables),
                   error = function(e) NULL)
    rng <- apply(E, 2, function(v) diff(range(v)))
    data.frame(
      n = n, n_share = n / N,
      env_coverage = mean(ref_bins %in% b),
      niche_kept = mean(all_bins %in% b),
      nn_km = nn, clark_evans = nn / expct,
      range_kept = min(ifelse(rng_all > 0, rng / rng_all, 1)),
      gap_max = if (is.null(tr)) NA_real_ else
        max(tr$gap_low, tr$gap_high),
      n_truncated = if (is.null(tr)) NA_integer_ else
        sum(tr$truncated_low) + sum(tr$truncated_high),
      shift = 1 - schoener(b, rec_bin))
  }
  null_coverage <- function(n) {
    if (n >= N) return(mean(ref_bins %in% rec_bin))
    mean(vapply(seq_len(n_null), function(i) {
      s <- with_seed(seed + i, sample.int(N, n))
      mean(ref_bins %in% rec_bin[s])
    }, numeric(1)))
  }
  thin_tab <- do.call(rbind, lapply(names(designs), function(k) {
    m <- rec_measures(designs[[k]])
    m$null_coverage <- null_coverage(m$n)
    m$coverage_gain <- m$env_coverage - m$null_coverage
    cbind(thinning = k,
          thin_method = sub("_.*", "", k),
          thin_level = ifelse(grepl("_", k), sub("^[^_]*_", "", k), ""),
          m, stringsAsFactors = FALSE)
  }))
  thin_tab$niche_intact <- thin_tab$range_kept >= intact_limit &
    thin_tab$niche_kept >= intact_limit

  # backgrounds
  say("Drawing backgrounds")
  bg_bins <- list()
  bg_rows <- list()
  unif_bins <- list()
  for (n in bg_n) {
    u <- quiet(sample_background(r, n = n, seed = seed + 7))
    unif_bins[[as.character(n)]] <- sp$bin_of(env_at_cells(r, u$cell))
  }
  for (nm in names(background)) for (n in bg_n) {
    args <- background[[nm]]
    args$n <- n
    if (is.null(args$seed)) args$seed <- seed + 7
    bg <- quiet(do.call(sample_background, c(list(raster = r), args)))
    B <- env_at_cells(r, bg$cell)
    B <- B[stats::complete.cases(B), , drop = FALSE]
    b <- sp$bin_of(B)
    key <- paste0(nm, "_n", n)
    bg_bins[[key]] <- b
    bg_rows[[key]] <- data.frame(
      background = key, bg_method = nm, bg_n = n, bg_drawn = nrow(B),
      bg_coverage = mean(ref_bins %in% b),
      bg_shift = 1 - schoener(b, unif_bins[[as.character(n)]]),
      stringsAsFactors = FALSE)
  }
  bg_tab <- do.call(rbind, bg_rows)
  rownames(bg_tab) <- NULL

  # every combination
  des <- merge(thin_tab, bg_tab, by = NULL)
  des$overlap <- mapply(function(t, b) {
    schoener(rec_bin[designs[[t]]], bg_bins[[b]])
  }, des$thinning, des$background)
  des$design <- paste(des$thinning, des$background, sep = " | ")
  des <- des[c("design", setdiff(names(des), "design"))]

  # how much each choice changes each measure: the range of the measure
  # over the levels of the choice, with the other choices held fixed,
  # relative to its range over all designs; each choice only over the
  # measures it can change
  rec_m <- c("n", "env_coverage", "niche_kept", "clark_evans", "range_kept",
             "shift", "overlap")
  bg_m <- c("bg_coverage", "bg_shift", "overlap")
  effect <- function(label, d, level, group, ms) {
    if (!nrow(d) || length(unique(level)) < 2) return(NULL)
    do.call(rbind, lapply(ms, function(m) {
      # mean of the measure per level within each group, then its range
      v <- tapply(d[[m]], list(group, level), mean, na.rm = TRUE)
      rg <- apply(v, 1, function(z) {
        z <- z[!is.na(z)]
        if (length(z) < 2) NA_real_ else diff(range(z))
      })
      sc <- diff(range(des[[m]], na.rm = TRUE))
      ch <- mean(rg, na.rm = TRUE)
      data.frame(choice = label, measure = m, change = ch,
                 relative = if (is.finite(sc) && sc > 0 && is.finite(ch))
                   ch / sc else 0, stringsAsFactors = FALSE)
    }))
  }
  bgkey <- paste(des$bg_method, des$bg_n)
  thk <- des$thin_level != ""
  eff <- rbind(
    effect("thinning method", des, des$thin_method, bgkey, rec_m),
    effect("thinning radius / width", des[thk, , drop = FALSE],
           des$thin_level[thk],
           paste(des$thin_method[thk], bgkey[thk]), rec_m),
    effect("background strategy", des, des$bg_method,
           paste(des$thinning, des$bg_n), bg_m),
    effect("background size", des, as.character(des$bg_n),
           paste(des$thinning, des$bg_method), bg_m))
  rownames(des) <- NULL
  rownames(thin_tab) <- NULL
  structure(list(designs = des, thinning = thin_tab, backgrounds = bg_tab,
                 effects = eff,
                 records = lapply(designs, function(i) base[i, ]),
                 settings = list(n_records = N, thin = thin,
                                 thin_km = thin_km, env_mult = env_mult,
                                 bg_n = bg_n, background = names(background),
                                 total_bins = total_bins,
                                 bins = length(ref_bins),
                                 intact_limit = intact_limit,
                                 n_null = n_null,
                                 m_area_km2 = m_area_km2)),
            class = "bias_sensitivity")
}

# Longitude/latitude as 3-D coordinates in km (chord distances, close to
# great-circle distances at the scale of thinning) (internal)
km_space <- function(lon, lat) {
  lo <- lon * pi / 180
  la <- lat * pi / 180
  cbind(cos(la) * cos(lo), cos(la) * sin(lo), sin(la)) * 6371.0088
}

# Spatial filtering (internal): records in random order, each kept unless a
# kept record lies within `d` km, so no two kept records are closer than d
thin_distance <- function(X, d) {
  n <- nrow(X)
  if (n < 2) return(seq_len(n))
  nb <- dbscan::frNN(X, eps = d)$id
  ord <- sample.int(n)
  removed <- logical(n)
  keep <- logical(n)
  for (i in ord) {
    if (removed[i]) next
    keep[i] <- TRUE
    removed[nb[[i]]] <- TRUE
  }
  which(keep)
}

#' @export
print.bias_sensitivity <- function(x, digits = 2, ...) {
  t <- x$thinning
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  pct <- function(v) share_label(v)
  cat("Bias-correction sensitivity: ", nrow(t), " thinning design(s) x ",
      nrow(x$backgrounds), " background(s); ", s$n_records,
      " records\n\nIN SHORT\n", sep = "")
  u <- t[t$thinning == "none", ]
  say(paste0("Uncorrected: ", u$n, " records covering ",
             pct(u$env_coverage), " of the environmental bins of M ",
             "(Clark-Evans ratio ", round(u$clark_evans, 2),
             if (u$clark_evans < 0.8) ", clustered" else
               if (u$clark_evans > 1.2) ", regular" else ", near random",
             ")."))
  ef <- x$effects
  if (!is.null(ef) && nrow(ef)) {
    agg <- tapply(ef$relative, ef$choice, mean, na.rm = TRUE)
    agg <- sort(agg, decreasing = TRUE)
    say(paste0("The choice that changes the design most is the ",
               names(agg)[1], if (length(agg) > 1) paste0(
                 "; the one that matters least, the ", names(agg)[length(agg)])
               else "", " (mean relative change over the measures: ",
               paste0(names(agg), " ", round(agg, 2), collapse = ", "),
               ")."))
    small <- names(agg)[agg < 0.1]
    if (length(small)) {
      say(paste0("Barely matters here: ", paste(small, collapse = ", "),
                 "."))
    }
  }
  th <- t[t$thinning != "none", , drop = FALSE]
  if (nrow(th)) {
    intact <- th$thinning[th$niche_intact]
    broken <- th$thinning[!th$niche_intact]
    say(paste0(if (length(intact)) paste0(
      "Leave the niche intact (range and occupied bins kept at least ",
      s$intact_limit * 100, "%): ", paste(intact, collapse = ", "), ". ")
      else "No thinning design leaves the niche fully intact. ",
      if (length(broken)) paste0("Cut the niche: ", paste0(
        broken, " (range ", pct(th$range_kept[!th$niche_intact]),
        ", bins ", pct(th$niche_kept[!th$niche_intact]), ")",
        collapse = "; "), ".") else ""))
    none_removed <- th$thinning[th$n >= s$n_records]
    th2 <- th[th$n < s$n_records, , drop = FALSE]
    better <- th2$thinning[th2$coverage_gain > 0.01]
    worse <- th2$thinning[th2$coverage_gain < -0.01]
    same <- setdiff(th2$thinning, c(better, worse))
    if (length(none_removed)) {
      say(paste0("Removed no record: ", paste(none_removed, collapse = ", "),
                 " (the records are already one per cell)."))
    }
    say(paste0("Against random removal of the same number of records: ",
               if (length(better)) paste0("more environment kept by ",
                                          paste(better, collapse = ", "))
               else "no design keeps more environment",
               if (length(same)) paste0("; no different from random: ",
                                        paste(same, collapse = ", "))
               else "",
               if (length(worse)) paste0("; less than random: ",
                                         paste(worse, collapse = ", "))
               else "", "."))
  }
  say(paste0("This describes the data, not model performance: no ",
             "correction can be recommended from first principles (no ",
             "optimal filtering distance, Lamboley & Fourcade 2024; ",
             "thinning often no better than random removal, Ten Caten & ",
             "Dallas 2023; the best correction depends on the bias, Inman ",
             "et al. 2021). Without independent test data, fit models with ",
             "and without correction and report both (Baker et al. 2024)."))
  cat("\nThinning designs:\n")
  show <- t[c("thinning", "n", "env_coverage", "niche_kept",
              "coverage_gain", "clark_evans", "range_kept", "gap_max",
              "n_truncated", "shift", "niche_intact")]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("\nBackgrounds:\n")
  b <- x$backgrounds
  b$bg_coverage <- round(b$bg_coverage, digits)
  b$bg_shift <- round(b$bg_shift, digits)
  print(b, row.names = FALSE)
  cat("", strwrap(paste0(
    "env_coverage = share of the ", s$bins, " environmental bins of M ",
    "holding a record; niche_kept = share of the bins of all records still ",
    "occupied; coverage_gain = env_coverage minus that of random subsets ",
    "of the same size (", s$n_null, " draws); clark_evans = mean ",
    "nearest-neighbour distance / its expectation for random points in M ",
    "(below 1 clustered); range_kept = smallest share of a layer's range ",
    "still spanned; gap_max, n_truncated as in niche_truncation(); shift = ",
    "1 - Schoener's D between kept and all records in environmental space; ",
    "bg_shift = 1 - D between the background and a uniform one of the same ",
    "size; overlap (in x$designs) = D between kept records and background. ",
    "The limits (", s$intact_limit, " for an intact niche, 0.01 for a gain, ",
    "0.1 for a choice that barely matters) are rules of thumb."),
    width = 76), sep = "\n")
  invisible(x)
}

#' Plot the designs of a bias-correction sensitivity analysis
#'
#' One point per design (thinning x background) on two chosen measures,
#' coloured by thinning method and shaped by background, with the
#' uncorrected design ringed: the designs side by side, as
#' [plot_m_area()] shows the scenarios of [m_scenarios()].
#'
#' @param x A `bias_sensitivity` from [bias_sensitivity()].
#' @param x_axis,y_axis Measures to plot (columns of `x$designs`), e.g.
#'   `"n"`, `"env_coverage"`, `"niche_kept"`, `"clark_evans"`,
#'   `"range_kept"`, `"shift"`, `"overlap"`, `"bg_coverage"`.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Baker, D.J., Maclean, I.M.D., Gaston, K.J. (2024) Effective
#'   strategies for correcting spatial sampling bias in species
#'   distribution models without independent test data. *Diversity and
#'   Distributions* 30, e13802. \doi{10.1111/ddi.13802}
#' @export
#' @examples
#' # see ?bias_sensitivity
plot_bias_sensitivity <- function(x, x_axis = "n", y_axis = "env_coverage",
                                  base_size = 11, file = NULL, width = 7.5,
                                  height = 5.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "bias_sensitivity")) {
    stop("`x` must come from bias_sensitivity().", call. = FALSE)
  }
  d <- x$designs
  for (a in c(x_axis, y_axis)) {
    if (!a %in% names(d) || !is.numeric(d[[a]])) {
      stop("Not a numeric measure of the designs: ", a, call. = FALSE)
    }
  }
  lab <- c(n = "Records kept", n_share = "Share of records kept",
           env_coverage = "Environmental bins of M with records (share)",
           niche_kept = "Bins of all records still occupied (share)",
           coverage_gain = "Coverage gain over random removal",
           clark_evans = "Clark-Evans ratio (below 1 clustered)",
           nn_km = "Mean nearest-neighbour distance (km)",
           range_kept = "Smallest share of a layer's range kept",
           shift = "Shift from all records (1 - D)",
           overlap = "Overlap of records and background (D)",
           bg_coverage = "Environmental bins of M with background (share)",
           bg_shift = "Background shift from uniform (1 - D)",
           gap_max = "Largest gap to the limit of M")
  nice <- function(a) if (a %in% names(lab)) lab[[a]] else a
  meth <- unique(x$thinning$thin_method)
  bg_dep0 <- c("overlap", "bg_coverage", "bg_shift", "bg_n", "bg_drawn")
  dd <- if (any(c(x_axis, y_axis) %in% bg_dep0)) d else
    d[!duplicated(d$thinning), , drop = FALSE]
  cnt <- table(factor(dd$thin_method, meth))
  # ColorBrewer Dark2
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02")[seq_along(meth)]
  ml <- paste0(meth, " (", as.integer(cnt), " design",
               ifelse(as.integer(cnt) == 1, "", "s"), ")")
  d$method <- factor(d$thin_method, meth, ml)
  d$bg <- factor(d$background, unique(d$background))
  shp <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0)[seq_along(levels(d$bg))]
  # measures of the records alone are the same for every background: one
  # point per thinning design then
  bg_dep <- c("overlap", "bg_coverage", "bg_shift", "bg_n", "bg_drawn")
  by_bg <- any(c(x_axis, y_axis) %in% bg_dep)
  if (!by_bg) d <- d[!duplicated(d$thinning), , drop = FALSE]
  u <- d[d$thin_method == "none", , drop = FALSE]
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[x_axis]], .data[[y_axis]])) +
    ggplot2::geom_point(data = u, shape = 21, size = 5, stroke = 0.8,
                        colour = "#3a3a38", fill = NA) +
    (if (by_bg) ggplot2::geom_point(ggplot2::aes(colour = .data$method,
                                                 shape = .data$bg),
                                    size = 2.6) else
      ggplot2::geom_point(ggplot2::aes(colour = .data$method), size = 2.8)) +
    labels_layer(d[d$thin_level != "", , drop = FALSE], base_size) +
    ggplot2::scale_colour_manual(values = stats::setNames(pal, ml),
                                 name = "Thinning") +
    ggplot2::labs(x = nice(x_axis), y = nice(y_axis),
                  title = "Data designs under sampling-bias corrections",
                  subtitle = paste0(
                    nrow(d), if (by_bg) " designs" else " thinning designs",
                    "; ringed: uncorrected records; labels: radius or ",
                    "width")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
  if (by_bg) {
    p <- p + ggplot2::scale_shape_manual(values = stats::setNames(
      shp, levels(d$bg)), name = "Background")
  }
  p
}

# Labels of the radius or width of each design, kept apart with ggrepel
# when it is installed (internal)
labels_layer <- function(d, base_size) {
  if (!nrow(d)) return(NULL)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = d, ggplot2::aes(label = .data$thin_level),
                             size = base_size / 3.6, colour = "#3a3a38",
                             min.segment.length = 0.2, seed = 1,
                             max.overlaps = Inf)
  } else {
    ggplot2::geom_text(data = d, ggplot2::aes(label = .data$thin_level),
                       size = base_size / 3.6, vjust = -1.1,
                       colour = "#3a3a38", check_overlap = TRUE)
  }
}
