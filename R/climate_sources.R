# Comparing climate products before choosing one

#' Compare climate products at the records and across M
#'
#' Climate products (WorldClim, CHELSA, regional atlases, reanalyses) are
#' built with different station networks, interpolation methods and
#' periods, and they disagree most where stations are sparse: mountains,
#' deserts, high latitudes. Models fitted with different products can
#' differ, and some algorithms are more sensitive to the choice than others
#' (Rodriguez-Rey & Jimenez-Valverde 2024). This function shows that
#' exposure before a product is chosen: it samples two or more products at
#' the occurrence records and at a random sample of cells of the
#' calibration area M, and reports, for every variable and every product
#' against a reference:
#' * **agreement**: Lin's concordance correlation (CCC; Lin 1989), which is
#'   1 only when the values are identical. It falls with scatter and with
#'   any offset or change of scale.
#' * **correlation**: Pearson and Spearman, as in the comparisons of
#'   Rodriguez-Rey & Jimenez-Valverde (2024). A high correlation with a
#'   low CCC means the products rank sites alike but differ in level or
#'   scale.
#' * **systematic offset**: the mean difference (product minus reference),
#'   in the variable's units and in standard deviations of the reference
#'   across M, and the RMSE.
#' * **distributions**: the distribution-free overlapping index of Pastore
#'   & Calcagni (2019), the area shared by the two densities of values (1 =
#'   the same distribution), and the share of the reference's range that
#'   the product's range covers (at the records and in M). Both were used by Rodriguez-Rey &
#'   Jimenez-Valverde (2024).
#' * **collinearity structure**: the Mantel correlation between the two
#'   products' matrices of correlations among the variables, with a
#'   permutation test (as in Jimenez-Valverde et al. 2009 and
#'   Rodriguez-Rey & Jimenez-Valverde 2024). A low value means variable
#'   selection on one product does not carry over to the other.
#' * **geographic pattern**: the standardised difference at every point,
#'   its rank correlation with longitude and latitude, and whether the
#'   records sit where the products disagree more than they do on average
#'   across M. [plot_climate_sources()] maps it.
#'
#' The function is generic: any rasters with matching variables can be
#' compared. Layers named after bioclimatic variables in the usual ways
#' (`bio1`, `bio_01`, `wc2.1_30s_bio_1`, `CHELSA_bio1_1981-2010_V.2.1`) are
#' matched automatically when `harmonize = TRUE`; other layers must have the
#' same names in every product. Products need not share a grid or CRS: each
#' is read at the same points with `terra::extract()`, so only those cells
#' are read, even from global 30 arc-second layers.
#'
#' Check that the products use the same units and scale (degrees C vs
#' kelvin, or values multiplied by 10) before comparing: an offset larger
#' than three standard deviations is reported as a likely mismatch.
#'
#' @param sources A named list of two or more climate products, each a
#'   `SpatRaster` (or a vector of file paths) with one layer per variable.
#' @param records Occurrence records: data frame with coordinates, CSV path,
#'   `sf` or `SpatVector`. Optional; without them only M is compared.
#' @param m Calibration area (`sf`, `SpatVector`, file path or an
#'   `m_area`). If `NULL`, the extent of the reference product.
#' @param reference Name or position of the reference product in
#'   `sources` (e.g. a regional atlas or the product used so far).
#' @param variables Variables to compare (after harmonising names). Default:
#'   all variables present in every product.
#' @param harmonize Match bioclimatic layer names across naming schemes?
#' @param n Number of cells of M sampled.
#' @param n_perm Permutations for the Mantel test (0 to skip).
#' @param coords Coordinate columns of `records` (guessed if `NULL`).
#' @param points_crs CRS of the record coordinates.
#' @param seed Random seed for the sample of M and the permutations.
#' @return An object of class `climate_source_compare`: a list with
#'   `metrics` (one row per variable, product and domain, `"records"` or
#'   `"M"`), `exposure` (one row per variable: the lowest agreement among
#'   the products and its grade), `collinearity` (Mantel r and p per
#'   product), `geography` (per variable and product: rank correlations of
#'   the standardised difference with longitude and latitude, and the ratio
#'   of the mean absolute difference at the records to that across M),
#'   `points` (the values of every product at every point, with `domain`,
#'   `lon` and `lat`) and `settings`.
#' @seealso [plot_climate_sources()], [write_climate_sources()];
#'   [climate_indices()] to compute the same indices from each product's
#'   monthly layers.
#' @references Jimenez-Valverde, A., Nakazawa, Y., Lira-Noriega, A.,
#'   Peterson, A.T. (2009) Environmental correlation structure and
#'   ecological niche model projections. *Biodiversity Informatics* 6,
#'   28-35. \doi{10.17161/bi.v6i1.1634}
#'
#'   Lin, L.I-K. (1989) A concordance correlation coefficient to evaluate
#'   reproducibility. *Biometrics* 45, 255-268. \doi{10.2307/2532051}
#'
#'   McBride, G.B. (2005) A proposal for strength-of-agreement criteria for
#'   Lin's concordance correlation coefficient. NIWA Client Report
#'   HAM2005-062. National Institute of Water & Atmospheric Research,
#'   Hamilton, New Zealand.
#'
#'   Pastore, M., Calcagni, A. (2019) Measuring distribution similarities
#'   between samples: a distribution-free overlapping index. *Frontiers in
#'   Psychology* 10, 1089. \doi{10.3389/fpsyg.2019.01089}
#'
#'   Rodriguez-Rey, M., Jimenez-Valverde, A. (2024) Differing sensitivity of
#'   species distribution modelling algorithms to climate data source.
#'   *Ecological Informatics* 79, 102387.
#'   \doi{10.1016/j.ecoinf.2023.102387}
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 80)
#' a <- sim$env
#' # a second "product": a regional bias in the north-east, plus noise
#' b <- a
#' x <- terra::init(a[[1]], "x")
#' y <- terra::init(a[[1]], "y")
#' bump <- exp(-((x - 22)^2 + (y - 54)^2) / 2)
#' set.seed(1)
#' noise <- terra::init(a[[1]], fun = stats::rnorm)
#' b[[1]] <- a[[1]] + 0.3 * bump + 0.03 * noise
#' b[[2]] <- a[[2]] * 1.1 + 0.02 * noise
#' cs <- climate_source_compare(list(A = a, B = b), sim$presence, n = 2000)
#' cs
#' plot_climate_sources(cs)
#' plot_climate_sources(cs, "map")
#' @export
climate_source_compare <- function(sources, records = NULL, m = NULL,
                                   reference = 1, variables = NULL,
                                   harmonize = TRUE, n = 5000, n_perm = 999,
                                   coords = NULL, points_crs = "EPSG:4326",
                                   seed = 1) {
  if (!is.list(sources) || length(sources) < 2) {
    stop("`sources` must be a list of two or more climate products.",
         call. = FALSE)
  }
  if (is.null(names(sources)) || any(!nzchar(names(sources)))) {
    names(sources) <- paste0("source", seq_along(sources))
  }
  if (anyDuplicated(names(sources))) {
    stop("The products in `sources` need distinct names.", call. = FALSE)
  }
  sources <- lapply(sources, function(s) {
    if (is.character(s)) s <- terra::rast(s)
    if (!inherits(s, "SpatRaster")) {
      stop("Each product must be a SpatRaster or file paths.", call. = FALSE)
    }
    if (harmonize) names(s) <- harmonize_bio(names(s))
    s
  })
  ref <- if (is.numeric(reference)) names(sources)[reference] else reference
  if (length(ref) != 1 || !ref %in% names(sources)) {
    stop("`reference` must name one product of `sources`.", call. = FALSE)
  }
  common <- Reduce(intersect, lapply(sources, names))
  if (is.null(variables)) variables <- common
  miss <- setdiff(variables, common)
  if (length(miss) || !length(variables)) {
    stop("Variables not found in every product: ",
         if (length(miss)) paste(miss, collapse = ", ") else "none shared",
         ". Layer names: ",
         paste(vapply(names(sources), function(k)
           paste0(k, " (", paste(utils::head(names(sources[[k]]), 6),
                                 collapse = ", "), ")"), ""),
           collapse = "; "),
         ". Rename the layers so they match.", call. = FALSE)
  }
  r0 <- sources[[ref]][[variables[1]]]
  # cells of M, on the grid of the reference
  cells <- suppressWarnings(suppressMessages(
    sample_background(r0, n = n, seed = seed)))
  if (!is.null(m)) {
    g <- if (inherits(m, "m_area")) m$m else as_sfc(m, "m")
    pts <- sf::st_as_sf(data.frame(x = cells$lon, y = cells$lat),
                        coords = c("x", "y"), crs = 4326)
    pts <- sf::st_transform(pts, sf::st_crs(g))
    cells <- cells[lengths(sf::st_intersects(pts, g)) > 0, , drop = FALSE]
  }
  if (!nrow(cells)) stop("No cell of the reference lies inside `m`.",
                         call. = FALSE)
  xy <- data.frame(lon = cells$lon, lat = cells$lat, domain = "M")
  n_rec <- 0
  if (!is.null(records)) {
    p <- as_points(records, coords, points_crs)
    if (!terra::is.lonlat(p)) p <- terra::project(p, "EPSG:4326")
    pc <- terra::crds(p)
    n_rec <- nrow(pc)
    xy <- rbind(data.frame(lon = pc[, 1], lat = pc[, 2], domain = "records"),
                xy)
  }
  vals <- lapply(names(sources), function(k) {
    v <- values_xy(sources[[k]][[variables]], as.matrix(xy[, 1:2]),
                   "EPSG:4326")
    colnames(v) <- variables
    v
  })
  names(vals) <- names(sources)
  pts <- xy
  for (k in names(sources)) {
    for (v in variables) pts[[paste0(v, "__", k)]] <- vals[[k]][, v]
  }
  others <- setdiff(names(sources), ref)
  domains <- if (n_rec) c("records", "M") else "M"
  inM <- xy$domain == "M"
  rows <- list()
  geo <- list()
  for (v in variables) {
    a <- vals[[ref]][, v]
    sd_m <- stats::sd(a[inM], na.rm = TRUE)
    for (k in others) {
      b <- vals[[k]][, v]
      d <- (b - a) / if (is.finite(sd_m) && sd_m > 0) sd_m else 1
      for (dm in domains) {
        i <- xy$domain == dm & !is.na(a) & !is.na(b)
        rows[[length(rows) + 1]] <- cbind(
          data.frame(variable = v, source = k, reference = ref, domain = dm,
                     stringsAsFactors = FALSE),
          agreement_metrics(a[i], b[i], sd_m,
                            if (any(i)) range(a[i]) else c(NA, NA)))
      }
      ok <- inM & is.finite(d)
      okr <- xy$domain == "records" & is.finite(d)
      cs <- function(z) if (sum(ok) > 2 && stats::sd(d[ok]) > 0 &&
                            stats::sd(z[ok]) > 0)
        suppressWarnings(stats::cor(d[ok], z[ok], method = "spearman")) else
          NA_real_
      mm <- mean(abs(d[ok]))
      geo[[length(geo) + 1]] <- data.frame(
        variable = v, source = k,
        rho_lon = cs(xy$lon), rho_lat = cs(xy$lat),
        mean_abs_M = mm,
        mean_abs_records = if (any(okr)) mean(abs(d[okr])) else NA_real_,
        records_ratio = if (any(okr) && mm > 0) mean(abs(d[okr])) / mm else
          NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  metrics <- do.call(rbind, rows)
  geography <- do.call(rbind, geo)
  # collinearity structure across M
  coll <- do.call(rbind, lapply(others, function(k) {
    out <- data.frame(source = k, reference = ref, mantel_r = NA_real_,
                      p = NA_real_)
    if (length(variables) < 3) return(out)
    A <- vals[[ref]][inM, variables, drop = FALSE]
    B <- vals[[k]][inM, variables, drop = FALSE]
    ok <- stats::complete.cases(A, B)
    if (sum(ok) < 5) return(out)
    ca <- suppressWarnings(stats::cor(A[ok, ]))
    cb <- suppressWarnings(stats::cor(B[ok, ]))
    mt <- mantel_cor(ca, cb, n_perm, seed)
    out$mantel_r <- mt[1]
    out$p <- mt[2]
    out
  }))
  # exposure per variable: the worst product, at the records if given
  dom <- if (n_rec) "records" else "M"
  ex <- do.call(rbind, lapply(variables, function(v) {
    mv <- metrics[metrics$variable == v & metrics$domain == dom, ,
                  drop = FALSE]
    w <- which.min(replace(mv$ccc, is.na(mv$ccc), Inf))
    if (!length(w)) w <- 1
    data.frame(variable = v, worst_source = mv$source[w], ccc = mv$ccc[w],
               r = mv$r[w], offset_sd = mv$offset_sd[w],
               grade = ccc_grade(mv$ccc[w]),
               unit_mismatch = any(abs(mv$offset_sd) > 3, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  ex <- ex[order(ex$ccc), , drop = FALSE]
  rownames(ex) <- NULL
  structure(list(
    metrics = metrics, exposure = ex, collinearity = coll,
    geography = geography, points = pts,
    settings = list(sources = names(sources), reference = ref,
                    variables = variables, n_records = n_rec,
                    n_cells = sum(inM), domain = dom, n_perm = n_perm,
                    seed = seed)),
    class = "climate_source_compare")
}

# Agreement of b with a (internal)
agreement_metrics <- function(a, b, sd_m, rng) {
  n <- length(a)
  na <- data.frame(n = n, r = NA_real_, rho = NA_real_, ccc = NA_real_,
                   offset = NA_real_, offset_sd = NA_real_, rmse = NA_real_,
                   rmse_sd = NA_real_, large = NA_real_, overlap = NA_real_,
                   coverage = NA_real_)
  if (n < 2) return(na)
  s <- if (is.finite(sd_m) && sd_m > 0) sd_m else NA_real_
  d <- b - a
  ma <- mean(a)
  mb <- mean(b)
  va <- mean((a - ma)^2)
  vb <- mean((b - mb)^2)
  cv <- mean((a - ma) * (b - mb))
  den <- va + vb + (ma - mb)^2
  out <- na
  out$ccc <- if (den > 0) 2 * cv / den else 1
  if (va > 0 && vb > 0) {
    out$r <- stats::cor(a, b)
    out$rho <- suppressWarnings(stats::cor(a, b, method = "spearman"))
  } else if (all(d == 0)) {
    out$r <- out$rho <- 1
  }
  out$offset <- mean(d)
  out$offset_sd <- mean(d) / s
  out$rmse <- sqrt(mean(d^2))
  out$rmse_sd <- out$rmse / s
  out$large <- mean(abs(d) / s > 0.5)
  out$overlap <- overlap_index(a, b)
  rb <- range(b)
  out$coverage <- if (diff(rng) > 0)
    max(0, min(rng[2], rb[2]) - max(rng[1], rb[1])) / diff(rng) else
      as.numeric(all(rb == rng))
  out
}

# Distribution-free overlapping index of two samples (Pastore & Calcagni
# 2019): the area under the minimum of the two kernel densities (internal)
overlap_index <- function(a, b, n_grid = 512) {
  if (identical(a, b)) return(1)
  lo <- min(a, b)
  hi <- max(a, b)
  if (hi == lo) return(1)
  bw <- function(z) {
    h <- stats::bw.nrd0(z)
    if (!is.finite(h) || h <= 0) (hi - lo) / 100 else h
  }
  pad <- 3 * max(bw(a), bw(b))
  fa <- stats::density(a, bw = bw(a), from = lo - pad, to = hi + pad,
                       n = n_grid)
  fb <- stats::density(b, bw = bw(b), from = lo - pad, to = hi + pad,
                       n = n_grid)
  dx <- fa$x[2] - fa$x[1]
  min(1, sum(pmin(fa$y, fb$y)) * dx)
}

# Mantel correlation of two correlation matrices with a permutation test
# (internal): r of the lower triangles, p from permuting the variables
mantel_cor <- function(ca, cb, n_perm, seed) {
  lt <- lower.tri(ca)
  r <- suppressWarnings(stats::cor(ca[lt], cb[lt]))
  if (!is.finite(r)) return(c(NA_real_, NA_real_))
  if (!n_perm) return(c(r, NA_real_))
  k <- ncol(ca)
  perm <- with_seed(seed, vapply(seq_len(n_perm), function(i) {
    o <- sample.int(k)
    suppressWarnings(stats::cor(ca[lt], cb[o, o][lt]))
  }, numeric(1)))
  c(r, (sum(perm >= r - 1e-12, na.rm = TRUE) + 1) / (n_perm + 1))
}

# Strength of agreement for Lin's CCC (McBride 2005) (internal)
ccc_grade <- function(x) {
  as.character(cut(x, c(-Inf, 0.90, 0.95, 0.99, Inf),
                   labels = c("poor", "moderate", "substantial",
                              "almost perfect"), right = FALSE))
}

# Common names for bioclimatic layers (internal): bio_01, wc2.1_30s_bio_1,
# CHELSA_bio1_1981-2010_V.2.1 -> bio1
harmonize_bio <- function(nm) {
  m <- regmatches(nm, regexec("(^|[_.-])bio_?0*([0-9]{1,2})([_.-]|$)", nm,
                              ignore.case = TRUE))
  vapply(seq_along(nm), function(i) {
    if (length(m[[i]]) == 4 && as.integer(m[[i]][3]) %in% 1:19) {
      paste0("bio", as.integer(m[[i]][3]))
    } else {
      nm[i]
    }
  }, character(1))
}

#' @export
print.climate_source_compare <- function(x, digits = 2, ...) {
  s <- x$settings
  e <- x$exposure
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  where <- if (s$domain == "records") "at the records" else "across M"
  cat("Comparison of ", length(s$sources), " climate products against ",
      s$reference, ": ", length(s$variables), " variable(s), ", s$n_records,
      " records, ", s$n_cells, " cells of M\n\nIN SHORT\n", sep = "")
  m <- x$metrics
  if (all(m$rmse == 0, na.rm = TRUE)) {
    say(paste0("The products give identical values at every point: there ",
               "is no exposure to the choice of product here. Check that ",
               "they are really different products."))
  } else {
    poor <- e$variable[e$grade %in% c("poor", "moderate")]
    good <- e$variable[e$grade %in% c("substantial", "almost perfect")]
    if (length(poor)) {
      say(paste0("Models are exposed to the choice of product through ",
                 length(poor), " variable(s) whose values disagree ", where,
                 ": ", paste0(poor, " (", e$worst_source[match(poor,
                   e$variable)], ", CCC ", formatC(e$ccc[match(poor,
                   e$variable)], format = "f", digits = digits), ")",
                   collapse = ", "), "."))
      say(paste0("Fit the model with each product and compare, prefer the ",
                 "product validated for the region, or drop these ",
                 "variables if the others carry the same information."))
    } else {
      say(paste0("The products agree well ", where, " for every variable ",
                 "(CCC at least 0.95): the choice of product should matter ",
                 "little for these variables."))
    }
    if (length(good) && length(poor)) {
      say(paste0("Good agreement: ", paste(good, collapse = ", "), "."))
    }
    off <- e[abs(e$offset_sd) > 0.25 & !e$unit_mismatch, , drop = FALSE]
    if (nrow(off)) {
      say(paste0("Systematic offsets: ", paste0(off$variable, " (",
                 sprintf("%+.2f", off$offset_sd), " SD in ",
                 off$worst_source, ")", collapse = ", "),
                 ". A model fitted with one product and projected with ",
                 "another shifts by this much."))
    }
    mM <- m[m$domain == "M" & !is.na(m$large) & m$large > 0.05, ,
            drop = FALSE]
    if (nrow(mM)) {
      say(paste0("Part of M differs by more than 0.5 SD: ",
                 paste0(mM$variable, " in ", mM$source, " (",
                        share_label(mM$large), " of cells)",
                        collapse = ", "),
                 ". Where these cells lie is shown by ",
                 "plot_climate_sources(x, \"map\")."))
    }
    if (any(e$unit_mismatch)) {
      say(paste0("Offsets above 3 SD in ",
                 paste(e$variable[e$unit_mismatch], collapse = ", "),
                 ": probably different units or scale (degrees C vs kelvin, ",
                 "values x10). Convert before comparing."))
    }
    co <- x$collinearity
    lowc <- co[!is.na(co$mantel_r) & co$mantel_r < 0.8, , drop = FALSE]
    if (nrow(lowc)) {
      say(paste0("The correlations among the variables differ in ",
                 paste0(lowc$source, " (Mantel r ",
                        formatC(lowc$mantel_r, format = "f",
                                digits = digits), ")", collapse = ", "),
                 ": a variable set chosen for low collinearity in one ",
                 "product may not hold in the other."))
    }
    g <- x$geography
    hot <- g[!is.na(g$records_ratio) & g$records_ratio > 1.5, , drop = FALSE]
    if (nrow(hot)) {
      say(paste0("The records sit where the products disagree more than ",
                 "on average across M for ",
                 paste0(hot$variable, " in ", hot$source, " (x",
                        formatC(hot$records_ratio, format = "f", digits = 1),
                        ")", collapse = ", "), "."))
    }
    gr <- g[pmax(abs(g$rho_lon), abs(g$rho_lat), na.rm = TRUE) > 0.5 &
              !is.na(pmax(abs(g$rho_lon), abs(g$rho_lat), na.rm = TRUE)), ,
            drop = FALSE]
    if (nrow(gr)) {
      say(paste0("The disagreement follows a geographic gradient for ",
                 paste0(gr$variable, " in ", gr$source, collapse = ", "),
                 ": see plot_climate_sources(x, \"map\")."))
    }
  }
  cat("\n")
  show <- e
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  if (any(!is.na(x$collinearity$mantel_r))) {
    cat("\nCollinearity structure across M (Mantel r, permutation p):\n")
    cl <- x$collinearity
    cl$mantel_r <- round(cl$mantel_r, digits)
    cl$p <- signif(cl$p, 2)
    print(cl, row.names = FALSE)
  }
  cat("", strwrap(paste0(
    "CCC = Lin's concordance correlation of the product with ", s$reference,
    " ", where, " (1 = identical values; it falls with scatter, offset and ",
    "change of scale); grade after McBride (2005): poor below 0.90, ",
    "moderate 0.90-0.95, substantial 0.95-0.99, almost perfect above 0.99, ",
    "limits from the literature. offset_sd = mean difference (product - ",
    "reference) in standard deviations of the reference across M; offsets ",
    "above 0.25 SD, differences above 0.5 SD in more than 5% of M, Mantel ",
    "r below 0.8 and a records-to-M disagreement ratio above 1.5 are rules ",
    "of thumb. All metrics, per domain, in ",
    "x$metrics; the geographic pattern in x$geography."), width = 76),
    sep = "\n")
  invisible(x)
}

#' Plot a comparison of climate products
#'
#' Four views of a [climate_source_compare()] result:
#' * `"summary"`: per variable and product, Lin's concordance (CCC) with
#'   the reference, the Pearson correlation and the standardised offset (mean
#'   difference +/- its standard deviation, in SD of the reference across
#'   M), at the records and across M. Dashed lines mark the CCC grades of
#'   McBride (2005).
#' * `"map"`: the geographic pattern of disagreement for one variable: the
#'   standardised difference (product minus reference) at the sampled cells
#'   of M, on a diverging scale, with the records outlined and the 5% of
#'   cells with the largest disagreement ringed.
#' * `"scatter"`: product against reference, with the 1:1 line.
#' * `"density"`: the distributions of the variable in each product, at
#'   the records and across M (their shared area is the overlapping index).
#'
#' @param x A `climate_source_compare` object.
#' @param type `"summary"`, `"map"`, `"scatter"` or `"density"`.
#' @param variable Variable for the map, scatter and density views
#'   (default: the one with the lowest agreement).
#' @param source Products shown in the map and scatter views (default: all
#'   but the reference).
#' @param base_size Base font size.
#' @param file Optional file name(s) to save the figure (see
#'   [save_figure()]).
#' @param width,height Figure size in inches when saving.
#' @return A ggplot object.
#' @seealso [climate_source_compare()]
#' @references Lin, L.I-K. (1989) A concordance correlation coefficient to
#'   evaluate reproducibility. *Biometrics* 45, 255-268.
#'   \doi{10.2307/2532051}
#'
#'   McBride, G.B. (2005) A proposal for strength-of-agreement criteria for
#'   Lin's concordance correlation coefficient. NIWA Client Report
#'   HAM2005-062.
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 80)
#' b <- sim$env
#' b[[1]] <- b[[1]] + 0.01 * terra::init(b[[1]], "y")
#' cs <- climate_source_compare(list(A = sim$env, B = b), sim$presence,
#'                              n = 1500)
#' plot_climate_sources(cs)
#' plot_climate_sources(cs, "map")
#' plot_climate_sources(cs, "scatter")
#' plot_climate_sources(cs, "density")
#' @export
plot_climate_sources <- function(x, type = c("summary", "map", "scatter",
                                             "density"),
                                 variable = NULL, source = NULL,
                                 base_size = 11, file = NULL, width = NULL,
                                 height = NULL) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "climate_source_compare")) {
    stop("`x` must come from climate_source_compare().", call. = FALSE)
  }
  type <- match.arg(type)
  s <- x$settings
  ink <- "#2b2b2b"
  muted <- "#6b6b6b"
  others <- setdiff(s$sources, s$reference)
  if (is.null(source)) source <- others
  bad <- setdiff(source, others)
  if (length(bad)) {
    stop("`source` must name product(s) other than the reference: ",
         paste(others, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(variable)) variable <- x$exposure$variable[1]
  if (!variable %in% s$variables) {
    stop("`variable` must be one of: ", paste(s$variables, collapse = ", "),
         ".", call. = FALSE)
  }
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02", "#a6761d", "#666666")
  src_col <- stats::setNames(rep(pal, length.out = length(others)), others)
  dom_lab <- c(
    records = sprintf("records (n = %s)", format(s$n_records, big.mark = ",")),
    M = sprintf("cells of M (n = %s)", format(s$n_cells, big.mark = ",")))
  pts <- x$points
  if (type == "summary") {
    m <- x$metrics
    m$var <- factor(m$variable, levels = rev(x$exposure$variable))
    m$dom <- factor(dom_lab[m$domain], levels = dom_lab)
    pdat <- rbind(
      data.frame(m[c("var", "source", "dom")], panel = "Lin's CCC",
                 value = m$ccc, lo = NA_real_, hi = NA_real_),
      data.frame(m[c("var", "source", "dom")], panel = "Pearson r",
                 value = m$r, lo = NA_real_, hi = NA_real_),
      data.frame(m[c("var", "source", "dom")],
                 panel = "Offset (SD of reference)", value = m$offset_sd,
                 lo = m$offset_sd - sqrt(pmax(m$rmse_sd^2 - m$offset_sd^2,
                                              0)),
                 hi = m$offset_sd + sqrt(pmax(m$rmse_sd^2 - m$offset_sd^2,
                                              0))))
    pdat$panel <- factor(pdat$panel, levels = c("Lin's CCC", "Pearson r",
                                                "Offset (SD of reference)"))
    ref_lines <- rbind(
      data.frame(panel = "Lin's CCC", x = c(0.90, 0.95, 0.99)),
      data.frame(panel = "Offset (SD of reference)", x = 0))
    ref_lines$panel <- factor(ref_lines$panel, levels = levels(pdat$panel))
    pd <- ggplot2::position_dodge(width = 0.6)
    p <- ggplot2::ggplot(pdat, ggplot2::aes(x = .data$value, y = .data$var,
                                            colour = .data$source,
                                            shape = .data$dom,
                                            group = interaction(.data$source,
                                                                .data$dom))) +
      ggplot2::geom_vline(data = ref_lines,
                          ggplot2::aes(xintercept = .data$x),
                          linetype = "dashed", colour = "#9e9e9e",
                          linewidth = 0.4) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi),
                              height = 0, colour = "#7f7f7f",
                              linewidth = 0.5, position = pd,
                              na.rm = TRUE) +
      ggplot2::geom_point(size = 2.4, position = pd, na.rm = TRUE) +
      ggplot2::facet_wrap(~panel, nrow = 1, scales = "free_x") +
      ggplot2::scale_colour_manual(values = src_col, name = "Product") +
      ggplot2::scale_shape_manual(values = c(16, 1), name = NULL) +
      ggplot2::labs(
        x = NULL, y = NULL,
        title = paste0("Agreement of climate products with ", s$reference),
        subtitle = paste0("Dashed: CCC grades of McBride (2005) at 0.90, ",
                          "0.95, 0.99; offset bars: +/- SD of the\n",
                          "standardised difference")) +
      theme_enm(base_size, ink, muted) +
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                     panel.spacing.x = ggplot2::unit(1.5, "lines"))
    return(p)
  }
  ref_col <- paste0(variable, "__", s$reference)
  if (type == "map") {
    sdm <- stats::sd(pts[[ref_col]][pts$domain == "M"], na.rm = TRUE)
    if (!is.finite(sdm) || sdm == 0) sdm <- 1
    long <- do.call(rbind, lapply(source, function(k) {
      d <- (pts[[paste0(variable, "__", k)]] - pts[[ref_col]]) / sdm
      data.frame(lon = pts$lon, lat = pts$lat, domain = pts$domain,
                 source = k, d = d)
    }))
    long <- long[is.finite(long$d), , drop = FALSE]
    lim <- max(stats::quantile(abs(long$d), 0.99, na.rm = TRUE), 1e-9)
    long$dc <- pmax(-lim, pmin(lim, long$d))
    mm <- long[long$domain == "M", , drop = FALSE]
    hot <- do.call(rbind, lapply(split(mm, mm$source), function(z) {
      z[abs(z$d) >= stats::quantile(abs(z$d), 0.95), , drop = FALSE]
    }))
    rec <- long[long$domain == "records", , drop = FALSE]
    mark <- c(sprintf("records (n = %s)", format(s$n_records,
                                                  big.mark = ",")),
              "largest 5% of differences in M")
    p <- ggplot2::ggplot(mm, ggplot2::aes(x = .data$lon, y = .data$lat)) +
      ggplot2::geom_point(ggplot2::aes(colour = .data$dc), size = 1.1) +
      ggplot2::scale_colour_distiller(
        palette = "RdBu", limits = c(-lim, lim), direction = -1,
        name = paste0("(product - ", s$reference, ")\n/ SD of ",
                      s$reference)) +
      ggplot2::facet_wrap(~source) +
      ggplot2::coord_quickmap() +
      ggplot2::labs(
        x = "Longitude", y = "Latitude",
        title = paste0("Disagreement in ", variable, " with ", s$reference),
        subtitle = sprintf(paste0("%s cells of M; colours capped at the ",
                                  "99th percentile (%.2f SD)"),
                           format(s$n_cells, big.mark = ","), lim)) +
      theme_enm(base_size, ink, muted)
    if (!is.null(hot) && nrow(hot)) {
      p <- p + ggplot2::geom_point(
        data = hot, ggplot2::aes(shape = mark[2]), colour = "#1a1a1a",
        size = 2.6, stroke = 0.5)
    }
    if (nrow(rec)) {
      p <- p + ggplot2::geom_point(
        data = rec, ggplot2::aes(fill = .data$dc, shape = mark[1]),
        colour = "#1a1a1a", size = 2.4, stroke = 0.6) +
        ggplot2::scale_fill_distiller(palette = "RdBu",
                                      limits = c(-lim, lim), direction = -1,
                                      guide = "none")
    }
    p <- p + ggplot2::scale_shape_manual(
      values = stats::setNames(c(24, 1), mark), name = NULL,
      guide = ggplot2::guide_legend(override.aes = list(
        fill = "#bdbdbd", size = 2.6)))
    return(p)
  }
  if (type == "scatter") {
    long <- do.call(rbind, lapply(source, function(k) {
      mt <- x$metrics[x$metrics$variable == variable &
                        x$metrics$source == k &
                        x$metrics$domain == s$domain, ]
      data.frame(ref = pts[[ref_col]], val = pts[[paste0(variable, "__", k)]],
                 dom = factor(dom_lab[pts$domain], levels = rev(dom_lab)),
                 source = sprintf("%s: CCC %.2f, r %.2f", k, mt$ccc, mt$r))
    }))
    long <- long[stats::complete.cases(long), , drop = FALSE]
    long <- long[order(long$dom), , drop = FALSE]
    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$ref, y = .data$val,
                                            colour = .data$dom)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, colour = "#9e9e9e",
                           linetype = "dashed") +
      ggplot2::geom_point(size = 1, alpha = 0.6) +
      ggplot2::facet_wrap(~source) +
      ggplot2::scale_colour_manual(
        values = stats::setNames(c("#bdbdbd", "#d95f02"),
                                 dom_lab[c("M", "records")]),
        name = NULL) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        x = paste0(variable, " in ", s$reference),
        y = paste0(variable, " in the product"),
        title = paste0(variable, ": products against ", s$reference),
        subtitle = paste0("Dashed: 1:1 line. CCC and r ",
                          if (s$domain == "records") "at the records" else
                            "across M")) +
      theme_enm(base_size, ink, muted) +
      ggplot2::theme(legend.position = "bottom")
    return(p)
  }
  # density
  long <- do.call(rbind, lapply(s$sources, function(k) {
    data.frame(val = pts[[paste0(variable, "__", k)]], source = k,
               dom = factor(dom_lab[pts$domain], levels = dom_lab))
  }))
  long <- long[is.finite(long$val), , drop = FALSE]
  all_col <- c(stats::setNames("#252525", s$reference), src_col)
  ov <- x$metrics[x$metrics$variable == variable, ]
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$val,
                                          colour = .data$source)) +
    ggplot2::geom_density(linewidth = 0.7) +
    ggplot2::facet_wrap(~dom, scales = "free_y") +
    ggplot2::scale_colour_manual(values = all_col, name = "Product") +
    ggplot2::labs(
      x = variable, y = "Density",
      title = paste0("Distribution of ", variable, " in each product"),
      subtitle = paste0("Overlapping index with ", s$reference, ": ",
                        paste0(ov$source, " ", ov$domain, " ",
                               formatC(ov$overlap, format = "f", digits = 2),
                               collapse = ", "))) +
    theme_enm(base_size, ink, muted) +
    ggplot2::theme(legend.position = "bottom")
  p
}

#' Write a comparison of climate products to CSV files
#'
#' Writes the tables of a [climate_source_compare()] result:
#' `<prefix>_metrics.csv`, `<prefix>_exposure.csv`,
#' `<prefix>_collinearity.csv`, `<prefix>_geography.csv` and
#' `<prefix>_points.csv` (the values of every product at every point, for
#' mapping in a GIS).
#'
#' @param x A `climate_source_compare` object.
#' @param dir Folder (created if needed).
#' @param prefix File name prefix.
#' @return The paths of the files, invisibly.
#' @seealso [climate_source_compare()]
#' @references Rodriguez-Rey, M., Jimenez-Valverde, A. (2024) Differing
#'   sensitivity of species distribution modelling algorithms to climate
#'   data source. *Ecological Informatics* 79, 102387.
#'   \doi{10.1016/j.ecoinf.2023.102387}
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 50)
#' b <- sim$env * 1.05
#' cs <- climate_source_compare(list(A = sim$env, B = b), sim$presence,
#'                              n = 500, n_perm = 0)
#' f <- write_climate_sources(cs, tempdir())
#' basename(f)
#' @export
write_climate_sources <- function(x, dir, prefix = "climate_sources") {
  if (!inherits(x, "climate_source_compare")) {
    stop("`x` must come from climate_source_compare().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  parts <- c("metrics", "exposure", "collinearity", "geography", "points")
  f <- file.path(dir, paste0(prefix, "_", parts, ".csv"))
  for (k in seq_along(parts)) {
    utils::write.csv(x[[parts[k]]], f[k], row.names = FALSE)
  }
  invisible(f)
}
