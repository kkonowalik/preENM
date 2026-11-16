# Are the predictors more informative than random variables? ---------------------

#' Does the species respond to each predictor, or could it be chance?
#'
#' A test to run before modelling. A set of uncorrelated predictors can
#' still contain layers that mean nothing to the species (collinearity
#' screening looks only at the predictors, not at the species), and models
#' will use them anyway: records are few and spatially clustered, so any
#' smooth layer separates them from the background to some degree, and
#' gains some "importance". This function measures, for every layer, how
#' much the conditions at the records differ from those available in the
#' calibration area, and compares it with the same measure for random
#' layers that have the spatial structure of real ones but no relation to
#' the species.
#'
#' **Statistic.** The Kolmogorov-Smirnov distance `D` between the values at
#' the records and the values in the calibration area: the largest
#' difference between the two cumulative distributions (0 = the species
#' uses the layer as it is available, 1 = it uses values the area hardly
#' offers). It responds to a shift (the species prefers high or low
#' values) and to a narrowing (it avoids both extremes).
#'
#' **Null model.** `D` is recomputed `n_null` times at the same record
#' cells for a random layer:
#'
#' * `"field"` (default): a smooth random surface with the same range of
#'   spatial autocorrelation as the layer (estimated as in
#'   [suggest_block_size()]). Works for a calibration area of any shape.
#' * `"shift"`: the layer itself, moved by a random distance in both
#'   directions (wrapping around the edges) and randomly mirrored, which
#'   keeps its exact spatial structure (the torus-translation test of Harms
#'   et al. 2001). Needs a calibration area that fills most of the raster's
#'   rectangle: records that land on cells without values are left out, and
#'   a shift is redrawn if fewer than half remain.
#'
#' Random values without spatial structure (white noise, as in
#' [random_variables()] with `method = "shuffle"`) would be a much weaker
#' null: real layers are smooth and records are clustered, so even a
#' meaningless smooth layer, such as a painting, differs from the
#' background at the records far more than noise does (Fourcade et al.
#' 2018).
#'
#' **Reading the result.** `p_value` is the share of random layers with a
#' `D` at least as large as the real one; `relevant` marks layers with
#' `p_value <= alpha`. Layers that are not distinguishable from random ones
#' are candidates for removal. The test is a screen, with three limits: it
#' looks at one layer at a time, so a layer that matters only in
#' combination with another can be missed; a layer correlated with a true
#' driver is relevant too; and recording bias along a layer (e.g. records
#' mostly from lowlands) looks like a response, so use thinned records.
#' The test is also deliberately cautious about layers that are a single
#' gradient across the whole area (e.g. temperature in a small flat
#' region): records of a species with a compact range differ from the
#' background on any such gradient, so the layer cannot be told apart from
#' geography. Then keep it on ecological grounds, or enlarge the area. To
#' test the layers inside your models, see [random_variables()] and
#' [compare_to_random()].
#'
#' @param points Records: a data.frame, CSV path, `sf` or
#'   [terra::SpatVector].
#' @param env Predictor layers of the calibration area
#'   ([terra::SpatRaster] or file paths).
#' @param null The null model: `"field"` or `"shift"`.
#' @param n_null Number of random layers per predictor.
#' @param alpha Largest `p_value` of a relevant layer.
#' @param n_bg Number of cells sampled to describe the calibration area.
#' @param min_shift With `null = "shift"`: smallest shift, as a share of
#'   the raster's rows and columns.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Seed.
#' @return A `variable_relevance` object: a data.frame with one row per
#'   layer, ordered from the most to the least distinct: `D`, `null_mean`,
#'   `null_95` (the 95th percentile of the null), `p_value`, `relevant`,
#'   and two descriptions of the response, `shift_sd` (mean at the records
#'   minus mean of the area, in standard deviations of the area) and
#'   `breadth_ratio` (standard deviation at the records / in the area).
#' @references Harms, K.E., Condit, R., Hubbell, S.P., Foster, R.B. (2001)
#'   Habitat associations of trees and shrubs in a 50-ha neotropical forest
#'   plot. *Journal of Ecology* 89, 947-959. Fourcade, Y., Besnard, A.G.,
#'   Secondi, J. (2018) Paintings predict the distribution of species, or
#'   the challenge of selecting environmental predictors and evaluation
#'   statistics. *Global Ecology and Biogeography* 27, 245-256. Konowalik,
#'   K. (2022) Phylogeography and colonization pattern of subendemic
#'   round-leaved oxeye daisy from the Dinarides to the Carpathians.
#'   *Scientific Reports* 12, 16443.
#' @seealso [random_variables()], [compare_to_random()]
#' @export
#' @examples
#' # three patchy layers; the species lives where the first one is high
#' sim <- simulate_enm_data(seed = 1)
#' env <- random_variables(sim$env, n = 3, range_km = 100, seed = 3)
#' names(env) <- c("temperature", "rainfall", "soil")
#' set.seed(1)
#' cells <- sample(terra::ncell(env), 150,
#'                 prob = terra::values(env$temperature)[, 1]^6)
#' records <- as.data.frame(terra::xyFromCell(env, cells))
#' names(records) <- c("lon", "lat")
#' variable_relevance(records, env, n_null = 199)
variable_relevance <- function(points, env, null = c("field", "shift"),
                               n_null = 499, alpha = 0.05, n_bg = 10000,
                               min_shift = 0.1, coords = NULL,
                               points_crs = "EPSG:4326", seed = 1) {
  null <- match.arg(null)
  r <- as_grid(env)
  vars <- names(r)
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
  cell <- unique(terra::cellFromXY(r, terra::crds(p)))
  E <- env_at_cells(r, cell)
  keep <- stats::complete.cases(E)
  cell <- cell[keep]
  E <- E[keep, , drop = FALSE]
  if (length(cell) < 10) {
    stop("Fewer than 10 occupied cells with values in all layers.",
         call. = FALSE)
  }
  valid <- valid_cells(r)
  bg <- with_seed(seed, sort(valid[sample.int(length(valid),
                                              min(length(valid), n_bg))]))
  M <- read_cells(r, bg)
  colnames(M) <- vars

  observed <- vapply(vars, function(v) ks_distance(E[, v], M[, v]),
                     numeric(1))
  nulls <- with_seed(seed, if (null == "field") {
    ranges <- suggest_block_size(r, seed = seed)$ranges
    rg <- stats::setNames(ranges$range_km, ranges$layer)
    rg[is.na(rg)] <- stats::median(rg, na.rm = TRUE)
    xy_rec <- terra::xyFromCell(r, cell)
    xy_bg <- terra::xyFromCell(r, bg)
    vapply(vars, function(v) {
      vapply(seq_len(n_null), function(b) {
        z <- noise_at(r, rg[[v]], rbind(xy_rec, xy_bg))
        i <- seq_len(nrow(xy_rec))
        ks_distance(z[i], z[-i])
      }, numeric(1))
    }, numeric(n_null))
  } else {
    shift_null(r, cell, M, n_null, min_shift)
  })
  nulls <- matrix(nulls, nrow = n_null, dimnames = list(NULL, vars))

  out <- data.frame(
    layer = vars, D = observed, null_mean = colMeans(nulls, na.rm = TRUE),
    null_95 = apply(nulls, 2, stats::quantile, 0.95, na.rm = TRUE,
                    names = FALSE),
    p_value = vapply(vars, function(v) {
      (1 + sum(nulls[, v] >= observed[[v]], na.rm = TRUE)) /
        (1 + sum(!is.na(nulls[, v])))
    }, numeric(1)),
    shift_sd = (colMeans(E) - colMeans(M)) / apply(M, 2, stats::sd),
    breadth_ratio = apply(E, 2, stats::sd) / apply(M, 2, stats::sd),
    row.names = NULL)
  out$relevant <- out$p_value <= alpha
  out <- out[order(out$p_value, -out$D), c("layer", "D", "null_mean",
                                           "null_95", "p_value", "relevant",
                                           "shift_sd", "breadth_ratio")]
  rownames(out) <- NULL
  attr(out, "settings") <- list(null = null, n_null = n_null, alpha = alpha,
                                n_records = length(cell))
  class(out) <- c("variable_relevance", "data.frame")
  out
}

#' @export
print.variable_relevance <- function(x, digits = 3, ...) {
  st <- attr(x, "settings")
  df <- as.data.frame(unclass(x), check.names = FALSE)
  if (!is.null(st)) {
    cat("Relevance of the predictors: ", st$n_records, " occupied cells, ",
        st$n_null, " random layers each (null: ", st$null, ")\n", sep = "")
    cat("  Not distinguishable from random (p > ", st$alpha, "): ",
        if (any(!df$relevant)) paste(df$layer[!df$relevant], collapse = ", ")
        else "none", "\n\n", sep = "")
  }
  num <- vapply(df, is.double, logical(1))
  df[num] <- lapply(df[num], signif, digits = digits)
  print(df, row.names = FALSE, ...)
  invisible(x)
}

# Largest difference between the cumulative distributions of x and of the
# reference sample (Kolmogorov-Smirnov distance) (internal)
ks_distance <- function(x, ref) {
  x <- sort(x[!is.na(x)])
  ref <- sort(ref[!is.na(ref)])
  n <- length(x)
  if (!n || !length(ref)) return(NA_real_)
  at <- unique(x)
  fx <- findInterval(at, x) / n
  fr <- findInterval(at, ref) / length(ref)
  # just below each value, where the step of x has not happened yet
  fx0 <- findInterval(at, x, left.open = TRUE) / n
  fr0 <- findInterval(at, ref, left.open = TRUE) / length(ref)
  max(abs(fx - fr), abs(fx0 - fr0))
}

# Values at `xy` of a smooth random surface with about the given range of
# autocorrelation: white noise on a coarse grid of that spacing, read by
# bilinear interpolation; three grids with random origins are summed so
# that the lattice does not show (internal)
noise_at <- function(r, range_km, xy) {
  e <- terra::ext(r)
  if (terra::is.lonlat(r)) {
    mid <- (e$ymin + e$ymax) / 2
    hx <- range_km / (111.32 * max(cos(mid * pi / 180), 0.05))
    hy <- range_km / 111.32
  } else {
    hx <- hy <- range_km * 1000
  }
  # at least a few grid cells, at most a few hundred along each side
  rs <- terra::res(r)
  hx <- min(max(hx, 2 * rs[1]), (e$xmax - e$xmin))
  hy <- min(max(hy, 2 * rs[2]), (e$ymax - e$ymin))
  z <- numeric(nrow(xy))
  nx <- ceiling((e$xmax - e$xmin) / hx) + 3
  ny <- ceiling((e$ymax - e$ymin) / hy) + 3
  for (k in 1:3) {
    # position in units of the coarse grid, with a random origin
    u <- (xy[, 1] - e$xmin) / hx + stats::runif(1)
    v <- (xy[, 2] - e$ymin) / hy + stats::runif(1)
    i <- pmin(pmax(floor(u), 0), nx - 2)
    j <- pmin(pmax(floor(v), 0), ny - 2)
    fu <- u - i
    fv <- v - j
    g <- stats::rnorm(nx * ny)
    at <- i + j * nx + 1
    z <- z + g[at] * (1 - fu) * (1 - fv) + g[at + 1] * fu * (1 - fv) +
      g[at + nx] * (1 - fu) * fv + g[at + nx + 1] * fu * fv
  }
  z
}

# Null distances from the layers moved around the torus (internal)
shift_null <- function(r, cell, M, n_null, min_shift) {
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  rc <- terra::rowColFromCell(r, cell)
  vars <- names(r)
  out <- matrix(NA_real_, n_null, length(vars))
  pick <- function(n) {
    lo <- ceiling(min_shift * n)
    if (n - lo <= lo) return(sample.int(n, 1) - 1L)
    sample(lo:(n - lo), 1)
  }
  failed <- 0L
  for (b in seq_len(n_null)) {
    for (try in 1:20) {
      row <- (rc[, 1] - 1 + pick(nr)) %% nr + 1
      col <- (rc[, 2] - 1 + pick(nc)) %% nc + 1
      if (stats::runif(1) < 0.5) row <- nr + 1 - row
      if (stats::runif(1) < 0.5) col <- nc + 1 - col
      X <- env_at_cells(r, terra::cellFromRowCol(r, row, col))
      ok <- stats::complete.cases(X)
      if (sum(ok) >= length(cell) / 2) break
    }
    if (sum(ok) < length(cell) / 2) {
      failed <- failed + 1L
      next
    }
    out[b, ] <- vapply(seq_along(vars), function(j) {
      ks_distance(X[ok, j], M[, j])
    }, numeric(1))
  }
  if (failed > n_null / 2) {
    stop("Most shifts move the records onto cells without values: the ",
         "calibration area does not fill the raster. Use null = \"field\".",
         call. = FALSE)
  }
  if (failed) warning(failed, " of ", n_null, " shifts left fewer than half ",
                      "of the records on cells with values and were ",
                      "dropped.", call. = FALSE)
  out
}

#' Random variables with the values and the spatial structure of real ones
#'
#' Makes layers that mean nothing to the species, to be modelled together
#' with the real predictors in any software: a real predictor that gets no
#' more importance than these is not informative (see
#' [compare_to_random()]). Each random layer takes the values of one real
#' layer, so it has exactly the same distribution of values (range, mean,
#' skewness, ties) and the same cells without values, and arranges them
#' anew:
#'
#' * `"field"` (default): in smooth patches with the same range of spatial
#'   autocorrelation as the real layer (estimated as in
#'   [suggest_block_size()]). This is the fair comparison: real layers are
#'   smooth and records are clustered, so a smooth layer gains importance
#'   by chance far more easily than noise does.
#' * `"shuffle"`: cell by cell at random (white noise), the classic random
#'   variable. A real predictor will almost always beat it; use it only to
#'   compare with earlier work.
#'
#' @param env Predictor layers ([terra::SpatRaster] or file paths).
#' @param n Number of random layers. They are made from the real layers in
#'   turn (the first from the first layer, and so on, starting again after
#'   the last), or from those named in `layers`.
#' @param method `"field"` or `"shuffle"`.
#' @param layers Optional names of the real layers to imitate.
#' @param range_km With `method = "field"`: size of the patches (range of
#'   spatial autocorrelation, km). `NULL` = that of each imitated layer;
#'   one number = the same for all.
#' @param seed Seed.
#' @param filename Optional output file.
#' @param overwrite Overwrite `filename`.
#' @return A [terra::SpatRaster] with layers named `random_1`, `random_2`,
#'   ...; `attr(, "source")` says which real layer each one imitates.
#' @seealso [variable_relevance()], [compare_to_random()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' rnd <- random_variables(sim$env, n = 2, seed = 1)
#' terra::plot(c(sim$env, rnd))
#' # same values as the real layers, arranged differently
#' range(terra::values(rnd[[1]])) == range(terra::values(sim$env[[1]]))
random_variables <- function(env, n = 5, method = c("field", "shuffle"),
                             layers = NULL, range_km = NULL, seed = NULL,
                             filename = "", overwrite = FALSE) {
  method <- match.arg(method)
  r <- as_grid(env)
  if (is.null(layers)) layers <- names(r)
  miss <- setdiff(layers, names(r))
  if (length(miss)) stop("Layer(s) not found: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  src <- rep_len(layers, n)
  out <- with_seed(seed, {
    rg <- if (method == "field" && !is.null(range_km)) {
      if (!is.numeric(range_km) || length(range_km) != 1 || is.na(range_km) ||
          range_km <= 0) {
        stop("`range_km` must be one positive number.", call. = FALSE)
      }
      stats::setNames(rep(range_km, length(unique(src))), unique(src))
    } else if (method == "field") {
      x <- suggest_block_size(r[[unique(src)]],
                              seed = if (is.null(seed)) 1 else seed)$ranges
      z <- stats::setNames(x$range_km, x$layer)
      z[is.na(z)] <- stats::median(z, na.rm = TRUE)
      z
    }
    terra::rast(lapply(seq_len(n), function(i) {
      real <- r[[src[i]]]
      v <- terra::values(real, mat = FALSE)
      ok <- which(!is.na(v))
      key <- if (method == "field") {
        noise_at(real, rg[[src[i]]], terra::xyFromCell(real, ok))
      } else stats::runif(length(ok))
      # the real values, placed in the order of the random surface
      v[ok][order(key)] <- sort(v[ok])
      res <- terra::rast(real)
      terra::values(res) <- v
      res
    }))
  })
  names(out) <- paste0("random_", seq_len(n))
  if (nzchar(filename)) {
    out <- terra::writeRaster(out, filename, overwrite = overwrite)
  }
  attr(out, "source") <- stats::setNames(src, names(out))
  out
}

#' Which predictors are more important than random variables in the models?
#'
#' Summarises variable importance from models fitted with the real
#' predictors together with random ones (see [random_variables()]), in any
#' software and for any number of algorithms and runs. In each run, a real
#' predictor *beats random* when its importance is above that of the random
#' variables. A predictor is kept unless it fails to beat random in most
#' runs of most algorithms, the consensus rule of Konowalik (2022):
#' algorithms measure importance differently, and a predictor that only one
#' of them finds useless may still matter.
#'
#' After removing the predictors that are not kept, fit the models again
#' (with new random variables) and repeat until all remaining predictors
#' are kept: importance is shared among predictors, so it changes when one
#' is removed.
#'
#' @param importance A data.frame in long format: one row per variable,
#'   algorithm and run.
#' @param variable,value Columns with the variable name and its importance
#'   (any measure where larger means more important: contribution,
#'   permutation importance, gain, ...).
#' @param algorithm,run Optional columns with the algorithm and the run
#'   (replicate, fold or background set). Without them all rows are one
#'   run of one algorithm.
#' @param random Regular expression that recognises the random variables by
#'   name.
#' @param threshold What a predictor must exceed in a run: `"max"`, the
#'   most important random variable (strict; the more random variables,
#'   the stricter), or a number between 0 and 1, that quantile of the
#'   random variables (e.g. `0.5` for their median).
#' @param min_runs Share of the runs of an algorithm in which a predictor
#'   must beat random to count as supported by that algorithm.
#' @return A data.frame with one row per real predictor: the share of runs
#'   in which it beats random, one column per algorithm; `algorithms` (how
#'   many support it); `mean_importance`; and `keep` (`FALSE` when most
#'   algorithms do not support it), ordered from the best supported.
#'   `attr(, "random")` holds the threshold of every run.
#' @references Konowalik, K. (2022) Phylogeography and colonization pattern
#'   of subendemic round-leaved oxeye daisy from the Dinarides to the
#'   Carpathians. *Scientific Reports* 12, 16443.
#' @seealso [random_variables()], [variable_relevance()]
#' @export
#' @examples
#' imp <- data.frame(
#'   algorithm = rep(c("maxent", "rf", "gbm"), each = 5),
#'   variable = rep(c("bio1", "bio12", "wind", "random_1", "random_2"), 3),
#'   importance = c(52, 30, 4, 8, 6,   40, 35, 12, 9, 4,   61, 25, 2, 7, 5))
#' compare_to_random(imp)
compare_to_random <- function(importance, variable = "variable",
                              value = "importance", algorithm = NULL,
                              run = NULL, random = "^random", threshold = "max",
                              min_runs = 0.5) {
  d <- as.data.frame(importance, stringsAsFactors = FALSE)
  if (is.null(algorithm) && "algorithm" %in% names(d)) algorithm <- "algorithm"
  if (is.null(run) && "run" %in% names(d)) run <- "run"
  miss <- setdiff(c(variable, value, algorithm, run), names(d))
  if (length(miss)) {
    stop("Column(s) not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  if (!identical(threshold, "max") &&
      !(is.numeric(threshold) && length(threshold) == 1 &&
        threshold >= 0 && threshold <= 1)) {
    stop("`threshold` must be \"max\" or a quantile between 0 and 1.",
         call. = FALSE)
  }
  v <- as.character(d[[variable]])
  imp <- as.numeric(d[[value]])
  alg <- if (is.null(algorithm)) rep("all", nrow(d)) else
    as.character(d[[algorithm]])
  rn <- if (is.null(run)) rep("1", nrow(d)) else as.character(d[[run]])
  is_rnd <- grepl(random, v)
  if (!any(is_rnd)) {
    stop("No random variables found (names matching \"", random, "\").",
         call. = FALSE)
  }
  if (all(is_rnd)) stop("Only random variables found.", call. = FALSE)
  key <- paste(alg, rn, sep = "\r")
  limit <- tapply(imp[is_rnd], key[is_rnd], function(z) {
    z <- z[!is.na(z)]
    if (!length(z)) NA_real_ else if (identical(threshold, "max")) max(z) else
      stats::quantile(z, threshold, names = FALSE)
  })
  no_rnd <- setdiff(unique(key), names(limit)[!is.na(limit)])
  if (length(no_rnd)) {
    stop(length(no_rnd), " run(s) have no random variable (e.g. ",
         sub("\r", ", run ", no_rnd[1]), ").", call. = FALSE)
  }
  real <- !is_rnd & !is.na(imp)
  beats <- imp[real] > limit[key[real]]
  algs <- unique(alg)
  vars <- unique(v[real])
  share <- vapply(algs, function(a) {
    vapply(vars, function(x) {
      i <- v[real] == x & alg[real] == a
      if (any(i)) mean(beats[i]) else NA_real_
    }, numeric(1))
  }, numeric(length(vars)))
  share <- matrix(share, nrow = length(vars), dimnames = list(vars, algs))
  supported <- share >= min_runs
  n_alg <- rowSums(!is.na(supported))
  out <- data.frame(variable = vars, share, check.names = FALSE,
                    row.names = NULL)
  out$algorithms <- rowSums(supported, na.rm = TRUE)
  out$mean_importance <- as.vector(tapply(imp[real], factor(v[real], vars),
                                          mean))
  # dropped only when most of the algorithms do not support it
  out$keep <- (n_alg - out$algorithms) <= n_alg / 2
  out <- out[order(-out$algorithms, -rowMeans(share, na.rm = TRUE),
                   -out$mean_importance), ]
  rownames(out) <- NULL
  lim <- data.frame(algorithm = sub("\r.*", "", names(limit)),
                    run = sub(".*\r", "", names(limit)),
                    threshold = as.vector(limit))
  attr(out, "random") <- lim
  out
}
