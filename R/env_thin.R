# Environmental thinning of presence records -------------------------------------

#' Thin presence records in environmental space, keeping the niche intact
#'
#' Dense sampling (e.g. the Alps, Scandinavia) puts many records into the same
#' part of environmental space, so models overweight those conditions.
#' `thin_environment()` keeps a set of records that represents the environment
#' of *all* records, without an arbitrary target number, and without
#' truncating the niche environmentally or geographically:
#'
#' 1. **Environmental space.** A PCA is fitted on a random sample of the whole
#'    raster area (the environment available to the species, e.g. M) and the
#'    records are placed in it; the first PCs explaining at least
#'    `pc_cumvar` of the variance are used (`axes = "raw"` uses the
#'    standardised layers instead). Precipitation amounts (`log_vars`) are
#'    log-transformed for this space only, so that differences count in
#'    relative terms.
#' 2. **Bin width from the data.** Every record is accurate to `radius_km`,
#'    so its true position may be any cell within that distance of its cell
#'    (the window spans more columns at high latitudes, where cells are
#'    narrow). The range of each axis within the window is the record's
#'    *environmental uncertainty*; the bin width of each axis is its median
#'    (`width_quantile`) over all records. Two records in one bin differ less
#'    than a single record's own uncertainty: they are environmental
#'    duplicates.
#' 3. **One record per bin** (`max_per_bin`), the one closest to the bin's
#'    centre (`pick = "central"`) or a random one. Every occupied bin keeps a
#'    record, so every removed record has a kept one less than one bin width
#'    away on every axis: only the density is evened out.
#' 4. **Always kept:** records you protect (e.g. genetic samples), the minimum
#'    and maximum of every layer (optionally the outer tails), and the convex
#'    hull vertices of every pair of axes (extreme *combinations* of
#'    conditions).
#' 5. **Geographic quota:** the area is cut into `block_km` equal-area blocks
#'    and every block keeps at least `min(block_min, records it had)`;
#'    missing records are added choosing the least represented environments
#'    first. Small or isolated areas and lone records are kept.
#'
#' Run it on records already reduced to one per raster cell (see
#' [thin_by_cell()]), with the rasters you will model with, clipped to M.
#'
#' @param points Presence records: a data.frame with coordinates (e.g. from
#'   [thin_by_cell()]), a CSV path, an `sf` object or a [terra::SpatVector].
#' @param env Predictor rasters ([terra::SpatRaster] or file paths), clipped
#'   to the accessible area.
#' @param axes `"pca"` (default) or `"raw"` (standardised layers).
#' @param log_vars Layers log-transformed (`log1p`) for the thinning space:
#'   names, or `NULL` to detect bioclim precipitation amounts (bio12, bio13,
#'   bio14, bio16, bio17, bio18, bio19) by name; `character(0)` for none.
#' @param pca_n Random cells of the raster area used to fit the PCA.
#' @param pc_cumvar,pc_min Keep the first PCs explaining at least `pc_cumvar`
#'   of the variance, but at least `pc_min` of them.
#' @param radius_km Positional accuracy of a record (km); the window of a
#'   record is every cell within this distance of its cell. `NA` uses a fixed
#'   window of `window` cells around the record's cell.
#' @param window Window half-width in cells when `radius_km = NA` (1 = 3 x 3).
#' @param width_quantile Quantile of the window ranges used as bin width.
#' @param width_mult Multiplier of the bin widths (sensitivity runs only).
#' @param max_per_bin Records kept per occupied bin (protected ones count).
#' @param pick `"central"` or `"random"`.
#' @param protect Records always kept: a logical vector, or a named list of
#'   values per column, e.g. `list(source = "sNMF")`.
#' @param protect_minmax Keep the minimum and maximum of every layer.
#' @param protect_tail Also keep this share of each layer's lowest and highest
#'   values (e.g. 0.005); 0 = off.
#' @param protect_hull Keep the convex-hull vertices of every pair of axes.
#' @param block_km,block_min Geographic quota: block size (km) and records
#'   kept per block (`block_min = 0` switches it off).
#' @param block_crs Equal-area CRS for the blocks (default: Lambert azimuthal
#'   equal-area centred on the records).
#' @param template Optional earlier `env_thinning` whose PCA, axes and bin
#'   widths are re-used (e.g. to thin lineages or subsets in the same space).
#' @param coords,points_crs Coordinate columns and their CRS when `points` is
#'   a table.
#' @param seed Optional random seed.
#' @param verbose Print progress.
#' @return An object of class `env_thinning`: `data` (all records with their
#'   layer values `env_*`, axis scores `ax_*`, `bin`, `block`, `kept` and
#'   `keep_reason`: protected, extreme_minmax, extreme_tail, hull, env_bin,
#'   block_quota or `NA` = removed), `pca`, `bin_widths`, `report` (checks)
#'   and `settings`. Use [kept()], [plot_env_thinning()] and
#'   [write_env_thinning()].
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 300)
#' th <- thin_environment(sim$presence, sim$env, radius_km = 20,
#'                        block_km = 100, seed = 1)
#' th
#' plot_env_thinning(th)
thin_environment <- function(points, env, axes = c("pca", "raw"),
                             log_vars = NULL, pca_n = 2e5, pc_cumvar = 0.95,
                             pc_min = 2, radius_km = 1, window = 1,
                             width_quantile = 0.5, width_mult = 1,
                             max_per_bin = 1, pick = c("central", "random"),
                             protect = NULL, protect_minmax = TRUE,
                             protect_tail = 0, protect_hull = TRUE,
                             block_km = 50, block_min = 2, block_crs = NULL,
                             template = NULL, coords = NULL,
                             points_crs = "EPSG:4326", seed = NULL,
                             verbose = TRUE) {
  axes <- match.arg(axes)
  pick <- match.arg(pick)
  say <- function(...) if (verbose) message(...)
  if (max_per_bin < 1) stop("`max_per_bin` must be at least 1.", call. = FALSE)
  r <- as_grid(env)
  vars <- names(r)
  geo <- terra::is.lonlat(r)

  # records, their cells and values
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
    lon <- as_number(df[[xy_cols[1]]])
    lat <- as_number(df[[xy_cols[2]]])
    ok <- !is.na(lon) & !is.na(lat)
    df <- df[ok, , drop = FALSE]
    p <- terra::vect(cbind(lon[ok], lat[ok]), crs = points_crs)
  } else {
    p <- as_points(points, coords, points_crs)
    df <- as.data.frame(p)
  }
  ll <- terra::crds(terra::project(p, "EPSG:4326"))
  df$lon <- ll[, 1]
  df$lat <- ll[, 2]
  pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
  cell <- terra::cellFromXY(r, terra::crds(pr))
  E <- env_at_cells(r, cell)
  ok <- stats::complete.cases(E)
  if (any(!ok)) say("Records without values in all layers (dropped): ",
                    sum(!ok))
  prot <- protect_flag(protect, df)
  df <- df[ok, , drop = FALSE]
  E <- E[ok, , drop = FALSE]
  cell <- cell[ok]
  prot <- prot[ok]
  n <- nrow(df)
  if (n < 2) stop("Fewer than 2 records with values in all layers.",
                  call. = FALSE)

  # environmental space
  if (!is.null(template)) {
    if (!inherits(template, "env_thinning")) {
      stop("`template` must come from thin_environment().", call. = FALSE)
    }
    if (!identical(template$settings$vars, vars)) {
      stop("`template` was made with other layers.", call. = FALSE)
    }
    pca <- template$pca
    lv <- template$settings$log_vars
    axes <- template$settings$axes
    n_ax <- template$settings$n_axes
  } else {
    lv <- if (is.null(log_vars)) default_log_vars(vars) else
      intersect(log_vars, vars)
    valid <- valid_cells(r)
    smp <- with_seed(seed, sort(valid[sample.int(length(valid),
                                                 min(length(valid), pca_n))]))
    M <- read_cells(r, smp)
    colnames(M) <- vars
    M <- M[stats::complete.cases(M), , drop = FALSE]
    if (length(lv) && any(M[, lv] < 0)) {
      stop("Negative values in a layer chosen for the log transformation: ",
           paste(lv, collapse = ", "), call. = FALSE)
    }
    pca <- stats::prcomp(log_tf(M, lv), center = TRUE, scale. = TRUE)
    ve <- pca$sdev^2 / sum(pca$sdev^2)
    n_ax <- if (axes == "pca") {
      min(length(vars), max(pc_min, which(cumsum(ve) >= pc_cumvar - 1e-9)[1]))
    } else length(vars)
    say("Environmental space: ", if (axes == "pca") paste0(
      "PC1-PC", n_ax, " (", round(100 * sum(ve[seq_len(n_ax)]), 1),
      "% of the variance of ", nrow(M), " cells)") else
        "the standardised layers",
      if (length(lv)) paste0("; log-transformed: ", paste(lv, collapse = ", ")))
  }
  ax_names <- if (axes == "pca") paste0("PC", seq_len(n_ax)) else vars
  to_axes <- function(X) {
    X <- log_tf(X[, vars, drop = FALSE], lv)
    A <- if (axes == "pca") {
      stats::predict(pca, newdata = X)[, seq_len(n_ax), drop = FALSE]
    } else {
      sweep(sweep(X, 2, pca$center, "-"), 2, pca$scale, "/")
    }
    A <- unname(as.matrix(A))
    colnames(A) <- ax_names
    A
  }
  Z <- to_axes(E)

  # bin widths: environmental uncertainty of one record
  if (!is.null(template)) {
    w <- template$bin_widths * width_mult / template$settings$width_mult
    win_tab <- NULL
  } else {
    win <- window_ranges(r, cell, df$lat, radius_km, window, to_axes, geo)
    w <- apply(win$axes, 2, stats::quantile, probs = width_quantile,
               na.rm = TRUE) * width_mult
    if (any(!is.finite(w) | w <= 0)) {
      stop("A bin width is 0 or not finite: the layers barely vary within ",
           "the window. Increase `radius_km` (or `window`) or ",
           "`width_quantile`.", call. = FALSE)
    }
    names(w) <- ax_names
    win_tab <- data.frame(layer = vars,
                          window_range = apply(win$raw, 2, stats::quantile,
                                               probs = width_quantile,
                                               na.rm = TRUE),
                          range_all_records = apply(E, 2, function(v)
                            diff(range(v))), row.names = NULL)
    say("Bin widths (", if (is.na(radius_km)) paste0((2 * window + 1),
                                                     "x", (2 * window + 1),
                                                     " cells") else
                                                       paste(radius_km, "km"),
        " window): ", paste(ax_names, signif(w, 3), sep = " = ",
                            collapse = ", "))
  }
  lo <- if (!is.null(template)) template$bin_origin else apply(Z, 2, min)
  B <- sweep(sweep(Z, 2, lo, "-"), 2, w, "/")
  BIN <- apply(floor(B), 1, paste, collapse = "_")

  # always kept
  KEEP <- rep(NA_character_, n)
  mark <- function(idx, why) {
    idx <- unique(idx[is.na(KEEP[idx])])
    KEEP[idx] <<- why
  }
  mark(which(prot), "protected")
  if (protect_minmax) {
    mark(unlist(lapply(vars, function(v) c(which.min(E[, v]),
                                           which.max(E[, v])))),
         "extreme_minmax")
  }
  if (protect_tail > 0) {
    mark(unlist(lapply(vars, function(v) {
      q <- stats::quantile(E[, v], c(protect_tail, 1 - protect_tail))
      which(E[, v] <= q[1] | E[, v] >= q[2])
    })), "extreme_tail")
  }
  if (protect_hull && n_ax >= 2 && n >= 3) {
    pairs <- utils::combn(n_ax, 2)
    mark(unlist(lapply(seq_len(ncol(pairs)), function(k) {
      grDevices::chull(Z[, pairs[1, k]], Z[, pairs[2, k]])
    })), "hull")
  }

  # one record per bin
  mu <- apply(B, 2, function(v) stats::ave(v, BIN))
  mu <- matrix(mu, nrow = n)
  score <- if (pick == "central") sqrt(rowSums((B - mu)^2)) else
    with_seed(seed, stats::runif(n))
  prot_in_bin <- stats::ave(as.numeric(!is.na(KEEP)), BIN, FUN = sum)
  need <- pmax(0, max_per_bin - prot_in_bin)
  cand <- which(is.na(KEEP))
  if (length(cand)) {
    o <- cand[order(BIN[cand], score[cand])]
    rk <- stats::ave(seq_along(o), BIN[o], FUN = seq_along)
    mark(o[rk <= need[o]], "env_bin")
  }

  # geographic quota
  if (is.null(block_crs)) {
    block_crs <- sprintf(paste0("+proj=laea +lat_0=%.4f +lon_0=%.4f +x_0=0 ",
                                "+y_0=0 +datum=WGS84 +units=m +no_defs"),
                         mean(df$lat), mean(df$lon))
  }
  P <- terra::crds(terra::project(terra::vect(cbind(df$lon, df$lat),
                                              crs = "EPSG:4326"), block_crs))
  bs <- block_km * 1000
  BLOCK <- paste(floor(P[, 1] / bs), floor(P[, 2] / bs), sep = "_")
  n_blk <- stats::ave(rep(1, n), BLOCK, FUN = sum)
  quota <- pmin(n_blk, block_min)
  added <- 0
  if (block_min > 0) {
    kept_idx <- which(!is.na(KEEP))
    d_nn <- rep(Inf, n)
    if (length(kept_idx) && length(kept_idx) < n) {
      d_nn[-kept_idx] <- nn_dist(B[-kept_idx, , drop = FALSE],
                                 B[kept_idx, , drop = FALSE])
    }
    for (ii in split(seq_len(n), BLOCK)) {
      miss <- quota[ii[1]] - sum(!is.na(KEEP[ii]))
      if (miss <= 0) next
      cc <- ii[is.na(KEEP[ii])]
      d <- d_nn[cc]
      for (step in seq_len(min(miss, length(cc)))) {
        j <- which.max(d)                     # least represented environment
        KEEP[cc[j]] <- "block_quota"
        added <- added + 1
        d <- pmin(d, sqrt(colSums((t(B[cc, , drop = FALSE]) - B[cc[j], ])^2)))
        d[j] <- -Inf
      }
    }
  }
  kept_flag <- !is.na(KEEP)
  say("Kept ", sum(kept_flag), " of ", n, " records (",
      round(100 * mean(kept_flag), 1), "%); occupied bins: ",
      length(unique(BIN)), "; added by the block quota: ", added)

  # checks
  Et <- E[kept_flag, , drop = FALSE]
  ranges <- data.frame(layer = vars, min_all = apply(E, 2, min),
                       min_kept = apply(Et, 2, min),
                       max_all = apply(E, 2, max),
                       max_kept = apply(Et, 2, max),
                       sd_all = apply(E, 2, stats::sd),
                       sd_kept = apply(Et, 2, stats::sd), row.names = NULL)
  intervals <- do.call(rbind, lapply(vars, function(v) {
    if (diff(range(E[, v])) == 0) return(NULL)
    br <- seq(min(E[, v]), max(E[, v]), length.out = 11)
    fa <- cut(E[, v], br, include.lowest = TRUE)
    data.frame(layer = v, interval = levels(fa),
               n_all = as.vector(table(fa)),
               n_kept = as.vector(table(fa[kept_flag])), row.names = NULL)
  }))
  d_rm <- if (any(!kept_flag)) nn_dist(B[!kept_flag, , drop = FALSE],
                                       B[kept_flag, , drop = FALSE]) else
                                         numeric(0)
  w0 <- w / width_mult
  sens <- data.frame(width_mult = c(0.5, 1, 1.5, 2, 3, 4))
  sens$occupied_bins <- vapply(sens$width_mult, function(m) {
    nrow(unique(floor(sweep(sweep(Z, 2, lo, "-"), 2, w0 * m, "/"))))
  }, integer(1))
  blocks <- data.frame(block = names(split(seq_len(n), BLOCK)))
  blocks$n_all <- as.vector(tapply(rep(1, n), BLOCK, sum)[blocks$block])
  blocks$n_kept <- as.vector(tapply(kept_flag, BLOCK, sum)[blocks$block])
  blocks$quota <- pmin(blocks$n_all, block_min)

  out <- df
  Etab <- as.data.frame(E)
  names(Etab) <- paste0("env_", vars)
  Ztab <- as.data.frame(Z)
  names(Ztab) <- paste0("ax_", ax_names)
  out <- cbind(out, Etab, Ztab)
  out$bin <- BIN
  out$block <- BLOCK
  out$kept <- kept_flag
  out$keep_reason <- KEEP
  rownames(out) <- NULL
  res <- list(
    data = out, pca = pca, bin_widths = w, bin_origin = lo,
    report = list(reasons = as.data.frame(table(keep_reason = KEEP,
                                                useNA = "ifany")),
                  window = win_tab, ranges = ranges, intervals = intervals,
                  removed_distance = stats::quantile(d_rm, c(0.5, 0.9, 0.95,
                                                             0.99, 1)),
                  width_sensitivity = sens, blocks = blocks,
                  blocks_below_quota = sum(blocks$n_kept < blocks$quota)),
    settings = list(vars = vars, axes = axes, n_axes = n_ax, log_vars = lv,
                    ax_names = ax_names, radius_km = radius_km,
                    window = window, width_quantile = width_quantile,
                    width_mult = width_mult, max_per_bin = max_per_bin,
                    pick = pick, block_km = block_km, block_min = block_min,
                    block_crs = block_crs, n_input = n))
  class(res) <- "env_thinning"
  res
}

#' @export
print.env_thinning <- function(x, ...) {
  s <- x$settings
  k <- sum(x$data$kept)
  cat("Environmental thinning: kept ", k, " of ", s$n_input, " records (",
      round(100 * k / s$n_input, 1), "%)\n", sep = "")
  cat("Space: ", if (s$axes == "pca") paste0("PC1-PC", s$n_axes) else
    "standardised layers", "; bin widths: ",
    paste(names(x$bin_widths), signif(x$bin_widths, 3), sep = " = ",
          collapse = ", "), "\n", sep = "")
  cat("Geographic quota: ", s$block_min, " per ", s$block_km,
      " km block; blocks below quota: ", x$report$blocks_below_quota, "\n\n",
      sep = "")
  r <- x$report$reasons
  r$keep_reason <- ifelse(is.na(r$keep_reason), "removed",
                          as.character(r$keep_reason))
  print(r, row.names = FALSE)
  cat("\nRanges, all vs kept (must be identical):\n")
  rg <- x$report$ranges
  rg[-1] <- lapply(rg[-1], signif, 4)
  print(rg, row.names = FALSE)
  invisible(x)
}

#' @rdname kept
#' @export
kept.env_thinning <- function(x) {
  d <- x$data[x$data$kept, , drop = FALSE]
  rownames(d) <- NULL
  d
}

#' @rdname kept
#' @export
removed.env_thinning <- function(x) {
  d <- x$data[!x$data$kept, , drop = FALSE]
  rownames(d) <- NULL
  d
}

# --- helpers -------------------------------------------------------------------------

# Bioclim precipitation amounts, recognised by name (internal)
default_log_vars <- function(vars) {
  vars[grepl("(^|[_.-])bio_?(12|13|14|16|17|18|19)([_.-]|$)", vars,
             ignore.case = TRUE)]
}

log_tf <- function(X, lv) {
  if (length(lv)) X[, lv] <- log1p(X[, lv])
  X
}

# Layer values at cell numbers (NA allowed) as a matrix (internal)
env_at_cells <- function(r, cells) {
  X <- matrix(NA_real_, length(cells), terra::nlyr(r),
              dimnames = list(NULL, names(r)))
  ok <- !is.na(cells)
  if (any(ok)) {
    X[ok, ] <- as.matrix(terra::extract(r, cells[ok]))[, names(r),
                                                       drop = FALSE]
  }
  X
}

# Range of the axes and layers within each record's window (internal)
window_ranges <- function(r, cell, lat, radius_km, window, to_axes, geo) {
  n <- length(cell)
  rc <- terra::rowColFromCell(r, cell)
  nl <- terra::nlyr(r)
  use_km <- !is.na(radius_km)
  if (use_km) {
    if (geo) {
      ns <- terra::res(r)[2] * 111.32
      ew <- terra::res(r)[1] * 111.32 * cos(lat * pi / 180)
    } else {
      ns <- terra::res(r)[2] / 1000
      ew <- rep(terra::res(r)[1] / 1000, n)
    }
    Dr <- ceiling(radius_km / ns) + 1
    Dc <- ceiling(radius_km / min(ew)) + 1
  } else {
    Dr <- Dc <- window
  }
  amin <- amax <- NULL
  rmin <- matrix(Inf, n, nl)
  rmax <- matrix(-Inf, n, nl)
  for (dr in -Dr:Dr) for (dc in -Dc:Dc) {
    inc <- if (use_km) {
      gap <- sqrt(pmax(0, abs(dr) * ns - ns / 2)^2 +
                    pmax(0, abs(dc) * ew - ew / 2)^2)
      gap <= radius_km
    } else rep(TRUE, n)
    if (!any(inc)) next
    cl <- rep(NA_real_, n)
    cl[inc] <- terra::cellFromRowCol(r, rc[inc, 1] + dr, rc[inc, 2] + dc)
    X <- env_at_cells(r, cl)
    A <- to_axes(X)
    if (is.null(amin)) {
      amin <- matrix(Inf, n, ncol(A))
      amax <- matrix(-Inf, n, ncol(A))
    }
    amin <- pmin(amin, A, na.rm = TRUE)
    amax <- pmax(amax, A, na.rm = TRUE)
    rmin <- pmin(rmin, X, na.rm = TRUE)
    rmax <- pmax(rmax, X, na.rm = TRUE)
  }
  colnames(rmax) <- colnames(rmin) <- names(r)
  list(axes = amax - amin, raw = rmax - rmin)
}

# Nearest-neighbour distance from each row of Q to the rows of X (internal)
nn_dist <- function(Q, X, chunk = 1000) {
  out <- numeric(nrow(Q))
  if (!nrow(Q)) return(out)
  x2 <- rowSums(X^2)
  for (s in seq(1, nrow(Q), by = chunk)) {
    i <- s:min(s + chunk - 1, nrow(Q))
    q <- Q[i, , drop = FALSE]
    d2 <- outer(rowSums(q^2), x2, "+") - 2 * q %*% t(X)
    out[i] <- sqrt(pmax(apply(d2, 1, min), 0))
  }
  out
}

protect_flag <- function(protect, df) {
  n <- nrow(df)
  if (is.null(protect)) return(rep(FALSE, n))
  if (is.logical(protect)) {
    if (length(protect) != n) {
      stop("A logical `protect` needs one value per record.", call. = FALSE)
    }
    return(!is.na(protect) & protect)
  }
  if (!is.list(protect) || is.null(names(protect))) {
    stop("`protect` must be a logical vector or a named list such as ",
         "list(source = \"sNMF\").", call. = FALSE)
  }
  miss <- setdiff(names(protect), names(df))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  Reduce(`|`, lapply(names(protect), function(k) df[[k]] %in% protect[[k]]))
}

#' Plot an environmental thinning
#'
#' Kept and removed records in environmental space (one or two pairs of
#' axes) and on the map, coloured by why each was kept.
#'
#' @param x An `env_thinning` from [thin_environment()].
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?thin_environment
plot_env_thinning <- function(x, base_size = 10) {
  if (!inherits(x, "env_thinning")) {
    stop("`x` must come from thin_environment().", call. = FALSE)
  }
  d <- x$data
  ax <- x$settings$ax_names
  pairs <- list(c(1, 2))
  if (length(ax) >= 4) pairs[[2]] <- c(3, 4) else
    if (length(ax) == 3) pairs[[2]] <- c(1, 3)
  grp <- ifelse(is.na(d$keep_reason), "removed",
                ifelse(d$keep_reason %in% c("extreme_minmax", "extreme_tail",
                                            "hull"), "extreme / hull",
                       d$keep_reason))
  lev <- c("removed", "env_bin", "block_quota", "protected", "extreme / hull")
  grp <- factor(grp, levels = lev)
  lab <- c(removed = "removed", env_bin = "kept: one per bin",
           block_quota = "kept: block quota", protected = "kept: protected",
           `extreme / hull` = "kept: extreme / hull")
  panels <- lapply(pairs, function(pp) {
    data.frame(panel = paste(ax[pp[1]], "vs", ax[pp[2]]),
               x = d[[paste0("ax_", ax[pp[1]])]],
               y = d[[paste0("ax_", ax[pp[2]])]], group = grp)
  })
  panels[[length(panels) + 1]] <- data.frame(panel = "Geography", x = d$lon,
                                             y = d$lat, group = grp)
  pd <- do.call(rbind, panels)
  pd$panel <- factor(pd$panel, levels = unique(pd$panel))
  pd <- pd[order(pd$group != "removed", decreasing = FALSE), ]
  ggplot2::ggplot(pd, ggplot2::aes(.data$x, .data$y, colour = .data$group,
                                   shape = .data$group)) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::scale_colour_manual(values = c(removed = "#c9c8c3",
                                            env_bin = "#1c5cab",
                                            block_quota = "#e5822e",
                                            protected = "#0f7b5f",
                                            `extreme / hull` = "#1f1f1d"),
                                 labels = lab, drop = FALSE, name = NULL) +
    ggplot2::scale_shape_manual(values = c(removed = 16, env_bin = 16,
                                           block_quota = 16, protected = 17,
                                           `extreme / hull` = 1),
                                labels = lab, drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~panel, scales = "free") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = paste0("Environmental thinning: kept ",
                                 sum(d$kept), " of ", nrow(d), " records")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom")
}

#' Write an environmental thinning to files
#'
#' Writes `<prefix>.csv` (kept records: lon, lat and, if present, species),
#' `<prefix>_all_data.csv` (kept records with all columns),
#' `<prefix>_env_table.csv` (all records with values, axes, bin, block,
#' `kept` and `keep_reason`), `<prefix>_pca_loadings.csv`,
#' `<prefix>_checks.csv` (ranges all vs kept), `<prefix>_pca.rds` (to re-use
#' as `template`) and `<prefix>_check.png`.
#'
#' @param x An `env_thinning`.
#' @param dir Output folder.
#' @param prefix File name prefix.
#' @return The paths of the written files, invisibly.
#' @export
#' @examples
#' # see ?thin_environment
write_env_thinning <- function(x, dir, prefix = "presences_envThinned") {
  if (!inherits(x, "env_thinning")) {
    stop("`x` must come from thin_environment().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, s))
  k <- kept(x)
  cols <- intersect(c("species", "lon", "lat"), names(k))
  utils::write.csv(k[cols], f(".csv"), row.names = FALSE)
  utils::write.csv(k, f("_all_data.csv"), row.names = FALSE)
  utils::write.csv(x$data, f("_env_table.csv"), row.names = FALSE)
  ld <- as.data.frame(x$pca$rotation)
  ld <- cbind(layer = rownames(ld), ld)
  ve <- x$pca$sdev^2 / sum(x$pca$sdev^2)
  ld <- rbind(ld, cbind(layer = "explained_variance",
                        as.data.frame(t(stats::setNames(ve, colnames(
                          x$pca$rotation))))))
  utils::write.csv(ld, f("_pca_loadings.csv"), row.names = FALSE)
  utils::write.csv(x$report$ranges, f("_checks.csv"), row.names = FALSE)
  tmpl <- x
  tmpl$data <- NULL
  saveRDS(tmpl, f("_pca.rds"))
  ggplot2::ggsave(f("_check.png"), plot_env_thinning(x), width = 11,
                  height = 4.5, dpi = 150)
  invisible(c(points = f(".csv"), all_data = f("_all_data.csv"),
              env_table = f("_env_table.csv"),
              loadings = f("_pca_loadings.csv"), checks = f("_checks.csv"),
              pca = f("_pca.rds"), plot = f("_check.png")))
}
