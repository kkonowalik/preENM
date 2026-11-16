# Extrapolation risk of a projection ---------------------------------------------

#' Where does a projection leave the conditions of the calibration data?
#'
#' Maps how far the conditions of a projection (another period, scenario or
#' region) are from those the model was calibrated on, with the two
#' standard measures, so that predictions there can be masked or read with
#' caution:
#'
#' | Layer | Meaning |
#' |---|---|
#' | `nt1` | ExDet novelty of single variables: 0 inside the reference range of every variable, negative outside (the sum over variables of the distance beyond the range, in units of the range) |
#' | `nt2` | ExDet novelty of combinations: Mahalanobis distance to the reference centre divided by the largest such distance in the reference; above 1 = a combination of values that the reference lacks, even if each value is within range |
#' | `mop` | mobility-oriented parity: mean distance, in standardised units, to the nearest `mop_share` of the reference; larger = less similar |
#' | `n_outside` | number of variables outside the reference range |
#' | `most_novel` | the variable farthest outside its reference range (its position in `names(reference)`; `NA` where none is outside) |
#'
#' `nt1` and `n_outside` show strict extrapolation; `nt2` detects new
#' combinations (e.g. warm and wet where the reference has warm and dry),
#' which MESS and range checks miss; `mop` grades similarity also inside
#' the reference range. [niche_truncation()] summarises the same question
#' per variable in a table; this function maps it.
#'
#' @param reference The calibration conditions: the rasters of the
#'   calibration area (M), from which `n_ref` cells are sampled, or a
#'   matrix or data.frame of values (e.g. at the records and background
#'   points used for fitting).
#' @param projection Rasters to project to, with the same layer names.
#' @param mop_share Share of the reference used as nearest neighbours of
#'   `mop` (Owens et al. used 0.1).
#' @param n_ref Number of reference cells sampled from a raster.
#' @param seed Seed for the sample.
#' @param filename Optional output file.
#' @param overwrite Overwrite `filename`.
#' @return A [terra::SpatRaster] with the five layers.
#' @references Mesgaran, M.B., Cousens, R.D., Webber, B.L. (2014) Here be
#'   dragons: a tool for quantifying novelty due to covariate range and
#'   correlation change when projecting species distribution models.
#'   *Diversity and Distributions* 20, 1147-1159. Owens, H.L. et al. (2013)
#'   Constraints on interpretation of ecological niche models by limited
#'   environmental ranges on calibration areas. *Ecological Modelling* 263,
#'   10-18.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' m <- terra::crop(sim$env, terra::ext(14, 21, 49, 55))    # calibration area
#' ex <- extrapolation_risk(m, sim$env)
#' terra::plot(ex[[c("nt1", "nt2", "mop")]])
extrapolation_risk <- function(reference, projection, mop_share = 0.1,
                               n_ref = 5000, seed = 1, filename = "",
                               overwrite = FALSE) {
  proj <- as_grid(projection)
  if (is.matrix(reference) || is.data.frame(reference)) {
    R <- as.matrix(reference[, vapply(as.data.frame(reference), is.numeric,
                                      logical(1)), drop = FALSE])
  } else {
    ref <- as_grid(reference)
    valid <- valid_cells(ref)
    s <- with_seed(seed, sort(valid[sample.int(length(valid),
                                               min(length(valid), n_ref))]))
    R <- read_cells(ref, s)
    colnames(R) <- names(ref)
  }
  R <- R[stats::complete.cases(R), , drop = FALSE]
  vars <- colnames(R)
  miss <- setdiff(vars, names(proj))
  if (is.null(vars) || length(miss)) {
    stop("`projection` lacks layer(s) of the reference: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (nrow(R) < 10) stop("Fewer than 10 reference rows with values.",
                         call. = FALSE)
  if (mop_share <= 0 || mop_share > 1) {
    stop("`mop_share` must be between 0 and 1.", call. = FALSE)
  }
  lo <- apply(R, 2, min)
  hi <- apply(R, 2, max)
  span <- ifelse(hi > lo, hi - lo, 1)
  mu <- colMeans(R)
  sdv <- apply(R, 2, stats::sd)
  sdv[sdv == 0] <- 1
  Rs <- sweep(sweep(R, 2, mu, "-"), 2, sdv, "/")
  k <- max(1, ceiling(mop_share * nrow(R)))
  S <- tryCatch(solve(stats::cov(R)), error = function(e) NULL)
  if (is.null(S)) {
    warning("The reference variables are perfectly collinear: `nt2` is NA.",
            call. = FALSE)
  } else {
    d_max <- sqrt(max(stats::mahalanobis(R, mu, S, inverted = TRUE)))
  }
  measure <- function(m) {
    out <- matrix(NA_real_, nrow(m), 5)
    ok <- which(stats::complete.cases(m))
    if (!length(ok)) return(out)
    P <- m[ok, , drop = FALSE]
    below <- sweep(sweep(P, 2, lo, "-"), 2, span, "/")       # < 0 = below
    above <- sweep(-sweep(P, 2, hi, "-"), 2, span, "/")      # < 0 = above
    beyond <- pmin(below, above, 0)
    out[ok, 1] <- rowSums(beyond)
    if (!is.null(S)) {
      out[ok, 2] <- sqrt(stats::mahalanobis(P, mu, S, inverted = TRUE)) / d_max
    }
    Ps <- sweep(sweep(P, 2, mu, "-"), 2, sdv, "/")
    # in blocks, so that the neighbour table stays small
    size <- max(1L, floor(2e6 / k))
    mop <- numeric(length(ok))
    for (s in seq(1L, length(ok), by = size)) {
      i <- s:min(s + size - 1L, length(ok))
      mop[i] <- rowMeans(dbscan::kNN(Rs, k = k,
                                     query = Ps[i, , drop = FALSE])$dist)
    }
    out[ok, 3] <- mop
    out[ok, 4] <- rowSums(beyond < 0)
    worst <- max.col(-beyond, ties.method = "first")
    out[ok, 5] <- ifelse(out[ok, 4] > 0, worst, NA)
    out
  }
  res <- terra::app(proj[[vars]], measure)
  names(res) <- c("nt1", "nt2", "mop", "n_outside", "most_novel")
  if (nzchar(filename)) {
    res <- terra::writeRaster(res, filename, overwrite = overwrite)
  }
  res
}
