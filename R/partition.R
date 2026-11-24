# Partitioning the uncertainty of projections --------------------------------------

#' Which source of uncertainty dominates a set of projections, and where?
#'
#' A realistic set of projections varies at the same time in climate model
#' (GCM), emission scenario, modelling algorithm and cross-validation run
#' (or any other factors). `partition_uncertainty()` splits, in every cell,
#' the variance of the projected suitability among these factors, so you
#' can see which choice matters most, and where.
#'
#' **Method.** For each cell, the projections are the response of a linear
#' model with the factors as main effects, as in an analysis of variance;
#' the variance explained by each factor (its sum of squares divided by
#' the total sum of squares) is its share. What the main effects do not
#' explain (interactions between factors, and replicate noise) is reported
#' as `interactions/residual`. The shares sum to 1 in every cell. This is
#' the variance partitioning used for climate projections (Hawkins and
#' Sutton 2009) and their ecological consequences (Wilson et al. 2024).
#'
#' * With a **balanced design** (every combination of levels equally often)
#'   the factors are independent and the shares do not depend on their
#'   order. Otherwise sequential (type I) sums of squares are used, in the
#'   order of the columns of `design`, and a warning says that the shares
#'   then depend on that order.
#' * The share of **scenario** depends on which scenarios are included: it
#'   describes the spread between the chosen futures, not an error, so it
#'   cannot be compared directly with the share of climate models (Wilson
#'   et al. 2024).
#' * Shares are relative: where the projections hardly differ, a large
#'   share of a tiny variance means little. The `variance` layer and the
#'   pooled shares (each cell weighted by its variance) show where the
#'   variance is.
#'
#' The rasters are read and written block by block, so stacks larger than
#' memory work.
#'
#' @param projections Projection rasters: a [terra::SpatRaster] with one
#'   layer per projection, or raster files and/or folders (as in
#'   [ensemble_models()]). All on the same grid.
#' @param design A data.frame with one row per projection (in the order of
#'   the layers, or with a column `layer` naming them) and one column per
#'   factor, e.g. `gcm`, `scenario`, `algorithm`, `run`.
#' @param region Optional part of the area for a second set of summaries:
#'   polygons (`sf`, `SpatVector`) or a raster (cells that are not `NA` and
#'   not 0), as in [transfer_bias_test()].
#' @param filename Optional GeoTIFF for the result.
#' @param overwrite Overwrite `filename`.
#' @param rescale Bring every projection to the 0-1 scale (layers on 0-100
#'   or 0-1000 are divided by 100 or 1000, with a message); see
#'   [extract_predictions()].
#' @param pattern,recursive File search in folders.
#' @param seed Seed for the cell sample used to recognise the scales.
#' @return A list of class `uncertainty_partition`: `raster` (a
#'   [terra::SpatRaster] with one share layer per factor, `share_residual`
#'   for interactions and residual, `variance` (the total variance between
#'   projections) and `dominant` (the factor with the largest share,
#'   categorical)), `shares` (the shares pooled over the area and, with
#'   `region`, over the region, each cell weighted by its variance, and the
#'   share of cells each factor dominates), `factors`, `balanced`, `design`
#'   and `scales`.
#' @references Wilson, R.J., Kay, S., Ciavatta, S. (2024) Partitioning
#'   climate uncertainty in ecological projections: Pacific oysters in a
#'   hotter Europe. *Ecological Informatics* 80, 102537.
#'   \doi{10.1016/j.ecoinf.2024.102537}
#'
#'   Hawkins, E., Sutton, R. (2009) The potential to narrow uncertainty in
#'   regional climate predictions. *Bulletin of the American Meteorological
#'   Society* 90, 1095-1108. \doi{10.1175/2009BAMS2607.1}
#' @seealso [plot_partition()], [ensemble_models()],
#'   [transfer_bias_test()], [extrapolation_risk()]
#' @export
#' @examples
#' # 2 GCMs x 2 scenarios x 3 algorithms
#' sim <- simulate_enm_data(seed = 1)
#' base <- sim$models[[c("GLM", "GBM", "Maxent_LQ")]]
#' design <- expand.grid(algorithm = names(base), gcm = c("A", "B"),
#'                       scenario = c("ssp1", "ssp5"),
#'                       stringsAsFactors = FALSE)
#' proj <- terra::rast(lapply(seq_len(nrow(design)), function(i) {
#'   d <- design[i, ]
#'   shift <- 0.05 * (d$gcm == "B") + 0.15 * (d$scenario == "ssp5")
#'   terra::clamp(base[[d$algorithm]] - shift, 0, 1)
#' }))
#' pu <- partition_uncertainty(proj, design)
#' pu
#' plot_partition(pu)
partition_uncertainty <- function(projections, design, region = NULL,
                                  filename = NULL, overwrite = FALSE,
                                  rescale = TRUE,
                                  pattern = "\\.tif{1,2}$",
                                  recursive = FALSE, seed = 1) {
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  x <- ensemble_stack(projections, pattern, recursive,
                      unique_names = "layer" %in% names(design))
  if ("layer" %in% names(design)) {
    miss <- setdiff(design$layer, names(x))
    if (length(miss)) {
      stop("No layer for: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    x <- x[[design$layer]]
    design$layer <- NULL
  } else if (nrow(design) != terra::nlyr(x)) {
    stop("`design` must have one row per layer (", terra::nlyr(x), ").",
         call. = FALSE)
  }
  facs <- names(design)
  if (!length(facs)) stop("`design` needs at least one factor column.",
                          call. = FALSE)
  design[] <- lapply(design, function(v) factor(v, unique(v)))
  one <- facs[vapply(design, nlevels, integer(1)) < 2]
  if (length(one)) {
    stop("Factor(s) with a single level: ", paste(one, collapse = ", "),
         ". Remove them from `design`.", call. = FALSE)
  }
  bad_names <- intersect(facs, c("residual", "variance", "dominant"))
  if (length(bad_names)) {
    stop("Rename the factor(s) ", paste(bad_names, collapse = ", "), ".",
         call. = FALSE)
  }
  n <- nrow(design)
  balanced <- length(unique(as.vector(table(design)))) == 1 &&
    all(table(design) > 0)
  if (!balanced) {
    warning("The design is not balanced (not every combination of levels ",
            "appears equally often): the shares depend on the order of the ",
            "factors (", paste(facs, collapse = ", "), "; sequential sums of ",
            "squares).", call. = FALSE)
  }

  # design matrix, orthogonalised once for all cells
  X <- stats::model.matrix(stats::reformulate(facs), design)
  qx <- qr(X)
  r <- qx$rank
  if (r >= n) {
    warning("As many parameters as projections: nothing is left for ",
            "interactions or residual, which are therefore 0.",
            call. = FALSE)
  }
  Q <- qr.Q(qx)[, seq_len(r), drop = FALSE]
  term <- attr(X, "assign")[qx$pivot[seq_len(r)]]

  # scales
  smp <- with_seed(seed, sample_complete(x, 5000))
  h <- harmonise_scales(as.data.frame(smp), rescale, "projections")
  div <- unname(h$divisor[names(x)])

  shares_of <- function(v) {
    v <- sweep(v, 2, div, "/")
    ok <- stats::complete.cases(v)
    out <- matrix(NA_real_, nrow(v), length(facs) + 3)
    if (!any(ok)) return(out)
    y <- v[ok, , drop = FALSE]
    z <- y %*% Q
    tot <- rowSums((y - rowMeans(y))^2)
    ss <- vapply(seq_along(facs), function(k) {
      rowSums(z[, term == k, drop = FALSE]^2)
    }, numeric(nrow(y)))
    ss <- matrix(ss, nrow = nrow(y))
    flat <- tot <= 1e-12 * pmax(1, rowSums(y^2))
    sh <- ss / tot
    res <- pmax(0, 1 - rowSums(sh))
    sh <- cbind(sh, res)
    sh[flat, ] <- NA
    dom <- max.col(sh, ties.method = "first")
    dom[flat] <- NA
    out[ok, ] <- cbind(sh, tot / (n - 1), dom)
    out
  }

  lay <- c(paste0("share_", facs), "share_residual", "variance", "dominant")
  out <- terra::rast(x, nlyrs = length(lay))
  names(out) <- lay
  f <- if (is.null(filename)) tempfile("preENM_partition", fileext = ".tif")
    else filename
  b <- terra::writeStart(out, f, overwrite = overwrite || is.null(filename),
                         names = lay)
  terra::readStart(x)
  on.exit(try(terra::readStop(x), silent = TRUE), add = TRUE)
  for (i in seq_len(b$n)) {
    v <- terra::readValues(x, b$row[i], b$nrows[i], 1, terra::ncol(x),
                           mat = TRUE)
    terra::writeValues(out, shares_of(v), b$row[i], b$nrows[i])
  }
  out <- terra::writeStop(out)
  terra::readStop(x)
  names(out) <- lay
  src <- c(facs, "interactions/residual")
  dom <- out[["dominant"]]
  levels(dom) <- data.frame(id = seq_along(src), source = src)
  out <- c(out[[setdiff(lay, "dominant")]], dom)
  names(out) <- lay

  # pooled shares: each cell weighted by its variance
  pooled <- function(mask = NULL) {
    sh <- out[[seq_along(src)]]
    vv <- out[["variance"]]
    if (!is.null(mask)) {
      sh <- terra::mask(sh, mask)
      vv <- terra::mask(vv, mask)
    }
    num <- terra::global(sh * vv, "sum", na.rm = TRUE)[, 1]
    den <- terra::global(vv, "sum", na.rm = TRUE)[, 1]
    dom <- terra::freq(if (is.null(mask)) out[["dominant"]] else
      terra::mask(out[["dominant"]], mask))
    cnt <- stats::setNames(rep(0, length(src)), src)
    if (!is.null(dom) && nrow(dom)) {
      lv <- as.character(dom$value)
      cnt[lv] <- dom$count
    }
    data.frame(source = src, share = if (den > 0) num / den else NA_real_,
               dominant_cells = as.vector(cnt),
               dominant_share = as.vector(cnt) / max(1, sum(cnt)),
               median_sd = NA_real_,
               stringsAsFactors = FALSE)
  }
  tab <- cbind(scope = "area", pooled())
  vsmp <- with_seed(seed, sample_complete(out[["variance"]], 1e5))[, 1]
  tab$median_sd <- sqrt(stats::median(vsmp, na.rm = TRUE))
  if (!is.null(region)) {
    m <- region_mask(region, out)
    rt <- cbind(scope = "region", pooled(m))
    rv <- with_seed(seed, sample_complete(terra::mask(out[["variance"]], m),
                                          1e5))[, 1]
    rt$median_sd <- sqrt(stats::median(rv, na.rm = TRUE))
    tab <- rbind(tab, rt)
  }
  rownames(tab) <- NULL
  structure(list(raster = out, shares = tab, factors = facs,
                 balanced = balanced, design = design, scales = h$table,
                 n = n),
            class = "uncertainty_partition")
}

# A region (polygons or raster) as a mask on the grid of r (internal):
# cells inside are 1, others NA; built by terra on disk when large
region_mask <- function(region, r) {
  g <- terra::rast(r[[1]])
  if (is.character(region)) {
    region <- if (grepl("\\.(tif|tiff|asc|nc|vrt|img)$", region[1],
                        ignore.case = TRUE)) terra::rast(region) else
      terra::vect(region)
  }
  if (inherits(region, c("sf", "sfc"))) region <- terra::vect(region)
  if (inherits(region, "SpatVector")) {
    if (!terra::same.crs(region, g)) region <- terra::project(region, g)
    m <- terra::rasterize(region, g, field = 1)
  } else if (inherits(region, "SpatRaster")) {
    m <- region[[1]]
    if (!isTRUE(terra::compareGeom(m, g, stopOnError = FALSE))) {
      m <- terra::project(m, g, method = "near")
    }
    m <- terra::ifel(is.na(m) | m == 0, NA, 1)
  } else {
    stop("`region` must be polygons or a raster.", call. = FALSE)
  }
  if (is.na(terra::global(m, "max", na.rm = TRUE)[1, 1])) {
    stop("`region` covers no cell.", call. = FALSE)
  }
  m
}

#' @export
print.uncertainty_partition <- function(x, digits = 3, ...) {
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  pc <- function(v) share_label(v)
  cat("Uncertainty partition: ", x$n, " projections; factors ",
      paste(x$factors, collapse = ", "), if (!x$balanced)
        " (unbalanced design)", "\n\nIN SHORT\n", sep = "")
  for (sc in unique(x$shares$scope)) {
    t <- x$shares[x$shares$scope == sc, , drop = FALSE]
    where <- if (sc == "area") "Across the area" else "In the region"
    if (all(is.na(t$share))) {
      say(paste0(where, ", the projections do not differ: there is no ",
                 "variance to partition."))
      next
    }
    o <- order(-t$share)
    lead <- t$source[o[1]]
    say(paste0(where, ", ", if (lead == "interactions/residual")
      "interactions between the factors and residual variation explain"
      else paste0("the choice of ", lead, " explains"),
      " most of the variance in predicted suitability (",
      pc(t$share[o[1]]), "), followed by ", paste0(
        t$source[o[-1]], " (", pc(t$share[o[-1]]), ")", collapse = ", "),
      ". ", {
        od <- which.max(t$dominant_share)
        if (t$source[od] == lead) {
          paste0(lead, " is the largest source in ",
                 pc(t$dominant_share[od]), " of the cells.")
        } else {
          paste0("Cell by cell, ", t$source[od], " is the largest source ",
                 "most often (", pc(t$dominant_share[od]), " of the cells), ",
                 "but ", lead, " leads where the variance is larger.")
        }
      }))
    if (is.finite(t$median_sd[1]) && t$median_sd[1] < 0.02) {
      say(paste0(where, ", the projections hardly differ (median standard ",
                 "deviation ", formatC(t$median_sd[1], digits = 3,
                                        format = "f"), "): the shares ",
                 "describe small differences."))
    }
  }
  if (!x$balanced) {
    say("The design is unbalanced: the shares depend on the order of the factors in `design`.")
  }
  if (any(grepl("scen|ssp|rcp", x$factors, ignore.case = TRUE))) {
    say("The share of scenarios depends on which scenarios were chosen: it is the spread between the futures considered, not an error, and is not directly comparable with the share of climate models (Wilson et al. 2024).")
  }
  cat("\nShares (pooled, each cell weighted by its variance):\n")
  t <- x$shares
  t$share <- pc(t$share)
  t$dominant_share <- pc(t$dominant_share)
  t$median_sd <- round(t$median_sd, digits)
  print(t, row.names = FALSE)
  cat("", strwrap(paste0(
    "share = sum of squares of a factor / total sum of squares, from a ",
    "linear model with the factors as main effects (as in an analysis of ",
    "variance); interactions/residual = the rest. dominant = the source ",
    "with the largest share in a cell. The 0.02 limit of the standard ",
    "deviation is a rule of thumb."), width = 78), sep = "\n")
  invisible(x)
}

#' Plot the partition of projection uncertainty
#'
#' The map of the dominant source of variance in every cell, next to a
#' stacked bar of the shares pooled over the area (and the region). Legend
#' entries give the number and share of cells each source dominates.
#'
#' @param x A result of [partition_uncertainty()].
#' @param type `"both"`, `"map"` or `"bar"`.
#' @param max_cells Largest number of cells drawn.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches (default 12 x 5.5
#'   for both panels).
#' @return A ggplot, or for `type = "both"` invisibly a list of two.
#' @references Wilson, R.J., Kay, S., Ciavatta, S. (2024) Partitioning
#'   climate uncertainty in ecological projections: Pacific oysters in a
#'   hotter Europe. *Ecological Informatics* 80, 102537.
#'   \doi{10.1016/j.ecoinf.2024.102537}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' base <- sim$models[[c("GLM", "GBM", "Maxent_LQ")]]
#' design <- expand.grid(algorithm = names(base), gcm = c("A", "B"),
#'                       stringsAsFactors = FALSE)
#' proj <- terra::rast(lapply(seq_len(nrow(design)), function(i) {
#'   terra::clamp(base[[design$algorithm[i]]] -
#'                  0.1 * (design$gcm[i] == "B"), 0, 1)
#' }))
#' plot_partition(partition_uncertainty(proj, design))
plot_partition <- function(x, type = c("both", "map", "bar"),
                           max_cells = 5e5, base_size = 11, file = NULL,
                           width = NULL, height = NULL) {
  if (!inherits(x, "uncertainty_partition")) {
    stop("`x` must come from partition_uncertainty().", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.null(width)) width <- if (type == "both") 12 else 6.5
  if (is.null(height)) height <- 5.5
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  src <- c(x$factors, "interactions/residual")
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02", "#a6761d")[seq_along(x$factors)]
  pal <- c(pal, "#bdbdbd")
  names(pal) <- src
  t <- x$shares[x$shares$scope == "area", , drop = FALSE]
  lab <- paste0(src, ": ", t$dominant_cells[match(src, t$source)], " cells (",
                share_label(t$dominant_share[match(src, t$source)]), ")")
  names(lab) <- src

  r <- x$raster[["dominant"]]
  if (terra::ncell(r) > max_cells) {
    r <- terra::spatSample(r, max_cells, method = "regular",
                           as.raster = TRUE)
  }
  d <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(d)[3] <- "source"
  d$source <- factor(as.character(d$source), src)
  p_map <- ggplot2::ggplot(d, ggplot2::aes(.data$x, .data$y,
                                           fill = .data$source)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_manual(values = pal, labels = lab, drop = FALSE,
                               name = "Dominant source") +
    (if (terra::is.lonlat(r)) ggplot2::coord_quickmap() else
      ggplot2::coord_equal()) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "Dominant source of variance",
                  subtitle = paste0("Largest share of the variance between ",
                                    x$n, " projections, per cell.")) +
    theme_enm(base_size, ink, muted)

  b <- x$shares
  b$source <- factor(b$source, rev(src))
  b$scope <- factor(ifelse(b$scope == "area", "Whole area", "Region"),
                    c("Whole area", "Region"))
  p_bar <- ggplot2::ggplot(b, ggplot2::aes(.data$scope, .data$share,
                                           fill = .data$source)) +
    ggplot2::geom_col(width = 0.6, colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(.data$share >= 0.04,
                                                   share_label(.data$share),
                                                   "")),
                       position = ggplot2::position_stack(vjust = 0.5),
                       colour = "white", size = base_size / 3.2) +
    ggplot2::scale_fill_manual(values = pal, name = "Source",
                               breaks = src) +
    ggplot2::scale_y_continuous(labels = function(v) paste0(100 * v, "%"),
                                expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = "Share of the variance",
                  title = "Shares of the variance",
                  subtitle = "Pooled over cells, each weighted by its variance.") +
    theme_enm(base_size, ink, muted)
  p <- switch(type, map = p_map, bar = p_bar, list(map = p_map, bar = p_bar))
  if (!is.null(file)) save_figure(p, file, width, height)
  if (type == "both") {
    print_side_by_side(p)
    return(invisible(p))
  }
  p
}
