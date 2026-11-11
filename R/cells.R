# One point per raster cell and background points ------------------------------

#' Keep one record per raster cell, centred in the cell
#'
#' Reduces presence records to one per cell of the rasters you will model with
#' (per species), and moves each kept record to the centre of its cell. This
#' removes the sampling bias of many records in the same cell and makes every
#' record match exactly one cell of the predictors. Records outside the
#' rasters, or in cells where *any* layer has no value, are dropped.
#'
#' When several records fall in one cell, the first one is kept. With
#' `priority`, records from preferred sources win (e.g. your own verified
#' samples over GBIF records); otherwise the input order decides.
#'
#' This is geographic thinning to the raster resolution. It is not the same
#' as environmental thinning (thinning in climate space).
#'
#' @param x Occurrence records: a data.frame with coordinates (e.g. from
#'   [merge_occurrences()] or [kept()]), a CSV path, an `sf` object or a
#'   [terra::SpatVector].
#' @param raster The predictor layers: a [terra::SpatRaster] (several layers
#'   allowed) or file paths. They define the grid; all layers must share
#'   extent, resolution and CRS.
#' @param by Column with species names: one record per cell *and species*.
#'   Ignored if absent.
#' @param priority Optional character vector of values of `priority_col`,
#'   from most to least preferred, e.g. `c("sNMF", "GBIF")`.
#' @param priority_col Column holding the source (default `"source"`, as
#'   created by [merge_occurrences()]).
#' @param centre Move kept records to the centre of their cell (default
#'   `TRUE`); the original coordinates are kept in `lon_orig`, `lat_orig`.
#' @param crs CRS of the output coordinates (default WGS84 longitude/latitude).
#' @param coords Optional names of the x/longitude and y/latitude columns of
#'   `x`; detected automatically if `NULL`.
#' @param points_crs CRS of the coordinates in `x` (default WGS84
#'   longitude/latitude).
#' @return A data.frame of the kept records with `lon`, `lat` (cell centres),
#'   `cell` (cell number), `lon_orig`, `lat_orig` and the other columns of
#'   `x`.
#' @export
#' @examples
#' r <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 40,
#'                  ymax = 50, vals = 1)
#' pts <- data.frame(lon = c(0.2, 0.7, 5.5, 20), lat = c(40.2, 40.6, 45.5, 45),
#'                   source = c("GBIF", "own", "GBIF", "GBIF"))
#' thin_by_cell(pts, r, priority = c("own", "GBIF"))
thin_by_cell <- function(x, raster, by = "species", priority = NULL,
                         priority_col = "source", centre = TRUE,
                         crs = "EPSG:4326", coords = NULL,
                         points_crs = "EPSG:4326") {
  r <- as_grid(raster)
  if (is.character(x) && length(x) == 1) x <- read_any_csv(x)
  if (is.data.frame(x)) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    xy_cols <- if (is.null(coords)) guess_coords(names(df)) else coords
    lon_in <- as_number(df[[xy_cols[1]]])
    lat_in <- as_number(df[[xy_cols[2]]])
    has_xy <- !is.na(lon_in) & !is.na(lat_in)
    p <- terra::vect(cbind(lon_in[has_xy], lat_in[has_xy]), crs = points_crs)
  } else {
    p <- as_points(x, coords, points_crs)
    df <- as.data.frame(p)
    has_xy <- rep(TRUE, nrow(df))
  }
  # every vector below has one element per row of df
  n <- nrow(df)
  xy0 <- matrix(NA_real_, n, 2)
  cell <- rep(NA_real_, n)
  if (any(has_xy)) {
    xy0[has_xy, ] <- terra::crds(terra::project(p, "EPSG:4326"))
    pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
    cell[has_xy] <- terra::cellFromXY(r, terra::crds(pr))
  }
  inside <- !is.na(cell)
  ok <- rep(FALSE, n)
  if (any(inside)) {
    v <- terra::extract(r, cell[inside])
    ok[inside] <- stats::complete.cases(v)
  }
  message(if (any(!has_xy)) paste0("Without coordinates: ", sum(!has_xy),
                                   "; "),
          "Outside the raster: ", sum(has_xy & !inside),
          "; inside but without values in all layers: ", sum(inside & !ok),
          "; remaining: ", sum(ok))
  d <- df[ok, , drop = FALSE]
  d$cell <- cell[ok]
  d$lon_orig <- xy0[ok, 1]
  d$lat_orig <- xy0[ok, 2]
  ord <- seq_len(nrow(d))
  if (!is.null(priority)) {
    if (!priority_col %in% names(d)) {
      stop("`priority` given but column `", priority_col, "` not found.",
           call. = FALSE)
    }
    pr_rank <- match(d[[priority_col]], priority)
    pr_rank[is.na(pr_rank)] <- length(priority) + 1
    ord <- order(pr_rank, ord)
  }
  d <- d[ord, , drop = FALSE]
  grp <- if (!is.null(by) && by %in% names(d)) d[[by]] else ""
  keep <- !duplicated(paste(grp, d$cell))
  message("Records in a cell already taken (dropped): ", sum(!keep))
  d <- d[keep, , drop = FALSE]
  if (centre) {
    ctr <- terra::xyFromCell(r, d$cell)
    ctr <- cell_coords(ctr, r, crs)
    d$lon <- ctr[, 1]
    d$lat <- ctr[, 2]
    if (nrow(d)) {
      shift <- sqrt(((d$lon - d$lon_orig) * cos(d$lat * pi / 180))^2 +
                      (d$lat - d$lat_orig)^2) * 111.2
      if (terra::is.lonlat(terra::vect(cbind(0, 0), crs = crs))) {
        message("Largest shift to a cell centre: ~", round(max(shift), 2),
                " km")
      }
    }
  } else {
    d$lon <- d$lon_orig
    d$lat <- d$lat_orig
  }
  first <- intersect(c(by, "lon", "lat"), names(d))
  d <- d[c(first, setdiff(names(d), first))]
  rownames(d) <- NULL
  message("Thinned records: ", nrow(d))
  d
}

#' Random background (pseudo-absence) points
#'
#' Draws background points from random cells of the predictor rasters, one
#' point per cell (without replacement), only from cells where *all* layers
#' have values. Several independent sets can be drawn, e.g. to repeat models
#' with different backgrounds (see the `set` argument of [evaluate_cv()]).
#' Each set uses its own seed (`seed + set number`), so a set can be
#' reproduced on its own.
#'
#' Restrict the rasters to the accessible area (M) beforehand: the background
#' should represent the area the species could have reached.
#'
#' @param raster The predictor layers ([terra::SpatRaster] or file paths),
#'   clipped to the accessible area.
#' @param n Number of points per set (default 10000). If fewer valid cells
#'   exist, all of them are used, with a warning.
#' @param n_sets Number of independent sets.
#' @param exclude Optional presence records whose cells are excluded: a
#'   data.frame with a `cell` column (as returned by [thin_by_cell()]) or
#'   with coordinates, or a vector of cell numbers.
#' @param seed Optional random seed.
#' @param crs CRS of the output coordinates (default WGS84 longitude/latitude).
#' @return A data.frame with `set`, `cell`, `lon`, `lat` (cell centres).
#' @export
#' @examples
#' r <- terra::rast(nrows = 20, ncols = 20, xmin = 0, xmax = 10, ymin = 40,
#'                  ymax = 50, vals = 1)
#' bg <- sample_background(r, n = 50, n_sets = 3, seed = 1)
#' table(bg$set)
sample_background <- function(raster, n = 10000, n_sets = 1, exclude = NULL,
                              seed = NULL, crs = "EPSG:4326") {
  r <- as_grid(raster)
  ok <- if (terra::nlyr(r) > 1) sum(r) else r
  valid <- terra::cells(ok)
  if (is.list(valid)) valid <- valid[[1]]
  if (!is.null(exclude)) {
    ex <- if (is.numeric(exclude) && is.null(dim(exclude))) {
      exclude
    } else if (is.data.frame(exclude) && "cell" %in% names(exclude)) {
      exclude$cell
    } else {
      p <- as_points(exclude)
      if (!terra::same.crs(p, r)) p <- terra::project(p, r)
      terra::cellFromXY(r, terra::crds(p))
    }
    valid <- setdiff(valid, ex)
  }
  message("Cells available for background: ", length(valid))
  n_take <- min(n, length(valid))
  if (n_take < n) {
    warning("Fewer valid cells (", length(valid), ") than n; using all.",
            call. = FALSE)
  }
  sets <- lapply(seq_len(n_sets), function(k) {
    cl <- with_seed(if (is.null(seed)) NULL else seed + k,
                    valid[sample.int(length(valid), n_take)])
    xy <- cell_coords(terra::xyFromCell(r, cl), r, crs)
    data.frame(set = k, cell = cl, lon = xy[, 1], lat = xy[, 2])
  })
  out <- do.call(rbind, sets)
  rownames(out) <- NULL
  out
}

# Read predictor layers as one SpatRaster (internal)
as_grid <- function(raster) {
  if (is.character(raster)) raster <- terra::rast(raster)
  if (inherits(raster, c("RasterLayer", "RasterStack", "RasterBrick"))) {
    raster <- terra::rast(raster)
  }
  if (!inherits(raster, "SpatRaster")) {
    stop("`raster` must be a SpatRaster or raster file path(s).",
         call. = FALSE)
  }
  raster
}

# Cell-centre coordinates in the output CRS (internal)
cell_coords <- function(xy, r, crs) {
  if (!nrow(xy)) return(xy)
  if (terra::same.crs(r, crs)) return(xy)
  terra::crds(terra::project(terra::vect(xy, crs = terra::crs(r)), crs))
}
