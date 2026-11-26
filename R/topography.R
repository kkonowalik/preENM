# Terrain-derived predictors from an elevation raster

#' Terrain indices from an elevation raster
#'
#' Derives terrain predictors from a digital elevation model (DEM): slope,
#' northness and eastness, terrain ruggedness, topographic position at one
#' or more radii, its scale-free form (deviation from mean elevation), the
#' topographic wetness index and the heat load index. They complement the
#' climatic indices of [climate_indices()]: terrain controls the local
#' microclimate, moisture and soil that coarse climate grids miss, and
#' matters most for plants and for fine grains (Schwager & Berg 2021). They
#' are listed in `climate_index_names(table = TRUE)` under the group `"terrain"`, so
#' `climate_index_names("terrain")` gives their names and the layers can
#' join any predictor stack for [sample_environment()] and
#' [select_variables()].
#'
#' Beyond what `terra::terrain()` offers, the function gives:
#' * **TPI at user-chosen radii** (one layer per radius, `tpi_<r>`), in
#'   cells or kilometres; kilometres are converted to cells separately in x
#'   and y, so the window stays round on longitude/latitude grids.
#' * **DEV** (deviation from mean elevation; De Reu et al. 2013): TPI
#'   divided by the standard deviation of elevation in the same window, so
#'   a ridge in gentle hills and one in high mountains get comparable
#'   values.
#' * **TWI** with its own D8 flow accumulation (terra has none). The
#'   upslope area is summed as real cell areas, so it is correct on
#'   longitude/latitude grids where cells shrink towards the poles.
#' * **HLI** (McCune & Keon 2002, equation 3) using each cell's own
#'   latitude, with aspect folded about the north-east/south-west axis in
#'   the northern hemisphere and about the north-west/south-east axis in the
#'   southern one, so south-west slopes are warmest in the north and
#'   north-west slopes in the south.
#' * Northness and eastness weighted by the sine of slope, so they are 0 on
#'   flat ground, where aspect has no meaning.
#'
#' **Definitions.**
#' * `slope`: steepest gradient (Horn 1981), degrees.
#' * `northness`, `eastness`: \eqn{\cos(aspect) \sin(slope)} and
#'   \eqn{\sin(aspect) \sin(slope)}, -1 to 1.
#' * `tri`: terrain ruggedness index of Riley et al. (1999), the root of
#'   the summed squared elevation differences to the 8 neighbours, in
#'   elevation units.
#' * `tpi_<r>`: elevation minus the mean elevation of the cells within
#'   radius `r` (centre excluded; Weiss 2001). Positive on ridges and convex
#'   slopes, negative in valleys and hollows, 0 on a plane.
#' * `dev_<r>`: `tpi_<r>` divided by the standard deviation of elevation in
#'   the same window; NA where the window is flat.
#' * `twi`: \eqn{\ln(a / \tan\beta)} (Beven & Kirkby 1979), with `a` the
#'   upslope area draining through the cell (the cell included) per unit
#'   contour width (m), and \eqn{\beta} the slope; slopes below `min_slope`
#'   (degrees) are raised to it so flats stay finite.
#' * `hli`: potential annual direct incident radiation, MJ cm\eqn{^{-2}}
#'   yr\eqn{^{-1}}, from slope, folded aspect and latitude. The equation
#'   was fitted for slopes up to 60 degrees and latitudes 30 to 60 degrees;
#'   outside them it extrapolates.
#'
#' **Memory.** Slope, aspect, TRI, TPI, DEV and HLI are computed block-wise
#' by terra and can be written straight to disk with `path`. TWI needs the
#' flow network of the whole DEM in memory, so it stops above `max_cells`
#' cells: clip the DEM to the study area or aggregate it first.
#'
#' @param dem Elevation raster (`SpatRaster` or file path), one layer, in
#'   metres. A projected CRS in metres or longitude/latitude.
#' @param indices Which indices: any of `"slope"`, `"northness"`,
#'   `"eastness"`, `"tri"`, `"tpi"`, `"dev"`, `"twi"`, `"hli"`, or `"all"`.
#' @param radius Radii of the TPI and DEV windows (one layer each).
#' @param radius_unit `"cells"` (default) or `"km"`.
#' @param min_slope Smallest slope used in the TWI, degrees.
#' @param max_cells Largest DEM for the TWI (whole raster in memory).
#' @param path Optional folder: each layer is written to
#'   `<prefix><name>.tif`.
#' @param prefix File name prefix.
#' @param overwrite Overwrite existing files?
#' @return A `SpatRaster` with one layer per index (and per radius for TPI
#'   and DEV), named as in `climate_index_names(table = TRUE)` (`tpi_<r>`, `dev_<r>`,
#'   with `r` followed by `km` when `radius_unit = "km"`).
#' @seealso [climate_indices()] for climatic indices; `climate_index_names(table = TRUE)`
#'   for meanings and references.
#' @references Beven, K.J., Kirkby, M.J. (1979) A physically based, variable
#'   contributing area model of basin hydrology. *Hydrological Sciences
#'   Bulletin* 24, 43-69. \doi{10.1080/02626667909491834}
#'
#'   De Reu, J., Bourgeois, J., Bats, M., et al. (2013) Application of the
#'   topographic position index to heterogeneous landscapes.
#'   *Geomorphology* 186, 39-49. \doi{10.1016/j.geomorph.2012.12.015}
#'
#'   Guisan, A., Weiss, S.B., Weiss, A.D. (1999) GLM versus CCA spatial
#'   modeling of plant species distribution. *Plant Ecology* 143, 107-122.
#'   \doi{10.1023/A:1009841519580}
#'
#'   Horn, B.K.P. (1981) Hill shading and the reflectance map. *Proceedings
#'   of the IEEE* 69, 14-47. \doi{10.1109/PROC.1981.11918}
#'
#'   McCune, B., Keon, D. (2002) Equations for potential annual direct
#'   incident radiation and heat load. *Journal of Vegetation Science* 13,
#'   603-606. \doi{10.1111/j.1654-1103.2002.tb02087.x}
#'
#'   Riley, S.J., DeGloria, S.D., Elliot, R. (1999) A terrain ruggedness
#'   index that quantifies topographic heterogeneity. *Intermountain Journal
#'   of Sciences* 5, 23-27.
#'
#'   Schwager, P., Berg, C. (2021) Remote sensing variables improve species
#'   distribution models for alpine plant species. *Basic and Applied
#'   Ecology* 54, 1-13. \doi{10.1016/j.baae.2021.04.002}
#'
#'   Weiss, A.D. (2001) Topographic position and landforms analysis. Poster,
#'   ESRI User Conference, San Diego.
#' @examples
#' # a small synthetic hill
#' dem <- terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 4000,
#'                    ymin = 0, ymax = 4000, crs = "EPSG:32633")
#' xy <- terra::xyFromCell(dem, 1:terra::ncell(dem))
#' terra::values(dem) <- 500 + 300 * exp(-((xy[, 1] - 2000)^2 +
#'                                         (xy[, 2] - 2000)^2) / 1e6)
#' ti <- topo_indices(dem, radius = c(3, 8))
#' names(ti)
#' terra::plot(ti[[c("northness", "tpi_8", "twi", "hli")]])
#' @export
topo_indices <- function(dem, indices = "all", radius = c(3, 10),
                         radius_unit = c("cells", "km"), min_slope = 0.1,
                         max_cells = 2e7, path = NULL, prefix = "",
                         overwrite = FALSE) {
  radius_unit <- match.arg(radius_unit)
  all_ix <- c("slope", "northness", "eastness", "tri", "tpi", "dev", "twi",
              "hli")
  if (identical(indices, "all") || identical(indices, "terrain")) {
    indices <- all_ix
  }
  bad <- setdiff(indices, all_ix)
  if (length(bad)) {
    stop("Unknown terrain index: ", paste(bad, collapse = ", "),
         ". Choose from ", paste(all_ix, collapse = ", "), ".", call. = FALSE)
  }
  if (is.character(dem)) dem <- terra::rast(dem)
  if (!inherits(dem, "SpatRaster") || terra::nlyr(dem) != 1) {
    stop("`dem` must be a single-layer SpatRaster or a file path.",
         call. = FALSE)
  }
  if (is.na(terra::crs(dem)) || !nzchar(terra::crs(dem))) {
    stop("`dem` has no CRS: slope and areas need one.", call. = FALSE)
  }
  if (any(c("tpi", "dev") %in% indices) &&
      (!length(radius) || any(!is.finite(radius)) || any(radius <= 0))) {
    stop("`radius` must be positive numbers.", call. = FALSE)
  }
  if ("twi" %in% indices && terra::ncell(dem) > max_cells) {
    stop("The TWI needs the whole DEM in memory and it has ",
         format(terra::ncell(dem), big.mark = ","), " cells (max_cells = ",
         format(max_cells, big.mark = ","), "). Clip it to the study area, ",
         "aggregate it, or raise max_cells.", call. = FALSE)
  }
  if (!is.null(path) && !dir.exists(path)) dir.create(path, recursive = TRUE)
  out <- list()
  make <- function(name, r) {
    names(r) <- name
    if (!is.null(path)) {
      f <- file.path(path, paste0(prefix, name, ".tif"))
      if (file.exists(f) && !overwrite) {
        stop("File exists: ", f, " (use overwrite = TRUE).", call. = FALSE)
      }
      r <- terra::writeRaster(r, f, overwrite = TRUE,
                              names = name)
    }
    out[[name]] <<- r
  }
  deg <- pi / 180
  need_slope <- any(c("slope", "northness", "eastness", "twi", "hli") %in%
                      indices)
  if (need_slope) slp <- terra::terrain(dem, "slope", unit = "degrees")
  if (any(c("northness", "eastness", "hli") %in% indices)) {
    asp <- terra::terrain(dem, "aspect", unit = "degrees")
    asp <- terra::ifel(is.na(asp) & !is.na(slp), 0, asp)
    # flat ground (slope 0 up to rounding): no aspect
    sin_s <- terra::ifel(slp < 1e-6, 0, sin(slp * deg))
  }
  if ("slope" %in% indices) make("slope", slp)
  if ("northness" %in% indices) {
    make("northness", cos(asp * deg) * sin_s)
  }
  if ("eastness" %in% indices) {
    make("eastness", sin(asp * deg) * sin_s)
  }
  if ("tri" %in% indices) make("tri", terra::terrain(dem, "TRIriley"))
  if (any(c("tpi", "dev") %in% indices)) {
    lab <- if (radius_unit == "km") paste0(format_radius(radius), "km") else
      format_radius(radius)
    for (k in seq_along(radius)) {
      w <- topo_window(dem, radius[k], radius_unit)
      mu <- terra::focal(dem, w = w, fun = "mean", na.rm = TRUE)
      tpi <- dem - mu
      if ("tpi" %in% indices) make(paste0("tpi_", lab[k]), tpi)
      if ("dev" %in% indices) {
        s <- terra::focal(dem, w = w, fun = "sd", na.rm = TRUE)
        make(paste0("dev_", lab[k]),
             terra::ifel(s > 0, tpi / s, NA))
      }
    }
  }
  if ("twi" %in% indices) make("twi", topo_twi(dem, slp, min_slope))
  if ("hli" %in% indices) {
    lat <- topo_latitude(dem)
    north <- lat >= 0
    fold <- terra::ifel(north, abs(180 - abs(asp - 225)),
                        abs(180 - abs(asp - 315)))
    L <- abs(lat) * deg
    S <- slp * deg
    A <- fold * deg
    make("hli", exp(-1.467 + 1.582 * cos(L) * cos(S) -
                      1.5 * cos(A) * sin(S) * sin(L) -
                      0.262 * sin(L) * sin(S) +
                      0.607 * sin(A) * sin(S)))
  }
  terra::rast(out)
}

# Radius labels without trailing zeros (internal)
format_radius <- function(r) {
  sub("\\.$", "", sub("0+$", "", formatC(r, format = "f", digits = 3)))
}

# Circular focal window, centre excluded; NA outside the circle (internal)
topo_window <- function(dem, r, unit) {
  if (unit == "cells") {
    rx <- ry <- r
  } else {
    rs <- terra::res(dem)
    if (terra::is.lonlat(dem)) {
      lat <- mean(as.vector(terra::ext(dem))[3:4]) * pi / 180
      ry <- r / (rs[2] * 111.32)
      rx <- r / (rs[1] * 111.32 * max(cos(lat), 0.01))
    } else {
      ry <- r * 1000 / rs[2]
      rx <- r * 1000 / rs[1]
    }
  }
  nx <- max(1, floor(rx))
  ny <- max(1, floor(ry))
  gx <- matrix(rep(-nx:nx, each = 2 * ny + 1), 2 * ny + 1)
  gy <- matrix(rep(-ny:ny, times = 2 * nx + 1), 2 * ny + 1)
  w <- ifelse((gx / max(rx, 1))^2 + (gy / max(ry, 1))^2 <= 1 + 1e-9, 1, NA)
  w[ny + 1, nx + 1] <- NA
  w
}

# Latitude of each cell, from a coarse grid for projected DEMs (internal)
topo_latitude <- function(dem) {
  if (terra::is.lonlat(dem)) return(terra::init(dem, "y"))
  fac <- ceiling(max(dim(dem)[1:2]) / 200)
  g <- if (fac > 1) terra::aggregate(terra::rast(dem), fac) else
    terra::rast(dem)
  xy <- terra::xyFromCell(g, seq_len(terra::ncell(g)))
  ll <- terra::project(xy, from = terra::crs(dem), to = "EPSG:4326")
  terra::values(g) <- ll[, 2]
  if (fac > 1) terra::resample(g, dem, method = "bilinear") else g
}

# Topographic wetness index with D8 flow accumulation (internal)
topo_twi <- function(dem, slp, min_slope) {
  fd <- terra::values(terra::terrain(dem, "flowdir"), mat = FALSE)
  z <- terra::values(dem, mat = FALSE)
  nr <- terra::nrow(dem)
  nc <- terra::ncol(dem)
  n <- nr * nc
  # true areas on longitude/latitude; map units (as slope) when projected
  area <- if (terra::is.lonlat(dem)) {
    terra::values(terra::cellSize(dem, unit = "m"), mat = FALSE)
  } else {
    rep(prod(terra::res(dem)), n)
  }
  # ESRI D8 codes: 1 E, 2 SE, 4 S, 8 SW, 16 W, 32 NW, 64 N, 128 NE
  codes <- c(1, 2, 4, 8, 16, 32, 64, 128)
  dr <- c(0, 1, 1, 1, 0, -1, -1, -1)
  dc <- c(1, 1, 0, -1, -1, -1, 0, 1)
  k <- match(fd, codes)
  id <- seq_len(n)
  row <- (id - 1) %/% nc + 1
  col <- (id - 1) %% nc + 1
  tr <- row + dr[k]
  tc <- col + dc[k]
  to <- (tr - 1) * nc + tc
  ok <- !is.na(k) & !is.na(z) & tr >= 1 & tr <= nr & tc >= 1 & tc <= nc
  to[!ok] <- NA
  to[ok][is.na(z[to[ok]])] <- NA
  acc <- ifelse(is.na(z), NA, area)
  indeg <- tabulate(to[!is.na(to)], nbins = n)
  todo <- which(indeg == 0 & !is.na(z))
  while (length(todo)) {
    t <- to[todo]
    keep <- !is.na(t)
    if (!any(keep)) break
    s <- rowsum(acc[todo[keep]], t[keep], reorder = FALSE)
    tgt <- as.integer(rownames(s))
    acc[tgt] <- acc[tgt] + s[, 1]
    indeg[tgt] <- indeg[tgt] - tabulate(match(t[keep], tgt),
                                        nbins = length(tgt))
    todo <- tgt[indeg[tgt] == 0]
  }
  sca <- acc / sqrt(area)
  b <- pmax(terra::values(slp, mat = FALSE), min_slope) * pi / 180
  out <- terra::rast(dem)
  terra::values(out) <- log(sca / tan(b))
  out
}
