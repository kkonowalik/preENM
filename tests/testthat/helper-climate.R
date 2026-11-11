# Synthetic bio-like climate for a present, a second present source and an
# LGM on one grid (Europe-like extent), built from three latent gradients
make_climate <- function(seed = 1, nr = 60, nc = 90) {
  set.seed(seed)
  r <- terra::rast(nrows = nr, ncols = nc, xmin = -10, xmax = 35, ymin = 40,
                   ymax = 70, crs = "EPSG:4326")
  xy <- terra::crds(r, na.rm = FALSE)
  lon <- (xy[, 1] + 10) / 45
  lat <- (xy[, 2] - 40) / 30
  topo <- sin(lon * 9) * cos(lat * 7) + stats::rnorm(nrow(xy), 0, 0.3)
  build <- function(Tm, Co, Pr) {
    v <- cbind(
      bio01 = Tm,
      bio04 = Co,
      bio05 = Tm + 0.6 * Co + stats::rnorm(nrow(xy), 0, 0.05),
      bio06 = Tm - 0.8 * Co + stats::rnorm(nrow(xy), 0, 0.05),
      bio10 = Tm + 0.4 * Co + stats::rnorm(nrow(xy), 0, 0.02),
      bio11 = Tm - 0.5 * Co + stats::rnorm(nrow(xy), 0, 0.02),
      bio12 = Pr,
      bio14 = 0.4 * Pr - 0.6 * Co + stats::rnorm(nrow(xy), 0, 0.1),
      bio15 = Co + 0.5 * stats::rnorm(nrow(xy)),
      scd = -Tm + 0.3 * topo + stats::rnorm(nrow(xy), 0, 0.1))
    out <- terra::rast(r, nlyrs = ncol(v))
    terra::values(out) <- v
    names(out) <- colnames(v)
    out
  }
  Tm <- 1 - lat - 0.5 * topo
  Co <- lon + 0.2 * stats::rnorm(nrow(xy))
  Pr <- 1 - lon + 0.5 * topo + 0.2 * stats::rnorm(nrow(xy))
  present <- build(Tm, Co, Pr)
  trace <- build(Tm + 0.1, Co * 1.05, Pr + stats::rnorm(nrow(xy), 0, 0.1))
  lgm <- build(Tm - 1.2 - 0.6 * lat, Co * 1.6, Pr * 0.6)
  # ice sheet: no values in the north-west in the LGM
  ice <- xy[, 2] > 62 & xy[, 1] < 20
  lgm[ice] <- NA
  regions <- sf::st_as_sf(data.frame(
    name = c("Alps", "Scandinavia"),
    wkt = c("POLYGON ((5 44, 17 44, 17 49, 5 49, 5 44))",
            "POLYGON ((5 58, 30 58, 30 70, 5 70, 5 58))")),
    wkt = "wkt", crs = 4326)
  list(present = present, trace = trace, lgm = lgm, regions = regions)
}
