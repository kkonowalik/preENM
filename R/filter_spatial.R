# Spatial helpers for filter_occurrences() (internal) -------------------------

# Values of a raster at lon/lat points (points projected to the raster's CRS)
extract_at <- function(r, lon, lat) {
  if (!length(lon)) return(numeric(0))
  p <- terra::vect(cbind(lon, lat), crs = "EPSG:4326")
  if (!terra::same.crs(p, r)) p <- terra::project(p, r)
  terra::extract(r, p, ID = FALSE)[[1]]
}

# Stated-country check with CoordinateCleaner::cc_coun()
run_country_check <- function(country, lon, lat, use, st, flag, skip) {
  if (is.null(country)) return(skip("country_mismatch", "no country column"))
  for (pk in c("CoordinateCleaner", "countrycode")) {
    if (!requireNamespace(pk, quietly = TRUE)) {
      return(skip("country_mismatch", paste("install the", pk, "package")))
    }
  }
  iso3 <- suppressWarnings(countrycode::countrycode(country, "iso2c", "iso3c",
                                                    warn = FALSE))
  ok <- use & !is.na(iso3) & !(country %in% st$country_check_skip)
  if (!any(ok)) return(flag(rep(FALSE, length(lon)), "country_mismatch"))
  d <- data.frame(decimalLongitude = lon[ok], decimalLatitude = lat[ok],
                  iso3 = iso3[ok])
  # Natural Earth codes France and Norway as -99 in iso_a3; iso_a3_eh is
  # complete
  args <- list(x = d, lon = "decimalLongitude", lat = "decimalLatitude",
               iso3 = "iso3", value = "flagged", ref_col = "iso_a3_eh",
               verbose = FALSE)
  if ("buffer" %in% names(formals(CoordinateCleaner::cc_coun)) &&
      st$country_buffer_m > 0) {
    args$buffer <- st$country_buffer_m
  }
  pass <- tryCatch(do.call(CoordinateCleaner::cc_coun, args),
                   error = function(e) e)
  if (inherits(pass, "error")) {
    return(skip("country_mismatch", paste("CoordinateCleaner failed:",
                                          conditionMessage(pass))))
  }
  bad <- rep(FALSE, length(lon))
  bad[which(ok)[!pass]] <- TRUE
  flag(bad, "country_mismatch")
}

# CoordinateCleaner::clean_coordinates() on the records still alive
run_cc <- function(lon, lat, use, st, urban, flag, skip, say) {
  if (!requireNamespace("CoordinateCleaner", quietly = TRUE)) {
    return(skip("coordinatecleaner", "install the CoordinateCleaner package"))
  }
  tests <- st$cc_tests
  if (!length(tests)) return(skip("coordinatecleaner", "cc_tests is empty"))
  if (!any(use)) return(flag(rep(FALSE, length(lon)), "coordinatecleaner"))
  if ("urban" %in% tests && is.null(urban)) {
    urban <- tryCatch({
      if (!requireNamespace("rnaturalearth", quietly = TRUE)) stop("no pkg")
      rnaturalearth::ne_download(scale = 10, type = "urban_areas",
                                 category = "cultural", returnclass = "sf")
    }, error = function(e) NULL)
    if (is.null(urban)) {
      warning("Urban areas could not be downloaded (offline?): the ",
              "CoordinateCleaner 'urban' test was skipped. Supply them with ",
              "filter_layers(urban = ...).", call. = FALSE)
      tests <- setdiff(tests, "urban")
    }
  }
  d <- data.frame(species = "sp", decimalLongitude = lon[use],
                  decimalLatitude = lat[use])
  cc <- function(urb) {
    CoordinateCleaner::clean_coordinates(
      x = d, lon = "decimalLongitude", lat = "decimalLatitude",
      species = "species", tests = tests, urban_ref = urb,
      centroids_rad = st$centroids_rad, value = "flagged", verbose = FALSE)
  }
  pass <- tryCatch(cc(urban), error = function(e) {
    if (inherits(urban, "sf")) {
      tryCatch(cc(terra::vect(urban)), error = function(e2) e2)
    } else e
  })
  if (inherits(pass, "error")) {
    return(skip("coordinatecleaner", paste("CoordinateCleaner failed:",
                                           conditionMessage(pass))))
  }
  say("  (CoordinateCleaner tests: ", paste(tests, collapse = ", "), ")")
  bad <- rep(FALSE, length(lon))
  bad[which(use)[!pass]] <- TRUE
  flag(bad, "coordinatecleaner")
}

# Projection for buffering: polar for circumpolar data, else centred on data
buffer_crs <- function(lon, lat, rng, setting) {
  if (!identical(setting, "auto")) return(setting)
  e <- terra::ext(terra::project(rng, "EPSG:4326"))
  span <- max(diff(range(lon, na.rm = TRUE)), e$xmax - e$xmin)
  if (span > 180) {
    pole <- if (mean(lat, na.rm = TRUE) >= 0) 90 else -90
    return(sprintf(paste0("+proj=aeqd +lat_0=%d +lon_0=0 +x_0=0 +y_0=0 ",
                          "+datum=WGS84 +units=m +no_defs"), pole))
  }
  sprintf(paste0("+proj=aeqd +lat_0=%.4f +lon_0=%.4f +x_0=0 +y_0=0 ",
                 "+datum=WGS84 +units=m +no_defs"),
          mean(lat, na.rm = TRUE), mean(lon, na.rm = TRUE))
}

# Inside the range polygon (core) and inside polygon + buffer (zone)
range_membership <- function(lon, lat, valid, rng, st) {
  n <- length(lon)
  core <- zone <- rep(NA, n)
  crs <- buffer_crs(lon[valid], lat[valid], rng, st$range_crs)
  rp <- terra::makeValid(terra::project(terra::makeValid(rng), crs))
  rp <- terra::aggregate(rp)
  zp <- if (st$range_buffer_km > 0) {
    terra::aggregate(terra::buffer(rp, st$range_buffer_km * 1000))
  } else rp
  p <- terra::project(terra::vect(cbind(lon[valid], lat[valid]),
                                  crs = "EPSG:4326"), crs)
  core[valid] <- as.vector(terra::is.related(p, rp, "intersects"))
  zone[valid] <- as.vector(terra::is.related(p, zp, "intersects"))
  if (any(core & !zone, na.rm = TRUE)) {
    stop(sum(core & !zone, na.rm = TRUE), " records are inside the range ",
         "polygon but outside polygon + buffer: the geometry is broken. ",
         "Try another `range_crs` or repair the polygon.", call. = FALSE)
  }
  list(core = core, zone = zone)
}

# Sea / open water: distance to land polygons, with a water-raster fallback
sea_test <- function(lon, lat, unc, use, layers, st, say) {
  n <- length(lon)
  land_max <- water_max <- dist <- rep(NA_real_, n)
  if (!is.null(layers$water)) {
    w <- water_window(layers$water, lon, lat, unc, use, st)
    land_max <- w$land_max
    water_max <- w$water_max
  }
  if (!is.null(layers$land)) {
    cand <- which(use & (is.na(water_max) | water_max >= 1))
    say("  records checked against the land polygons: ", length(cand))
    if (length(cand)) dist[cand] <- distance_to_land(lon[cand], lat[cand],
                                                     layers$land)
  }
  list(land_max = land_max, water_max = water_max, dist = dist)
}

# Distance (m) from points to land polygons, processed in ~1 degree blocks;
# 0 on land, Inf when no land lies within the `pad_km` search window
distance_to_land <- function(lon, lat, land, pad_km = 25) {
  if (!terra::is.lonlat(land)) land <- terra::project(land, "EPSG:4326")
  d <- rep(NA_real_, length(lon))
  blk <- paste(floor(lon), floor(lat))
  pad_lat <- pad_km / 111.2
  for (b in unique(blk)) {
    s <- which(blk == b)
    pad_lon <- pad_km / (111.2 * max(cos(mean(lat[s]) * pi / 180), 0.05))
    e <- terra::ext(max(-180, min(lon[s]) - pad_lon),
                    min(180, max(lon[s]) + pad_lon),
                    max(-90, min(lat[s]) - pad_lat),
                    min(90, max(lat[s]) + pad_lat))
    lb <- tryCatch(terra::crop(land, e), error = function(err) NULL)
    if (is.null(lb)) next
    if (nrow(lb) == 0) {
      d[s] <- Inf
      next
    }
    # distances in a local equidistant projection centred on the block:
    # exact enough within the search window, and independent of how a terra
    # version measures distances in longitude/latitude
    loc <- sprintf(paste0("+proj=aeqd +lat_0=%.5f +lon_0=%.5f +x_0=0 +y_0=0 ",
                          "+datum=WGS84 +units=m +no_defs"),
                   mean(lat[s]), mean(lon[s]))
    p <- terra::project(terra::vect(cbind(lon[s], lat[s]),
                                    crs = "EPSG:4326"), loc)
    lp <- terra::project(lb, loc)
    dm <- terra::distance(p, lp)
    d[s] <- apply(matrix(dm, nrow = length(s)), 1, min)
  }
  d
}

# Water raster: most land (%) in a window sized by the uncertainty, and the
# most water (%) in the 5x5 cells around (coastal / lake candidates)
water_window <- function(wat, lon, lat, unc, use, st) {
  n <- length(lon)
  land_max <- water_max <- rep(NA_real_, n)
  idx <- which(use)
  if (!length(idx)) return(list(land_max = land_max, water_max = water_max))
  p <- terra::vect(cbind(lon[idx], lat[idx]), crs = "EPSG:4326")
  if (!terra::same.crs(p, wat)) p <- terra::project(p, wat)
  cells <- terra::cellFromXY(wat, terra::crds(p))
  rc <- terra::rowColFromCell(wat, cells)
  vals <- function(cl) {
    out <- rep(NA_real_, length(cl))
    ok <- !is.na(cl)
    if (any(ok)) out[ok] <- terra::extract(wat, cl[ok])[[1]]
    out
  }
  own <- vals(cells)
  mx <- suppressWarnings(max(own, na.rm = TRUE))
  sc <- if (is.finite(mx) && mx <= 1) 100 else 1
  u <- unc[idx]
  k <- ifelse(is.na(u), 1L, ifelse(u <= st$sea_ring_unc_m[1], 0L,
              ifelse(u <= st$sea_ring_unc_m[2], 1L, 2L)))
  lm <- wm <- rep(NA_real_, length(idx))
  for (dr in -2:2) for (dc in -2:2) {
    cl <- terra::cellFromRowCol(wat, rc[, 1] + dr, rc[, 2] + dc)
    w <- vals(cl) * sc
    wm <- pmax(wm, w, na.rm = TRUE)
    inside <- !is.na(cells) & pmax(abs(dr), abs(dc)) <= k
    lm[inside] <- pmax(lm[inside], 100 - w[inside], na.rm = TRUE)
  }
  land_max[idx] <- lm
  water_max[idx] <- wm
  list(land_max = land_max, water_max = water_max)
}

# Lowest (or highest) raster value within radius_km of a point
window_extreme <- function(r, lon, lat, radius_km, fun = min) {
  p <- terra::vect(cbind(lon, lat), crs = "EPSG:4326")
  if (terra::is.lonlat(r)) {
    cl <- cos(lat * pi / 180)
    dlat <- radius_km / 111.2
    dlon <- radius_km / (111.2 * max(cl, 0.05))
    e <- terra::ext(max(lon - dlon, -180), min(lon + dlon, 180),
                    max(lat - dlat, -90), min(lat + dlat, 90))
    dist_fun <- function(cc) {
      sqrt(((cc[, 1] - lon) * cl)^2 + (cc[, 2] - lat)^2) * 111.2
    }
  } else {
    xy <- terra::crds(terra::project(p, r))
    rm <- radius_km * 1000
    e <- terra::ext(xy[1] - rm, xy[1] + rm, xy[2] - rm, xy[2] + rm)
    dist_fun <- function(cc) {
      sqrt((cc[, 1] - xy[1])^2 + (cc[, 2] - xy[2])^2) / 1000
    }
  }
  w <- tryCatch(terra::crop(r, e), error = function(err) NULL)
  if (is.null(w)) return(NA_real_)
  v <- terra::values(w, mat = FALSE)
  cc <- terra::xyFromCell(w, seq_len(terra::ncell(w)))
  v <- v[dist_fun(cc) <= radius_km & !is.na(v)]
  if (length(v)) fun(v) else NA_real_
}

# Points on a sphere as 3-D coordinates in km (chord distances, no projection)
to_xyz <- function(lon, lat, R = 6371) {
  lon <- lon * pi / 180
  lat <- lat * pi / 180
  cbind(R * cos(lat) * cos(lon), R * cos(lat) * sin(lon), R * sin(lat))
}

# Number of other points within eps_km
neighbour_counts <- function(lon, lat, eps_km) {
  if (length(lon) < 2) return(rep(0L, length(lon)))
  lengths(dbscan::frNN(to_xyz(lon, lat), eps = eps_km)$id)
}

# Size of the cluster each point belongs to (points chained closer than eps)
cluster_sizes <- function(lon, lat, eps_km) {
  if (length(lon) < 2) return(rep(1L, length(lon)))
  cl <- dbscan::dbscan(to_xyz(lon, lat), eps = eps_km, minPts = 2)$cluster
  id <- ifelse(cl == 0, -seq_along(cl), cl)     # noise = own cluster of 1
  as.vector(table(id)[as.character(id)])
}
