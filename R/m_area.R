# Accessible area (M) delimitation ---------------------------------------------

#' Delimit the accessible area (M) from occurrence records
#'
#' Draws the area the species could have reached (the M of the BAM diagram),
#' which defines where background points are drawn and to which the
#' predictors are clipped. Instead of choosing buffer radii by hand, one
#' number is derived from the data: the **linking distance L**.
#'
#' **Linking distance.** After thinning (one record per `thin_km` cell), the
#' records are joined by single-linkage clustering, i.e. the minimum spanning
#' tree. Distances between neighbouring records within a continuous range are
#' sampling gaps; the gaps between separate land masses or isolated parts of
#' the range are much longer. L is the last link before the largest relative
#' jump in the sorted spanning-tree links between `link_search[1]` and
#' `link_search[2]` km. A second distance L2 is used for corridors: the link
#' just after the largest jump among the links longer than L (up to
#' `link2_max_km`), or the longest of them when there are only one or two.
#' See [linking_distance()].
#'
#' Whatever the method, compare it with others in [m_scenarios()]: the
#' calibration area strongly affects model accuracy (Luna et al. 2024).
#'
#' **Methods** (radius r = L/2 by default: two buffers of L/2 touch when
#' records are L apart, so every chain of records linked at L becomes one
#' patch):
#' * `"delaunay"` (default): Delaunay triangles whose edges are all at most L,
#'   plus buffers of r around every record. Follows concave range outlines.
#' * `"buffer"`: buffers of radius `buffer_km` (default r) around every
#'   record.
#' * `"cluster_hull"`: for each cluster of records linked at L, the convex
#'   hull (`hull_ratio = 1`) or a concave hull (`hull_ratio < 1`) of the
#'   buffers. Fills gaps inside each cluster automatically.
#' * `"convex_hull"`: one convex hull around all buffers (the simplest and
#'   most generous M).
#' * `"simulation"`: M simulated by dispersal from the records over the
#'   landscape in `env` ([m_simulate()], Machado-Stredel et al. 2021); give
#'   its settings in `simulation`, e.g. `simulation = list(kernel_sd = 2,
#'   n_events = 30, barrier = mountains)`. L is still computed and
#'   reported but does not shape this M.
#' * `"adaptive"`: spatially adaptive buffers for unevenly surveyed areas
#'   (Barbosa 2026, "cluster mean distance"): around each cluster of records
#'   linked at L, a buffer whose radius is the mean distance between the
#'   records of that cluster (times `adaptive_mult`), so well-surveyed,
#'   extensive clusters reach further and isolated records keep only a
#'   small circle (`adaptive_min_km`) instead of pulling unsurveyed land
#'   into M. With `adaptive_weight = TRUE` each radius is further scaled by
#'   the cluster's number of records relative to the largest cluster.
#'   One large, extensive cluster gets a wide buffer (its mean distance can
#'   be several times L), so this method is generous for well-surveyed
#'   ranges and conservative for scattered records.
#'
#' **Steps after the method**, in order: corridors between clusters
#' (`corridors = TRUE`: the spanning-tree links between L and L2, buffered by
#' r), a closing (buffer out and back in by `closing_km`; fills bays and gaps
#' narrower than twice that), filling of interior holes, the expert additions
#' (`add` polygons and `add_sites`), sealing slivers between them
#' (`seal_km`), filling holes again, and finally clipping to `extent`, removing
#' `exclude` and trimming to `land`.
#'
#' **Expert intervention.** Everything automatic can be overridden: set
#' `link_km` by hand, add polygons drawn in a GIS (e.g. a valley or migration
#' route the automatic M leaves out) with `add`, add sites known from the
#' literature with `add_sites`, and remove areas with `exclude`. The report
#' states how much area each addition contributes, so it can be reported as
#' expert-based in the Methods.
#'
#' The default (`"delaunay"`, closed by L/2, holes filled) is the M selected
#' for *Dryas octopetala* in Europe, where L was 340 km.
#'
#' @param points Occurrence records: a data.frame with coordinates (lon/lat or
#'   detected column names), a CSV path, an `sf` object or a
#'   [terra::SpatVector]. Combine modern and fossil records beforehand if both
#'   should shape M (e.g. with [merge_occurrences()]).
#' @param method `"delaunay"`, `"buffer"`, `"cluster_hull"`,
#'   `"convex_hull"`, `"adaptive"` or `"simulation"`.
#' @param env,simulation For `"simulation"`: the predictor layers and a
#'   list of further arguments for [m_simulate()].
#' @param link_km Linking distance L in km; `NULL` (default) = derived from the
#'   data.
#' @param link_search Range (km) in which the automatic L is searched.
#' @param link2_max_km Upper limit (km) for the second linking distance L2.
#' @param link_mult Multiplier of L for the method (sensitivity analyses, e.g.
#'   0.5 or 2).
#' @param buffer_km Buffer radius for `method = "buffer"` (default L/2).
#' @param adaptive_mult,adaptive_weight,adaptive_min_km For `"adaptive"`:
#'   multiplier of the cluster radii, scaling by cluster size, and the
#'   radius (km) of single-record clusters.
#' @param hull_ratio For `"cluster_hull"`: 1 = convex, smaller = more concave
#'   (needs GEOS >= 3.11).
#' @param corridors Join clusters lying within L2 of each other by corridors.
#' @param closing_km Closing radius (km). Default: L/2 for `"delaunay"` (or
#'   L2/2 with corridors), 0 otherwise.
#' @param close_km Small closing (km) applied to the raw buffers to smooth
#'   them (default 25).
#' @param fill_holes Fill interior holes.
#' @param add Expert polygons to add (`sf`, [terra::SpatVector] or a path,
#'   any CRS).
#' @param add_sites Expert sites to add: a data.frame with `lon`, `lat` and
#'   `buffer_km` (and optionally `name`).
#' @param seal_km Closing (km) after the additions, sealing slivers where an
#'   added polygon almost touches the automatic M.
#' @param exclude Polygons to remove from M (e.g. another mountain system).
#'   Records inside them are dropped before M is built.
#' @param extent Optional `c(xmin, xmax, ymin, ymax)` in longitude/latitude:
#'   records outside are dropped and M is clipped to it.
#' @param land Optional land polygons to trim the sea from M. Not needed for
#'   clipping rasters, whose sea cells are usually `NA` anyway.
#' @param thin_km Thinning cell size (km) before linking; one record per cell
#'   (and per `group`).
#' @param group Optional column (e.g. modern / fossil) thinned separately and
#'   reported separately.
#' @param priority,priority_col Sources preferred when thinning (see
#'   [thin_by_cell()]).
#' @param crs Equal-area projection used for all geometry. `"auto"`: Lambert
#'   azimuthal equal-area centred on the data (polar for circumpolar data).
#' @param coords,points_crs Coordinate columns and their CRS when `points` is
#'   a table.
#' @return An object of class `m_area`: a list with `m` (the final M, `sfc`
#'   in `crs`), `automatic` (before expert additions), `corridors`, `points`
#'   (thinned records with their cluster), `all_points` (all records used),
#'   `linking` (see [linking_distance()]), `report` and `settings`. Use
#'   [plot_m_area()], [m_report()], [write_m_area()] and [crop_to_m()].
#' @references Barbosa, A.M. (2026) Spatially adaptive background
#'   delimitation strategies for species distribution models with unevenly
#'   surveyed occurrences. *Biodiversity Informatics* 20, 81-88.
#'   \doi{10.17161/bi.v20i1.25045}
#'
#'   Machado-Stredel, F., Cobos, M.E., Peterson, A.T. (2021) A
#'   simulation-based method for selecting calibration areas for ecological
#'   niche models and species distribution models. *Frontiers of
#'   Biogeography* 13, e48814. \doi{10.21425/F5FBG48814}
#'
#'   Luna, S., Pena-Peniche, A., Mendoza-Alfaro, R. (2024) Species
#'   distribution model accuracy is strongly influenced by the choice of
#'   calibration area. *Biodiversity Informatics* 18, 43-55.
#'   \doi{10.17161/bi.v18i.22655}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' m
#' plot_m_area(m)
#' # an expert polygon added to the automatic M
#' extra <- sf::st_sfc(sf::st_polygon(list(rbind(c(22, 50), c(24, 50),
#'   c(24, 51), c(22, 51), c(22, 50)))), crs = 4326)
#' m2 <- delimit_m(sim$presence, link_km = 120, thin_km = 5, add = extra)
#' m_report(m2)
delimit_m <- function(points, method = c("delaunay", "buffer", "cluster_hull",
                                         "convex_hull", "adaptive",
                                         "simulation"),
                      link_km = NULL, link_search = c(100, 600),
                      link2_max_km = 1000, link_mult = 1, buffer_km = NULL,
                      adaptive_mult = 1, adaptive_weight = FALSE,
                      adaptive_min_km = 1, env = NULL, simulation = list(),
                      hull_ratio = 1, corridors = FALSE, closing_km = NULL,
                      close_km = 25, fill_holes = TRUE, add = NULL,
                      add_sites = NULL, seal_km = 5, exclude = NULL,
                      extent = NULL, land = NULL, thin_km = 10, group = NULL,
                      priority = NULL, priority_col = "source",
                      crs = "auto", coords = NULL,
                      points_crs = "EPSG:4326") {
  method <- match.arg(method)
  prep <- m_prepare(points, extent = extent, exclude = exclude,
                    thin_km = thin_km, group = group, priority = priority,
                    priority_col = priority_col, crs = crs, coords = coords,
                    points_crs = points_crs)
  lk <- m_linking(prep, link_km = link_km, link_search = link_search,
                  link2_max_km = link2_max_km)
  sim <- NULL
  if (method == "simulation") {
    if (is.null(env)) {
      stop("method = \"simulation\" needs `env` (see m_simulate()).",
           call. = FALSE)
    }
    sim <- do.call(m_simulate, c(list(points = prep$all[c("lon", "lat")],
                                      env = env), simulation))
  }
  out <- m_build(prep, lk, method = method, link_mult = link_mult, sim = sim,
          buffer_km = buffer_km, adaptive = list(
            mult = adaptive_mult, weight = adaptive_weight,
            min_km = adaptive_min_km), hull_ratio = hull_ratio,
          corridors = corridors, closing_km = closing_km,
          close_km = close_km, fill_holes = fill_holes, add = add,
          add_sites = add_sites, seal_km = seal_km, land = land)
  if (!is.null(sim)) out$simulation <- sim
  out
}

#' Linking distance of occurrence records
#'
#' The connectivity analysis behind [delimit_m()]: thins the records,
#' computes their minimum spanning tree (single-linkage clustering) and
#' derives the linking distance L and the second distance L2. Use it, with
#' [plot_linking()], to see how the number of separate clusters changes with
#' the linking distance before choosing or accepting L.
#'
#' The spanning tree is computed from the Delaunay triangulation, which
#' contains it, so memory use grows linearly with the number of records
#' (a full distance matrix would not fit for tens of thousands of records).
#'
#' @inheritParams delimit_m
#' @return An object of class `m_linking`: a list with `L` and `L2` (km),
#'   `auto` (whether L was derived from the data), `n_clusters` (at L),
#'   `table` (number of clusters, largest cluster and singletons for a range
#'   of linking distances), `links` (the spanning-tree edges, km) and the
#'   thinned points.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' lk <- linking_distance(sim$presence, thin_km = 5, link_search = c(20, 200))
#' lk
#' plot_linking(lk)
linking_distance <- function(points, link_km = NULL,
                             link_search = c(100, 600), link2_max_km = 1000,
                             thin_km = 10, group = NULL, priority = NULL,
                             priority_col = "source", extent = NULL,
                             exclude = NULL, crs = "auto", coords = NULL,
                             points_crs = "EPSG:4326") {
  prep <- m_prepare(points, extent = extent, exclude = exclude,
                    thin_km = thin_km, group = group, priority = priority,
                    priority_col = priority_col, crs = crs, coords = coords,
                    points_crs = points_crs)
  m_linking(prep, link_km = link_km, link_search = link_search,
            link2_max_km = link2_max_km)
}

#' Compare several M scenarios
#'
#' Builds several M variants from the same records (thinning and linking
#' distance are computed once) and compares them: area, share of records
#' inside (per `group`), and overlap with a reference scenario (share of this
#' scenario inside the reference, and the Jaccard index).
#'
#' **Use it by default, not as an optional check.** The choice of
#' calibration area was the most influential factor in simulations of model
#' accuracy, explaining half of its variance, more than regularisation,
#' collinearity treatment and spatial thinning together (Luna et al. 2024).
#' Comparing several delimitations (geometric, simulated with
#' `method = "simulation"`, and the minimal extent of [m_minimal()]) and
#' reporting how much they differ is therefore part of choosing M, not a
#' supplement to it.
#'
#' @inheritParams delimit_m
#' @param scenarios Named list; each element is a list of arguments of
#'   [delimit_m()] (method, link_mult, corridors, buffer_km, ...). The default
#'   compares the Delaunay M at 0.5, 1 and 2 times L, with corridors, the
#'   cluster and overall convex hulls, and fixed 100 km buffers.
#' @param reference Name of the reference scenario (default: the first).
#' @param ... Arguments passed to [delimit_m()] for every scenario (e.g.
#'   `add`, `exclude`, `extent`, `link_km`). A scenario with
#'   `method = "simulation"` needs `env` (and optionally `simulation`) in its
#'   own list.
#' @return An object of class `m_scenarios`: a list with `areas` (the
#'   `m_area` objects) and `comparison` (a data.frame).
#' @references Luna, S., Pena-Peniche, A., Mendoza-Alfaro, R. (2024)
#'   Species distribution model accuracy is strongly influenced by the
#'   choice of calibration area. *Biodiversity Informatics* 18, 43-55.
#'   \doi{10.17161/bi.v18i.22655}
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' sc <- m_scenarios(sim$presence, link_km = 120, thin_km = 5,
#'                   scenarios = list(delaunay = list(),
#'                                    hull = list(method = "cluster_hull"),
#'                                    buffer_50 = list(method = "buffer",
#'                                                     buffer_km = 50)))
#' sc$comparison
#' plot_m_area(sc)
m_scenarios <- function(points, scenarios = NULL, reference = NULL,
                        link_km = NULL, link_search = c(100, 600),
                        link2_max_km = 1000, thin_km = 10, group = NULL,
                        priority = NULL, priority_col = "source",
                        extent = NULL, exclude = NULL, crs = "auto",
                        coords = NULL, points_crs = "EPSG:4326", ...) {
  if (is.null(scenarios)) {
    scenarios <- list(
      delaunay = list(),
      delaunay_0.5L = list(link_mult = 0.5),
      delaunay_2L = list(link_mult = 2),
      delaunay_corridors = list(corridors = TRUE),
      cluster_hull = list(method = "cluster_hull"),
      convex_hull = list(method = "convex_hull"),
      buffer_100km = list(method = "buffer", buffer_km = 100,
                          fill_holes = FALSE))
  }
  if (is.null(names(scenarios)) || any(names(scenarios) == "")) {
    stop("`scenarios` must be a named list.", call. = FALSE)
  }
  if (is.null(reference)) reference <- names(scenarios)[1]
  if (!reference %in% names(scenarios)) {
    stop("`reference` must be one of the scenario names.", call. = FALSE)
  }
  prep <- m_prepare(points, extent = extent, exclude = exclude,
                    thin_km = thin_km, group = group, priority = priority,
                    priority_col = priority_col, crs = crs, coords = coords,
                    points_crs = points_crs)
  lk <- m_linking(prep, link_km = link_km, link_search = link_search,
                  link2_max_km = link2_max_km)
  common <- list(...)
  areas <- lapply(names(scenarios), function(nm) {
    args <- utils::modifyList(common, scenarios[[nm]])
    # a simulated M: run the simulation for this scenario
    if (identical(args$method, "simulation")) {
      if (is.null(args$env)) {
        stop("Scenario '", nm, "': method \"simulation\" needs `env`.",
             call. = FALSE)
      }
      args$sim <- do.call(m_simulate, c(list(
        points = prep$all[c("lon", "lat")], env = args$env),
        if (is.null(args$simulation)) list() else args$simulation))
    }
    args$env <- NULL
    args$simulation <- NULL
    m <- tryCatch(do.call(m_build, c(list(prep = prep, lk = lk), args)),
                  error = function(e) {
                    warning("Scenario '", nm, "' failed: ", conditionMessage(e),
                            call. = FALSE)
                    NULL
                  })
    if (!is.null(m)) m$settings$scenario <- nm
    m
  })
  names(areas) <- names(scenarios)
  areas <- Filter(Negate(is.null), areas)
  if (!reference %in% names(areas)) {
    stop("The reference scenario failed.", call. = FALSE)
  }
  ref <- areas[[reference]]$m
  a_ref <- geom_area(ref)
  cmp <- do.call(rbind, lapply(names(areas), function(nm) {
    m <- areas[[nm]]
    a <- geom_area(m$m)
    ai <- geom_area(safe_intersection(m$m, ref))
    d <- data.frame(scenario = nm, method = m$settings$method,
                    area_km2 = round(a / 1e6),
                    stringsAsFactors = FALSE)
    d <- cbind(d, inside_shares(m))
    d$share_in_reference_pct <- round(100 * ai / a, 1)
    d$jaccard_with_reference <- round(ai / (a + a_ref - ai), 3)
    d
  }))
  rownames(cmp) <- NULL
  out <- list(areas = areas, comparison = cmp, reference = reference,
              linking = lk)
  class(out) <- "m_scenarios"
  out
}

# --- preparation -----------------------------------------------------------------

m_prepare <- function(points, extent, exclude, thin_km, group, priority,
                      priority_col, crs, coords, points_crs) {
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy <- if (is.null(coords)) guess_coords(names(df)) else coords
    lon <- as_number(df[[xy[1]]])
    lat <- as_number(df[[xy[2]]])
    ok <- !is.na(lon) & !is.na(lat)
    v <- terra::vect(cbind(lon[ok], lat[ok]), crs = points_crs)
    df <- df[ok, , drop = FALSE]
  } else {
    v <- as_points(points, coords, points_crs)
    df <- as.data.frame(v)
  }
  ll <- terra::crds(terra::project(v, "EPSG:4326"))
  df$lon <- ll[, 1]
  df$lat <- ll[, 2]
  n_input <- nrow(df)
  if (!is.null(group) && !group %in% names(df)) {
    stop("Column `", group, "` (group) not found.", call. = FALSE)
  }
  df$.group <- if (is.null(group)) "all" else as.character(df[[group]])
  keep <- rep(TRUE, nrow(df))
  if (!is.null(extent)) {
    e <- unname(extent)
    keep <- keep & df$lon >= e[1] & df$lon <= e[2] & df$lat >= e[3] &
      df$lat <= e[4]
  }
  if (nrow(df) == 0 || !any(keep)) {
    stop("No records left (check coordinates and `extent`).", call. = FALSE)
  }
  if (identical(crs, "auto")) crs <- equal_area_crs(df$lon[keep],
                                                    df$lat[keep])
  ex <- NULL
  if (!is.null(exclude)) {
    ex <- union_geom(sf::st_transform(as_sfc(exclude, "exclude"), crs))
    inside_ex <- lengths(sf::st_intersects(
      sf::st_transform(sf::st_as_sf(df, coords = c("lon", "lat"), crs = 4326,
                                    remove = FALSE), crs), ex)) > 0
    keep <- keep & !inside_ex
  }
  dropped <- sum(!keep)
  df <- df[keep, , drop = FALSE]
  if (!nrow(df)) stop("No records left after `exclude`.", call. = FALSE)
  pts <- sf::st_transform(sf::st_as_sf(df, coords = c("lon", "lat"),
                                       crs = 4326, remove = FALSE), crs)
  xy <- sf::st_coordinates(pts)
  df$x <- xy[, 1]
  df$y <- xy[, 2]
  # thinning: one record per thin_km cell (and group), preferred sources first
  ord <- seq_len(nrow(df))
  if (!is.null(priority)) {
    if (!priority_col %in% names(df)) {
      stop("`priority` given but column `", priority_col, "` not found.",
           call. = FALSE)
    }
    pr <- match(df[[priority_col]], priority)
    pr[is.na(pr)] <- length(priority) + 1
    ord <- order(pr, ord)
  }
  res <- thin_km * 1000
  cell <- paste(df$.group, floor(df$x / res), floor(df$y / res))[ord]
  th <- df[ord, , drop = FALSE][!duplicated(cell), , drop = FALSE]
  # identical coordinates in different groups: keep one for the geometry
  th <- th[!duplicated(paste(th$x, th$y)), , drop = FALSE]
  rownames(th) <- NULL
  ext_geom <- NULL
  if (!is.null(extent)) {
    ext_geom <- sf::st_make_valid(sf::st_transform(box_ll(extent), crs))
  }
  list(all = df, thin = th, crs = crs, n_input = n_input,
       n_dropped = dropped, exclude = ex, extent = ext_geom,
       extent_ll = extent, thin_km = thin_km, group = group)
}

# Lambert azimuthal equal-area centred on the data (polar if circumpolar)
equal_area_crs <- function(lon, lat) {
  if (diff(range(lon)) > 180) {
    lat0 <- if (mean(lat) >= 0) 90 else -90
    lon0 <- 0
  } else {
    lat0 <- round(mean(range(lat)), 1)
    lon0 <- round(mean(range(lon)), 1)
  }
  sprintf(paste0("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 ",
                 "+datum=WGS84 +units=m +no_defs"), lat0, lon0)
}

# Densified lon/lat box as an sfc polygon
box_ll <- function(e, n = 100) {
  e <- unname(e)
  xs <- seq(e[1], e[2], length.out = n)
  ys <- seq(e[3], e[4], length.out = n)
  ring <- rbind(cbind(xs, e[3]), cbind(e[2], ys), cbind(rev(xs), e[4]),
                cbind(e[1], rev(ys)))
  ring <- rbind(ring, ring[1, ])
  sf::st_sfc(sf::st_polygon(list(unname(ring))), crs = 4326)
}

as_sfc <- function(x, what) {
  if (is.character(x)) x <- sf::st_read(x, quiet = TRUE)
  if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  if (inherits(x, "sf")) x <- sf::st_geometry(x)
  if (!inherits(x, "sfc")) {
    stop("`", what, "` must be an sf object, a SpatVector or a file path.",
         call. = FALSE)
  }
  if (is.na(sf::st_crs(x))) {
    warning("`", what, "` has no CRS; assuming EPSG:4326.", call. = FALSE)
    sf::st_crs(x) <- 4326
  }
  sf::st_make_valid(x)
}

# --- minimum spanning tree and linking distance -----------------------------------

# Unique edges of the Delaunay triangulation (indices into x, y) and the
# triangles as index triples
delaunay_graph <- function(x, y) {
  n <- length(x)
  key <- paste(x, y)
  tri <- tryCatch({
    mp <- sf::st_sfc(sf::st_multipoint(cbind(x, y)))
    suppressWarnings(sf::st_collection_extract(sf::st_triangulate(mp),
                                               "POLYGON"))
  }, error = function(e) NULL)
  if (is.null(tri) || !length(tri)) {
    # collinear points: the spanning tree joins them in order
    o <- order(x, y)
    return(list(edges = cbind(o[-n], o[-1]), tri = matrix(integer(0), 0, 3)))
  }
  cc <- sf::st_coordinates(tri)
  first3 <- stats::ave(seq_len(nrow(cc)), cc[, "L2"], FUN = seq_along) <= 3
  idx <- match(paste(cc[first3, "X"], cc[first3, "Y"]), key)
  t3 <- matrix(idx, ncol = 3, byrow = TRUE)
  t3 <- t3[stats::complete.cases(t3), , drop = FALSE]
  e <- rbind(t3[, 1:2], t3[, 2:3], t3[, c(3, 1)])
  e <- t(apply(e, 1, sort))
  e <- unique(e)
  list(edges = e, tri = t3)
}

# Kruskal on the Delaunay edges
mst_edges <- function(x, y) {
  n <- length(x)
  empty <- data.frame(from = integer(0), to = integer(0), km = numeric(0))
  if (n < 2) return(list(mst = empty, tri = matrix(integer(0), 0, 3)))
  g <- delaunay_graph(x, y)
  e <- g$edges
  len <- sqrt((x[e[, 1]] - x[e[, 2]])^2 + (y[e[, 1]] - y[e[, 2]])^2) / 1000
  o <- order(len)
  parent <- seq_len(n)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  keep <- logical(length(o))
  k <- 0L
  for (j in o) {
    a <- find(e[j, 1])
    b <- find(e[j, 2])
    if (a != b) {
      parent[a] <- b
      keep[j] <- TRUE
      k <- k + 1L
      if (k == n - 1L) break
    }
  }
  mst <- data.frame(from = e[keep, 1], to = e[keep, 2], km = len[keep])
  list(mst = mst[order(mst$km), , drop = FALSE], tri = g$tri)
}

# Cluster membership when records are linked by edges up to `km`
link_components <- function(n, mst, km) {
  parent <- seq_len(n)
  find <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  use <- mst[mst$km <= km, , drop = FALSE]
  for (j in seq_len(nrow(use))) {
    a <- find(use$from[j])
    b <- find(use$to[j])
    if (a != b) parent[a] <- b
  }
  root <- vapply(seq_len(n), find, integer(1))
  as.integer(factor(root, levels = unique(root)))
}

# Last link before the largest relative jump within [lo, hi]
auto_link <- function(h, lo, hi) {
  cand <- sort(h[h >= lo & h <= hi])
  if (length(cand) < 3) return(NA_real_)
  r <- cand[-1] / cand[-length(cand)]
  cand[which.max(r)]
}

m_linking <- function(prep, link_km, link_search, link2_max_km) {
  th <- prep$thin
  ms <- mst_edges(th$x, th$y)
  h <- ms$mst$km
  auto <- is.null(link_km)
  if (auto) {
    L <- auto_link(h, link_search[1], link_search[2])
    if (is.na(L)) {
      L <- link_search[2]
      warning("Fewer than 3 spanning-tree links between ", link_search[1],
              " and ", link_search[2], " km: L set to ", L, " km. Widen ",
              "`link_search` or set `link_km`.", call. = FALSE)
    }
  } else {
    L <- link_km
  }
  # L2: the next jump above L; with only one or two longer links within
  # link2_max_km, the longest of them (so that isolated parts within reach
  # can still be joined by corridors)
  cand2 <- sort(h[h > L & h <= link2_max_km])
  L2 <- if (length(cand2) >= 3) {
    c2 <- c(L, cand2)
    c2[which.max(c2[-1] / c2[-length(c2)]) + 1]
  } else if (length(cand2)) max(cand2) else L
  n <- nrow(th)
  grid <- sort(unique(c(10, 15, 20, 30, 50, 75, 100, 150, 200, 250, 300,
                        400, 500, 600, 800, 1000, 1500, L)))
  tab <- do.call(rbind, lapply(grid, function(s) {
    cl <- link_components(n, ms$mst, s)
    sz <- tabulate(cl)
    data.frame(link_km = round(s, 1), n_clusters = max(cl), largest = max(sz),
               singletons = sum(sz == 1))
  }))
  cl <- link_components(n, ms$mst, L)
  th$cluster <- cl
  out <- list(L = L, L2 = L2, auto = auto, link_search = link_search,
              n_clusters = max(cl), table = tab, links = ms$mst,
              tri = ms$tri, points = th, crs = prep$crs)
  class(out) <- "m_linking"
  out
}

#' @export
print.m_linking <- function(x, ...) {
  cat("Linking distance L = ", round(x$L), " km",
      if (x$auto) paste0(" (derived from the data, searched in ",
                         x$link_search[1], "-", x$link_search[2], " km)")
      else " (set by hand)", "\n", sep = "")
  cat("Second linking distance L2 = ", round(x$L2), " km\n", sep = "")
  cat(nrow(x$points), " thinned records form ", x$n_clusters,
      " cluster(s) at L\n", sep = "")
  cat("Longest spanning-tree links (km): ",
      paste(round(rev(utils::tail(x$links$km, 12))), collapse = " "),
      "\n\n", sep = "")
  print(x$table, row.names = FALSE)
  invisible(x)
}

# --- building M -------------------------------------------------------------------

m_build <- function(prep, lk, method = "delaunay", link_mult = 1, sim = NULL,
                    buffer_km = NULL, adaptive = list(mult = 1, weight = FALSE,
                                                      min_km = 1),
                    hull_ratio = 1, corridors = FALSE,
                    closing_km = NULL, close_km = 25, fill_holes = TRUE,
                    add = NULL, add_sites = NULL, seal_km = 5, land = NULL) {
  method <- match.arg(method, c("delaunay", "buffer", "cluster_hull",
                                "convex_hull", "adaptive", "simulation"))
  crs <- prep$crs
  th <- lk$points
  L <- lk$L * link_mult
  r <- L / 2 * 1000
  pts_of <- function(i) {
    sf::st_sfc(lapply(i, function(j) sf::st_point(c(th$x[j], th$y[j]))),
               crs = crs)
  }
  # union of buffers around many points: each point is buffered on its own
  # and the circles are merged tile by tile (exact; buffering one multipoint
  # makes GEOS intersect all circles at once and can need many GB of memory)
  buffer_points <- function(i, rad) {
    if (!length(i)) return(sf::st_sfc(crs = crs))
    tile <- 5 * rad
    key <- paste(floor(th$x[i] / tile), floor(th$y[i] / tile))
    parts <- lapply(split(i, key), function(j) {
      sf::st_union(sf::st_buffer(pts_of(j), rad, nQuadSegs = 16))
    })
    union_geom(do.call(c, unname(parts)))
  }
  all_i <- seq_len(nrow(th))
  close_up <- function(g, km) {
    if (km > 0) buffer_geom(buffer_geom(g, km * 1000), -km * 1000) else g
  }
  base <- switch(method,
    delaunay = {
      t3 <- lk$tri
      parts <- list()
      used <- integer(0)
      if (nrow(t3)) {
        el <- function(a, b) sqrt((th$x[a] - th$x[b])^2 +
                                    (th$y[a] - th$y[b])^2)
        mx <- pmax(el(t3[, 1], t3[, 2]), el(t3[, 2], t3[, 3]),
                   el(t3[, 3], t3[, 1]))
        k <- t3[mx <= L * 1000, , drop = FALSE]
        if (nrow(k)) {
          tri <- sf::st_sfc(lapply(seq_len(nrow(k)), function(i) {
            id <- k[i, c(1, 2, 3, 1)]
            sf::st_polygon(list(cbind(th$x[id], th$y[id])))
          }), crs = crs)
          # the buffer of the triangles also covers r around their corners
          parts[[1]] <- buffer_geom(union_geom(tri), r)
          used <- unique(as.vector(k))
        }
      }
      parts[[length(parts) + 1]] <- buffer_points(setdiff(all_i, used), r)
      close_up(union_geom(do.call(c, parts)), close_km)
    },
    buffer = {
      rb <- if (is.null(buffer_km)) r else buffer_km * 1000
      close_up(buffer_points(all_i, rb), close_km)
    },
    cluster_hull = {
      cl <- if (link_mult == 1) th$cluster else
        link_components(nrow(th), lk$links, L)
      hulls <- lapply(sort(unique(cl)), function(i) {
        idx <- which(cl == i)
        if (hull_ratio >= 1) {
          # hull of the buffers = buffer of the hull of the points (exact)
          buffer_geom(sf::st_convex_hull(union_geom(pts_of(idx))), r)
        } else {
          hull_of(buffer_points(idx, r), hull_ratio)
        }
      })
      union_geom(do.call(c, hulls))
    },
    convex_hull = buffer_geom(sf::st_convex_hull(union_geom(pts_of(all_i))),
                              r),
    simulation = union_geom(sf::st_transform(sim$m, crs)),
    adaptive = {
      cl <- if (link_mult == 1) th$cluster else
        link_components(nrow(th), lk$links, L)
      ids <- sort(unique(cl))
      size <- vapply(ids, function(i) sum(cl == i), numeric(1))
      rad <- vapply(ids, function(i) {
        idx <- which(cl == i)
        if (length(idx) < 2) return(0)
        # mean pairwise distance (a sample of 2000 records for large ones)
        if (length(idx) > 2000) idx <- idx[seq(1, length(idx),
                                               length.out = 2000)]
        mean(stats::dist(cbind(th$x[idx], th$y[idx])))
      }, numeric(1)) * adaptive$mult
      if (isTRUE(adaptive$weight)) rad <- rad * size / max(size)
      rad <- pmax(rad, adaptive$min_km * 1000)
      union_geom(do.call(c, lapply(seq_along(ids), function(k) {
        buffer_points(which(cl == ids[k]), rad[k])
      })))
    })
  base <- union_geom(base)

  corr <- sf::st_sfc(crs = crs)
  L2 <- lk$L2
  if (corridors) {
    ce <- lk$links[lk$links$km > L & lk$links$km <= L2, , drop = FALSE]
    if (nrow(ce)) {
      corr <- sf::st_sfc(lapply(seq_len(nrow(ce)), function(i) {
        sf::st_linestring(rbind(c(th$x[ce$from[i]], th$y[ce$from[i]]),
                                c(th$x[ce$to[i]], th$y[ce$to[i]])))
      }), crs = crs)
      base <- union_geom(c(base, buffer_geom(corr, r)))
    }
  }
  if (is.null(closing_km)) {
    closing_km <- if (corridors) L2 / 2 else if (method == "delaunay") L / 2
    else 0
  }
  auto_m <- close_up(base, closing_km)
  if (fill_holes) auto_m <- drop_holes(auto_m)
  auto_m <- clip_m(auto_m, prep)

  # expert additions
  extra <- list()
  if (!is.null(add)) {
    extra[["add"]] <- union_geom(sf::st_transform(as_sfc(add, "add"), crs))
  }
  if (!is.null(add_sites)) {
    if (!all(c("lon", "lat", "buffer_km") %in% names(add_sites))) {
      stop("`add_sites` needs columns lon, lat and buffer_km.", call. = FALSE)
    }
    s <- sf::st_transform(sf::st_as_sf(add_sites, coords = c("lon", "lat"),
                                       crs = 4326), crs)
    extra[["sites"]] <- union_geom(sf::st_buffer(sf::st_geometry(s),
                                                 add_sites$buffer_km * 1000,
                                                 nQuadSegs = 16))
  }
  final <- auto_m
  if (length(extra)) {
    final <- union_geom(c(auto_m, do.call(c, unname(extra))))
    final <- close_up(final, seal_km)
    if (fill_holes) final <- drop_holes(final)
  }
  final <- clip_m(final, prep)
  if (!is.null(land)) {
    lp <- union_geom(sf::st_transform(as_sfc(land, "land"), crs))
    final <- union_geom(safe_intersection(final, lp))
  }
  added <- if (length(extra)) safe_difference(final, auto_m) else
    sf::st_sfc(crs = crs)

  m <- list(m = final, automatic = auto_m, added = added, corridors = corr,
            points = th, all_points = prep$all, linking = lk,
            settings = list(method = method, L = L, L2 = L2,
                            link_mult = link_mult, buffer_km = buffer_km,
                            hull_ratio = hull_ratio, corridors = corridors,
                            closing_km = closing_km, close_km = close_km,
                            fill_holes = fill_holes, seal_km = seal_km,
                            thin_km = prep$thin_km, group = prep$group,
                            n_input = prep$n_input,
                            n_dropped = prep$n_dropped,
                            expert = names(extra), crs = crs,
                            extent = prep$extent_ll))
  class(m) <- "m_area"
  m$report <- m_report(m)
  m
}

clip_m <- function(g, prep) {
  if (!is.null(prep$extent)) g <- safe_intersection(g, prep$extent)
  if (!is.null(prep$exclude)) g <- safe_difference(g, prep$exclude)
  union_geom(g)
}

hull_of <- function(g, ratio) {
  if (ratio >= 1) return(sf::st_convex_hull(g))
  tryCatch(sf::st_concave_hull(g, ratio = ratio, allow_holes = FALSE),
           error = function(e) {
             warning("Concave hull unavailable (", conditionMessage(e),
                     "); convex hull used.", call. = FALSE)
             sf::st_convex_hull(g)
           })
}

# --- geometry helpers ---------------------------------------------------------------

union_geom <- function(g) {
  if (!length(g)) return(g)
  sf::st_make_valid(sf::st_union(sf::st_make_valid(g)))
}

buffer_geom <- function(g, d) {
  sf::st_make_valid(sf::st_buffer(g, d, nQuadSegs = 16))
}

safe_intersection <- function(a, b) {
  if (!length(a) || !length(b)) return(sf::st_sfc(crs = sf::st_crs(a)))
  out <- suppressWarnings(sf::st_intersection(sf::st_make_valid(a),
                                              sf::st_make_valid(b)))
  polys_only(out)
}

safe_difference <- function(a, b) {
  if (!length(a)) return(a)
  if (!length(b)) return(a)
  out <- suppressWarnings(sf::st_difference(sf::st_make_valid(a),
                                            sf::st_make_valid(b)))
  polys_only(out)
}

# Keep only the polygon parts of a geometry, as single polygons (internal)
as_polygons <- function(g) {
  if (!length(g)) return(g)
  g <- sf::st_geometry(g)
  tp <- as.character(sf::st_geometry_type(g))
  if (any(tp == "GEOMETRYCOLLECTION")) {
    g <- sf::st_collection_extract(g, "POLYGON", warn = FALSE)
    tp <- as.character(sf::st_geometry_type(g))
  }
  g <- g[tp %in% c("POLYGON", "MULTIPOLYGON")]
  if (!length(g)) return(g)
  suppressWarnings(sf::st_cast(g, "POLYGON"))
}

polys_only <- function(g) {
  g <- as_polygons(g)
  if (!length(g)) return(g)
  union_geom(g)
}

drop_holes <- function(g) {
  p <- as_polygons(g)
  if (!length(p)) return(g)
  union_geom(sf::st_sfc(lapply(p, function(x) sf::st_polygon(list(x[[1]]))),
                        crs = sf::st_crs(g)))
}

geom_area <- function(g) {
  if (!length(g)) return(0)
  as.numeric(sum(sf::st_area(g)))
}

# share (%) of all records (per group) inside M
inside_shares <- function(m) {
  ap <- m$all_points
  p <- sf::st_as_sf(ap, coords = c("x", "y"), crs = m$settings$crs)
  inside <- if (length(m$m)) lengths(sf::st_intersects(p, m$m)) > 0 else
    rep(FALSE, nrow(ap))
  grp <- unique(ap$.group)
  vals <- vapply(grp, function(g) round(100 * mean(inside[ap$.group == g]), 1),
                 numeric(1))
  names(vals) <- if (length(grp) == 1 && grp == "all") {
    "records_inside_pct"
  } else paste0("records_inside_pct_", grp)
  as.data.frame(as.list(vals), check.names = FALSE)
}

# Simulated accessible area -------------------------------------------------------

#' Accessible area (M) by simulating dispersal from the records
#'
#' M is the area that has been accessible to the species over a relevant
#' period. Instead of drawing a hull or buffers around the records,
#' `m_simulate()` simulates that history, following the method of
#' Machado-Stredel et al. (2021), in a stable climate or through a glacial
#' cycle (`past`), which is conceptually closer to the definition of M
#' than any geometry:
#'
#' 1. **Niche.** An ellipsoid in the space of the first principal components
#'    of the layers (those explaining `pc_cumvar` of the variance, at most
#'    `max_pc`), with the centroid and covariance of the records. A cell is
#'    suitable when its Mahalanobis distance from the centroid is within the
#'    `niche_level` quantile of the chi-squared distribution (default 0.95,
#'    leaving out the most extreme 5%), with suitability
#'    `exp(-0.5 * MD^2)`.
#' 2. **Dispersal.** In each of `n_rep` replicates, a random half of the
#'    records inside the niche are the starting populations. In each of
#'    `n_events` dispersal events, every colonised cell sends dispersers in
#'    random directions over distances drawn from a normal kernel with
#'    standard deviation `kernel_sd` (in cells): 1 to `max_dispersers`
#'    of them, more from more suitable cells. Every cell reached is
#'    *accessed*; it is *colonised*, and sends dispersers itself, only if
#'    suitable. Cells without values (the sea, cells outside the layers) and
#'    `barrier` cells can be neither accessed nor colonised, so dispersal
#'    can only jump a barrier narrower than the kernel allows.
#' 3. **M.** The mean number of accesses per cell over the replicates; the
#'    accessed cells, without the `access_threshold` share (default 5%)
#'    accessed least often, form M. The cells of all records are always
#'    added, and their number is reported.
#'
#' **Past climate** (their "changing framework"). With `past` (the same
#' layers for a glacial period, e.g. the Last Glacial Maximum from
#' CHELSA-TraCE21k), the simulation runs through a cycle of climatic
#' scenarios of `span_ka` thousand years each: a transition from
#' current-like to glacial conditions, a glacial period, a transition back
#' and a current period (`periods`, in thousands of years; default 5, 10,
#' 5 and 10, as in their paper). Transitional climates are weighted means
#' of the two layers, projected with the PCA of the present, and each
#' scenario lasts `n_events` events. Populations go extinct where a new
#' scenario is unsuitable, so the records' areas may have been reached
#' from refugia. Cells with values only in the past (e.g. shelves exposed
#' at low sea level) are land while the glacial weight dominates and can
#' carry dispersal; M itself is the accessed cells that have values today
#' (the others are counted in the report).
#'
#' The simulation works on the whole grid, so it needs the layers in
#' memory: grids with more than `max_cells` cells are first aggregated on
#' disk (as in [grain_sensitivity()]). The kernel is in cells of the grid
#' used, so choose `kernel_sd` and `n_events` together with the resolution:
#' a disperser moves about `kernel_sd` cells per event, and M grows with
#' both.
#'
#' The result can be compared with the other delimitations in
#' [m_scenarios()] (`method = "simulation"` in [delimit_m()], with `env`).
#'
#' @param points Records: a data.frame with coordinates, a CSV path, `sf`
#'   or [terra::SpatVector].
#' @param env Predictor layers covering a generous study region (larger
#'   than any plausible M), [terra::SpatRaster] or files.
#' @param barrier Optional barriers to dispersal: polygons (`sf`,
#'   SpatVector, path) or a raster whose cells other than `NA` and 0 are
#'   barriers.
#' @param n_rep,n_events Replicates and dispersal events per replicate.
#' @param kernel_sd Standard deviation of the dispersal kernel, in cells.
#' @param max_dispersers Largest number of dispersers per colonised cell
#'   and event.
#' @param niche_level Chi-squared quantile that bounds the niche ellipsoid.
#' @param access_threshold Share of the least accessed cells left out of M.
#' @param past Optional layers of a glacial climate with the names of
#'   `env` (on its grid, or resampled to it), for a glacial cycle.
#' @param periods Durations (thousands of years) of the transition to
#'   glacial conditions, the glacial period, the transition back and the
#'   current period.
#' @param span_ka Duration of one climatic scenario (thousands of years).
#' @param pc_cumvar,max_pc Principal components used for the niche.
#' @param max_cells Largest grid simulated; larger grids are aggregated.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Random seed.
#' @return An object of class `m_simulation`: `m` (M as polygons, `sfc` in
#'   the CRS of `env`), `access` (raster of the mean number of accesses),
#'   `suitability` (the niche suitability today), `m_raster` (1 in M),
#'   `report` (records, starting records inside the niche, records added,
#'   cells and area, cells accessed that have no values today, replicates
#'   in which the species went extinct) and `settings`.
#' @references Machado-Stredel, F., Cobos, M.E., Peterson, A.T. (2021) A
#'   simulation-based method for selecting calibration areas for ecological
#'   niche models and species distribution models. *Frontiers of
#'   Biogeography* 13, e48814. \doi{10.21425/F5FBG48814}
#'
#'   Barve, N., Barve, V., Jimenez-Valverde, A., Lira-Noriega, A., Maher,
#'   S.P., Peterson, A.T., Soberon, J., Villalobos, F. (2011) The crucial
#'   role of the accessible area in ecological niche modeling and species
#'   distribution modeling. *Ecological Modelling* 222, 1810-1819.
#'   \doi{10.1016/j.ecolmodel.2011.02.011}
#' @seealso [delimit_m()], [m_scenarios()], [m_minimal()]
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_enm_data(seed = 1, n_presence = 60)
#' ms <- m_simulate(sim$presence, sim$env, n_rep = 5, n_events = 10,
#'                  kernel_sd = 2, seed = 1)
#' ms
#' plot(ms)
#' }
m_simulate <- function(points, env, barrier = NULL, n_rep = 10,
                       n_events = 20, kernel_sd = 1, max_dispersers = 4,
                       niche_level = 0.95, access_threshold = 0.05,
                       past = NULL,
                       periods = c(to_glacial = 5, glacial = 10,
                                   to_current = 5, current = 10),
                       span_ka = 1, pc_cumvar = 0.95, max_pc = 3,
                       max_cells = 2e6, coords = NULL,
                       points_crs = "EPSG:4326", seed = 1) {
  r <- as_grid(env)
  fact <- 1
  if (terra::ncell(r) > max_cells) {
    fact <- ceiling(sqrt(terra::ncell(r) / max_cells))
    tmp <- tempfile("preENM_msim")
    r <- aggregate_on_disk(r, fact, "average", tmp)
    message("Grid aggregated by ", fact, " to ", terra::ncell(r), " cells.")
  }
  nr <- terra::nrow(r)
  nc <- terra::ncol(r)
  V <- terra::values(r)
  valid <- stats::complete.cases(V)
  Pv <- NULL
  if (!is.null(past)) {
    P <- as_grid(past)
    miss <- setdiff(names(r), names(P))
    if (length(miss)) {
      stop("`past` lacks layer(s) of `env`: ", paste(miss, collapse = ", "),
           call. = FALSE)
    }
    P <- P[[names(r)]]
    if (fact > 1) P <- aggregate_on_disk(P, fact, "average",
                                         tempfile("preENM_msim_past"))
    if (!terra::compareGeom(P, r, stopOnError = FALSE)) {
      P <- terra::resample(P, r, method = "bilinear")
    }
    Pv <- terra::values(P)
  }
  bar <- rep(FALSE, nr * nc)
  if (!is.null(barrier)) {
    b <- if (inherits(barrier, "SpatRaster") || (is.character(barrier) &&
           grepl("\\.(tif|tiff|grd|asc|nc)$", barrier[1],
                 ignore.case = TRUE))) {
      br <- as_grid(barrier)[[1]]
      if (!terra::compareGeom(br, r, stopOnError = FALSE)) {
        br <- terra::resample(br, r, method = "near")
      }
      bv <- terra::values(br, mat = FALSE)
      !is.na(bv) & bv != 0
    } else {
      bp <- terra::vect(sf::st_transform(as_sfc(barrier, "barrier"),
                                         terra::crs(r)))
      bv <- terra::values(terra::rasterize(bp, r, touches = TRUE),
                          mat = FALSE)
      !is.na(bv)
    }
    bar <- b
    valid <- valid & !b
  }
  # records
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  p <- if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy <- if (is.null(coords)) guess_coords(names(df)) else coords
    x <- as_number(df[[xy[1]]])
    y <- as_number(df[[xy[2]]])
    ok <- !is.na(x) & !is.na(y)
    terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else as_points(points, coords, points_crs)
  pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
  rec_cell <- terra::cellFromXY(r, terra::crds(pr))
  rec_cell <- rec_cell[!is.na(rec_cell)]
  rec_ok <- unique(rec_cell[valid[rec_cell]])
  if (length(rec_ok) < 5) {
    stop("Fewer than 5 records in cells with values (outside barriers).",
         call. = FALSE)
  }
  # niche: ellipsoid in principal-component space
  vc <- which(valid)
  smp <- if (length(vc) > 1e5) with_seed(seed, sample(vc, 1e5)) else vc
  pc <- stats::prcomp(V[smp, , drop = FALSE], scale. = TRUE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  k <- min(max_pc, ncol(V), which(cumsum(ve) >= pc_cumvar - 1e-9)[1])
  scores <- function(i) stats::predict(pc, V[i, , drop = FALSE])[, seq_len(k),
                                                                 drop = FALSE]
  Sr <- scores(rec_ok)
  ctr <- colMeans(Sr)
  cv <- stats::cov(Sr)
  if (k == 1) cv <- matrix(cv, 1, 1)
  icv <- tryCatch(solve(cv), error = function(e) MASS_ginv(cv))
  limit <- stats::qchisq(niche_level, df = k)
  # suitability of every valid cell of a climate X (rows = cells); the
  # past is projected with the current PCA, as in Machado-Stredel et al.
  suit_of <- function(X, ok) {
    out <- numeric(nr * nc)
    idx <- which(ok)
    for (s0 in seq(1, length(idx), by = 2e5)) {
      i <- idx[s0:min(length(idx), s0 + 2e5 - 1)]
      d <- sweep(stats::predict(pc, X[i, , drop = FALSE])[, seq_len(k),
                                                          drop = FALSE],
                 2, ctr)
      md2 <- rowSums((d %*% icv) * d)
      out[i] <- ifelse(md2 <= limit, exp(-0.5 * md2), 0)
    }
    out
  }
  S <- suit_of(V, valid)
  start <- rec_ok[S[rec_ok] > 0]
  if (!length(start)) stop("No record lies inside the niche ellipsoid.",
                           call. = FALSE)

  # climatic scenarios: the present only, or a glacial cycle interpolated
  # between the present and `past` (weight of the present: 1 in current
  # conditions, 0 in glacial ones), each kept for `n_events` events; a cell
  # with values only in the past (e.g. an exposed shelf) is land when the
  # weight of the past dominates, and vice versa
  scen <- if (is.null(Pv)) {
    list(list(w = 1, S = S, valid = valid))
  } else {
    need <- c("to_glacial", "glacial", "to_current", "current")
    if (!all(need %in% names(periods))) {
      stop("`periods` must be named ", paste(need, collapse = ", "), ".",
           call. = FALSE)
    }
    nper <- round(periods[need] / span_ka)
    nper[nper < 0] <- 0
    wd <- c(rev(seq_len(nper[["to_glacial"]])) / (nper[["to_glacial"]] + 1),
            rep(0, nper[["glacial"]]),
            seq_len(nper[["to_current"]]) / (nper[["to_current"]] + 1),
            rep(1, nper[["current"]]))
    okc <- stats::complete.cases(V)
    okp <- stats::complete.cases(Pv)
    lapply(wd, function(w) {
      X <- V
      both <- okc & okp
      X[both, ] <- w * V[both, ] + (1 - w) * Pv[both, ]
      only_p <- okp & !okc
      X[only_p, ] <- Pv[only_p, ]
      ok <- (both | (only_p & w < 0.5) | (okc & !okp & w >= 0.5)) & !bar
      St <- suit_of(X, ok)
      list(w = w, idx = which(St > 0), val = St[St > 0], valid = ok)
    })
  }
  s_of <- function(sc) {
    if (!is.null(sc$S)) return(sc$S)
    out <- numeric(nr * nc)
    out[sc$idx] <- sc$val
    out
  }

  # dispersal
  move <- function(cells, n, valid) {
    rc <- cbind((cells - 1) %/% nc + 1, (cells - 1) %% nc + 1)
    d <- stats::rnorm(n, 0, kernel_sd)
    a <- stats::runif(n, 0, 2 * pi)
    rr <- rc[, 1] + round(d * sin(a))
    cc <- rc[, 2] + round(d * cos(a))
    inside <- rr >= 1 & rr <= nr & cc >= 1 & cc <= nc
    out <- rep(NA_integer_, n)
    out[inside] <- as.integer((rr[inside] - 1) * nc + cc[inside])
    out[!is.na(out)][!valid[out[!is.na(out)]]] <- NA_integer_
    out
  }
  Dsum <- numeric(nr * nc)
  extinct <- 0L
  with_seed(seed, {
    for (rep_i in seq_len(n_rep)) {
      st <- start[sample.int(length(start), max(1, ceiling(length(start) /
                                                             2)))]
      col <- logical(nr * nc)
      col[st] <- TRUE
      D <- numeric(nr * nc)
      D[st] <- 1
      for (sc in scen) {
        St <- s_of(sc)
        col[col & St <= 0] <- FALSE        # extinction where unsuitable
        for (ev in seq_len(n_events)) {
          cl <- which(col)
          if (!length(cl)) break
          nd <- pmin(max_dispersers, pmax(1, ceiling(St[cl] *
                                                       max_dispersers)))
          from <- rep(cl, nd)
          to <- move(from, length(from), sc$valid)
          to <- to[!is.na(to)]
          if (!length(to)) next
          D <- D + tabulate(to, nbins = nr * nc)
          col[to[St[to] > 0]] <- TRUE
        }
      }
      if (!any(col)) extinct <- extinct + 1L
      Dsum <- Dsum + D
    }
  })
  Dm <- Dsum / n_rep
  # M is a calibration area for today: cells with present-day values
  not_today <- sum(Dm > 0 & !valid)
  acc <- Dm > 0 & valid
  thr <- stats::quantile(Dm[acc], access_threshold, names = FALSE)
  inM <- acc & Dm >= thr
  outside <- setdiff(unique(rec_cell), which(inM))
  inM[unique(rec_cell)] <- TRUE

  mk <- function(v, nm) {
    x <- terra::rast(r[[1]])
    terra::values(x) <- v
    names(x) <- nm
    x
  }
  access <- mk(ifelse(valid, Dm, NA), "access")
  m_r <- mk(ifelse(inM, 1, NA), "M")
  suit <- mk(ifelse(valid, S, NA), "suitability")
  poly <- sf::st_as_sfc(sf::st_as_sf(terra::as.polygons(m_r, dissolve = TRUE)))
  poly <- union_geom(poly)
  area <- terra::global(terra::mask(terra::cellSize(m_r, unit = "km"), m_r),
                        "sum", na.rm = TRUE)[[1]]
  structure(list(
    m = poly, access = access, suitability = suit, m_raster = m_r,
    report = data.frame(records = length(unique(rec_cell)),
                        in_niche = length(start),
                        added_to_m = length(outside),
                        cells_accessed = sum(acc), cells_in_m = sum(inM),
                        area_km2 = round(area),
                        accessed_not_today = not_today,
                        extinct_replicates = extinct),
    settings = list(n_rep = n_rep, n_events = n_events,
                    n_scenarios = length(scen), past = !is.null(Pv),
                    periods = if (is.null(Pv)) NULL else periods,
                    span_ka = span_ka,
                    kernel_sd = kernel_sd, max_dispersers = max_dispersers,
                    niche_level = niche_level,
                    access_threshold = access_threshold, n_pc = k,
                    variance = sum(ve[seq_len(k)]), aggregated = fact,
                    barrier = !is.null(barrier))),
    class = "m_simulation")
}

# Generalised inverse for a singular covariance (internal)
MASS_ginv <- function(X) {
  s <- svd(X)
  pos <- s$d > max(sqrt(.Machine$double.eps) * s$d[1], 0)
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

#' @export
print.m_simulation <- function(x, ...) {
  s <- x$settings
  rp <- x$report
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  cat("Simulated accessible area (M): ", s$n_rep, " replicates x ",
      if (isTRUE(s$past)) paste0(s$n_scenarios, " climatic scenarios x "),
      s$n_events, " dispersal events, kernel SD ", s$kernel_sd,
      " cell(s)\n\nIN SHORT\n", sep = "")
  say(paste0("M covers ", formatC(rp$area_km2, format = "d", big.mark = " "),
             " km2 (", rp$cells_in_m, " cells): the cells reached by ",
             "dispersal from the records, without the ",
             s$access_threshold * 100, "% least accessed."))
  say(paste0(rp$in_niche, " of the ", rp$records, " record cells lie inside ",
             "the niche ellipsoid and served as starting populations",
             if (rp$records > rp$in_niche) paste0("; the other ",
               rp$records - rp$in_niche, " are outside it (extreme ",
               "conditions) and could not colonise") else "", "."))
  if (isTRUE(s$past)) {
    say(paste0("Simulated through a glacial cycle of ",
               sum(s$periods), " thousand years (", s$n_scenarios,
               " scenarios of ", s$span_ka, " ka).",
               if (rp$extinct_replicates > 0) paste0(
                 " The species went extinct in ", rp$extinct_replicates,
                 " of ", s$n_rep, " replicates: the niche ellipsoid finds ",
                 "no glacial refugium from those starting points.") else "",
               if (rp$accessed_not_today > 0) paste0(
                 " ", rp$accessed_not_today, " accessed cells have no ",
                 "values today (e.g. shelves exposed at low sea level) and ",
                 "are not part of M.") else ""))
  }
  if (rp$added_to_m > 0) {
    say(paste0(rp$added_to_m, " record cell(s) were not reached by the ",
               "simulation and were added to M: consider more events, a ",
               "wider kernel, or whether these records are reliable."))
  }
  say(paste0("M grows with the kernel and the number of events: compare ",
             "a few settings, and other delimitations, with m_scenarios()."))
  cat("", strwrap(paste0(
    "Niche: ellipsoid on PC1-PC", s$n_pc, " (", round(100 * s$variance),
    "% of the variance) bounded at the ", s$niche_level, " chi-squared ",
    "quantile; up to ", s$max_dispersers, " dispersers per cell and event",
    if (s$aggregated > 1) paste0("; grid aggregated by ", s$aggregated)
    else "", if (s$barrier) "; with barriers" else "",
    ". Method of Machado-Stredel et al. (2021), ",
    if (isTRUE(s$past)) "glacial cycle." else "stable climate."),
    width = 76), sep = "\n")
  invisible(x)
}

#' Plot a simulated accessible area
#'
#' The mean number of accesses per cell over the replicates, the outline of
#' M and the records' cells.
#'
#' @param x An `m_simulation` from [m_simulate()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @param ... Not used.
#' @return A ggplot object.
#' @references Machado-Stredel, F., Cobos, M.E., Peterson, A.T. (2021) A
#'   simulation-based method for selecting calibration areas for ecological
#'   niche models and species distribution models. *Frontiers of
#'   Biogeography* 13, e48814. \doi{10.21425/F5FBG48814}
#' @export
#' @examples
#' # see ?m_simulate
plot.m_simulation <- function(x, base_size = 11, file = NULL, width = 7,
                              height = 6, ...) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  a <- x$access
  if (terra::ncell(a) > 5e5) {
    a <- terra::spatSample(a, 5e5, method = "regular", as.raster = TRUE)
  }
  d <- terra::as.data.frame(a, xy = TRUE, na.rm = TRUE)
  d <- d[d$access > 0, , drop = FALSE]
  lonlat <- terra::is.lonlat(a)
  ggplot2::ggplot() +
    ggplot2::geom_raster(data = d, ggplot2::aes(.data$x, .data$y,
                                                fill = .data$access)) +
    ggplot2::scale_fill_distiller(palette = "YlGnBu", direction = 1,
                                  trans = "log10",
                                  name = "Mean accesses\nper cell") +
    ggplot2::geom_sf(data = sf::st_sf(geometry = x$m), fill = NA,
                     colour = "#d95f02", linewidth = 0.6,
                     inherit.aes = FALSE) +
    ggplot2::labs(x = NULL, y = NULL, title = "Simulated accessible area",
                  subtitle = paste0(x$settings$n_rep, " replicates x ",
                                    x$settings$n_events, " events; ",
                                    "outline: M (",
                                    formatC(x$report$area_km2, format = "d",
                                            big.mark = " "), " km2)")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
}

# Minimal background area ---------------------------------------------------------

#' Smallest calibration area whose environment has stabilised
#'
#' How large must the calibration area be? Rotllan-Puig & Traveset (2021)
#' grew the background in steps from the centre of the records and found
#' the smallest area beyond which models stopped improving (their MinBAR
#' method). preENM fits no model, so `m_minimal()` keeps the geometry of
#' their procedure and asks a model-free question instead: from which
#' extent on does the *environment available* stop changing as the area
#' grows?
#'
#' * **Extents.** Circles around the centre of the records (their median
#'   coordinates, or the mean), whose radii contain growing shares of the
#'   records (`steps` equal steps of records, as in MinBAR, so that
#'   scattered records do not leave large empty steps), then larger
#'   circles beyond the farthest record (`beyond`, multiples of its
#'   distance).
#' * **Environment of each extent.** `n` random points within the circle,
#'   with values in every layer.
#' * **Change.** The largest Kolmogorov-Smirnov distance, over the layers,
#'   between the environment of an extent and of the previous one
#'   (`ks_step`), and between each extent and the largest (`ks_largest`).
#'   The environment has stabilised when `ks_step` is below the 5%
#'   critical value of the two-sample test (`1.36 * sqrt(2 / n)`, about
#'   0.019 for `n = 10000`; consecutive extents cannot be told apart) for
#'   `k_stable` consecutive steps.
#'
#' The stable extent is a lower bound for a calibration area that
#' represents the conditions around the records; it is not a delimitation
#' of M by dispersal (see [delimit_m()] and [m_simulate()]). The extent of
#' the calibration area strongly affects model accuracy (Luna et al. 2024),
#' so compare it with the other delimitations in [m_scenarios()].
#'
#' @param points Records: a data.frame with coordinates, a CSV path, `sf`
#'   or [terra::SpatVector].
#' @param env Predictor layers covering a generous region around the
#'   records.
#' @param steps Number of extents up to the farthest record.
#' @param beyond Multiples of the farthest record's distance added as
#'   further extents.
#' @param centre `"median"` or `"mean"` of the record coordinates.
#' @param n Random points per extent.
#' @param k_stable Consecutive stable steps required.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Random seed.
#' @return An object of class `m_minimal`: `extents` (one row per extent:
#'   `radius_km`, `area_km2`, `records_pct`, `n_points`, `ks_step`,
#'   `ks_largest`, `stable`), `minimal` (the row of the first stable
#'   extent, or `NULL`), `m` (that circle as an `sfc` polygon in WGS84),
#'   `centre` and `settings`.
#' @references Rotllan-Puig, X., Traveset, A. (2021) Determining the minimal
#'   background area for species distribution models: MinBAR package.
#'   *Ecological Modelling* 439, 109353.
#'   \doi{10.1016/j.ecolmodel.2020.109353}
#'
#'   Luna, S., Pena-Peniche, A., Mendoza-Alfaro, R. (2024) Species
#'   distribution model accuracy is strongly influenced by the choice of
#'   calibration area. *Biodiversity Informatics* 18, 43-55.
#'   \doi{10.17161/bi.v18i.22655}
#' @seealso [delimit_m()], [m_simulate()], [m_scenarios()]
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_enm_data(seed = 1, n_presence = 60)
#' mm <- m_minimal(sim$presence, sim$env, steps = 6, n = 3000, seed = 1)
#' mm
#' plot(mm)
#' }
m_minimal <- function(points, env, steps = 10, beyond = c(1.25, 1.5, 2),
                      centre = c("median", "mean"), n = 10000, k_stable = 2,
                      coords = NULL, points_crs = "EPSG:4326", seed = 1) {
  centre <- match.arg(centre)
  r <- as_grid(env)
  if (is.character(points) && length(points) == 1) {
    points <- read_any_csv(points)
  }
  p <- if (is.data.frame(points) && !inherits(points, "sf")) {
    df <- as.data.frame(points, stringsAsFactors = FALSE)
    xy <- if (is.null(coords)) guess_coords(names(df)) else coords
    x <- as_number(df[[xy[1]]])
    y <- as_number(df[[xy[2]]])
    ok <- !is.na(x) & !is.na(y)
    terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else as_points(points, coords, points_crs)
  ll <- terra::crds(terra::project(p, "EPSG:4326"))
  if (nrow(ll) < 3) stop("At least 3 records are needed.", call. = FALSE)
  ctr <- if (centre == "median") apply(ll, 2, stats::median) else
    colMeans(ll)
  crs <- sprintf(paste0("+proj=laea +lat_0=%s +lon_0=%s +x_0=0 +y_0=0 ",
                        "+datum=WGS84 +units=m +no_defs"),
                 num_arg(ctr[2]), num_arg(ctr[1]))
  pe <- terra::crds(terra::project(terra::vect(ll, crs = "EPSG:4326"), crs))
  dist <- sqrt(rowSums(pe^2)) / 1000
  far <- max(dist)
  radii <- c(stats::quantile(dist, seq_len(steps) / steps, names = FALSE),
             far * beyond)
  radii <- sort(unique(pmax(radii, 1e-3)))
  # environment of each extent: random points within the circle
  draw <- function(R, k) {
    with_seed(seed + k, {
      rad <- R * 1000 * sqrt(stats::runif(2 * n))
      a <- stats::runif(2 * n, 0, 2 * pi)
      xy <- cbind(rad * cos(a), rad * sin(a))
      pts <- terra::project(terra::vect(xy, crs = crs), terra::crs(r))
      v <- as.matrix(terra::extract(r, pts, ID = FALSE))
      v <- v[stats::complete.cases(v), , drop = FALSE]
      v[seq_len(min(n, nrow(v))), , drop = FALSE]
    })
  }
  envs <- lapply(seq_along(radii), function(k) draw(radii[k], k))
  vars <- names(r)
  ks <- function(A, B) {
    if (!nrow(A) || !nrow(B)) return(NA_real_)
    max(vapply(seq_along(vars), function(j) {
      u <- sort(unique(c(A[, j], B[, j])))
      max(abs(stats::ecdf(A[, j])(u) - stats::ecdf(B[, j])(u)))
    }, numeric(1)))
  }
  m <- length(radii)
  ks_step <- c(NA, vapply(2:m, function(k) ks(envs[[k]], envs[[k - 1]]),
                          numeric(1)))
  ks_largest <- vapply(seq_len(m), function(k) ks(envs[[k]], envs[[m]]),
                       numeric(1))
  npts <- vapply(envs, nrow, integer(1))
  crit <- 1.36 * sqrt((npts + c(NA, npts[-m])) / (npts * c(NA, npts[-m])))
  below <- !is.na(ks_step) & ks_step <= crit
  # stable from step k when this and the next k_stable - 1 steps are below
  stable <- vapply(seq_len(m), function(k) {
    j <- k:min(m, k + k_stable - 1)
    length(j) == k_stable && all(below[j])
  }, logical(1))
  tab <- data.frame(step = seq_len(m), radius_km = radii,
                    area_km2 = pi * radii^2,
                    records_pct = vapply(radii, function(R)
                      100 * mean(dist <= R + 1e-9), numeric(1)),
                    n_points = npts, ks_step = ks_step, critical = crit,
                    ks_largest = ks_largest, stable = stable)
  first <- which(stable)[1]
  # the extent before the first indistinguishable step is the minimal one
  minimal <- if (is.na(first)) NULL else tab[max(1, first - 1), ]
  circle <- if (is.null(minimal)) NULL else {
    g <- sf::st_buffer(sf::st_sfc(sf::st_point(c(0, 0)), crs = crs),
                       minimal$radius_km * 1000, nQuadSegs = 32)
    sf::st_transform(g, 4326)
  }
  structure(list(extents = tab, minimal = minimal, m = circle,
                 centre = ctr,
                 settings = list(steps = steps, beyond = beyond, n = n,
                                 k_stable = k_stable, centre = centre,
                                 farthest_km = far)),
            class = "m_minimal")
}

#' @export
print.m_minimal <- function(x, digits = 3, ...) {
  t <- x$extents
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  f <- function(v) formatC(round(v), format = "d", big.mark = " ")
  cat("Minimal background area: ", nrow(t), " extents around the ", s$centre,
      " centre of the records, up to ", f(max(t$radius_km)),
      " km\n\nIN SHORT\n", sep = "")
  if (is.null(x$minimal)) {
    say(paste0("The environment did not stabilise within the extents ",
               "tested: each larger circle still adds conditions that ",
               "differ from the smaller one. Test larger extents (`beyond`), ",
               "or accept that the environment around the records keeps ",
               "changing (e.g. a continental gradient), and choose M on ",
               "other grounds (dispersal, m_simulate())."))
  } else {
    mm <- x$minimal
    say(paste0("The available environment stabilises from a radius of about ",
               f(mm$radius_km), " km (", f(mm$area_km2), " km2, ",
               round(mm$records_pct), "% of the records inside): larger ",
               "circles add conditions that cannot be told apart from it."))
    if (mm$records_pct < 100) {
      say(paste0("That circle leaves ", round(100 - mm$records_pct),
                 "% of the records outside: a calibration area should still ",
                 "contain all records, so use it as a lower bound."))
    }
  }
  say(paste0("This concerns the environment only, not model performance ",
             "(MinBAR used models). The extent of the calibration area ",
             "strongly affects model accuracy (Luna et al. 2024): compare ",
             "it with other delimitations in m_scenarios()."))
  cat("\n")
  show <- t[c("step", "radius_km", "area_km2", "records_pct", "n_points",
              "ks_step", "critical", "ks_largest", "stable")]
  show$radius_km <- round(show$radius_km)
  show$area_km2 <- round(show$area_km2)
  show$records_pct <- round(show$records_pct, 1)
  for (k in c("ks_step", "critical", "ks_largest")) {
    show[[k]] <- round(show[[k]], digits)
  }
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "ks_step = largest Kolmogorov-Smirnov distance, over the layers, ",
    "between the environment of an extent and the previous one; critical ",
    "= its 5% critical value; ks_largest = the same against the largest ",
    "extent. Stable = ks_step below the critical value for ", s$k_stable,
    " consecutive steps; the minimal extent is the one before the first ",
    "stable step."), width = 76), sep = "\n")
  invisible(x)
}

#' Plot the stabilisation of the environment with the extent
#'
#' The change of the environment between consecutive extents and against
#' the largest extent, with the 5% critical value and the minimal extent,
#' from [m_minimal()].
#'
#' @param x An `m_minimal` object.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @param ... Not used.
#' @return A ggplot object.
#' @references Rotllan-Puig, X., Traveset, A. (2021) Determining the minimal
#'   background area for species distribution models: MinBAR package.
#'   *Ecological Modelling* 439, 109353.
#'   \doi{10.1016/j.ecolmodel.2020.109353}
#' @export
#' @examples
#' # see ?m_minimal
plot.m_minimal <- function(x, base_size = 11, file = NULL, width = 7,
                           height = 4.5, ...) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  t <- x$extents
  lv <- c("Change from the previous extent", "Difference from the largest")
  d <- rbind(data.frame(radius = t$radius_km, value = t$ks_step,
                        series = lv[1]),
             data.frame(radius = t$radius_km, value = t$ks_largest,
                        series = lv[2]))
  d <- d[!is.na(d$value), , drop = FALSE]
  d$series <- factor(d$series, lv)
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$radius, .data$value,
                                       colour = .data$series)) +
    ggplot2::geom_line(data = t[!is.na(t$critical), , drop = FALSE],
                       ggplot2::aes(.data$radius_km, .data$critical),
                       inherit.aes = FALSE, linetype = "dashed",
                       colour = "#8a8a86") +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c("#1b9e77", "#7570b3"), lv), name = NULL) +
    ggplot2::labs(x = "Radius around the centre of the records (km)",
                  y = "Largest Kolmogorov-Smirnov distance",
                  title = "Environment of growing extents",
                  subtitle = paste0(
                    "Dashed: 5% critical value",
                    if (!is.null(x$minimal)) paste0(
                      "; vertical: minimal extent (",
                      round(x$minimal$radius_km), " km)") else "")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom")
  if (!is.null(x$minimal)) {
    p <- p + ggplot2::geom_vline(xintercept = x$minimal$radius_km,
                                 colour = "#d95f02")
  }
  p
}
