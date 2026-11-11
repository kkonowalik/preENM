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
#' @param method `"delaunay"`, `"buffer"`, `"cluster_hull"` or
#'   `"convex_hull"`.
#' @param link_km Linking distance L in km; `NULL` (default) = derived from the
#'   data.
#' @param link_search Range (km) in which the automatic L is searched.
#' @param link2_max_km Upper limit (km) for the second linking distance L2.
#' @param link_mult Multiplier of L for the method (sensitivity analyses, e.g.
#'   0.5 or 2).
#' @param buffer_km Buffer radius for `method = "buffer"` (default L/2).
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
                                         "convex_hull"),
                      link_km = NULL, link_search = c(100, 600),
                      link2_max_km = 1000, link_mult = 1, buffer_km = NULL,
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
  m_build(prep, lk, method = method, link_mult = link_mult,
          buffer_km = buffer_km, hull_ratio = hull_ratio,
          corridors = corridors, closing_km = closing_km,
          close_km = close_km, fill_holes = fill_holes, add = add,
          add_sites = add_sites, seal_km = seal_km, land = land)
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
#' scenario inside the reference, and the Jaccard index). Use it to show how
#' sensitive M is to the method, as a supplement to the chosen M.
#'
#' @inheritParams delimit_m
#' @param scenarios Named list; each element is a list of arguments of
#'   [delimit_m()] (method, link_mult, corridors, buffer_km, ...). The default
#'   compares the Delaunay M at 0.5, 1 and 2 times L, with corridors, the
#'   cluster and overall convex hulls, and fixed 100 km buffers.
#' @param reference Name of the reference scenario (default: the first).
#' @param ... Arguments passed to [delimit_m()] for every scenario (e.g.
#'   `add`, `exclude`, `extent`, `link_km`).
#' @return An object of class `m_scenarios`: a list with `areas` (the
#'   `m_area` objects) and `comparison` (a data.frame).
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

m_build <- function(prep, lk, method = "delaunay", link_mult = 1,
                    buffer_km = NULL, hull_ratio = 1, corridors = FALSE,
                    closing_km = NULL, close_km = 25, fill_holes = TRUE,
                    add = NULL, add_sites = NULL, seal_km = 5, land = NULL) {
  method <- match.arg(method, c("delaunay", "buffer", "cluster_hull",
                                "convex_hull"))
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
                              r))
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
