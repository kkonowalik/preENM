# Reports, previews and files for the accessible area (M) -----------------------

#' @export
print.m_area <- function(x, ...) {
  s <- x$settings
  cat("Accessible area (M): method ", s$method,
      if (s$link_mult != 1) paste0(" at ", s$link_mult, " x L"), "\n",
      sep = "")
  cat("Linking distance L = ", round(x$linking$L), " km",
      if (x$linking$auto) " (from the data)" else " (set by hand)",
      "; L2 = ", round(s$L2), " km; ", x$linking$n_clusters,
      " cluster(s) at L\n", sep = "")
  cat("Records: ", s$n_input, " read, ", s$n_dropped,
      " outside the extent or excluded, ", nrow(x$points),
      " after thinning to ", s$thin_km, " km\n", sep = "")
  cat("Closing ", round(s$closing_km), " km; holes filled: ", s$fill_holes,
      "; corridors: ", s$corridors,
      if (length(s$expert)) paste0("; expert additions: ",
                                  paste(s$expert, collapse = ", ")),
      "\n\n", sep = "")
  print(x$report, row.names = FALSE)
  invisible(x)
}

#' @export
print.m_scenarios <- function(x, ...) {
  cat("M scenarios (reference: ", x$reference, "); L = ", round(x$linking$L),
      " km, L2 = ", round(x$linking$L2), " km\n\n", sep = "")
  print(x$comparison, row.names = FALSE)
  invisible(x)
}

#' Report of an accessible area (M)
#'
#' Areas of the automatic M, of what the expert additions contribute (area
#' outside the automatic M) and of the final M, with the share of records
#' inside. Areas include the sea unless `land` was used; give a `raster`
#' (any predictor, which is `NA` over the sea) to also get the land area.
#'
#' @param m An `m_area` object from [delimit_m()].
#' @param raster Optional raster; adds the area of its non-`NA` cells.
#' @return A data.frame with one row per part.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' m_report(m, raster = sim$env)
m_report <- function(m, raster = NULL) {
  check_m(m)
  parts <- list(automatic = m$automatic, added_by_expert = m$added,
                M = m$m)
  out <- data.frame(part = names(parts),
                    area_km2 = vapply(parts, function(g) round(geom_area(g) /
                                                                 1e6),
                                      numeric(1)),
                    stringsAsFactors = FALSE)
  if (!is.null(raster)) {
    r <- as_grid(raster)[[1]]
    out$land_area_km2 <- vapply(parts, function(g) {
      if (!length(g)) return(0)
      v <- terra::vect(sf::st_transform(g, terra::crs(r)))
      mk <- terra::mask(terra::crop(r, v, snap = "out"), v)
      cs <- terra::mask(terra::cellSize(mk, unit = "km"), mk)
      round(terra::global(cs, "sum", na.rm = TRUE)[[1]])
    }, numeric(1))
  }
  ins <- inside_shares(m)
  auto_m <- m
  auto_m$m <- m$automatic
  ins_auto <- inside_shares(auto_m)
  for (k in names(ins)) {
    out[[k]] <- c(ins_auto[[k]], NA, ins[[k]])
  }
  rownames(out) <- NULL
  out
}

check_m <- function(m) {
  if (!inherits(m, "m_area")) {
    stop("`m` must be the output of delimit_m().", call. = FALSE)
  }
}

#' Shares of raster classes within M (e.g. Koppen-Geiger climates)
#'
#' Percentage of the land area of the automatic and the final M in each class
#' of a categorical raster, such as the CHELSA Koppen-Geiger map. Useful to
#' describe M and to see what expert additions or alternative scenarios
#' change.
#'
#' @param m An `m_area` object.
#' @param raster Categorical raster (classes as integers).
#' @param labels Optional class labels, indexed by class value (e.g.
#'   [koppen_geiger_labels()]).
#' @param fact Aggregation factor (modal) to speed up large rasters.
#' @return A data.frame with `class`, `label`, `pct_automatic` and `pct_M`.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' # two classes (1 = west, 2 = east) from a continuous layer
#' cls <- terra::classify(sim$env[["env1"]],
#'                        rbind(c(-Inf, 0.5, 1), c(0.5, Inf, 2)))
#' m_class_shares(m, cls, labels = c("west", "east"))
m_class_shares <- function(m, raster, labels = NULL, fact = 1) {
  check_m(m)
  r <- as_grid(raster)[[1]]
  share <- function(g) {
    if (!length(g)) return(data.frame(class = numeric(0), km2 = numeric(0)))
    v <- terra::vect(sf::st_transform(g, terra::crs(r)))
    x <- terra::crop(r, v, snap = "out")
    if (fact > 1) x <- terra::aggregate(x, fact, fun = "modal", na.rm = TRUE)
    x <- terra::mask(x, v)
    z <- terra::zonal(terra::cellSize(x, unit = "km"), x, fun = "sum")
    names(z) <- c("class", "km2")
    z
  }
  a <- share(m$automatic)
  f <- share(m$m)
  out <- merge(data.frame(class = a$class,
                          pct_automatic = round(100 * a$km2 / sum(a$km2), 2)),
               data.frame(class = f$class,
                          pct_M = round(100 * f$km2 / sum(f$km2), 2)),
               all = TRUE)
  out$label <- if (!is.null(labels)) labels[out$class] else
    as.character(out$class)
  out <- out[order(out$class), c("class", "label", "pct_automatic", "pct_M")]
  rownames(out) <- NULL
  out
}

#' Koppen-Geiger class labels
#'
#' Labels of the 31 Koppen-Geiger classes in the order used by the CHELSA
#' `kg0` map (class 1 = "Af", ..., 31 = "EF"), for [m_class_shares()] and
#' [class_polygons()].
#'
#' @return Character vector of 31 labels.
#' @export
#' @examples
#' koppen_geiger_labels()[c(12, 30)]
koppen_geiger_labels <- function() {
  c("Af", "Am", "As", "Aw", "BWk", "BWh", "BSk", "BSh", "Cfa", "Cfb", "Cfc",
    "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Dfa", "Dfb", "Dfc", "Dfd",
    "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "ET", "EF")
}

#' Polygons of selected raster classes
#'
#' Turns the cells of chosen classes of a categorical raster into polygons,
#' e.g. Mediterranean and arid Koppen-Geiger climates to test how much M
#' changes without them (pass the result to `exclude` in [delimit_m()] or
#' [m_scenarios()]).
#'
#' @param raster Categorical raster.
#' @param classes Class values, or labels if `labels` is given.
#' @param labels Optional class labels indexed by class value.
#' @param fact Aggregation factor (modal) to simplify the polygons.
#' @return An `sf` object with one (multi)polygon.
#' @export
#' @examples
#' r <- terra::rast(nrows = 10, ncols = 10, vals = rep(1:2, 50))
#' class_polygons(r, 2)
class_polygons <- function(raster, classes, labels = NULL, fact = 1) {
  r <- as_grid(raster)[[1]]
  if (!is.null(labels) && is.character(classes)) {
    classes <- match(classes, labels)
    if (anyNA(classes)) stop("Unknown class labels.", call. = FALSE)
  }
  if (fact > 1) r <- terra::aggregate(r, fact, fun = "modal", na.rm = TRUE)
  sel <- terra::classify(r, cbind(classes, 1), others = NA)
  p <- terra::as.polygons(sel, dissolve = TRUE)
  sf::st_as_sf(terra::aggregate(p))
}

#' Clip rasters to the accessible area (M)
#'
#' Crops the rasters to M and sets cells outside M to `NA`, ready for
#' [sample_background()] and model fitting.
#'
#' @param raster A [terra::SpatRaster] or file path(s).
#' @param m An `m_area` object, or polygons (`sf`, [terra::SpatVector] or a
#'   path; for a GeoPackage written by [write_m_area()] its `M` layer is
#'   used).
#' @param filename Optional output file.
#' @return The clipped [terra::SpatRaster].
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' env_m <- crop_to_m(sim$env, m)
#' terra::plot(env_m)
crop_to_m <- function(raster, m, filename = "") {
  r <- as_grid(raster)
  if (is.character(m) && length(m) == 1 &&
      "M" %in% sf::st_layers(m)$name) {
    m <- sf::st_read(m, layer = "M", quiet = TRUE)   # from write_m_area()
  }
  g <- if (inherits(m, "m_area")) m$m else as_sfc(m, "m")
  v <- terra::vect(sf::st_transform(g, terra::crs(r)))
  terra::mask(terra::crop(r, v, snap = "out"), v, filename = filename)
}

# --- previews -----------------------------------------------------------------------

#' Preview of an accessible area (M) or of M scenarios
#'
#' Map of M with the records used, the expert additions highlighted and the
#' corridors, on a light land background (Natural Earth countries when the
#' `rnaturalearth` package is installed, or your own `land` polygons). For
#' [m_scenarios()] each scenario is shown in its own panel.
#'
#' @param x An `m_area` or `m_scenarios` object.
#' @param land Optional land polygons for the background.
#' @param points Show the records (thinned).
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' plot_m_area(m)
plot_m_area <- function(x, land = NULL, points = TRUE, base_size = 11) {
  if (inherits(x, "m_scenarios")) {
    areas <- x$areas
  } else {
    check_m(x)
    areas <- list(M = x)
  }
  first <- areas[[1]]
  crs <- first$settings$crs
  polys <- do.call(rbind, lapply(names(areas), function(nm) {
    a <- areas[[nm]]
    rows <- list(sf::st_sf(scenario = nm, part = "automatic M",
                           geometry = a$automatic))
    if (length(a$added)) {
      rows[[2]] <- sf::st_sf(scenario = nm, part = "expert addition",
                             geometry = a$added)
    }
    do.call(rbind, rows)
  }))
  polys$scenario <- factor(polys$scenario, levels = names(areas))
  pt <- sf::st_as_sf(first$points, coords = c("x", "y"), crs = crs)
  bb <- sf::st_bbox(c(sf::st_geometry(polys), sf::st_geometry(pt)))
  pad <- 0.05 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]])
  xl <- c(bb[["xmin"]] - pad, bb[["xmax"]] + pad)
  yl <- c(bb[["ymin"]] - pad, bb[["ymax"]] + pad)
  bg <- basemap(land, xl, yl, crs)

  ink <- "#3a3a38"
  p <- ggplot2::ggplot()
  if (!is.null(bg)) {
    p <- p + ggplot2::geom_sf(data = bg, fill = "#ecebe7", colour = "#c4c3be",
                              linewidth = 0.2)
  }
  p <- p + ggplot2::geom_sf(data = polys,
                            ggplot2::aes(fill = .data$part),
                            colour = NA, alpha = 0.55) +
    ggplot2::scale_fill_manual(values = c("automatic M" = "#3987e5",
                                          "expert addition" = "#e5822e"),
                               name = NULL)
  cor <- do.call(c, lapply(areas, function(a) a$corridors))
  if (length(cor)) {
    p <- p + ggplot2::geom_sf(data = sf::st_sf(geometry = cor),
                              colour = "#104281", linewidth = 0.5,
                              linetype = "dashed")
  }
  if (points) {
    # the most numerous group is drawn in dark ink, the others highlighted
    grp <- names(sort(table(first$points$.group), decreasing = TRUE))
    if (length(grp) > 1) {
      pt$group <- factor(first$points$.group, levels = grp)
      pt <- pt[order(-as.integer(pt$group)), ]        # rare groups on top
      cols <- stats::setNames(rep(c("#1f1f1d", "#c9480c", "#6b3fa0",
                                    "#0f7b5f"), length.out = length(grp)), grp)
      shp <- stats::setNames(rep(c(16, 17, 15, 18), length.out = length(grp)),
                             grp)
      p <- p + ggplot2::geom_sf(data = pt, ggplot2::aes(shape = .data$group,
                                                        colour = .data$group),
                                size = 0.9) +
        ggplot2::scale_colour_manual(values = cols, name = NULL) +
        ggplot2::scale_shape_manual(values = shp, name = NULL)
    } else {
      p <- p + ggplot2::geom_sf(data = pt, colour = "#1f1f1d", size = 0.6)
    }
  }
  s <- first$settings
  p <- p + ggplot2::coord_sf(crs = crs, xlim = xl, ylim = yl, expand = FALSE) +
    ggplot2::labs(
      title = if (length(areas) > 1) "Accessible area (M): scenarios" else
        paste0("Accessible area (M): ", s$method),
      subtitle = paste0("L = ", round(first$linking$L), " km (",
                        if (first$linking$auto) "from the data" else "by hand",
                        "); ", nrow(first$points), " thinned records")) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(text = ggplot2::element_text(colour = ink),
                   plot.subtitle = ggplot2::element_text(colour = "#8a8a86"),
                   panel.grid = ggplot2::element_line(colour = "#e6e5e1",
                                                      linewidth = 0.2),
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = NA),
                   legend.position = "bottom")
  if (length(areas) > 1) {
    p <- p + ggplot2::facet_wrap(~scenario)
  }
  p
}

# Land background cropped to the plotted area (NULL when unavailable)
basemap <- function(land, xl, yl, crs) {
  tryCatch({
    if (is.null(land)) {
      if (!requireNamespace("rnaturalearth", quietly = TRUE)) return(NULL)
      land <- suppressMessages(rnaturalearth::ne_countries(
        scale = "small", returnclass = "sf"))
    }
    g <- as_sfc(land, "land")
    box <- sf::st_as_sfc(sf::st_bbox(c(xmin = xl[1], xmax = xl[2],
                                       ymin = yl[1], ymax = yl[2]),
                                     crs = sf::st_crs(crs)))
    box_ll_poly <- sf::st_transform(sf::st_segmentize(box, diff(xl) / 50),
                                    4326)
    old <- suppressMessages(sf::sf_use_s2(FALSE))
    on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
    g <- sf::st_transform(g, 4326)
    bbx <- sf::st_bbox(box_ll_poly)
    crop1 <- function(x) suppressMessages(suppressWarnings(sf::st_crop(x, bbx)))
    g <- tryCatch(crop1(g), error = function(e) {
      # an invalid polygon (e.g. one crossing the 180th meridian) must not
      # remove the whole background: crop polygon by polygon, skip failures
      pp <- suppressWarnings(sf::st_cast(g, "POLYGON"))
      parts <- lapply(seq_along(pp), function(i) {
        tryCatch(crop1(pp[i]), error = function(e2) NULL)
      })
      do.call(c, Filter(function(z) !is.null(z) && length(z), parts))
    })
    g <- g[!sf::st_is_empty(g)]
    if (!length(g)) return(NULL)
    sf::st_sf(geometry = sf::st_transform(g, crs))
  }, error = function(e) NULL)
}

#' Plot the linking-distance diagnostics
#'
#' `type = "links"` shows the longest links of the minimum spanning tree (the
#' jump at L separates sampling gaps from gaps between separate parts of the
#' range); `type = "clusters"` shows the number of clusters as a function of
#' the linking distance. Dashed lines mark L (and L2).
#'
#' @param x An `m_linking` object from [linking_distance()] or an `m_area`.
#' @param type `"links"` or `"clusters"`.
#' @param n_links Number of longest links shown.
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' lk <- linking_distance(sim$presence, thin_km = 5, link_search = c(20, 200))
#' plot_linking(lk)
#' plot_linking(lk, "clusters")
plot_linking <- function(x, type = c("links", "clusters"), n_links = 40,
                         base_size = 11) {
  type <- match.arg(type)
  if (inherits(x, "m_area")) x <- x$linking
  if (!inherits(x, "m_linking")) {
    stop("`x` must come from linking_distance() or delimit_m().",
         call. = FALSE)
  }
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  marks <- data.frame(what = c("L", "L2"), v = c(x$L, x$L2))
  marks <- marks[!duplicated(round(marks$v, 6)), , drop = FALSE]
  if (type == "links") {
    h <- utils::tail(x$links$km, n_links)
    d <- data.frame(rank = seq_along(h), km = h)
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$rank, .data$km)) +
      ggplot2::geom_hline(data = marks, ggplot2::aes(yintercept = .data$v),
                          linetype = "dashed", colour = "#c9480c",
                          linewidth = 0.4) +
      ggplot2::geom_text(data = marks, ggplot2::aes(x = 1, y = .data$v,
                                                    label = paste0(.data$what, " = ", round(.data$v), " km")),
                         hjust = 0, vjust = -0.5, colour = "#c9480c",
                         size = base_size / 3.4) +
      ggplot2::geom_point(colour = "#1c5cab", size = 2) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(x = paste0("Rank (1 = shortest of the ", length(h),
                               " longest links)"),
                    y = "Link length (km, log scale)",
                    title = "Longest links of the minimum spanning tree")
  } else {
    d <- x$table
    p <- ggplot2::ggplot(d, ggplot2::aes(.data$link_km, .data$n_clusters)) +
      ggplot2::geom_vline(data = marks, ggplot2::aes(xintercept = .data$v),
                          linetype = "dashed", colour = "#c9480c",
                          linewidth = 0.4) +
      ggplot2::geom_line(colour = "#1c5cab", linewidth = 0.6) +
      ggplot2::geom_point(colour = "#1c5cab", size = 2) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::labs(x = "Linking distance (km, log scale)",
                    y = "Number of clusters (log scale)",
                    title = "Single-linkage clusters")
  }
  p + ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(text = ggplot2::element_text(colour = ink),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = NA),
                   axis.title = ggplot2::element_text(colour = muted))
}

# --- files ---------------------------------------------------------------------------

#' Write an accessible area (M) and its reports to files
#'
#' Writes a GeoPackage (`<prefix>.gpkg`, EPSG:4326) with the layers `M`,
#' `M_automatic`, `expert_additions` and `corridors` (when present) and
#' `records` (thinned records with their cluster), plus `<prefix>_report.csv`,
#' `<prefix>_connectivity.csv` (clusters by linking distance),
#' `<prefix>_links.csv` (spanning-tree links), and the previews
#' `<prefix>_map.png` and `<prefix>_links.png`. For [m_scenarios()], each
#' scenario is a layer and the comparison table is written.
#'
#' @param x An `m_area` or `m_scenarios` object.
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @param raster Optional raster for land areas in the report (see
#'   [m_report()]).
#' @param land Optional land polygons for the map background.
#' @return The paths of the written files, invisibly.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' f <- write_m_area(m, tempdir(), "M_example")
#' basename(f)
write_m_area <- function(x, dir, prefix = "M", raster = NULL, land = NULL) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, s))
  gpkg <- f(".gpkg")
  if (file.exists(gpkg)) unlink(gpkg)
  lay <- function(g, name, ...) {
    if (!length(g)) return(invisible())
    sf::st_write(sf::st_sf(name = name, ...,
                           geometry = sf::st_transform(g, 4326)),
                 gpkg, layer = name, quiet = TRUE)
  }
  if (inherits(x, "m_scenarios")) {
    for (nm in names(x$areas)) {
      lay(x$areas[[nm]]$m, nm,
          area_km2 = round(geom_area(x$areas[[nm]]$m) / 1e6))
    }
    lk <- x$linking
    first <- x$areas[[1]]
    utils::write.csv(x$comparison, f("_comparison.csv"), row.names = FALSE)
    paths <- c(gpkg = gpkg, comparison = f("_comparison.csv"))
  } else {
    check_m(x)
    lay(x$m, "M", area_km2 = round(geom_area(x$m) / 1e6))
    lay(x$automatic, "M_automatic",
        area_km2 = round(geom_area(x$automatic) / 1e6))
    lay(x$added, "expert_additions")
    lay(x$corridors, "corridors")
    lk <- x$linking
    first <- x
    rep <- m_report(x, raster = raster)
    utils::write.csv(rep, f("_report.csv"), row.names = FALSE)
    paths <- c(gpkg = gpkg, report = f("_report.csv"))
  }
  pts <- first$points
  pts_sf <- sf::st_transform(sf::st_as_sf(pts[setdiff(names(pts), ".group")],
                                          coords = c("x", "y"),
                                          crs = first$settings$crs), 4326)
  pts_sf$group <- pts$.group
  sf::st_write(pts_sf, gpkg, layer = "records", quiet = TRUE)
  utils::write.csv(lk$table, f("_connectivity.csv"), row.names = FALSE)
  utils::write.csv(data.frame(rank_from_longest = rev(seq_len(nrow(lk$links))),
                              link_km = round(lk$links$km, 2)),
                   f("_links.csv"), row.names = FALSE)
  ggplot2::ggsave(f("_map.png"), plot_m_area(x, land = land), width = 9,
                  height = 7, dpi = 150)
  ggplot2::ggsave(f("_links.png"), plot_linking(lk), width = 7, height = 4.5,
                  dpi = 150)
  invisible(c(paths, connectivity = f("_connectivity.csv"),
              links = f("_links.csv"), map = f("_map.png"),
              links_plot = f("_links.png")))
}
