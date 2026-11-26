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
#' @param env Optional predictor layers: adds, as attribute
#'   `"availability"`, how evenly M offers the records' conditions of each
#'   predictor ([m_availability()]).
#' @return A data.frame with one row per part.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' m_report(m, raster = sim$env)
m_report <- function(m, raster = NULL, env = NULL) {
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
  if (!is.null(env)) {
    attr(out, "availability") <- m_availability(
      env, m$all_points[c("lon", "lat")], m = m)
  }
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
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' plot_m_area(m)
plot_m_area <- function(x, land = NULL, points = TRUE, base_size = 11,
                        file = NULL, width = 8, height = 7) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
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
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF. The plot is also
#'   returned.
#' @param width,height Size of the saved plot in inches.
#' @return A `ggplot` object.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' lk <- linking_distance(sim$presence, thin_km = 5, link_search = c(20, 200))
#' plot_linking(lk)
#' plot_linking(lk, "clusters")
plot_linking <- function(x, type = c("links", "clusters"), n_links = 40,
                         base_size = 11,
                         file = NULL, width = 8, height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width,
                     height))
  }
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

# Availability of conditions within M -------------------------------------------

#' How evenly does M offer the conditions the species occupies?
#'
#' Presence records come from the conditions that M offers, and M rarely
#' offers them evenly: a climate common in M will be sampled often whatever
#' the species prefers, and one that M offers in a few cells will be
#' sampled rarely even if the species thrives there. Uneven availability of
#' climates inside the calibration area therefore biases the estimated
#' niche (Jimenez & Soberon 2022). `m_availability()` is the continuous
#' analogue of [m_class_shares()]: for each predictor, it divides the
#' range of the records' values into `bins` equal bins and compares the
#' share of M's cells in each bin with the share an even spread over M's
#' range would give.
#'
#' * `ratio` (per bin, in `bins`): observed / even share; 1 = as available
#'   as an even spread, 0.1 = a tenth of it.
#' * `well_represented`: the share of the records' range in bins with
#'   `ratio` of at least `thin_ratio` (default 0.25); the rest is thinly
#'   represented.
#' * `complete`: the share of the records' range in which M has any cell
#'   (records outside M's range, or gaps).
#' * `evenness`: Pielou's evenness of M's cells over the bins of the
#'   records' range (1 = even).
#' * `distorted`: `well_represented` below `min_well` (default 0.8),
#'   `complete` below 1, or `evenness` below `min_even` (default 0.7); M
#'   then gives a distorted view of the availability of that variable.
#'
#' For a flagged variable, read its response curve with care: the records
#' are thin where M is thin. Jimenez & Soberon (2022) weight the likelihood
#' of the records by the availability of their conditions to estimate the
#' fundamental niche; other remedies are a larger or differently shaped M,
#' or a background stratified in environmental space ([sample_background()]
#' with `method = "envstrat"`).
#'
#' @param env Predictor layers ([terra::SpatRaster] or files), clipped to M
#'   or covering it.
#' @param records Records: a data.frame with coordinates, a CSV path, `sf`
#'   or [terra::SpatVector].
#' @param m Optional M: an `m_area` from [delimit_m()] or polygons; cells
#'   outside it are ignored. `NULL` when `env` is already clipped to M.
#' @param bins Number of bins over the records' range of each predictor.
#' @param thin_ratio Ratio below which a bin is thinly represented.
#' @param min_well,min_even Limits of `well_represented` and `evenness`
#'   below which a variable is flagged.
#' @param n Cells of M sampled.
#' @param coords,points_crs Coordinate columns and their CRS.
#' @param seed Random seed.
#' @return An object of class `m_availability`: `variables` (one row per
#'   predictor: records' range, `well_represented`, `thin`, `complete`,
#'   `evenness`, `distorted` and `reason`), `bins` (per predictor and bin:
#'   `lower`, `upper`, `m_share`, `even_share`, `ratio`, `records`) and
#'   `settings`. See [plot_m_availability()].
#' @references Jimenez, L., Soberon, J. (2022) Estimating the fundamental
#'   niche: accounting for the uneven availability of existing climates in
#'   the calibration area. *Ecological Modelling* 464, 109823.
#'   \doi{10.1016/j.ecolmodel.2021.109823}
#' @seealso [m_report()], [m_class_shares()], [niche_truncation()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 100)
#' av <- m_availability(sim$env, sim$presence)
#' av
#' plot_m_availability(av)
m_availability <- function(env, records, m = NULL, bins = 20,
                           thin_ratio = 0.25, min_well = 0.8,
                           min_even = 0.7, n = 1e5, coords = NULL,
                           points_crs = "EPSG:4326", seed = 1) {
  r <- as_grid(env)
  vars <- names(r)
  # cells of M
  cells <- suppressWarnings(suppressMessages(sample_background(r, n = n,
                                                              seed = seed)))
  if (!is.null(m)) {
    g <- if (inherits(m, "m_area")) m$m else as_sfc(m, "m")
    pts <- sf::st_as_sf(data.frame(x = cells$lon, y = cells$lat),
                        coords = c("x", "y"), crs = 4326)
    pts <- sf::st_transform(pts, sf::st_crs(g))
    cells <- cells[lengths(sf::st_intersects(pts, g)) > 0, , drop = FALSE]
    if (!nrow(cells)) stop("No cell of `env` lies inside `m`.",
                           call. = FALSE)
  }
  M <- env_at_cells(r, cells$cell)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  # records
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  p <- if (is.data.frame(records) && !inherits(records, "sf")) {
    df <- as.data.frame(records, stringsAsFactors = FALSE)
    xy <- if (is.null(coords)) guess_coords(names(df)) else coords
    x <- as_number(df[[xy[1]]])
    y <- as_number(df[[xy[2]]])
    ok <- !is.na(x) & !is.na(y)
    terra::vect(cbind(x[ok], y[ok]), crs = points_crs)
  } else as_points(records, coords, points_crs)
  pr <- if (terra::same.crs(p, r)) p else terra::project(p, r)
  E <- env_at_cells(r, terra::cellFromXY(r, terra::crds(pr)))
  E <- E[stats::complete.cases(E), , drop = FALSE]
  if (nrow(E) < 5) stop("Fewer than 5 records with values in all layers.",
                        call. = FALSE)
  rows <- list()
  bin_rows <- list()
  for (v in vars) {
    e <- E[, v]
    mv <- M[, v]
    lo <- min(e)
    hi <- max(e)
    if (hi <= lo) {
      rows[[v]] <- data.frame(layer = v, rec_min = lo, rec_max = hi,
                              well_represented = NA_real_, thin = NA_real_,
                              complete = NA_real_, evenness = NA_real_,
                              distorted = NA, reason = "records constant",
                              stringsAsFactors = FALSE)
      next
    }
    br <- seq(lo, hi, length.out = bins + 1)
    w <- diff(br)
    mrange <- diff(range(mv))
    cnt <- tabulate(findInterval(mv[mv >= lo & mv <= hi], br,
                                 rightmost.closed = TRUE), bins)
    share <- cnt / length(mv)
    even <- if (mrange > 0) w / mrange else rep(1 / bins, bins)
    ratio <- share / even
    rcnt <- tabulate(findInterval(e, br, rightmost.closed = TRUE), bins)
    pin <- cnt / max(sum(cnt), 1)
    ev <- if (sum(cnt) > 0) {
      pp <- pin[pin > 0]
      -sum(pp * log(pp)) / log(bins)
    } else 0
    well <- mean(ratio >= thin_ratio)
    comp <- mean(cnt > 0)
    why <- c(if (comp < 1) paste0("no cells of M in ",
                                  round(100 * (1 - comp)),
                                  "% of the records' range"),
             if (well < min_well) paste0(round(100 * (1 - well)),
                                         "% of the records' range thinly ",
                                         "represented"),
             if (ev < min_even) paste0("uneven (evenness ", round(ev, 2),
                                       ")"))
    rows[[v]] <- data.frame(layer = v, rec_min = lo, rec_max = hi,
                            well_represented = well, thin = 1 - well,
                            complete = comp, evenness = ev,
                            distorted = length(why) > 0,
                            reason = paste(why, collapse = "; "),
                            stringsAsFactors = FALSE)
    bin_rows[[v]] <- data.frame(layer = v, bin = seq_len(bins),
                                lower = br[-(bins + 1)], upper = br[-1],
                                m_share = share, even_share = even,
                                ratio = ratio, records = rcnt,
                                stringsAsFactors = FALSE)
  }
  out <- list(variables = do.call(rbind, rows),
              bins = do.call(rbind, bin_rows),
              settings = list(bins = bins, thin_ratio = thin_ratio,
                              min_well = min_well, min_even = min_even,
                              n_cells = nrow(M), n_records = nrow(E)))
  rownames(out$variables) <- NULL
  rownames(out$bins) <- NULL
  class(out) <- "m_availability"
  out
}

#' @export
print.m_availability <- function(x, digits = 2, ...) {
  v <- x$variables
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  cat("Availability of the records' conditions in M: ", nrow(v),
      " predictor(s), ", s$n_records, " records, ", s$n_cells,
      " cells of M\n\nIN SHORT\n", sep = "")
  bad <- v[v$distorted %in% TRUE, , drop = FALSE]
  if (!nrow(bad)) {
    say(paste0("M offers the records' range of every predictor fairly ",
               "evenly: at least ", s$min_well * 100, "% of each range is ",
               "well represented."))
  } else {
    say(paste0("M gives a distorted view of the availability of ",
               nrow(bad), " predictor(s): ",
               paste0(bad$layer, " (", bad$reason, ")", collapse = "; "),
               "."))
    say(paste0("Where M is thin, records are thin whatever the species ",
               "prefers: read these response curves with care, and ",
               "consider weighting by availability (Jimenez & Soberon ",
               "2022), a different M, or a background stratified in ",
               "environmental space."))
  }
  ok <- v$layer[v$distorted %in% FALSE]
  if (length(ok) && nrow(bad)) {
    say(paste0("Fairly even: ", paste(ok, collapse = ", "), "."))
  }
  cat("\n")
  show <- v[c("layer", "rec_min", "rec_max", "well_represented",
              "complete", "evenness", "distorted")]
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], signif, 3)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "The records' range of each predictor is cut into ", s$bins, " equal ",
    "bins; ratio = share of M's cells in a bin / the share an even spread ",
    "over M's range would give. well_represented = share of the range in ",
    "bins with ratio of at least ", s$thin_ratio, "; complete = share with ",
    "any cell of M; evenness = Pielou's evenness of M's cells over the ",
    "bins. Flagged when well_represented < ", s$min_well, ", complete < 1 ",
    "or evenness < ", s$min_even, ": rules of thumb, not values from the ",
    "literature."), width = 76), sep = "\n")
  invisible(x)
}

#' Plot the availability of the records' conditions in M
#'
#' For each predictor, the availability ratio of every bin of the records'
#' range (share of M's cells relative to an even spread), with the limit of
#' thin representation, and the number of records per bin.
#'
#' @param x An `m_availability` from [m_availability()].
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Jimenez, L., Soberon, J. (2022) Estimating the fundamental
#'   niche: accounting for the uneven availability of existing climates in
#'   the calibration area. *Ecological Modelling* 464, 109823.
#'   \doi{10.1016/j.ecolmodel.2021.109823}
#' @export
#' @examples
#' # see ?m_availability
plot_m_availability <- function(x, base_size = 11, file = NULL, width = 8,
                                height = 5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "m_availability")) {
    stop("`x` must come from m_availability().", call. = FALSE)
  }
  b <- x$bins
  th <- x$settings$thin_ratio
  b$mid <- (b$lower + b$upper) / 2
  b$wid <- b$upper - b$lower
  lab <- c("Well represented", "Thinly represented")
  b$cls <- factor(ifelse(b$ratio >= th, lab[1], lab[2]), lab)
  cnt <- table(b$cls)
  ll <- paste0(lab, ": ", as.integer(cnt), " bins (",
               share_label(as.integer(cnt) / nrow(b)), ")")
  b$cls <- factor(b$cls, lab, ll)
  fl <- x$variables$layer[x$variables$distorted %in% TRUE]
  b$panel <- ifelse(b$layer %in% fl, paste0(b$layer, " (flagged)"), b$layer)
  b$panel <- factor(b$panel, unique(b$panel))
  b$ratio_plot <- pmax(b$ratio, 1e-3)
  base_y <- min(1e-3, min(b$ratio_plot)) / 2
  ggplot2::ggplot(b) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .data$lower, xmax = .data$upper,
                                    ymin = base_y, ymax = .data$ratio_plot,
                                    fill = .data$cls)) +
    ggplot2::geom_hline(yintercept = 1, colour = "#8a8a86") +
    ggplot2::geom_hline(yintercept = th, linetype = "dashed",
                        colour = "#d95f02") +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = stats::setNames(
      c("#7570b3", "#d95f02"), ll), name = NULL) +
    ggplot2::labs(x = "Value (range of the records)",
                  y = "Availability in M relative to an even spread",
                  title = "Availability of the records' conditions in M",
                  subtitle = paste0("Solid: even availability; dashed: ",
                                    th, " of it (thin below)")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom")
}
