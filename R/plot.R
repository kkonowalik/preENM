# Visualisation ---------------------------------------------------------------

#' Plot models in the space of two metrics with the Pareto front
#'
#' Each model is a point; by default x = MAE and y = AUC, the comparison of
#' Konowalik & Nosol (2021), but any two metrics can be shown. Axes are
#' oriented so that **better models are towards the top-right** ("lower is
#' better" axes are reversed). The dashed line connects the models on the
#' Pareto front. With a bootstrapped evaluation, bars show the confidence
#' intervals and the fill shows `p_pareto`, the share of resamples in which
#' the model was on the front.
#'
#' By default the axes are the first two `criteria` used in
#' [evaluate_models()]. If more than two criteria were used, the front was
#' computed on all of them, so front models may look dominated in a 2-D plot
#' and no front line is drawn.
#'
#' @param x An `enm_evaluation` object from [evaluate_models()].
#' @param x_metric,y_metric Metrics (columns of `x`) for the axes. Default:
#'   the second and first criterion (MAE and AUC for the default criteria).
#' @param labels Which models to label: `"all"`, `"front"` (Pareto models only)
#'   or `"none"`. Default: all when there are 15 models or fewer, otherwise
#'   front only.
#' @param error_bars Show bootstrap confidence intervals (if available).
#' @param zoom `TRUE` zooms onto the region of the best models (see
#'   [plot_selection()]); or a list with the axis limits,
#'   `list(x = c(0.19, 0.24), y = c(0.89, 0.93))`. Models outside are named
#'   in the subtitle.
#' @param base_size Base font size passed to [ggplot2::theme_minimal()].
#' @return A `ggplot` object that can be customised further or saved with
#'   [ggplot2::ggsave()].
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' res <- evaluate_models(sim$models, sim$presence, sim$absence,
#'                        n_boot = 100, seed = 1)
#' plot_pareto(res)
#' plot_pareto(res, x_metric = "TSS")   # any other metric on an axis
#' # ggplot2::ggsave("pareto.png", plot_pareto(res), width = 7, height = 5,
#' #                 dpi = 300)
plot_pareto <- function(x, x_metric = NULL, y_metric = NULL, labels = NULL,
                        error_bars = TRUE, zoom = FALSE, base_size = 12) {
  if (!inherits(x, "enm_evaluation")) {
    stop("`x` must be the output of evaluate_models().", call. = FALSE)
  }
  crit <- attr(x, "criteria")
  if (is.null(crit)) crit <- c(AUC = "max", MAE = "min")
  ax <- default_axes(crit, x_metric, y_metric)
  x_metric <- ax[["x"]]
  y_metric <- ax[["y"]]
  df <- as.data.frame(x)
  miss <- setdiff(c(x_metric, y_metric), names(df))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  dirs <- c(x = dir_of(x_metric, crit), y = dir_of(y_metric, crit))
  dirs[is.na(dirs)] <- "none"
  if (is.null(labels)) labels <- if (nrow(df) <= 15) "all" else "front"
  labels <- match.arg(labels, c("all", "front", "none"))
  has_boot <- "p_pareto" %in% names(df)

  df$.x <- df[[x_metric]]
  df$.y <- df[[y_metric]]
  df$status <- factor(ifelse(df$pareto_front == 1, "Pareto front",
                             "Dominated"),
                      levels = c("Pareto front", "Dominated"))
  front <- df[df$pareto_front == 1, , drop = FALSE]
  ink <- "#3a3a38"
  muted <- "#8a8a86"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y))

  # staircase of the front: only meaningful when the plotted metrics are
  # exactly the criteria and both are "max"/"min"
  if (nrow(front) > 1 && setequal(names(crit), c(x_metric, y_metric)) &&
      all(dirs %in% c("max", "min"))) {
    p <- p + ggplot2::geom_path(data = front_path(front, dirs),
                                ggplot2::aes(x = .data$x, y = .data$y),
                                colour = muted, linewidth = 0.6,
                                linetype = "dashed")
  }
  bar_col <- if (!isFALSE(zoom) && !is.null(zoom)) "#cfcec9" else muted
  lo <- function(m) paste0(m, "_lower")
  up <- function(m) paste0(m, "_upper")
  if (has_boot && error_bars) {
    if (all(c(lo(y_metric), up(y_metric)) %in% names(df))) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[lo(y_metric)]],
                     ymax = .data[[up(y_metric)]]),
        width = 0, colour = bar_col, linewidth = 0.4)
    }
    if (all(c(lo(x_metric), up(x_metric)) %in% names(df))) {
      p <- p + ggplot2::geom_errorbar(
        ggplot2::aes(xmin = .data[[lo(x_metric)]],
                     xmax = .data[[up(x_metric)]]),
        orientation = "y", width = 0, colour = bar_col, linewidth = 0.4)
    }
  }
  if (has_boot) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(fill = .data$p_pareto,
                                       shape = .data$status),
                          size = 3.6, colour = "white", stroke = 0.8) +
      ggplot2::scale_fill_gradient(
        low = "#cde2fb", high = "#104281", limits = c(0, 1),
        name = "P(on Pareto front)\nacross bootstraps",
        guide = ggplot2::guide_colourbar(order = 2))
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(fill = .data$status,
                                       shape = .data$status),
                          size = 3.6, colour = "white", stroke = 0.8) +
      ggplot2::scale_fill_manual(values = c("Pareto front" = "#1c5cab",
                                            "Dominated" = "#a3a39e"),
                                 name = NULL, drop = FALSE)
  }
  # with a bootstrap, fill shows p_pareto, so the shape legend gets a neutral
  # fill; without, fill and shape both show the status and share one legend
  p <- p + ggplot2::scale_shape_manual(
    values = c("Pareto front" = 21, "Dominated" = 22), name = NULL,
    drop = FALSE, guide = if (has_boot) {
      ggplot2::guide_legend(order = 1, override.aes = list(
        fill = "#5c5c58", colour = "#5c5c58"))
    } else "legend")
  zero_ref <- c(x = "x", y = "y")[dirs == "zero"]
  if ("x" %in% zero_ref) {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = muted,
                                 linewidth = 0.4, linetype = "dotted")
  }
  if ("y" %in% zero_ref) {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = muted,
                                 linewidth = 0.4, linetype = "dotted")
  }

  hidden <- character(0)
  zoomed <- !isFALSE(zoom) && !is.null(zoom)
  if (zoomed) {
    zl <- zoom_limits(df, dirs, zoom)
    hidden <- zl$hidden
  }
  if (labels != "none") {
    lab <- if (labels == "front") front else df
    p <- p + label_layer(lab[!lab$model %in% hidden, , drop = FALSE], ink,
                         muted, base_size)
  }

  if (dirs[["x"]] == "min") p <- p + ggplot2::scale_x_reverse()
  if (dirs[["y"]] == "min") p <- p + ggplot2::scale_y_reverse()
  if (zoomed) {
    p <- p + ggplot2::coord_cartesian(xlim = zl$xlim, ylim = zl$ylim)
  }
  sub <- c(
    if (length(crit) > 2) {
      paste0("Pareto front on ", length(crit), " criteria: ",
             paste(names(crit), collapse = ", "), ".")
    },
    if (has_boot) {
      sprintf("Bars: %d%% bootstrap intervals (%d resamples).",
              round(100 * attr(x, "conf")), attr(x, "n_boot"))
    },
    better_text(dirs, x_metric, y_metric),
    hidden_text(hidden))
  p +
    ggplot2::labs(x = axis_label(x_metric, dirs[["x"]]),
                  y = axis_label(y_metric, dirs[["y"]]),
                  title = paste0("Model comparison: ", y_metric, " vs ",
                                 x_metric, if (zoomed) " (zoom)"),
                  subtitle = paste(sub, collapse = "\n")) +
    theme_enm(base_size, ink, muted)
}

# Direction of a metric: from the criteria, else the package default
dir_of <- function(m, crit) {
  if (m %in% names(crit)) return(crit[[m]])
  if (m %in% names(default_directions)) return(default_directions[[m]])
  NA_character_
}

axis_label <- function(m, dir) {
  if (is.na(dir) || dir == "none") m else paste0(m, " (", direction_label(dir), ")")
}

better_text <- function(dirs, x_metric, y_metric) {
  if (all(dirs %in% c("max", "min"))) {
    return("Better models are top-right.")
  }
  part <- function(m, d) switch(d, max = paste("higher", m),
                                min = paste("lower", m),
                                zero = paste(m, "closer to 0"), NULL)
  bits <- c(part(y_metric, dirs[["y"]]), part(x_metric, dirs[["x"]]))
  if (length(bits)) paste0("Better: ", paste(bits, collapse = ", "), ".")
}

# Staircase through the front, drawn in data units (internal)
front_path <- function(front, dirs) {
  g <- function(v, d) if (d == "max") v else -v      # "higher is better"
  gx <- g(front$.x, dirs[["x"]])
  gy <- g(front$.y, dirs[["y"]])
  o <- order(gx, -gy)
  fx <- front$.x[o]
  fy <- front$.y[o]
  n <- length(fx)
  # from each front point: move along y to the next point's level, then x
  px <- c(rbind(fx[-n], fx[-n]), fx[n])
  py <- c(rbind(fy[-n], fy[-1]), fy[n])
  data.frame(x = px, y = py)
}

#' @importFrom rlang .data
NULL

#' Plot a model selection for any pair of metrics
#'
#' Shows the models from [select_models()] in the space of two metrics.
#' Selected models are filled, Pareto-front models are circles (others are
#' squares), dashed lines mark the minimum values and bars show +/- 1 SD
#' across cross-validation runs when `<metric>_sd` columns are present. Call
#' it several times for several panels, e.g. AUC vs MAE and AUC vs fossil AUC.
#'
#' **Zoom.** When many models are close together the plot gets cluttered.
#' `zoom = TRUE` shows only the region of the best models: the 2-D Pareto
#' fronts of the two plotted metrics are added one after another until at
#' least half of the models (and at least three) are included, and the view
#' is their range. Models outside the view are named in the subtitle. Give
#' your own limits with `zoom = list(x = c(0.19, 0.24), y = c(0.89, 0.93))`.
#' `plot(sel)` draws the full plot and, next to it, the zoom (when the zoom
#' leaves out at least one model), and returns both.
#'
#' Selected models are drawn on top, so they are never hidden behind other
#' points. Pareto-front membership is computed on all `criteria` of
#' [select_models()]; with more than the two plotted criteria, a front model
#' can look dominated in the plot (noted in the subtitle). See
#' [plot_models()] for interactive and 3-D plots.
#'
#' @param x An `enm_selection` object from [select_models()].
#' @param x_metric,y_metric Columns to plot on the x and y axes. Default:
#'   the second and the first selection criterion.
#' @param labels Label the points with model names.
#' @param zoom `FALSE`, `TRUE` (region of the best models) or a list of
#'   axis limits (see Details). For `plot()`: `"auto"` adds the zoomed
#'   plot when it leaves out at least one model.
#' @param base_size Base font size.
#' @return A `ggplot` object; `plot()` returns a list with `all` and `zoom`
#'   invisibly.
#' @export
#' @examples
#' df <- data.frame(model = c("GLM", "GBM", "RF", "SRE"),
#'                  AUC = c(0.86, 0.88, 0.87, 0.71),
#'                  TSS = c(0.62, 0.65, 0.66, 0.38),
#'                  MAE = c(0.30, 0.27, 0.24, 0.41))
#' sel <- select_models(df)
#' plot_selection(sel)
#' plot_selection(sel, x_metric = "TSS")
#' plot_selection(sel, zoom = TRUE)        # region of the best models
#' plot(sel)                               # full view and zoom side by side
plot_selection <- function(x, x_metric = NULL, y_metric = NULL,
                           labels = TRUE, zoom = FALSE, base_size = 12) {
  if (!inherits(x, "enm_selection")) {
    stop("`x` must be the output of select_models().", call. = FALSE)
  }
  ax <- default_axes(attr(x, "criteria"), x_metric, y_metric)
  x_metric <- ax[["x"]]
  y_metric <- ax[["y"]]
  df <- as.data.frame(x)
  miss <- setdiff(c(x_metric, y_metric), names(df))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  df$.x <- df[[x_metric]]
  df$.y <- df[[y_metric]]
  df$status <- factor(ifelse(df$selected, "Selected", "Not selected"),
                      levels = c("Selected", "Not selected"))
  df$front <- factor(ifelse(df$pareto, "Pareto front", "Dominated"),
                     levels = c("Pareto front", "Dominated"))
  # selected models are drawn last, so other points never hide them
  df <- df[order(df$selected, df$pareto), , drop = FALSE]
  crit <- attr(x, "criteria")
  mins <- attr(x, "min_metrics")
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  dirs <- c(x = dir_of(x_metric, crit), y = dir_of(y_metric, crit))
  dirs[is.na(dirs)] <- "none"
  zoomed <- !isFALSE(zoom) && !is.null(zoom)
  # in a zoom the bars are often wider than the view: draw them lighter
  bar_col <- if (zoomed) "#cfcec9" else muted

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$.x, y = .data$.y))
  if (dirs[["x"]] == "zero") {
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = muted,
                                 linewidth = 0.4, linetype = "dotted")
  }
  if (dirs[["y"]] == "zero") {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = muted,
                                 linewidth = 0.4, linetype = "dotted")
  }
  if (x_metric %in% names(mins)) {
    p <- p + ggplot2::geom_vline(xintercept = mins[[x_metric]],
                                 linetype = "dashed", colour = muted,
                                 linewidth = 0.4)
  }
  if (y_metric %in% names(mins)) {
    p <- p + ggplot2::geom_hline(yintercept = mins[[y_metric]],
                                 linetype = "dashed", colour = muted,
                                 linewidth = 0.4)
  }
  xs <- paste0(x_metric, "_sd")
  ys <- paste0(y_metric, "_sd")
  if (ys %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$.y - .data[[ys]],
                   ymax = .data$.y + .data[[ys]]),
      width = 0, colour = bar_col, linewidth = 0.4)
  }
  if (xs %in% names(df)) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$.x - .data[[xs]],
                   xmax = .data$.x + .data[[xs]]),
      orientation = "y", width = 0, colour = bar_col, linewidth = 0.4)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(fill = .data$status,
                                     shape = .data$front),
                        size = 3.6, colour = "white", stroke = 0.8) +
    ggplot2::scale_fill_manual(values = c("Selected" = "#1c5cab",
                                          "Not selected" = "#b4b3ae"),
                               name = NULL, drop = FALSE,
                               guide = ggplot2::guide_legend(
                                 override.aes = list(shape = 21))) +
    ggplot2::scale_shape_manual(values = c("Pareto front" = 21,
                                           "Dominated" = 22),
                                name = NULL, drop = FALSE,
                                guide = ggplot2::guide_legend(
                                  override.aes = list(fill = "#5c5c58")))
  hidden <- character(0)
  if (zoomed) {
    zl <- zoom_limits(df, dirs, zoom)
    hidden <- zl$hidden
  }
  if (labels) {
    p <- p + label_layer(if (zoomed) df[!df$model %in% hidden, ] else df,
                         ink, muted, base_size)
  }
  if (dirs[["x"]] == "min") p <- p + ggplot2::scale_x_reverse()
  if (dirs[["y"]] == "min") p <- p + ggplot2::scale_y_reverse()
  if (zoomed) {
    p <- p + ggplot2::coord_cartesian(xlim = zl$xlim, ylim = zl$ylim)
  }
  extra_crit <- setdiff(names(crit), c(x_metric, y_metric))
  sub <- c(
    paste0("Selection rule: ", attr(x, "rule"), ". ",
           better_text(dirs, x_metric, y_metric)),
    paste0(if (length(mins)) "Dashed lines: minimum values. ",
           if (ys %in% names(df) || xs %in% names(df))
             "Bars: \u00b11 SD across runs."),
    if (length(extra_crit)) {
      paste0("Pareto front on ", length(crit), " criteria (",
             paste(names(crit), collapse = ", "), ").")
    },
    hidden_text(hidden))
  p +
    ggplot2::labs(x = axis_label(x_metric, dirs[["x"]]),
                  y = axis_label(y_metric, dirs[["y"]]),
                  title = paste0("Model selection: ", y_metric, " vs ",
                                 x_metric, if (zoomed) " (zoom)"),
                  subtitle = paste(sub[nzchar(sub)], collapse = "\n")) +
    theme_enm(base_size, ink, muted)
}

#' @rdname plot_selection
#' @param y Ignored.
#' @param ... Passed to [plot_selection()].
#' @export
plot.enm_selection <- function(x, y, x_metric = NULL, y_metric = NULL,
                               zoom = "auto", ...) {
  p_all <- plot_selection(x, x_metric, y_metric, ...)
  out <- list(all = p_all)
  if (!isFALSE(zoom) && !is.null(zoom)) {
    p_zoom <- plot_selection(x, x_metric, y_metric, zoom = TRUE, ...)
    # draw the zoom only when it leaves out at least one model
    if (length(zoom_hidden(x, x_metric, y_metric)) ||
        !identical(zoom, "auto")) {
      out$zoom <- p_zoom
    }
  }
  print_side_by_side(out)
  invisible(out)
}

#' @rdname plot_pareto
#' @param y Ignored.
#' @param ... Passed to [plot_pareto()].
#' @export
plot.enm_evaluation <- function(x, y, x_metric = NULL, y_metric = NULL,
                                zoom = "auto", ...) {
  out <- list(all = plot_pareto(x, x_metric, y_metric, ...))
  if (!isFALSE(zoom) && !is.null(zoom)) {
    p_zoom <- plot_pareto(x, x_metric, y_metric, zoom = TRUE, ...)
    if (length(zoom_hidden(x, x_metric, y_metric)) ||
        !identical(zoom, "auto")) {
      out$zoom <- p_zoom
    }
  }
  print_side_by_side(out)
  invisible(out)
}

# Print one or two ggplots next to each other (internal)
print_side_by_side <- function(plots) {
  if (length(plots) == 1) {
    print(plots[[1]])
    return(invisible())
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  for (i in 1:2) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = 1,
                                          layout.pos.col = i))
  }
  invisible()
}

# Default axes: the first criterion on y, the second on x (internal)
default_axes <- function(crit, x_metric = NULL, y_metric = NULL) {
  if (is.null(crit)) crit <- c(AUC = "max", MAE = "min")
  if (is.null(y_metric)) y_metric <- names(crit)[1]
  if (is.null(x_metric)) {
    x_metric <- if (length(crit) > 1) names(crit)[2] else
      setdiff(c("MAE", "AUC"), y_metric)[1]
  }
  c(x = x_metric, y = y_metric)
}

# Models a zoomed view would leave out (internal)
zoom_hidden <- function(x, x_metric, y_metric) {
  crit <- attr(x, "criteria")
  ax <- default_axes(crit, x_metric, y_metric)
  df <- as.data.frame(x)
  df$.x <- df[[ax[["x"]]]]
  df$.y <- df[[ax[["y"]]]]
  dirs <- c(x = dir_of(ax[["x"]], crit), y = dir_of(ax[["y"]], crit))
  dirs[is.na(dirs)] <- "none"
  zoom_limits(df, dirs, TRUE)$hidden
}

hidden_text <- function(hidden) {
  if (!length(hidden)) return(character(0))
  paste0("Outside this view: ", paste(hidden, collapse = ", "), ".")
}

# Axis limits of the region of the best models (internal). The 2-D Pareto
# fronts of the plotted metrics are added one by one until at least half
# of the models (and at least three) are included; the view is their range,
# padded. `zoom` may instead be list(x = , y = ) limits.
zoom_limits <- function(df, dirs, zoom) {
  if (is.list(zoom)) {
    xl <- zoom$x
    yl <- zoom$y
    if (is.null(xl)) xl <- range(df$.x, na.rm = TRUE)
    if (is.null(yl)) yl <- range(df$.y, na.rm = TRUE)
  } else {
    ok <- !is.na(df$.x) & !is.na(df$.y)
    d <- df[ok, , drop = FALSE]
    cr <- c(x = dirs[["x"]], y = dirs[["y"]])
    cr <- cr[cr %in% c("max", "min", "zero")]
    inside <- rep(TRUE, nrow(d))
    if (length(cr) && nrow(d) > 3) {
      m <- to_minimise(data.frame(x = d$.x, y = d$.y), cr)
      fr <- pareto_fronts(m)
      need <- max(3, ceiling(nrow(d) / 2))
      k <- 1
      while (sum(fr <= k) < need) k <- k + 1
      inside <- fr <= k
    }
    pad <- function(r) r + c(-1, 1) * max(diff(r), 1e-3) * 0.12
    xl <- pad(range(d$.x[inside]))
    yl <- pad(range(d$.y[inside]))
  }
  xl <- sort(xl)
  yl <- sort(yl)
  hidden <- df$model[is.na(df$.x) | is.na(df$.y) |
                       df$.x < xl[1] | df$.x > xl[2] |
                       df$.y < yl[1] | df$.y > yl[2]]
  # with a reversed axis the limits are given in reversed order
  list(xlim = if (dirs[["x"]] == "min") rev(xl) else xl,
       ylim = if (dirs[["y"]] == "min") rev(yl) else yl,
       hidden = as.character(hidden))
}

label_layer <- function(df, ink, muted, base_size) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = df, ggplot2::aes(label = .data$model),
                             colour = ink, size = base_size / 3.6,
                             box.padding = 0.5, min.segment.length = 0.2,
                             segment.colour = muted, segment.size = 0.3,
                             seed = 1, max.overlaps = Inf)
  } else {
    ggplot2::geom_text(data = df, ggplot2::aes(label = .data$model),
                       colour = ink, size = base_size / 3.6, vjust = -1.1)
  }
}

theme_enm <- function(base_size, ink, muted) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#e6e5e1",
                                               linewidth = 0.3),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      text = ggplot2::element_text(colour = ink),
      plot.subtitle = ggplot2::element_text(colour = muted),
      legend.position = "right")
}
