# Flexible 2-D / 3-D model plots ---------------------------------------------------

#' Plot models on any two or three metrics, static or interactive
#'
#' One function for every view of a model comparison: choose the metrics
#' for the axes (e.g. bias against AUC, MAEbal against AUC, or AUC, MAEbal
#' and bias in 3-D), and whether the plot is a static `ggplot` or an
#' interactive one (package `plotly`) that can be rotated, zoomed and
#' hovered over in the RStudio viewer or a browser, and saved as an HTML
#' file.
#'
#' * With two metrics and `interactive = FALSE`, the plot is the same as
#'   [plot_selection()] (for [select_models()] output) or [plot_pareto()]
#'   (for [evaluate_models()] output), including `zoom`.
#' * With `z`, the plot is 3-D and always interactive.
#'
#' Axes are oriented so that better is always away from the origin
#' ("lower is better" axes are reversed). Colours show the selection
#' (selected / not selected) or, for a bootstrapped evaluation, the share
#' of resamples in which the model was on the Pareto front; circles are
#' Pareto-front models and squares dominated ones. Bars show the SD across
#' runs (`<metric>_sd`) or the bootstrap intervals (`<metric>_lower`,
#' `<metric>_upper`) when present.
#'
#' @param object Output of [select_models()], [evaluate_models()] or
#'   [evaluate_cv()], or any data.frame with a `model` column and one
#'   column per metric.
#' @param x,y,z Metrics (column names) for the axes. Default: the second and
#'   first criterion of `object` (e.g. MAE and AUC). `z` adds a third axis
#'   (3-D, interactive).
#' @param interactive Draw with `plotly` (default `TRUE` for 3-D, `FALSE`
#'   for 2-D).
#' @param zoom For static 2-D plots: `TRUE` to show only the region of the
#'   best models, or axis limits; see [plot_selection()].
#' @param labels Show model names.
#' @param error_bars Draw error bars (bootstrap intervals, else +/- 1 SD).
#' @param bar_colour Colour of the error bars.
#' @param directions Optional named vector of directions for metrics the
#'   package does not know, e.g. `c(fossil_AUC_deep = "max")` (`"max"`,
#'   `"min"` or `"zero"`). Known: the criteria of `object`, AUC, MAE,
#'   MAEbal, bias, TSS, Boyce and any column ending in `_AUC*`, `_sens`,
#'   `_lift` (higher is better) or `_area` (lower is better).
#' @param file Optional `.html` file to save an interactive plot to (for a
#'   static plot use [ggplot2::ggsave()]).
#' @param browse Open the saved file in the web browser.
#' @param base_size Base font size of static plots.
#' @return A `ggplot` (static) or `plotly` object (interactive). Printing an
#'   interactive plot opens it in the RStudio viewer or the browser.
#' @export
#' @examples
#' df <- data.frame(model = c("GLM", "GBM", "RF", "MAXNET", "SRE"),
#'                  AUC = c(0.918, 0.919, 0.912, 0.918, 0.701),
#'                  MAEbal = c(0.218, 0.236, 0.204, 0.332, 0.299),
#'                  bias = c(-0.04, -0.05, -0.02, 0.10, -0.18),
#'                  TSS = c(0.70, 0.71, 0.57, 0.71, 0))
#' sel <- select_models(df, criteria = c(AUC = "max", MAEbal = "min"),
#'                      plot = FALSE)
#' plot_models(sel, x = "bias", y = "AUC")
#' plot_models(sel, x = "MAEbal", y = "AUC", zoom = TRUE)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   plot_models(sel, x = "MAEbal", y = "AUC", z = "bias")    # 3-D
#' }
plot_models <- function(object, x = NULL, y = NULL, z = NULL,
                        interactive = !is.null(z), zoom = FALSE,
                        labels = TRUE, error_bars = TRUE,
                        bar_colour = "#c8553d", directions = NULL,
                        file = NULL, browse = interactive(),
                        base_size = 12) {
  if (!is.data.frame(object) || !"model" %in% names(object)) {
    stop("`object` must be a data.frame with a `model` column, e.g. from ",
         "select_models() or evaluate_models().", call. = FALSE)
  }
  crit <- attr(object, "criteria")
  if (inherits(object, "enm_cv") || is.null(crit)) {
    crit <- NULL
  }
  ax <- default_axes(crit, x, y)
  mets <- c(x = ax[["x"]], y = ax[["y"]], z = z)
  miss <- setdiff(mets, names(object))
  if (length(miss)) {
    stop("Columns not found: ", paste(miss, collapse = ", "), ". Available: ",
         paste(setdiff(names(object), "model"), collapse = ", "),
         call. = FALSE)
  }
  dirs <- vapply(mets, metric_direction, character(1), crit = crit,
                 directions = directions)

  # a plain table (or evaluate_cv output): Pareto front of the plotted
  # metrics, so points are still classified
  if (!inherits(object, c("enm_selection", "enm_evaluation"))) {
    cr <- dirs[dirs != "none"]
    names(cr) <- mets[names(cr)]
    if (!length(cr)) {
      stop("Unknown direction of ", paste(mets, collapse = ", "),
           "; give `directions`, e.g. c(", mets[1], " = \"max\").",
           call. = FALSE)
    }
    object <- select_models(as.data.frame(object), rule = "pareto",
                            min_metrics = numeric(0), criteria = cr,
                            plot = FALSE)
  }
  if (is.null(z) && !interactive) {
    if (inherits(object, "enm_evaluation")) {
      return(plot_pareto(object, mets[["x"]], mets[["y"]],
                         labels = if (labels) "all" else "none",
                         error_bars = error_bars, zoom = zoom,
                         bar_colour = bar_colour, base_size = base_size))
    }
    return(plot_selection(object, mets[["x"]], mets[["y"]], labels = labels,
                          zoom = zoom, error_bars = error_bars,
                          bar_colour = bar_colour, base_size = base_size))
  }
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Interactive and 3-D plots need the plotly package: ",
         "install.packages(\"plotly\")", call. = FALSE)
  }
  p <- plotly_models(as.data.frame(object), mets, dirs, labels,
                     is_eval = inherits(object, "enm_evaluation"),
                     bars = error_bars, bar_colour = bar_colour)
  if (!is.null(file)) {
    file <- normalizePath(file, mustWork = FALSE)
    lib <- file.path(dirname(file), paste0(tools::file_path_sans_ext(
      basename(file)), "_files"))
    had_lib <- dir.exists(lib)
    ok <- tryCatch({
      htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (ok && !had_lib) {
      unlink(lib, recursive = TRUE)     # not needed by a self-contained file
    }
    if (!ok) htmlwidgets::saveWidget(p, file, selfcontained = FALSE)
    message("Saved: ", file)
    if (browse) utils::browseURL(file)
  }
  p
}

# Direction of a metric for plotting (internal)
metric_direction <- function(m, crit = NULL, directions = NULL) {
  if (!is.null(directions) && m %in% names(directions)) {
    return(directions[[m]])
  }
  d <- dir_of(m, crit)
  if (!is.na(d)) return(d)
  if (grepl("_AUC(_|$)|_sens$|_lift$|_Boyce$|_TSS$", m)) return("max")
  if (grepl("_area$|_MAE(bal)?$", m)) return("min")
  "none"
}

# Interactive 2-D / 3-D scatter with plotly (internal)
plotly_models <- function(df, mets, dirs, labels, is_eval, bars = TRUE,
                          bar_colour = "#c8553d") {
  three_d <- "z" %in% names(mets)
  front <- if (is_eval) df$pareto_front == 1 else df$pareto
  df$.front <- ifelse(front, "Pareto front", "Dominated")
  df$.group <- if (is_eval) df$.front else
    paste(ifelse(df$selected, "Selected", "Not selected"), df$.front,
          sep = ", ")
  colours <- c("Selected, Pareto front" = "#1c5cab",
               "Selected, Dominated" = "#1c5cab",
               "Not selected, Pareto front" = "#a3a39e",
               "Not selected, Dominated" = "#a3a39e",
               "Pareto front" = "#1c5cab", "Dominated" = "#a3a39e")
  symbols <- c("Pareto front" = "circle", "Dominated" = "square")
  hover <- function(d) {
    vals <- vapply(mets, function(m) {
      paste0(m, ": ", formatC(d[[m]], digits = 3, format = "f"))
    }, character(nrow(d)))
    vals <- matrix(vals, nrow = nrow(d))
    paste0("<b>", d$model, "</b><br>",
           apply(vals, 1, paste, collapse = "<br>"),
           if (is_eval && "p_pareto" %in% names(d)) {
             paste0("<br>P(front): ", round(d$p_pareto, 2))
           },
           if (!is_eval) paste0("<br>", d$.group))
  }
  err <- function(d, m) {
    if (!bars) return(list(visible = FALSE))
    if (paste0(m, "_lower") %in% names(d)) {
      return(list(type = "data", symmetric = FALSE,
                  array = d[[paste0(m, "_upper")]] - d[[m]],
                  arrayminus = d[[m]] - d[[paste0(m, "_lower")]],
                  color = bar_colour, thickness = 1.5, width = 0))
    }
    if (paste0(m, "_sd") %in% names(d)) {
      return(list(type = "data", array = d[[paste0(m, "_sd")]],
                  color = bar_colour, thickness = 1.5, width = 0))
    }
    list(visible = FALSE)
  }
  axis <- function(k) {
    a <- list(title = axis_label(mets[[k]], dirs[[k]]), zeroline = FALSE)
    if (dirs[[k]] == "min") a$autorange <- "reversed"
    a
  }
  p <- plotly::plot_ly()
  # selected models last, so they are drawn on top
  grp <- unique(df$.group[order(grepl("^Selected", df$.group))])
  for (g in grp) {
    d <- df[df$.group == g, , drop = FALSE]
    mk <- list(symbol = symbols[[d$.front[1]]], size = if (three_d) 6 else 11,
               color = colours[[g]], line = list(color = "white", width = 1))
    if (is_eval && "p_pareto" %in% names(d)) {
      mk$color <- d$p_pareto
      mk$colorscale <- list(c(0, "#cde2fb"), c(1, "#104281"))
      mk$cmin <- 0
      mk$cmax <- 1
      mk$showscale <- g == grp[length(grp)]
      mk$colorbar <- list(title = "P(front)")
    }
    args <- list(p, data = d, x = d[[mets[["x"]]]], y = d[[mets[["y"]]]],
                 type = if (three_d) "scatter3d" else "scatter",
                 mode = if (labels) "markers+text" else "markers",
                 text = d$model, textposition = "top center",
                 hovertext = hover(d), hoverinfo = "text", name = g,
                 marker = mk, error_x = err(d, mets[["x"]]),
                 error_y = err(d, mets[["y"]]))
    if (three_d) {
      args$z <- d[[mets[["z"]]]]
      args$error_z <- err(d, mets[["z"]])
    }
    p <- do.call(plotly::add_trace, args)
  }
  title <- paste("Models:", paste(mets, collapse = " / "))
  if (three_d) {
    plotly::layout(p, title = title,
                   scene = list(xaxis = axis("x"), yaxis = axis("y"),
                                zaxis = axis("z")))
  } else {
    plotly::layout(p, title = title, xaxis = axis("x"), yaxis = axis("y"))
  }
}
