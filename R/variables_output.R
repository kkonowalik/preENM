# Output of the variable selection -----------------------------------------------

#' @export
print.var_selection <- function(x, n = 5, ...) {
  s <- x$settings
  cat("Variable selection: ", nrow(x$variables), " variables; screen |rho| < ",
      s$thr_r, ", VIF <= ", s$thr_vif, " in ", length(x$rho), " data sets",
      if (length(x$regional_rho)) paste0(" (regions: |rho| < ",
                                         s$thr_r_regional, ", VIF <= ",
                                         s$thr_vif_regional, " in ",
                                         length(x$regional_rho), " matrices)"),
      "\n", sep = "")
  cl <- split(x$clusters$var, x$clusters$cluster)
  cat("Clusters of the union tree (|rho| >= ", s$thr_r, "): ",
      paste(vapply(cl, paste, character(1), collapse = "+"), collapse = " | "),
      "\n", sep = "")
  cat("Criteria: ", paste(names(s$criteria), s$criteria, sep = ":",
                          collapse = ", "), "\n\n", sep = "")
  cat("Sets surviving each filter:\n")
  print(x$counts, row.names = FALSE)
  if (!is.null(x$sets)) {
    cat("\nBest sets per size (by mean rank; * = Pareto front):\n")
    b <- best_sets(x, n = n)
    b$set <- paste0(ifelse(b$pareto, "* ", "  "), b$set)
    b$pareto <- NULL
    num <- vapply(b, is.double, logical(1))
    b[num] <- lapply(b[num], round, 3)
    print(b, row.names = FALSE)
  }
  invisible(x)
}

#' Dendrogram of the variables
#'
#' The union tree of [select_variables()]: complete linkage on 1 - the largest
#' |Spearman rho| over all data sets, so variables joined below the dashed
#' line (|rho| >= `thr_r`) are collinear in at least one data set. Give a
#' `set` to highlight its variables (e.g. one of [best_sets()]).
#'
#' @param x A `var_selection`.
#' @param set Optional variables to highlight (a character vector, or one
#'   string with names separated by spaces).
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?select_variables; e.g. plot_variable_tree(vs, best_sets(vs)$set[1])
plot_variable_tree <- function(x, set = NULL, base_size = 11) {
  check_vs(x)
  hc <- x$tree
  if (length(set) == 1) set <- strsplit(set, " ")[[1]]
  seg <- dendro_segments(hc)
  lab <- data.frame(x = seq_along(hc$order), label = hc$labels[hc$order])
  lab$chosen <- lab$label %in% set
  thr <- x$settings$thr_r
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, y = .data$y,
                                                   xend = .data$xend,
                                                   yend = .data$yend),
                          colour = "#5c5c58", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 1 - thr, linetype = "dashed",
                        colour = "#c9480c", linewidth = 0.4) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = -0.015,
                                                label = .data$label,
                                                colour = .data$chosen),
                       angle = 90, hjust = 1, size = base_size / 3.2,
                       fontface = ifelse(lab$chosen, "bold", "plain")) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#1c5cab",
                                            `FALSE` = "#5c5c58"),
                                 guide = "none") +
    ggplot2::scale_y_continuous(limits = c(-0.025 * max(nchar(lab$label)) -
                                             0.03, max(hc$height) * 1.05),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = NULL, y = "1 - max |rho| over the data sets",
                  title = "Variables: union tree (complete linkage)",
                  subtitle = paste0(
                    "Dashed: |rho| = ", thr, "; variables joined below it are ",
                    "collinear in at least one data set",
                    if (length(set)) paste0("\nHighlighted: ",
                                            paste(set, collapse = ", ")))) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(colour = "#8a8a86"),
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = NA))
}

# Segments of a dendrogram drawn from an hclust object (internal)
dendro_segments <- function(hc) {
  n <- length(hc$order)
  leaf_x <- match(seq_len(n), hc$order)
  node_x <- numeric(n - 1)
  segs <- vector("list", n - 1)
  pos <- function(j) if (j < 0) c(leaf_x[-j], 0) else c(node_x[j],
                                                       hc$height[j])
  for (i in seq_len(n - 1)) {
    a <- pos(hc$merge[i, 1])
    b <- pos(hc$merge[i, 2])
    h <- hc$height[i]
    node_x[i] <- (a[1] + b[1]) / 2
    segs[[i]] <- data.frame(x = c(a[1], b[1], a[1]), y = c(a[2], b[2], h),
                            xend = c(a[1], b[1], b[1]), yend = c(h, h, h))
  }
  do.call(rbind, segs)
}

#' Heatmap of the correlations between variables
#'
#' `which = "union"` shows the largest |rho| over all data sets (what the
#' screen uses); a data-set name shows its own signed rho. Variables are
#' ordered by the union tree.
#'
#' @param x A `var_selection`.
#' @param which `"union"` or a data-set name (see `names(x$rho)`).
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?select_variables; e.g. plot_correlation(vs)
plot_correlation <- function(x, which = "union", base_size = 10) {
  check_vs(x)
  M <- if (which == "union") x$union else {
    if (!which %in% names(x$rho)) {
      stop("`which` must be \"union\" or one of: ",
           paste(names(x$rho), collapse = ", "), call. = FALSE)
    }
    x$rho[[which]]
  }
  o <- x$tree$labels[x$tree$order]
  d <- expand.grid(v1 = o, v2 = o, stringsAsFactors = FALSE)
  d$value <- M[cbind(d$v1, d$v2)]
  d$v1 <- factor(d$v1, levels = o)
  d$v2 <- factor(d$v2, levels = rev(o))
  thr <- x$settings$thr_r
  diag_cell <- as.character(d$v1) == as.character(d$v2)
  d$over <- abs(d$value) >= thr & !diag_cell
  d$dark <- abs(d$value) > 0.55                 # text contrast on the fill
  d$label <- ifelse(diag_cell, "", sprintf("%.2f", d$value))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$v1, .data$v2,
                                       fill = .data$value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4)
  if (length(o) <= 25) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label,
                                             colour = .data$dark,
                                             fontface = ifelse(.data$over,
                                                               "bold",
                                                               "plain")),
                                size = base_size / 3.6) +
      ggplot2::scale_colour_manual(values = c(`TRUE` = "white",
                                              `FALSE` = "#3a3a38"),
                                   guide = "none")
  }
  p <- if (which == "union") {
    p + ggplot2::scale_fill_gradient(low = "#f4f3ef", high = "#104281",
                                     limits = c(0, 1), name = "max |rho|")
  } else {
    p + ggplot2::scale_fill_gradient2(low = "#b2182b", mid = "#f4f3ef",
                                      high = "#104281", limits = c(-1, 1),
                                      name = "rho")
  }
  p + ggplot2::labs(x = NULL, y = NULL,
                    title = if (which == "union")
                      "Largest |Spearman rho| over the data sets" else
                        paste("Spearman rho:", which),
                    subtitle = paste0("Bold: |rho| >= ", thr,
                                      " (collinear)")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5),
                   panel.grid = ggplot2::element_blank(),
                   plot.subtitle = ggplot2::element_text(colour = "#8a8a86"),
                   plot.background = ggplot2::element_rect(fill = "white",
                                                           colour = NA))
}

#' Plot the scored variable sets
#'
#' Every passing set in the space of two criteria, one panel per set size;
#' Pareto-front sets are filled and the best-ranked sets are labelled.
#'
#' @param x A `var_selection`.
#' @param x_metric,y_metric Scores on the axes (default: the first two
#'   criteria).
#' @param n_labels Number of best-ranked sets labelled per size.
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?select_variables; e.g. plot_variable_sets(vs)
plot_variable_sets <- function(x, x_metric = NULL, y_metric = NULL,
                               n_labels = 3, base_size = 10) {
  check_vs(x)
  r <- x$sets
  if (is.null(r) || !nrow(r)) stop("No passing sets to plot.", call. = FALSE)
  crit <- x$settings$criteria
  ax <- c(names(crit), "worst_vif", "worst_maxr")
  if (is.null(x_metric)) x_metric <- ax[1]
  if (is.null(y_metric)) y_metric <- setdiff(ax, x_metric)[1]
  miss <- setdiff(c(x_metric, y_metric), names(r))
  if (length(miss)) stop("Unknown score(s): ", paste(miss, collapse = ", "),
                         call. = FALSE)
  r$.x <- r[[x_metric]]
  r$.y <- r[[y_metric]]
  r$size <- factor(paste(r$n_vars, "variables"),
                   levels = paste(sort(unique(r$n_vars)), "variables"))
  r$front <- factor(ifelse(r$pareto, "Pareto front", "other"),
                    levels = c("Pareto front", "other"))
  lab <- r[r$rank <= n_labels, , drop = FALSE]
  dl <- function(m) {
    d <- if (m %in% names(crit)) crit[[m]] else
      if (m %in% c("worst_vif", "worst_maxr")) "min" else "none"
    axis_label(m, d)
  }
  p <- ggplot2::ggplot(r, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$front), shape = 21,
                        colour = "white", size = 2.6, stroke = 0.5) +
    ggplot2::scale_fill_manual(values = c("Pareto front" = "#1c5cab",
                                          other = "#b4b3ae"), name = NULL) +
    ggplot2::facet_wrap(~size, scales = "free")
  if (nrow(lab)) p <- p + label_layer(cbind(lab, model = lab$set), "#3a3a38",
                                      "#8a8a86", base_size)
  p + ggplot2::labs(x = dl(x_metric), y = dl(y_metric),
                    title = "Variable sets passing the collinearity screen",
                    subtitle = paste0("Labelled: the ", n_labels,
                                      " best sets by mean rank per size")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86")
}

#' Write a variable selection to files
#'
#' Writes CSV tables (`<prefix>_sets.csv` with every passing set,
#' `_pareto.csv`, `_counts.csv`, `_pairs.csv`, `_variables.csv`,
#' `_frequency.csv`, `_reference_sets.csv`, `_consistency.csv`,
#' `_structure_agreement.csv`, `_union_rho.csv`, `_tree.nwk`) and PNG figures
#' (`_tree.png` with the best set highlighted, `_union_heatmap.png`,
#' `_sets.png`).
#'
#' @param x A `var_selection`.
#' @param dir Output folder.
#' @param prefix File name prefix.
#' @return The paths of the written files, invisibly.
#' @export
#' @examples
#' # see ?select_variables; e.g. write_variable_selection(vs, tempdir())
write_variable_selection <- function(x, dir, prefix = "variables") {
  check_vs(x)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f <- function(s) file.path(dir, paste0(prefix, s))
  paths <- character(0)
  csv <- function(d, s, rn = FALSE) {
    if (is.null(d)) return(invisible())
    utils::write.csv(d, f(s), row.names = rn)
    paths[[s]] <<- f(s)
  }
  csv(x$sets, "_sets.csv")
  if (!is.null(x$sets)) csv(x$sets[x$sets$pareto, ], "_pareto.csv")
  csv(x$counts, "_counts.csv")
  csv(x$pairs, "_pairs.csv")
  csv(x$variables, "_variables.csv")
  csv(x$frequency, "_frequency.csv")
  csv(x$reference_sets, "_reference_sets.csv")
  csv(x$consistency, "_consistency.csv")
  csv(x$structure_agreement, "_structure_agreement.csv", rn = TRUE)
  csv(x$union, "_union_rho.csv", rn = TRUE)
  writeLines(hclust_newick(x$tree), f("_tree.nwk"))
  paths[["_tree.nwk"]] <- f("_tree.nwk")
  best <- if (!is.null(x$sets)) x$sets$set[1] else NULL
  ggplot2::ggsave(f("_tree.png"), plot_variable_tree(x, best), width = 8,
                  height = 5.5, dpi = 150)
  ggplot2::ggsave(f("_union_heatmap.png"), plot_correlation(x), width = 7.5,
                  height = 7, dpi = 150)
  paths[["_tree.png"]] <- f("_tree.png")
  paths[["_union_heatmap.png"]] <- f("_union_heatmap.png")
  if (!is.null(x$sets) && nrow(x$sets)) {
    ggplot2::ggsave(f("_sets.png"), plot_variable_sets(x), width = 10,
                    height = 6, dpi = 150)
    paths[["_sets.png"]] <- f("_sets.png")
  }
  invisible(paths)
}

# Newick string of an hclust tree (internal; no extra package needed)
hclust_newick <- function(hc) {
  lab <- gsub("[(),:; ]", "_", hc$labels)
  node <- function(i) {
    ch <- hc$merge[i, ]
    s <- vapply(ch, function(j) {
      if (j < 0) {
        paste0(lab[-j], ":", format(hc$height[i] / 2, digits = 4))
      } else {
        paste0(node(j), ":", format((hc$height[i] - hc$height[j]) / 2,
                                    digits = 4))
      }
    }, character(1))
    paste0("(", paste(s, collapse = ","), ")")
  }
  paste0(node(nrow(hc$merge)), ";")
}

check_vs <- function(x) {
  if (!inherits(x, "var_selection")) {
    stop("`x` must come from select_variables().", call. = FALSE)
  }
}
