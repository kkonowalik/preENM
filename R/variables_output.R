# Output of the variable selection -----------------------------------------------

#' @export
print.var_selection <- function(x, n = 5, advice = TRUE, ...) {
  s <- x$settings
  cat("Variable selection: ", nrow(x$variables), " variables; screen |rho| < ",
      s$thr_r, ", VIF <= ", s$thr_vif, " in ", length(x$rho), " data sets",
      if (length(x$regional_rho)) paste0(" (regions: |rho| < ",
                                         s$thr_r_regional, ", VIF <= ",
                                         s$thr_vif_regional, " in ",
                                         length(x$regional_rho), " matrices)"),
      "\n", sep = "")
  cat("Criteria: ", paste(names(s$criteria), s$criteria, sep = ":",
                          collapse = ", "), "\n\n", sep = "")
  # the main result first
  if (!is.null(x$sets)) {
    cat("BEST SETS per size (by mean rank; * = Pareto front):\n")
    b <- best_sets(x, n = n)
    b$set <- paste0(ifelse(b$pareto, "* ", "  "), b$set)
    b$pareto <- NULL
    num <- vapply(b, is.double, logical(1))
    b[num] <- lapply(b[num], round, 3)
    print(b, row.names = FALSE)
  } else {
    cat("NO SET passes the screen.\n")
  }
  if (advice) {
    cat("\nIN SHORT\n")
    cat(bullets(recommend_sets(x, print = FALSE)), sep = "\n")
  }
  cat("\nDETAILS\n")
  cl <- split(x$clusters$var, x$clusters$cluster)
  cat("Clusters of the union tree (|rho| >= ", s$thr_r, "): ",
      paste(vapply(cl, paste, character(1), collapse = "+"), collapse = " | "),
      "\n\n", sep = "")
  cat("Sets surviving each filter:\n")
  cn <- x$counts
  dead <- cn$pairs_compatible == 0
  print(cn[!dead, , drop = FALSE], row.names = FALSE)
  if (any(dead)) {
    cat("  (no compatible combination of ", range_text(cn$size[dead]),
        " variables)\n", sep = "")
  }
  if (!is.null(x$regions)) {
    cat("\nRegions, all variables (limits |rho| < ", s$thr_r_regional,
        ", VIF <= ", s$thr_vif_regional, "; see correlation_by_region()):\n",
        sep = "")
    rg <- x$regions[x$regions$screened,
                    c("data_set", "region", "cells", "max_abs_rho", "pair",
                      "pairs_above")]
    rg$max_abs_rho <- round(rg$max_abs_rho, 3)
    print(rg, row.names = FALSE)
    off <- x$regions[!x$regions$screened, , drop = FALSE]
    if (nrow(off)) {
      cat("  Not screened: ",
          paste0(off$data_set, " | ", off$region, " (", off$cells, " cells",
                 if ("why_not" %in% names(off)) paste0(", ", off$why_not),
                 ")", collapse = "; "), "\n", sep = "")
    }
  }
  invisible(x)
}

# "4 to 13", "5 and 6" or "4, 6, 9" (internal)
range_text <- function(v) {
  v <- sort(unique(v))
  if (length(v) > 2 && all(diff(v) == 1)) paste(v[1], "to", v[length(v)]) else
    if (length(v) == 2) paste(v, collapse = " and ") else
      paste(v, collapse = ", ")
}

# Sentences as a wrapped bullet list (internal)
bullets <- function(x) {
  vapply(x, function(z) {
    paste(strwrap(z, width = 78, initial = "- ", prefix = "  "),
          collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

#' The outcome of a variable selection in plain words
#'
#' Reads the tables of [select_variables()] for you: which set to start
#' from and how good it is, the close alternatives, the sizes for which
#' nothing passes and why, the variables that the good sets share and those
#' that never pass, the correlations that change between data sets, and
#' what the regional screen excluded. It is printed with the selection
#' (under "IN SHORT"). The suggestion is the best-ranked set of the largest
#' size that passes; it is a starting point, and the choice among the top
#' sets should be made on ecological grounds.
#'
#' @param x A `var_selection` from [select_variables()].
#' @param print Print the sentences (otherwise only return them).
#' @return The sentences, as a character vector (invisibly when printed).
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' e1 <- sim$env[["env1"]]
#' e2 <- sim$env[["env2"]]
#' now <- c(e1, e2, e1 + e2 + 0.05 * e1 * e2, e1 * e2, sqrt(e2))
#' names(now) <- c("t", "p", "tp", "txp", "sp")
#' past <- now - 0.2 * now[["t"]]
#' names(past) <- names(now)
#' es <- sample_environment(list(present = now, past = past), n = 3000,
#'                          seed = 1, verbose = FALSE)
#' vs <- select_variables(es, sizes = 2:3)
#' recommend_sets(vs)
recommend_sets <- function(x, print = TRUE) {
  check_vs(x)
  s <- x$settings
  f <- function(v) format(round(v, 2), nsmall = 2)
  out <- character(0)
  r <- x$sets
  if (is.null(r) || !nrow(r)) {
    out <- c(out, paste0("No set passes the screen. ",
                         why_no_sets(x, x$counts$size),
                         " Relax thr_r / thr_vif, use fewer or other ",
                         "variables, or smaller sizes."))
  } else {
    k <- max(r$n_vars)
    top <- r[r$n_vars == k, , drop = FALSE]
    b <- top[1, ]
    reg <- if (!is.na(b$regional_maxr)) paste0(
      ", ", f(b$regional_maxr), " inside a region (", b$worst_region, ")")
    out <- c(out, paste0(
      "Suggested set: ", b$set, " (", k, " variables, the largest size that ",
      "passes", if (!b$pareto) "; not on the Pareto front", "). Largest |rho| ",
      f(b$worst_maxr), " over the data sets", reg, "; VIF ", f(b$worst_vif),
      if (!is.na(b$novelty)) paste0("; ", round(100 * b$novelty, 1),
                                    "% of the target cells outside the ",
                                    "reference conditions"), "."))
    alt <- top[-1, , drop = FALSE]
    alt <- utils::head(alt[alt$pareto, , drop = FALSE], 3)
    if (nrow(alt)) {
      out <- c(out, paste0("Close alternatives of that size (Pareto front): ",
                           paste(alt$set, collapse = "; "), "."))
    }
    if (nrow(top) == 1) {
      out <- c(out, paste0("It is the only set of ", k, " variables that ",
                           "passes."))
    }
    smaller <- setdiff(sort(unique(r$n_vars), decreasing = TRUE), k)
    if (length(smaller)) {
      out <- c(out, paste0(
        "Best of the smaller sizes: ",
        paste(vapply(smaller, function(j) {
          paste0(r$set[r$n_vars == j][1], " (", j, ")")
        }, character(1)), collapse = "; "), "."))
    }
    none <- setdiff(x$counts$size, r$n_vars)
    if (length(none)) {
      out <- c(out, paste0("No set of ", range_text(none), " variables ",
                           "passes. ", why_no_sets(x, none)))
    }
    fq <- x$frequency
    often <- fq$var[fq$in_pareto_sets > 0]
    never <- fq$var[fq$in_passing_sets == 0]
    out <- c(out, paste0(
      "Most often in the Pareto sets: ",
      paste(utils::head(often, 4), collapse = ", "),
      if (length(never)) paste0(". Never in a passing set: ",
                                paste(never, collapse = ", ")), "."))
  }
  # correlations that change between the data sets
  p <- x$pairs
  own <- grep("^rho_", names(p), value = TRUE)
  own <- own[!grepl("@common$", own)]
  ch <- p[p$collinear_only_in_some, , drop = FALSE]
  if (nrow(ch) && length(own) > 1) {
    spread <- apply(abs(as.matrix(ch[own])), 1, function(z) max(z) - min(z))
    e <- ch[which.max(spread), ]
    v <- abs(unlist(e[own]))
    out <- c(out, paste0(
      nrow(ch), " pair", if (nrow(ch) != 1) "s are" else " is",
      " collinear in some data sets and not in others, e.g. ", e$var1, " ~ ",
      e$var2, ": |rho| ", f(max(v)), " in ", sub("^rho_", "", names(v)[which.max(v)]),
      ", ", f(min(v)), " in ", sub("^rho_", "", names(v)[which.min(v)]),
      ". A set chosen on one data set alone could contain such a pair."))
  }
  # what the regional screen adds
  if (length(x$regional_rho)) {
    only <- p[p$max_abs < s$thr_r & !is.na(p$max_abs_regional) &
                p$max_abs_regional >= s$thr_r_regional, , drop = FALSE]
    if (nrow(only)) {
      tb <- sort(table(only$regional_driver), decreasing = TRUE)
      out <- c(out, paste0(
        nrow(only), " pair", if (nrow(only) != 1) "s are" else " is",
        " excluded only by the regional screen (below ", s$thr_r,
        " over the whole area, at or above ", s$thr_r_regional,
        " inside a region), most often in ",
        paste0(names(tb)[seq_len(min(2, length(tb)))], " (",
               tb[seq_len(min(2, length(tb)))], ")", collapse = " and "),
        ". Check that these regions are meaningful in those data sets (for ",
        "example not under an ice sheet in that period); otherwise mask ",
        "them in sample_environment() or use `skip_regions`."))
    }
  }
  off <- if (is.null(x$regions)) NULL else
    x$regions[!x$regions$screened, , drop = FALSE]
  if (!is.null(off) && nrow(off)) {
    small <- unique(off$region[!"why_not" %in% names(off) |
                                 off$why_not == "too few cells"])
    if (length(small)) {
      out <- c(out, paste0(
        "Regions with too few sampled cells to be screened in at least one ",
        "data set: ", paste(small, collapse = ", "), ". Draw more cells ",
        "(`n`) or lower `min_region_n` to include them."))
    }
  }
  out <- c(out, paste0("The ranking is a guide: choose among the top sets ",
                       "on ecological grounds."))
  if (print) {
    cat(bullets(out), sep = "\n")
    return(invisible(out))
  }
  out
}

#' Dendrogram of the variables, with the best sets below it
#'
#' The union tree of [select_variables()]: complete linkage on 1 - the largest
#' |Spearman rho| over all data sets, so variables joined below the dashed
#' line (|rho| >= `thr_r`) are collinear in at least one data set; these
#' clusters are shaded alternately. Below the tree, each row is one variable
#' set and the dots mark its variables: a passing set takes at most one
#' variable from each cluster. By default the rows are the best-ranked sets
#' of each size (largest size first), labelled with their size and rank,
#' and the variables of the first row are printed in bold.
#'
#' @param x A `var_selection`.
#' @param set Sets to show as rows: `NULL` (default) = the `n_sets`
#'   best-ranked sets of each size; one set as a character vector of
#'   variable names; or several sets, each one string with names separated
#'   by spaces (as in [best_sets()]), optionally named to label the rows.
#' @param n_sets Number of best sets per size shown when `set` is `NULL`
#'   (at most 10 rows in total).
#' @param base_size Base font size.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?select_variables; e.g. plot_variable_tree(vs), or with your own
#' # choice next to the best set:
#' # plot_variable_tree(vs, c(best = best_sets(vs)$set[1],
#' #                          mine = "bio05 bio15 bio17"))
plot_variable_tree <- function(x, set = NULL, n_sets = 2, base_size = 11) {
  check_vs(x)
  hc <- x$tree
  thr <- x$settings$thr_r
  # the sets to show, as a named list of variable vectors
  if (is.null(set)) {
    b <- x$sets
    sets <- list()
    if (!is.null(b) && nrow(b)) {
      b <- b[b$rank <= n_sets, , drop = FALSE]
      b <- utils::head(b[order(-b$n_vars, b$rank, b$mean_rank), ], 10)
      sets <- strsplit(b$set, " ")
      names(sets) <- paste0(b$n_vars, " var., #", b$rank)
    }
    what <- "the best-ranked sets of each size (size, #rank)"
  } else {
    sets <- if (length(set) > 1 && !any(grepl(" ", set))) {
      list(set)
    } else strsplit(set, " ")
    names(sets) <- if (!is.null(names(set)) && length(set) == length(sets)) {
      names(set)
    } else if (length(sets) == 1) "set" else paste("set", seq_along(sets))
    bad <- setdiff(unlist(sets), hc$labels)
    if (length(bad)) stop("Unknown variable(s): ",
                          paste(bad, collapse = ", "), call. = FALSE)
    what <- "the given sets"
  }
  seg <- dendro_segments(hc)
  o <- hc$labels[hc$order]
  lab <- data.frame(x = seq_along(o), label = o)
  lab$chosen <- if (length(sets)) lab$label %in% sets[[1]] else FALSE
  # layout below the tree: variable names, then one row per set
  top <- max(hc$height)
  band <- (0.03 * max(nchar(o)) + 0.04) * max(top, 0.5)
  step <- 0.07 * max(top, 0.5)
  row_y <- -band - step * seq_along(sets)
  bottom <- if (length(sets)) min(row_y) - step / 2 else -band
  # clusters of collinear variables, as alternating bands
  cl <- stats::cutree(hc, h = 1 - thr)[o]
  run <- cumsum(c(1, diff(as.integer(factor(cl, unique(cl)))) != 0))
  bands <- data.frame(xmin = tapply(lab$x, run, min) - 0.5,
                      xmax = tapply(lab$x, run, max) + 0.5,
                      shade = seq_along(unique(run)) %% 2 == 0)
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = bands[bands$shade, , drop = FALSE],
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax),
                       ymin = bottom, ymax = 1 - thr, fill = "#eeede8") +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = .data$x, y = .data$y,
                                                   xend = .data$xend,
                                                   yend = .data$yend),
                          colour = "#5c5c58", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = 1 - thr, linetype = "dashed",
                        colour = "#c9480c", linewidth = 0.4) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = .data$x, y = -0.015 * top,
                                                label = .data$label,
                                                colour = .data$chosen),
                       angle = 90, hjust = 1, size = base_size / 3.2,
                       fontface = ifelse(lab$chosen, "bold", "plain"))
  if (length(sets)) {
    grid <- do.call(rbind, lapply(seq_along(sets), function(i) {
      data.frame(x = lab$x, y = row_y[i], inset = lab$label %in% sets[[i]],
                 first = i == 1)
    }))
    rows <- data.frame(y = row_y, label = names(sets),
                       first = seq_along(sets) == 1)
    p <- p +
      ggplot2::geom_segment(data = rows,
                            ggplot2::aes(y = .data$y, yend = .data$y),
                            x = 0.6, xend = length(o) + 0.4,
                            colour = "#d8d7d2", linewidth = 0.3) +
      ggplot2::geom_point(data = grid[!grid$inset, , drop = FALSE],
                          ggplot2::aes(.data$x, .data$y), shape = 21,
                          fill = "white", colour = "#c4c3be", size = 1.6) +
      ggplot2::geom_point(data = grid[grid$inset, , drop = FALSE],
                          ggplot2::aes(.data$x, .data$y,
                                       colour = .data$first), size = 3) +
      ggplot2::geom_text(data = rows,
                         ggplot2::aes(y = .data$y, label = .data$label,
                                      colour = .data$first),
                         x = 0.35, hjust = 1, size = base_size / 3.4)
  }
  left <- if (length(sets)) 0.35 - 0.2 * max(nchar(names(sets))) *
    length(o) / 13 else 0.4
  p +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#1c5cab",
                                            `FALSE` = "#5c5c58"),
                                 guide = "none") +
    ggplot2::scale_x_continuous(limits = c(left, length(o) + 0.6),
                                breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(bottom, top * 1.05),
                                breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = NULL, y = "1 - max |rho| over the data sets",
                  title = "Variables: union tree and the best sets",
                  subtitle = paste0(
                    "Dashed: |rho| = ", thr, "; variables joined below it are ",
                    "collinear in at least one data set (bands)",
                    if (length(sets)) paste0(
                      "\nRows: ", what, "; dots: their variables. Bold ",
                      "and blue: the first row"))) +
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
#' Two views of the passing sets, one panel per set size:
#'
#' * `type = "scores"` (default): a score card. Rows are the best-ranked
#'   sets, columns the criteria and the largest correlations; each cell
#'   shows the value, shaded by how good it is among all passing sets of
#'   that size (darker = better). Sets on the Pareto front are marked `*`.
#' * `type = "scatter"`: every passing set in the space of two criteria.
#'   Pareto-front sets are filled blue, and the best-ranked sets carry
#'   their rank, with the key below the plot. Scores that take few
#'   distinct values (e.g. `gradient_min`, which is decided by one
#'   variable) put the sets on lines; choose other axes with `x_metric`
#'   and `y_metric`, or use the score card.
#'
#' @param x A `var_selection`.
#' @param x_metric,y_metric Scores on the axes of the scatter (default: the
#'   two criteria with the most distinct values). Giving one switches to
#'   the scatter.
#' @param n_labels Number of best-ranked sets labelled per size in the
#'   scatter.
#' @param base_size Base font size.
#' @param type `"scores"` or `"scatter"`.
#' @param n Number of best-ranked sets per size in the score card.
#' @return A `ggplot` object.
#' @export
#' @examples
#' # see ?select_variables; e.g. plot_variable_sets(vs) and
#' # plot_variable_sets(vs, type = "scatter")
plot_variable_sets <- function(x, x_metric = NULL, y_metric = NULL,
                               n_labels = 5, base_size = 10,
                               type = c("scores", "scatter"), n = 8) {
  check_vs(x)
  type <- if (missing(type) && (!is.null(x_metric) || !is.null(y_metric))) {
    "scatter"
  } else match.arg(type)
  r <- x$sets
  if (is.null(r) || !nrow(r)) stop("No passing sets to plot.", call. = FALSE)
  crit <- x$settings$criteria
  ax <- c(names(crit), "worst_vif", "worst_maxr")
  miss <- setdiff(c(x_metric, y_metric), names(r))
  if (length(miss)) stop("Unknown score(s): ", paste(miss, collapse = ", "),
                         call. = FALSE)
  direction <- function(m) {
    if (m %in% names(crit)) crit[[m]] else
      if (m %in% c("worst_vif", "worst_maxr", "regional_maxr",
                   "regional_vif")) "min" else "none"
  }
  sizes <- paste(sort(unique(r$n_vars)), "variables")
  r$size <- factor(paste(r$n_vars, "variables"), levels = sizes)

  if (type == "scores") {
    cols <- unique(c(names(crit), "worst_maxr",
                     if (!all(is.na(r$regional_maxr))) "regional_maxr"))
    cols <- cols[vapply(cols, function(m) !all(is.na(r[[m]])), logical(1))]
    d <- do.call(rbind, lapply(cols, function(m) {
      # 1 = the best value among the passing sets of that size
      good <- stats::ave(r[[m]], r$n_vars, FUN = function(v) {
        rg <- range(v, na.rm = TRUE)
        if (diff(rg) == 0) return(rep(0.5, length(v)))
        g <- (v - rg[1]) / diff(rg)
        if (direction(m) == "min") 1 - g else g
      })
      data.frame(set = r$set, size = r$size, rank = r$rank,
                 mean_rank = r$mean_rank, pareto = r$pareto, metric = m,
                 value = r[[m]], good = good)
    }))
    d <- d[d$rank <= n, , drop = FALSE]
    # at most n rows per size also when ranks are tied
    keep <- unlist(lapply(split(r, r$size), function(z) utils::head(z$set, n)))
    d <- d[d$set %in% keep, , drop = FALSE]
    d$row <- paste0(ifelse(d$pareto, "* ", "  "), d$set, "   #", d$rank)
    ord <- unique(d[order(d$size, -d$mean_rank), "row"])
    d$row <- factor(d$row, levels = ord)
    head_lab <- vapply(cols, function(m) {
      paste0(m, "\n", switch(direction(m), min = "(lower is better)",
                             max = "(higher is better)", ""))
    }, character(1))
    d$metric <- factor(head_lab[d$metric], levels = head_lab)
    d$label <- formatC(d$value, digits = 3, format = "fg", flag = "#")
    d$dark <- d$good > 0.6
    return(
      ggplot2::ggplot(d, ggplot2::aes(.data$metric, .data$row)) +
        ggplot2::geom_tile(ggplot2::aes(fill = .data$good), colour = "white",
                           linewidth = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = .data$label,
                                        colour = .data$dark),
                           size = base_size / 3.4) +
        ggplot2::scale_fill_gradient(low = "#f4f3ef", high = "#1c5cab",
                                     limits = c(0, 1), breaks = c(0, 1),
                                     labels = c("worst", "best"),
                                     name = "Among the\npassing sets\nof the size") +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "white",
                                                `FALSE` = "#3a3a38"),
                                     guide = "none") +
        ggplot2::scale_x_discrete(position = "top") +
        ggplot2::facet_grid(size ~ ., scales = "free_y", space = "free_y") +
        ggplot2::labs(x = NULL, y = NULL,
                      title = "Best variable sets and their scores",
                      subtitle = paste0(
                        "The ", n, " best sets per size by mean rank (#); ",
                        "* = on the Pareto front")) +
        theme_enm(base_size, "#3a3a38", "#8a8a86") +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(
                         hjust = 0, family = "mono", colour = "#3a3a38"),
                       strip.text.y = ggplot2::element_text(angle = 0))
    )
  }

  # scatter: axes with the most distinct values, unless chosen
  spread <- vapply(ax, function(m) {
    if (all(is.na(r[[m]]))) 0 else length(unique(signif(r[[m]], 4)))
  }, numeric(1))
  pick <- names(sort(-spread))
  if (is.null(x_metric)) x_metric <- setdiff(pick, y_metric)[1]
  if (is.null(y_metric)) y_metric <- setdiff(pick, x_metric)[1]
  r$.x <- r[[x_metric]]
  r$.y <- r[[y_metric]]
  r$front <- factor(ifelse(r$pareto, "Pareto front", "other"),
                    levels = c("Pareto front", "other"))
  lab <- do.call(rbind, lapply(split(r, r$size), utils::head, n_labels))
  lab$model <- as.character(lab$rank)
  key <- vapply(split(lab, lab$size, drop = TRUE), function(z) {
    paste0(z$size[1], ": ", paste0(z$rank, " = ", z$set, collapse = ";  "))
  }, character(1))
  dl <- function(m) axis_label(m, direction(m))
  p <- ggplot2::ggplot(r, ggplot2::aes(.data$.x, .data$.y)) +
    ggplot2::geom_point(ggplot2::aes(fill = .data$front), shape = 21,
                        colour = "white", size = 2.8, stroke = 0.5,
                        alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("Pareto front" = "#1c5cab",
                                          other = "#b4b3ae"), name = NULL,
                               drop = FALSE) +
    ggplot2::facet_wrap(~size, scales = "free")
  if (nrow(lab)) p <- p + label_layer(lab, "#3a3a38", "#8a8a86",
                                      base_size * 1.15)
  p + ggplot2::labs(x = dl(x_metric), y = dl(y_metric),
                    title = "Variable sets passing the collinearity screen",
                    subtitle = paste0("Numbers: rank by mean rank within ",
                                      "the size (key below)"),
                    caption = paste(key, collapse = "\n")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "top",
                   plot.caption = ggplot2::element_text(hjust = 0,
                                                        colour = "#3a3a38"),
                   plot.caption.position = "plot")
}

#' Write a variable selection to files
#'
#' Writes CSV tables (`<prefix>_best_sets.csv` with the ten best sets per
#' size, `_sets.csv` with every passing set, `_summary.txt` with the
#' sentences of [recommend_sets()], `_regions.csv`, `_regions_best_set.csv`,
#' `_pareto.csv`, `_counts.csv`, `_pairs.csv`, `_variables.csv`,
#' `_frequency.csv`, `_reference_sets.csv`, `_consistency.csv`,
#' `_structure_agreement.csv`, `_union_rho.csv`, `_tree.nwk`) and PNG figures
#' (`_tree.png` with the best sets below the tree, `_union_heatmap.png`,
#' `_sets.png` with the score card and `_sets_scatter.png`).
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
  if (!is.null(x$sets)) csv(best_sets(x, n = 10), "_best_sets.csv")
  writeLines(recommend_sets(x, print = FALSE), f("_summary.txt"))
  paths[["_summary.txt"]] <- f("_summary.txt")
  if (!is.null(x$sets)) csv(x$sets[x$sets$pareto, ], "_pareto.csv")
  csv(x$counts, "_counts.csv")
  csv(x$pairs, "_pairs.csv")
  csv(x$regions, "_regions.csv")
  if (!is.null(x$regions) && !is.null(x$sets) && nrow(x$sets)) {
    csv(correlation_by_region(x, x$sets$set[1]), "_regions_best_set.csv")
  }
  csv(x$variables, "_variables.csv")
  csv(x$frequency, "_frequency.csv")
  csv(x$reference_sets, "_reference_sets.csv")
  csv(x$consistency, "_consistency.csv")
  csv(x$structure_agreement, "_structure_agreement.csv", rn = TRUE)
  csv(x$union, "_union_rho.csv", rn = TRUE)
  writeLines(hclust_newick(x$tree), f("_tree.nwk"))
  paths[["_tree.nwk"]] <- f("_tree.nwk")
  ggplot2::ggsave(f("_tree.png"), plot_variable_tree(x), width = 8.5,
                  height = 6.5, dpi = 150)
  ggplot2::ggsave(f("_union_heatmap.png"), plot_correlation(x), width = 7.5,
                  height = 7, dpi = 150)
  paths[["_tree.png"]] <- f("_tree.png")
  paths[["_union_heatmap.png"]] <- f("_union_heatmap.png")
  if (!is.null(x$sets) && nrow(x$sets)) {
    nrow_card <- sum(pmin(table(x$sets$n_vars), 8))
    ggplot2::ggsave(f("_sets.png"), plot_variable_sets(x), width = 10,
                    height = max(3.5, 1.6 + 0.3 * nrow_card), dpi = 150,
                    limitsize = FALSE)
    ggplot2::ggsave(f("_sets_scatter.png"),
                    plot_variable_sets(x, type = "scatter"), width = 10,
                    height = 6.5, dpi = 150)
    paths[["_sets.png"]] <- f("_sets.png")
    paths[["_sets_scatter.png"]] <- f("_sets_scatter.png")
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
