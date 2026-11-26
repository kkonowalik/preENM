# How many background points? ----------------------------------------------------

#' How many background points are enough to represent the environment of M?
#'
#' [sample_background()] draws `n = 10000` points by default, yet the
#' number and the spatial spread of background points change
#' presence-background models (Whitford et al. 2024), and the best number
#' depends on the species and the study area (Rausell-Moreno et al. 2025).
#' `background_sufficiency()` checks the part of the question that does not
#' depend on the model: how many points are needed for the background to
#' represent the environment available in the calibration area. It draws
#' background at an increasing series of `n`, several times each, and
#' measures:
#'
#' * **Coverage of environmental space**: the environment of M is a PCA of
#'   a large reference sample of its cells (`ref_n`), cut into about
#'   `total_bins` bins on the first PCs (as in [sample_background()]'s
#'   environmental methods). Coverage is the share of the occupied bins
#'   that hold at least one background point (`coverage`), and the share
#'   of the area of M in covered bins (`coverage_area`). Rare
#'   environments, in bins with few cells, are the last to be covered.
#' * **Stability of the marginal distributions**: the largest
#'   Kolmogorov-Smirnov distance between the background and the reference
#'   values of any layer (`ks_max`). It falls roughly as one over the
#'   square root of `n` for a uniform background; for a biased (effort or
#'   target-group) background it levels off at the size of the bias.
#'
#' **Recommended n**: the smallest `n` whose mean coverage is within `tol`
#' (default 1 percentage point) of the coverage at the largest `n` tested.
#' When coverage still rises by more than `tol` between the two largest
#' `n`, the curve has not saturated within the range tested, and the
#' printout says so.
#'
#' This is about representing M, not about the accuracy of a particular
#' model: in simulations, the best background was about 0.5-1 times the
#' number of presences or about 5% of the cells of the area, depending on
#' prevalence, and the number of background points mattered little for
#' projections to future climates (Rausell-Moreno et al. 2025). Use the
#' saturation point as a defensible lower bound for representing the
#' environment, and test the models themselves if they are sensitive.
#'
#' @param raster Predictor layers clipped to M ([terra::SpatRaster] or
#'   files), as for [sample_background()].
#' @param ns Numbers of background points to try.
#' @param n_rep Replicate draws per `n` (different seeds).
#' @param ref_n Cells in the reference sample of M (all cells if fewer).
#' @param total_bins Approximate number of environmental bins.
#' @param env_cumvar,env_axes PCs used for the bins: those explaining
#'   `env_cumvar` of the variance (at most 3), or exactly `env_axes`.
#' @param tol Tolerance of the saturation (coverage, 0-1).
#' @param n_presences Optional number of presences, to report the ratio of
#'   background to presences.
#' @param log_vars Layers log-transformed for the PCA (default: bioclim
#'   precipitation amounts).
#' @param seed Random seed of the reference sample and the draws.
#' @param ... Further arguments for [sample_background()] (e.g. `method`,
#'   `effort`, `target`, `exclude`), so that the background you will use is
#'   the one assessed.
#' @return An object of class `background_sufficiency`: `draws` (one row
#'   per `n` and replicate), `summary` (per `n`: mean and SD of the
#'   measures), `recommended_n`, `saturated`, `reference` (bins, PCs) and
#'   `settings`. See [plot_background_sufficiency()].
#' @references Whitford, A.M., Shipley, B.R., McGuire, J.L. (2024) The
#'   influence of the number and distribution of background points in
#'   presence-background species distribution models. *Ecological
#'   Modelling* 488, 110604. \doi{10.1016/j.ecolmodel.2023.110604}
#'
#'   Rausell-Moreno, A., Galiana, N., Naimi, B., Araujo, M.B. (2025)
#'   Improving species distribution models by optimising background
#'   points: impacts on current and future climate projections.
#'   *Ecological Modelling* 507, 111177.
#'   \doi{10.1016/j.ecolmodel.2025.111177}
#' @seealso [sample_background()]
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1)
#' bs <- background_sufficiency(sim$env, ns = c(100, 250, 500, 1000, 2500),
#'                              n_rep = 3, total_bins = 400, seed = 1)
#' bs
#' plot_background_sufficiency(bs)
background_sufficiency <- function(raster,
                                   ns = c(250, 500, 1000, 2500, 5000, 10000,
                                          25000, 50000),
                                   n_rep = 5, ref_n = 1e5, total_bins = 1000,
                                   env_cumvar = 0.95, env_axes = NULL,
                                   tol = 0.01, n_presences = NULL,
                                   log_vars = NULL, seed = 1, ...) {
  r <- as_grid(raster)
  ns <- sort(unique(as.integer(ns)))
  if (length(ns) < 2) stop("Give at least two values in `ns`.",
                           call. = FALSE)
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  # reference sample of the environment of M
  ref <- quiet(sample_background(r, n = ref_n, seed = seed))
  avail <- nrow(ref)
  if (avail < max(ns)) {
    message("Only ", avail, " cells are available: n is capped there.")
    ns <- unique(pmin(ns, avail))
  }
  R <- env_at_cells(r, ref$cell)
  vars <- colnames(R)
  R <- R[stats::complete.cases(R), , drop = FALSE]
  sp <- env_space(R, total_bins, env_cumvar, env_axes, log_vars)
  k <- sp$k
  nb <- sp$nb
  ve <- sp$ve
  ref_bin <- sp$ref_bin
  ref_tab <- table(ref_bin)
  ref_share <- ref_tab / sum(ref_tab)
  ref_ecdf <- lapply(vars, function(v) stats::ecdf(R[, v]))
  names(ref_ecdf) <- vars
  ref_u <- lapply(vars, function(v) sort(unique(R[, v])))
  names(ref_u) <- vars
  # two-sample Kolmogorov-Smirnov distance, exact with ties: both
  # step functions compared at every value of either sample
  ks <- function(X) {
    max(vapply(vars, function(v) {
      x <- X[, v]
      if (!length(x)) return(NA_real_)
      u <- sort(unique(c(ref_u[[v]], x)))
      max(abs(stats::ecdf(x)(u) - ref_ecdf[[v]](u)))
    }, numeric(1)))
  }

  rows <- list()
  for (n in ns) for (i in seq_len(n_rep)) {
    bg <- quiet(sample_background(r, n = n, seed = seed * 1000 + i, ...))
    B <- env_at_cells(r, bg$cell)
    B <- B[stats::complete.cases(B), , drop = FALSE]
    covered <- unique(sp$bin_of(B))
    inb <- names(ref_tab) %in% covered
    rows[[length(rows) + 1]] <- data.frame(
      n = n, replicate = i, n_drawn = nrow(bg),
      coverage = mean(inb), coverage_area = sum(ref_share[inb]),
      ks_max = ks(B))
  }
  draws <- do.call(rbind, rows)
  agg <- function(v, f) as.vector(tapply(draws[[v]], draws$n, f))
  summ <- data.frame(n = ns,
                     coverage = agg("coverage", mean),
                     coverage_sd = agg("coverage", stats::sd),
                     coverage_area = agg("coverage_area", mean),
                     ks_max = agg("ks_max", mean),
                     ks_sd = agg("ks_max", stats::sd))
  m <- length(ns)
  plateau <- summ$coverage[m]
  every_cell <- max(ns) >= avail
  rec <- ns[which(summ$coverage >= plateau - tol)[1]]
  saturated <- (summ$coverage[m] - summ$coverage[m - 1]) <= tol
  method <- list(...)$method
  structure(list(draws = draws, summary = summ, recommended_n = rec,
                 saturated = saturated, every_cell = every_cell,
                 reference = list(cells = nrow(R), bins = length(ref_tab),
                                  axes = k, bins_per_axis = nb,
                                  variance = sum(ve[seq_len(k)]),
                                  rare_bins = sum(ref_tab <= 2)),
                 settings = list(ns = ns, n_rep = n_rep, tol = tol,
                                 method = if (is.null(method)) "uniform"
                                 else method,
                                 n_presences = n_presences)),
            class = "background_sufficiency")
}

#' @export
print.background_sufficiency <- function(x, digits = 3, ...) {
  s <- x$summary
  st <- x$settings
  rf <- x$reference
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  pct <- function(v) paste0(formatC(100 * v, format = "f", digits = 1), "%")
  f <- function(v) formatC(v, format = "d", big.mark = " ")
  rec <- x$recommended_n
  i <- match(rec, s$n)
  cat("Background sufficiency (", st$method, " background): n from ",
      f(min(s$n)), " to ", f(max(s$n)), ", ", st$n_rep,
      " draws each\n\nIN SHORT\n", sep = "")
  if (x$saturated) {
    say(paste0("Coverage of the environmental space of M saturates at about ",
               "n = ", f(rec), ": ", pct(s$coverage[i]), " of its ",
               rf$bins, " environmental bins (", pct(s$coverage_area[i]),
               " of its area) hold a background point, within ",
               st$tol * 100, " percentage point(s) of the ", f(max(s$n)),
               "-point background (", pct(s$coverage[nrow(s)]), ")."))
    if (10000 %in% s$n || max(s$n) >= 10000) {
      say(if (rec <= 10000) {
        paste0("The usual default of 10 000 is ",
               if (rec < 10000) "more than " else "", "enough to ",
               "represent the environment.")
      } else {
        paste0("The usual default of 10 000 is not enough here: use at ",
               "least ", f(rec), ".")
      })
    }
  } else {
    say(paste0("Coverage has not saturated within the range tested: it ",
               "still rises from ", pct(s$coverage[nrow(s) - 1]), " to ",
               pct(s$coverage[nrow(s)]), " between n = ",
               f(s$n[nrow(s) - 1]), " and ", f(max(s$n)), ". Test larger ",
               "n, or accept that the rarest environments (", rf$rare_bins,
               " bins with at most 2 of the reference cells) stay ",
               "partly uncovered."))
  }
  if (isTRUE(x$every_cell)) {
    say(paste0("The largest n uses every available cell of M, so its ",
               "coverage is complete by construction; the recommended n is ",
               "relative to the whole area."))
  }
  say(paste0("The marginal distributions differ from those of M by at most ",
             formatC(s$ks_max[i], format = "f", digits = 3),
             " (Kolmogorov-Smirnov distance, worst layer) at that n, ",
             formatC(s$ks_max[nrow(s)], format = "f", digits = 3),
             " at the largest.",
             if (st$method != "uniform") paste0(
               " A ", st$method, " background is biased on purpose, so ",
               "this distance does not fall to 0.") else ""))
  if (!is.null(st$n_presences)) {
    say(paste0("With ", f(st$n_presences), " presences this is ",
               formatC(rec / st$n_presences, format = "f", digits = 1),
               " background points per presence; in simulations, about ",
               "0.5-1 per presence gave the most accurate current ",
               "predictions (Rausell-Moreno et al. 2025), a different ",
               "question from representing M."))
  }
  say(paste0("The number and spread of background points change model ",
             "predictions in a species-dependent way (Whitford et al. ",
             "2024); test your models at the chosen n and at a larger ",
             "one if the result matters."))
  cat("\n")
  show <- s
  show[, -1] <- lapply(show[, -1], round, digits)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "coverage = share of the environmental bins of M (PC1-PC", rf$axes,
    ", ", round(100 * rf$variance), "% of the variance, ", rf$bins_per_axis,
    " bins per axis, from ", f(rf$cells), " reference cells) holding at ",
    "least one background point; coverage_area = share of the area of M in ",
    "covered bins; ks_max = largest Kolmogorov-Smirnov distance between ",
    "background and reference values of any layer; _sd over ", st$n_rep,
    " draws. Saturation = within ", st$tol, " of the coverage at the largest ",
    "n: a rule of thumb, not a value from the literature."), width = 76),
    sep = "\n")
  invisible(x)
}

#' Plot the saturation curve of background coverage
#'
#' Coverage of the environmental space of M (bins and area) and the
#' largest Kolmogorov-Smirnov distance of the marginal distributions
#' against the number of background points, with every draw and the mean,
#' from [background_sufficiency()].
#'
#' @param x A `background_sufficiency` object.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to; the format follows the
#'   extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Whitford, A.M., Shipley, B.R., McGuire, J.L. (2024) The
#'   influence of the number and distribution of background points in
#'   presence-background species distribution models. *Ecological
#'   Modelling* 488, 110604. \doi{10.1016/j.ecolmodel.2023.110604}
#' @export
#' @examples
#' # see ?background_sufficiency
plot_background_sufficiency <- function(x, base_size = 11, file = NULL,
                                        width = 7, height = 5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "background_sufficiency")) {
    stop("`x` must come from background_sufficiency().", call. = FALSE)
  }
  d <- x$draws
  s <- x$summary
  panels <- c("Coverage of the environment of M (%)",
              "Largest Kolmogorov-Smirnov distance")
  lv <- c("Environmental bins", "Area", "Marginal distributions")
  long <- function(df, v, panel, series, scale = 1) {
    data.frame(n = df$n, value = df[[v]] * scale, panel = panel,
               series = series, stringsAsFactors = FALSE)
  }
  pts <- rbind(long(d, "coverage", panels[1], lv[1], 100),
               long(d, "coverage_area", panels[1], lv[2], 100),
               long(d, "ks_max", panels[2], lv[3]))
  lin <- rbind(long(s, "coverage", panels[1], lv[1], 100),
               long(s, "coverage_area", panels[1], lv[2], 100),
               long(s, "ks_max", panels[2], lv[3]))
  for (q in c("pts", "lin")) {
    z <- get(q)
    z$panel <- factor(z$panel, panels)
    z$series <- factor(z$series, lv)
    assign(q, z)
  }
  # ColorBrewer Dark2
  pal <- stats::setNames(c("#1b9e77", "#7570b3", "#d95f02"), lv)
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = x$recommended_n, linetype = "dashed",
                        colour = "#8a8a86") +
    ggplot2::geom_point(data = pts, ggplot2::aes(.data$n, .data$value,
                                                 colour = .data$series),
                        alpha = 0.35, size = 1.3) +
    ggplot2::geom_line(data = lin, ggplot2::aes(.data$n, .data$value,
                                                colour = .data$series,
                                                linetype = .data$series),
                       linewidth = 0.7) +
    ggplot2::scale_linetype_manual(values = stats::setNames(
      c("solid", "22", "solid"), lv), name = NULL) +
    ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_log10(labels = function(v) formatC(v, format = "d",
                                                         big.mark = " ")) +
    ggplot2::scale_colour_manual(values = pal, name = NULL) +
    ggplot2::labs(x = "Number of background points", y = NULL,
                  title = "Background coverage against its size",
                  subtitle = paste0(
                    x$settings$method, " background, ", x$settings$n_rep,
                    " draws per n (points) and their mean (lines);\n",
                    "dashed: n = ", formatC(x$recommended_n, format = "d",
                                            big.mark = " "),
                    ", where coverage is within ", x$settings$tol * 100,
                    " percentage point(s) of the largest n")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom")
}

# Environmental space of a reference sample (internal): a PCA of the layers
# (precipitation amounts log-transformed) cut into about `total_bins`
# equal-width bins on the first PCs (those explaining `cumvar`, at most 3,
# or exactly `axes`). `bin_of()` gives the bin of any matrix of layer
# values (outside the reference range: the outermost bin).
env_space <- function(R, total_bins = 1000, cumvar = 0.95, axes = NULL,
                      log_vars = NULL) {
  vars <- colnames(R)
  lv <- if (is.null(log_vars)) default_log_vars(vars) else
    intersect(log_vars, vars)
  tf <- function(M) log_tf(pmax_cols(M, lv), lv)
  pc <- stats::prcomp(tf(R), scale. = TRUE)
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  k <- if (!is.null(axes)) min(axes, ncol(R)) else
    min(3, which(cumsum(ve) >= cumvar - 1e-9)[1])
  nb <- max(4, round(total_bins^(1 / k)))
  S <- pc$x[, seq_len(k), drop = FALSE]
  brk <- lapply(seq_len(k), function(j) {
    seq(min(S[, j]), max(S[, j]), length.out = nb + 1)
  })
  bin_scores <- function(Sc) {
    do.call(paste, lapply(seq_len(k), function(j) {
      pmin(pmax(findInterval(Sc[, j], brk[[j]], rightmost.closed = TRUE),
                1), nb)
    }))
  }
  bin_of <- function(X) {
    X <- X[, vars, drop = FALSE]
    bin_scores(stats::predict(pc, tf(X))[, seq_len(k), drop = FALSE])
  }
  list(bin_of = bin_of, ref_bin = bin_scores(S), k = k, nb = nb, ve = ve,
       pca = pc)
}
