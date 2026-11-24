# Plausibility of response curves -------------------------------------------------

#' Check the shape of response curves
#'
#' Performance metrics do not tell whether a model responds to the
#' environment in an ecologically plausible way. A model can score well and
#' still have an implausible response, which then extrapolates badly under
#' climate change (Fiorentino et al. 2025); machine-learning models in
#' particular overfit into irregular responses that conventional metrics do
#' not reveal (Chollet Ramampiandra et al. 2023). The shape of the response
#' also decides how a model behaves under novel conditions: predictions hold
#' up where new conditions fall near the optimum and fail in the tails
#' (Allyn et al. 2025), and models trained on a truncated niche extrapolate
#' worse (Yau et al. 2026). `response_check()` describes every response
#' curve so that implausible ones can be found before the model is used.
#'
#' **Input.** The response curves from your modelling software: a long
#' table with the predictor (`variable`), its value (`value`) and the
#' predicted suitability (`suitability`), optionally per `model` and `run`
#' (as exported by biomod2, ENMeval, maxnet, `predicts`, `pdp`, ...).
#' Alternatively, raw predictions with the predictor values at the same
#' points (one column per predictor and a `pred` column): the response is
#' then the mean prediction in `n_bins` classes of each predictor. That is
#' a marginal response, which mixes in the effects of correlated
#' predictors; partial-dependence or "evaluation strip" curves from the
#' modelling software are better. Bin means carry sampling noise, so for
#' them a turning point must also exceed 3 x sqrt(2) median standard errors
#' of the bin means.
#'
#' **Per curve** (variable x model x run):
#' * `turning_points`: reversals of the curve (a rise followed by a fall or
#'   the reverse) larger than `tolerance` of the curve's range, so that
#'   small numerical noise is not counted.
#' * `wiggliness`: the total variation of the curve divided by its range: 1
#'   for a monotonic curve, at most 2 for a curve with one peak, larger the
#'   more the curve goes up and down.
#' * `shape`: `"flat"` (range below `flat`), `"increasing"`,
#'   `"decreasing"` (no turning point), `"unimodal"` (one peak),
#'   `"U-shaped"` (one trough; rarely plausible for suitability) or
#'   `"irregular"` (more turning points, or wiggliness above `max_wiggle`).
#' * `optimum`: the predictor value with the highest suitability, and
#'   `optimum_position`, its position in the calibrated range (0 = lower
#'   edge, 1 = upper edge).
#' * `rising_low`, `rising_high`: suitability still rising towards the lower
#'   or upper edge of the calibrated range (by more than `tolerance` of the
#'   range over the outermost `edge` share of the range). Beyond such an
#'   edge the prediction depends on how the algorithm extrapolates
#'   (clamping, a straight continuation or a drop); this is where
#'   [extrapolation_risk()] and [environment_shift()] matter most.
#'
#' **Records.** With `occurrences` (predictor values at the presences), the
#' density of the records along each gradient is compared with the curve
#' (`agreement`, the Spearman correlation between the curve and the
#' density). With `background` too, the ratio of presence to background
#' density is used, which follows the response of the species rather than
#' the availability of conditions; without it, the density of records also
#' reflects which conditions are common. The calibrated range is the range
#' of `calibration` (default: `background` and `occurrences`, else the
#' curve itself).
#'
#' **Reading the shapes: which curve is plausible?** There is no single
#' correct shape; what to expect depends on the gradient and on how much
#' of it the data cover.
#' * Along a **direct gradient** (one acting on the physiology, such as
#'   temperature or pH) a species is expected to respond with a bell: the
#'   highest suitability at its optimum, falling towards zero beyond its
#'   tolerance; narrow or broad, symmetric or skewed, with a peak or a
#'   plateau at the optimum (Austin 2002; Fiorentino et al. 2025).
#' * Along a **resource gradient** (water, nutrients) a rise that levels off
#'   (a plateau) is plausible: more of the resource helps until it no longer
#'   limits.
#' * A **monotonic** curve is plausible when the data cover only one side of
#'   the niche, e.g. a cold-adapted species whose warm limit lies inside the
#'   area but whose cold limit lies beyond it. Then it is still rising at an
#'   edge (`rising_low`/`rising_high`), and the response beyond that edge is
#'   unknown: it is a truncated niche, not a wrong model (see
#'   [niche_truncation()]).
#' * Along an **indirect gradient** (elevation, latitude), which acts through
#'   others, any shape can arise, and the curve says little about the
#'   species.
#' * **Irregular** curves (several peaks, steps) and **U-shapes** are
#'   rarely plausible; they usually mean the model fits noise or follows
#'   where records happen to lie (Chollet Ramampiandra et al. 2023).
#'
#' The curve should not follow how common the values of a predictor are in
#' the area (the distribution of the raster, which is often unimodal
#' itself), but how much more often the species is found at a value than
#' that value is available. That is why the agreement uses the ratio of
#' presence to background density when `background` is given: without it,
#' a bell-shaped density of records can simply mirror a bell-shaped
#' availability of the conditions.
#'
#' **Projection.** With `projection` (predictor values of the area or
#' period the model is projected to: a data.frame or rasters), the share of
#' it beyond each edge where the curve is still rising is reported.
#'
#' Background drawn evenly in environmental space, as
#' `sample_background(method = "envstrat")` does, makes response curves
#' closer to the true response of the species (Broussin et al. 2024).
#'
#' @param x Response curves (long table) or raw predictions (wide table);
#'   see Details.
#' @param variable,value,suitability,model,run Column names of the long
#'   table (`model` and `run` are optional).
#' @param pred For raw predictions, the column with the predictions; the
#'   other numeric columns (except `model`, `run`) are the predictors.
#' @param occurrences Optional predictor values at the presences
#'   (data.frame with one column per variable).
#' @param background Optional predictor values at the background points.
#' @param calibration Optional predictor values defining the calibrated
#'   range (default: `background` and `occurrences`).
#' @param projection Optional predictor values of a projection: a
#'   data.frame, or rasters ([terra::SpatRaster]; a sample of cells is
#'   read).
#' @param n_bins Number of classes for raw predictions.
#' @param tolerance Smallest reversal counted as a turning point, and
#'   smallest rise at an edge, as a share of the curve's range.
#' @param flat A curve whose suitability varies less than this is flat.
#' @param max_wiggle Wiggliness above which a curve is irregular.
#' @param edge Share of the calibrated range at each end used to judge
#'   whether the curve is still rising.
#' @return A list of class `response_check`: `shapes` (one row per curve:
#'   `model`, `run`, `variable`, `shape`, `turning_points`, `wiggliness`,
#'   `range`, `optimum`, `optimum_position`, `rising_low`, `rising_high`,
#'   `agreement`, `beyond_low`, `beyond_high`), `curves` (the curves
#'   used), `density` (the record density along each gradient), `limits`
#'   (calibrated range per variable) and `settings`.
#' @references Fiorentino, D., Nunez-Riboni, I., Pierce, M.E., Oesterwind,
#'   D., Akimova, A. (2025) Improving species distribution models for
#'   climate change studies: ecological plausibility and performance
#'   metrics. *Ecological Modelling* 508, 111207.
#'   \doi{10.1016/j.ecolmodel.2025.111207}
#'
#'   Austin, M.P. (2002) Spatial prediction of species distribution: an
#'   interface between ecological theory and statistical modelling.
#'   *Ecological Modelling* 157, 101-118.
#'   \doi{10.1016/S0304-3800(02)00205-3}
#'
#'   Chollet Ramampiandra, E., Scheidegger, A., Wydler, J., Schuwirth, N.
#'   (2023) A comparison of machine learning and statistical species
#'   distribution models: quantifying overfitting supports model
#'   interpretation. *Ecological Modelling* 481, 110353.
#'   \doi{10.1016/j.ecolmodel.2023.110353}
#'
#'   Allyn, A.J. et al. (2025) Contrasting species distribution model
#'   predictability under novel temperature conditions. *Diversity and
#'   Distributions* 31, e70036. \doi{10.1111/ddi.70036}
#'
#'   Yau, E.Y.H., Bonebrake, T.C., Hughes, A.C. (2026) Optimizing species
#'   distribution models for extrapolation - guidance through algorithm,
#'   variable, and ensemble modelling decisions. *Ecological Informatics*
#'   99, 104008. \doi{10.1016/j.ecoinf.2026.104008}
#'
#'   Broussin, J., Mouchet, M., Goberville, E. (2024) Generating
#'   pseudo-absences in the ecological space improves the biological
#'   relevance of response curves in species distribution models.
#'   *Ecological Modelling* 498, 110865.
#'   \doi{10.1016/j.ecolmodel.2024.110865}
#' @seealso [plot_response_check()], [extrapolation_risk()],
#'   [environment_shift()], [sample_background()]
#' @export
#' @examples
#' v <- seq(0, 30, length.out = 60)
#' curves <- rbind(
#'   data.frame(model = "GLM", variable = "temp", value = v,
#'              suitability = exp(-(v - 15)^2 / 40)),
#'   data.frame(model = "RF", variable = "temp", value = v,
#'              suitability = pmin(1, exp(-(v - 15)^2 / 40) +
#'                                   0.3 * (sin(v) > 0.6))),
#'   data.frame(model = "GLM", variable = "rain", value = v,
#'              suitability = plogis((v - 10) / 3)))
#' set.seed(1)
#' occ <- data.frame(temp = rnorm(200, 15, 4), rain = runif(200, 5, 30))
#' rc <- response_check(curves, occurrences = occ)
#' rc
#' plot_response_check(rc)
response_check <- function(x, variable = "variable", value = "value",
                           suitability = "suitability", model = "model",
                           run = "run", pred = "pred", occurrences = NULL,
                           background = NULL, calibration = NULL,
                           projection = NULL, n_bins = 30, tolerance = 0.05,
                           flat = 0.05, max_wiggle = 3, edge = 0.1) {
  if (!is.data.frame(x)) stop("`x` must be a data.frame.", call. = FALSE)
  long <- all(c(variable, value, suitability) %in% names(x))
  if (long) {
    cur <- data.frame(
      model = if (model %in% names(x)) as.character(x[[model]]) else "model",
      run = if (run %in% names(x)) as.character(x[[run]]) else "1",
      variable = as.character(x[[variable]]),
      value = as.numeric(x[[value]]),
      suitability = as.numeric(x[[suitability]]),
      stringsAsFactors = FALSE)
  } else if (pred %in% names(x)) {
    cur <- binned_responses(x, pred, model, run, n_bins)
  } else {
    stop("`x` needs the columns ", variable, ", ", value, " and ",
         suitability, " (response curves) or a `", pred, "` column with ",
         "predictor columns (raw predictions).", call. = FALSE)
  }
  if (is.null(cur$se)) cur$se <- NA_real_
  cur <- cur[stats::complete.cases(cur[c("model", "run", "variable",
                                         "value", "suitability")]), ,
             drop = FALSE]
  if (!nrow(cur)) stop("No complete rows in `x`.", call. = FALSE)
  vars <- unique(cur$variable)

  # calibrated range per variable
  calib <- calibration
  if (is.null(calib)) {
    calib <- if (!is.null(background) || !is.null(occurrences)) {
      rbind_common(background, occurrences)
    } else NULL
  }
  limits <- do.call(rbind, lapply(vars, function(v) {
    r <- if (!is.null(calib) && v %in% names(calib)) {
      range(calib[[v]], na.rm = TRUE)
    } else range(cur$value[cur$variable == v])
    data.frame(variable = v, low = r[1], high = r[2],
               from = if (!is.null(calib) && v %in% names(calib))
                 "data" else "curve", stringsAsFactors = FALSE)
  }))

  # record density along each gradient
  dens <- NULL
  if (!is.null(occurrences)) {
    dens <- do.call(rbind, lapply(vars, function(v) {
      if (!v %in% names(occurrences)) return(NULL)
      p <- stats::na.omit(occurrences[[v]])
      if (length(p) < 3) return(NULL)
      lim <- range(c(p, cur$value[cur$variable == v]))
      dp <- stats::density(p, from = lim[1], to = lim[2], n = 128)
      y <- dp$y
      ratio <- FALSE
      if (!is.null(background) && v %in% names(background)) {
        b <- stats::na.omit(background[[v]])
        if (length(b) >= 3) {
          db <- stats::density(b, from = lim[1], to = lim[2], n = 128,
                               bw = dp$bw)
          y <- ifelse(db$y > max(db$y) * 1e-3, dp$y / db$y, NA)
          ratio <- TRUE
        }
      }
      data.frame(variable = v, value = dp$x, density = y, ratio = ratio)
    }))
  }

  # projection values
  proj <- NULL
  if (!is.null(projection)) {
    proj <- if (inherits(projection, "SpatRaster")) {
      as.data.frame(sample_complete(projection[[intersect(vars,
                                                          names(projection))]],
                                    1e5))
    } else as.data.frame(projection)
  }

  keys <- unique(cur[c("model", "run", "variable")])
  shapes <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    d <- cur[cur$model == k$model & cur$run == k$run &
               cur$variable == k$variable, ]
    lim <- unlist(limits[limits$variable == k$variable, c("low", "high")])
    # binned curves: a reversal must also exceed 3 x sqrt(2) typical
    # standard errors of the bin means (noise of a difference of two bins)
    noise <- if (all(is.na(d$se))) 0 else
      3 * sqrt(2) * stats::median(d$se, na.rm = TRUE)
    s <- curve_shape(d$value, d$suitability, lim, tolerance, flat,
                     max_wiggle, edge, noise)
    ag <- NA_real_
    if (!is.null(dens) && k$variable %in% dens$variable) {
      dd <- dens[dens$variable == k$variable & !is.na(dens$density), ]
      inside <- dd$value >= min(d$value) & dd$value <= max(d$value)
      if (sum(inside) > 5) {
        ag <- suppressWarnings(stats::cor(
          stats::approx(s$x, s$y, dd$value[inside])$y,
          dd$density[inside], method = "spearman"))
      }
    }
    bl <- bh <- NA_real_
    if (!is.null(proj) && k$variable %in% names(proj)) {
      pv <- stats::na.omit(proj[[k$variable]])
      if (length(pv)) {
        bl <- if (isTRUE(s$rising_low)) mean(pv < lim[1]) else 0
        bh <- if (isTRUE(s$rising_high)) mean(pv > lim[2]) else 0
      }
    }
    data.frame(model = k$model, run = k$run, variable = k$variable,
               shape = s$shape, turning_points = s$turning_points,
               wiggliness = s$wiggliness, range = s$range,
               optimum = s$optimum, optimum_position = s$optimum_position,
               rising_low = s$rising_low, rising_high = s$rising_high,
               agreement = ag, beyond_low = bl, beyond_high = bh,
               stringsAsFactors = FALSE)
  }))
  rownames(shapes) <- NULL
  structure(list(shapes = shapes, curves = cur, density = dens,
                 limits = limits,
                 settings = list(tolerance = tolerance, flat = flat,
                                 max_wiggle = max_wiggle, edge = edge,
                                 ratio = !is.null(background),
                                 binned = !long,
                                 projection = !is.null(projection))),
            class = "response_check")
}

# Rows of two predictor tables with their common columns (internal)
rbind_common <- function(a, b) {
  if (is.null(a)) return(as.data.frame(b))
  if (is.null(b)) return(as.data.frame(a))
  v <- intersect(names(a), names(b))
  rbind(as.data.frame(a)[v], as.data.frame(b)[v])
}

# Mean prediction in classes of each predictor (internal)
binned_responses <- function(x, pred, model, run, n_bins) {
  groups <- c(model, run)[c(model, run) %in% names(x)]
  vars <- setdiff(names(x)[vapply(x, is.numeric, logical(1))],
                  c(pred, groups))
  if (!length(vars)) stop("No predictor columns in `x`.", call. = FALSE)
  key <- if (length(groups)) interaction(x[groups], drop = TRUE) else
    factor(rep("model", nrow(x)))
  do.call(rbind, lapply(split(seq_len(nrow(x)), key), function(i) {
    do.call(rbind, lapply(vars, function(v) {
      z <- x[[v]][i]
      p <- x[[pred]][i]
      ok <- !is.na(z) & !is.na(p)
      if (sum(ok) < 5) return(NULL)
      br <- seq(min(z[ok]), max(z[ok]), length.out = n_bins + 1)
      b <- findInterval(z[ok], br, rightmost.closed = TRUE)
      mz <- tapply(z[ok], b, mean)
      mp <- tapply(p[ok], b, mean)
      se <- tapply(p[ok], b, function(q) {
        if (length(q) > 1) stats::sd(q) / sqrt(length(q)) else NA_real_
      })
      data.frame(
        model = if (model %in% names(x)) as.character(x[[model]][i[1]])
          else "model",
        run = if (run %in% names(x)) as.character(x[[run]][i[1]]) else "1",
        variable = v, value = as.vector(mz), suitability = as.vector(mp),
        se = as.vector(se), stringsAsFactors = FALSE)
    }))
  }))
}

# Shape of one curve (internal)
curve_shape <- function(x, y, lim, tolerance, flat, max_wiggle, edge,
                        noise = 0) {
  o <- order(x)
  x <- x[o]
  y <- y[o]
  # average repeated x
  if (anyDuplicated(x)) {
    y <- as.vector(tapply(y, x, mean))
    x <- sort(unique(x))
  }
  rg <- diff(range(y))
  out <- list(x = x, y = y, range = rg, turning_points = 0L,
              wiggliness = 1, shape = "flat",
              optimum = x[which.max(y)], optimum_position = NA_real_,
              rising_low = FALSE, rising_high = FALSE)
  span <- lim[2] - lim[1]
  if (span > 0) out$optimum_position <- (out$optimum - lim[1]) / span
  if (length(x) < 3 || rg < flat) return(out)
  tol <- max(tolerance * rg, noise)
  # turning points: reversals larger than tol (zig-zag with hysteresis)
  ext <- integer(0)
  dir <- 0
  piv <- 1
  for (i in seq_along(y)[-1]) {
    if (dir >= 0 && y[i] > y[piv]) {
      piv <- i
      dir <- if (dir == 0 && y[i] - y[1] > tol) 1 else dir
    } else if (dir <= 0 && y[i] < y[piv]) {
      piv <- i
      dir <- if (dir == 0 && y[1] - y[i] > tol) -1 else dir
    }
    if (dir == 1 && y[piv] - y[i] > tol) {
      ext <- c(ext, piv)                 # a peak
      dir <- -1
      piv <- i
    } else if (dir == -1 && y[i] - y[piv] > tol) {
      ext <- c(ext, -piv)                # a trough
      dir <- 1
      piv <- i
    } else if (dir == 0) {
      if (y[i] > y[1] + tol) {
        dir <- 1
        piv <- i
      } else if (y[i] < y[1] - tol) {
        dir <- -1
        piv <- i
      }
    }
  }
  tp <- length(ext)
  wig <- sum(abs(diff(y))) / rg
  shape <- if (tp == 0) {
    if (y[length(y)] > y[1]) "increasing" else "decreasing"
  } else if (tp == 1) {
    if (ext > 0) "unimodal" else "U-shaped"
  } else "irregular"
  if (wig > max_wiggle) shape <- "irregular"
  # still rising towards an edge of the calibrated range?
  at <- function(v) stats::approx(x, y, v, rule = 2)$y
  w <- edge * span
  rising_low <- span > 0 && at(lim[1]) - at(lim[1] + w) > tol
  rising_high <- span > 0 && at(lim[2]) - at(lim[2] - w) > tol
  out$turning_points <- tp
  out$wiggliness <- wig
  out$shape <- shape
  out$rising_low <- rising_low
  out$rising_high <- rising_high
  out
}

#' @export
print.response_check <- function(x, digits = 2, ...) {
  s <- x$shapes
  st <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  n_curves <- nrow(s)
  cat("Response check: ", n_curves, " curve(s), ",
      length(unique(s$variable)), " variable(s), ",
      length(unique(s$model)), " model(s)", if (st$binned)
        " (binned from raw predictions)", "\n\nIN SHORT\n", sep = "")
  tab <- table(factor(s$shape, c("unimodal", "increasing", "decreasing",
                                 "flat", "U-shaped", "irregular")))
  tab <- tab[tab > 0]
  say(paste0("Shapes: ", paste0(as.vector(tab), " ", names(tab),
                                collapse = ", "), "."))
  lab <- function(r) paste0(r$variable, if (length(unique(s$model)) > 1)
    paste0(" (", r$model, ")"), if (length(unique(s$run)) > 1)
      paste0(" run ", r$run))
  irr <- s[s$shape %in% c("irregular", "U-shaped"), , drop = FALSE]
  if (nrow(irr)) {
    say(paste0("Implausible shapes: ", paste(unique(lab(irr)),
                                              collapse = "; "),
               ". Irregular responses (several turning points or a curve ",
               "that goes up and down) usually mean the model fits noise, ",
               "and U-shaped responses rarely make ecological sense; check ",
               "these models before using their projections."))
  } else {
    say("No curve is irregular or U-shaped.")
  }
  edge_r <- s[s$rising_low | s$rising_high, , drop = FALSE]
  if (nrow(edge_r)) {
    side <- ifelse(edge_r$rising_low & edge_r$rising_high, "both edges",
                   ifelse(edge_r$rising_low, "lower edge", "upper edge"))
    say(paste0("Still rising at the edge of the calibrated range: ",
               paste(unique(paste0(lab(edge_r), " - ", side)),
                     collapse = "; "),
               ". The optimum may lie beyond the data (a truncated niche), ",
               "and predictions beyond that edge depend on how the ",
               "algorithm extrapolates."))
  }
  if (st$projection) {
    b <- s[(s$beyond_low > 0 | s$beyond_high > 0) %in% TRUE, , drop = FALSE]
    if (nrow(b)) {
      say(paste0("Part of the projection lies beyond such an edge: ",
                 paste(unique(paste0(lab(b), " ", share_label(
                   pmax(b$beyond_low, b$beyond_high, na.rm = TRUE)))),
                   collapse = "; "),
                 " of the projected values; see extrapolation_risk()."))
    } else {
      say("No part of the projection lies beyond an edge where a curve is still rising.")
    }
  }
  if (any(!is.na(s$agreement))) {
    low <- s[!is.na(s$agreement) & s$agreement < 0.3, , drop = FALSE]
    if (nrow(low)) {
      say(paste0("Curves that do not follow the records (agreement below ",
                 "0.3): ", paste(unique(lab(low)), collapse = "; "), ".",
                 if (!st$ratio) paste0(" Without background points the ",
                                       "record density also reflects which ",
                                       "conditions are common.")))
    } else {
      say("Every curve follows the records along its gradient (agreement 0.3 or more).")
    }
  }
  cat("\n")
  show <- s[c("model", "run", "variable", "shape", "turning_points",
              "wiggliness", "optimum", "optimum_position", "rising_low",
              "rising_high", "agreement")]
  if (length(unique(show$run)) == 1) show$run <- NULL
  if (length(unique(show$model)) == 1) show$model <- NULL
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "turning_points = reversals larger than ", 100 * st$tolerance,
    "% of the curve's range; wiggliness = total variation / range (1 = ",
    "monotonic, up to 2 = one peak); irregular = two or more turning ",
    "points or wiggliness above ", st$max_wiggle, "; flat = range below ",
    st$flat, ". rising_low/high = still rising over the outer ",
    100 * st$edge, "% of the calibrated range. optimum_position: 0 = lower ",
    "edge, 1 = upper edge. agreement = Spearman correlation between the ",
    "curve and the ", if (st$ratio) "presence/background density ratio"
    else "density of records", ". All limits are rules of thumb."),
    width = 78), sep = "\n")
  invisible(x)
}

#' Plot response curves with their shape and the density of records
#'
#' One panel per variable: the response curves of the models (thin lines
#' for runs), labelled with their shape, the calibrated range (dotted
#' lines) and, beneath the curves, the density of the records along the
#' gradient (grey), so you can see whether the modelled response follows
#' the data or the algorithm.
#'
#' @param x A result of [response_check()].
#' @param variables Optional subset of variables.
#' @param models Optional subset of models.
#' @param base_size Base font size.
#' @param file Optional file(s) to save the plot to, e.g. `"fig.pdf"` or
#'   `c("fig.png", "fig.pdf")`; the format follows the extension, and a
#'   name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved plot in inches.
#' @return A ggplot object.
#' @references Fiorentino, D., Nunez-Riboni, I., Pierce, M.E., Oesterwind,
#'   D., Akimova, A. (2025) Improving species distribution models for
#'   climate change studies: ecological plausibility and performance
#'   metrics. *Ecological Modelling* 508, 111207.
#'   \doi{10.1016/j.ecolmodel.2025.111207}
#' @export
#' @examples
#' v <- seq(0, 30, length.out = 60)
#' curves <- data.frame(variable = "temp", value = v,
#'                      suitability = exp(-(v - 15)^2 / 40))
#' plot_response_check(response_check(curves))
plot_response_check <- function(x, variables = NULL, models = NULL,
                                base_size = 11, file = NULL, width = 10,
                                height = 6) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "response_check")) {
    stop("`x` must come from response_check().", call. = FALSE)
  }
  cur <- x$curves
  s <- x$shapes
  if (!is.null(variables)) {
    cur <- cur[cur$variable %in% variables, ]
    s <- s[s$variable %in% variables, ]
  }
  if (!is.null(models)) {
    cur <- cur[cur$model %in% models, ]
    s <- s[s$model %in% models, ]
  }
  if (!nrow(cur)) stop("Nothing to plot.", call. = FALSE)
  ink <- "#3a3a38"
  muted <- "#8a8a86"
  mods <- unique(cur$model)
  pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
           "#e6ab02", "#a6761d", "#666666")
  pal <- if (length(mods) <= 8) pal[seq_along(mods)] else
    grDevices::hcl.colors(length(mods), "Dark 3")
  names(pal) <- mods
  cur$group <- paste(cur$model, cur$run)
  ymax <- max(cur$suitability, na.rm = TRUE)
  p <- ggplot2::ggplot()
  if (!is.null(x$density)) {
    dd <- x$density[x$density$variable %in% unique(cur$variable) &
                      !is.na(x$density$density), ]
    if (nrow(dd)) {
      # scaled to the lower 30% of the panel, per variable
      dd$y <- stats::ave(dd$density, dd$variable,
                         FUN = function(v) 0.3 * ymax * v / max(v))
      p <- p + ggplot2::geom_area(
        data = dd, ggplot2::aes(.data$value, .data$y),
        fill = "#bdbdbd", alpha = 0.7)
    }
  }
  lim <- x$limits[x$limits$variable %in% unique(cur$variable), ]
  p <- p + ggplot2::geom_vline(data = lim,
                               ggplot2::aes(xintercept = .data$low),
                               linetype = "dotted", colour = muted) +
    ggplot2::geom_vline(data = lim, ggplot2::aes(xintercept = .data$high),
                        linetype = "dotted", colour = muted)
  multi_run <- length(unique(cur$run)) > 1
  p <- p + ggplot2::geom_line(
    data = cur, ggplot2::aes(.data$value, .data$suitability,
                             colour = .data$model, group = .data$group),
    linewidth = if (multi_run) 0.35 else 0.8,
    alpha = if (multi_run) 0.6 else 1)
  # shape labels: per model, the most common shape over runs
  lab <- do.call(rbind, lapply(split(s, list(s$variable, s$model),
                                     drop = TRUE), function(r) {
    tb <- sort(table(r$shape), decreasing = TRUE)
    data.frame(variable = r$variable[1], model = r$model[1],
               label = paste0(r$model[1], ": ", names(tb)[1],
                              if (length(tb) > 1)
                                paste0(" (", tb[1], "/", nrow(r), " runs)"),
                              if (any(r$rising_low | r$rising_high))
                                ", rising at edge" else ""),
               stringsAsFactors = FALSE)
  }))
  lab$rank <- stats::ave(seq_len(nrow(lab)), lab$variable, FUN = seq_along)
  xr <- stats::aggregate(value ~ variable, cur, min)
  lab$x <- xr$value[match(lab$variable, xr$variable)]
  lab$y <- ymax * (1.12 - 0.08 * (lab$rank - 1))
  p <- p + ggplot2::geom_text(
    data = lab, ggplot2::aes(.data$x, .data$y, label = .data$label,
                             colour = .data$model),
    hjust = 0, size = base_size / 3.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = pal, name = "Model") +
    ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::labs(x = "Predictor value", y = "Predicted suitability",
                  title = "Response curves",
                  subtitle = paste0(
                    "Dotted: calibrated range",
                    if (!is.null(x$density)) paste0(
                      "; grey: ", if (x$settings$ratio)
                        "presence/background density ratio"
                      else "density of records", " (scaled)"), ".")) +
    theme_enm(base_size, ink, muted)
  if (length(mods) == 1) p <- p + ggplot2::guides(colour = "none")
  p
}
