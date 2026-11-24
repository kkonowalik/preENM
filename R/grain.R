# Sensitivity of the modelling design to grain --------------------------------

#' How the modelling design changes with grain
#'
#' The grain (cell size) of an ENM is usually fixed by whatever resolution
#' the predictor rasters happen to have, yet spatial scale and resolution
#' can change the conclusions of species distribution models (Henderson et
#' al. 2023), and how well suitability predicts local abundance depends on
#' grain and on the scale of the covariates (Brambilla et al. 2023;
#' Ferrer-Ferrando et al. 2025). `grain_sensitivity()` repeats the
#' preparation steps of preENM at a series of coarser grains and reports
#' how the design changes, so that the grain can be chosen, or at least
#' reported, knowingly. It fits no model.
#'
#' **At each grain** (the native resolution times each of `factors`):
#' 1. The layers are aggregated on disk with GDAL, file to file, one grain
#'    at a time (as in [crop_rasters()]): every coarse cell gets the mean
#'    (`method = "average"`) of the native cells it covers, written as
#'    32-bit floating point so that averages of integer layers are not
#'    rounded. Each grain's files are deleted before the next is made,
#'    unless `out_dir` is given.
#' 2. [thin_by_cell()]: one record per cell (and species).
#' 3. [thin_environment()] on those records (`env_thin = TRUE`). As the
#'    records now sit at cell centres, their positional accuracy is taken
#'    as one cell of that grain (`radius_km`), unless `env_thin_args` says
#'    otherwise.
#' 4. [uncertainty_environment()] on all records, if they have a coordinate
#'    uncertainty: the share of records whose uncertainty circle spans too
#'    much of the niche. The share can go either way: a coarse cell
#'    contains the circle more often, but a circle that crosses the edge of
#'    a coarse cell spans two cells whose values differ more.
#' 5. [sample_environment()] and [select_variables()] on every data set in
#'    `env`: the passing variable sets, and the *selected* set, i.e. the
#'    top-ranked set of `size` variables. With a single data set the only
#'    criterion is the VIF; give the periods or scenarios the model will be
#'    projected to as further data sets so that novelty and gradient
#'    stability count too, as in your own selection.
#'
#' **Agreement.** The selected set at each grain is compared with the one at
#' the finest grain by the Jaccard index (shared variables / all variables
#' of the two sets; 1 = the same set). `set_agreement` is the mean Jaccard
#' over all pairs of grains (1 when the same set is selected at every
#' grain; a grain without a passing set counts as 0 against one with a
#' set). The collections of passing sets are compared in the same way
#' (`passing_jaccard`), which is less sensitive to near ties in the
#' ranking. The correlation structure is summarised as the number of pairs
#' with |Spearman rho| of at least `thr_r` in any data set, and the mean
#' absolute change of rho from the finest grain.
#'
#' **Stable grain.** The coarsest grain up to which the selected set agrees
#' with the finest grain at every grain, with a Jaccard of at least
#' `stable` (default 1, the same set). This is the practical answer: the
#' grain can be coarsened that far without changing the predictors.
#'
#' Memory use stays low: only the values at the records and random samples
#' of cells are read (via [terra::extract()] and the samplers of preENM).
#'
#' @param records Presence records: a data.frame with coordinates, a CSV
#'   path, an `sf` object or a [terra::SpatVector]. They may carry a
#'   coordinate-uncertainty column.
#' @param env The predictor layers at their native grain: a
#'   [terra::SpatRaster] (or file paths), or a named list of them, one per
#'   data set (e.g. present, LGM, 2070), as for [sample_environment()]. The
#'   first data set is used for the thinning and the uncertainty.
#' @param factors Aggregation factors (whole numbers, 1 = native grain).
#' @param method GDAL resampling method for aggregation: `"average"`
#'   (default), `"med"` (median), `"mode"` (for categorical layers) or
#'   `"near"`.
#' @param sizes,size Set sizes searched by [select_variables()], and the
#'   size of the selected set compared between grains (default: the largest
#'   size with a passing set at the finest grain).
#' @param thr_r,thr_vif Collinearity limits of [select_variables()].
#' @param n_sample Cells drawn per data set by [sample_environment()].
#' @param env_thin Run [thin_environment()] at every grain.
#' @param uncertainty Coordinate uncertainty for [uncertainty_environment()]:
#'   a column name (default: detected), one value per record or a single
#'   value; `NA` to skip this step.
#' @param max_rel,max_radius_km Passed to [uncertainty_environment()].
#' @param stable Lowest Jaccard index with the finest grain that still
#'   counts as the same selection.
#' @param thin_args,env_thin_args,sample_args,select_args Further arguments
#'   (named lists) for [thin_by_cell()], [thin_environment()],
#'   [sample_environment()] and [select_variables()].
#' @param out_dir Optional folder to keep the aggregated layers (one
#'   sub-folder per grain); by default they are temporary.
#' @param coords,points_crs Coordinate columns and their CRS when `records`
#'   is a table.
#' @param seed Optional random seed, used in the same way at every grain so
#'   that differences come from the grain, not from the sampling.
#' @param verbose Print progress.
#' @return An object of class `grain_sensitivity`: `table` (one row per
#'   grain: `factor`, `res_x`, `res_y`, `cell_km` (side of a square cell of
#'   the same area at the records' median latitude), `records_cell`,
#'   `records_env`, `unc_examined`, `unc_lost`, `unc_lost_share`,
#'   `n_passing`, `selected`, `tied` (passing sets sharing the top rank),
#'   `jaccard`, `passing_jaccard`,
#'   `collinear_pairs`, `rho_change`, `stable`), `sets` (the passing sets of
#'   each grain), `rho` (the Spearman matrices of the first data set),
#'   `set_agreement`, `stable_factor`, `notes` (steps that failed at a
#'   grain) and `settings`. See
#'   [plot_grain_sensitivity()] and [write_grain_sensitivity()].
#' @references Henderson, A.F., Santoro, J.A., Kremer, P. (2023) Impacts of
#'   spatial scale and resolution on species distribution models of American
#'   chestnut (*Castanea dentata*) in Pennsylvania, USA. *Forest Ecology and
#'   Management* 529, 120741. \doi{10.1016/j.foreco.2022.120741}
#'
#'   Brambilla, M., Bazzi, G., Ilahiane, L. (2023) The effectiveness of
#'   species distribution models in predicting local abundance depends on
#'   model grain size. *Ecology* 105, e4224. \doi{10.1002/ecy.4224}
#'
#'   Ferrer-Ferrando, D., Tarroso, P., Telleria, J.L., Acevedo, P.,
#'   Fernandez-Lopez, J. (2025) Disentangling the effect of the spatial
#'   scale and species spatial pattern on the abundance-suitability
#'   relationship. *Ecography* 2025, e07766. \doi{10.1111/ecog.07766}
#' @seealso [thin_by_cell()], [thin_environment()],
#'   [uncertainty_environment()], [select_variables()]
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_enm_data(seed = 1, n_presence = 200)
#' e <- sim$env
#' env <- c(e, e[[1]] * e[[2]], sqrt(e[[2]]))
#' names(env) <- c("t", "p", "txp", "sp")
#' g <- grain_sensitivity(sim$presence, env, factors = c(1, 2, 4),
#'                        sizes = 2:3, n_sample = 2000, env_thin = FALSE,
#'                        uncertainty = 5000, seed = 1, verbose = FALSE)
#' g
#' plot_grain_sensitivity(g)
#' }
grain_sensitivity <- function(records, env, factors = c(1, 2, 4, 8),
                              method = c("average", "med", "mode", "near"),
                              sizes = 3:6, size = NULL, thr_r = 0.8,
                              thr_vif = 5, n_sample = 1e4, env_thin = TRUE,
                              uncertainty = NULL, max_rel = 0.1,
                              max_radius_km = 20, stable = 1,
                              thin_args = list(), env_thin_args = list(),
                              sample_args = list(), select_args = list(),
                              out_dir = NULL, coords = NULL,
                              points_crs = "EPSG:4326", seed = NULL,
                              verbose = TRUE) {
  method <- match.arg(method)
  say <- function(...) if (verbose) message(...)
  factors <- sort(unique(as.integer(round(factors))))
  if (!length(factors) || any(factors < 1)) {
    stop("`factors` must be whole numbers of at least 1.", call. = FALSE)
  }
  if (!is.list(env) || inherits(env, "SpatRaster")) env <- list(env)
  if (is.null(names(env)) || any(!nzchar(names(env)))) {
    names(env) <- if (length(env) == 1) "present" else
      paste0("set", seq_along(env))
  }
  env <- lapply(env, as_grid)
  nV <- terra::nlyr(env[[1]])
  sizes <- sizes[sizes <= nV]
  if (!length(sizes)) stop("Every set size exceeds the number of layers.",
                           call. = FALSE)
  if (is.character(records) && length(records) == 1) {
    records <- read_any_csv(records)
  }
  if (!is.data.frame(records) || inherits(records, "sf")) {
    p <- as_points(records, coords, points_crs)
    records <- data.frame(as.data.frame(p),
                          lon = terra::crds(p)[, 1],
                          lat = terra::crds(p)[, 2])
    points_crs <- terra::crs(p)
    coords <- c("lon", "lat")
  }
  if (is.null(coords)) coords <- guess_coords(names(records))
  if (identical(uncertainty, NA)) {
    do_unc <- FALSE
  } else if (is.null(uncertainty)) {
    uncertainty <- find_column(names(records), c(
      "coordinateUncertaintyInMeters", "uncertainty", "unc_m"))
    do_unc <- !is.na(uncertainty)
    if (!do_unc) say("No coordinate-uncertainty column: that step is ",
                     "skipped.")
  } else {
    do_unc <- TRUE
  }

  # median latitude of the records, for cell sizes in km
  r0 <- env[[1]]
  xy <- cbind(as_number(records[[coords[1]]]), as_number(records[[coords[2]]]))
  xy <- xy[stats::complete.cases(xy), , drop = FALSE]
  lat0 <- if (nrow(xy)) stats::median(terra::crds(terra::project(
    terra::vect(xy, crs = points_crs), "EPSG:4326"))[, 2]) else 0
  cell_km <- function(r) {
    rs <- terra::res(r)
    if (terra::is.lonlat(r)) {
      sqrt(rs[1] * 111.32 * cos(lat0 * pi / 180) * rs[2] * 111.32)
    } else {
      sqrt(prod(rs)) / 1000
    }
  }

  tmp <- tempfile("preENM_grain")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  rows <- list()
  sets <- list()
  rho <- list()
  unions <- list()
  notes <- list()
  for (f in factors) {
    say(sprintf("Grain x%d ...", f))
    gdir <- if (is.null(out_dir)) file.path(tmp, paste0("x", f)) else
      file.path(out_dir, paste0("x", f))
    eg <- if (f == 1) env else lapply(names(env), function(k) {
      aggregate_on_disk(env[[k]], f, method, file.path(gdir, k))
    })
    names(eg) <- names(env)
    ref <- eg[[1]]

    # cell thinning, environmental thinning
    ct <- quiet(do.call(thin_by_cell, c(list(x = records, raster = ref,
                                             coords = coords,
                                             points_crs = points_crs),
                                        thin_args)))
    n_env <- NA_integer_
    if (env_thin && nrow(ct) >= 3) {
      # records now sit at cell centres: accurate to about one cell
      ea <- env_thin_args
      if (is.null(ea$radius_km) && is.null(ea$window)) {
        ea$radius_km <- cell_km(ref)
      }
      th <- tryCatch(quiet(do.call(thin_environment, c(list(
        points = ct, env = ref, coords = c("lon", "lat"), seed = seed,
        verbose = FALSE), ea))), error = function(e) {
          notes[[length(notes) + 1]] <<- paste0(
            "x", f, ": thin_environment() failed (", conditionMessage(e),
            ")")
          NULL
        })
      if (!is.null(th)) n_env <- sum(th$data$kept)
    }
    # records lost to coordinate uncertainty
    u_exam <- u_lost <- NA_integer_
    if (do_unc) {
      u <- quiet(uncertainty_environment(records, ref,
                                         uncertainty = uncertainty,
                                         max_rel = max_rel,
                                         max_radius_km = max_radius_km,
                                         coords = coords,
                                         points_crs = points_crs,
                                         verbose = FALSE))
      u_exam <- sum(!is.na(u$unc_ok))
      u_lost <- sum(!u$unc_ok, na.rm = TRUE)
    }
    # variable selection
    es <- quiet(do.call(sample_environment, c(list(
      datasets = eg, n = n_sample, seed = seed, verbose = FALSE),
      sample_args)))
    vs <- quiet(do.call(select_variables, c(list(
      x = es, sizes = sizes, thr_r = thr_r, thr_vif = thr_vif),
      select_args)))
    if (is.null(vs$sets)) {
      notes[[length(notes) + 1]] <- paste0("x", f, ": no variable set ",
                                           "passes the screen")
    }
    sets[[as.character(f)]] <- if (is.null(vs$sets)) data.frame(
      set = character(0), n_vars = integer(0), rank = integer(0)) else
        vs$sets[c("set", "n_vars", "rank", "mean_rank")]
    rho[[as.character(f)]] <- vs$rho[[1]]
    unions[[as.character(f)]] <- vs$union
    rows[[length(rows) + 1]] <- data.frame(
      factor = f, res_x = terra::res(ref)[1], res_y = terra::res(ref)[2],
      cell_km = cell_km(ref), records_cell = nrow(ct),
      records_env = n_env, unc_examined = u_exam, unc_lost = u_lost,
      unc_lost_share = if (is.na(u_exam) || u_exam == 0) NA_real_ else
        u_lost / u_exam)
    if (f > 1 && is.null(out_dir)) unlink(gdir, recursive = TRUE)
  }
  tab <- do.call(rbind, rows)

  # the selected set at each grain and its agreement with the finest
  k1 <- as.character(factors[1])
  if (is.null(size)) {
    have <- unique(sets[[k1]]$n_vars)
    size <- if (length(have)) max(have) else max(sizes)
  }
  pick <- function(s) {
    s <- s[s$n_vars == size, , drop = FALSE]
    if (!nrow(s)) NA_character_ else s$set[order(s$rank, s$mean_rank)][1]
  }
  sel <- vapply(sets, pick, character(1))
  # sets sharing the top rank: the selected one is then the first of them
  tied <- vapply(sets, function(s) {
    s <- s[s$n_vars == size, , drop = FALSE]
    if (!nrow(s)) NA_integer_ else sum(s$mean_rank == min(s$mean_rank))
  }, integer(1))
  passing <- lapply(sets, function(s) s$set[s$n_vars == size])
  tab$n_passing <- vapply(passing, length, integer(1))
  tab$selected <- unname(sel)
  tab$tied <- unname(tied)
  tab$jaccard <- vapply(sel, jaccard_sets, numeric(1), b = sel[[k1]])
  tab$passing_jaccard <- vapply(passing, function(p) {
    a <- passing[[k1]]
    if (!length(a) && !length(p)) NA_real_ else
      length(intersect(a, p)) / length(union(a, p))
  }, numeric(1))
  lt <- lower.tri(rho[[k1]])
  tab$collinear_pairs <- vapply(unions, function(A) sum(A[lt] >= thr_r),
                                integer(1))
  tab$rho_change <- vapply(rho, function(R) mean(abs(R[lt] - rho[[k1]][lt])),
                           numeric(1))
  tab$stable <- cumprod(!is.na(tab$jaccard) &
                          tab$jaccard >= stable - 1e-12) == 1
  pairs <- utils::combn(seq_along(sel), 2)
  # a grain without a passing set disagrees with one that has one
  agreement <- if (length(sel) < 2) NA_real_ else
    mean(apply(pairs, 2, function(j) {
      a <- sel[j[1]]
      b <- sel[j[2]]
      if (is.na(a) && is.na(b)) NA_real_ else if (is.na(a) || is.na(b)) 0
      else jaccard_sets(a, b)
    }), na.rm = TRUE)
  stable_factor <- if (any(tab$stable)) max(tab$factor[tab$stable]) else
    NA_integer_
  rownames(tab) <- NULL

  out <- list(table = tab, sets = sets, rho = rho,
              set_agreement = agreement, stable_factor = stable_factor,
              notes = unlist(notes),
              settings = list(factors = factors, method = method,
                              size = size, sizes = sizes, thr_r = thr_r,
                              thr_vif = thr_vif, stable = stable,
                              max_rel = max_rel, env_thin = env_thin,
                              uncertainty = do_unc, data_sets = names(env),
                              native_res = terra::res(r0),
                              lonlat = terra::is.lonlat(r0),
                              out_dir = out_dir))
  class(out) <- "grain_sensitivity"
  out
}

# Jaccard index of two variable sets written as space-separated names
# (internal); NA when either is missing
jaccard_sets <- function(a, b) {
  if (is.na(a) || is.na(b)) return(NA_real_)
  a <- strsplit(a, " ", fixed = TRUE)[[1]]
  b <- strsplit(b, " ", fixed = TRUE)[[1]]
  length(intersect(a, b)) / length(union(a, b))
}

# Aggregate a raster stack by `fact` on disk with GDAL, one warp per source
# file (internal). Layers held in memory are written to a file first. The
# coarse grid starts at the top-left corner of the input, as in
# terra::aggregate(), and the output is Float32 so averages are not rounded.
aggregate_on_disk <- function(r, fact, method, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  src <- terra::sources(r, bands = TRUE)
  e <- terra::ext(r)
  rs <- terra::res(r) * fact
  nc <- ceiling(terra::ncol(r) / fact)
  nr <- ceiling(terra::nrow(r) / fact)
  te <- c(terra::xmin(e), terra::ymax(e) - nr * rs[2],
          terra::xmin(e) + nc * rs[1], terra::ymax(e))
  layers <- vector("list", terra::nlyr(r))
  for (s in unique(src$sid)) {
    i <- which(src$sid == s)
    f <- src$source[i[1]]
    band <- src$bands[i]
    if (!nzchar(f) || !file.exists(f)) {
      # in memory (or not a plain file): write these layers first
      f <- file.path(dir, sprintf("in_%d.tif", s))
      terra::writeRaster(r[[i]], f, overwrite = TRUE)
      band <- seq_along(i)
    }
    out <- file.path(dir, sprintf("agg_%d.tif", s))
    opts <- c("-te", num_arg(te), "-tr", num_arg(rs), "-r", method,
              "-ot", "Float32", "-overwrite", "-co", "COMPRESS=LZW",
              "-co", "TILED=YES", "-co", "BIGTIFF=IF_SAFER")
    sf::gdal_utils("warp", source = f, destination = out, options = opts,
                   quiet = TRUE)
    if (!file.exists(out)) stop("GDAL did not write ", out, call. = FALSE)
    a <- terra::rast(out)
    for (j in seq_along(i)) layers[[i[j]]] <- a[[band[j]]]
  }
  out <- terra::rast(layers)
  names(out) <- names(r)
  out
}

#' @export
print.grain_sensitivity <- function(x, digits = 2, ...) {
  t <- x$table
  s <- x$settings
  say <- function(txt) cat(strwrap(txt, width = 76, initial = "- ",
                                   exdent = 2), sep = "\n")
  km <- function(i) {
    v <- t$cell_km[i]
    paste0("about ", if (v < 10) signif(v, 2) else round(v), " km")
  }
  g <- function(i) paste0("x", t$factor[i], " (", km(i), ")")
  n <- nrow(t)
  cat("Grain sensitivity: ", n, " grain(s), x", t$factor[1], " to x",
      t$factor[n], " the native ", paste(signif(s$native_res, 3),
                                         collapse = " x "),
      if (s$lonlat) " degrees" else " map units", "; ",
      length(s$data_sets), " data set(s)\n\nIN SHORT\n", sep = "")
  if (n < 2) {
    say("Only one grain was examined: nothing to compare. Give several `factors`.")
  }
  sel1 <- t$selected[1]
  if (is.na(sel1)) {
    say(paste0("No set of ", s$size, " variables passes the screen at the ",
               "finest grain, so the stability of the selection cannot be ",
               "judged. Lower `size` or relax `thr_r` / `thr_vif`."))
  } else if (n >= 2) {
    k <- max(which(t$stable))
    if (k == n) {
      say(paste0("The selected set of ", s$size, " variables (", sel1,
                 ") is the same at every grain tested, up to ", g(n),
                 ", the coarsest tested; coarser grains were not ",
                 "examined."))
    } else {
      chg <- if (is.na(t$selected[k + 1])) "no set of this size passes" else
        paste0("the selected set changes to ", t$selected[k + 1],
               " (Jaccard ", round(t$jaccard[k + 1], 2), ")")
      say(paste0(if (k == 1) paste0("The selected set of ", s$size,
                                    " variables (", sel1, ") already changes",
                                    " at the first coarser grain: at ",
                                    g(2), ", ", chg, ". The selection ",
                                    "depends on grain; choose the grain ",
                                    "first, on ecological grounds.")
                 else paste0("The selected set of ", s$size, " variables (",
                             sel1, ") stays stable up to ", g(k),
                             "; at ", g(k + 1), ", ", chg, ".")))
    }
    say(paste0("Over all pairs of grains, the selected sets agree with a ",
               "mean Jaccard index of ", round(x$set_agreement, 2),
               if (any(!is.na(t$passing_jaccard[-1]))) paste0(
                 "; the collections of passing sets with the finest ",
                 "grain, ", paste(round(t$passing_jaccard[-1], 2),
                                  collapse = ", ")), "."))
  }
  if (!is.na(sel1) && any(t$tied > 1, na.rm = TRUE)) {
    w <- which(t$tied > 1)
    say(paste0("At ", paste0("x", t$factor[w], collapse = ", "), ", ",
               if (length(unique(t$tied[w])) == 1) t$tied[w][1] else
                 "several", " passing sets share the top rank, so the ",
               "selected set is only the first of equals there and its ",
               "Jaccard index may change by chance; the agreement of the ",
               "passing sets is the safer measure. Giving the data sets of ",
               "the projection periods adds criteria that separate the ",
               "sets."))
  }
  if (n >= 2) {
    say(paste0("Cell thinning keeps ", t$records_cell[1], " records at ",
               g(1), " and ", t$records_cell[n], ", ",
               share_label(t$records_cell[n] / max(t$records_cell[1], 1)),
               " of them, at ", g(n), if (!all(is.na(t$records_env))) paste0(
                 "; environmental thinning keeps ", t$records_env[1],
                 " and ", t$records_env[n]), "."))
    if (s$uncertainty && !all(is.na(t$unc_lost_share))) {
      say(paste0("Records too uncertain for the grain: ",
                 share_label(t$unc_lost_share[1]), " at x", t$factor[1],
                 ", ", share_label(t$unc_lost_share[n]), " at x",
                 t$factor[n], "."))
    }
    if (any(t$collinear_pairs != t$collinear_pairs[1])) {
      say(paste0("The number of collinear pairs changes from ",
                 t$collinear_pairs[1], " to ", t$collinear_pairs[n],
                 ": aggregation changes the correlation structure."))
    }
  }
  for (nt in x$notes) say(paste0("Note: ", nt, "."))
  say(paste0("Grain is a choice about the ecology (the scale at which the ",
             "species responds) and the data (the positional accuracy of the ",
             "records), not only about the variables: use the stable range ",
             "as the room within which to choose, and report the grain."))
  cat("\n")
  show <- t[c("factor", "cell_km", "records_cell", "records_env",
              "unc_lost_share", "n_passing", "selected", "tied", "jaccard",
              "passing_jaccard", "collinear_pairs", "rho_change")]
  if (all(is.na(show$records_env))) show$records_env <- NULL
  if (all(is.na(show$unc_lost_share))) show$unc_lost_share <- NULL
  num <- vapply(show, is.double, logical(1))
  show[num] <- lapply(show[num], round, digits)
  print(show, row.names = FALSE)
  cat("", strwrap(paste0(
    "factor = aggregation factor of the native cells (", s$method,
    "); cell_km = side of a square cell of the same area at the records' ",
    "median latitude; records_cell / records_env = records kept by ",
    "thin_by_cell() / thin_environment(); unc_lost_share = share of the ",
    "examined records whose uncertainty circle spans more than ", s$max_rel,
    " of the niche breadth (uncertainty_environment()); selected = ",
    "top-ranked passing set of ", s$size, " variables (select_variables(), ",
    "|rho| < ", s$thr_r, ", VIF <= ", s$thr_vif, "); tied = passing sets ",
    "sharing the top rank; jaccard = shared / all ",
    "variables of this and the finest grain's set; passing_jaccard = the ",
    "same for the collections of passing sets; collinear_pairs = pairs with ",
    "|rho| >= ", s$thr_r, " in any data set; rho_change = mean |change| of ",
    "Spearman rho from the finest grain. Stable = Jaccard >= ", s$stable,
    " at this and every finer grain; this limit and the reading of the ",
    "Jaccard index (1 the same set, 0.5 or less mostly different) are ",
    "rules of thumb, not values from the literature."), width = 76),
    sep = "\n")
  invisible(x)
}

#' Plot how the modelling design changes with grain
#'
#' Four panels against cell size: records kept after cell and environmental
#' thinning, the share of records too uncertain for the grain, the
#' agreement of the selected and passing variable sets with the finest
#' grain (Jaccard index), and the number of collinear variable pairs.
#'
#' @param x A `grain_sensitivity` from [grain_sensitivity()].
#' @param base_size Base font size.
#' @param file Optional file name(s) to save the figure; the format follows
#'   the extension, a name without extension is saved as PNG and PDF.
#' @param width,height Size of the saved figure in inches.
#' @return A ggplot object (invisibly the saved file names when `file` is
#'   given).
#' @references Henderson, A.F., Santoro, J.A., Kremer, P. (2023) Impacts of
#'   spatial scale and resolution on species distribution models of American
#'   chestnut (*Castanea dentata*) in Pennsylvania, USA. *Forest Ecology and
#'   Management* 529, 120741. \doi{10.1016/j.foreco.2022.120741}
#' @export
#' @examples
#' \donttest{
#' # see ?grain_sensitivity
#' }
plot_grain_sensitivity <- function(x, base_size = 11, file = NULL,
                                   width = 7, height = 5.5) {
  if (!is.null(file)) {
    return(save_call(match.call(), parent.frame(), file, width, height))
  }
  if (!inherits(x, "grain_sensitivity")) {
    stop("`x` must come from grain_sensitivity().", call. = FALSE)
  }
  t <- x$table
  panels <- c("Records kept", "Records too uncertain (%)",
              "Agreement with the finest grain (Jaccard)",
              "Collinear variable pairs")
  long <- function(v, panel, series) {
    data.frame(cell_km = t$cell_km, value = v,
               panel = panel, series = series, stringsAsFactors = FALSE)
  }
  d <- rbind(long(t$records_cell, panels[1], "After cell thinning"),
             long(t$records_env, panels[1], "After environmental thinning"),
             long(100 * t$unc_lost_share, panels[2], "Too uncertain"),
             long(t$jaccard, panels[3], "Selected set"),
             long(t$passing_jaccard, panels[3], "Passing sets"),
             long(t$collinear_pairs, panels[4], "Collinear pairs"))
  d <- d[!is.na(d$value), , drop = FALSE]
  d$panel <- factor(d$panel, panels)
  lv <- c("After cell thinning", "After environmental thinning",
          "Too uncertain", "Selected set", "Passing sets", "Collinear pairs")
  d$series <- factor(d$series, intersect(lv, d$series))
  # ColorBrewer Dark2
  pal <- stats::setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                           "#66a61e", "#a6761d"), lv)
  lab <- paste0("x", t$factor, "\n", ifelse(t$cell_km < 10,
                                            signif(t$cell_km, 2),
                                            round(t$cell_km)))
  shp <- stats::setNames(c(16, 17, 16, 16, 1, 16), lv)
  ltp <- stats::setNames(c("solid", "solid", "solid", "solid", "dashed",
                           "solid"), lv)
  # the Jaccard panel always spans 0 to 1
  lim <- data.frame(cell_km = t$cell_km[1], value = c(0, 1),
                    panel = factor(panels[3], panels))
  # whole-number breaks for counts and percentages; the Jaccard panel
  # (0 to 1) keeps the decimal ones
  brk <- function(l) {
    b <- pretty(l)
    if (max(l) <= 1.05) b else b[b == round(b)]
  }
  sv <- levels(d$series)
  ggplot2::ggplot(d, ggplot2::aes(.data$cell_km, .data$value,
                                  colour = .data$series)) +
    ggplot2::geom_blank(data = lim, inherit.aes = FALSE,
                        ggplot2::aes(.data$cell_km, .data$value)) +
    ggplot2::geom_line(ggplot2::aes(linetype = .data$series),
                       linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$series), size = 2.2,
                        stroke = 0.9) +
    ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_y") +
    ggplot2::scale_x_continuous(trans = "log2", breaks = t$cell_km,
                                labels = lab) +
    ggplot2::scale_y_continuous(breaks = brk) +
    ggplot2::scale_colour_manual(values = pal[sv], name = NULL) +
    ggplot2::scale_shape_manual(values = shp[sv], name = NULL) +
    ggplot2::scale_linetype_manual(values = ltp[sv], name = NULL) +
    ggplot2::labs(x = "Aggregation factor and cell size (km)", y = NULL,
                  title = "Modelling design at each grain",
                  subtitle = paste0("Selected set: top-ranked passing set ",
                                    "of ", x$settings$size, " variables")) +
    theme_enm(base_size, "#3a3a38", "#8a8a86") +
    ggplot2::theme(legend.position = "bottom")
}

#' Write the grain-sensitivity results to CSV files
#'
#' @param x A `grain_sensitivity` from [grain_sensitivity()].
#' @param dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @return Invisibly, the paths of the two files written: the per-grain
#'   table and the passing variable sets of every grain.
#' @references Henderson, A.F., Santoro, J.A., Kremer, P. (2023) Impacts of
#'   spatial scale and resolution on species distribution models of American
#'   chestnut (*Castanea dentata*) in Pennsylvania, USA. *Forest Ecology and
#'   Management* 529, 120741. \doi{10.1016/j.foreco.2022.120741}
#' @export
#' @examples
#' \donttest{
#' # see ?grain_sensitivity
#' }
write_grain_sensitivity <- function(x, dir, prefix = "grain_sensitivity") {
  if (!inherits(x, "grain_sensitivity")) {
    stop("`x` must come from grain_sensitivity().", call. = FALSE)
  }
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f1 <- file.path(dir, paste0(prefix, "_table.csv"))
  f2 <- file.path(dir, paste0(prefix, "_sets.csv"))
  utils::write.csv(x$table, f1, row.names = FALSE)
  s <- do.call(rbind, lapply(names(x$sets), function(k) {
    d <- x$sets[[k]]
    if (!nrow(d)) return(NULL)
    data.frame(factor = as.integer(k), d, row.names = NULL)
  }))
  if (is.null(s)) s <- data.frame(factor = integer(0), set = character(0))
  utils::write.csv(s, f2, row.names = FALSE)
  invisible(c(f1, f2))
}
