#' Search all variable sets that are uncorrelated in every period and region
#'
#' A systematic, transparent alternative to picking predictors by hand from a
#' correlation dendrogram. Every subset of the candidate variables is screened
#' for collinearity in *all* data sets (e.g. the present and the past or
#' future periods the model will be projected to) and, more leniently, inside
#' regions; the sets that pass are scored on transferability, and the best
#' trade-offs are reported, ranked, for you to choose from on ecological
#' grounds.
#'
#' **Screen.** A pair of variables is *compatible* when its |Spearman rho| is
#' below `thr_r` in every data set (own cells and, for the non-reference data
#' sets, the common cells) and below `thr_r_regional` in every region. A set
#' passes when all its pairs are compatible, its rank-based VIF is at most
#' `thr_vif` in every data set and `thr_vif_regional` in every region, and it
#' contains at least one variable of each group in `require_one_of` (e.g. a
#' temperature and a moisture variable). Because the dendrogram ("union tree"
#' of the largest |rho| over the data sets, complete linkage) puts every pair
#' with |rho| >= `thr_r` in one cluster, a passing set never holds two
#' variables of one cluster: choosing "one variable per cluster" is contained
#' in this search, but the pair check that cutting the tree alone does not
#' guarantee is enforced.
#'
#' **Scores** (for each passing set):
#' * `novelty`: the share of cells of the target data sets (default: all but
#'   the reference) outside the reference conditions, by the multivariate
#'   environmental similarity surface (MESS < 0; Elith et al. 2010); the
#'   worst target counts, the others are in `novelty_<target>`. Lower is
#'   better: less extrapolation.
#' * `gradient_min`: the weakest rank correlation, between the paired data
#'   sets on the same cells, of any member: does each variable keep its
#'   spatial pattern through time? Higher is better.
#' * `agreement_min`: with `agreement = c(a, b)` (two data sources for the
#'   same period, e.g. CHELSA and TraCE21k present), the weakest agreement of
#'   any member. Higher is better.
#' * `worst_vif`, `regional_vif`: the largest rank-based VIF over the data
#'   sets and over the regions. Lower is better.
#'
#' There is no hand-made total score: within each set size the Pareto front
#' (sets no other set beats on all criteria) is marked and a transparent mean
#' rank over the criteria orders the sets.
#'
#' @param x An `env_sample` from [sample_environment()].
#' @param sizes Set sizes to search.
#' @param thr_r,thr_vif Collinearity limits in the data sets.
#' @param regional Use the regions (if the sample has them) in the screen.
#' @param thr_r_regional,thr_vif_regional Softer limits inside regions.
#' @param min_region_n Regions with fewer sampled cells are ignored.
#' @param screen_sets Data sets used in the screen (default: all).
#' @param require_one_of Optional named list of variable groups; every set
#'   must contain at least one variable of each group.
#' @param targets Data sets whose novelty is scored (default: all but the
#'   reference).
#' @param gradient_pairs List of data-set pairs compared on common cells for
#'   gradient stability (default: the reference with each target).
#' @param agreement Optional pair of data sets of the same period from
#'   different sources.
#' @param criteria Named vector of the scores used for the Pareto front and
#'   the mean rank, with `"min"` or `"max"`. Default: `novelty`,
#'   `gradient_min`, `agreement_min`, `worst_vif`, `regional_vif` (scores
#'   that cannot be computed are left out).
#' @param weights Optional named weights of the criteria in the mean rank.
#' @param reference_sets Optional named list of sets (e.g. your earlier
#'   choice) scored for comparison even if they fail the screen.
#' @param n_mess Target cells used for the novelty scores.
#' @param use_common Also screen the non-reference data sets on the common
#'   cells.
#' @param max_subsets Safety limit on the number of subsets of one size.
#' @return An object of class `var_selection`; see [best_sets()],
#'   [plot_variable_tree()], [plot_correlation()], [plot_variable_sets()] and
#'   [write_variable_selection()]. Its elements: `sets` (all passing sets with
#'   their scores, `pareto`, `mean_rank` and `rank` within the size),
#'   `counts` (sets surviving each filter), `pairs` (each pair: rho per data
#'   set, largest and smallest |rho|, the data set driving it, sign changes,
#'   regional maximum, `compatible`), `variables` (single-variable
#'   properties and cluster), `frequency` (how often each variable appears in
#'   passing, Pareto and top sets), `reference_sets`, `consistency` (the two
#'   `agreement` data sets compared per variable), `structure_agreement`
#'   (similarity of the correlation structure between data sets), `tree`,
#'   `rho` and `regional_rho`.
#' @references Elith, J., Kearney, M., Phillips, S. (2010) The art of
#'   modelling range-shifting species. *Methods in Ecology and Evolution* 1,
#'   330-342. \doi{10.1111/j.2041-210X.2010.00036.x}
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
#' best_sets(vs)
select_variables <- function(x, sizes = 3:6, thr_r = 0.8, thr_vif = 5,
                             regional = TRUE, thr_r_regional = 0.85,
                             thr_vif_regional = 10, min_region_n = 1000,
                             screen_sets = NULL, require_one_of = NULL,
                             targets = NULL, gradient_pairs = NULL,
                             agreement = NULL, criteria = NULL,
                             weights = NULL, reference_sets = NULL,
                             n_mess = 1e5, use_common = TRUE,
                             max_subsets = 2e6) {
  if (!inherits(x, "env_sample")) {
    stop("`x` must come from sample_environment().", call. = FALSE)
  }
  V <- x$vars
  nV <- length(V)
  ref <- x$reference
  if (is.null(screen_sets)) screen_sets <- names(x$sets)
  bad <- setdiff(c(screen_sets, targets, unlist(gradient_pairs), agreement),
                 names(x$sets))
  if (length(bad)) {
    stop("Unknown data set(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  if (max(sizes) > nV) stop("A set size exceeds the number of variables.",
                            call. = FALSE)
  for (g in require_one_of) {
    if (!all(g %in% V)) stop("`require_one_of` contains unknown variables.",
                             call. = FALSE)
  }
  spear <- function(d) {
    R <- suppressWarnings(stats::cor(d[V], method = "spearman"))
    dimnames(R) <- list(V, V)
    R
  }

  # correlation matrices: data sets (own cells, + common cells) and regions
  rho <- lapply(x$sets[screen_sets], spear)
  if (use_common && !is.null(x$common)) {
    extra <- setdiff(intersect(screen_sets, names(x$common)),
                     c(ref, x$anomalies))
    for (k in extra) rho[[paste0(k, "@common")]] <- spear(x$common[[k]])
  }
  rho <- rho[vapply(rho, function(R) !anyNA(R), logical(1))]
  if (!length(rho)) stop("Correlations could not be computed (constant ",
                         "variables?).", call. = FALSE)
  reg_rho <- list()
  if (regional && x$regions) {
    for (k in screen_sets) {
      d <- x$sets[[k]]
      for (rg in sort(unique(stats::na.omit(d$.region)))) {
        i <- which(d$.region == rg)
        if (length(i) >= min_region_n) {
          R <- spear(d[i, , drop = FALSE])
          if (!anyNA(R)) reg_rho[[paste(k, rg, sep = " | ")]] <- R
        }
      }
    }
    if (!length(reg_rho)) {
      message("No region has ", min_region_n, " sampled cells: the regional ",
              "screen is skipped.")
    }
  }
  absmax <- function(Rl) {
    m <- Reduce(pmax, lapply(Rl, abs))
    dimnames(m) <- list(V, V)
    m
  }
  A <- absmax(rho)
  A_reg <- if (length(reg_rho)) absmax(reg_rho) else NULL

  # union tree: one cluster = variables collinear in at least one data set
  tree <- stats::hclust(stats::as.dist(1 - A), method = "complete")
  cl <- stats::cutree(tree, h = 1 - thr_r)

  # pairs
  pr <- t(utils::combn(V, 2))
  pairs <- data.frame(var1 = pr[, 1], var2 = pr[, 2], stringsAsFactors = FALSE)
  rr <- vapply(rho, function(R) R[pr], numeric(nrow(pr)))
  rr <- matrix(rr, nrow = nrow(pr), dimnames = list(NULL, names(rho)))
  for (k in names(rho)) pairs[[paste0("rho_", k)]] <- round(rr[, k], 3)
  pairs$max_abs <- apply(abs(rr), 1, max)
  pairs$min_abs <- apply(abs(rr), 1, min)
  pairs$driver <- names(rho)[apply(abs(rr), 1, which.max)]
  pairs$sign_change <- apply(rr, 1, function(z) any(z > 0.3) && any(z < -0.3))
  pairs$collinear_only_in_some <- pairs$max_abs >= thr_r &
    pairs$min_abs < thr_r
  pairs$max_abs_regional <- if (is.null(A_reg)) NA_real_ else A_reg[pr]
  pairs$compatible <- pairs$max_abs < thr_r &
    (is.null(A_reg) | is.na(pairs$max_abs_regional) |
       pairs$max_abs_regional < thr_r_regional)
  pairs <- pairs[order(-pairs$max_abs), ]
  rownames(pairs) <- NULL
  Cmp <- matrix(TRUE, nV, nV, dimnames = list(V, V))
  Cmp[cbind(pairs$var1, pairs$var2)] <- pairs$compatible
  Cmp[cbind(pairs$var2, pairs$var1)] <- pairs$compatible

  # per-variable ingredients: novelty, gradients, agreement
  if (is.null(targets)) targets <- setdiff(names(x$sets), ref)
  S <- list()
  for (t in targets) {
    d <- x$sets[[t]]
    i <- if (nrow(d) > n_mess) sample.int(nrow(d), n_mess) else seq_len(nrow(d))
    S[[t]] <- lapply(stats::setNames(V, V), function(v) {
      mess_var(x$sets[[ref]][[v]], d[[v]][i])
    })
  }
  if (is.null(gradient_pairs)) {
    gradient_pairs <- lapply(intersect(targets, names(x$common)),
                             function(t) c(ref, t))
  }
  grad <- NULL
  if (length(gradient_pairs) && !is.null(x$common)) {
    grad <- vapply(gradient_pairs, function(p) {
      vapply(V, function(v) suppressWarnings(stats::cor(
        x$common[[p[1]]][[v]], x$common[[p[2]]][[v]], method = "spearman")),
        numeric(1))
    }, numeric(nV))
    grad <- matrix(grad, nrow = nV, dimnames = list(V, vapply(
      gradient_pairs, paste, character(1), collapse = "~")))
  }
  consistency <- NULL
  agr <- NULL
  if (!is.null(agreement)) {
    if (is.null(x$common)) {
      warning("`agreement` needs common cells; skipped.", call. = FALSE)
    } else {
      a <- x$common[[agreement[1]]]
      b <- x$common[[agreement[2]]]
      consistency <- data.frame(
        var = V,
        spearman = vapply(V, function(v) stats::cor(a[[v]], b[[v]],
                                                    method = "spearman"), 1),
        mean_a = vapply(V, function(v) mean(a[[v]]), 1),
        mean_b = vapply(V, function(v) mean(b[[v]]), 1),
        bias_b_minus_a = vapply(V, function(v) mean(b[[v]] - a[[v]]), 1),
        rmse = vapply(V, function(v) sqrt(mean((b[[v]] - a[[v]])^2)), 1),
        sd_ratio = vapply(V, function(v) stats::sd(b[[v]]) / stats::sd(a[[v]]),
                          1),
        std_bias = vapply(V, function(v) abs(mean(b[[v]] - a[[v]])) /
                            stats::sd(a[[v]]), 1),
        row.names = NULL)
      names(consistency)[3:4] <- paste0("mean_", agreement)
      agr <- stats::setNames(consistency$spearman, V)
    }
  }
  variables <- data.frame(var = V, cluster = as.integer(cl[V]),
                          vif_in_pool = vapply(V, function(v) {
                            max(vapply(rho, function(R) {
                              tryCatch(diag(solve(R))[[v]],
                                       error = function(e) Inf)
                            }, 1))
                          }, 1), row.names = NULL)
  for (t in names(S)) {
    variables[[paste0("novelty_alone_", t)]] <- vapply(S[[t]], function(s)
      mean(s < 0), 1)
  }
  if (!is.null(grad)) {
    for (k in colnames(grad)) variables[[paste0("gradient_", k)]] <- grad[, k]
  }
  if (!is.null(agr)) variables$agreement <- agr

  # enumerate all compatible subsets
  grp <- lapply(require_one_of, function(g) match(g, V))
  vif_max <- function(idx, Rl) {
    m <- 0
    for (R in Rl) {
      m <- max(m, tryCatch(max(diag(solve(R[idx, idx, drop = FALSE]))),
                           error = function(e) Inf))
    }
    m
  }
  maxr <- function(idx, Rl) {
    max(vapply(Rl, function(R) {
      m <- abs(R[idx, idx, drop = FALSE])
      diag(m) <- 0
      max(m)
    }, 1))
  }
  cand <- list()
  counts <- list()
  for (k in sizes) {
    nk <- choose(nV, k)
    if (nk > max_subsets) {
      stop(format(nk, big.mark = ","), " subsets of size ", k, ": too many. ",
           "Use fewer variables or smaller sizes, or raise `max_subsets`.",
           call. = FALSE)
    }
    cm <- utils::combn(nV, k)
    ok <- rep(TRUE, ncol(cm))
    if (k > 1) {
      for (a in 1:(k - 1)) for (b in (a + 1):k) {
        ok <- ok & Cmp[cbind(cm[a, ], cm[b, ])]
      }
    }
    n_pair <- sum(ok)
    for (g in grp) ok <- ok & colSums(matrix(cm %in% g, nrow = k)) >= 1
    n_role <- sum(ok)
    keep <- which(ok)
    keep <- keep[vapply(keep, function(j) vif_max(cm[, j], rho) <= thr_vif,
                        logical(1))]
    n_vif <- length(keep)
    if (length(reg_rho)) {
      keep <- keep[vapply(keep, function(j) {
        vif_max(cm[, j], reg_rho) <= thr_vif_regional
      }, logical(1))]
    }
    counts[[length(counts) + 1]] <- data.frame(
      size = k, all_subsets = ncol(cm), pairs_compatible = n_pair,
      with_required_groups = n_role, vif_ok = n_vif,
      regional_vif_ok = length(keep))
    for (j in keep) cand[[length(cand) + 1]] <- cm[, j]
  }
  counts <- do.call(rbind, counts)

  metrics <- function(idx) {
    m <- c(n_vars = length(idx), worst_maxr = maxr(idx, rho),
           worst_vif = vif_max(idx, rho),
           regional_maxr = if (length(reg_rho)) maxr(idx, reg_rho) else NA,
           regional_vif = if (length(reg_rho)) vif_max(idx, reg_rho) else NA)
    nov <- vapply(S, function(s) mean(Reduce(pmin, s[idx]) < 0), 1)
    m["novelty"] <- if (length(nov)) max(nov) else NA
    if (length(nov)) m[paste0("novelty_", names(nov))] <- nov
    m["gradient_min"] <- if (is.null(grad)) NA else min(grad[idx, ])
    m["agreement_min"] <- if (is.null(agr)) NA else min(agr[idx])
    m["clusters_covered"] <- length(unique(cl[idx]))
    m
  }
  if (is.null(criteria)) {
    criteria <- c(novelty = "min", gradient_min = "max",
                  agreement_min = "max", worst_vif = "min",
                  regional_vif = "min")
  }
  check_criteria(criteria, names(metrics(1)))
  res <- NULL
  if (length(cand)) {
    res <- do.call(rbind, lapply(cand, metrics))
    res <- data.frame(set = vapply(cand, function(i) paste(V[i],
                                                           collapse = " "),
                                   character(1)),
                      res, row.names = NULL, check.names = FALSE)
    res$n_vars <- as.integer(res$n_vars)
    usable <- names(criteria)[vapply(names(criteria), function(k)
      !all(is.na(res[[k]])), logical(1))]
    dropped <- setdiff(names(criteria), usable)
    if (length(dropped)) {
      message("Criteria not available for these data (left out): ",
              paste(dropped, collapse = ", "))
    }
    criteria <- criteria[usable]
    w <- stats::setNames(rep(1, length(criteria)), names(criteria))
    if (!is.null(weights)) w[intersect(names(weights), names(w))] <-
      weights[intersect(names(weights), names(w))]
    res$pareto <- FALSE
    res$mean_rank <- NA_real_
    res$rank <- NA_integer_
    for (k in unique(res$n_vars)) {
      i <- which(res$n_vars == k)
      if (length(criteria)) {
        M <- to_minimise(res[i, , drop = FALSE], criteria)
        res$pareto[i] <- pareto_fronts(M) == 1
        rk <- apply(M, 2, rank, ties.method = "average")
        rk <- matrix(rk, nrow = length(i))
        res$mean_rank[i] <- as.vector(rk %*% w) / sum(w)
      } else {
        res$pareto[i] <- TRUE
        res$mean_rank[i] <- 1
      }
      res$rank[i] <- rank(res$mean_rank[i], ties.method = "min")
    }
    res <- res[order(res$n_vars, res$mean_rank), ]
    rownames(res) <- NULL
  } else {
    warning("No set passes the screen: relax thr_r / thr_vif, the sizes or ",
            "require_one_of.", call. = FALSE)
  }

  # how often each variable appears
  has <- function(r, v) if (is.null(r) || !nrow(r)) logical(0) else
    grepl(paste0("(^| )", v, "( |$)"), r$set)
  top <- if (is.null(res)) NULL else res[res$rank <= 10, , drop = FALSE]
  freq <- data.frame(var = V, cluster = as.integer(cl[V]),
                     in_passing_sets = vapply(V, function(v)
                       mean(has(res, v)), 1),
                     in_pareto_sets = vapply(V, function(v)
                       mean(has(res[res$pareto, , drop = FALSE], v)), 1),
                     in_top_ranked = vapply(V, function(v)
                       mean(has(top, v)), 1), row.names = NULL)
  freq[is.na(freq)] <- 0
  freq <- freq[order(-freq$in_pareto_sets, -freq$in_passing_sets), ]
  rownames(freq) <- NULL

  # reference sets under the same metrics
  refs <- NULL
  if (length(reference_sets)) {
    refs <- do.call(rbind, lapply(names(reference_sets), function(nm) {
      v <- reference_sets[[nm]]
      if (!all(v %in% V)) stop("Reference set '", nm, "' has unknown ",
                               "variables.", call. = FALSE)
      idx <- match(v, V)
      m <- metrics(idx)
      data.frame(name = nm, set = paste(v, collapse = " "), t(m),
                 check.names = FALSE, row.names = NULL)
    }))
    refs$passes <- refs$worst_maxr < thr_r & refs$worst_vif <= thr_vif &
      (is.na(refs$regional_vif) | refs$regional_vif <= thr_vif_regional)
  }

  # similarity of the correlation structure between data sets
  lt <- lower.tri(rho[[1]])
  agree <- outer(names(rho), names(rho), Vectorize(function(a, b) {
    stats::cor(abs(rho[[a]][lt]), abs(rho[[b]][lt]))
  }))
  dimnames(agree) <- list(names(rho), names(rho))

  out <- list(sets = res, counts = counts, pairs = pairs,
              variables = variables, frequency = freq,
              reference_sets = refs, consistency = consistency,
              structure_agreement = agree, tree = tree,
              clusters = data.frame(var = V, cluster = as.integer(cl[V])),
              rho = rho, regional_rho = reg_rho, union = A,
              settings = list(sizes = sizes, thr_r = thr_r, thr_vif = thr_vif,
                              thr_r_regional = thr_r_regional,
                              thr_vif_regional = thr_vif_regional,
                              criteria = criteria, targets = targets,
                              reference = ref, agreement = agreement,
                              require_one_of = require_one_of))
  class(out) <- "var_selection"
  out
}

# MESS similarity of one variable (Elith et al. 2010): negative = outside
# the reference range
mess_var <- function(ref, p) {
  ref <- ref[!is.na(ref)]
  mn <- min(ref)
  mx <- max(ref)
  if (mx == mn) return(ifelse(p == mn, 100, -100))
  f <- stats::ecdf(ref)(p) * 100
  ifelse(p < mn, (p - mn) / (mx - mn) * 100,
         ifelse(p > mx, (mx - p) / (mx - mn) * 100,
                ifelse(f <= 50, 2 * f, 2 * (100 - f))))
}

#' The best variable sets of a selection
#'
#' @param x A `var_selection` from [select_variables()].
#' @param n Number of sets per size.
#' @param size Optional set size(s).
#' @param pareto_only Only sets on the Pareto front.
#' @return A data.frame of sets ordered by size and mean rank, with their
#'   scores.
#' @export
#' @examples
#' # see ?select_variables
best_sets <- function(x, n = 5, size = NULL, pareto_only = FALSE) {
  if (!inherits(x, "var_selection")) {
    stop("`x` must come from select_variables().", call. = FALSE)
  }
  r <- x$sets
  if (is.null(r)) return(r)
  if (!is.null(size)) r <- r[r$n_vars %in% size, , drop = FALSE]
  if (pareto_only) r <- r[r$pareto, , drop = FALSE]
  r <- do.call(rbind, lapply(split(r, r$n_vars), utils::head, n))
  cols <- c("set", "n_vars", "rank", "mean_rank", "pareto",
            names(x$settings$criteria), "worst_maxr")
  r <- r[intersect(unique(cols), names(r))]
  rownames(r) <- NULL
  r
}
