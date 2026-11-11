cl <- make_climate()

test_that("block-wise reading equals terra", {
  r <- cl$lgm
  v <- terra::values(r)
  expect_equal(valid_cells(r), which(stats::complete.cases(v)))
  cells <- sort(sample(terra::ncell(r), 300))
  expect_equal(unname(read_cells(r, cells)), unname(v[cells, ]))
})

test_that("layers are matched by name tokens", {
  r <- cl$present[[1:2]]
  names(r) <- c("CHELSA_TraCE21k_bio01_-100_V1.0", "CHELSA_TraCE21k_bio04_-100_V1.0")
  m <- match_layers(r, c("bio04", "bio01"), "x")
  expect_equal(names(m), c("bio04", "bio01"))
  names(r) <- c("bio01", "bio010")
  expect_equal(names(match_layers(r, "bio01", "x")), "bio01")   # exact first
  expect_error(match_layers(r, "bio12", "x"), "matches 0")
})

test_that("MESS of one variable follows Elith et al. (2010)", {
  ref <- 1:100
  expect_equal(mess_var(ref, c(50, 1, 100)), c(100, 2, 0))
  expect_equal(mess_var(ref, c(-9.9, 109.9)), c(-11, -10), tolerance = 0.01)
  expect_true(all(mess_var(ref, c(-5, 200)) < 0))
})

test_that("sampling: own and common cells, anomalies, regions, weighting", {
  es <- sample_environment(
    list(CHELSA = cl$present, TraCE = cl$trace, LGM = cl$lgm), n = 1500,
    regions = cl$regions, seed = 2, verbose = FALSE,
    anomalies = list(LGM_delta = c(current = "CHELSA", past = "LGM",
                                   model_present = "TraCE")),
    ratio_vars = "bio12", bounds = list(scd = c(-1, 1)))
  expect_s3_class(es, "env_sample")
  expect_equal(names(es$sets), c("CHELSA", "TraCE", "LGM", "LGM_delta"))
  expect_true(all(vapply(es$sets[1:3], nrow, 1L) == 1500))
  # common cells hold the same places in every data set
  expect_equal(es$common$CHELSA[c(".x", ".y")], es$common$LGM[c(".x", ".y")])
  # LGM ice sheet: no own LGM cell there
  expect_false(any(es$sets$LGM$.y > 62 & es$sets$LGM$.x < 20))
  # anomaly = anomaly_climate() at the common cells, with bounds
  V <- es$vars
  ref <- anomaly_climate(es$common$CHELSA[V], es$common$LGM[V],
                         es$common$TraCE[V], ratio_vars = "bio12")
  expect_equal(es$common$LGM_delta$bio12, ref$bio12)
  expect_equal(es$common$LGM_delta$scd, pmin(pmax(ref$scd, -1), 1))
  # regions assigned from the polygons
  d <- es$sets$CHELSA
  alps <- d$.region %in% "Alps"
  expect_true(all(d$.x[alps] >= 5 & d$.x[alps] <= 17 & d$.y[alps] <= 49))
  expect_output(print(es), "Regions")
  # area weighting: fewer high-latitude cells than an unweighted draw
  uw <- sample_environment(list(a = cl$present), n = 3000, seed = 3,
                           area_weighted = FALSE, verbose = FALSE)
  aw <- sample_environment(list(a = cl$present), n = 3000, seed = 3,
                           verbose = FALSE)
  expect_lt(mean(aw$sets$a$.y), mean(uw$sets$a$.y))
  # tiles as regions
  et <- sample_environment(list(a = cl$present), n = 500, seed = 1,
                           regions = list(x = c(-10, 10, 35), y = c(40, 55, 70)),
                           verbose = FALSE)
  expect_setequal(unique(et$sets$a$.region),
                  c("x1_y1", "x1_y2", "x2_y1", "x2_y2"))
})

test_that("scale check warns and fixes a variable stored x 10", {
  # a copy of the present with bio12 stored x 10 (spread x 10)
  cur <- cl$present
  past <- cl$present
  past[["bio12"]] <- past[["bio12"]] * 10
  expect_warning(sample_environment(list(now = cur, past = past), n = 500,
                                    seed = 1, verbose = FALSE), "x 10")
  es <- suppressMessages(sample_environment(list(now = cur, past = past),
                                            n = 500, seed = 1,
                                            fix_scale = TRUE))
  expect_equal(stats::IQR(es$sets$past$bio12), stats::IQR(es$sets$now$bio12),
               tolerance = 0.2)
  # a shifted climate with the same spread raises no alarm
  shifted <- cl$present
  shifted[["bio01"]] <- shifted[["bio01"]] - 5
  expect_silent(sample_environment(list(now = cur, past = shifted), n = 500,
                                   seed = 1, verbose = FALSE))
})

test_that("the screen keeps exactly the sets a brute-force check keeps", {
  es <- sample_environment(list(CHELSA = cl$present, LGM = cl$lgm), n = 1500,
                           regions = cl$regions, seed = 4, verbose = FALSE)
  grp <- list(t = c("bio01", "bio05", "bio06", "bio10", "bio11", "scd"),
              m = c("bio12", "bio14"))
  vs <- select_variables(es, sizes = 2:4, min_region_n = 100,
                         require_one_of = grp)
  V <- es$vars
  mats <- c(lapply(es$sets, function(d) cor(d[V], method = "spearman")),
            list(cor(es$common$LGM[V], method = "spearman")))
  regm <- list()
  for (k in names(es$sets)) {
    d <- es$sets[[k]]
    for (rg in c("Alps", "Scandinavia")) {
      i <- which(d$.region == rg)
      if (length(i) >= 100) regm[[length(regm) + 1]] <- cor(d[i, V],
                                                            method = "spearman")
    }
  }
  passes <- function(s) {
    for (R in mats) {
      m <- abs(R[s, s]); diag(m) <- 0
      if (max(m) >= 0.8) return(FALSE)
      if (max(diag(solve(R[s, s]))) > 5) return(FALSE)
    }
    for (R in regm) {
      m <- abs(R[s, s]); diag(m) <- 0
      if (max(m) >= 0.85) return(FALSE)
      if (max(diag(solve(R[s, s]))) > 10) return(FALSE)
    }
    all(vapply(grp, function(g) any(s %in% g), logical(1)))
  }
  brute <- unlist(lapply(2:4, function(k) {
    cm <- utils::combn(V, k)
    keep <- apply(cm, 2, passes)
    apply(cm[, keep, drop = FALSE], 2, paste, collapse = " ")
  }))
  expect_setequal(vs$sets$set, brute)
  # no passing set holds two variables of one cluster of the union tree
  expect_true(all(vs$sets$clusters_covered == vs$sets$n_vars))
  expect_equal(vs$counts$regional_vif_ok, as.vector(table(
    factor(vs$sets$n_vars, levels = 2:4))))
})

test_that("scores, Pareto front and ranking", {
  es <- sample_environment(list(CHELSA = cl$present, TraCE = cl$trace,
                                LGM = cl$lgm), n = 1500, seed = 5,
                           verbose = FALSE)
  vs <- select_variables(es, sizes = 3, targets = "LGM",
                         agreement = c("CHELSA", "TraCE"),
                         reference_sets = list(bad = c("bio01", "bio10",
                                                       "bio12")))
  s <- vs$sets
  V <- es$vars
  # novelty of the first set, computed by hand
  mem <- strsplit(s$set[1], " ")[[1]]
  m <- Reduce(pmin, lapply(mem, function(v) mess_var(es$sets$CHELSA[[v]],
                                                     es$sets$LGM[[v]])))
  expect_equal(s$novelty[1], mean(m < 0))
  # gradient: rank correlation between CHELSA and LGM on common cells
  g <- vapply(mem, function(v) cor(es$common$CHELSA[[v]], es$common$LGM[[v]],
                                   method = "spearman"), 1)
  expect_equal(s$gradient_min[1], min(g))
  expect_equal(vs$consistency$var, V)
  expect_equal(s$agreement_min[1], min(vs$consistency$spearman[
    match(mem, V)]))
  # Pareto front checked by brute force
  crit <- vs$settings$criteria
  M <- to_minimise(s, crit)
  dominated <- vapply(seq_len(nrow(M)), function(i) any(apply(M, 1, function(o)
    all(o <= M[i, ]) && any(o < M[i, ]))), logical(1))
  expect_equal(s$pareto, !dominated)
  expect_true(all(diff(s$mean_rank) >= 0))
  expect_equal(min(s$rank), 1)
  # reference set flagged as failing (bio01 and bio10 are collinear)
  expect_false(vs$reference_sets$passes)
  b <- best_sets(vs, n = 2)
  expect_equal(nrow(b), 2)
  expect_equal(nrow(best_sets(vs, pareto_only = TRUE)),
               min(5, sum(s$pareto)))
  expect_output(print(vs), "Best sets")
  expect_equal(diag(vs$structure_agreement), rep(1, nrow(
    vs$structure_agreement)), ignore_attr = TRUE)
})

test_that("pair table and frequencies", {
  es <- sample_environment(list(CHELSA = cl$present, LGM = cl$lgm), n = 1500,
                           seed = 6, verbose = FALSE)
  vs <- select_variables(es, sizes = 3)
  p <- vs$pairs
  expect_equal(nrow(p), choose(length(es$vars), 2))
  expect_true(all(p$compatible == (p$max_abs < 0.8)))
  expect_true(all(p$max_abs >= p$min_abs))
  expect_true(all(p$driver %in% names(vs$rho)))
  f <- vs$frequency
  expect_true(all(f$in_passing_sets >= 0 & f$in_passing_sets <= 1))
  expect_equal(sum(f$in_passing_sets), 3)        # each set has 3 variables
})

test_that("works with one data set, different grids and impossible screens", {
  one <- sample_environment(list(now = cl$present), n = 1000, seed = 1,
                            verbose = FALSE)
  expect_message(vs1 <- select_variables(one, sizes = 2:3), "left out")
  expect_equal(names(vs1$settings$criteria), "worst_vif")
  coarse <- terra::aggregate(cl$lgm, 2)
  es2 <- sample_environment(list(now = cl$present, past = coarse), n = 800,
                            seed = 1, verbose = FALSE)
  expect_null(es2$common)
  vs2 <- suppressMessages(select_variables(es2, sizes = 3))
  expect_true(all(is.na(vs2$sets$gradient_min)))
  expect_false(all(is.na(vs2$sets$novelty)))
  expect_warning(vs3 <- select_variables(es2, sizes = 3, thr_r = 0.05),
                 "No set passes")
  expect_null(vs3$sets)
  expect_error(select_variables(es2, sizes = 3, max_subsets = 10), "too many")
  expect_error(select_variables(es2, sizes = 20), "exceeds")
  expect_error(select_variables(es2, targets = "zz"), "Unknown data set")
  expect_error(sample_environment(list(cl$present, cl$lgm)), "named list")
})

test_that("plots and files", {
  es <- sample_environment(list(CHELSA = cl$present, LGM = cl$lgm), n = 1000,
                           seed = 7, verbose = FALSE)
  vs <- select_variables(es, sizes = 2:3)
  expect_s3_class(plot_variable_tree(vs, best_sets(vs)$set[1]), "ggplot")
  expect_s3_class(plot_correlation(vs), "ggplot")
  expect_s3_class(plot_correlation(vs, "LGM"), "ggplot")
  expect_error(plot_correlation(vs, "nope"), "union")
  expect_s3_class(plot_variable_sets(vs), "ggplot")
  expect_error(plot_variable_sets(vs, "nope"), "Unknown score")
  f <- write_variable_selection(vs, tempfile(), "t")
  expect_true(all(file.exists(f)))
  nwk <- readLines(f[["_tree.nwk"]])
  expect_match(nwk, ";$")
  expect_true(all(vapply(es$vars, grepl, logical(1), x = nwk, fixed = TRUE)))
  expect_error(best_sets(1), "select_variables")
})
