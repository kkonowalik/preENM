# ODMAP protocol filled from the objects of a workflow ---------------------------

#' Fill the ODMAP reporting protocol from your preENM workflow
#'
#' ODMAP (Overview, Data, Model, Assessment, Prediction; Zurell et al.
#' 2020) is the standard protocol for reporting species distribution
#' models, required or recommended by several journals. Much of its Data
#' section, and parts of the others, describe exactly the steps this
#' package performs. The function writes those entries from the objects
#' you already have (which tests removed how many records, how M was
#' delimited, how the records were thinned, how the predictors were
#' screened, how the models were evaluated and selected), so that they are
#' reported completely and identically to what was run. The remaining
#' elements (the study, the taxon, the algorithms and their settings) are
#' left empty for you to fill, or supplied through `fields`.
#'
#' Every argument is optional: pass what your workflow used.
#'
#' @param filter An `occ_filter` from [filter_occurrences()].
#' @param uncertainty The result of [uncertainty_environment()].
#' @param m An `m_area` from [delimit_m()].
#' @param thinning An `env_thinning` from [thin_environment()].
#' @param records The final presence records used for modelling (a
#'   data.frame; only their number is reported).
#' @param background Background points from [sample_background()].
#' @param variables A `var_selection` from [select_variables()].
#' @param predictors Names of the predictors finally used.
#' @param truncation A `niche_truncation` from [niche_truncation()].
#' @param block_size A `block_size` from [suggest_block_size()].
#' @param folds Fold numbers from [spatial_folds()].
#' @param evaluation An `enm_evaluation` from [evaluate_models()] or an
#'   `enm_cv` from [evaluate_cv()].
#' @param selection An `enm_selection` from [select_models()].
#' @param null_test The result of [null_model_test()].
#' @param fields Named list of your own entries, named by element, e.g.
#'   `list("Study title" = "...", "Focal Taxon" = "...")`. Where an element
#'   name occurs in several subsections, write
#'   `"Subsection: Element"`, e.g. `"Transfer data: Data sources"`. Your
#'   entry replaces the automatic one.
#' @param file Optional path of a CSV file to write.
#' @return A data.frame with `section`, `subsection`, `element`, `value` and
#'   `filled_by` (`"preENM"`, `"user"` or empty), invisibly if written to a
#'   file. The sections, subsections and elements follow ODMAP version 1.0;
#'   compare them with the template of the journal you submit to.
#' @references Zurell, D. et al. (2020) A standard protocol for reporting
#'   species distribution models. *Ecography* 43, 1261-1277.
#' @export
#' @examples
#' sim <- simulate_enm_data(seed = 1, n_presence = 150)
#' m <- delimit_m(sim$presence, link_km = 120, thin_km = 5)
#' bg <- sample_background(sim$env, 500, seed = 1)
#' ev <- evaluate_models(sim$models, sim$presence, sim$absence)
#' od <- odmap_table(m = m, records = sim$presence, background = bg,
#'                   evaluation = ev, predictors = names(sim$env),
#'                   fields = list("Focal Taxon" = "a simulated species"))
#' od[od$filled_by != "", c("subsection", "element", "value")]
odmap_table <- function(filter = NULL, uncertainty = NULL, m = NULL,
                        thinning = NULL, records = NULL, background = NULL,
                        variables = NULL, predictors = NULL,
                        truncation = NULL, block_size = NULL, folds = NULL,
                        evaluation = NULL, selection = NULL,
                        null_test = NULL, fields = list(), file = NULL) {
  od <- odmap_elements()
  od$value <- ""
  od$filled_by <- ""
  add <- function(subsection, element, ...) {
    i <- which(od$subsection == subsection & od$element == element)
    if (length(i) != 1) stop("Internal: unknown ODMAP element ", subsection,
                             ": ", element, call. = FALSE)
    txt <- paste0(...)
    od$value[i] <<- if (od$value[i] == "") txt else paste(od$value[i], txt)
    od$filled_by[i] <<- "preENM"
  }
  kv <- function(x) paste(names(x), unlist(x), sep = " = ", collapse = ", ")

  add("Software", "Software", R.version.string, "; preENM ",
      as.character(utils::packageVersion("preENM")), ".")

  if (!is.null(filter)) {
    check_filter(filter)
    d <- filter$data
    lg <- filter$log
    run <- lg[lg$status == "run" & !is.na(lg$removed) & lg$removed > 0, ]
    add("Biodiversity data", "Sample size", nrow(d),
        " records before cleaning, ", sum(is.na(d$reason)), " after.")
    add("Biodiversity data", "Cleaning",
        "Rule-based filtering (preENM::filter_occurrences); records removed ",
        "per test: ", if (nrow(run)) paste0(run$reason, " (", run$removed, ")",
                                            collapse = ", ") else "none",
        ". Tests run without removals: ",
        paste(lg$test[lg$status == "run" & !is.na(lg$removed) &
                        lg$removed == 0], collapse = ", "),
        ". Tests skipped: ",
        if (any(lg$status == "skipped")) paste0(
          lg$test[lg$status == "skipped"], " [",
          lg$note[lg$status == "skipped"], "]", collapse = ", ") else "none",
        ".")
  }
  if (!is.null(uncertainty)) {
    if (!"unc_ok" %in% names(uncertainty)) {
      stop("`uncertainty` must come from uncertainty_environment().",
           call. = FALSE)
    }
    add("Biodiversity data", "Errors and biases",
        "Coordinate uncertainty was judged by the range of each predictor ",
        "among the cells within the uncertainty radius, relative to the ",
        "niche breadth (preENM::uncertainty_environment): ",
        sum(uncertainty$unc_ok, na.rm = TRUE), " records usable, ",
        sum(!uncertainty$unc_ok, na.rm = TRUE), " too uncertain, ",
        sum(is.na(uncertainty$unc_ok)), " undecided.")
  }
  if (!is.null(m)) {
    if (!inherits(m, "m_area")) stop("`m` must come from delimit_m().",
                                     call. = FALSE)
    st <- m$settings
    area <- m$report$area_km2[m$report$part == "M"]
    add("Biodiversity data", "Clipping",
        "Accessible area (M) delimited with preENM::delimit_m, method ",
        st$method, if (!is.null(st$L)) paste0(", linking distance ",
                                              signif(st$L, 4), " km"),
        if (!is.null(st$buffer_km)) paste0(", buffer ", st$buffer_km, " km"),
        "; records thinned to ", st$thin_km, " km for the delimitation; ",
        "area ", format(round(area), big.mark = ","), " km2",
        if (length(st$expert)) paste0("; expert additions: ",
                                      paste(st$expert, collapse = ", ")),
        ".")
    add("Scale of Analysis", "Boundary", "Accessible area (M) of the species.")
  }
  if (!is.null(thinning)) {
    if (!inherits(thinning, "env_thinning")) {
      stop("`thinning` must come from thin_environment().", call. = FALSE)
    }
    st <- thinning$settings
    add("Biodiversity data", "Scaling",
        "Records thinned in environmental space (preENM::thin_environment): ",
        sum(thinning$data$kept), " of ", st$n_input, " kept; ",
        if (st$axes == "pca") paste0("PC1-PC", st$n_axes) else
          "standardised layers", " of ", paste(st$vars, collapse = ", "),
        "; at most ", st$max_per_bin, " per bin; at least ", st$block_min,
        " per ", st$block_km, " km block; extremes of every layer kept.")
  }
  if (!is.null(records)) {
    add("Biodiversity data", "Sample size", NROW(records),
        " presence records used for modelling.")
  }
  if (!is.null(background)) {
    info <- attr(background, "background")
    sets <- if ("set" %in% names(background)) length(unique(background$set))
            else 1
    add("Biodiversity data", "Background data",
        NROW(background) / sets, " background points",
        if (sets > 1) paste0(" in each of ", sets, " sets"),
        " (preENM::sample_background, method ",
        if (is.null(info$method)) "unknown" else info$method,
        if (!is.null(info$n_bins)) paste0(", ", info$n_bins,
                                          " environmental bins on PC1-PC",
                                          info$n_axes),
        if (!is.null(info$effort)) paste0(", recording effort with ",
                                          kv(as.list(info$effort))), ").")
    if (!is.null(info$effort)) {
      add("Biodiversity data", "Errors and biases",
          "Sampling bias was addressed by drawing background points in ",
          "proportion to recording effort.")
    }
  }
  if (!is.null(predictors)) {
    add("Predictor variables", "Predictor variables",
        paste(predictors, collapse = ", "), ".")
  }
  if (!is.null(variables)) {
    if (!inherits(variables, "var_selection")) {
      stop("`variables` must come from select_variables().", call. = FALSE)
    }
    st <- variables$settings
    add("Variable pre-selection", "Variable pre-selection",
        "All subsets of ", paste(range(st$sizes), collapse = " to "),
        " variables among ", nrow(variables$variables),
        " candidates were enumerated (preENM::select_variables); ",
        if (is.null(variables$sets)) 0 else nrow(variables$sets),
        " passed the collinearity screen and were ranked on ",
        kv(as.list(st$criteria)), ".")
    add("Multicollinearity", "Multicollinearity",
        "Sets were required to have |Spearman r| < ", st$thr_r,
        " and VIF <= ", st$thr_vif, " in every data set (",
        paste(names(variables$rho), collapse = ", "), ")",
        if (length(variables$regional_rho)) paste0(
          " and |r| < ", st$thr_r_regional, ", VIF <= ", st$thr_vif_regional,
          " within ", length(variables$regional_rho), " regions"), ".")
    if (length(st$targets)) {
      add("Transfer data", "Quantification of Novelty",
          "Share of cells outside the reference range of ", st$reference,
          " (MESS < 0) computed for ", paste(st$targets, collapse = ", "),
          " and used to rank the predictor sets.")
    }
  }
  if (!is.null(truncation)) {
    if (!inherits(truncation, "niche_truncation")) {
      stop("`truncation` must come from niche_truncation().", call. = FALSE)
    }
    v <- truncation$variables
    ends <- c(if (any(v$truncated_low))
      paste(v$layer[v$truncated_low], "(low end)"),
      if (any(v$truncated_high))
        paste(v$layer[v$truncated_high], "(high end)"))
    add("Biodiversity data", "Errors and biases",
        "Niche truncation by the calibration area was checked ",
        "(preENM::niche_truncation): ",
        if (length(ends)) paste0("truncated on ", paste(ends, collapse = ", "))
        else "no predictor truncated", ".")
    p <- truncation$projections
    if (!is.null(p)) {
      worst <- tapply(pmax(p$below_m, p$above_m), p$projection, max)
      risky <- p[p$risk_low | p$risk_high, , drop = FALSE]
      txt <- paste0(
        "Largest share of cells outside the calibration range of any ",
        "predictor: ", paste0(names(worst), " ", round(100 * worst, 1), "%",
                              collapse = ", "), ". Truncated responses ",
        "extrapolated: ", if (nrow(risky)) paste0(
          risky$layer, " in ", risky$projection, collapse = ", ") else
            "none", ".")
      add("Transfer data", "Quantification of Novelty", txt)
      add("Uncertainty quantification", "Novel environments", txt)
    }
  }
  if (!is.null(block_size)) {
    if (!inherits(block_size, "block_size")) {
      stop("`block_size` must come from suggest_block_size().", call. = FALSE)
    }
    add("Analysis and Correction of non-independence",
        "Spatial autocorrelation",
        "Range of spatial autocorrelation of the predictors (median over ",
        "layers, semivariance reaching ",
        100 * block_size$settings$sill_share, "% of the variance): ",
        block_size$block_size_km, " km.")
  }
  if (!is.null(folds)) {
    add("Data partitioning", "Validation data",
        "Spatial block cross-validation (preENM::spatial_folds): ",
        length(unique(folds)), " folds",
        if (!is.null(attr(folds, "block_size_km"))) paste0(
          ", blocks of ", attr(folds, "block_size_km"), " km"),
        ", blocks assigned to folds so that presences and absences are ",
        "balanced.")
    add("Analysis and Correction of non-independence",
        "Spatial autocorrelation",
        "Models were evaluated on spatially separated blocks.")
  }
  if (!is.null(evaluation)) {
    mets <- intersect(names(default_directions), names(evaluation))
    crit <- attr(evaluation, "criteria")
    nb <- attr(evaluation, "n_boot")
    cv <- inherits(evaluation, "enm_cv")
    add("Performance statistics",
        if (cv) "Performance on validation data" else
          "Performance on test data",
        paste(mets, collapse = ", "), " for ", nrow(evaluation), " models",
        if (cv) " (means over folds)" else "",
        if (!is.null(nb) && nb > 0) paste0("; ", nb, " stratified bootstrap ",
                                           "replicates for intervals"),
        if (!is.null(crit)) paste0("; ranked by Pareto fronts on ",
                                   kv(as.list(crit))),
        " (preENM).")
  }
  if (!is.null(selection)) {
    if (!inherits(selection, "enm_selection")) {
      stop("`selection` must come from select_models().", call. = FALSE)
    }
    mins <- attr(selection, "min_metrics")
    maxs <- attr(selection, "max_metrics")
    add("Model selection - model averaging - ensembles", "Model selection",
        "Rule: ", attr(selection, "rule"),
        if (length(mins)) paste0("; required ", paste(names(mins), mins,
                                                      sep = " >= ",
                                                      collapse = ", ")),
        if (length(maxs)) paste0("; required ", paste(names(maxs), maxs,
                                                      sep = " <= ",
                                                      collapse = ", ")),
        "; Pareto front on ", kv(as.list(attr(selection, "criteria"))),
        "; selected: ", paste(selection$model[selection$selected],
                              collapse = ", "), " (preENM::select_models).")
  }
  if (!is.null(null_test)) {
    if (!all(c("model", "metric", "p_value") %in% names(null_test))) {
      stop("`null_test` must come from null_model_test().", call. = FALSE)
    }
    add("Performance statistics", "Performance on test data",
        "Significance against random locations",
        if (isTRUE(attr(null_test, "weighted"))) " drawn by recording effort",
        " (", attr(null_test, "n_null"), " replicates, ",
        "preENM::null_model_test): ",
        paste0(null_test$model, " ", null_test$metric, " p = ",
               signif(null_test$p_value, 2), collapse = ", "), ".")
  }

  # the user's entries replace the automatic ones
  for (nm in names(fields)) {
    key <- tolower(trimws(strsplit(nm, ":", fixed = TRUE)[[1]]))
    i <- if (length(key) == 2) {
      which(tolower(od$subsection) == key[1] & tolower(od$element) == key[2])
    } else which(tolower(od$element) == key[1])
    if (!length(i)) stop("Unknown ODMAP element: ", nm, call. = FALSE)
    if (length(i) > 1) {
      stop("'", nm, "' occurs in several subsections (",
           paste(od$subsection[i], collapse = ", "),
           "): write \"Subsection: ", od$element[i[1]], "\".", call. = FALSE)
    }
    od$value[i] <- paste(as.character(fields[[nm]]), collapse = "; ")
    od$filled_by[i] <- "user"
  }
  if (!is.null(file)) {
    utils::write.csv(od, file, row.names = FALSE)
    return(invisible(od))
  }
  od
}

# Sections, subsections and elements of ODMAP 1.0 (internal)
odmap_elements <- function() {
  x <- list(
    Overview = list(
      "Authorship" = c("Study title", "Author names", "Contact",
                       "Study link"),
      "Model objective" = c("Model objective", "Target output"),
      "Focal Taxon" = "Focal Taxon",
      "Location" = "Location",
      "Scale of Analysis" = c("Spatial extent", "Spatial resolution",
                              "Temporal extent", "Temporal resolution",
                              "Boundary"),
      "Biodiversity data overview" = c("Observation type",
                                       "Response data type"),
      "Predictors" = "Predictor types",
      "Hypotheses" = "Hypotheses",
      "Assumptions" = "Model assumptions",
      "Algorithms" = c("Modelling techniques", "Model complexity",
                       "Model averaging"),
      "Workflow" = "Model workflow",
      "Software" = c("Software", "Code availability", "Data availability")),
    Data = list(
      "Biodiversity data" = c("Taxon names", "Taxonomic reference system",
                              "Ecological level", "Data sources",
                              "Sampling design", "Sample size", "Clipping",
                              "Scaling", "Cleaning", "Absence data",
                              "Background data", "Errors and biases"),
      "Data partitioning" = c("Training data", "Validation data",
                              "Test data"),
      "Predictor variables" = c("Predictor variables", "Data sources",
                                "Spatial extent", "Spatial resolution",
                                "Coordinate reference system",
                                "Temporal extent", "Temporal resolution",
                                "Data processing", "Errors and biases",
                                "Dimension reduction"),
      "Transfer data" = c("Data sources", "Spatial extent",
                          "Spatial resolution", "Temporal extent",
                          "Temporal resolution", "Models and scenarios",
                          "Data processing", "Quantification of Novelty")),
    Model = list(
      "Variable pre-selection" = "Variable pre-selection",
      "Multicollinearity" = "Multicollinearity",
      "Model settings" = c("Model settings (fitting)",
                           "Model settings (extrapolation)"),
      "Model estimates" = c("Coefficients", "Parameter uncertainty",
                            "Variable importance"),
      "Model selection - model averaging - ensembles" =
        c("Model selection", "Model averaging", "Model ensembles"),
      "Analysis and Correction of non-independence" =
        c("Spatial autocorrelation", "Temporal autocorrelation",
          "Nested data"),
      "Threshold selection" = "Threshold selection"),
    Assessment = list(
      "Performance statistics" = c("Performance on training data",
                                   "Performance on validation data",
                                   "Performance on test data"),
      "Plausibility check" = c("Response shapes", "Expert judgement")),
    Prediction = list(
      "Prediction output" = c("Prediction unit", "Post-processing"),
      "Uncertainty quantification" = c("Algorithmic uncertainty",
                                       "Input data uncertainty",
                                       "Parameter uncertainty",
                                       "Scenario uncertainty",
                                       "Novel environments")))
  do.call(rbind, lapply(names(x), function(s) {
    do.call(rbind, lapply(names(x[[s]]), function(sub) {
      data.frame(section = s, subsection = sub, element = x[[s]][[sub]],
                 stringsAsFactors = FALSE)
    }))
  }))
}
