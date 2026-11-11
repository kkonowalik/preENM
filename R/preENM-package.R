#' preENM: data preparation and model selection for ENM/SDM
#'
#' Data preparation:
#' * [read_gbif()], [filter_occurrences()], [filter_settings()],
#'   [filter_layers()]: rule-based filtering of occurrence records.
#' * [delimit_m()], [linking_distance()], [m_scenarios()], [plot_m_area()],
#'   [write_m_area()], [crop_to_m()]: the accessible area (M).
#' * [merge_occurrences()], [thin_by_cell()], [sample_background()]: one
#'   record per raster cell and background points.
#'
#' Model evaluation and selection:
#' * [evaluate_models()]: AUC, MAE and bias for one or many models, Pareto
#'   ranking and bootstrap robustness.
#' * [plot_pareto()]: publication-ready AUC-MAE plot with the Pareto front.
#' * [extract_predictions()]: model predictions at evaluation points.
#' * [pareto_rank()]: generic multi-criteria Pareto ranking.
#' * [simulate_enm_data()]: example data.
#'
#' Please cite Konowalik & Nosol (2021) \doi{10.1038/s41598-020-80062-1}
#' when using this package (see `citation("preENM")`).
#'
#' @keywords internal
"_PACKAGE"
