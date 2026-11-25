#' preENM: data preparation and model selection for ENM/SDM
#'
#' Data preparation:
#' * [match_names()]: scientific names matched to GBIF, ChecklistBank or a
#'   local backbone (accepted name, rank, synonyms, fuzzy matches).
#' * [read_gbif()], [filter_occurrences()], [filter_settings()],
#'   [filter_layers()]: rule-based filtering of occurrence records,
#'   including tests of the names.
#' * [delimit_m()], [linking_distance()], [m_scenarios()], [plot_m_area()],
#'   [write_m_area()], [crop_to_m()]: the accessible area (M).
#' * [crop_rasters()]: cut folders of raster files to M, keeping values exact.
#' * [download_chelsa_monthly()], [download_chelsa_trace21k()],
#'   [monthly_files()], [climate_indices()]: twenty-five derived climate
#'   variables (water balance, heat, cold, season and light; see
#'   [climate_index_names()]) from monthly climatologies of the present,
#'   the future or the past.
#' * [uncertainty_environment()]: coordinate uncertainty judged by how much
#'   the predictors vary within it (keep, propagate or drop).
#' * [jitter_uncertainty()], [plot_jitter()], [write_jitter()]: replicate
#'   datasets with records moved within their uncertainty, to propagate it
#'   into the models.
#' * [niche_truncation()], [plot_niche_truncation()]: is the niche cut off
#'   by the calibration area, and do projections go beyond it?
#' * [niche_coverage()], [plot_niche_coverage()], [classify_points()],
#'   [unfilled_records()]: how much of a widespread species' niche a region
#'   of interest contains; is a nested (hierarchical) model needed, and
#'   would it help the points you predict for?
#' * [extrapolation_risk()]: maps of novel conditions in a projection
#'   (ExDet and MOP).
#' * [check_sample_size()], [match_climate_period()]: enough records, and
#'   from the period of the climate data?
#' * [sample_environment()], [select_variables()], [best_sets()],
#'   [recommend_sets()]: predictor
#'   sets uncorrelated in every period, scenario and region, ranked.
#' * [variable_relevance()], [random_variables()], [compare_to_random()]:
#'   are the predictors more informative than random variables, before
#'   modelling and inside the models?
#' * [grain_sensitivity()], [plot_grain_sensitivity()],
#'   [write_grain_sensitivity()]: how records kept, records too uncertain,
#'   correlation structure and the selected variable set change at coarser
#'   grains, and up to which grain the selection is stable.
#' * [merge_occurrences()], [thin_by_cell()]: one record per raster cell;
#'   the data type of each source (presence-only, presence-absence,
#'   abundance) is kept through filtering.
#' * [export_integrated()], [write_integrated()]: input for integrated SDMs
#'   (intSDM, PointedSDMs, biomod2), one dataset per data type.
#' * [thin_environment()], [plot_env_thinning()], [write_env_thinning()]:
#'   thinning in environmental space without truncating the niche.
#' * [sample_background()]: background points, uniform, stratified in
#'   environmental space, by recording effort, or both; effort layers from
#'   GBIF with [gbif_effort_sql()], [download_gbif_effort()] and
#'   [gbif_effort_rasters()].
#' * [download_neotoma()], [neotoma_records()], [fossil_site_bins()],
#'   [site_bin_slices()]: dated fossil records for [temporal_validation()].
#'
#' Model evaluation and selection:
#' * [suggest_block_size()], [spatial_folds()], [evaluate_cv()]: spatial
#'   block cross-validation with a block size from the autocorrelation of
#'   the predictors.
#' * [fold_report()], [plot_folds()], [environmental_folds()],
#'   [spatial_sorting_bias()]: audit of cross-validation folds for leakage,
#'   spatial and environmental separation and spatial sorting bias.
#' * [ensemble_models()], [plot_ensemble()]: ensemble of the selected
#'   models and the disagreement between them.
#' * [partition_uncertainty()], [plot_partition()]: which source (GCM,
#'   scenario, algorithm, run) dominates the variance of projections, and
#'   where.
#' * [coordinate_spread()], [plot_coordinate_spread()]: spread of the
#'   predictions due to coordinate error, from models fitted to jittered
#'   replicates.
#' * [ignorance_map()], [record_ignorance()], [taxon_complexity()],
#'   [plot_ignorance()]:
#'   biogeographical ignorance from recording effort or from the records
#'   (completeness, age, reliability, precision), and the suitability where
#'   the data are sparse.
#' * [range_metrics()], [plot_range_metrics()], [write_range_metrics()]:
#'   EOO, AOO, locations and the patch structure of the range under each
#'   threshold rule, and their change between projections.
#' * [niche_overlap()], [overlap_scenarios()], [plot_niche_overlap()],
#'   [prediction_overlap()]: niche overlap between groups of records, its
#'   tests and its sensitivity to methodological choices.
#' * [null_model_test()]: are the models better than expected from where
#'   records are collected?
#' * [null_occurrences()], [transfer_bias_test()], [environment_shift()]:
#'   is a projected change driven by the records or by the study design?
#' * [odmap_table()]: the ODMAP reporting protocol filled from the objects
#'   of your workflow.
#' * [evaluate_models()]: AUC, MAE and bias for one or many models, Pareto
#'   ranking and bootstrap robustness.
#' * [plot_pareto()]: publication-ready AUC-MAE plot with the Pareto front.
#' * [boyce_index()], [plot_boyce()]: the continuous Boyce index by moving
#'   window or smoothed P/E curve, with bootstrap interval.
#' * [response_check()], [plot_response_check()]: are the response curves
#'   ecologically plausible (shape, turning points, rising edges)?
#' * [calibration_curve()], [brier_decompose()], [plot_calibration()]:
#'   calibration of the predictions (reliability diagram, calibration slope
#'   and intercept, Brier score decomposition).
#' * [thresholds()], [threshold_agreement()], [plot_thresholds()]: the
#'   threshold of a prediction under several rules, and the cells on which
#'   they agree.
#' * [select_models()], [plot_selection()], [plot_models()]: selection by
#'   minimum values or Pareto front; any two or three metrics, zoom, static
#'   or interactive 3-D plots.
#' * [extract_predictions()]: model predictions at evaluation points.
#' * [pareto_rank()]: generic multi-criteria Pareto ranking.
#' * [simulate_enm_data()]: example data.
#'
#' Please cite Konowalik & Nosol (2021) \doi{10.1038/s41598-020-80062-1}
#' when using this package (see `citation("preENM")`).
#'
#' @keywords internal
"_PACKAGE"
