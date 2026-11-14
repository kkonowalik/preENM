# preENM 2.7.0

* `temporal_validation()`: **paired record bootstrap** (`n_boot`, `conf`,
  `seed`). The records are resampled and the same resampled records are
  scored for every model. Every metric gets a percentile interval, and the
  `compare` metric (default: the first `combine` group) gets the gap to the
  best model with a paired interval, `_p_best`, `_p_worse` and `_worse`
  (the gap exceeds `margin`). Replicates in `attr(x, "bootstrap")`.
* `select_models(max_metrics =)`: upper limits for "lower is better"
  metrics, e.g. `c(MAEbal = 0.3)`, so one rule can require discrimination
  and calibration. Shown in the printout and as dashed lines in plots.
* Documentation: the background strategies are presented as choices to
  compare, not with a recommended winner.

# preENM 2.6.0

Improvements to model selection, prompted by a comparison with the
*Dryas octopetala* algorithm selection:

* `temporal_validation()`: out-of-fold records are now **pooled across
  runs** (new `pool`, default: pooled when no record is scored by more than
  one run of a model). Before, the run values were averaged, so a few
  records in a small fold weighed as much as many records in a large one;
  with unevenly sized folds this lowered e.g. the deep fossil AUC by about
  0.04. `pool = FALSE` gives the previous values. Per-record scores are in
  `attr(x, "records")`.
* New built-in metric **`MAEbal`** (`mae_balanced()`): the MAE averaged
  over presences and background, independent of the number of background
  points. A user metric with a built-in name now replaces it (with a
  message) instead of being an error.
* `select_models()`: `rule = "both"` (minimum values and Pareto front);
  documentation on when the Pareto front is the right tool (best single
  models) and when minimum values are (ensembles); `plot = interactive()`
  draws the selection.
* **Zoom**: `plot_selection(zoom = TRUE)` and `plot_pareto(zoom = TRUE)` show
  the region of the best models (or given limits) and name the models
  outside it; `plot(sel)` / `plot(evaluation)` draw the full view and the
  zoom side by side. Selected models are drawn on top of the others, bars
  are lighter in the zoom, and a `"zero"` metric (bias) gets a reference
  line.
* New `plot_models()`: any two or three metrics of a selection, evaluation
  or plain table; static (`ggplot`) or interactive and 3-D (`plotly`,
  opens in the viewer or browser, saved as one HTML file with `file`).
* `evaluate_models()`: models with missing predictions are evaluated on
  their own points (also in the bootstrap); the printout lists models with
  fewer points. Rows with an `NA` in one model should not be dropped.

# preENM 2.5.0

Environmental thinning and background strategies, generalised from the
*Dryas octopetala* workflow and its background experiments:

* `thin_environment()`: thins presence records in a PCA (or standardised)
  environmental space. The bin width is not arbitrary: it is the median
  environmental uncertainty of a record (the range of each axis within the
  record's positional accuracy, `radius_km`). One record per occupied bin is
  kept, so density is evened out while every removed record has a kept
  neighbour within one bin. Always kept: protected records (e.g. genetic
  samples), the minimum and maximum of every layer, convex-hull vertices of
  every pair of axes, and a geographic quota per equal-area block, so the
  niche is not truncated environmentally or geographically. Templates
  re-use the same space for subsets. `plot_env_thinning()`,
  `write_env_thinning()`, `kept()`/`removed()` (now generics).
* `sample_background()` gains `method`: `"uniform"` (default, unchanged),
  `"envstrat"` (even coverage of environmental bins), `"effort"`
  (proportional to recording effort, a target-group-like background) and
  `"enveffort"` (environmental bins with effort-weighted cells inside each
  bin; the best combination for *Dryas*, reducing the sampling bias).
* `gbif_effort_sql()`, `download_gbif_effort()`, `gbif_effort_rasters()`:
  recording-effort layers (records, species, target group per cell) from one
  aggregated GBIF SQL download. Credentials are read from `~/.Renviron`
  (`GBIF_USER`, `GBIF_PWD`, `GBIF_EMAIL`), never from scripts.
* `download_neotoma()`, `neotoma_records()`, `calibrate_radiocarbon()`
  (IntCal20 via 'rintcal'), `trace21k_id()`, `fossil_site_bins()`,
  `site_bin_slices()`: dated fossil records with their own age windows,
  flags and time slices, ready for `temporal_validation()`.
* `temporal_validation()` with `thresholds` also gives the predicted area
  (`_area`) and `_lift` = sensitivity / area, so a model cannot score well
  by predicting everything as suitable.

# preENM 2.4.1

* Fix: `crop_rasters()` failed for every file with recent GDAL versions
  ("Failed to parse ' 2' as number"), because numbers passed to GDAL were
  padded with spaces. Affected 2.3.0 and 2.4.0 on Linux and Windows.

# preENM 2.4.0

Selection of predictor variables across periods, scenarios and regions,
generalised from the correlation workflow for *Dryas octopetala*:

* `sample_environment()`: area-proportional random cells from several data
  sets (present, past periods, future scenarios), cells common to all,
  anomaly ("delta") data sets built at those cells, regions from polygons or
  tiles, layers matched by name, rasters read block by block, and a check for
  variables stored x 10 in one data set.
* `select_variables()`: systematic search of every variable subset: a pair
  and VIF screen in every data set and region, optional required groups
  (e.g. temperature and moisture), and scores for novelty (MESS),
  gradient stability, agreement between data sources and VIF. Results: the
  Pareto front, mean rank, variable frequencies, reference sets, a pair
  table (rho per data set, driver, sign changes), consistency between data
  sources, and the union dendrogram.
* `best_sets()`, `plot_variable_tree()`, `plot_correlation()`,
  `plot_variable_sets()`, `write_variable_selection()` (CSV, Newick tree,
  PNG).

# preENM 2.3.0

* `crop_rasters()`: batch cropping and masking of raster files (e.g. all
  CHELSA-TraCE21k time slices) to M or any polygon, file to file with GDAL
  (low memory), skipping finished files on re-runs, selecting variables by
  name, and continuing past failing files. Raw values, data type,
  scale/offset and NoData are kept exactly on the input grid; reading scaled
  integer rasters into R and writing them back in their integer type can
  silently round them and drop the scale and offset.

# preENM 2.2.0

Delimitation of the accessible area (M), generalised from the M-area
workflow for *Dryas octopetala*:

* `delimit_m()`: M from the records with a data-derived linking distance L
  (minimum spanning tree / single-linkage clustering). Methods `"delaunay"`
  (default: Delaunay triangles with edges up to L plus L/2 buffers, closed by
  L/2, holes filled), `"cluster_hull"` (convex or concave), `"convex_hull"`
  and `"buffer"`. Expert intervention: added polygons and sites, corridors up
  to a second linking distance L2, exclusions, extent, land trimming, manual
  L.
* `linking_distance()` and `plot_linking()`: the connectivity diagnostics.
* `m_scenarios()`: several M variants compared (area, records inside per
  group, overlap and Jaccard with a reference).
* `plot_m_area()`, `m_report()`, `m_class_shares()` (e.g. Koppen-Geiger
  shares, with `koppen_geiger_labels()`), `class_polygons()`,
  `write_m_area()` (GeoPackage, CSV, PNG) and `crop_to_m()`.
* The spanning tree is computed from the Delaunay triangulation (identical
  to single linkage on the full distance matrix, but with memory growing
  linearly), and buffers of many points are merged tile by tile, so tens of
  thousands of records are handled in seconds with modest memory.
* `merge_occurrences()` repairs text that is not valid UTF-8 (`encoding`).
* `sf` is now imported.

# preENM 2.1.0

Point filtering, one point per cell and background points, generalised from
a complete GBIF workflow for *Dryas octopetala*:

* `filter_occurrences()`: rule-based filtering of occurrence records,
  designed especially for GBIF downloads. 29 tests in a fixed order (each
  record keeps the first one it fails), each optional, with adjustable
  thresholds: duplicates, atlas grid cell centres, GBIF issues, basis of
  record, establishment means, coordinate uncertainty and decimals, garden
  and settlement text, country names as localities, record year, manual
  lists, gridded datasets (uncertainty blocks, lattices, stacked datasets),
  CoordinateCleaner tests, range polygon, distance to the shoreline, water,
  climate, built-up, cropland and GHS settlement layers, isolation and
  climate limits with a refuge rule.
* `filter_settings()` with the presets `"strict"` and `"lenient"`;
  `filter_layers()` for the optional layers; `filter_tests()`.
* Results: `kept()`, `removed()`, `filter_report()`, `filter_sensitivity()`,
  `record_fate()`, `write_filter_results()`.
* `read_gbif()` reads GBIF downloads (`.csv` or `.zip`) keeping coordinates
  as written.
* `merge_occurrences()`, `thin_by_cell()` (one record per raster cell, with
  source priority, centred in the cell) and `sample_background()` (random
  background sets from cells valid in all layers).
* Optional download helpers: `download_shoreline()` (GSHHG),
  `download_landcover()` (via geodata) and `aggregate_ghs()` (GHS-BUILT-S).
* `dbscan` is now imported; CoordinateCleaner, countrycode, geodata, readr,
  rnaturalearth and stringi are suggested.

# preENM 2.0.0

* The package is renamed from `ENMselect` to **`preENM`** and the repository
  moved to <https://github.com/kkonowalik/preENM>. The package will grow into a
  toolkit for data preparation for ENM/SDM (point filtering, background point
  selection, environmental thinning, correlation analysis); the existing
  evaluation and model-selection functions are unchanged. Replace
  `library(ENMselect)` with `library(preENM)`.
* New hex logo.
* Horizontal error bars use `geom_errorbar(orientation = "y")` instead of
  `geom_errorbarh()`, which is deprecated in ggplot2 4.0.
* `Metrics` is declared in Suggests (used only in a test that checks AUC, MAE
  and bias against it).

# ENMselect 1.2.0

* Choose which metrics to compare: `evaluate_models(criteria = ...)` accepts
  any metrics with a direction, `"max"`, `"min"` or the new `"zero"` (closer
  to 0 is better, e.g. `bias`). The default is still AUC vs MAE. The Pareto
  front, compromise distance, bootstrap intervals and `p_pareto` all follow
  the chosen criteria; `tol` is also available.
* Add your own metrics with `metrics = list(name = function(obs, pred) ...)`
  in `evaluate_models()` and `evaluate_cv()`.
* `plot_pareto()` gains `x_metric` / `y_metric`, defaults to the chosen
  criteria and orients axes so that better models are top-right (the MAE axis
  is now reversed, matching `plot_selection()`).
* `select_models()` and `pareto_rank()` support "closer to 0 is better"
  criteria (`"zero"` / `zero =`).
* `evaluate_models(weights = )` now defaults to equal weights for whichever
  criteria are used.

# ENMselect 1.1.0

New tools generalised from a full biomod2 workflow (spatial cross-validation,
fossil validation and selection of algorithms for an ensemble):

* `spatial_folds()`: spatial block cross-validation folds, balanced for
  presences and background points.
* `evaluate_cv()`: hold-out metrics per cross-validation run, summarised per
  model across folds and background sets, with diagnostics (SD between sets
  vs between folds, AUC per fold) and optional fixed thresholds for TSS.
* `temporal_validation()`: validation of hindcasts against dated records
  (e.g. fossils) ranked against background cells of the same time slice, by
  period, for pooled periods, for analog climates only, and as sensitivity.
* `select_models()`: selection by conventional minimum values or by the
  Pareto front (with tolerance for ties), mean rank, manual override.
* `plot_selection()`: selection figure for any pair of metrics.
* New metrics `tss_score()`, `boyce_index()` (identical to
  `ecospat::ecospat.boyce()`), `threshold_p10()`.
* Helpers `anomaly_climate()` (delta method for past climate),
  `outside_range()` (non-analog climate), `age_to_period()`.
* `pareto_rank()` gains `tol`.
* `evaluate_models()` now also reports TSS and the Boyce index.
* `simulate_enm_data()` also returns environmental layers.

# ENMselect 1.0.0

* First release: `evaluate_models()`, `pareto_rank()`, `plot_pareto()`,
  `extract_predictions()`, AUC/MAE/bias, bootstrap robustness of selection.
