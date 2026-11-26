# preENM 2.38.0

Bias-correction decision support (Lamboley & Fourcade 2024, *Journal of
Biogeography* 51, 1783-1794; Ten Caten & Dallas 2023, *Ecosphere* 14,
e4703; Inman et al. 2021, *Ecosphere* 12, e03422; Baker et al. 2024,
*Diversity and Distributions* 30, e13802):

* `bias_sensitivity()`: applies a grid of corrections (no thinning, one
  record per cell, spatial filtering at several distances, environmental
  thinning at several bin widths; background strategies and sizes) and
  reports how each changes the data, without fitting a model: records kept,
  coverage of the environmental bins of M and of the records' own bins,
  the coverage of random subsets of the same size (is thinning better than
  random removal?), Clark-Evans ratio, range kept and gaps to the limits of
  M (as `niche_truncation()`), shift from the uncorrected records, and for
  backgrounds their coverage, shift from uniform and overlap with the
  records. The printout names the choices that change the design most and
  least, and the thinning designs that leave the niche intact.
* `plot_bias_sensitivity()`: the designs on any two of these measures.
* Internal: the environmental binning of `background_sufficiency()` is
  shared with `bias_sensitivity()`.

# preENM 2.37.0

Target-group background, background sufficiency and adaptive M (Barber et
al. 2021, *Diversity and Distributions* 28, 128-141; Balej et al. 2025,
*Ecography* e08202; Whitford et al. 2024, *Ecological Modelling* 488,
110604; Rausell-Moreno et al. 2025, *Ecological Modelling* 507, 111177;
Barbosa 2026, *Biodiversity Informatics* 20, 81-88):

* `sample_background(method = "target_group")`: background weighted by the
  records of a target group (`target`: records or a raster of counts),
  per cell (`target_bandwidth_km = 0`, the classic target-group
  background) or smoothed by a Gaussian kernel (TGOB+ of Balej et al.);
  with `target_observer` and `presences`, each observer's normalised
  density weighted by their share of the focal species' records (PW-OOA);
  `target_stratify = TRUE` combines it with the environmental bins of
  `"enveffort"`. Densities are computed only at a random pool of cells, so
  global layers stay within memory.
* `gbif_target_group_sql()`: the GBIF SQL download of a genus, family or
  other group without the focal species, record by record (with
  `recordedBy`) or counted per cell; run with `download_gbif_effort()`.
* `background_sufficiency()`: draws the background (any method) at an
  increasing series of `n` and returns the saturation curve of the
  coverage of the environmental space of M (bins and area) and of the
  largest Kolmogorov-Smirnov distance of the marginal distributions; the
  printout recommends the `n` at which coverage saturates and says whether
  the curve saturated within the range tested.
  `plot_background_sufficiency()`.
* `delimit_m(method = "adaptive")`: spatially adaptive buffers for unevenly
  surveyed areas (Barbosa's cluster mean distance), with `adaptive_mult`,
  `adaptive_weight` and `adaptive_min_km`.

# preENM 2.36.0

Preparing data for integrated species distribution models (Zulian et al.
2021, *Diversity and Distributions* 27, 2498-2509; Morera-Pujol et al.
2023, *Ecography* e06451; Mostert et al. 2025, *Ecology and Evolution* 15,
e71029; Gueguen et al. 2026, *Ecography* e08588; Makinen et al. 2026,
*Diversity and Distributions* e70157):

* `merge_occurrences()` gains `data_type` (per source: `"presence_only"`,
  `"presence_absence"`, `"abundance"` or a label of your own) and
  `response`; the new columns `data_type` and `response` are read from
  `occurrenceStatus`, `present`, `individualCount` or a named column.
  Presence-only records with status `ABSENT` are reported.
* `filter_occurrences()` accepts a list of settings named by data type (and
  `default`), filtering each type with its own tests and thresholds.
  Absences are exempt from the presence-plausibility tests (range,
  isolation, climate limits) and are not duplicates of presences at the
  same place.
* `filter_report()` gives one column per data type, `record_fate()` the
  data type and response, the printout a summary per type, and
  `odmap_table()` the removals per type.
* `export_integrated()`: one dataset per data type (or source) in the input
  structure of intSDM (`$addStructured()` arguments), PointedSDMs (`sf`
  datasets with `present` / `counts`, and the matching `startISDM()`
  arguments) or biomod2 (`BIOMOD_FormatingData()` arguments with
  `data.type`), without depending on those packages; warns about
  presence-absence data without absences, extreme prevalence and small
  datasets. `write_integrated()` writes them to CSV.

# preENM 2.35.0

Propagating coordinate uncertainty (Gabor et al. 2023, *Ecography*
e06358; Wu et al. 2026, *Environmental Modelling & Software* 204, 107096):

* `jitter_uncertainty()`: `n` replicate datasets in which every record is
  moved to a random point within its own uncertainty circle (uniform over
  its area), optionally only inside M, on land and in cells with values in
  every layer; records with an uncertainty of 0 do not move, records
  without one stay fixed or get `na_uncertainty`. Reports how much the
  environment at each record changes across the replicates, relative to
  the niche breadth (the scale of `uncertainty_environment()`), and says
  whether propagating is worth the extra model runs. `propagate =
  "variable"` moves only records whose environment varies within the
  radius. With `plot_jitter()` and `write_jitter()`.
* `coordinate_spread()`: from one prediction per replicate, the per-cell
  mean, standard deviation and coefficient of variation due to coordinate
  error alone, computed as the disagreement layer of `ensemble_models()`
  and optionally compared with it; `plot_coordinate_spread()`.
* `uncertainty_environment()` adds `unc_action`: `"keep"` (uniform
  environment), `"propagate"` (it varies) or `"drop"` (radius above
  `max_radius_km`).

# preENM 2.34.0

Taxonomic name matching (Fuehrding-Potschkat et al. 2022, *Ecology and
Evolution* 12, e9168; Van Eupen et al. 2021, *Ecological Modelling* 444,
109453; Lipinski 2026, *Conservation* 6, 40):

* `match_names()`: resolves the scientific names of occurrence records
  against GBIF (rgbif), a ChecklistBank dataset (by default the Catalogue
  of Life) or a local Darwin Core table (e.g. the World Flora Online
  backbone, without internet). Adds the matched and accepted name, the
  accepted species, rank, status, match type (exact, fuzzy, higher rank,
  none), confidence and flags for synonyms, hybrids, aggregates,
  higher-rank and fuzzy matches. Each distinct name is looked up once;
  results are cached for the session and optionally in a file.
* Three name tests run first in `filter_occurrences()`: `name_unmatched`
  (not resolved, or only above the genus; optionally fuzzy matches and low
  confidence), `name_synonym_conflict` (resolves to an accepted name other
  than `target_taxon`; hybrid names) and `rank_conflict` (identified only
  to a coarser rank than `target_rank`, an infraspecific taxon out of
  scope, or optionally an aggregate). New settings: `target_taxon`,
  `target_rank`, `infraspecific_include`, `infraspecific_exclude`,
  `aggregates_in_scope`, `hybrids_in_scope`, `accept_fuzzy`,
  `min_name_confidence`. Without the columns of `match_names()` the tests
  are skipped and reported as such.
* `record_fate()` shows the queried and accepted name, match type and rank.
* jsonlite added to Suggests (for ChecklistBank).

# preENM 2.33.0

* `niche_overlap()`: new `truncate` option computes D, I and breadth from
  densities set to 0 outside the occupied envelope (`marginal_records`), as
  ecospat does with `th.sp`; the tests use the same choice. It is off by
  default: in simulations with known overlap it lowered D in every case
  (true D 0.5 estimated as 0.40 instead of 0.48 with 200 records per
  group) and only brought fully separate niches from 0.007 to exactly 0.
* `plot_niche_overlap()`: `limits = "available"` (default, all available
  environments), `"niches"` (zoom to the occupied envelopes of the pair)
  or `c(xmin, xmax, ymin, ymax)`, with a `margin` around the automatic
  windows.

# preENM 2.32.0

Grain sensitivity (Henderson et al. 2023, *Forest Ecology and Management*
529, 120741; Brambilla et al. 2023, *Ecology* 105, e4224; Ferrer-Ferrando
et al. 2025, *Ecography* e07766):

* `grain_sensitivity()`: aggregates the predictor layers to a series of
  coarser grains, on disk with GDAL and one grain at a time, and repeats
  cell thinning, environmental thinning, the coordinate-uncertainty check
  and the variable selection at each. It reports the records kept, the
  share of records too uncertain for the grain, the collinear pairs and the
  change of Spearman rho, the selected variable set and its Jaccard
  agreement with the finest grain, and says up to which grain the selected
  set stays the same. Ties for the top rank are flagged, with the agreement
  of all passing sets as the safer measure there.
* `plot_grain_sensitivity()` and `write_grain_sensitivity()`.

# preENM 2.31.1

* `plot_niche_overlap()` uses the ColorBrewer diverging palette RdYlBu: the
  niche of the first group only in blue, of the second only in red, the
  shared niche in yellow, and the legend in that order.

# preENM 2.31.0

Niche overlap between groups (Broennimann et al. 2012, *Global Ecology and
Biogeography* 21, 481-497; Warren et al. 2008, *Evolution* 62, 2868-2883;
Warren et al. 2021, *Ecography* 44, 504-511; Rodder and Engler 2011, *Global
Ecology and Biogeography* 20, 915-927; Guisan et al. 2014, *Trends in
Ecology & Evolution* 29, 260-269; Lo Parrino et al. 2023, *Ecography*
e06432):

* `niche_overlap()`: two or more groups of records (a grouping column, e.g.
  lineages or cytotypes, or a named list) in a shared PCA (or raw)
  environmental space built from their pooled backgrounds; kernel densities
  of records and environments; Schoener's D and the corrected Hellinger I,
  graded with the classes of Rodder and Engler; Levins' breadth;
  expansion, stability and unfilling; equivalency and similarity tests (in
  both directions) with null distributions, p-values and standardised
  effect sizes; with `predictions`, D, I and rank correlation of model
  predictions in geographic space.
* `overlap_scenarios()`: the same over a grid of choices (background,
  space, thinning, grid size), saying whether each verdict holds in every
  scenario or which choice changes it.
* `plot_niche_overlap()`: the niches of a pair in the shared space with
  the available environments; `prediction_overlap()`: overlap of
  suitability rasters (ENMTools).

Ignorance maps (requested):

* `taxon_complexity()`: the number of synonyms and of subspecies and
  varieties of each species in World Flora Online (WorldFlora package, now
  in Suggests) or another Darwin Core backbone, as a reliability weight for
  `record_ignorance()`; a rule of thumb, not a published weighting.
* `record_ignorance(env = )`: empty cells interpolated with the
  spatio-environmental distance of Tessarolo et al. (2021); cells with
  records of the focal species carry half of the weight, as in their
  method.

# preENM 2.30.0

Range metrics for conservation assessment (Curd et al. 2023, *Global
Change Biology* 29, 631-647; Dauby et al. 2017, *Ecology and Evolution* 7,
11292-11303; Omar and Elgamal 2021, *Kew Bulletin* 76, 477-496; Auffret et
al. 2024, *Diversity and Distributions* 30, e13834):

* `range_metrics()`: from binary maps, or continuous maps with the result
  of `thresholds()`, under every threshold rule: suitable area, extent of
  occurrence (minimum convex polygon), area of occupancy on the IUCN 2 x 2
  km grid and the number of locations on a grid of the size of one threat
  (both in an equal-area projection, with the grid placed to give the
  fewest cells, as in ConR), and the landscape structure of the range:
  number and sizes of patches, largest patch index, edge density, distance
  to the nearest patch and mean suitability per patch. With several
  projections, the change of every metric from the baseline. The printout
  says when the threshold rule changes the Red List threshold met (the
  category would then depend on the rule), notes that projected changes
  are uncertain in both directions (Auffret et al. 2024) and that the
  suitable area is larger than the occupied area; `occupied = TRUE` keeps
  only the patches with records.
* `plot_range_metrics()`: a metric under each rule with the criterion-B
  thresholds; `write_range_metrics()` writes the tables and the plot.

Requested additions:

* The scale check (0-1, 0-100, 0-1000) now also covers `thresholds()`,
  `null_model_test()` (pool divided like the data), `temporal_validation()`
  and `transfer_bias_test()` (thresholds given on the model's scale divided
  alike).
* `record_ignorance()`: maps of biogeographical ignorance from occurrence
  records (e.g. a GBIF download of the reference group), with the
  components of Tessarolo et al. (2021): completeness of each cell's
  inventory (from the species accumulation curve), age of the records,
  optional taxonomic reliability weights and the precision of the
  coordinates; the best record of the focal species per cell, and
  interpolation of empty cells.
* `response_check()` and the manual explain which response shapes are
  plausible for direct, resource and indirect gradients and for a niche
  covered on one side only, and that a curve should follow how much more
  often the species occurs at a value than the value is available, not the
  distribution of the raster.

# preENM 2.29.0

Plausibility of response curves (Fiorentino et al. 2025, *Ecological
Modelling* 508, 111207; Chollet Ramampiandra et al. 2023, *Ecological
Modelling* 481, 110353; Allyn et al. 2025, *Diversity and Distributions* 31,
e70036; Yau et al. 2026, *Ecological Informatics* 99, 104008; Broussin et
al. 2024, *Ecological Modelling* 498, 110865):

* `response_check()`: for the response curves exported by any modelling
  software (or binned from raw predictions), per variable, model and run:
  the shape (increasing, decreasing, unimodal, U-shaped, irregular, flat),
  the number of turning points, a wiggliness measure, the optimum and its
  position in the calibrated range, and whether the curve is still rising
  at an edge of that range. With records (and background) it says whether
  each curve follows the records along the gradient; with a projection, how
  much of it lies beyond a rising edge. For binned curves, reversals within
  the sampling noise of the bins are not counted. The printout names the
  implausible curves first.
* `plot_response_check()`: the curves with their shapes, the calibrated
  range and the density of records beneath; `file` saves it.

# preENM 2.28.0

Projection uncertainty and maps of ignorance (Wilson et al. 2024,
*Ecological Informatics* 80, 102537; Tessarolo et al. 2021, *Ecography* 44,
1743-1755):

* `partition_uncertainty()`: splits, in every cell, the variance of a set of
  projections among factors such as GCM, scenario, algorithm and run (an
  analysis of variance with the factors as main effects; interactions and
  residual as the rest). Returns a share map per factor, the total variance
  and the dominant source of each cell, and the shares pooled over the area
  and an optional `region`, each cell weighted by its variance. Warns when
  the design is unbalanced (shares then depend on the order of the
  factors) and notes that the share of scenarios depends on which
  scenarios were chosen. Rasters are processed block by block on disk.
* `plot_partition()`: the dominant-source map next to a stacked bar of the
  pooled shares; `file` saves it.
* `ignorance_map()`: turns a recording-effort raster (e.g. from
  `gbif_effort_rasters()`) into a score of biogeographical ignorance from 0
  (well recorded) to 1 (nothing recorded), with the half-ignorance,
  log-normalised or normalised algorithms of Ruete (2015), optionally from
  records per species and with effort counting within a radius. With a
  suitability raster it adds the suitability where ignorance is low and a
  joint map (suitable or not, known or ignorant). `plot_ignorance()` draws
  both.

Prediction scales and ensembles:

* Algorithms write predictions on different scales (0-1, 0-100, 0-1000).
  `extract_predictions()`, `evaluate_models()`, `evaluate_cv()`,
  `calibration_curve()`, `ensemble_models()` and `partition_uncertainty()`
  now recognise the scale of each model from its largest value and, with
  `rescale = TRUE` (default), divide models on 0-100 or 0-1000 by 100 or
  1000, with a message naming them; `rescale = FALSE` keeps the values and
  warns. Unrecognised scales (negative values or above 1000, e.g. raw
  Maxent output) are left as they are, with a warning. Before, only a
  pooled warning about values outside 0-1 was given.
* `ensemble_models()` now uses PCAm by default (the median of the models
  closest to their common trend; Zhu et al. 2021); with fewer than three
  models it falls back to the median, with a message. **Results of
  `ensemble_models()` with the default method change.**
* The figures of the manual (README) are no longer built into the package,
  which kept the installed size above 5 MB; they stay on GitHub.

# preENM 2.27.1

The manual is reordered to follow the workflow:

* "How the selection works" is now section 18, before the full workflow
  (now section 19), which uses it.
* The subsections of the full workflow follow its steps: auditing the
  folds, fossil (temporal) validation, plotting, ensembles, and thresholds
  (moved from section 14, as a binary map is usually made from the final
  ensemble; a new step 8 says so).
* The overview table follows the order of the manual, and every internal
  link is checked (one link to "Choosing which metrics to compare" was
  broken).
* A test of `ensemble_models()` compared file paths literally and failed on
  macOS, where temporary files lie under `/private/var`; it now compares
  normalised paths.

# preENM 2.27.0

Ensembles and the disagreement between models (Zhu et al. 2021,
*Ecological Modelling* 448, 109502):

* `ensemble_models()`: combines the prediction rasters of the selected
  models (`select_models()`) or of the Pareto front (`evaluate_models()`);
  dominated or unselected models are left out. Methods: median (default),
  unweighted mean, PCAm (the median of the half of the models closest to
  the first principal component of the predictions, recommended by Zhu et
  al.) and the mean weighted by the compromise score of the evaluation.
  The help page explains why the weighted mean is not the default. Returns
  the ensemble with the standard deviation and coefficient of variation
  between models per cell; the printout says how much and where the models
  disagree (by level of the ensemble prediction, the area above levels of
  disagreement, the model that departs most). Rasters are read from files
  or folders and combined block by block; large stacks keep every
  intermediate layer on disk.
* `plot_ensemble()`: the ensemble and the disagreement side by side;
  `file` saves it.

# preENM 2.26.1

Completing the audit against Table 1 of Stock et al. (2023):

* `fold_report()` gains `time`: the share of training records dated after
  the earliest test record, per fold. For a model to be projected forward
  in time, folds should be ordered in time; for a model of one period it
  does not matter, so it is reported, not called leakage.
  `evaluate_cv(time = )` passes a column of dates or years to the audit.
* `fold_report(env = )` names predictors that look like coordinates or
  time (`lon`, `lat`, `x`, `y`, `year`, `doy`, ...), which let a model
  learn where and when records were made rather than the response to the
  environment. Only the names are checked.
* `variable_relevance()` and the manual say that, when models are
  cross-validated, the test should be run on the training records of each
  fold: selecting predictors on all records is feature selection on the
  response before cross-validation, a form of data leakage.
* The year of Yoo and Koo is 2025 (*Cartography and Geographic Information
  Science* 52, 596-609), not 2024.

# preENM 2.26.0

Auditing cross-validation folds (Stock et al. 2023, *Nature Ecology &
Evolution* 7, 1743-1745; Yoo and Koo 2025; Koldasbayeva and Zaytsev 2025;
Hijmans 2012, *Ecology* 93, 679-688):

* `fold_report()`: checks a fold assignment for leakage (records with the
  same coordinates, or in the same raster cell, in more than one fold) and
  reports per fold the numbers of training and testing records, the
  distance from test to nearest training records against the block size
  and the spacing of the records, the overlap of the predictor
  distributions between training and testing, the share of test records
  outside the training range, and the spatial sorting bias. The printout
  says plainly whether leakage was found and in which folds.
* `evaluate_cv()` gains `coords` (and `env`, `block_size`, `audit`): the
  folds are audited alongside the metrics, the audit is printed with them,
  and a warning is given when leakage is found. The metrics themselves do
  not change.
* `environmental_folds()`: folds by k-means clustering in the environmental
  space of `thin_environment()`, as a complement to `spatial_folds()`.
* `spatial_sorting_bias()`: Hijmans' SSB (mean distance from test presences
  to the nearest training presence divided by that from test background)
  and a background sample matched by pairwise distance, so that a corrected
  evaluation can be run.
* `plot_folds()`: map of the folds with the leaking records marked; `file`
  saves it.
* Distances are great-circle distances for longitude/latitude, computed by
  nearest-neighbour search (fast for tens of thousands of records).
* Internal: the PCA of `thin_environment()` is now a shared helper; its
  results do not change.

# preENM 2.25.0

Smoothed Boyce index (Liu et al. 2024, *Ecography* 2025, e07218):

* `boyce_index()` gains `smooth`: besides the moving window (`"none"`, the
  default, unchanged, so existing results do not change), the P/E curve can
  be estimated by a penalised logistic regression spline (`"spline"`,
  without extra packages), the thin plate spline of Liu et al. (`"gam"`,
  needs mgcv, now in Suggests) or loess (`"loess"`). With few presences the
  smoothed index is much less variable and not underestimated.
* `n_boot` gives a bootstrap confidence interval of the index (presences
  and background resampled separately) and, for smoothed curves, a band
  around the P/E curve. With `n_boot > 0` or `details = TRUE` the result is
  a `boyce_index` object with the P/E curve, the interval and a test of a
  flat curve; the printout says in plain words how the frequency of
  presences changes with the prediction, and warns when a smoothed curve is
  not distinguishable from flat (its index can then be anywhere between -1
  and 1) and when the moving window has fewer than 50 presences.
* `plot_boyce()`: the P/E curve with its bootstrap band; `file` saves it.
* `Boyce_smooth`: optional metric for `evaluate_models()`,
  `evaluate_cv()` and `null_model_test()` (`metrics = "Boyce_smooth"` or a
  criterion), next to the moving-window `Boyce`.
* Why the default stays: in simulated random models with 50 presences the
  smoothed index was 0.9 or more in a third of the runs (the moving window
  never), and good models often score exactly 1 and tie.

# preENM 2.24.0

Calibration diagnostics (Silva Filho et al. 2023, *Machine Learning* 112,
3211-3260), completing the evaluation with discrimination, calibration and
ranking:

* `calibration_curve()`: the reliability diagram of one model or of every
  model in an `extract_predictions()` table, with bins of equal width or
  equal count, the observed frequency against the mean prediction of each
  bin and binomial (Clopper-Pearson) intervals; plus the calibration slope
  and intercept with confidence intervals. The printout says in plain words
  whether a model is over- or under-confident and whether its predictions
  are too high or too low on average, and says when the slope is not
  different from 0, so calibration cannot be judged.
* `brier_decompose()`: the Brier score split into reliability, resolution
  and uncertainty, with the two within-bin terms of Stephenson et al.
  (2008) that make the sum exact. `brier_score()`: the Brier score.
* `plot_calibration()`: the reliability diagram with the 1:1 line, one
  panel per model; `file` saves it.
* `evaluate_models()`, `evaluate_cv()` and `null_model_test()` accept
  `calibration_slope`, `calibration_intercept` and `brier` as optional
  metrics, computed when named in `criteria` or in `metrics` (e.g.
  `metrics = "brier"`), so they can be Pareto criteria.
* New criterion direction `"one"` (closer to 1 is better), for the
  calibration slope; also in `select_models()`, `pareto_rank(one = )` and
  the plots (dotted reference line at 1).
* The help pages note that with background instead of absences these
  values compare models evaluated on the same background only
  (Jimenez-Valverde 2025, 2026).

# preENM 2.23.0

Thresholds under several rules (Hellegers et al. 2025, *Diversity and
Distributions* 31, e70019; Yu et al. 2024, *Ecology and Evolution* 14,
e11208):

* `thresholds()`: the threshold of a continuous prediction under nine
  rules (maxSSS, maxTSS, equal sensitivity and specificity, maxMCC, maxF,
  prevalence, fixed sensitivity, p10 and the mutation point of the P/E
  curve found by threshold regression, BTQR), with sensitivity,
  specificity, TSS, MCC, F and the predicted-suitable share of each. Takes two vectors or the table of
  `extract_predictions()` for several models at once; `rules` restricts the
  set. The printout says how much the mapped range varies between the most
  and the least inclusive rule.
* `threshold_agreement()`: for a prediction raster, the number of rules
  that call each cell suitable, and the area suitable under all rules,
  under most and under at least one.
* `plot_thresholds()`: the predicted-suitable share under each rule, or the
  map of the agreement.

# preENM 2.22.1

* `plot_niche_coverage()`: the subtitle of the map is shorter and wrapped,
  so it fits in the panel; shares below 1% in the legends get one decimal
  (e.g. 245 of 49501 records: 0.5%, not 0%).

# preENM 2.22.0

`niche_coverage()`:

* **Much less memory and time** with large rasters: values are read only at
  the records and at random locations in the region, and whether a record
  lies in the region is a point-in-polygon test. Before, the region was
  rasterized onto the whole grid and every valid cell was listed, which
  could fail with `std::bad_alloc` on global 30 arc-second layers. In a
  test with 200 million cells: peak memory 9.3 GB -> 0.6 GB, 67 s -> 2 s,
  identical result.
* **Graded advice**: small, moderate, large and very large gap
  (`limits = c(0.1, 0.25, 0.5)`; two limits as before are still accepted).
  Only the last limit is from Goicolea et al. (2025); the printout says
  so.
* **Automatic grid**: `n_bins = "auto"` (20 x 20, coarser with fewer than
  about 1,200 records). The unfilling is also given for a coarser and a
  finer grid, and the printout sums up the range over grids and
  `min_records` ("19% to 36%: from a moderate to a large gap").
* The note on what a nested model adds now explains its main use, regional
  projections (future conditions of the region that the species occupies
  elsewhere), and its assumption (no strong local adaptation).
* Figure for papers: no advice in the subtitle; legends with counts and
  percentages; the region in ColorBrewer green on grey land, its outline
  drawn above the records.

# preENM 2.21.0

`niche_coverage()` extended after its first use on a real circumpolar
species:

* **How serious is the gap?** A bootstrap 95% interval of `unfilling`
  (`n_boot`) and a null comparison: the unfilling that a random sample of
  the whole-range records of the size of the regional records would leave
  (`n_null`), with SES and p-value. A large SES means the region lacks
  conditions of the niche; a small one that the gap is about what the
  number of regional records alone would leave.
* **Robustness:** the unfilling for `min_records` = 1, 2, 3 and 5, with a
  note when the advice changes (sparse, remote records can decide it).
* `unfilled_records()`: the whole-range records in unfilled niche cells,
  with coordinates, to map and check them.
* `classify_points()` (and `predict()`): places new points (the region
  under another climate, evaluation or fossil sites) in the fitted niche
  space and grid, in seven categories (filled; occupied elsewhere with the
  conditions available in the region or absent from it; region only;
  available but unoccupied; unoccupied inside the space; outside it), with
  the regional record density of each point's cell. With `predicted` it
  counts the misses per category, i.e. whether a nested model could help.
  Unit guard: warns when points land far beyond the whole range in the
  niche space or a variable's median is far from that of the records
  (e.g. Kelvin against degrees Celsius).
* `log_vars`: variables log-transformed before the PCA (`"auto"` for
  bioclim precipitation amounts).
* Per niche cell (`bins`): records, regional records, regional density
  percentile and category; the advice text mentions the sensitivity to
  sparse records and when a nested model helps.
* `plot_niche_coverage()` draws the niche space **and a map** of the
  whole-range records coloured by their niche cell, with the region
  outlined (`type = "both"`, `"pca"` or `"map"`).

# preENM 2.20.0

* `niche_coverage()`: **`unfilling` and `not_available` now count cells of
  environmental space, not records.** Each cell occupied by the species
  counts once, however many records it holds. Before, they were shares of
  records, which can never exceed the share of records outside the region:
  with a well-sampled region (often 90% or more of the records) a large
  part of the niche known from few records elsewhere hardly counted, and
  the advice could wrongly say that a regional model is enough. The
  record-based shares are still given (`unfilling_records`,
  `not_available_records`) for comparison. New `min_records` (records a
  cell needs to count) and `bins` (one row per niche cell).
* `plot_niche_coverage()` frames the niche cells without regional records
  in red, with their number in the legend; regional records are drawn
  smaller so the records below stay visible.
* The printout explains why the per-variable shares are smaller than those
  of the niche (single variables against combinations of conditions).

# preENM 2.19.0

* **Every plotting function can save its figure**: `file = "fig.pdf"`
  (PNG, JPEG, TIFF, SVG by extension), several files at once
  (`c("fig.png", "fig.pdf")`), or a name without extension for PNG and
  PDF; `width` and `height` in inches. Added to `plot_pareto()`,
  `plot_selection()`, the `plot()` methods of selections and evaluations
  (both panels), `plot_models()` (static plots too), `plot_env_thinning()`,
  `plot_m_area()`, `plot_linking()`, `plot_niche_truncation()`,
  `plot_variable_tree()`, `plot_correlation()`, `plot_transfer_bias()`,
  `plot_niche_coverage()`; `plot_variable_sets()` accepts several files.
* `niche_coverage()`:
  * reports records that lie outside the layers (they were left out
    silently) and warns when all remaining records lie in the region: the
    region then trivially contains their niche, and the advice is now
    "cannot be judged" instead of "regional model sufficient";
  * `region_points`: a separate regional data set (e.g. a national
    survey), compared with the records of the whole range and drawn as
    crosses;
  * `file` saves the plot directly;
  * clearer printout (plain sentences, definitions at the end) and plot
    (legend with the numbers of records, "whole-range records inside /
    outside the region").

# preENM 2.18.0

* `niche_coverage()` and `plot_niche_coverage()`: for a widespread species
  modelled in one region (e.g. a country), how much of its niche the region
  contains: `unfilling` (share of the whole-range niche without regional
  records, as in Goicolea et al. 2025, *Ecography* e07328), the share in
  conditions the region lacks, per-variable coverage, and advice on whether
  a spatially nested (hierarchical) model is needed. Manual section 10
  explains nested models and when to use them.
* `plot_variable_sets()`: `sizes` shows only some set sizes (e.g.
  `sizes = 2:3`), `n = "all"` shows every passing set, and `file` saves the
  plot (PDF, PNG, SVG, ... by extension; height from the number of rows).
  The legend and subtitle of the score card now say what the shading means:
  darker = better, i.e. the lowest value where lower is better and the
  highest where higher is better.
* `write_variable_selection(format = )`: figures as PNG, PDF or both.

# preENM 2.17.0

How large is a design bias? `transfer_bias_test()` now measures how far
the change predicted from the real records is from that of random points,
and how much of it the design produces:

* `ses`: standardised effect size, (observed - null mean) / null SD. Not
  limited by the number of null models, unlike the p-value, so "just
  outside" and "far outside" the null interval can be told apart, and
  comparable between projections and metrics.
* `design_share`: the part of the predicted change that random points
  reproduce (1 = all of it; above 1 they predict more change than the
  records; below 0 the opposite direction).
* `prob_exceed`: probability that a model of the real records predicts a
  change beyond that of a model of random points. `observed` can now be a
  list of models of the real records (e.g. cross-validation runs), whose
  variation then enters this probability.
* A fourth verdict, **design-dominated**: outside the null interval, but
  random points reproduce at least `dominance` (default half) of the
  change, or more. **Changed:** such results were called "data-driven"
  before; a prediction whose change random points reproduce (or exceed)
  is not evidence about the species even when it differs from them.
* `plot_transfer_bias()` labels each projection with its SES and design
  share and colours the point by the verdict; the printed summary explains
  each projection in these terms.

# preENM 2.16.0

Design bias of projections, after Warren et al. (2021, *Methods in Ecology
and Evolution* 12, 985-995) and Warren et al. (2026, *Ecography* e08559):

* `transfer_bias_test()` asks whether a projected change (to another
  period, scenario or region) is driven by the records or by the study
  design. It compares the change predicted from the real records with the
  change predicted by the same models fitted to random points, for the
  mean change of suitability, the share of cells with falling suitability
  and the change of suitable area, for the whole area and an optional
  region. It returns a verdict per projection (design-driven, not
  distinguishable, data-driven), a summary in plain words and maps of the
  mean and spread of the change across the null models. The predictions
  can come from any software (rasters or files), or from your own fitting
  function. `plot_transfer_bias()` draws the result.
* `null_occurrences()` draws the random point sets: uniformly, in
  proportion to recording effort (`method = "bias"`, the remedy of Warren
  et al. 2026 for the inflated false positive rates of null models that
  ignore sampling bias) or by moving and turning the real point pattern
  (`method = "pattern"`), which keeps its clustering.
* `environment_shift()` is a check before modelling: per projection and
  variable, the share of present conditions that disappear and the share
  of projected conditions that are new (Williams et al. 2007), with the
  design bias to expect.

# preENM 2.15.0

Variable selection (`sample_environment()`, `select_variables()`):

* The printed selection leads with the result: the best sets of every
  size, then a summary in plain words ("IN SHORT"), then the details
  (clusters, filter counts, regions). `recommend_sets()` returns the
  summary: the suggested set and its scores, close alternatives, the sizes
  for which nothing passes and why, variables shared by the good sets or
  never passing, correlations that change between data sets, and the
  regions in which pairs are excluded only by the regional screen.
* `best_sets(vs, size = k)` no longer returns `NULL` when no set of that
  size passes: it returns a table without rows and says at which filter
  the sets were lost.
* Masks for single data sets: `sample_environment(exclude = , keep = )`
  leaves out cells (e.g. the ice sheet of a period) or restricts a data
  set (e.g. to the land of a period), with polygons or rasters in any CRS.
  `select_variables(skip_regions = )` leaves a region out of the regional
  screen, in all data sets or in one (`"LGM | Scandes"`). Data sets that
  look like past periods and have no mask give a warning, because
  palaeoclimate layers usually have values under the ice and over the sea
  (`options(preENM.past_warning = FALSE)` switches it off).
* A misplaced argument such as `bounds = list(fix_scale = TRUE)` is now
  reported instead of being ignored silently; `bounds` and `ratio_vars`
  without anomalies give a note.
* `plot_variable_sets()` draws a score card by default: the best sets as
  rows, the scores as columns, shaded by how good each value is among the
  passing sets. The scatter (`type = "scatter"`) labels the best sets with
  their rank and gives the key below the plot, and chooses the two scores
  with the most distinct values, so the points no longer fall on lines.
* `plot_variable_tree()` shows the best sets of each size as rows of dots
  below the tree and shades the clusters of collinear variables; several
  sets, e.g. your own next to the best, can be given.
* `write_variable_selection()` also writes `_best_sets.csv`,
  `_summary.txt` and `_sets_scatter.png`.

# preENM 2.14.0

* Five new day-length indices in `climate_indices()` (25 in total), all in
  the `light` group and computed from latitude and the temperature season,
  without extra data: `daylength_end` (photoperiod when the growing season
  ends), `daylength_gs` (mean day length of the season), `daylight_gs`
  (hours of daylight in the season), `daylength_max` (the longest day) and
  `days_to_solstice` (days from the start of the season to the summer
  solstice, with the start interpolated between monthly means).
* `climate_index_names(table = TRUE)` and the index table of the manual now
  give, for every index, its group, what it is useful for and its source
  in the literature; indices defined here are marked "this package". The
  references of `?climate_indices` are complete for all indices.
* `download_chelsa_trace21k()`, `trace21k_urls()`, `trace21k_time()`:
  CHELSA-TraCE21k (monthly layers and bioclimatic variables of any century
  of the last 22,000 years), for a region only if wanted, so the climate
  indices can be computed for the past. `monthly_files()` also finds months
  written without a leading zero.
* `select_variables()` reports the regions explicitly: the printed summary
  lists every region of every data set with its cells, strongest
  correlation and pair; `correlation_by_region()` gives the same table for
  all variables or for a chosen set; sets carry `regional_maxr` and
  `worst_region`, pairs `regional_driver`; regions too small to be screened
  are listed as such. `write_variable_selection()` writes both tables.
* `variable_relevance(null = "noise")`: the classic random variable
  without spatial structure, as an option. The default is unchanged.
* The manual (README) is reordered into the steps of a project: installation,
  an overview table with links, the steps before modelling, the steps
  after modelling. One index table instead of two, generic examples
  (three periods instead of names from one study), new figures.

# preENM 2.13.0

Are the selected predictors better than random variables? The elimination
procedure of Konowalik (2022, Scientific Reports 12, 16443), before
modelling and inside the models:

* `variable_relevance()`: a test that needs no model. For each layer, the
  Kolmogorov-Smirnov distance between the values at the records and in the
  calibration area is compared with the same distance for random layers
  with the same spatial structure: smooth random surfaces with the layer's
  autocorrelation range (`null = "field"`) or the layer itself shifted
  around the raster (`null = "shift"`, torus translation). White noise is
  not used as the null, because any smooth layer beats it at clustered
  records.
* `random_variables()`: random layers for any modelling software. Each has
  exactly the values and empty cells of a real layer, arranged in smooth
  random patches of the same size (`method = "field"`) or cell by cell
  (`method = "shuffle"`, the classic random variable).
* `compare_to_random()`: summary of variable importance from models with
  random variables, over algorithms and runs: the share of runs in which
  each predictor beats the random variables, and `keep`, which drops a
  predictor only when most algorithms do not support it.

# preENM 2.12.0

Six additions around the model, all working with any modelling software:

* `extrapolation_risk()`: maps where a projection (period, scenario or
  region) leaves the calibration conditions: ExDet novelty of single
  variables (`nt1`) and of their combinations (`nt2`, Mahalanobis
  distance), mobility-oriented parity (`mop`), the number of variables
  outside their range and the most novel variable.
* `null_model_test()`: significance of AUC, MAEbal or any metric against
  a null model in which the presences are random locations, drawn
  uniformly or in proportion to recording effort. With effort, a model
  that only predicts where people record is no longer significant.
* `suggest_block_size()`: a block size for `spatial_folds()` from the
  range of spatial autocorrelation of the predictors (empirical
  semivariograms, median range over the layers), without further
  packages.
* `check_sample_size()`: the number of records (as occupied cells)
  against the minimum numbers of van Proosdij et al. (2016), records per
  predictor and records per fold.
* `match_climate_period()`: labels records collected within, before or
  after the period of the climate data.
* `odmap_table()`: the ODMAP protocol (Zurell et al. 2020) with the
  entries on cleaning, clipping (M), thinning, background, predictors,
  collinearity, novelty, partitioning, performance and selection written
  from the objects of your workflow; your own entries through `fields`;
  the rest left empty to fill. Optionally written to CSV.

# preENM 2.11.0

Two new data-preparation checks and five more climate indices:

* `uncertainty_environment()`: judges coordinate uncertainty in
  environmental terms. For each record it measures how much every
  predictor varies among the cells that the uncertainty circle touches,
  relative to the breadth of the species' niche on that predictor, and
  marks as usable (`unc_ok`) the records whose true conditions cannot be
  far from the recorded ones. Imprecise records on flat terrain are kept;
  precise-looking ones on steep gradients can be dropped.
* `niche_truncation()` and `plot_niche_truncation()`: before modelling,
  reports for each predictor whether the records are still frequent where
  the conditions of the calibration area (M) end (a truncated response),
  how far the records stay from the limits of M, and which share of each
  projection (period, scenario or region) lies beyond M and beyond the
  records. Truncated ends that a projection goes beyond are flagged.
* `climate_indices()`: Kira's warmth and coldness indices
  (`kira_warmth`, `kira_coldness`), extraterrestrial radiation over the
  year and over the growing season (`ra`, `ra_gs`) and the day length at
  the start of the growing season (`daylength_start`).
* `climate_index_names()`: lists the twenty indices, by group (`water`,
  `heat`, `cold`, `season`, `light`), with `drop` to leave some out and
  `table = TRUE` for their meaning and units. `indices` in
  `climate_indices()` also accepts group names, e.g.
  `indices = c("cold", "gdd")`.
* **Changed default:** `climate_indices()` now computes all indices
  (before: the seven of version 2.9.0). For the previous result use
  `indices = climate_index_names()[1:7]`.

# preENM 2.10.0

New climate indices and faster computation (results unchanged):

* `climate_indices(indices = "all")` adds eight indices of limits shared
  by many organisms: growing degree days (`gdd`, base `gdd_base`), frost
  severity and duration (`cold_stress`, `frost_months`, below
  `frost_threshold`), `continentality`, the Bagnouls-Gaussen dry season
  (`dry_months`), the FAO length of the growing period (`lgp`), and the
  Emberger and Ellenberg quotients. The default remains the seven
  indices of 2.9.0.
* Faster, with identical results:
  * Pareto fronts (`evaluate_models()`, `select_models()`,
    `pareto_rank()`, `select_variables()` and their bootstraps):
    dominance is computed once in blocks instead of once per point and
    per front (about 10 times faster for 300 models on two criteria).
  * `boyce_index()`: windows counted by binary search on sorted values
    (about 14 times faster with 100,000 background points), which also
    speeds up bootstraps with `boot_metrics = "all"`.
  * `thin_environment()`: nearest kept record found with a k-d tree
    (`dbscan::kNN()`) instead of a full distance matrix, and bins named
    in one vectorised step (about 6 times faster for 30,000 records).
  * `select_variables()`: the novelty (MESS) of each candidate set is
    computed on the distinct patterns of novel variables instead of on
    every sampled cell (about 20 times faster for that step with
    thousands of candidate sets).
* Fixed: with `tol > 0`, three or more models could dominate each other
  in a circle (each better than the next by more than `tol` on one
  criterion), and the Pareto ranking then never finished. Such models
  now share a front.

# preENM 2.9.0

Derived climate variables, generalised from the C3/C4 workflow:

* `climate_indices()`: vapour pressure deficit, potential evapotranspiration
  (Hargreaves-Samani, FAO-56), climatic moisture deficit, aridity index,
  heat stress, summer rain index and heat-rain synchronicity from monthly
  Tmax, Tmin, (Tmean) and precipitation; any CRS (latitude per cell);
  written to files and resumable; temperature units checked (Kelvin
  converted, values divided by 10 by mistake detected).
* `download_chelsa_monthly()`: CHELSA V2.1 monthly climatologies for the
  present or any future period, GCM and scenario, optionally only for a
  region (read remotely with GDAL, values, scale and offset kept exactly),
  resumable. `chelsa_monthly_urls()` finds the files on the CHELSA server.
* `monthly_files()`: the 12 monthly files of a variable in a folder, in
  order, optionally filtered by GCM, scenario or period.

# preENM 2.8.0

* `select_models()` reports selected models that are far worse than the
  other selected models on a criterion without a required value (e.g. a
  model with good AUC and fossil AUC but much higher MAEbal) and suggests
  the `max_metrics`/`min_metrics` to add (`attr(x, "outliers")`). Criteria
  only decide the Pareto front and the ranking; models pass through the
  required values.
* Error bars in `plot_selection()`, `plot_pareto()`, `plot()` and
  `plot_models()`: drawn in a distinct colour (`bar_colour`, default red)
  so they no longer blend with grid and reference lines; taken from the
  bootstrap interval (`_lower`, `_upper`) or else +/- 1 SD (`_sd`) for any
  metric on either axis; `error_bars = FALSE` hides them. The subtitle says
  which kind each axis shows.
* `evaluate_models(boot_metrics = "all")`: bootstrap intervals for every
  metric (MAE, MAEbal, bias, TSS, Boyce), not only the criteria.

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
