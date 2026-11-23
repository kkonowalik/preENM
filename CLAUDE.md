# preENM: conventions for work on this package

These are the maintainer's standing requests; follow them in every change
and in every session.

## Before writing code
- Read DESCRIPTION (for the current version), NAMESPACE, NEWS.md and
  README.md, then the specific R/ files the task names.
- Repository: github.com/kkonowalik/preENM. Releases go to `main`; check
  that the GitHub Actions run (R-CMD-check) passes after every push.

## Design rule (must not be broken)
- **preENM never fits a model.** It covers the steps around fitting ENM/SDM
  models with any external software.
  - Before modelling it consumes occurrence tables and raster stacks and
    emits cleaned data, calibration areas, selected predictors and sampling
    designs.
  - After modelling it consumes prediction rasters, extracted prediction
    values or variable-importance tables from any software, and emits
    metrics, rankings, diagnostics and reports.
- If a feature seems to need preENM to call Maxent, a GLM or a random forest,
  the design is wrong: restructure it to take the user's model output as
  input instead (as `transfer_bias_test()` takes predictions, or an optional
  user-supplied `fit` function).

## Dependencies
- Imports are deliberately minimal: dbscan, sf, terra, ggplot2, grid, rlang,
  stats, tools, utils. **Do not add a package to Imports.**
- Anything else goes in Suggests and is used conditionally via
  `requireNamespace()`, with an informative error if absent (e.g.
  "Install plotly: install.packages(\"plotly\")"). Tests that need it use
  `skip_if_not_installed()`; network tests also `skip_on_cran()`,
  `skip_on_ci()` and `skip_if_offline()`.

## Documentation
- roxygen2 with markdown enabled (RoxygenNote 7.3.1). Document every
  exported function with `@param`, `@return`, `@examples` and `@references`.
  Regenerate NAMESPACE and man/ with roxygen (`devtools::document()`).
- Examples fast and runnable; wrap anything slow or network-dependent in
  `\donttest{}`.
- Update with every change: README (the manual; the overview table and the
  section of the step), NEWS.md (a new version section in the voice of the
  existing entries: minor version for a feature, patch for a fix), the
  Version and, for a feature, the Description field in DESCRIPTION, and the
  package help (`R/preENM-package.R`).

## Output objects and printed output
- Return a classed list with a print method that **leads with the result in
  plain words**: an "IN SHORT" summary, as the existing modules do (variable
  selection, M, niche coverage, transfer bias).
- Plain sentences first, definitions and details at the end (not in
  parentheses in the middle of a sentence).
- Warn when a result is trivial or cannot be judged (e.g. all records inside
  the region) instead of printing a confident verdict. Graded advice is
  better than a single cut-off; say which limits come from the literature
  and which are rules of thumb.
- Where a result is worth plotting, add a `plot_*` function in the existing
  style; where it is worth writing out, add a `write_*` function.

## Figures
- **Every function that draws a figure can save it.** Give it
  `file = NULL, width = <default>, height = <default>` and, at the top of the
  function, `if (!is.null(file)) return(save_call(match.call(),
  parent.frame(), file, width, height))` (helpers `save_figure()` and
  `save_call()` in `R/plot.R`). `file` accepts one or more paths; the format
  follows the extension, and a name without extension is saved as both PNG
  and PDF. Functions that build a result and usually show it (e.g.
  `niche_coverage()`) take `file` too and save the plot directly.
- Figures must be fit for a paper or supplement: no advice or verdicts in
  titles or subtitles, only what is shown; legends with counts and
  percentages (one decimal below 1%); subtitles that fit the panel (about
  80 characters per line in a 7-inch panel; wrap with `\n`); no labels that
  can be misread (e.g. "best" only when the shading is direction-aware, and
  say so). Colours distinct from each other and from the background
  (ColorBrewer where a palette is needed); error bars in a colour other
  than the points and grid.

## Generality
- The package is for any species, region and period. Examples from one
  study (e.g. *Dryas octopetala*) are illustrations only: no
  species-specific defaults, names or claims in code or manual.

## Code and checks
- R sources ASCII only (CRAN); `\uXXXX` escapes for symbols.
- Memory: never read or rasterize a whole raster when values at points or a
  sample suffice (users work with global 30 arc-second layers); read values
  with `terra::extract()` at the points needed.
- Tests: testthat edition 3, one file per R/ module,
  `tests/testthat/test-<module>.R`. Reuse or extend the fixtures in
  `tests/testthat/helper-climate.R` and `helper-occurrences.R`;
  `simulate_enm_data()` is available for examples. Every new function gets
  tests with **known answers** (a case where the right result is known in
  advance), not only "runs without error".
- Work in small commits. Before saying a task is done, run
  `devtools::document()`, `devtools::test()` and `R CMD check --as-cran`, and
  report every NOTE or WARNING instead of suppressing it (in the cloud
  container, the locale warning, the clock note and missing Suggests are
  environment-only; say so).
