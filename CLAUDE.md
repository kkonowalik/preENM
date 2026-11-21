# preENM: conventions for work on this package

These are the maintainer's standing requests; follow them in every change.

## Figures
- **Every function that draws a figure can save it.** Give it
  `file = NULL, width = <default>, height = <default>` and, at the top of the
  function, `if (!is.null(file)) return(save_call(match.call(),
  parent.frame(), file, width, height))` (helpers in `R/plot.R`).
  `file` accepts one or more paths; the format follows the extension, and a
  name without extension is saved as both PNG and PDF. Functions that build
  a result and usually show it (e.g. `niche_coverage()`) take `file` too and
  save the plot directly.
- Plot text must say what it shows: legends with counts, subtitles that fit
  a 7-inch-wide figure, no labels that can be misread (e.g. "best" only when
  the shading is direction-aware, and say so).

## Generality
- The package is for any species, region and period. Examples from one
  study (e.g. *Dryas octopetala*) are illustrations only: no
  species-specific defaults, names or claims in code or manual.

## Printed output
- Plain sentences first, definitions and details at the end (not in
  parentheses in the middle of a sentence). An "IN SHORT" summary in words.
- Warn when a result is trivial or cannot be judged (e.g. all records inside
  the region) instead of printing a confident verdict.

## Code and checks
- R sources ASCII only (CRAN); `\uXXXX` escapes for symbols.
- Tests with known answers for every new function; `R CMD check --as-cran`
  clean (apart from environment-only notes); verify the GitHub Actions run
  after pushing.
- Update README (manual), NEWS and the version in DESCRIPTION with every
  change.
