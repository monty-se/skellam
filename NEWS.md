# skellam (development version)

## Infrastructure

- Added a testthat test suite (`tests/testthat/`) replacing `tests/skellam.R`, which only drew plots.
- Added `dev/skellam_diagnostics.R`, an annotated health check of the package.
- Removed build artefacts (`doc/`, `Meta/`, `build/`, `.Rproj.user/`) from version control; installing from a clone failed because of them.
- `DESCRIPTION`: added Oguz Ersan as author and a `BugReports` field; removed the fields CRAN writes itself (`Date`, `Packaged`, `Repository`).

# skellam 0.2.4

---

## New Features and Improvements
- Added a comprehensive `README.md` with usage examples, installation instructions, and badges (e.g., CRAN status, license).
- Improved the **vignette** (`validate_skellam.Rmd`) with clearer explanations, reproducible examples, and better visualizations.

## Documentation Updates
- Revised documentation for `skellam.reg` to clarify parameter definitions, return values, and included practical examples.

## Internal Changes
- Updated vignette-building dependencies (now requires `rmarkdown` and `knitr`).
- Minor fixes to ensure compatibility with `R CMD check` across all tested platforms (see `cran-results.md`).

---