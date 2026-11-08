# skellam (development version)

## Bug fixes

- `skellam.mle()` and `skellam.reg()` no longer switch off warnings for the rest of the R session (they set `options(warn = -1)` and never restored it).

- `dskellam()` returns `NA` for missing `x` instead of stopping with "missing value where TRUE/FALSE needed"; one `NA` no longer breaks a whole vector.
- `dskellam()` returns 0 (or `-Inf` on the log scale) for non-integer `x`, as `dpois()` does, instead of the density at `trunc(x)`. Values within 1e-7 of an integer are treated as that integer.
- `dskellam()` warns when argument lengths are not multiples of each other (the check could never fire), keeps `names(x)`, and returns 0 for infinite `x`.
- `pskellam(q, 0, 0)` is 0 or 1 instead of `NaN`, and any zero rate is now computed exactly with `ppois()` (previously an exact 0 from `pchisq()` was sent to the saddlepoint approximation, which divides by the zero rate).
- `pskellam()` returns `NA` for missing `q` (was `NaN`), 0/1 for `q = -Inf`/`Inf` (was `NaN` with warnings), warns about recycling, and keeps `names(q)`.

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