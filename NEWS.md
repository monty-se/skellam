# skellam (development version)

## Bug fixes

- `dskellam()` is now accurate (relative error below 1e-12) for rates up to at least 1e7, and `log = TRUE` is finite in the deepest tails. Previously, for rates above about 5e4 the pmf was computed by a saddlepoint approximation that applied only half of its correction term (the probabilities summed to 1 - 1/(32 lambda)), `log = TRUE` returned natural-scale probabilities there, and `dskellam(500, 3, 2, log = TRUE)` was `-Inf`. The new core (`R/bessel.R`) evaluates the log Bessel function by the uniform asymptotic expansion for large orders, a stable downward recurrence for large arguments, and the power series where `besselI()` underflows.
- `dskellam.sp()` now honours `log = TRUE` and applies the full second-order saddlepoint correction.
- `skellam.mle()` and `skellam.reg()` no longer switch off warnings for the rest of the R session (they set `options(warn = -1)` and never restored it).

- `dskellam()` returns `NA` for missing `x` instead of stopping with "missing value where TRUE/FALSE needed"; one `NA` no longer breaks a whole vector.
- `dskellam()` returns 0 (or `-Inf` on the log scale) for non-integer `x`, as `dpois()` does, instead of the density at `trunc(x)`. Values within 1e-7 of an integer are treated as that integer.
- `dskellam()` warns when argument lengths are not multiples of each other (the check could never fire), keeps `names(x)`, and returns 0 for infinite `x`.
- `pskellam(q, 0, 0)` is 0 or 1 instead of `NaN`, and any zero rate is now computed exactly with `ppois()` (previously an exact 0 from `pchisq()` was sent to the saddlepoint approximation, which divides by the zero rate).
- `pskellam()` returns `NA` for missing `q` (was `NaN`), 0/1 for `q = -Inf`/`Inf` (was `NaN` with warnings), warns about recycling, and keeps `names(q)`.
- `qskellam()` with `lower.tail = FALSE`, or with `log.p = TRUE`, no longer returns one more than the correct quantile whenever `p` is exactly a cdf value. `qskellam(pskellam(x, ...), ...) == x` now holds on both tails and both scales.
- `qskellam(1, ...)` returns `Inf` (was `-Inf`), and `p = 0`/`p = 1` give the ends of the support on both tails and on the log scale (`qskellam(0, 3, 2, log.p = TRUE)` stopped with an error).
- `qskellam()` with `lambda1 = 0` no longer returns one more than the correct quantile at cdf values.
- `qskellam()` validates its input: `p` outside [0, 1] and negative rates give `NaN` with a warning (were `Inf` or an error), missing values give `NA` (were an error). It warns about recycling and keeps `names(p)`.
- `qskellam()` no longer gives different answers when `options(verbose = TRUE)` is set.
- `rskellam()` returns an integer vector, like `rpois()`; invalid rates give `NA` with "NAs produced" (was `NaN`). The draws themselves are unchanged for a given seed.

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