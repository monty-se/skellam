# Changelog

## skellam (development version)

### New features

- [`skellam.fit()`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  fits the Skellam distribution to a sample by maximum likelihood (or
  the method of moments) and returns a `"skellam_fit"` object with
  standard errors from the observed information. It works with
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) (estimates,
  standard errors and confidence intervals, sample versus fitted
  moments, AIC and BIC), [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html) (hence
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) and
  [`BIC()`](https://rdrr.io/r/stats/AIC.html)),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html) (pmf, cdf or
  quantiles of the fitted distribution) and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html). Boundary
  estimates (one rate 0) are flagged with a warning and `NA` standard
  errors; `na.rm = TRUE` drops missing values.

- [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  has a formula interface, `skellam.reg(y ~ x1 + x2 | z1, data)`: the
  part before `|` models `lambda1`, the part after it `lambda2` (one
  part means the same covariates for both). Factors, transformations and
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms work as in
  [`glm()`](https://rdrr.io/r/stats/glm.html); incomplete rows are
  removed with a warning. The result is a `"skellam_reg"` object with
  methods for `print`, `summary`, `coef`, `vcov`, `logLik` (hence `AIC`,
  `BIC`), `nobs`, `confint`, `fitted`, `residuals` (response, Pearson
  and randomized quantile residuals), `predict` (mean, rates, or
  probabilities and cdf at given values, on new data), `simulate`,
  `anova` (likelihood-ratio tests of nested models) and `update`. The
  original call `skellam.reg(y, x)` still works and its result still
  contains `loglik`, `param1` and `param2`.

- [`skellam.moments()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  returns the mean, variance, standard deviation, skewness and excess
  kurtosis, and the cumulants, raw and central moments of any order;
  [`skellam.mgf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md),
  [`skellam.cf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  and
  [`skellam.pgf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  are the moment, characteristic and probability generating functions.

- [`skellam.test()`](https://monty-se.github.io/skellam/reference/skellam.test.md)
  tests goodness of fit with the Pearson chi-squared or likelihood-ratio
  statistic on cells merged to a minimum expected count. With estimated
  rates the default p-value comes from a parametric bootstrap (refitting
  each simulated sample), because no chi-squared reference is exact in
  that case (Chernoff and Lehmann, 1954); an asymptotic chi-squared
  p-value is available. Returns an `"htest"` object.

- Diagnostic plots for
  [`skellam.fit()`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  and
  [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  fits: [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
  draw a hanging rootogram (default), observed versus fitted
  probabilities, a PIT histogram for count data, or a normal QQ plot of
  randomized quantile residuals. All are ggplot2 objects that can be
  modified with `+`; ggplot2 is suggested, not required.
  [`theme_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md),
  [`skellam_palette()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  (colour-blind-safe Okabe-Ito colours),
  [`scale_colour_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  and
  [`scale_fill_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  give other figures the same look.

### Deprecated

- [`skellam.mle()`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
  is deprecated in favour of
  [`skellam.fit()`](https://monty-se.github.io/skellam/reference/skellam.fit.md).
  It returns the same estimates and its old return value, warns once per
  session, and will be removed in version 1.0.0.

## skellam 0.2.5

### Breaking change

- [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md):
  `param1` now describes the rate of the counts that are **added** and
  `param2` the rate of the counts that are **subtracted**, as in
  `dskellam(y, lambda1, lambda2)` and as the documentation always said.
  The log-likelihood had the sign of the `y * log(lambda1 / lambda2)`
  term reversed, so the two blocks were swapped: for `y = Y1 - Y2`,
  `param1` estimated the rate of `Y2`. Results from earlier versions
  should be re-run; the old documented example compensated by writing
  `y <- y2 - y1`.

### Bug fixes

- [`dskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  is now accurate (relative error below 1e-12) for rates up to at least
  1e7, and `log = TRUE` is finite in the deepest tails. Previously, for
  rates above about 5e4 the pmf was computed by a saddlepoint
  approximation that applied only half of its correction term (the
  probabilities summed to 1 - 1/(32 lambda)), `log = TRUE` returned
  natural-scale probabilities there, and
  `dskellam(500, 3, 2, log = TRUE)` was `-Inf`. The new core
  (`R/bessel.R`) evaluates the log Bessel function by the uniform
  asymptotic expansion for large orders, a stable downward recurrence
  for large arguments, and the power series where
  [`besselI()`](https://rdrr.io/r/base/Bessel.html) underflows.

- [`pskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  is now accurate (relative error below 1e-11) at all rates. Above
  `lambda1 + lambda2 = 1e4` the non-central chi-square route through
  [`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) was slow, lost
  accuracy, and beyond about 5e5 stopped converging and returned a wrong
  value (`pskellam(0, 1e6, 1e6)` gave 0.875; the answer is 0.500). Those
  rates, and any tail where
  [`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) underflows, now
  use a direct log-space summation of the pmf, so far tails are exact on
  the log scale and never negative. Quantiles at large rates are
  correspondingly faster:
  [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  at rates of 1e6 takes 0.04 s instead of 0.7 s with warnings.

- [`dskellam.sp()`](https://monty-se.github.io/skellam/reference/skellam.md)
  now honours `log = TRUE` and applies the full second-order saddlepoint
  correction.

- [`skellam.mle()`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
  and
  [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  no longer switch off warnings for the rest of the R session (they set
  `options(warn = -1)` and never restored it).

- [`dskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  returns `NA` for missing `x` instead of stopping with “missing value
  where TRUE/FALSE needed”; one `NA` no longer breaks a whole vector.

- [`dskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  returns 0 (or `-Inf` on the log scale) for non-integer `x`, as
  [`dpois()`](https://rdrr.io/r/stats/Poisson.html) does, instead of the
  density at `trunc(x)`. Values within 1e-7 of an integer are treated as
  that integer.

- [`dskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  warns when argument lengths are not multiples of each other (the check
  could never fire), keeps `names(x)`, and returns 0 for infinite `x`.

- `pskellam(q, 0, 0)` is 0 or 1 instead of `NaN`, and any zero rate is
  now computed exactly with
  [`ppois()`](https://rdrr.io/r/stats/Poisson.html) (previously an exact
  0 from [`pchisq()`](https://rdrr.io/r/stats/Chisquare.html) was sent
  to the saddlepoint approximation, which divides by the zero rate).

- [`pskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  returns `NA` for missing `q` (was `NaN`), 0/1 for `q = -Inf`/`Inf`
  (was `NaN` with warnings), warns about recycling, and keeps
  `names(q)`.

- [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  with `lower.tail = FALSE`, or with `log.p = TRUE`, no longer returns
  one more than the correct quantile whenever `p` is exactly a cdf
  value. `qskellam(pskellam(x, ...), ...) == x` now holds on both tails
  and both scales.

- `qskellam(1, ...)` returns `Inf` (was `-Inf`), and `p = 0`/`p = 1`
  give the ends of the support on both tails and on the log scale
  (`qskellam(0, 3, 2, log.p = TRUE)` stopped with an error).

- [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  with `lambda1 = 0` no longer returns one more than the correct
  quantile at cdf values.

- [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  validates its input: `p` outside \[0, 1\] and negative rates give
  `NaN` with a warning (were `Inf` or an error), missing values give
  `NA` (were an error). It warns about recycling and keeps `names(p)`.

- [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  no longer gives different answers when `options(verbose = TRUE)` is
  set.

- [`rskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  returns an integer vector, like
  [`rpois()`](https://rdrr.io/r/stats/Poisson.html); invalid rates give
  `NA` with “NAs produced” (was `NaN`). The draws themselves are
  unchanged for a given seed.

- [`skellam.mle()`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
  no longer returns negative rates. For samples whose variance does not
  exceed the absolute mean (including constant samples) it returned a
  negative “Poisson rate” such as (1.5, -1.5); for an all-zero sample it
  stopped with “non-finite value supplied by ‘nlm’”. It now maximizes
  the profile likelihood over the admissible range with
  [`optimize()`](https://rdrr.io/r/stats/optimize.html) and returns the
  boundary estimate when that is the maximum: (0, 0) for zeros, (c, 0)
  for a constant c \> 0. Missing, non-integer or too few observations
  give clear errors. The return value keeps its components (`iters` now
  counts log-likelihood evaluations).

- [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  now works at any rate, gives identical results on repeated calls, and
  reports problems. It used
  [`besselI()`](https://rdrr.io/r/base/Bessel.html) without scaling,
  which overflowed for rates above about 360 and underflowed at its
  starting values when the data had large differences (“function cannot
  be evaluated at initial parameters”); it started from random values
  ([`rnorm()`](https://rdrr.io/r/stats/Normal.html)) without a seed, so
  two calls on the same data differed; and it chained three optimizers
  without checking convergence. It now uses `dskellam(log = TRUE)`,
  starts deterministically from intercept-only moment estimates, runs
  one BFGS optimization with the analytic gradient, and warns if the
  optimizer does not converge or the information matrix cannot be
  inverted (standard errors are then `NA`). Observations with missing
  values are removed with a warning (previously a missing covariate
  silently misaligned `y` and `x`).

### Documentation

- [`?skellam.reg`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  cites the published Karlis & Ntzoufras (2009) paper instead of a
  conference handout URL; the
  [`skellam.mle()`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
  examples set a seed.
- `citation("skellam")` lists all five authors and the GitHub URL (it
  pointed to an R-Forge project and omitted the maintainers).
- The
  [`?skellam`](https://monty-se.github.io/skellam/reference/skellam.md)
  help page now states the conventions the functions follow (recycling,
  `NA`, invalid input, the quantile definition) and describes the
  numerical methods and their accuracy; the hand-written usage block
  that duplicated the generated one is gone, and the note about R
  versions before 2.9 with it.

### Infrastructure

- The package has a website, <https://monty-se.github.io/skellam/>,
  built with pkgdown and published by GitHub Actions; its address is in
  `DESCRIPTION` and the README.
- Added a testthat test suite (`tests/testthat/`) replacing
  `tests/skellam.R`, which only drew plots.
- Added `dev/skellam_diagnostics.R`, an annotated health check of the
  package.
- Removed build artefacts (`doc/`, `Meta/`, `build/`, `.Rproj.user/`)
  from version control; installing from a clone failed because of them.
- `DESCRIPTION`: added Oguz Ersan as author and a `BugReports` field;
  removed the fields CRAN writes itself (`Date`, `Packaged`,
  `Repository`).

## skellam 0.2.4

------------------------------------------------------------------------

### New Features and Improvements

- Added a comprehensive `README.md` with usage examples, installation
  instructions, and badges (e.g., CRAN status, license).
- Improved the **vignette** (`validate_skellam.Rmd`) with clearer
  explanations, reproducible examples, and better visualizations.

### Documentation Updates

- Revised documentation for `skellam.reg` to clarify parameter
  definitions, return values, and included practical examples.

### Internal Changes

- Updated vignette-building dependencies (now requires `rmarkdown` and
  `knitr`).
- Minor fixes to ensure compatibility with `R CMD check` across all
  tested platforms (see `cran-comments.md`).

------------------------------------------------------------------------
