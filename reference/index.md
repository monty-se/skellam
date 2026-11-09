# Package index

## The Skellam distribution

Density, distribution function, quantile function and random generation,
with saddlepoint approximations to the density and distribution
function.

- [`dskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  [`dskellam.sp()`](https://monty-se.github.io/skellam/reference/skellam.md)
  [`pskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  [`pskellam.sp()`](https://monty-se.github.io/skellam/reference/skellam.md)
  [`qskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  [`rskellam()`](https://monty-se.github.io/skellam/reference/skellam.md)
  : The Skellam Distribution

## Moments and generating functions

Theoretical moments and cumulants of any order, and the moment,
characteristic and probability generating functions.

- [`skellam.moments()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  [`print(`*`<skellam_moments>`*`)`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  [`skellam.mgf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  [`skellam.cf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  [`skellam.pgf()`](https://monty-se.github.io/skellam/reference/skellam.moments.md)
  : Moments and Generating Functions of the Skellam Distribution

## Estimation

Fitting the distribution to a sample of differences, and regression with
covariates on both rates.

- [`skellam.fit()`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`print(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`summary(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`coef(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`vcov(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`logLik(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`nobs(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`confint(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`predict(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  [`simulate(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
  : Fit the Skellam Distribution to a Sample
- [`skellam.reg()`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
  : Skellam Regression
- [`print(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`summary(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`coef(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`vcov(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`logLik(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`nobs(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`fitted(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`residuals(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`predict(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`simulate(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  [`anova(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md)
  : Methods for Skellam Regression Fits

## Goodness of fit

Testing whether a sample comes from a Skellam distribution.

- [`skellam.test()`](https://monty-se.github.io/skellam/reference/skellam.test.md)
  : Goodness-of-Fit Test for the Skellam Distribution

## Diagnostic plots

Rootograms, PIT histograms and quantile-residual plots comparing data
with a fit, and the package’s ggplot2 theme and colours.

- [`plot(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam-plots.md)
  [`plot(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam-plots.md)
  [`autoplot(`*`<skellam_fit>`*`)`](https://monty-se.github.io/skellam/reference/skellam-plots.md)
  [`autoplot(`*`<skellam_reg>`*`)`](https://monty-se.github.io/skellam/reference/skellam-plots.md)
  : Diagnostic Plots for Skellam Fits
- [`skellam_palette()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  [`theme_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  [`scale_colour_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  [`scale_fill_skellam()`](https://monty-se.github.io/skellam/reference/skellam_palette.md)
  : Colours and Theme for skellam Plots

## Data

Example data used in the articles.

- [`icu_bavaria`](https://monty-se.github.io/skellam/reference/icu_bavaria.md)
  : COVID-19 Intensive Care Occupancy in Bavaria

## Deprecated

Kept for backward compatibility until version 1.0.0.

- [`skellam.mle()`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
  : Maximum Likelihood Estimation for the Skellam Distribution
  (deprecated)
