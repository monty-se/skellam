# Fit the Skellam Distribution to a Sample

Estimates the rates \\(\lambda_1, \lambda_2)\\ of a Skellam distribution
from a sample of differences of counts, by maximum likelihood (default)
or by the method of moments, and returns a fitted-model object that
works with the usual extractor functions.

## Usage

``` r
skellam.fit(x, method = c("mle", "moments"), na.rm = FALSE)

# S3 method for class 'skellam_fit'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'skellam_fit'
summary(object, level = 0.95, ...)

# S3 method for class 'skellam_fit'
coef(object, ...)

# S3 method for class 'skellam_fit'
vcov(object, ...)

# S3 method for class 'skellam_fit'
logLik(object, ...)

# S3 method for class 'skellam_fit'
nobs(object, ...)

# S3 method for class 'skellam_fit'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'skellam_fit'
predict(object, x, type = c("pmf", "cdf", "quantile"), p, ...)

# S3 method for class 'skellam_fit'
simulate(object, nsim = 1, seed = NULL, ...)
```

## Arguments

- x:

  A numeric vector of integers (differences of counts; any sign).

- method:

  `"mle"` (default) for maximum likelihood, or `"moments"` for the
  method-of-moments estimates \\\lambda_1 = (s^2 + \bar x)/2\\,
  \\\lambda_2 = (s^2 - \bar x)/2\\, truncated at 0.

- na.rm:

  Logical; if `TRUE`, missing values are removed before fitting. If
  `FALSE` (default), missing values give an error.

- object, digits, ..., level, parm, type, nsim, seed:

  Arguments of the methods: `object` is a `"skellam_fit"` object; `parm`
  and `level` as in [`confint`](https://rdrr.io/r/stats/confint.html);
  `type` in `predict` is `"pmf"` (default), `"cdf"` or `"quantile"`,
  evaluated at `x` (quantiles at probabilities `p`); `nsim` and `seed`
  as in [`simulate`](https://rdrr.io/r/stats/simulate.html).

- p:

  Probabilities, for `predict(type = "quantile")`.

## Value

An object of class `"skellam_fit"`, a list with components

- estimate:

  Named vector `c(lambda1 = , lambda2 = )`.

- vcov:

  2 x 2 covariance matrix of the estimates (`NA` on the boundary and for
  the moment method).

- loglik:

  Log-likelihood at the estimates.

- nobs:

  Number of observations used.

- method:

  `"mle"` or `"moments"`.

- boundary:

  `TRUE` if one rate is estimated as 0.

- evals:

  Number of log-likelihood evaluations (maximum likelihood).

- x:

  The data used.

- call:

  The matched call.

Methods are available for `print`, `summary`, `coef`, `vcov`, `logLik`
(hence `AIC` and `BIC`), `nobs`, `confint`, `predict` and `simulate`.

## Details

**Maximum likelihood.** At the maximum the fitted mean equals the sample
mean, \\\hat\lambda_1 - \hat\lambda_2 = \bar{x}\\, so the log-likelihood
is maximised over \\\lambda_2\\ alone with
[`optimize`](https://rdrr.io/r/stats/optimize.html), on the range that
keeps both rates non-negative. The covariance matrix is the inverse of
the observed information, computed from the analytic score.

**Boundary estimates.** A Skellam with two positive rates has variance
\\\lambda_1 + \lambda_2 \ge \|\lambda_1 - \lambda_2\|\\. When the sample
variance does not exceed the absolute sample mean, the maximum lies on
the boundary, with one rate equal to 0 (a Poisson or a reflected
Poisson): a sample of zeros gives \\(0, 0)\\, a constant sample \\c \>
0\\ gives \\(c, 0)\\. Standard errors do not exist there; the covariance
matrix is then `NA` and the fit is flagged with `boundary = TRUE` and a
warning.

**Moments.** The moment estimates are returned without a covariance
matrix and serve mainly as a quick check or a starting point.

Confidence intervals are Wald intervals on the log scale,
\\\exp(\log\hat\lambda \pm z\\\mathrm{se}/\hat\lambda)\\, which keeps
them positive; they are `NA` when the covariance matrix is.

## References

Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football
outcomes: using the Skellam's distribution for the goal difference. *IMA
Journal of Management Mathematics* **20**(2), 133-145.

## See also

[`skellam`](https://monty-se.github.io/skellam/reference/skellam.md) for
the distribution functions,
[`skellam.reg`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
for regression on covariates.

## Examples

``` r
set.seed(1)
x <- rpois(500, 6) - rpois(500, 4)
fit <- skellam.fit(x)
fit
#> Skellam distribution fitted by maximum likelihood (n = 500)
#> 
#>         Estimate Std. Error
#> lambda1    5.683     0.3126
#> lambda2    3.717     0.3063
#> 
#> log-likelihood -1268.513, AIC 2541.026
summary(fit)
#> 
#> Call:
#> skellam.fit(x = x)
#> 
#> Skellam distribution fitted by maximum likelihood (n = 500)
#> 
#>         Estimate Std. Error 2.5 % 97.5 %
#> lambda1    5.683     0.3126 5.102  6.330
#> lambda2    3.717     0.3063 3.163  4.368
#> 
#> Sample and fitted moments:
#>                   sample  fitted
#> mean             1.96600 1.96600
#> variance         9.36357 9.39958
#> skewness         0.04562 0.06822
#> excess kurtosis -0.21702 0.10639
#> 
#> var(x) - |mean(x)| = 7.398 (must be >= 0 for two positive rates)
#> log-likelihood -1268.513, AIC 2541.026, BIC 2549.455
coef(fit)
#>  lambda1  lambda2 
#> 5.682792 3.716792 
confint(fit)
#>            2.5 %   97.5 %
#> lambda1 5.101969 6.329738
#> lambda2 3.162510 4.368223
AIC(fit)
#> [1] 2541.026

# probabilities and quantities from the fitted distribution
predict(fit, x = -2:2)
#> [1] 0.05635807 0.08278187 0.10844110 0.12656940 0.13174771
predict(fit, x = 0, type = "cdf")
#> [1] 0.3166358

# the boundary: a constant sample is best fitted by a Poisson
skellam.fit(rep(3, 50))
#> Warning: the estimate lies on the boundary (a rate is 0): the sample variance does not exceed the absolute mean; standard errors are not available
#> Skellam distribution fitted by maximum likelihood (n = 50)
#> 
#>         Estimate Std. Error
#> lambda1        3         NA
#> lambda2        0         NA
#> 
#> log-likelihood -74.79613, AIC 153.5923
#> Estimate on the boundary: one rate is 0 (standard errors not available).
```
