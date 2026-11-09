# Methods for Skellam Regression Fits

Extractor, prediction and testing methods for objects returned by
[`skellam.reg`](https://monty-se.github.io/skellam/reference/skellam.reg.md).

## Usage

``` r
# S3 method for class 'skellam_reg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'skellam_reg'
summary(object, ...)

# S3 method for class 'skellam_reg'
coef(object, model = c("full", "lambda1", "lambda2"), ...)

# S3 method for class 'skellam_reg'
vcov(object, ...)

# S3 method for class 'skellam_reg'
logLik(object, ...)

# S3 method for class 'skellam_reg'
nobs(object, ...)

# S3 method for class 'skellam_reg'
fitted(object, ...)

# S3 method for class 'skellam_reg'
residuals(object, type = c("response", "pearson", "quantile"), ...)

# S3 method for class 'skellam_reg'
predict(
  object,
  newdata = NULL,
  type = c("response", "lambda1", "lambda2", "rates", "prob", "cdf"),
  at = NULL,
  ...
)

# S3 method for class 'skellam_reg'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'skellam_reg'
anova(object, ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further fitted models (for `anova`); otherwise unused.

- object, x:

  A `"skellam_reg"` object.

- model:

  For `coef`: `"full"` (all coefficients, named `lambda1:...` and
  `lambda2:...`), `"lambda1"` or `"lambda2"`.

- type:

  For `predict`: `"response"` (the fitted mean \\\lambda_1 -
  \lambda_2\\), `"lambda1"`, `"lambda2"`, `"rates"` (both), `"prob"`
  (the probabilities of the values in `at`) or `"cdf"`. For `residuals`:
  `"response"` (\\y - (\lambda_1 - \lambda_2)\\), `"pearson"` (divided
  by \\\sqrt{\lambda_1 + \lambda_2}\\) or `"quantile"` (randomized
  quantile residuals, Dunn and Smyth 1996, standard normal when the
  model is right; they use random numbers).

- newdata:

  Optional data frame for `predict`; if omitted, the fitted values are
  used.

- at:

  Integer values at which `type = "prob"` or `"cdf"` are evaluated; one
  column per value.

- nsim, seed:

  As in [`simulate`](https://rdrr.io/r/stats/simulate.html).

## Value

`coef`, `vcov`, `logLik`, `nobs`, `fitted` and `residuals` return what
their names say; `predict` returns a vector, or a matrix for
`type = "rates"`, `"prob"` and `"cdf"`; `simulate` a data frame of
simulated responses; `anova` a likelihood-ratio test table for nested
models.

## References

Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
*Journal of Computational and Graphical Statistics* **5**(3), 236-244.
