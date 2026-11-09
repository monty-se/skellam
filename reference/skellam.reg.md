# Skellam Regression

Fits a regression model in which a difference of counts \\y = Y_1 -
Y_2\\ follows a Skellam distribution and both rates depend on covariates
through log links.

## Usage

``` r
skellam.reg(formula, data = NULL, ...)
```

## Arguments

- formula:

  A formula `y ~ x1 + x2` (the same covariates for both rates) or a
  two-part formula `y ~ x1 + x2 | z1 + z2`, where the part before `|` is
  the model for \\\lambda_1\\ and the part after it the model for
  \\\lambda_2\\. Offsets can be given with
  [`offset()`](https://rdrr.io/r/stats/offset.html) in either part. For
  the original interface, a numeric vector `y`; see Details.

- data:

  A data frame containing the variables in the formula. For the original
  interface, the covariates `x`.

- ...:

  For the original interface, `y` and `x` may also be given by name.

## Value

An object of class `"skellam_reg"`, a list with components

- coefficients:

  A list with the coefficient vectors `lambda1` and `lambda2`.

- vcov:

  Covariance matrix of all coefficients.

- loglik:

  Maximised log-likelihood.

- param1, param2:

  Coefficient tables (estimate, standard error, Wald statistic, p-value)
  for \\\lambda_1\\ and \\\lambda_2\\, as returned by earlier versions.

- fitted.values:

  A data frame with the fitted rates `lambda1` and `lambda2`.

- y, nobs, df.residual, converged, call, formula, terms, levels:

  The response used, the number of observations, the residual degrees of
  freedom, whether the optimizer converged, and what is needed to
  predict on new data.

Methods are available for `print`, `summary`, `coef`, `vcov`, `logLik`
(hence `AIC`, `BIC`), `nobs`, `confint` (the default method), `fitted`,
`residuals`, `predict`, `simulate`, `anova` and `update`.

## Details

The response is modelled as \\y_i = Y\_{1i} - Y\_{2i}\\ with independent
\\Y\_{1i} \sim \mathrm{Poisson}(\lambda\_{1i})\\ and \\Y\_{2i} \sim
\mathrm{Poisson}(\lambda\_{2i})\\, and \\\log\lambda\_{1i} =
x\_{1i}^\top\beta_1\\, \\\log\lambda\_{2i} = x\_{2i}^\top\beta_2\\. So
the `lambda1` coefficients describe the rate of the counts that are
added and the `lambda2` coefficients the rate of the counts that are
subtracted, as in `dskellam(y, lambda1, lambda2)`.

The log-likelihood is computed with `dskellam(log = TRUE)` and maximised
by [`optim`](https://rdrr.io/r/stats/optim.html) (BFGS) with the
analytic gradient, from deterministic starting values (least-squares
fits of the log moment estimates of an intercept-only model), so
repeated calls give identical results. Standard errors come from the
inverse of the observed information
([`optimHess`](https://rdrr.io/r/stats/optim.html)). Observations with
missing values in any variable used are removed with a warning. A
warning is also given if the optimizer does not converge or the
information matrix is not invertible (standard errors are then `NA`).

**Original interface.** `skellam.reg(y, x)` with a numeric response `y`
and a vector, matrix or data frame of covariates `x` fits the model with
the same covariates (and an intercept) for both rates, exactly as in
earlier versions of the package.

## References

- Skellam, J. G. (1946) The frequency distribution of the difference
  between two Poisson variates belonging to different populations.
  *Journal of the Royal Statistical Society, Series A* **109**(3), 296.

- Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football
  outcomes: using the Skellam's distribution for the goal difference.
  *IMA Journal of Management Mathematics* **20**(2), 133-145.

## See also

[`predict.skellam_reg`](https://monty-se.github.io/skellam/reference/skellam_reg-methods.md),
[`skellam.fit`](https://monty-se.github.io/skellam/reference/skellam.fit.md)
for a single sample,
[`skellam`](https://monty-se.github.io/skellam/reference/skellam.md) for
the distribution.

## Author

Michail Tsagris; formula interface and methods by Montasser Ghachem and
Oguz Ersan

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(attack = rnorm(n), defence = rnorm(n))
d$goaldiff <- rpois(n, exp(0.4 + 0.3 * d$attack)) -
              rpois(n, exp(0.1 - 0.2 * d$defence))

# different covariates for the two rates
m <- skellam.reg(goaldiff ~ attack | defence, data = d)
m
#> Skellam regression
#> 
#> Call:
#> skellam.reg(formula = goaldiff ~ attack | defence, data = d)
#> 
#> Coefficients for log(lambda1), the rate of the counts that are added:
#> (Intercept)      attack 
#>      0.4848      0.2097 
#> 
#> Coefficients for log(lambda2), the rate of the counts that are subtracted:
#> (Intercept)     defence 
#>     0.22359    -0.07533 
#> 
#> n = 500, log-likelihood -974.2896, AIC 1956.579
summary(m)
#> 
#> Call:
#> skellam.reg(formula = goaldiff ~ attack | defence, data = d)
#> 
#> log(lambda1): rate of the counts that are added
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.48485    0.06562   7.389 1.48e-13 ***
#> attack       0.20972    0.04320   4.855 1.20e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> log(lambda2): rate of the counts that are subtracted
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)  0.22359    0.08114   2.755  0.00586 **
#> defence     -0.07533    0.05224  -1.442  0.14926   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> n = 500, residual df = 496, log-likelihood -974.2896, AIC 1956.579, BIC 1973.438

# nested models: likelihood-ratio test
m0 <- skellam.reg(goaldiff ~ 1 | 1, data = d)
anova(m0, m)
#> Likelihood ratio test
#> 
#> Model 1: goaldiff ~ 1 | 1
#> Model 2: goaldiff ~ attack | defence
#> 
#>         #Df  LogLik Df  Chisq Pr(>Chisq)    
#> Model 1   2 -987.07                         
#> Model 2   4 -974.29  2 25.564  2.811e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# probabilities of each goal difference for new data
predict(m, newdata = data.frame(attack = 1, defence = 0), type = "prob", at = -2:2)
#>             -2        -1         0         1         2
#> [1,] 0.0648087 0.1400983 0.2158241 0.2243766 0.1662352

# original interface: same covariates (and an intercept) for both rates
skellam.reg(d$goaldiff, d$attack)$param1
#>              Estimate Std. Error Wald value      p-value
#> (Intercept) 0.4854379 0.06578229   7.379462 1.589310e-13
#> x           0.2097604 0.06084470   3.447472 5.658591e-04
```
