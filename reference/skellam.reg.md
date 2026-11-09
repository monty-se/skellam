# Skellam Regression

Fits a regression model assuming a Skellam distribution for the response
variable.

## Usage

``` r
skellam.reg(y, x)
```

## Arguments

- y:

  A vector of integers (positive or negative)

- x:

  A matrix, vector or data.frame of covariates

## Value

A list with components:

- loglik:

  Maximized log-likelihood value

- param1:

  Matrix for \\\lambda_1\\ parameters:

  - Column 1: Estimated coefficients

  - Column 2: Standard errors

  - Column 3: Wald z statistics (coef/se)

  - Column 4: p-values (Wald test)

- param2:

  Matrix for \\\lambda_2\\ parameters (same structure as param1)

## Details

The response is modelled as \\y_i = Y\_{1i} - Y\_{2i}\\ with \\Y\_{1i}
\sim \mathrm{Poisson}(\lambda\_{1i})\\ and \\Y\_{2i} \sim
\mathrm{Poisson}(\lambda\_{2i})\\, so that `param1` describes the rate
of the counts that are added and `param2` the rate of the counts that
are subtracted, as in `dskellam(y, lambda1, lambda2)`. Both rates use a
log link, \\\log\lambda\_{ki} = x_i^\top\beta_k\\, with the same
covariates (an intercept is added).

The log-likelihood is computed with `dskellam(log = TRUE)` and maximized
by [`optim`](https://rdrr.io/r/stats/optim.html) (BFGS) with the
analytic gradient, starting deterministically from the moment estimates
of an intercept-only model, so repeated calls give identical results.
Standard errors come from the inverse of the observed information
([`optimHess`](https://rdrr.io/r/stats/optim.html)). Observations with a
missing response or covariate are removed with a warning. A warning is
also given if the optimizer does not converge or the information matrix
is not invertible (standard errors are then `NA`).

## References

- Skellam, J. G. (1946) The frequency distribution of the difference
  between two Poisson variates belonging to different populations.
  *Journal of the Royal Statistical Society, Series A* **109**(3), 296.

- Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution
  of the difference between two Poisson variates. *Statistica
  Neerlandica* **16**(1), 17-23.

- Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football
  outcomes: using the Skellam's distribution for the goal difference.
  *IMA Journal of Management Mathematics* **20**(2), 133-145.

## Author

Michail Tsagris

## Examples

``` r
set.seed(0)
x <- rnorm(100)
y1 <- rpois(100, exp(1 + 1 * x))
y2 <- rpois(100, exp(-1 + 1 * x))
y <- y1 - y2
skellam.reg(y, x)   # param1 near (1, 1), param2 near (-1, 1)
#> $loglik
#> [1] -195.8916
#> 
#> $param1
#>              Estimate Std. Error Wald value      p-value
#> (Intercept) 0.8944038 0.09865942   9.065569 1.239554e-19
#> x           0.9908116 0.06145987  16.121276 1.808368e-58
#> 
#> $param2
#>               Estimate Std. Error Wald value    p-value
#> (Intercept) -1.1752637  0.6720832 -1.7486879 0.08034498
#> x            0.1933462  0.7030906  0.2749947 0.78332035
#> 
```
