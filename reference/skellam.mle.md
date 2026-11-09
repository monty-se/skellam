# Maximum Likelihood Estimation for the Skellam Distribution (deprecated)

Estimates the parameters of a Skellam distribution using maximum
likelihood.

## Usage

``` r
skellam.mle(x)
```

## Arguments

- x:

  A vector of integers (positive or negative).

## Value

A list with components:

- iters:

  Number of log-likelihood evaluations used by the optimizer.

- loglik:

  Maximized log-likelihood value.

- param:

  Estimated parameters (\\\hat{\lambda}\_1\\, \\\hat{\lambda}\_2\\).

## Details

**Deprecated:** use
[`skellam.fit`](https://monty-se.github.io/skellam/reference/skellam.fit.md),
which returns a fitted-model object with standard errors and works with
`coef`, `vcov`, `logLik`, `AIC`, `confint` and `summary`. `skellam.mle`
gives the same estimates, keeps its old return value, and warns once per
session. It will be removed in version 1.0.0.

At the maximum likelihood estimate the fitted mean equals the sample
mean, \\\hat\lambda_1 - \hat\lambda_2 = \bar{x}\\, so the function
maximizes the profile log-likelihood over \\\lambda_2\\ alone, with
\\\lambda_1 = \lambda_2 + \bar{x}\\, using
[`optimize`](https://rdrr.io/r/stats/optimize.html) on the admissible
range \\\lambda_2 \ge \max(0, -\bar{x})\\, which keeps both rates
non-negative.

The estimate can lie on the boundary, with one rate equal to 0 (the data
are then best described by a Poisson or a reflected Poisson). This
happens when the sample variance does not exceed the absolute sample
mean, which a Skellam with two positive rates cannot produce (its
variance is \\\lambda_1 + \lambda_2 \ge \|\lambda_1 - \lambda_2\|\\);
for example, a sample of zeros gives \\(0, 0)\\ and a constant sample
\\c \> 0\\ gives \\(c, 0)\\.

## References

- Butler, R. (2007) *Saddlepoint Approximations with Applications*,
  Cambridge University Press.

- Johnson, N. L. (1959) On an extension of the connection between
  Poisson and \\\chi^2\\ distributions. *Biometrika* **46**, 352-362.

- Johnson, N. L.; Kotz, S.; Kemp, A. W. (1993) *Univariate Discrete
  Distributions*, 2nd ed., John Wiley and Sons.

- Skellam, J. G. (1946) The frequency distribution of the difference
  between two Poisson variates belonging to different populations.
  *Journal of the Royal Statistical Society, Series A* **109**(3), 296.

- Strackee, J.; van der Gon, J. J. D. (1962) The frequency distribution
  of the difference between two Poisson variates. *Statistica
  Neerlandica* **16**(1), 17-23.

- Abdulhamid, A. A.; Maha, A. O. (2010) On The Poisson Difference
  Distribution Inference and Applications. *Bulletin of the Malaysian
  Mathematical Sciences Society* **33**(1), 17-45.

- Wikipedia: Skellam distribution
  <https://en.wikipedia.org/wiki/Skellam_distribution>

## Author

Michail Tsagris

## Examples

``` r
# Basic example
set.seed(1)
x1 <- rpois(1000, 10)
x2 <- rpois(1000, 6)
x <- x1 - x2
skellam.mle(x)
#> Warning: skellam.mle() is deprecated as of skellam 0.3.0 and will be removed in 1.0.0; use skellam.fit() instead.
#> This warning is shown once per session.
#> $iters
#> [1] 25
#> 
#> $loglik
#> [1] -2869.714
#> 
#> $param
#>      mu 1      mu 2 
#> 11.001598  7.217598 
#> 

# Larger sample size
x1 <- rpois(10000, 10)
x2 <- rpois(10000, 6)
x <- x1 - x2
skellam.mle(x)
#> $iters
#> [1] 22
#> 
#> $loglik
#> [1] -28072.52
#> 
#> $param
#>      mu 1      mu 2 
#> 10.042114  6.021814 
#> 
```
