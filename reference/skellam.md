# The Skellam Distribution

Density, distribution function, quantile function and random generation
for the Skellam distribution, the distribution of the difference of two
independent Poisson variables, together with saddlepoint approximations
to the density and distribution function.

## Usage

``` r
dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)

dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)

pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

pskellam.sp(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

rskellam(n, lambda1, lambda2 = lambda1)
```

## Arguments

- x, q:

  Numeric vector of quantiles (integers for the pmf).

- lambda1, lambda2:

  Numeric vectors of non-negative rates; `lambda2` defaults to
  `lambda1`.

- log, log.p:

  Logical; if `TRUE`, probabilities are given as `log(p)`.

- lower.tail:

  Logical; if `TRUE` (default), probabilities are \\P(X \le x)\\,
  otherwise \\P(X \> x)\\.

- p:

  Numeric vector of probabilities.

- n:

  Number of observations; if `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dskellam` gives the (log) probability mass function, `pskellam` the
(log) distribution function, `qskellam` the quantile function and
`rskellam` generates random deviates (an integer vector). `dskellam.sp`
and `pskellam.sp` give the saddlepoint approximations. Names of the
first argument are kept.

Invalid `lambda` values return `NaN` with a warning; a zero rate is
allowed and reduces the distribution to a (reflected) Poisson.

## Details

If \\Y_1 \sim \mathrm{Poisson}(\lambda_1)\\ and \\Y_2 \sim
\mathrm{Poisson}(\lambda_2)\\ are independent, then \\X = Y_1 - Y_2\\
has the Skellam distribution with probability mass function \$\$P(X = x)
= e^{-(\lambda_1+\lambda_2)} (\lambda_1/\lambda_2)^{x/2}
I\_{\|x\|}(2\sqrt{\lambda_1 \lambda_2}), \qquad x = \ldots, -1, 0, 1,
\ldots\$\$ where \\I\_\nu\\ is the modified Bessel function of the first
kind. The mean is \\\lambda_1 - \lambda_2\\ and the variance
\\\lambda_1 + \lambda_2\\. (If \\Y_1\\ and \\Y_2\\ are correlated
Poisson variables with means \\\mu_i\\ and correlation \\\rho\\, the
difference is still Skellam, with \\\lambda_i = \mu_i - \rho\sqrt{\mu_1
\mu_2}\\.)

All functions follow the conventions of the stats distribution
functions: arguments are recycled to a common length (with a warning
when the lengths are not multiples), missing values propagate, invalid
rates and probabilities give `NaN` with a warning, a non-integer `x` has
density 0 (with a warning), and `qskellam` returns the smallest integer
\\x\\ with \\P(X \le x) \ge p\\ (for `lower.tail = TRUE`), so that
`qskellam(pskellam(x, ...), ...)` recovers `x`.

**Numerical methods.** `dskellam` evaluates the log of the Bessel
function directly, so that `log = TRUE` is finite and accurate for any
input: through [`besselI`](https://rdrr.io/r/base/Bessel.html) where
that is reliable, through the uniform asymptotic expansion for large
orders, a stable downward recurrence for large arguments, and the power
series where `besselI` underflows. The relative error is below
\\10^{-12}\\ for rates up to at least \\10^7\\, and below \\10^{-9}\\
for extremely unequal large rates such as \\(10^6, 10^3)\\. `pskellam`
uses the non-central chi-square identity of Johnson (1959) through
[`pchisq`](https://rdrr.io/r/stats/Chisquare.html) when \\\lambda_1 +
\lambda_2 \le 10^4\\, and otherwise, as well as in tails where `pchisq`
underflows, sums the log-pmf directly toward the nearer tail, so that
far tails are exact on the log scale. `qskellam` starts from a
Cornish-Fisher approximation and searches the integers with `pskellam`.
`rskellam` returns the difference of two
[`rpois`](https://rdrr.io/r/stats/Poisson.html) draws.

`dskellam.sp` and `pskellam.sp` are the second-order saddlepoint
approximation to the pmf and the Lugannani-Rice approximation to the cdf
with Butler's second continuity correction (Butler 2007). Their relative
error is of order \\1/(\lambda_1+\lambda_2)^2\\ and
\\1/(\lambda_1+\lambda_2)\\ respectively; they are provided for users
who want the approximations themselves and are not used by `dskellam` or
`pskellam`.

## References

- Butler, R. (2007) *Saddlepoint Approximations with Applications*,
  Cambridge University Press.

- Johnson, N. L. (1959) On an extension of the connection between
  Poisson and \\\chi^2\\ distributions. *Biometrika* 46, 352-362.

- Johnson, N. L., Kotz, S., & Kemp, A. W. (1993) *Univariate Discrete
  Distributions*, 2nd ed., John Wiley and Sons.

- Olver, F. W. J. et al. (eds.) *NIST Digital Library of Mathematical
  Functions*, Section 10.41 (uniform asymptotic expansions of Bessel
  functions). <https://dlmf.nist.gov/10.41>

- Skellam, J. G. (1946) The frequency distribution of the difference
  between two Poisson variates. *Journal of the Royal Statistical
  Society, Series A* 109(3), 296.

- Strackee, J. & van der Gon, J. J. D. (1962) The frequency distribution
  of the difference between two Poisson variates. *Statistica
  Neerlandica* 16(1), 17-23.

## See also

[`skellam.mle`](https://monty-se.github.io/skellam/reference/skellam.mle.md)
and
[`skellam.reg`](https://monty-se.github.io/skellam/reference/skellam.reg.md)
for estimation; [`dpois`](https://rdrr.io/r/stats/Poisson.html) for the
Poisson distribution.

## Examples

``` r
# Compare with Poisson when one lambda = 0
dskellam(0:10, 5, 0)
#>  [1] 0.006737947 0.033689735 0.084224337 0.140373896 0.175467370 0.175467370
#>  [7] 0.146222808 0.104444863 0.065278039 0.036265577 0.018132789
dpois(0:10, 5)
#>  [1] 0.006737947 0.033689735 0.084224337 0.140373896 0.175467370 0.175467370
#>  [7] 0.146222808 0.104444863 0.065278039 0.036265577 0.018132789

# Both lambdas non-zero
dskellam(c(-1,1), c(12,10), c(10,12))
#> [1] 0.0697968 0.0697968
pskellam(c(-1,0), c(12,10), c(10,12))
#> [1] 0.2965079 0.7034921

# Quantile function inverts the distribution function
qskellam(c(0.05, 0.95), 3, 4)
#> [1] -5  3
qskellam(pskellam(-3:3, 3, 4), 3, 4)
#> [1] -3 -2 -1  0  1  2  3

# Log-scale values stay finite far in the tails and at large rates
dskellam(500, 3, 2, log = TRUE)
#> [1] -2067.012
pskellam(0, 1e6, 1e6)
#> [1] 0.500141

# Random generation
set.seed(1)
rskellam(10, 8.5, 10.25)
#>  [1]  0 -2 -1 -5 -6  5  5  1  1 -6
```
