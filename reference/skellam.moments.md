# Moments and Generating Functions of the Skellam Distribution

Theoretical moments, cumulants and generating functions of the Skellam
distribution with rates \\\lambda_1\\ and \\\lambda_2\\.

## Usage

``` r
skellam.moments(lambda1, lambda2 = lambda1, order = 4)

# S3 method for class 'skellam_moments'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

skellam.mgf(t, lambda1, lambda2 = lambda1)

skellam.cf(t, lambda1, lambda2 = lambda1)

skellam.pgf(z, lambda1, lambda2 = lambda1)
```

## Arguments

- lambda1, lambda2:

  Non-negative rates (single numbers for `skellam.moments`, recycled
  with `t` or `z` for the generating functions); `lambda2` defaults to
  `lambda1`.

- order:

  Highest order of the cumulants and moments returned.

- x:

  A `"skellam_moments"` object.

- digits:

  Number of significant digits to print.

- ...:

  Unused.

- t:

  Real argument of the moment generating and characteristic functions.

- z:

  Non-zero argument of the probability generating function (real or
  complex).

## Value

`skellam.moments` returns an object of class `"skellam_moments"`: a list
with `mean`, `variance`, `sd`, `skewness`, `kurtosis` (excess), and the
vectors `cumulants`, `raw` and `central` of orders `1:order`. The
generating functions return numeric (complex for `skellam.cf`, and for
`skellam.pgf` with complex `z`) vectors.

## Details

The cumulant generating function of \\X = Y_1 - Y_2\\ is \\K(t) =
\lambda_1(e^t - 1) + \lambda_2(e^{-t} - 1)\\, so the cumulants are
\\\kappa_r = \lambda_1 + (-1)^r \lambda_2\\: the mean is \\\lambda_1 -
\lambda_2\\, the variance \\\lambda_1 + \lambda_2\\, the skewness
\\(\lambda_1 - \lambda_2)/(\lambda_1 + \lambda_2)^{3/2}\\ and the excess
kurtosis \\1/(\lambda_1 + \lambda_2)\\. Raw moments follow from the
cumulants by the recursion \\\mu'\_n = \sum\_{k=0}^{n-1} \binom{n-1}{k}
\kappa\_{n-k} \mu'\_k\\, and central moments by the same recursion with
\\\kappa_1 = 0\\.

The generating functions are \\M(t) = E\[e^{tX}\] = \exp\\K(t)\\\\,
\\\varphi(t) = E\[e^{itX}\] = \exp\\\lambda_1(e^{it}-1) +
\lambda_2(e^{-it}-1)\\\\ and \\G(z) = E\[z^X\] = \exp\\\lambda_1(z-1) +
\lambda_2(1/z-1)\\\\.

Sample moments of data are reported by
[`summary()`](https://rdrr.io/r/base/summary.html) of a fit from
[`skellam.fit`](https://monty-se.github.io/skellam/reference/skellam.fit.md),
next to the fitted ones.

## Examples

``` r
skellam.moments(3, 2)
#> Skellam distribution, lambda1 = 3, lambda2 = 2
#> 
#>            mean        variance              sd        skewness excess kurtosis 
#>         1.00000         5.00000         2.23607         0.08944         0.20000 
#> 
#> Cumulants, raw and central moments by order:
#>          1 2  3   4
#> cumulant 1 5  1   5
#> raw      1 6 17 115
#> central  0 5  1  80
skellam.moments(3, 2, order = 6)$raw
#>    1    2    3    4    5    6 
#>    1    6   17  115  512 3867 

# the mgf at 0 is 1 and its derivative is the mean
skellam.mgf(0, 3, 2)
#> [1] 1
(skellam.mgf(1e-6, 3, 2) - skellam.mgf(-1e-6, 3, 2)) / 2e-6
#> [1] 1

# E[z^X] equals the sum of z^x P(X = x)
x <- -60:60
c(skellam.pgf(0.7, 3, 2), sum(0.7^x * dskellam(x, 3, 2)))
#> [1] 0.9580482 0.9580482
```
