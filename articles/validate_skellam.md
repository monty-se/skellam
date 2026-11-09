# Validating the skellam package

## What this vignette checks

If \\Y_1 \sim \mathrm{Poisson}(\lambda_1)\\ and \\Y_2 \sim
\mathrm{Poisson}(\lambda_2)\\ are independent, their difference \\X =
Y_1 - Y_2\\ has the Skellam distribution,

\\ P(X = x) = e^{-(\lambda_1+\lambda_2)}
\left(\frac{\lambda_1}{\lambda_2}\right)^{x/2}
I\_{\|x\|}\\\left(2\sqrt{\lambda_1\lambda_2}\right), \qquad x \in
\mathbb{Z}, \\

with mean \\\lambda_1 - \lambda_2\\ and variance \\\lambda_1 +
\lambda_2\\. This vignette checks, with numbers rather than by eye, that
the functions of the package agree with that definition and with each
other:

1.  `dskellam` against an exact sum over Poisson probabilities;
2.  `pskellam` against cumulative sums of `dskellam`, and `qskellam` as
    its inverse;
3.  `rskellam` against `dskellam` with a goodness-of-fit test;
4.  `skellam.mle` and `skellam.reg` against known parameters in
    simulated data.

Every random number below is drawn after `set.seed(20261108)`, so the
output is the same each time the vignette is built.

## 1. The probability mass function

By definition, \\P(X = x) = \sum_k P(Y_1 = k + x)\\P(Y_2 = k)\\. Summing
that series directly is slow but exact, which makes it a good reference.

``` r

pmf_by_sum <- function(x, l1, l2) {
  k <- max(0, -x):(abs(x) + 200 + 20 * (l1 + l2))
  sum(dpois(k + x, l1) * dpois(k, l2))
}
x <- -10:15
for (L in list(c(1.5, 0.5), c(12, 8))) {
  ref <- vapply(x, pmf_by_sum, 0, l1 = L[1], l2 = L[2])
  cat(sprintf("lambda = (%g, %g): largest relative difference %.1e\n",
              L[1], L[2], max(abs(dskellam(x, L[1], L[2]) / ref - 1))))
}
#> lambda = (1.5, 0.5): largest relative difference 6.7e-15
#> lambda = (12, 8): largest relative difference 2.7e-15
```

The probabilities also add up to one, including at rates far beyond what
the series above can handle quickly:

``` r

for (L in list(c(1.5, 0.5), c(12, 8), c(1e4, 1e4), c(1e6, 1e6))) {
  s <- sqrt(sum(L)); m <- L[1] - L[2]
  x <- seq(floor(m - 30 * s), ceiling(m + 30 * s))
  cat(sprintf("lambda = (%g, %g): 1 - sum = %.1e\n", L[1], L[2],
              1 - sum(dskellam(x, L[1], L[2]))))
}
#> lambda = (1.5, 0.5): 1 - sum = 0.0e+00
#> lambda = (12, 8): 1 - sum = -8.9e-16
#> lambda = (10000, 10000): 1 - sum = 0.0e+00
#> lambda = (1e+06, 1e+06): 1 - sum = 0.0e+00
```

When one rate is zero the distribution is a Poisson, or a Poisson turned
upside down:

``` r

all.equal(dskellam(0:10, 5, 0), dpois(0:10, 5))
#> [1] TRUE
all.equal(dskellam(-(0:10), 0, 5), dpois(0:10, 5))
#> [1] TRUE
```

## 2. Distribution and quantile functions

`pskellam(x)` should equal the running sum of `dskellam` up to `x`, and
the two tails should add up to one.

``` r

l1 <- 12; l2 <- 8
x <- -40:50
max(abs(pskellam(x, l1, l2) - cumsum(dskellam(x, l1, l2))))
#> [1] 1.110223e-15
max(abs(pskellam(x, l1, l2) + pskellam(x, l1, l2, lower.tail = FALSE) - 1))
#> [1] 0
```

For a discrete distribution the quantile function returns the smallest
integer \\x\\ with \\P(X \le x) \ge p\\. Feeding a value of the
distribution function back in must therefore return the same \\x\\, on
both tails and on the log scale:

``` r

x <- -15:25
all(qskellam(pskellam(x, l1, l2), l1, l2) == x)
#> [1] TRUE
all(qskellam(pskellam(x, l1, l2, lower.tail = FALSE), l1, l2, lower.tail = FALSE) == x)
#> [1] TRUE
all(qskellam(pskellam(x, l1, l2, log.p = TRUE), l1, l2, log.p = TRUE) == x)
#> [1] TRUE
qskellam(c(0, 0.025, 0.5, 0.975, 1), l1, l2)
#> [1] -Inf   -5    4   13  Inf
```

## 3. Random generation

`rskellam` draws the difference of two Poisson variates. We compare the
observed frequencies of 20,000 draws with the probabilities from
`dskellam`, pooling the sparse tails so every cell expects at least 5
observations.

``` r

N <- 20000
lambda1 <- 1.5; lambda2 <- 0.5
z <- rskellam(N, lambda1, lambda2)
c(mean = mean(z), expected = lambda1 - lambda2)
#>     mean expected 
#>   1.0044   1.0000
c(variance = var(z), expected = lambda1 + lambda2)
#> variance expected 
#> 2.013481 2.000000

lo <- qskellam(0.001, lambda1, lambda2); hi <- qskellam(0.999, lambda1, lambda2)
cells <- lo:hi
prob <- dskellam(cells, lambda1, lambda2)
prob[1] <- pskellam(lo, lambda1, lambda2)                                 # all x <= lo
prob[length(prob)] <- pskellam(hi - 1, lambda1, lambda2, lower.tail = FALSE)  # all x >= hi
obs <- table(factor(pmin(pmax(z, lo), hi), levels = cells))
chisq.test(obs, p = prob)
#> 
#>  Chi-squared test for given probabilities
#> 
#> data:  obs
#> X-squared = 3.2351, df = 9, p-value = 0.9542
```

``` r

plot(cells, as.vector(obs), type = "h", lwd = 6, col = "grey70",
     xlab = "x", ylab = "Frequency", main = "")
points(cells, N * prob, pch = 19, col = "#0072B2")
legend("topright", c("rskellam", "N * dskellam"), lwd = c(6, NA), pch = c(NA, 19),
       col = c("grey70", "#0072B2"), bty = "n")
```

![Observed frequencies of 20,000 draws (bars) and the expected
frequencies from dskellam
(points).](validate_skellam_files/figure-html/random-plot-1.png)

Observed frequencies of 20,000 draws (bars) and the expected frequencies
from dskellam (points).

## 4. Estimation

### Maximum likelihood for a sample

``` r

x <- rpois(2000, 12) - rpois(2000, 8)
fit <- skellam.mle(x)
fit$param
#>      mu 1      mu 2 
#> 12.065749  8.159249
```

At the maximum the fitted mean equals the sample mean, and the estimates
are close to the true values (12, 8). Samples whose variance does not
exceed the absolute mean cannot come from a Skellam with two positive
rates; the estimate then puts one rate at zero instead of returning a
negative rate:

``` r

skellam.mle(rep(3, 50))$param
#> mu 1 mu 2 
#>    3    0
skellam.mle(c(rep(2, 45), rep(3, 5)))$param
#> mu 1 mu 2 
#>  2.1  0.0
```

### Regression

`skellam.reg(y, x)` lets both rates depend on covariates through a log
link. For \\y = Y_1 - Y_2\\, `param1` describes the rate of the counts
that are added (\\\lambda_1\\) and `param2` the rate of the counts that
are subtracted (\\\lambda_2\\).

``` r

n <- 1000
covariate <- rnorm(n)
y <- rpois(n, exp(1 + 0.5 * covariate)) - rpois(n, exp(-0.5 + 0.3 * covariate))
reg <- skellam.reg(y, covariate)
round(reg$param1, 3)   # true values 1 and 0.5
#>             Estimate Std. Error Wald value p-value
#> (Intercept)    0.982      0.031     31.615       0
#> x              0.520      0.026     19.876       0
round(reg$param2, 3)   # true values -0.5 and 0.3
#>             Estimate Std. Error Wald value p-value
#> (Intercept)   -0.550      0.118     -4.650   0.000
#> x              0.329      0.110      2.975   0.003
```

## Getting help

If you find a bug, please file an issue with a minimal reproducible
example on [GitHub](https://github.com/monty-se/skellam/issues).
