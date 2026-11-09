# Goodness-of-Fit Test for the Skellam Distribution

Tests whether a sample of differences of counts comes from a Skellam
distribution, with rates either given or estimated by maximum
likelihood.

## Usage

``` r
skellam.test(
  x,
  lambda1 = NULL,
  lambda2 = NULL,
  statistic = c("pearson", "lr"),
  method = c("bootstrap", "asymptotic"),
  B = 999,
  min.expected = 5
)
```

## Arguments

- x:

  A numeric vector of integers.

- lambda1, lambda2:

  Rates of the hypothesised distribution. If both are `NULL` (default)
  they are estimated from `x` with
  [`skellam.fit`](https://monty-se.github.io/skellam/reference/skellam.fit.md).

- statistic:

  `"pearson"` (default) for the Pearson chi-squared statistic or `"lr"`
  for the likelihood-ratio (G) statistic.

- method:

  `"bootstrap"` (default) for a parametric-bootstrap p-value, or
  `"asymptotic"` for the chi-squared approximation.

- B:

  Number of bootstrap samples.

- min.expected:

  Minimum expected count per cell; sparse tail cells are merged until
  every cell reaches it.

## Value

An object of class `"htest"` with the usual components (`statistic`,
`parameter`: degrees of freedom or `B`, `p.value`, `method`,
`data.name`, and `estimate`: the rates used), plus `observed`,
`expected` and `cells` (a data frame with the limits of each cell).

## Details

The integers are grouped into cells so that every cell has an expected
count of at least `min.expected` under the hypothesised (or fitted)
distribution: the cells start one per integer, the two end cells take
the tails, and the cell with the smallest expected count is merged with
its smaller neighbour until all reach the minimum. With \\k\\ cells, the
statistic is \\X^2 = \sum (O - E)^2 / E\\ or \\G = 2 \sum O \log(O/E)\\.

**p-values.** When the rates are given, both statistics are
approximately \\\chi^2\_{k-1}\\. When the rates are estimated from the
ungrouped data, the limiting distribution lies between \\\chi^2\_{k-3}\\
and \\\chi^2\_{k-1}\\ (Chernoff and Lehmann, 1954), so no chi-squared
reference is exact: `method = "asymptotic"` uses \\k - 3\\ degrees of
freedom, which tends to reject slightly too often. The default
parametric bootstrap avoids the problem: it simulates `B` samples from
the fitted distribution, re-estimates the rates and rebuilds the cells
for each, and reports \\(1 + \\\\T^\* \ge T\\)/(B + 1)\\. Set a seed for
reproducible bootstrap p-values.

## References

Chernoff, H. and Lehmann, E. L. (1954) The use of maximum likelihood
estimates in \\\chi^2\\ tests for goodness of fit. *Annals of
Mathematical Statistics* **25**(3), 579-586.

## See also

[`skellam.fit`](https://monty-se.github.io/skellam/reference/skellam.fit.md),
[`chisq.test`](https://rdrr.io/r/stats/chisq.test.html)

## Examples

``` r
set.seed(1)
x <- rpois(500, 4) - rpois(500, 2)
skellam.test(x, B = 199)
#> 
#>  Skellam goodness-of-fit test (Pearson, parametric bootstrap)
#> 
#> data:  x
#> X-squared = 2.3808, B = 199, p-value = 1
#> sample estimates:
#>  lambda1  lambda2 
#> 3.847458 1.877458 
#> 

# a sample with extra zeros is rejected
z <- x; z[sample(500, 60)] <- 0
skellam.test(z, B = 199)
#> 
#>  Skellam goodness-of-fit test (Pearson, parametric bootstrap)
#> 
#> data:  z
#> X-squared = 41.32, B = 199, p-value = 0.005
#> sample estimates:
#>  lambda1  lambda2 
#> 3.556766 1.826766 
#> 

# rates given in advance: asymptotic chi-squared with k - 1 df
skellam.test(x, lambda1 = 4, lambda2 = 2, method = "asymptotic")
#> 
#>  Skellam goodness-of-fit test (Pearson, asymptotic chi-squared)
#> 
#> data:  x
#> X-squared = 4.1428, df = 12, p-value = 0.9807
#> sample estimates:
#> lambda1 lambda2 
#>       4       2 
#> 
```
