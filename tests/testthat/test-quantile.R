# For a discrete distribution the quantile is defined by an inequality:
#   lower tail: smallest integer x with P(X <= x) >= p
#   upper tail: smallest integer x with P(X >  x) <= p
# so feeding a cdf value back in must return the same x ("round trip").

roundtrip <- function(x, l1, l2, lower.tail = TRUE, log.p = FALSE) {
  qskellam(pskellam(x, l1, l2, lower.tail = lower.tail, log.p = log.p),
           l1, l2, lower.tail = lower.tail, log.p = log.p)
}

test_that("round trip holds on both tails and both scales  [N3]", {
  for (L in list(c(3, 2), c(0.5, 4), c(1e3, 1e3), c(1e4, 1))) {
    x <- support_grid(L[1], L[2], w = 5)
    for (lt in c(TRUE, FALSE)) for (lg in c(FALSE, TRUE)) {
      expect_identical(roundtrip(x, L[1], L[2], lt, lg), as.numeric(x),
                       label = sprintf("lambda=(%g,%g) lower.tail=%s log.p=%s", L[1], L[2], lt, lg))
    }
  }
})

test_that("upper and lower tails agree: q(p, upper) == q(1 - p, lower)", {
  p <- c(0.01, 0.1, 0.3, 0.77)
  expect_identical(qskellam(p, 3, 2, lower.tail = FALSE), qskellam(1 - p, 3, 2))
})
