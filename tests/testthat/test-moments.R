test_that("cumulants and moments agree with sums over the pmf", {
  m <- skellam.moments(3, 2, order = 6)
  x <- -100:100; d <- dskellam(x, 3, 2)
  expect_equal(unname(m$cumulants), 3 + (-1)^(1:6) * 2)
  expect_equal(unname(m$raw), vapply(1:6, function(r) sum(x^r * d), 0), tolerance = 1e-10)
  expect_equal(unname(m$central), vapply(1:6, function(r) sum((x - 1)^r * d), 0), tolerance = 1e-10)
  expect_equal(m$mean, 1); expect_equal(m$variance, 5)
  expect_equal(m$skewness, 1 / 5^1.5); expect_equal(m$kurtosis, 1 / 5)
  expect_output(print(m), "excess kurtosis")
})

test_that("Poisson limit and input checks", {
  m <- skellam.moments(4, 0)
  expect_equal(unname(m$cumulants), rep(4, 4))                 # Poisson: all cumulants = lambda
  expect_error(skellam.moments(c(1, 2), 1), "single")
  expect_error(skellam.moments(-1, 1), "non-negative")
  expect_error(skellam.moments(1, 1, order = 2), "at least 4")
})

test_that("generating functions agree with expectations over the pmf", {
  x <- -100:100; d <- dskellam(x, 3, 2)
  expect_equal(skellam.mgf(0, 3, 2), 1)
  expect_equal(skellam.mgf(0.4, 3, 2), sum(exp(0.4 * x) * d), tolerance = 1e-12)
  expect_equal(skellam.cf(0.3, 3, 2), sum(exp(1i * 0.3 * x) * d), tolerance = 1e-12)
  expect_equal(skellam.pgf(0.7, 3, 2), sum(0.7^x * d), tolerance = 1e-12)
  expect_equal(skellam.pgf(1, 3, 2), 1)
  expect_error(skellam.pgf(0, 3, 2), "non-zero")
})
