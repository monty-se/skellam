test_that("skellam.mle() does not change options(warn)  [E1]", {
  old <- options(warn = 0); on.exit(options(old))
  set.seed(42); x <- rpois(200, 6) - rpois(200, 4)
  skellam.mle(x)
  expect_equal(getOption("warn"), 0)
})

test_that("skellam.reg() does not change options(warn), even if the fit fails  [E1]", {
  old <- options(warn = 0); on.exit(options(old))
  set.seed(42); y <- rpois(200, 6) - rpois(200, 4)
  try(suppressWarnings(skellam.reg(y, rnorm(200))), silent = TRUE)
  expect_equal(getOption("warn"), 0)
})

test_that("skellam.mle recovers the rates and keeps its return shape", {
  set.seed(11)
  x <- rpois(5000, 10) - rpois(5000, 6)
  f <- skellam.mle(x)
  expect_named(f, c("iters", "loglik", "param"))
  expect_named(f$param, c("mu 1", "mu 2"))
  expect_equal(unname(f$param[1] - f$param[2]), mean(x), tolerance = 1e-8)   # MLE matches the mean
  expect_lt(abs(f$param[[1]] - 10), 3 * sqrt(2 * 16 / 5000) * 4)
  expect_equal(f$loglik, sum(dskellam(x, f$param[[1]], f$param[[2]], log = TRUE)), tolerance = 1e-12)
  # the maximum: nudging lambda2 either way (with lambda1 tracking) lowers the log-likelihood
  ll <- function(d) sum(dskellam(x, f$param[[1]] + d, f$param[[2]] + d, log = TRUE))
  expect_lt(ll(0.01), f$loglik); expect_lt(ll(-0.01), f$loglik)
})

test_that("skellam.mle returns non-negative rates for under-dispersed samples  [E2]", {
  f <- skellam.mle(c(rep(2, 45), rep(3, 5)))     # var 0.09 < mean 2.1
  expect_true(all(f$param >= 0))
  expect_equal(unname(f$param), c(2.1, 0))
})

test_that("skellam.mle handles degenerate samples with the boundary estimate  [E3]", {
  expect_equal(unname(skellam.mle(rep(0, 50))$param), c(0, 0))
  expect_equal(skellam.mle(rep(0, 50))$loglik, 0)
  expect_equal(unname(skellam.mle(rep(3, 50))$param), c(3, 0))
  expect_equal(unname(skellam.mle(rep(-3, 50))$param), c(0, 3))
})

test_that("skellam.mle rejects invalid input with clear messages  [E3]", {
  expect_error(skellam.mle(c(1, NA, 2)), "missing values")
  expect_error(skellam.mle(c(1.5, 2)), "integers")
  expect_error(skellam.mle(3), "two observations")
  expect_error(skellam.mle("a"), "numeric")
})

test_that("skellam.mle works at very large and very small rates", {
  set.seed(9)
  x <- rpois(500, 1e5) - rpois(500, 1e5 - 50)
  f <- skellam.mle(x)
  expect_true(all(is.finite(f$param)) && all(f$param > 0))
  expect_equal(unname(f$param[1] - f$param[2]), mean(x), tolerance = 1e-8)
  x <- rpois(2000, 0.05) - rpois(2000, 0.03)
  expect_true(all(is.finite(skellam.mle(x)$param)))
})
