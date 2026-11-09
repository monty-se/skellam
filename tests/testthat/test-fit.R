test_that("skellam.fit recovers the rates, with standard errors that match the observed information", {
  set.seed(11)
  x <- rpois(5000, 10) - rpois(5000, 6)
  f <- skellam.fit(x)
  expect_s3_class(f, "skellam_fit")
  expect_named(coef(f), c("lambda1", "lambda2"))
  expect_equal(unname(coef(f)[1] - coef(f)[2]), mean(x), tolerance = 1e-8)
  se <- sqrt(diag(vcov(f)))
  expect_lt(max(abs(coef(f) - c(10, 6)) / se), 4)
  # covariance = inverse of a numerical Hessian of minus the log-likelihood
  nll <- function(p) -sum(dskellam(x, p[1], p[2], log = TRUE))
  p0 <- unname(coef(f)); h <- 1e-3
  H <- matrix(0, 2, 2)
  for (i in 1:2) for (j in 1:2) {
    ei <- replace(numeric(2), i, h); ej <- replace(numeric(2), j, h)
    H[i, j] <- (nll(p0 + ei + ej) - nll(p0 + ei - ej) - nll(p0 - ei + ej) + nll(p0 - ei - ej)) / (4 * h^2)
  }
  expect_equal(unname(vcov(f)), solve(H), tolerance = 1e-4)
})

test_that("extractor methods agree with the definitions", {
  set.seed(3)
  x <- rpois(400, 4) - rpois(400, 2)
  f <- skellam.fit(x)
  ll <- sum(dskellam(x, coef(f)[[1]], coef(f)[[2]], log = TRUE))
  expect_equal(as.numeric(logLik(f)), ll)
  expect_identical(attr(logLik(f), "df"), 2L)
  expect_identical(nobs(f), 400L)
  expect_equal(AIC(f), -2 * ll + 4)
  expect_equal(BIC(f), -2 * ll + 2 * log(400))
  ci <- confint(f)
  expect_identical(dim(ci), c(2L, 2L))
  expect_true(all(ci[, 1] < coef(f) & coef(f) < ci[, 2]))
  expect_identical(rownames(confint(f, "lambda2")), "lambda2")
  expect_equal(predict(f, x = -3:3), dskellam(-3:3, coef(f)[[1]], coef(f)[[2]]))
  expect_equal(predict(f, x = 0, type = "cdf"), pskellam(0, coef(f)[[1]], coef(f)[[2]]))
  expect_equal(predict(f, type = "quantile", p = 0.5), qskellam(0.5, coef(f)[[1]], coef(f)[[2]]))
  s <- simulate(f, nsim = 3, seed = 1)
  expect_identical(dim(s), c(400L, 3L))
  expect_type(s[[1]], "integer")
})

test_that("print and summary run and show the essentials", {
  set.seed(4)
  f <- skellam.fit(rpois(300, 5) - rpois(300, 3))
  expect_output(print(f), "maximum likelihood")
  s <- summary(f)
  expect_s3_class(s, "summary.skellam_fit")
  expect_output(print(s), "Sample and fitted moments")
  expect_equal(unname(s$moments["mean", "fitted"]), unname(coef(f)[1] - coef(f)[2]))
})

test_that("boundary estimates are flagged, with NA standard errors  [E3]", {
  expect_warning(f <- skellam.fit(rep(3, 50)), "boundary")
  expect_equal(unname(coef(f)), c(3, 0))
  expect_true(f$boundary)
  expect_true(all(is.na(vcov(f))))
  expect_true(all(is.na(confint(f))))
  expect_output(print(f), "boundary")
  expect_warning(f0 <- skellam.fit(rep(0, 20)), "boundary")
  expect_equal(as.numeric(logLik(f0)), 0)
})

test_that("method of moments", {
  set.seed(5)
  x <- rpois(2000, 7) - rpois(2000, 2)
  f <- skellam.fit(x, method = "moments")
  expect_equal(unname(coef(f)), c((var(x) + mean(x)) / 2, (var(x) - mean(x)) / 2))
  expect_true(all(is.na(vcov(f))))
  # truncated at zero for under-dispersed samples
  g <- skellam.fit(c(rep(2, 45), rep(3, 5)), method = "moments")
  expect_equal(unname(coef(g)), c(2.1, 0))
})

test_that("input checks, and na.rm", {
  expect_error(skellam.fit(c(1, NA, 2)), "na.rm")
  set.seed(6); x <- rpois(50, 3) - rpois(50, 2)
  expect_equal(coef(skellam.fit(c(x, NA), na.rm = TRUE)), coef(skellam.fit(x)))
  expect_error(skellam.fit(c(1.5, 2)), "integers")
  expect_error(skellam.fit(3), "two observations")
})

test_that("skellam.mle is deprecated: same estimates, old return value, warns once per session", {
  set.seed(7)
  x <- rpois(500, 6) - rpois(500, 4)
  state <- get(".skellam_state", envir = asNamespace("skellam")); state$warned <- character()
  expect_warning(m <- skellam.mle(x), "deprecated")
  expect_silent(skellam.mle(x))                                # second call: no repeat
  expect_named(m, c("iters", "loglik", "param"))
  expect_named(m$param, c("mu 1", "mu 2"))
  f <- suppressWarnings(skellam.fit(x))
  expect_equal(unname(m$param), unname(coef(f)))
  expect_equal(m$loglik, as.numeric(logLik(f)))
})
