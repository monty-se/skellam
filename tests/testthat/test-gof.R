test_that("cells have the minimum expected count and cover all observations", {
  set.seed(1)
  x <- rpois(500, 4) - rpois(500, 2)
  t <- skellam.test(x, method = "asymptotic")
  expect_s3_class(t, "htest")
  expect_true(all(t$expected >= 5))
  expect_equal(sum(t$expected), 500)
  expect_equal(sum(t$observed), 500)
  expect_identical(t$cells$lower[1], -Inf); expect_identical(t$cells$upper[nrow(t$cells)], Inf)
  expect_equal(unname(t$statistic), sum((t$observed - t$expected)^2 / t$expected))
  expect_equal(unname(t$parameter), nrow(t$cells) - 3)
})

test_that("given rates use k - 1 degrees of freedom; the LR statistic is G", {
  set.seed(2)
  x <- rpois(400, 4) - rpois(400, 2)
  t <- skellam.test(x, lambda1 = 4, lambda2 = 2, method = "asymptotic")
  expect_equal(unname(t$parameter), length(t$observed) - 1)
  expect_equal(unname(t$estimate), c(4, 2))
  g <- skellam.test(x, statistic = "lr", method = "asymptotic")
  nz <- g$observed > 0
  expect_equal(unname(g$statistic), 2 * sum(g$observed[nz] * log(g$observed[nz] / g$expected[nz])))
  expect_named(g$statistic, "G")
  expect_error(skellam.test(x, lambda1 = 4), "both")
})

test_that("size under the null is close to the nominal level (rates given)", {
  set.seed(3)
  p <- replicate(300, skellam.test(rpois(300, 4) - rpois(300, 2), lambda1 = 4, lambda2 = 2,
                                   method = "asymptotic")$p.value)
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.04)
})

test_that("bootstrap p-value: valid range, reproducible with a seed, and detects extra zeros", {
  set.seed(4)
  x <- rpois(400, 4) - rpois(400, 2)
  set.seed(10); a <- skellam.test(x, B = 99)
  set.seed(10); b <- skellam.test(x, B = 99)
  expect_identical(a$p.value, b$p.value)
  expect_gte(a$p.value, 1 / 100); expect_lte(a$p.value, 1)
  expect_identical(unname(a$parameter), 99L)
  z <- x; z[sample(400, 60)] <- 0
  set.seed(11)
  expect_lt(skellam.test(z, B = 99)$p.value, 0.05)
})
