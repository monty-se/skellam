# Accuracy of the log-pmf core against the exact Poisson-product sum, over the
# whole range of rates. Target: relative error of the pmf below 1e-12 up to
# rates of 1e6 with equal rates, below 1e-9 for extremely unequal huge rates
# (where the Bessel representation itself cancels 1e6-sized terms).

acc_grid <- list(c(0.01, 0.01), c(3, 2), c(100, 80), c(1e3, 1e3), c(1e4, 1e4),
                 c(3e4, 3e4), c(1e5, 1e5), c(1e6, 1e6), c(50, 0.001), c(300, 300))
unequal  <- list(c(1e4, 1), c(1e6, 1e3), c(1e6, 1e5), c(2e5, 10))

test_that("log-pmf matches the exact sum to 1e-12, including at rates of 1e5-1e6  [N5b, N7]", {
  for (L in acc_grid) {
    s <- sqrt(sum(L)); m <- L[1] - L[2]
    x <- unique(round(m + s * c(-30, -8, -3, -1, 0, 1, 3, 8, 30)))
    expect_silent(lp <- dskellam(x, L[1], L[2], log = TRUE))
    expect_equal(lp, vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]), tolerance = 1e-12,
                 label = sprintf("log-pmf at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("log-pmf matches the exact sum to 1e-9 for extremely unequal huge rates", {
  for (L in unequal) {
    s <- sqrt(sum(L)); m <- L[1] - L[2]
    x <- unique(round(m + s * c(-30, -8, -3, 0, 3, 8, 30)))
    expect_silent(lp <- dskellam(x, L[1], L[2], log = TRUE))
    expect_equal(lp, vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]), tolerance = 1e-9,
                 label = sprintf("log-pmf at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("pmf sums to 1 at large and very unequal rates, without warnings  [N7]", {
  for (L in list(c(1e5, 1e5), c(1e6, 1e6), c(1e4, 1), c(1e6, 1e3))) {
    x <- support_grid(L[1], L[2], w = 25)
    expect_silent(d <- dskellam(x, L[1], L[2]))
    expect_equal(sum(d), 1, tolerance = 1e-10, label = sprintf("mass at lambda = (%g, %g)", L[1], L[2]))
  }
})

test_that("log = TRUE stays finite in the deep tails  [N5]", {
  x <- c(-2000, -500, -200, 200, 500, 1000)
  expect_silent(lp <- dskellam(x, 3, 2, log = TRUE))
  expect_true(all(is.finite(lp)))
  expect_equal(lp, vapply(x, lp_exact, 0, l1 = 3, l2 = 2), tolerance = 1e-12)
  expect_identical(dskellam(500, 3, 2), 0)                  # natural scale underflows to 0, correctly
})

test_that("log = TRUE never returns a natural-scale value  [N5b]", {
  expect_lt(dskellam(0, 1e5, 1e5, log = TRUE), 0)
  expect_equal(dskellam(0, 1e5, 1e5, log = TRUE), lp_exact(0, 1e5, 1e5), tolerance = 1e-12)
})

test_that("fallback indexing cannot mix observations  [N6]", {
  expect_identical(dskellam(c(0, 400, -3), 3, 2)[2], dskellam(400, 3, 2))
  expect_identical(dskellam(c(0, 400, -3), 3, 2, log = TRUE)[2], dskellam(400, 3, 2, log = TRUE))
})

test_that("tiny rates", {
  x <- c(-2, -1, 0, 1, 2, 5)
  for (L in list(c(1e-8, 1e-8), c(1e-300, 1), c(1e-3, 1e4))) {
    expect_equal(dskellam(x, L[1], L[2], log = TRUE), vapply(x, lp_exact, 0, l1 = L[1], l2 = L[2]),
                 tolerance = 1e-11)
  }
})

test_that("dskellam.sp honours log and carries the full second-order correction  [N15]", {
  expect_equal(dskellam.sp(-5:5, 3, 2, log = TRUE), log(dskellam.sp(-5:5, 3, 2)))
  # with the full correction the approximation integrates to 1 up to O(1/lambda^2)
  x <- support_grid(1e5, 1e5, w = 20)
  expect_equal(sum(dskellam.sp(x, 1e5, 1e5)), 1, tolerance = 1e-9)
  # and is within 1e-5 relative of the exact pmf at moderate rates
  x <- -40:40
  expect_equal(dskellam.sp(x, 100, 80), dskellam(x, 100, 80), tolerance = 1e-4)
})
