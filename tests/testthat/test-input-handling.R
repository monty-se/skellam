test_that("missing values propagate instead of erroring  [N1]", {
  expect_identical(dskellam(NA, 3, 2), NA_real_)
  expect_identical(dskellam(NaN, 3, 2), NaN)
  out <- dskellam(c(0, NA, 1), 3, 2)
  expect_equal(out[c(1, 3)], dskellam(c(0, 1), 3, 2))
  expect_true(is.na(out[2]))
  expect_true(is.na(dskellam(0, NA, 2)))
})

test_that("non-integer x has density 0 with a warning  [N11]", {
  expect_warning(d <- dskellam(2.5, 3, 2), "non-integer")
  expect_identical(d, 0)
  expect_warning(d <- dskellam(-2.7, 3, 2, log = TRUE), "non-integer")
  expect_identical(d, -Inf)
  # values within floating-point noise of an integer count as that integer
  expect_identical(dskellam(3 - 1e-12, 3, 2), dskellam(3, 3, 2))
  expect_silent(dskellam(0.1 * 3 * 10, 3, 2))
})

test_that("infinite x has density 0", {
  expect_identical(dskellam(c(-Inf, Inf), 3, 2), c(0, 0))
})

test_that("invalid rates give NaN with a warning", {
  expect_warning(d <- dskellam(0:1, c(-1, 3), 2), "NaNs produced")
  expect_true(is.nan(d[1]))
  expect_equal(d[2], dskellam(1, 3, 2))
})

test_that("recycling warns when lengths are not multiples  [N13]", {
  expect_warning(dskellam(1:3, 1:2, 1), "multiple")
  expect_silent(dskellam(1:4, 1:2, 1))
  expect_silent(dskellam(1:6, 1:2, 1:3))
  expect_identical(dskellam(numeric(0), 3, 2), numeric(0))
})

test_that("names of x are kept  [A10]", {
  expect_named(dskellam(c(a = 0, b = 1), 3, 2), c("a", "b"))
})
