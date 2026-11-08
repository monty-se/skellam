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
