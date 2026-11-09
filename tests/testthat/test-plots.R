skip_if_not_installed("ggplot2")

fit <- local({ set.seed(1); suppressWarnings(skellam.fit(rpois(400, 6) - rpois(400, 4))) })
reg <- local({
  set.seed(2); n <- 400
  d <- data.frame(a = rnorm(n), b = rnorm(n))
  d$y <- rpois(n, exp(1 + 0.5 * d$a)) - rpois(n, exp(0.2 - 0.3 * d$b))
  skellam.reg(y ~ a | b, data = d)
})

test_that("every plot type returns a ggplot that builds", {
  for (obj in list(fit, reg)) for (t in c("rootogram", "pmf", "pit", "qq")) {
    p <- plot(obj, type = t)
    expect_s3_class(p, "ggplot")
    expect_silent(ggplot2::ggplot_build(p))
  }
  expect_s3_class(plot(fit, style = "standing"), "ggplot")
  expect_s3_class(ggplot2::autoplot(fit, type = "pit"), "ggplot")
  expect_s3_class(ggplot2::autoplot(reg), "ggplot")
})

test_that("rootogram data: observed counts add up to n, expected counts to about n", {
  dd <- skellam:::.skellam_diag_data(fit$x, coef(fit)[[1]], coef(fit)[[2]])
  expect_identical(sum(dd$counts$observed), 400L)
  expect_equal(sum(dd$counts$expected), 400, tolerance = 2e-3)
  expect_equal(mean(dd$pit$height), 1)                        # a density on [0, 1]
})

test_that("theme and palette", {
  expect_s3_class(theme_skellam(), "theme")
  expect_named(skellam_palette(), c("observed", "fitted", "reference", "secondary", "neutral", "band"))
  expect_s3_class(scale_colour_skellam(), "Scale")
})
