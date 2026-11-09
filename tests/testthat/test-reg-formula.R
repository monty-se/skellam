sim_goals <- function(n = 600, seed = 1) {
  set.seed(seed)
  d <- data.frame(attack = rnorm(n), defence = rnorm(n),
                  venue = factor(sample(c("home", "neutral", "away"), n, TRUE)),
                  expo = runif(n, 1, 2))
  d$y <- rpois(n, exp(0.4 + 0.3 * d$attack)) - rpois(n, exp(0.1 - 0.2 * d$defence))
  d
}

test_that("one-part formula reproduces the original (y, x) interface", {
  d <- sim_goals()
  a <- skellam.reg(y ~ attack, data = d)
  b <- skellam.reg(d$y, d$attack)
  expect_s3_class(a, "skellam_reg"); expect_s3_class(b, "skellam_reg")
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-8)
  expect_equal(a$loglik, b$loglik, tolerance = 1e-10)
  expect_identical(colnames(b$param1), c("Estimate", "Std. Error", "Wald value", "p-value"))
  expect_equal(unname(skellam.reg(y = d$y, x = d$attack)$loglik), b$loglik)
})

test_that("two-part formula: different covariates for the two rates recover the truth", {
  d <- sim_goals(3000, seed = 2)
  m <- skellam.reg(y ~ attack | defence, data = d)
  expect_named(coef(m, "lambda1"), c("(Intercept)", "attack"))
  expect_named(coef(m, "lambda2"), c("(Intercept)", "defence"))
  est <- coef(m); se <- sqrt(diag(vcov(m)))
  expect_lt(max(abs(est - c(0.4, 0.3, 0.1, -0.2)) / se), 4)
  expect_identical(names(coef(m))[1], "lambda1:(Intercept)")
  expect_identical(dim(vcov(m)), c(4L, 4L))
  expect_equal(attr(logLik(m), "df"), 4)
  expect_equal(AIC(m), -2 * m$loglik + 8)
  expect_identical(nobs(m), 3000L)
  expect_identical(dim(confint(m)), c(4L, 2L))
})

test_that("log-likelihood equals the sum of dskellam at the fitted rates", {
  d <- sim_goals()
  m <- skellam.reg(y ~ attack + venue | defence, data = d)
  expect_equal(m$loglik, sum(dskellam(d$y, m$fitted.values$lambda1, m$fitted.values$lambda2, log = TRUE)),
               tolerance = 1e-10)
})

test_that("factors, offsets and missing values", {
  d <- sim_goals()
  m <- skellam.reg(y ~ attack + venue | defence, data = d)
  expect_named(coef(m, "lambda1"), c("(Intercept)", "attack", "venuehome", "venueneutral"))
  # an offset shifts log(lambda1) exactly
  mo <- skellam.reg(y ~ attack + offset(log(expo)) | defence, data = d)
  nd <- data.frame(attack = 0.5, defence = 0, expo = 2)
  expect_equal(unname(predict(mo, nd, type = "lambda1")),
               exp(coef(mo, "lambda1")[[1]] + 0.5 * coef(mo, "lambda1")[[2]] + log(2)))
  d$attack[3] <- NA; d$defence[7] <- NA
  expect_warning(m2 <- skellam.reg(y ~ attack | defence, data = d), "2 observation")
  expect_identical(nobs(m2), nrow(d) - 2L)
})

test_that("predict: response, rates, probabilities on new data", {
  d <- sim_goals()
  m <- skellam.reg(y ~ attack | defence, data = d)
  expect_equal(predict(m), fitted(m))
  nd <- data.frame(attack = c(-1, 0, 1), defence = c(0, 0, 1))
  r <- predict(m, nd, type = "rates")
  expect_equal(unname(r[, "lambda1"]), exp(coef(m, "lambda1")[[1]] + coef(m, "lambda1")[[2]] * nd$attack))
  expect_equal(predict(m, nd), unname(r[, 1] - r[, 2]))
  pr <- predict(m, nd, type = "prob", at = -40:40)
  expect_identical(dim(pr), c(3L, 81L))
  expect_equal(rowSums(pr), rep(1, 3), tolerance = 1e-12)
  expect_equal(predict(m, nd, type = "cdf", at = 0)[, 1], pskellam(0, r[, 1], r[, 2]))
  expect_error(predict(m, nd, type = "prob"), "at")
  old <- skellam.reg(d$y, d$attack)
  expect_error(predict(old, nd), "formula interface")
})

test_that("residuals: response, Pearson and quantile residuals", {
  d <- sim_goals(2000, seed = 3)
  m <- skellam.reg(y ~ attack | defence, data = d)
  l1 <- m$fitted.values$lambda1; l2 <- m$fitted.values$lambda2
  expect_equal(residuals(m), d$y - (l1 - l2))
  expect_equal(residuals(m, "pearson"), (d$y - (l1 - l2)) / sqrt(l1 + l2))
  set.seed(4); q <- residuals(m, "quantile")
  expect_lt(abs(mean(q)), 0.1); expect_lt(abs(sd(q) - 1), 0.1)
})

test_that("anova gives the likelihood-ratio test for nested models; update refits", {
  d <- sim_goals()
  m1 <- skellam.reg(y ~ attack | defence, data = d)
  m0 <- skellam.reg(y ~ 1 | 1, data = d)
  a <- anova(m0, m1)
  expect_s3_class(a, "anova")
  expect_equal(a$Chisq[2], 2 * (m1$loglik - m0$loglik))
  expect_equal(a$Df[2], 2)
  expect_equal(anova(m1, m0)$Chisq[2], a$Chisq[2])     # order does not matter
  u <- update(m1, y ~ attack | defence + venue)
  expect_named(coef(u, "lambda2"), c("(Intercept)", "defence", "venuehome", "venueneutral"))
  expect_error(anova(m1), "two")
})

test_that("simulate, print and summary", {
  d <- sim_goals()
  m <- skellam.reg(y ~ attack | defence, data = d)
  s <- simulate(m, nsim = 2, seed = 1)
  expect_identical(dim(s), c(nrow(d), 2L))
  expect_type(s[[1]], "integer")
  expect_output(print(m), "counts that are added")
  expect_output(print(summary(m)), "Pr\\(>\\|z\\|\\)")
})

test_that("formula errors are clear", {
  d <- sim_goals()
  expect_error(skellam.reg(~ attack, data = d), "response")
  expect_error(skellam.reg(y ~ attack | defence | venue, data = d), "only one")
})
