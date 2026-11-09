## --------------------------------------------------------------------------
## skellam.fit.R -- fitting the Skellam distribution to a sample
##
## Public functions:
##   skellam.fit()   fit (lambda1, lambda2) by maximum likelihood or moments;
##                   returns a "skellam_fit" object
##   S3 methods for "skellam_fit": print, summary, coef, vcov, logLik, nobs,
##                   confint, predict, simulate
##
## Internal functions:
##   .skellam_profile_mle()  the profile-likelihood maximisation shared with
##                           skellam.mle()
##   .skellam_check_sample() input checks shared with skellam.mle()
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

# Validate a sample of differences; return it rounded, without NAs if na.rm.
.skellam_check_sample <- function(x, na.rm = FALSE) {
  if (!is.numeric(x)) stop("'x' must be a numeric vector of integers", call. = FALSE)
  if (anyNA(x)) {
    if (!na.rm) stop("'x' contains missing values; use na.rm = TRUE or remove them", call. = FALSE)
    x <- x[!is.na(x)]
  }
  if (any(!is.finite(x)) || any(.skellam_nonint(x)))
    stop("'x' must contain finite integers (differences of counts)", call. = FALSE)
  if (length(x) < 2L) stop("at least two observations are needed", call. = FALSE)
  round(x)
}

# Profile maximum likelihood. At the MLE lambda1 - lambda2 = mean(x), so the
# log-likelihood is maximised over lambda2 alone, on the admissible range
# lambda2 >= max(0, -mean(x)) that keeps both rates non-negative.
.skellam_profile_mle <- function(x) {
  mx <- mean(x)
  lo <- max(0, -mx)
  evals <- 0L
  prof <- function(theta) {
    evals <<- evals + 1L
    l2 <- lo + theta
    sum(dskellam(x, l2 + mx, l2, log = TRUE))
  }
  upper <- 10 * (stats::var(x) + abs(mx)) + 1
  opt <- stats::optimize(prof, c(0, upper), maximum = TRUE, tol = 1e-10 * (1 + upper))
  # optimize() never evaluates the end points; the MLE can sit on theta = 0
  at0 <- prof(0)
  if (at0 >= opt$objective) {
    theta <- 0; loglik <- at0
  } else {
    theta <- opt$maximum; loglik <- opt$objective
  }
  l2 <- lo + theta
  list(lambda = c(lambda1 = l2 + mx, lambda2 = l2), loglik = loglik, evals = evals)
}

# Gradient of the log-likelihood in (lambda1, lambda2), both > 0.
#   d log f / d lambda1 = (z/(2 lambda1)) (r + nu/z) - 1 + x/(2 lambda1)
#   d log f / d lambda2 = (z/(2 lambda2)) (r + nu/z) - 1 - x/(2 lambda2)
# with nu = |x|, z = 2 sqrt(lambda1 lambda2), r = I_{nu+1}(z) / I_nu(z).
.skellam_score <- function(lambda, x) {
  l1 <- lambda[1]; l2 <- lambda[2]
  nu <- abs(x)
  z  <- 2 * sqrt(l1) * sqrt(l2)
  zn <- rep(z, length(x))
  r  <- exp(.log_besselI_scaled(zn, nu + 1) - .log_besselI_scaled(zn, nu))
  a  <- 0.5 * z * r + 0.5 * nu
  c(sum(a / l1 - 1 + x / (2 * l1)), sum(a / l2 - 1 - x / (2 * l2)))
}

#' Fit the Skellam Distribution to a Sample
#'
#' Estimates the rates \eqn{(\lambda_1, \lambda_2)}{(lambda1, lambda2)} of a
#' Skellam distribution from a sample of differences of counts, by maximum
#' likelihood (default) or by the method of moments, and returns a fitted-model
#' object that works with the usual extractor functions.
#'
#' @param x A numeric vector of integers (differences of counts; any sign).
#' @param method \code{"mle"} (default) for maximum likelihood, or
#'   \code{"moments"} for the method-of-moments estimates
#'   \eqn{\lambda_1 = (s^2 + \bar x)/2}{lambda1 = (var + mean)/2},
#'   \eqn{\lambda_2 = (s^2 - \bar x)/2}{lambda2 = (var - mean)/2}, truncated at 0.
#' @param na.rm Logical; if \code{TRUE}, missing values are removed before
#'   fitting. If \code{FALSE} (default), missing values give an error.
#'
#' @details
#' \strong{Maximum likelihood.} At the maximum the fitted mean equals the sample
#' mean, \eqn{\hat\lambda_1 - \hat\lambda_2 = \bar{x}}{lambda1 - lambda2 = mean(x)},
#' so the log-likelihood is maximised over \eqn{\lambda_2}{lambda2} alone with
#' \code{\link[stats]{optimize}}, on the range that keeps both rates non-negative.
#' The covariance matrix is the inverse of the observed information, computed
#' from the analytic score.
#'
#' \strong{Boundary estimates.} A Skellam with two positive rates has variance
#' \eqn{\lambda_1 + \lambda_2 \ge |\lambda_1 - \lambda_2|}{lambda1 + lambda2 >= |lambda1 - lambda2|}.
#' When the sample variance does not exceed the absolute sample mean, the
#' maximum lies on the boundary, with one rate equal to 0 (a Poisson or a
#' reflected Poisson): a sample of zeros gives \eqn{(0, 0)}, a constant sample
#' \eqn{c > 0} gives \eqn{(c, 0)}. Standard errors do not exist there; the
#' covariance matrix is then \code{NA} and the fit is flagged with
#' \code{boundary = TRUE} and a warning.
#'
#' \strong{Moments.} The moment estimates are returned without a covariance
#' matrix and serve mainly as a quick check or a starting point.
#'
#' @return An object of class \code{"skellam_fit"}, a list with components
#' \describe{
#'   \item{estimate}{Named vector \code{c(lambda1 = , lambda2 = )}.}
#'   \item{vcov}{2 x 2 covariance matrix of the estimates (\code{NA} on the
#'     boundary and for the moment method).}
#'   \item{loglik}{Log-likelihood at the estimates.}
#'   \item{nobs}{Number of observations used.}
#'   \item{method}{\code{"mle"} or \code{"moments"}.}
#'   \item{boundary}{\code{TRUE} if one rate is estimated as 0.}
#'   \item{evals}{Number of log-likelihood evaluations (maximum likelihood).}
#'   \item{x}{The data used.}
#'   \item{call}{The matched call.}
#' }
#' Methods are available for \code{print}, \code{summary}, \code{coef},
#' \code{vcov}, \code{logLik} (hence \code{AIC} and \code{BIC}), \code{nobs},
#' \code{confint}, \code{predict} and \code{simulate}.
#'
#' @seealso \code{\link{skellam}} for the distribution functions,
#'   \code{\link{skellam.reg}} for regression on covariates.
#'
#' @references
#' Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football outcomes:
#' using the Skellam's distribution for the goal difference. \emph{IMA Journal of
#' Management Mathematics} \bold{20}(2), 133-145.
#'
#' @examples
#' set.seed(1)
#' x <- rpois(500, 6) - rpois(500, 4)
#' fit <- skellam.fit(x)
#' fit
#' summary(fit)
#' coef(fit)
#' confint(fit)
#' AIC(fit)
#'
#' # probabilities and quantities from the fitted distribution
#' predict(fit, x = -2:2)
#' predict(fit, x = 0, type = "cdf")
#'
#' # the boundary: a constant sample is best fitted by a Poisson
#' skellam.fit(rep(3, 50))
#'
#' @importFrom stats coef vcov logLik nobs confint predict simulate
#' @export
skellam.fit <- function(x, method = c("mle", "moments"), na.rm = FALSE) {
  cl <- match.call()
  method <- match.arg(method)
  x <- .skellam_check_sample(x, na.rm = na.rm)
  n <- length(x)
  m <- mean(x); v <- stats::var(x)
  V <- matrix(NA_real_, 2, 2, dimnames = list(c("lambda1", "lambda2"), c("lambda1", "lambda2")))
  evals <- NA_integer_

  if (method == "mle") {
    est <- .skellam_profile_mle(x)
    lambda <- est$lambda
    loglik <- est$loglik
    evals <- est$evals
  } else {
    l2 <- max((v - m) / 2, 0, -m)
    lambda <- c(lambda1 = l2 + m, lambda2 = l2)
    loglik <- sum(dskellam(x, lambda[1], lambda[2], log = TRUE))
  }
  boundary <- any(lambda <= 0)

  if (method == "mle") {
    if (boundary) {
      warning("the estimate lies on the boundary (a rate is 0): the sample variance does not ",
              "exceed the absolute mean; standard errors are not available", call. = FALSE)
    } else {
      nll  <- function(p) -sum(dskellam(x, p[1], p[2], log = TRUE))
      grad <- function(p) -.skellam_score(p, x)
      H <- stats::optimHess(unname(lambda), nll, grad)
      inv <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(inv) || any(diag(inv) <= 0)) {
        warning("the observed information is not positive definite; standard errors are not available",
                call. = FALSE)
      } else {
        V[] <- (inv + t(inv)) / 2
      }
    }
  }

  structure(list(estimate = lambda, vcov = V, loglik = loglik, nobs = n,
                 method = method, boundary = boundary, evals = evals,
                 x = x, call = cl),
            class = "skellam_fit")
}

#' @rdname skellam.fit
#' @param object,digits,...,level,parm,type,nsim,seed Arguments of the
#'   methods: \code{object} is a \code{"skellam_fit"} object; \code{parm} and
#'   \code{level} as in \code{\link[stats]{confint}}; \code{type} in
#'   \code{predict} is \code{"pmf"} (default), \code{"cdf"} or
#'   \code{"quantile"}, evaluated at \code{x} (quantiles at probabilities
#'   \code{p}); \code{nsim} and \code{seed} as in \code{\link[stats]{simulate}}.
#' @export
print.skellam_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Skellam distribution fitted by",
      if (x$method == "mle") "maximum likelihood" else "the method of moments",
      sprintf("(n = %d)\n\n", x$nobs))
  se <- sqrt(diag(x$vcov))
  tab <- cbind(Estimate = x$estimate, `Std. Error` = se)
  print(tab, digits = digits, na.print = "NA")
  cat(sprintf("\nlog-likelihood %s, AIC %s\n", format(x$loglik, digits = digits + 3),
              format(-2 * x$loglik + 4, digits = digits + 3)))
  if (x$boundary) cat("Estimate on the boundary: one rate is 0 (standard errors not available).\n")
  invisible(x)
}

#' @rdname skellam.fit
#' @export
summary.skellam_fit <- function(object, level = 0.95, ...) {
  x <- object$x
  l1 <- object$estimate[[1]]; l2 <- object$estimate[[2]]
  m3 <- function(z) mean((z - mean(z))^3) / mean((z - mean(z))^2)^1.5
  m4 <- function(z) mean((z - mean(z))^4) / mean((z - mean(z))^2)^2 - 3
  mom <- rbind(
    mean            = c(sample = mean(x),  fitted = l1 - l2),
    variance        = c(sample = stats::var(x), fitted = l1 + l2),
    skewness        = c(sample = m3(x), fitted = if (l1 + l2 > 0) (l1 - l2) / (l1 + l2)^1.5 else NA),
    `excess kurtosis` = c(sample = m4(x), fitted = if (l1 + l2 > 0) 1 / (l1 + l2) else NA))
  ci <- confint(object, level = level)
  coefs <- cbind(Estimate = object$estimate, `Std. Error` = sqrt(diag(object$vcov)), ci)
  structure(list(coefficients = coefs, moments = mom, loglik = object$loglik,
                 aic = -2 * object$loglik + 4, bic = -2 * object$loglik + 2 * log(object$nobs),
                 nobs = object$nobs, method = object$method, boundary = object$boundary,
                 dispersion = stats::var(x) - abs(mean(x)), call = object$call),
            class = "summary.skellam_fit")
}

#' @export
print.summary.skellam_fit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Skellam distribution fitted by",
      if (x$method == "mle") "maximum likelihood" else "the method of moments",
      sprintf("(n = %d)\n\n", x$nobs))
  print(x$coefficients, digits = digits, na.print = "NA")
  cat("\nSample and fitted moments:\n")
  print(x$moments, digits = digits)
  cat(sprintf("\nvar(x) - |mean(x)| = %s (must be >= 0 for two positive rates)\n",
              format(x$dispersion, digits = digits)))
  cat(sprintf("log-likelihood %s, AIC %s, BIC %s\n", format(x$loglik, digits = digits + 3),
              format(x$aic, digits = digits + 3), format(x$bic, digits = digits + 3)))
  if (x$boundary) cat("Estimate on the boundary: one rate is 0 (standard errors not available).\n")
  invisible(x)
}

#' @rdname skellam.fit
#' @export
coef.skellam_fit <- function(object, ...) object$estimate

#' @rdname skellam.fit
#' @export
vcov.skellam_fit <- function(object, ...) object$vcov

#' @rdname skellam.fit
#' @export
logLik.skellam_fit <- function(object, ...) {
  structure(object$loglik, df = 2L, nobs = object$nobs, class = "logLik")
}

#' @rdname skellam.fit
#' @export
nobs.skellam_fit <- function(object, ...) object$nobs

#' @rdname skellam.fit
#' @details Confidence intervals are Wald intervals on the log scale,
#'   \eqn{\exp(\log\hat\lambda \pm z\,\mathrm{se}/\hat\lambda)}{exp(log(lambda) +- z se / lambda)},
#'   which keeps them positive; they are \code{NA} when the covariance matrix is.
#' @export
confint.skellam_fit <- function(object, parm, level = 0.95, ...) {
  est <- object$estimate
  se <- sqrt(diag(object$vcov))
  z <- stats::qnorm(1 - (1 - level) / 2)
  ci <- cbind(est * exp(-z * se / est), est * exp(z * se / est))
  a <- (1 - level) / 2
  colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%")
  rownames(ci) <- names(est)
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}

#' @rdname skellam.fit
#' @param p Probabilities, for \code{predict(type = "quantile")}.
#' @export
predict.skellam_fit <- function(object, x, type = c("pmf", "cdf", "quantile"), p, ...) {
  type <- match.arg(type)
  l <- object$estimate
  switch(type,
         pmf      = dskellam(x, l[[1]], l[[2]]),
         cdf      = pskellam(x, l[[1]], l[[2]]),
         quantile = qskellam(p, l[[1]], l[[2]]))
}

#' @rdname skellam.fit
#' @export
simulate.skellam_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  l <- object$estimate
  n <- object$nobs
  out <- as.data.frame(matrix(rskellam(n * nsim, l[[1]], l[[2]]), n, nsim))
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}
