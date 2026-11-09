## --------------------------------------------------------------------------
## skellam.moments.R -- moments, cumulants and generating functions
##
## Public functions:
##   skellam.moments()  mean, variance, skewness, kurtosis, and cumulants,
##                      raw and central moments of any order
##   skellam.mgf()      moment generating function
##   skellam.cf()       characteristic function
##   skellam.pgf()      probability generating function
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' Moments and Generating Functions of the Skellam Distribution
#'
#' Theoretical moments, cumulants and generating functions of the Skellam
#' distribution with rates \eqn{\lambda_1}{lambda1} and \eqn{\lambda_2}{lambda2}.
#'
#' @param lambda1,lambda2 Non-negative rates (single numbers for
#'   \code{skellam.moments}, recycled with \code{t} or \code{z} for the
#'   generating functions); \code{lambda2} defaults to \code{lambda1}.
#' @param order Highest order of the cumulants and moments returned.
#' @param t Real argument of the moment generating and characteristic functions.
#' @param z Non-zero argument of the probability generating function (real or
#'   complex).
#' @param x A \code{"skellam_moments"} object.
#' @param digits Number of significant digits to print.
#' @param ... Unused.
#'
#' @details
#' The cumulant generating function of \eqn{X = Y_1 - Y_2}{X = Y1 - Y2} is
#' \eqn{K(t) = \lambda_1(e^t - 1) + \lambda_2(e^{-t} - 1)}{K(t) = lambda1 (exp(t) - 1) + lambda2 (exp(-t) - 1)},
#' so the cumulants are \eqn{\kappa_r = \lambda_1 + (-1)^r \lambda_2}{kappa_r = lambda1 + (-1)^r lambda2}:
#' the mean is \eqn{\lambda_1 - \lambda_2}{lambda1 - lambda2}, the variance
#' \eqn{\lambda_1 + \lambda_2}{lambda1 + lambda2}, the skewness
#' \eqn{(\lambda_1 - \lambda_2)/(\lambda_1 + \lambda_2)^{3/2}}{(lambda1 - lambda2)/(lambda1 + lambda2)^(3/2)}
#' and the excess kurtosis \eqn{1/(\lambda_1 + \lambda_2)}{1/(lambda1 + lambda2)}.
#' Raw moments follow from the cumulants by the recursion
#' \eqn{\mu'_n = \sum_{k=0}^{n-1} \binom{n-1}{k} \kappa_{n-k} \mu'_k}{mu'_n = sum_k choose(n-1, k) kappa_(n-k) mu'_k},
#' and central moments by the same recursion with \eqn{\kappa_1 = 0}{kappa_1 = 0}.
#'
#' The generating functions are
#' \eqn{M(t) = E[e^{tX}] = \exp\{K(t)\}}{M(t) = exp(K(t))},
#' \eqn{\varphi(t) = E[e^{itX}] = \exp\{\lambda_1(e^{it}-1) + \lambda_2(e^{-it}-1)\}}{phi(t) = exp(lambda1 (exp(it) - 1) + lambda2 (exp(-it) - 1))}
#' and \eqn{G(z) = E[z^X] = \exp\{\lambda_1(z-1) + \lambda_2(1/z-1)\}}{G(z) = exp(lambda1 (z - 1) + lambda2 (1/z - 1))}.
#'
#' Sample moments of data are reported by \code{summary()} of a fit from
#' \code{\link{skellam.fit}}, next to the fitted ones.
#'
#' @return \code{skellam.moments} returns an object of class
#'   \code{"skellam_moments"}: a list with \code{mean}, \code{variance},
#'   \code{sd}, \code{skewness}, \code{kurtosis} (excess), and the vectors
#'   \code{cumulants}, \code{raw} and \code{central} of orders \code{1:order}.
#'   The generating functions return numeric (complex for \code{skellam.cf},
#'   and for \code{skellam.pgf} with complex \code{z}) vectors.
#'
#' @examples
#' skellam.moments(3, 2)
#' skellam.moments(3, 2, order = 6)$raw
#'
#' # the mgf at 0 is 1 and its derivative is the mean
#' skellam.mgf(0, 3, 2)
#' (skellam.mgf(1e-6, 3, 2) - skellam.mgf(-1e-6, 3, 2)) / 2e-6
#'
#' # E[z^X] equals the sum of z^x P(X = x)
#' x <- -60:60
#' c(skellam.pgf(0.7, 3, 2), sum(0.7^x * dskellam(x, 3, 2)))
#'
#' @export
skellam.moments <- function(lambda1, lambda2 = lambda1, order = 4) {
  if (length(lambda1) != 1L || length(lambda2) != 1L)
    stop("'lambda1' and 'lambda2' must be single numbers", call. = FALSE)
  if (!is.finite(lambda1) || !is.finite(lambda2) || lambda1 < 0 || lambda2 < 0)
    stop("'lambda1' and 'lambda2' must be finite and non-negative", call. = FALSE)
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 4L) stop("'order' must be at least 4", call. = FALSE)
  r <- seq_len(order)
  kappa <- lambda1 + (-1)^r * lambda2
  from_cumulants <- function(k) {
    m <- numeric(order + 1L); m[1] <- 1                   # m[n + 1] = moment of order n
    for (n in r) {
      j <- 0:(n - 1)
      m[n + 1] <- sum(choose(n - 1, j) * k[n - j] * m[j + 1])
    }
    m[-1]
  }
  raw <- from_cumulants(kappa)
  central <- from_cumulants(replace(kappa, 1, 0))
  v <- lambda1 + lambda2
  names(kappa) <- names(raw) <- names(central) <- r
  structure(list(lambda1 = lambda1, lambda2 = lambda2,
                 mean = lambda1 - lambda2, variance = v, sd = sqrt(v),
                 skewness = if (v > 0) (lambda1 - lambda2) / v^1.5 else NaN,
                 kurtosis = if (v > 0) 1 / v else NaN,
                 cumulants = kappa, raw = raw, central = central),
            class = "skellam_moments")
}

#' @rdname skellam.moments
#' @export
print.skellam_moments <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(sprintf("Skellam distribution, lambda1 = %s, lambda2 = %s\n\n",
              format(x$lambda1, digits = digits), format(x$lambda2, digits = digits)))
  s <- c(mean = x$mean, variance = x$variance, sd = x$sd, skewness = x$skewness,
         `excess kurtosis` = x$kurtosis)
  print(s, digits = digits)
  cat("\nCumulants, raw and central moments by order:\n")
  print(rbind(cumulant = x$cumulants, raw = x$raw, central = x$central), digits = digits)
  invisible(x)
}

.skellam_gf_check <- function(lambda1, lambda2) {
  if (any(!is.finite(lambda1) | lambda1 < 0) || any(!is.finite(lambda2) | lambda2 < 0))
    stop("'lambda1' and 'lambda2' must be finite and non-negative", call. = FALSE)
}

#' @rdname skellam.moments
#' @export
skellam.mgf <- function(t, lambda1, lambda2 = lambda1) {
  .skellam_gf_check(lambda1, lambda2)
  exp(lambda1 * expm1(t) + lambda2 * expm1(-t))
}

#' @rdname skellam.moments
#' @export
skellam.cf <- function(t, lambda1, lambda2 = lambda1) {
  .skellam_gf_check(lambda1, lambda2)
  exp(lambda1 * (exp(1i * t) - 1) + lambda2 * (exp(-1i * t) - 1))
}

#' @rdname skellam.moments
#' @export
skellam.pgf <- function(z, lambda1, lambda2 = lambda1) {
  .skellam_gf_check(lambda1, lambda2)
  if (any(z == 0)) stop("'z' must be non-zero", call. = FALSE)
  exp(lambda1 * (z - 1) + lambda2 * (1 / z - 1))
}
