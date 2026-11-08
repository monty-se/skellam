## --------------------------------------------------------------------------
## dskellam.r -- Skellam probability mass function
##
## Public functions:
##   dskellam()      P(X = x), on the natural or log scale
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' The Skellam Distribution
#'
#' Density, distribution function, quantile function, and random generation for the Skellam distribution.
#'
#' The Skellam distribution describes the difference between two independent Poisson random variables.
#' This documentation covers:
#'
#' \strong{Density:} \cr
#' \preformatted{
#' dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)
#' }
#'
#' \strong{Distribution Function:} \cr
#' \preformatted{
#' pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#' \strong{Quantile Function:} \cr
#' \preformatted{
#' qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#' \strong{Random Generation:} \cr
#' \preformatted{
#' rskellam(n, lambda1, lambda2 = lambda1)
#' }
#'
#' \strong{Saddlepoint Approximations:} \cr
#' \preformatted{
#' dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)
#' pskellam.sp(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
#' }
#'
#'
#' @details
#' If \eqn{Y_1} and \eqn{Y_2} are Poisson variables with means \eqn{\mu_1} and \eqn{\mu_2} and correlation \eqn{\rho},
#' then \eqn{X = Y_1 - Y_2} is Skellam with parameters:
#'
#' \deqn{\lambda_1 = \mu_1 - \rho \sqrt{\mu_1 \mu_2}}
#'
#' \deqn{\lambda_2 = \mu_2 - \rho \sqrt{\mu_1 \mu_2}}
#'
#' The density is given by:
#'
#' \deqn{I(2 \sqrt{\lambda_1 \lambda_2}, |x|) (\lambda_1/\lambda_2)^{x/2} \exp(-\lambda_1-\lambda_2)}
#'
#' where \eqn{I(y,\nu)} is the modified Bessel function of the first kind.
#'
#' @param x,q For functions \code{dskellam}, \code{dskellam.sp}, and \code{pskellam.sp}: a numeric vector of quantiles.
#' @param p For \code{qskellam}: a numeric vector of probabilities.
#' @param n For \code{rskellam}: a non-negative integer specifying the number of observations.
#' @param lambda1,lambda2 Numeric vectors of (non-negative) means; \code{lambda2} defaults to \code{lambda1} if not provided.
#' @param log,log.p Logical; if TRUE, returns the logarithm of the computed value.
#' @param lower.tail Logical; if TRUE (default), returns \eqn{P(X \le x)}; otherwise, returns \eqn{P(X > x)}.
#'
#' @return
#' \itemize{
#'   \item \code{dskellam} returns the (log) density.
#'   \item \code{pskellam} returns the (log) cumulative distribution function.
#'   \item \code{qskellam} returns the quantile function.
#'   \item \code{rskellam} generates random deviates.
#' }
#'
#' Invalid \code{lambda} values will return \code{NaN} with a warning.
#'
#' @note
#' The \pkg{VGAM} package also provides Skellam functions. This implementation offers a broader working range,
#' correct handling when one rate parameter is zero, enhanced argument checking, and improved accuracy for \eqn{x < 0}
#' (in R versions prior to 2.9). Use \code{skellam::dskellam} or \code{VGAM::dskellam} to specify which implementation to use.
#'
#' @references
#' \itemize{
#'   \item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
#'         Cambridge University Press.
#'   \item Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{\chi^2} distributions.
#'         \emph{Biometrika} 46, 352-362.
#'   \item Johnson, N. L., Kotz, S., & Kemp, A. W. (1993) \emph{Univariate Discrete Distributions},
#'         2nd ed., John Wiley and Sons.
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates.
#'         \emph{Journal of the Royal Statistical Society, Series A} 109(3), 296.
#'   \item Strackee, J. & van der Gon, J. J. D. (1962) The frequency distribution of the difference
#'         between two Poisson variates. \emph{Statistica Neerlandica} 16(1), 17-23.
#'   \item Wikipedia: \url{https://en.wikipedia.org/wiki/Skellam_distribution}
#' }
#'
#' @examples
#' # Compare with Poisson when one lambda = 0
#' dskellam(0:10, 5, 0)
#' dpois(0:10, 5)
#'
#' # Both lambdas non-zero
#' dskellam(c(-1,1), c(12,10), c(10,12))
#' pskellam(c(-1,0), c(12,10), c(10,12))
#'
#' # Quantile function
#' qskellam(c(0.05, 0.95), 3, 4)
#'
#' # Random generation
#' rskellam(10, 8.5, 10.25)
#'
#' @name skellam
#' @aliases dskellam pskellam qskellam rskellam dskellam.sp pskellam.sp
NULL

NULL


#' @export
#' @rdname skellam
dskellam <- function(x, lambda1, lambda2=lambda1, log=FALSE){
 # density (PMF) of Skellam distriubition (difference of Poissons)
    log <- isTRUE(as.logical(log[1]))
    nms <- names(x)
    nx <- length(x)
   # recycle to a common length; warn when lengths are not multiples (N13)
    a <- .skellam_recycle(x, lambda1, lambda2)
    x <- as.numeric(a[[1]]); lambda1 <- as.numeric(a[[2]]); lambda2 <- as.numeric(a[[3]])
    len <- length(x)
    ret <- rep(NaN, length.out = len)
   # invalid rates -> NaN with a warning, as dpois does
    bad <- .skellam_bad_rate(lambda1) | .skellam_bad_rate(lambda2)
    if (any(bad)) warning("NaNs produced")
   # missing values propagate (N1): NA stays NA, NaN stays NaN
    miss <- !bad & (is.na(x) | is.na(lambda1) | is.na(lambda2))
    ret[miss] <- (x + lambda1 + lambda2)[miss]
    ok <- !bad & !miss
   # outside the integer support the pmf is 0 (N11); infinite x likewise
    nonint <- ok & .skellam_nonint(x)
    if (any(nonint)) {
        xn <- x[nonint]
        warning("non-integer x = ", paste(format(xn[seq_len(min(5L, length(xn)))]), collapse = ", "),
                if (length(xn) > 5) ", ..." else "", call. = FALSE)
    }
    off <- nonint | (ok & is.infinite(x))
    ret[off] <- if (log) -Inf else 0
    ok <- ok & !off
    x[ok] <- round(x[ok])
   # handle a zero lambda separately (avoids division by zero & accuracy issues for large values of lambda or x)
    z1 <- ok & lambda1 == 0
    z2 <- ok & lambda2 == 0
    ret[z1] <- stats::dpois(-x[z1], lambda2[z1], log = log)
    ret[z2] <- stats::dpois( x[z2], lambda1[z2], log = log)    # corrects bug in VGAM 0.7-9
   # both rates positive: log-space core (see R/bessel.R), finite for every input
    nz <- ok & (lambda1 > 0) & (lambda2 > 0)
    if (any(nz)) {
        lp <- .skellam_logpmf(x[nz], lambda1[nz], lambda2[nz])
        ret[nz] <- if (log) lp else exp(lp)
    }
    if (length(nms) && nx == len) names(ret) <- nms
    ret
}
