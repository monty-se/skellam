% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dskellam.r, R/dskellam.sp.r, R/pskellam.r,
%   R/pskellam.sp.r, R/qskellam.r, R/rskellam.r
\name{skellam}
\alias{skellam}
\alias{dskellam}
\alias{pskellam}
\alias{qskellam}
\alias{rskellam}
\alias{dskellam.sp}
\alias{pskellam.sp}
\title{The Skellam Distribution}
\usage{
dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)

dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)

pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

pskellam.sp(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)

rskellam(n, lambda1, lambda2 = lambda1)
}
\arguments{
\item{x, q}{For functions \code{dskellam}, \code{dskellam.sp}, and \code{pskellam.sp}: a numeric vector of quantiles.}

\item{lambda1, lambda2}{Numeric vectors of (non-negative) means; \code{lambda2} defaults to \code{lambda1} if not provided.}

\item{log, log.p}{Logical; if TRUE, returns the logarithm of the computed value.}

\item{lower.tail}{Logical; if TRUE (default), returns \eqn{P(X \le x)}; otherwise, returns \eqn{P(X > x)}.}

\item{p}{For \code{qskellam}: a numeric vector of probabilities.}

\item{n}{For \code{rskellam}: a non-negative integer specifying the number of observations.}
}
\value{
\itemize{
\item \code{dskellam} returns the (log) density.
\item \code{pskellam} returns the (log) cumulative distribution function.
\item \code{qskellam} returns the quantile function.
\item \code{rskellam} generates random deviates.
}

Invalid \code{lambda} values will return \code{NaN} with a warning.
}
\description{
Density, distribution function, quantile function, and random generation for the Skellam distribution.
}
\details{
The Skellam distribution describes the difference between two independent Poisson random variables.
This documentation covers:

\strong{Density:} \cr
\preformatted{
dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)
}

\strong{Distribution Function:} \cr
\preformatted{
pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
}

\strong{Quantile Function:} \cr
\preformatted{
qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
}

\strong{Random Generation:} \cr
\preformatted{
rskellam(n, lambda1, lambda2 = lambda1)
}

\strong{Saddlepoint Approximations:} \cr
\preformatted{
dskellam.sp(x, lambda1, lambda2 = lambda1, log = FALSE)
pskellam.sp(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
}

If \eqn{Y_1} and \eqn{Y_2} are Poisson variables with means \eqn{\mu_1} and \eqn{\mu_2} and correlation \eqn{\rho},
then \eqn{X = Y_1 - Y_2} is Skellam with parameters:

\deqn{\lambda_1 = \mu_1 - \rho \sqrt{\mu_1 \mu_2}}

\deqn{\lambda_2 = \mu_2 - \rho \sqrt{\mu_1 \mu_2}}

The density is given by:

\deqn{I(2 \sqrt{\lambda_1 \lambda_2}, |x|) (\lambda_1/\lambda_2)^{x/2} \exp(-\lambda_1-\lambda_2)}

where \eqn{I(y,\nu)} is the modified Bessel function of the first kind.
}
\note{
The \pkg{VGAM} package also provides Skellam functions. This implementation offers a broader working range,
correct handling when one rate parameter is zero, enhanced argument checking, and improved accuracy for \eqn{x < 0}
(in R versions prior to 2.9). Use \code{skellam::dskellam} or \code{VGAM::dskellam} to specify which implementation to use.
}
\examples{
# Compare with Poisson when one lambda = 0
dskellam(0:10, 5, 0)
dpois(0:10, 5)

# Both lambdas non-zero
dskellam(c(-1,1), c(12,10), c(10,12))
pskellam(c(-1,0), c(12,10), c(10,12))

# Quantile function
qskellam(c(0.05, 0.95), 3, 4)

# Random generation
rskellam(10, 8.5, 10.25)

}
\references{
\itemize{
\item Butler, R. (2007) \emph{Saddlepoint Approximations with Applications},
Cambridge University Press.
\item Johnson, N. L. (1959) On an extension of the connection between Poisson and \eqn{\chi^2} distributions.
\emph{Biometrika} 46, 352-362.
\item Johnson, N. L., Kotz, S., & Kemp, A. W. (1993) \emph{Univariate Discrete Distributions},
2nd ed., John Wiley and Sons.
\item Skellam, J. G. (1946) The frequency distribution of the difference between two Poisson variates.
\emph{Journal of the Royal Statistical Society, Series A} 109(3), 296.
\item Strackee, J. & van der Gon, J. J. D. (1962) The frequency distribution of the difference
between two Poisson variates. \emph{Statistica Neerlandica} 16(1), 17-23.
\item Wikipedia: \url{https://en.wikipedia.org/wiki/Skellam_distribution}
}
}
