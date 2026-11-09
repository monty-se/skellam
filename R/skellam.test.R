## --------------------------------------------------------------------------
## skellam.test.R -- goodness-of-fit test for the Skellam distribution
##
## Public functions:
##   skellam.test()  Pearson chi-squared or likelihood-ratio statistic on
##                   binned counts, with a parametric-bootstrap or an
##                   asymptotic chi-squared p-value; returns an "htest"
##
## Internal functions:
##   .skellam_bins()     cells with expected counts >= min.expected
##   .skellam_gof_stat() the statistic for given data and rates
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

# Cells for the binned test. Start from one cell per integer between the
# 1e-6 and 1 - 1e-6 quantiles, with the two end cells absorbing the tails, then
# merge cells (always the smallest expected count with its smaller neighbour)
# until every expected count is at least min.expected. Returns the lower and
# upper limits of each cell (-Inf / Inf at the ends) and the probabilities.
.skellam_bins <- function(l1, l2, n, min.expected) {
  lo <- qskellam(1e-6, l1, l2); hi <- qskellam(1 - 1e-6, l1, l2)
  if (!is.finite(lo)) lo <- 0
  if (!is.finite(hi)) hi <- 0
  if (hi <= lo) return(data.frame(lower = -Inf, upper = Inf, prob = 1))
  v <- lo:hi
  p <- dskellam(v, l1, l2)
  p[1] <- pskellam(lo, l1, l2)
  p[length(p)] <- pskellam(hi - 1, l1, l2, lower.tail = FALSE)
  lower <- c(-Inf, v[-1]); upper <- c(v[-length(v)], Inf)
  repeat {
    e <- n * p
    k <- length(p)
    if (k == 1L || min(e) >= min.expected) break
    i <- which.min(e)
    j <- if (i == 1L) 2L else if (i == k) k - 1L else if (e[i - 1L] <= e[i + 1L]) i - 1L else i + 1L
    a <- min(i, j); b <- max(i, j)
    p[a] <- p[a] + p[b]; upper[a] <- upper[b]
    p <- p[-b]; lower <- lower[-b]; upper <- upper[-b]
  }
  data.frame(lower = lower, upper = upper, prob = p)
}

.skellam_gof_stat <- function(x, l1, l2, min.expected, type) {
  bins <- .skellam_bins(l1, l2, length(x), min.expected)
  cell <- findInterval(x, c(bins$lower[-1] - 0.5))         # 0-based cell index
  obs <- tabulate(cell + 1L, nbins = nrow(bins))
  exp <- length(x) * bins$prob
  stat <- if (type == "pearson") sum((obs - exp)^2 / exp) else {
    nz <- obs > 0
    2 * sum(obs[nz] * log(obs[nz] / exp[nz]))
  }
  list(stat = stat, obs = obs, exp = exp, bins = bins)
}

#' Goodness-of-Fit Test for the Skellam Distribution
#'
#' Tests whether a sample of differences of counts comes from a Skellam
#' distribution, with rates either given or estimated by maximum likelihood.
#'
#' @param x A numeric vector of integers.
#' @param lambda1,lambda2 Rates of the hypothesised distribution. If both are
#'   \code{NULL} (default) they are estimated from \code{x} with
#'   \code{\link{skellam.fit}}.
#' @param statistic \code{"pearson"} (default) for the Pearson chi-squared
#'   statistic or \code{"lr"} for the likelihood-ratio (G) statistic.
#' @param method \code{"bootstrap"} (default) for a parametric-bootstrap
#'   p-value, or \code{"asymptotic"} for the chi-squared approximation.
#' @param B Number of bootstrap samples.
#' @param min.expected Minimum expected count per cell; sparse tail cells are
#'   merged until every cell reaches it.
#'
#' @details
#' The integers are grouped into cells so that every cell has an expected count
#' of at least \code{min.expected} under the hypothesised (or fitted)
#' distribution: the cells start one per integer, the two end cells take the
#' tails, and the cell with the smallest expected count is merged with its
#' smaller neighbour until all reach the minimum. With \eqn{k} cells, the
#' statistic is \eqn{X^2 = \sum (O - E)^2 / E}{X^2 = sum (O - E)^2 / E} or
#' \eqn{G = 2 \sum O \log(O/E)}{G = 2 sum O log(O/E)}.
#'
#' \strong{p-values.} When the rates are given, both statistics are
#' approximately \eqn{\chi^2_{k-1}}{chi-squared with k - 1 df}. When the rates
#' are estimated from the ungrouped data, the limiting distribution lies between
#' \eqn{\chi^2_{k-3}}{chi-squared(k - 3)} and \eqn{\chi^2_{k-1}}{chi-squared(k - 1)}
#' (Chernoff and Lehmann, 1954), so no chi-squared reference is exact:
#' \code{method = "asymptotic"} uses \eqn{k - 3} degrees of freedom, which tends
#' to reject slightly too often. The default parametric bootstrap avoids the
#' problem: it simulates \code{B} samples from the fitted distribution,
#' re-estimates the rates and rebuilds the cells for each, and reports
#' \eqn{(1 + \#\{T^* \ge T\})/(B + 1)}{(1 + #(T* >= T))/(B + 1)}. Set a seed
#' for reproducible bootstrap p-values.
#'
#' @return An object of class \code{"htest"} with the usual components
#'   (\code{statistic}, \code{parameter}: degrees of freedom or \code{B},
#'   \code{p.value}, \code{method}, \code{data.name}, and \code{estimate}: the
#'   rates used), plus \code{observed}, \code{expected} and \code{cells} (a data
#'   frame with the limits of each cell).
#'
#' @references Chernoff, H. and Lehmann, E. L. (1954) The use of maximum
#'   likelihood estimates in \eqn{\chi^2} tests for goodness of fit.
#'   \emph{Annals of Mathematical Statistics} \bold{25}(3), 579-586.
#'
#' @seealso \code{\link{skellam.fit}}, \code{\link[stats]{chisq.test}}
#'
#' @examples
#' set.seed(1)
#' x <- rpois(500, 4) - rpois(500, 2)
#' skellam.test(x, B = 199)
#'
#' # a sample with extra zeros is rejected
#' z <- x; z[sample(500, 60)] <- 0
#' skellam.test(z, B = 199)
#'
#' # rates given in advance: asymptotic chi-squared with k - 1 df
#' skellam.test(x, lambda1 = 4, lambda2 = 2, method = "asymptotic")
#'
#' @export
skellam.test <- function(x, lambda1 = NULL, lambda2 = NULL,
                         statistic = c("pearson", "lr"), method = c("bootstrap", "asymptotic"),
                         B = 999, min.expected = 5) {
  dname <- paste(deparse(substitute(x)), collapse = " ")
  statistic <- match.arg(statistic); method <- match.arg(method)
  x <- .skellam_check_sample(x)
  n <- length(x)
  given <- !is.null(lambda1) || !is.null(lambda2)
  if (given) {
    if (is.null(lambda1) || is.null(lambda2)) stop("give both 'lambda1' and 'lambda2', or neither", call. = FALSE)
    if (length(lambda1) != 1L || length(lambda2) != 1L || !is.finite(lambda1) || !is.finite(lambda2) ||
        lambda1 < 0 || lambda2 < 0) stop("'lambda1' and 'lambda2' must be single non-negative numbers", call. = FALSE)
    lam <- c(lambda1 = lambda1, lambda2 = lambda2)
  } else {
    lam <- .skellam_profile_mle(x)$lambda
  }
  obs <- .skellam_gof_stat(x, lam[[1]], lam[[2]], min.expected, statistic)
  k <- length(obs$obs)
  sname <- if (statistic == "pearson") "X-squared" else "G"

  if (method == "asymptotic") {
    df <- k - if (given) 1L else 3L
    if (df < 1L) stop(sprintf("only %d cells; too few for an asymptotic test (use method = \"bootstrap\")", k),
                      call. = FALSE)
    p <- stats::pchisq(obs$stat, df, lower.tail = FALSE)
    param <- c(df = df)
    mtext <- sprintf("Skellam goodness-of-fit test (%s, asymptotic chi-squared%s)",
                     if (statistic == "pearson") "Pearson" else "likelihood ratio",
                     if (given) "" else ", df = cells - 3")
  } else {
    B <- as.integer(B)
    if (B < 1L) stop("'B' must be a positive integer", call. = FALSE)
    tstar <- numeric(B)
    for (b in seq_len(B)) {
      xs <- rskellam(n, lam[[1]], lam[[2]])
      lb <- if (given || stats::var(xs) == 0) lam else .skellam_profile_mle(xs)$lambda
      tstar[b] <- .skellam_gof_stat(xs, lb[[1]], lb[[2]], min.expected, statistic)$stat
    }
    p <- (1 + sum(tstar >= obs$stat)) / (B + 1)
    param <- c(B = B)
    mtext <- sprintf("Skellam goodness-of-fit test (%s, parametric bootstrap)",
                     if (statistic == "pearson") "Pearson" else "likelihood ratio")
  }
  structure(list(statistic = stats::setNames(obs$stat, sname), parameter = param, p.value = p,
                 method = mtext, data.name = dname, estimate = lam,
                 observed = obs$obs, expected = obs$exp,
                 cells = data.frame(lower = obs$bins$lower, upper = obs$bins$upper,
                                    observed = obs$obs, expected = obs$exp)),
            class = "htest")
}
