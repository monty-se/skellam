## --------------------------------------------------------------------------
## pskellam.r -- Skellam cumulative distribution function
##
## Public functions:
##   pskellam()      P(X <= q) or P(X > q), via the non-central chi-square identity
##
## Original author: Jerry W. Lewis
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @export
#' @rdname skellam
pskellam <- function(q, lambda1, lambda2=lambda1, lower.tail=TRUE, log.p=FALSE){
 # CDF of Skellam distribution (difference of Poissons)
    lower.tail <- isTRUE(as.logical(lower.tail[1]))
    log.p <- isTRUE(as.logical(log.p[1]))
    nms <- names(q)
    nq <- length(q)
    a <- .skellam_recycle(q, lambda1, lambda2)
    x <- as.numeric(a[[1]]); lambda1 <- as.numeric(a[[2]]); lambda2 <- as.numeric(a[[3]])
    len <- length(x)
    ret <- rep(NaN, length.out = len)
    bad <- .skellam_bad_rate(lambda1) | .skellam_bad_rate(lambda2)
    if (any(bad)) warning("NaNs produced")
    miss <- !bad & (is.na(x) | is.na(lambda1) | is.na(lambda2))
    ret[miss] <- (x + lambda1 + lambda2)[miss]
    ok <- !bad & !miss
   # CDF is a step function, so convert to integer values without warning
   # (the small tolerance keeps 2.9999999999996 at 3, as ppois does)
    fin <- ok & is.finite(x)
    x[fin] <- floor(x[fin] + 1e-7 * pmax(1, abs(x[fin])))
   # q = +-Inf: the cdf is 0 or 1
    p0 <- if (log.p) -Inf else 0
    p1 <- if (log.p) 0 else 1
    ninf <- ok & x == -Inf
    pinf <- ok & x == Inf
    ret[ninf] <- if (lower.tail) p0 else p1
    ret[pinf] <- if (lower.tail) p1 else p0
    ok <- ok & is.finite(x)
   # a zero rate reduces to a (reflected) Poisson; handled exactly (N10)
   #   lambda2 = 0: X = Y1,  P(X <= x) = P(Y1 <= x)
   #   lambda1 = 0: X = -Y2, P(X <= x) = P(Y2 >= -x) = P(Y2 > -x - 1)
    z2 <- ok & lambda2 == 0
    z1 <- ok & lambda1 == 0 & !z2
    ret[z2] <- stats::ppois(x[z2], lambda1[z2], lower.tail = lower.tail, log.p = log.p)
    ret[z1] <- stats::ppois(-x[z1] - 1, lambda2[z1], lower.tail = !lower.tail, log.p = log.p)
   # both rates positive: different formulas for negative & nonnegative x
    nz <- ok & !z1 & !z2
    neg <- nz & (x <  0)
    pos <- nz & (x >= 0)
    if (lower.tail) {
        ret[neg] <- stats::pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],log.p=log.p)
        ret[pos] <- stats::pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],lower.tail=FALSE,log.p=log.p)
    } else {
        ret[neg] <- stats::pchisq(2*lambda2[neg],-2*x[neg],2*lambda1[neg],lower.tail=FALSE,log.p=log.p)
        ret[pos] <- stats::pchisq(2*lambda1[pos],2*(x[pos]+1),2*lambda2[pos],log.p=log.p)
    }
    chk <- nz & (!is.finite(ret) | (!log.p)&(ret<1e-308))    # use saddlepoint approximation if outside the working range of pchisq
    if (any(chk)) ret[chk] <- pskellam.sp(x[chk],lambda1[chk],lambda2[chk],lower.tail,log.p)
    if (length(nms) && nq == len) names(ret) <- nms
    ret
}
