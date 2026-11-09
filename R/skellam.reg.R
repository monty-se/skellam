## --------------------------------------------------------------------------
## skellam.reg.R -- Skellam regression
##
## Public functions:
##   skellam.reg()   log-linear regression on both rates; formula interface
##                   (y ~ x1 + x2 | z1, data) and the original (y, x) form;
##                   returns a "skellam_reg" object
##   S3 methods for "skellam_reg": print, summary, coef, vcov, logLik, nobs,
##                   fitted, residuals, predict, simulate, anova
##
## Internal functions:
##   .skellam_reg_fit()   the maximum likelihood engine (BFGS, analytic score)
##   .skellam_split_rhs() split a two-part formula at "|"
##
## Original author: Michail Tsagris
## Corrections:     Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' Skellam Regression
#'
#' Fits a regression model in which a difference of counts \eqn{y = Y_1 - Y_2}{y = Y1 - Y2}
#' follows a Skellam distribution and both rates depend on covariates through
#' log links.
#'
#' @param formula A formula \code{y ~ x1 + x2} (the same covariates for both
#'   rates) or a two-part formula \code{y ~ x1 + x2 | z1 + z2}, where the part
#'   before \code{|} is the model for \eqn{\lambda_1}{lambda1} and the part after
#'   it the model for \eqn{\lambda_2}{lambda2}. Offsets can be given with
#'   \code{offset()} in either part. For the original interface, a numeric
#'   vector \code{y}; see Details.
#' @param data A data frame containing the variables in the formula. For the
#'   original interface, the covariates \code{x}.
#' @param ... For the original interface, \code{y} and \code{x} may also be
#'   given by name.
#'
#' @details
#' The response is modelled as \eqn{y_i = Y_{1i} - Y_{2i}}{y = Y1 - Y2} with
#' independent \eqn{Y_{1i} \sim \mathrm{Poisson}(\lambda_{1i})}{Y1 ~ Poisson(lambda1)} and
#' \eqn{Y_{2i} \sim \mathrm{Poisson}(\lambda_{2i})}{Y2 ~ Poisson(lambda2)}, and
#' \eqn{\log\lambda_{1i} = x_{1i}^\top\beta_1}{log(lambda1) = x1'beta1},
#' \eqn{\log\lambda_{2i} = x_{2i}^\top\beta_2}{log(lambda2) = x2'beta2}. So the
#' \code{lambda1} coefficients describe the rate of the counts that are added
#' and the \code{lambda2} coefficients the rate of the counts that are
#' subtracted, as in \code{dskellam(y, lambda1, lambda2)}.
#'
#' The log-likelihood is computed with \code{dskellam(log = TRUE)} and
#' maximised by \code{\link[stats]{optim}} (BFGS) with the analytic gradient,
#' from deterministic starting values (least-squares fits of the log moment
#' estimates of an intercept-only model), so repeated calls give identical
#' results. Standard errors come from the inverse of the observed information
#' (\code{\link[stats]{optimHess}}). Observations with missing values in any
#' variable used are removed with a warning. A warning is also given if the
#' optimizer does not converge or the information matrix is not invertible
#' (standard errors are then \code{NA}).
#'
#' \strong{Original interface.} \code{skellam.reg(y, x)} with a numeric response
#' \code{y} and a vector, matrix or data frame of covariates \code{x} fits the
#' model with the same covariates (and an intercept) for both rates, exactly as
#' in earlier versions of the package.
#'
#' @return An object of class \code{"skellam_reg"}, a list with components
#' \describe{
#'   \item{coefficients}{A list with the coefficient vectors \code{lambda1} and
#'     \code{lambda2}.}
#'   \item{vcov}{Covariance matrix of all coefficients.}
#'   \item{loglik}{Maximised log-likelihood.}
#'   \item{param1, param2}{Coefficient tables (estimate, standard error, Wald
#'     statistic, p-value) for \eqn{\lambda_1}{lambda1} and
#'     \eqn{\lambda_2}{lambda2}, as returned by earlier versions.}
#'   \item{fitted.values}{A data frame with the fitted rates \code{lambda1} and
#'     \code{lambda2}.}
#'   \item{y, nobs, df.residual, converged, call, formula, terms, levels}{The
#'     response used, the number of observations, the residual degrees of
#'     freedom, whether the optimizer converged, and what is needed to predict
#'     on new data.}
#' }
#' Methods are available for \code{print}, \code{summary}, \code{coef},
#' \code{vcov}, \code{logLik} (hence \code{AIC}, \code{BIC}), \code{nobs},
#' \code{confint} (the default method), \code{fitted}, \code{residuals},
#' \code{predict}, \code{simulate}, \code{anova} and \code{update}.
#'
#' @seealso \code{\link{predict.skellam_reg}}, \code{\link{skellam.fit}} for a
#'   single sample, \code{\link{skellam}} for the distribution.
#'
#' @references
#' \itemize{
#'   \item Skellam, J. G. (1946) The frequency distribution of the difference between
#'   two Poisson variates belonging to different populations. \emph{Journal of the
#'   Royal Statistical Society, Series A} \bold{109}(3), 296.
#'
#'   \item Karlis, D. and Ntzoufras, I. (2009) Bayesian modelling of football outcomes:
#'   using the Skellam's distribution for the goal difference. \emph{IMA Journal of
#'   Management Mathematics} \bold{20}(2), 133-145.
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(attack = rnorm(n), defence = rnorm(n))
#' d$goaldiff <- rpois(n, exp(0.4 + 0.3 * d$attack)) -
#'               rpois(n, exp(0.1 - 0.2 * d$defence))
#'
#' # different covariates for the two rates
#' m <- skellam.reg(goaldiff ~ attack | defence, data = d)
#' m
#' summary(m)
#'
#' # nested models: likelihood-ratio test
#' m0 <- skellam.reg(goaldiff ~ 1 | 1, data = d)
#' anova(m0, m)
#'
#' # probabilities of each goal difference for new data
#' predict(m, newdata = data.frame(attack = 1, defence = 0), type = "prob", at = -2:2)
#'
#' # original interface: same covariates (and an intercept) for both rates
#' skellam.reg(d$goaldiff, d$attack)$param1
#'
#' @author Michail Tsagris; formula interface and methods by Montasser Ghachem
#'   and Oguz Ersan
#' @keywords models regression
#' @importFrom stats fitted residuals anova update
#' @export
skellam.reg <- function(formula, data = NULL, ...) {
  cl <- match.call()
  dots <- list(...)
  if (missing(formula) && !is.null(dots$y)) {
    return(.skellam_reg_xy(dots$y, if (!is.null(dots$x)) dots$x else data, cl))
  }
  if (!inherits(formula, "formula")) {
    xx <- if (!is.null(dots$x)) dots$x else data
    if (is.null(xx)) stop("covariates 'x' are missing", call. = FALSE)
    return(.skellam_reg_xy(formula, xx, cl))
  }

  if (is.null(data)) data <- environment(formula)
  parts <- .skellam_split_rhs(formula)
  f1 <- parts$f1; f2 <- parts$f2
  mf1 <- stats::model.frame(f1, data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  mf2 <- stats::model.frame(f2, data, na.action = stats::na.pass, drop.unused.levels = TRUE)
  t1 <- stats::terms(mf1); t2 <- stats::terms(mf2)
  y  <- stats::model.response(mf1)
  X1 <- stats::model.matrix(t1, mf1)
  X2 <- stats::model.matrix(t2, mf2)
  o1 <- stats::model.offset(mf1); if (is.null(o1)) o1 <- rep(0, length(y))
  o2 <- stats::model.offset(mf2); if (is.null(o2)) o2 <- rep(0, length(y))

  ok <- !is.na(y) & stats::complete.cases(X1) & stats::complete.cases(X2) & !is.na(o1) & !is.na(o2)
  if (!all(ok)) {
    warning(sprintf("%d observation(s) with missing values removed", sum(!ok)), call. = FALSE)
    y <- y[ok]; X1 <- X1[ok, , drop = FALSE]; X2 <- X2[ok, , drop = FALSE]; o1 <- o1[ok]; o2 <- o2[ok]
  }
  fit <- .skellam_reg_fit(y, X1, X2, o1, o2)
  fit$call <- cl
  fit$formula <- formula
  fit$terms <- list(lambda1 = stats::delete.response(t1), lambda2 = stats::delete.response(t2))
  fit$levels <- list(lambda1 = stats::.getXlevels(t1, mf1), lambda2 = stats::.getXlevels(t2, mf2))
  fit$contrasts <- list(lambda1 = attr(X1, "contrasts"), lambda2 = attr(X2, "contrasts"))
  fit
}

# Original interface: same covariates and an intercept for both rates.
.skellam_reg_xy <- function(y, x, cl) {
  xdf <- data.frame(x)
  X <- stats::model.matrix(~., stats::model.frame(~., xdf, na.action = stats::na.pass))
  if (length(y) != nrow(X)) stop("'y' and 'x' must have the same number of observations", call. = FALSE)
  if (!is.numeric(y)) stop("'y' must be a numeric vector of integers", call. = FALSE)
  ok <- !is.na(y) & stats::complete.cases(X)
  if (!all(ok)) {
    warning(sprintf("%d observation(s) with missing values removed", sum(!ok)), call. = FALSE)
    y <- y[ok]; X <- X[ok, , drop = FALSE]
  }
  fit <- .skellam_reg_fit(y, X, X, rep(0, length(y)), rep(0, length(y)))
  fit$call <- cl
  fit$formula <- NULL
  fit
}

# Split  y ~ a + b | c  into  y ~ a + b  and  y ~ c ; without "|" both parts
# are the whole right-hand side.
.skellam_split_rhs <- function(formula) {
  if (length(formula) != 3L) stop("the formula needs a response: y ~ ...", call. = FALSE)
  lhs <- formula[[2]]; rhs <- formula[[3]]
  # update() can wrap the right-hand side in parentheses
  while (is.call(rhs) && identical(rhs[[1]], as.name("("))) rhs <- rhs[[2]]
  is_bar <- function(e) is.call(e) && identical(e[[1]], as.name("|"))
  if (is_bar(rhs)) {
    r1 <- rhs[[2]]; r2 <- rhs[[3]]
    if (is_bar(r1) || is_bar(r2))
      stop("only one '|' is allowed: y ~ (model for lambda1) | (model for lambda2)", call. = FALSE)
  } else {
    r1 <- r2 <- rhs
  }
  env <- environment(formula)
  f1 <- stats::as.formula(call("~", lhs, r1), env = env)
  f2 <- stats::as.formula(call("~", lhs, r2), env = env)
  list(f1 = f1, f2 = f2)
}

# Maximum likelihood for y with log(lambda_k) = X_k beta_k + o_k.
.skellam_reg_fit <- function(y, X1, X2, o1, o2) {
  if (any(!is.finite(y)) || any(.skellam_nonint(y)))
    stop("'y' must contain finite integers (differences of counts)", call. = FALSE)
  y <- round(as.numeric(y))
  n <- length(y); p1 <- ncol(X1); p2 <- ncol(X2); p <- p1 + p2
  if (n <= p) stop("too few observations for the number of parameters", call. = FALSE)
  i1 <- seq_len(p1); i2 <- p1 + seq_len(p2)

  # d log f / d log lambda1 = (z/2) I'_nu(z)/I_nu(z) - lambda1 + y/2,
  # d log f / d log lambda2 = (z/2) I'_nu(z)/I_nu(z) - lambda2 - y/2,
  # with nu = |y|, z = 2 sqrt(lambda1 lambda2), I'_nu/I_nu = I_{nu+1}/I_nu + nu/z.
  rates <- function(pa) list(l1 = exp(drop(X1 %*% pa[i1]) + o1), l2 = exp(drop(X2 %*% pa[i2]) + o2))
  # the line search can try steps where a rate overflows or underflows;
  # report those points as "infinitely bad" instead of evaluating there
  usable <- function(r) all(is.finite(r$l1) & is.finite(r$l2) & r$l1 > 0 & r$l2 > 0)
  nll <- function(pa) {
    r <- rates(pa)
    if (!usable(r)) return(1e300)
    v <- -sum(dskellam(y, r$l1, r$l2, log = TRUE))
    if (is.finite(v)) v else 1e300
  }
  grad <- function(pa) {
    r <- rates(pa)
    if (!usable(r)) return(rep(0, p))
    nu <- abs(y)
    z  <- 2 * sqrt(r$l1) * sqrt(r$l2)
    ratio <- exp(.log_besselI_scaled(z, nu + 1) - .log_besselI_scaled(z, nu))
    common <- 0.5 * z * ratio + 0.5 * nu
    -c(crossprod(X1, common - r$l1 + 0.5 * y), crossprod(X2, common - r$l2 - 0.5 * y))
  }

  # deterministic start: moment estimates of an intercept-only model, turned
  # into coefficients by least squares (an intercept gets the log rate,
  # slopes start near 0)
  m  <- mean(y)
  l2 <- max((stats::var(y) - m) / 2, 0.5)
  l1 <- max(l2 + m, 0.5)
  ls_start <- function(X, o, l) {
    b <- tryCatch(stats::lm.fit(X, rep(log(l), n) - o)$coefficients, error = function(e) rep(0, ncol(X)))
    b[is.na(b)] <- 0
    b
  }
  start <- c(ls_start(X1, o1, l1), ls_start(X2, o2, l2))

  mod <- stats::optim(start, nll, grad, method = "BFGS",
                      control = list(maxit = 1000, reltol = 1e-12))
  converged <- mod$convergence == 0
  if (!converged)
    warning("optimizer did not converge (code ", mod$convergence, "); estimates may be unreliable",
            call. = FALSE)

  b <- mod$par
  nm1 <- colnames(X1); nm2 <- colnames(X2)
  H <- stats::optimHess(b, nll, grad)
  V <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(diag(V))) || any(diag(V) <= 0)) {
    warning("the information matrix is not positive definite; standard errors are not available",
            call. = FALSE)
    V <- matrix(NA_real_, p, p)
  } else {
    V <- (V + t(V)) / 2
  }
  vnames <- c(paste0("lambda1:", nm1), paste0("lambda2:", nm2))
  dimnames(V) <- list(vnames, vnames)
  s <- sqrt(diag(V))
  tab <- function(est, se, nm) {
    out <- cbind(est, se, est / se, stats::pchisq((est / se)^2, 1, lower.tail = FALSE))
    dimnames(out) <- list(nm, c("Estimate", "Std. Error", "Wald value", "p-value"))
    out
  }
  r <- rates(b)
  structure(list(
    loglik = -mod$value,
    param1 = tab(b[i1], s[i1], nm1),
    param2 = tab(b[i2], s[i2], nm2),
    coefficients = list(lambda1 = stats::setNames(b[i1], nm1), lambda2 = stats::setNames(b[i2], nm2)),
    vcov = V,
    fitted.values = data.frame(lambda1 = r$l1, lambda2 = r$l2),
    y = y, nobs = n, df.residual = n - p, converged = converged,
    offset = list(lambda1 = o1, lambda2 = o2)),
    class = "skellam_reg")
}

#' Methods for Skellam Regression Fits
#'
#' Extractor, prediction and testing methods for objects returned by
#' \code{\link{skellam.reg}}.
#'
#' @param object,x A \code{"skellam_reg"} object.
#' @param model For \code{coef}: \code{"full"} (all coefficients, named
#'   \code{lambda1:...} and \code{lambda2:...}), \code{"lambda1"} or
#'   \code{"lambda2"}.
#' @param newdata Optional data frame for \code{predict}; if omitted, the
#'   fitted values are used.
#' @param type For \code{predict}: \code{"response"} (the fitted mean
#'   \eqn{\lambda_1 - \lambda_2}{lambda1 - lambda2}), \code{"lambda1"},
#'   \code{"lambda2"}, \code{"rates"} (both), \code{"prob"} (the probabilities
#'   of the values in \code{at}) or \code{"cdf"}. For \code{residuals}:
#'   \code{"response"} (\eqn{y - (\lambda_1 - \lambda_2)}{y - (lambda1 - lambda2)}),
#'   \code{"pearson"} (divided by \eqn{\sqrt{\lambda_1 + \lambda_2}}{sqrt(lambda1 + lambda2)})
#'   or \code{"quantile"} (randomized quantile residuals, Dunn and Smyth 1996,
#'   standard normal when the model is right; they use random numbers).
#' @param at Integer values at which \code{type = "prob"} or \code{"cdf"} are
#'   evaluated; one column per value.
#' @param nsim,seed As in \code{\link[stats]{simulate}}.
#' @param digits Number of significant digits to print.
#' @param ... Further fitted models (for \code{anova}); otherwise unused.
#'
#' @return \code{coef}, \code{vcov}, \code{logLik}, \code{nobs}, \code{fitted}
#'   and \code{residuals} return what their names say; \code{predict} returns a
#'   vector, or a matrix for \code{type = "rates"}, \code{"prob"} and
#'   \code{"cdf"}; \code{simulate} a data frame of simulated responses;
#'   \code{anova} a likelihood-ratio test table for nested models.
#'
#' @references Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile
#'   residuals. \emph{Journal of Computational and Graphical Statistics}
#'   \bold{5}(3), 236-244.
#'
#' @name skellam_reg-methods
#' @aliases predict.skellam_reg
NULL

#' @rdname skellam_reg-methods
#' @export
print.skellam_reg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Skellam regression")
  if (!is.null(x$call)) cat("\n\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nCoefficients for log(lambda1), the rate of the counts that are added:\n")
  print(x$coefficients$lambda1, digits = digits)
  cat("\nCoefficients for log(lambda2), the rate of the counts that are subtracted:\n")
  print(x$coefficients$lambda2, digits = digits)
  cat(sprintf("\nn = %d, log-likelihood %s, AIC %s\n", x$nobs, format(x$loglik, digits = digits + 3),
              format(stats::AIC(x), digits = digits + 3)))
  invisible(x)
}

#' @rdname skellam_reg-methods
#' @export
summary.skellam_reg <- function(object, ...) {
  mk <- function(tab) {
    colnames(tab) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    tab
  }
  structure(list(call = object$call, lambda1 = mk(object$param1), lambda2 = mk(object$param2),
                 loglik = object$loglik, aic = stats::AIC(object), bic = stats::BIC(object),
                 nobs = object$nobs, df.residual = object$df.residual, converged = object$converged),
            class = "summary.skellam_reg")
}

#' @export
print.summary.skellam_reg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!is.null(x$call)) cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n", sep = "")
  cat("\nlog(lambda1): rate of the counts that are added\n")
  stats::printCoefmat(x$lambda1, digits = digits, has.Pvalue = TRUE)
  cat("\nlog(lambda2): rate of the counts that are subtracted\n")
  stats::printCoefmat(x$lambda2, digits = digits, has.Pvalue = TRUE)
  cat(sprintf("\nn = %d, residual df = %d, log-likelihood %s, AIC %s, BIC %s\n", x$nobs, x$df.residual,
              format(x$loglik, digits = digits + 3), format(x$aic, digits = digits + 3),
              format(x$bic, digits = digits + 3)))
  if (!isTRUE(x$converged)) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}

#' @rdname skellam_reg-methods
#' @export
coef.skellam_reg <- function(object, model = c("full", "lambda1", "lambda2"), ...) {
  model <- match.arg(model)
  if (model == "full") {
    c(stats::setNames(object$coefficients$lambda1, paste0("lambda1:", names(object$coefficients$lambda1))),
      stats::setNames(object$coefficients$lambda2, paste0("lambda2:", names(object$coefficients$lambda2))))
  } else object$coefficients[[model]]
}

#' @rdname skellam_reg-methods
#' @export
vcov.skellam_reg <- function(object, ...) object$vcov

#' @rdname skellam_reg-methods
#' @export
logLik.skellam_reg <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients$lambda1) + length(object$coefficients$lambda2),
            nobs = object$nobs, class = "logLik")
}

#' @rdname skellam_reg-methods
#' @export
nobs.skellam_reg <- function(object, ...) object$nobs

#' @rdname skellam_reg-methods
#' @export
fitted.skellam_reg <- function(object, ...) object$fitted.values$lambda1 - object$fitted.values$lambda2

#' @rdname skellam_reg-methods
#' @export
residuals.skellam_reg <- function(object, type = c("response", "pearson", "quantile"), ...) {
  type <- match.arg(type)
  l1 <- object$fitted.values$lambda1; l2 <- object$fitted.values$lambda2; y <- object$y
  switch(type,
    response = y - (l1 - l2),
    pearson  = (y - (l1 - l2)) / sqrt(l1 + l2),
    quantile = {
      lo <- pskellam(y - 1, l1, l2); hi <- pskellam(y, l1, l2)
      stats::qnorm(lo + stats::runif(length(y)) * (hi - lo))
    })
}

#' @rdname skellam_reg-methods
#' @export
predict.skellam_reg <- function(object, newdata = NULL,
                                type = c("response", "lambda1", "lambda2", "rates", "prob", "cdf"),
                                at = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    l1 <- object$fitted.values$lambda1; l2 <- object$fitted.values$lambda2
  } else {
    if (is.null(object$terms)) stop("predicting on new data needs a model fitted with the formula interface",
                                     call. = FALSE)
    mm <- function(k) {
      mf <- stats::model.frame(object$terms[[k]], newdata, na.action = stats::na.pass,
                               xlev = object$levels[[k]])
      X  <- stats::model.matrix(object$terms[[k]], mf, contrasts.arg = object$contrasts[[k]])
      off <- stats::model.offset(mf); if (is.null(off)) off <- 0
      unname(exp(drop(X %*% object$coefficients[[k]]) + off))
    }
    l1 <- mm("lambda1"); l2 <- mm("lambda2")
  }
  if (type %in% c("prob", "cdf")) {
    if (is.null(at)) stop("give the values of the response in 'at'", call. = FALSE)
    f <- if (type == "prob") dskellam else pskellam
    out <- vapply(at, function(a) f(rep(a, length(l1)), l1, l2), numeric(length(l1)))
    out <- matrix(out, nrow = length(l1))
    colnames(out) <- as.character(at)
    return(out)
  }
  switch(type,
         response = l1 - l2,
         lambda1  = l1,
         lambda2  = l2,
         rates    = cbind(lambda1 = l1, lambda2 = l2))
}

#' @rdname skellam_reg-methods
#' @export
simulate.skellam_reg <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  l1 <- object$fitted.values$lambda1; l2 <- object$fitted.values$lambda2
  n <- length(l1)
  out <- as.data.frame(matrix(rskellam(n * nsim, rep(l1, nsim), rep(l2, nsim)), n, nsim))
  names(out) <- paste0("sim_", seq_len(nsim))
  out
}

#' @rdname skellam_reg-methods
#' @export
anova.skellam_reg <- function(object, ...) {
  models <- c(list(object), list(...))
  if (length(models) < 2L) stop("give at least two nested models to compare", call. = FALSE)
  if (!all(vapply(models, inherits, NA, "skellam_reg"))) stop("all models must be skellam.reg fits", call. = FALSE)
  n <- vapply(models, function(m) m$nobs, 0L)
  if (length(unique(n)) > 1L) stop("models were fitted to different numbers of observations", call. = FALSE)
  ll <- vapply(models, function(m) as.numeric(logLik(m)), 0)
  df <- vapply(models, function(m) attr(logLik(m), "df"), 0)
  o <- order(df); ll <- ll[o]; df <- df[o]
  chisq <- c(NA, 2 * diff(ll)); ddf <- c(NA, diff(df))
  p <- stats::pchisq(chisq, ddf, lower.tail = FALSE)
  tab <- data.frame(`#Df` = df, LogLik = ll, Df = ddf, Chisq = chisq, `Pr(>Chisq)` = p,
                    check.names = FALSE)
  rownames(tab) <- paste("Model", seq_along(models))
  forms <- vapply(models[o], function(m) if (is.null(m$formula)) deparse(m$call) else
    paste(deparse(m$formula), collapse = ""), "")
  structure(tab, heading = c("Likelihood ratio test\n", paste0("Model ", seq_along(forms), ": ", forms, collapse = "\n"), ""),
            class = c("anova", "data.frame"))
}
