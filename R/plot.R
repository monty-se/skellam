## --------------------------------------------------------------------------
## plot.R -- diagnostic plots and the package's ggplot2 theme
##
## Public functions:
##   theme_skellam()          ggplot2 theme for the package's figures
##   skellam_palette()        the colour roles used in all plots
##   scale_colour_skellam(), scale_fill_skellam()  scales using the palette
##   plot()/autoplot() methods for "skellam_fit" and "skellam_reg":
##                            rootogram, pmf overlay, PIT histogram, QQ plot
##
## Internal functions:
##   .skellam_diag_data()     observed and expected counts, PIT, residuals
##
## ggplot2 is suggested, not required: the functions stop with a clear
## message if it is not installed. Every plot is returned as a ggplot object
## so it can be modified with + like any other ggplot.
##
## Original author: Montasser Ghachem and Oguz Ersan, 2026
## License:         GPL (>= 2)
## Package:         skellam, https://github.com/monty-se/skellam
## --------------------------------------------------------------------------

#' @importFrom utils globalVariables
NULL

utils::globalVariables(c("value", "observed", "expected", "ymin", "ymax", "role",
                         "u", "height", "theoretical", "sample", "lo", "hi"))

.need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting; install it with install.packages('ggplot2')",
         call. = FALSE)
}

#' Colours and Theme for skellam Plots
#'
#' The colour roles and the \pkg{ggplot2} theme used by the package's plots.
#' The colours come from the Okabe-Ito palette, which stays distinguishable
#' for the common forms of colour blindness and in greyscale; colour is never
#' the only cue (observed and fitted values also differ in shape and line type).
#'
#' @param base_size Base font size, in points.
#' @param base_family Base font family.
#' @param ... Passed to \code{\link[ggplot2]{scale_colour_manual}} or
#'   \code{\link[ggplot2]{scale_fill_manual}}.
#'
#' @return \code{skellam_palette()} returns a named character vector of
#'   colours (\code{observed}, \code{fitted}, \code{reference},
#'   \code{secondary}, \code{neutral}, and \code{band}, a light tint of the
#'   reference colour for reference bands; it is opaque so figures print and
#'   export the same on every graphics device); \code{theme_skellam()} a ggplot2 theme;
#'   the scale functions ggplot2 scales.
#'
#' @examples
#' skellam_palette()
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   library(ggplot2)
#'   d <- data.frame(x = -6:6, p = dskellam(-6:6, 3, 2))
#'   ggplot(d, aes(x, p)) +
#'     geom_col(fill = skellam_palette()["fitted"]) +
#'     theme_skellam()
#' }
#' @export
skellam_palette <- function() {
  c(observed = "#4D4D4D", fitted = "#0072B2", reference = "#D55E00",
    secondary = "#009E73", neutral = "#999999", band = "#F7DCCB")
}

#' @rdname skellam_palette
#' @export
theme_skellam <- function(base_size = 11, base_family = "") {
  .need_ggplot2()
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.3),
      axis.line.x        = ggplot2::element_line(colour = "#4D4D4D", linewidth = 0.3),
      axis.ticks         = ggplot2::element_line(colour = "#4D4D4D", linewidth = 0.3),
      axis.title         = ggplot2::element_text(size = ggplot2::rel(0.95)),
      plot.title         = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.05), hjust = 0),
      plot.subtitle      = ggplot2::element_text(colour = "#666666", size = ggplot2::rel(0.9)),
      plot.caption       = ggplot2::element_text(colour = "#999999", size = ggplot2::rel(0.75), hjust = 1),
      legend.position    = "top",
      legend.justification = "left",
      legend.title       = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(6, 10, 6, 6),
      strip.text         = ggplot2::element_text(hjust = 0, face = "plain"),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

#' @rdname skellam_palette
#' @export
scale_colour_skellam <- function(...) {
  .need_ggplot2()
  ggplot2::scale_colour_manual(values = skellam_palette(), ...)
}

#' @rdname skellam_palette
#' @export
scale_fill_skellam <- function(...) {
  .need_ggplot2()
  ggplot2::scale_fill_manual(values = skellam_palette(), ...)
}

# Observed and expected counts on the integer support, the non-randomised PIT
# histogram and randomised quantile residuals, for data y with rates l1, l2
# (scalars, or one pair per observation).
.skellam_diag_data <- function(y, l1, l2, bins = 10) {
  n <- length(y)
  l1 <- rep_len(l1, n); l2 <- rep_len(l2, n)
  lo <- min(min(y), min(qskellam(5e-4, l1, l2)))
  hi <- max(max(y), max(qskellam(1 - 5e-4, l1, l2)))
  k <- lo:hi
  expected <- vapply(k, function(v) sum(dskellam(rep(v, n), l1, l2)), 0)
  observed <- tabulate(y - lo + 1L, nbins = length(k))
  # non-randomised PIT (Czado, Gneiting and Held 2009)
  Fm <- pskellam(y - 1, l1, l2); Fy <- pskellam(y, l1, l2)
  grid <- seq(0, 1, length.out = bins + 1)
  Fbar <- vapply(grid, function(g) mean(pmin(pmax((g - Fm) / pmax(Fy - Fm, 1e-300), 0), 1)), 0)
  pit <- data.frame(u = (grid[-1] + grid[-(bins + 1)]) / 2, height = diff(Fbar) * bins)
  band <- stats::qbinom(c(0.025, 0.975), n, 1 / bins) / (n / bins)
  # randomised quantile residuals (Dunn and Smyth 1996)
  qres <- stats::qnorm(Fm + stats::runif(n) * (Fy - Fm))
  list(counts = data.frame(value = k, observed = observed, expected = expected),
       pit = pit, band = band, qres = qres, n = n)
}

.skellam_diag_plot <- function(dd, type, style, subtitle, envelope) {
  .need_ggplot2()
  pal <- skellam_palette()
  cnt <- dd$counts
  int_breaks <- function(lim) {
    b <- pretty(lim, n = 7)
    b[b == round(b)]
  }
  if (type == "rootogram") {
    cnt$ymax <- sqrt(cnt$expected)
    cnt$ymin <- if (style == "hanging") sqrt(cnt$expected) - sqrt(cnt$observed) else 0
    if (style == "standing") cnt$ymax <- sqrt(cnt$observed)
    p <- ggplot2::ggplot(cnt) +
      ggplot2::geom_hline(yintercept = 0, colour = pal[["observed"]], linewidth = 0.3) +
      ggplot2::geom_rect(ggplot2::aes(xmin = value - 0.4, xmax = value + 0.4, ymin = ymin, ymax = ymax),
                         fill = "grey70", colour = NA) +
      ggplot2::geom_line(ggplot2::aes(value, sqrt(expected)), colour = pal[["fitted"]], linewidth = 0.6) +
      ggplot2::geom_point(ggplot2::aes(value, sqrt(expected)), colour = pal[["fitted"]], size = 1.6) +
      ggplot2::scale_x_continuous(breaks = int_breaks) +
      ggplot2::labs(x = "Value", y = "sqrt(Frequency)", subtitle = subtitle)
    if (style == "hanging")
      p <- p + ggplot2::geom_hline(yintercept = c(-1, 1), colour = pal[["reference"]],
                                   linetype = "dashed", linewidth = 0.4)
    return(p + theme_skellam())
  }
  if (type == "pmf") {
    n <- dd$n
    long <- rbind(data.frame(value = cnt$value, p = cnt$observed / n, role = "Observed"),
                  data.frame(value = cnt$value, p = cnt$expected / n, role = "Skellam fit"))
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data = long[long$role == "Observed", ],
                        ggplot2::aes(value, p, fill = role), width = 0.8) +
      ggplot2::geom_line(data = long[long$role == "Skellam fit", ],
                         ggplot2::aes(value, p, colour = role), linetype = "dashed", linewidth = 0.5) +
      ggplot2::geom_point(data = long[long$role == "Skellam fit", ],
                          ggplot2::aes(value, p, colour = role), size = 1.8) +
      ggplot2::scale_fill_manual(values = c(Observed = "grey70")) +
      ggplot2::scale_colour_manual(values = c(`Skellam fit` = pal[["fitted"]])) +
      ggplot2::guides(fill = ggplot2::guide_legend(order = 1), colour = ggplot2::guide_legend(order = 2)) +
      ggplot2::scale_x_continuous(breaks = int_breaks) +
      ggplot2::labs(x = "Value", y = "Probability", subtitle = subtitle)
    return(p + theme_skellam())
  }
  if (type == "pit") {
    p <- ggplot2::ggplot(dd$pit) +
      ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = dd$band[1], ymax = dd$band[2],
                        fill = pal[["band"]]) +
      ggplot2::geom_col(ggplot2::aes(u, height), width = 1 / nrow(dd$pit) * 0.92, fill = "grey60") +
      ggplot2::geom_hline(yintercept = 1, colour = pal[["fitted"]], linewidth = 0.6) +
      ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      ggplot2::labs(x = "PIT", y = "Density", subtitle = subtitle)
    return(p + theme_skellam())
  }
  # qq
  q <- sort(dd$qres)
  n <- length(q)
  th <- stats::qnorm(stats::ppoints(n))
  df <- data.frame(theoretical = th, sample = q)
  p <- ggplot2::ggplot(df)
  if (envelope > 0) {
    sims <- apply(matrix(stats::rnorm(n * envelope), n), 2, sort)
    df$lo <- apply(sims, 1, stats::quantile, 0.025); df$hi <- apply(sims, 1, stats::quantile, 0.975)
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_ribbon(ggplot2::aes(theoretical, ymin = lo, ymax = hi), fill = pal[["band"]])
  }
  p <- p + ggplot2::geom_abline(slope = 1, intercept = 0, colour = pal[["fitted"]], linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(theoretical, sample), colour = pal[["observed"]], size = 1.2) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Theoretical quantiles", y = "Quantile residuals", subtitle = subtitle)
  p + theme_skellam()
}

#' Diagnostic Plots for Skellam Fits
#'
#' Plots comparing data with a fitted Skellam distribution or Skellam
#' regression. All plots are \pkg{ggplot2} objects, so they can be modified
#' with \code{+}; \pkg{ggplot2} must be installed.
#'
#' @param x,object A fit from \code{\link{skellam.fit}} or
#'   \code{\link{skellam.reg}}.
#' @param type \describe{
#'   \item{\code{"rootogram"}}{(default) Hanging rootogram (Kleiber and Zeileis
#'     2016): bars of \eqn{\sqrt{\mathrm{observed}}}{sqrt(observed)} hang from the
#'     curve of \eqn{\sqrt{\mathrm{expected}}}{sqrt(expected)}; a bar ending
#'     below 0 means too few expected, above 0 too many. Dashed lines at
#'     \eqn{\pm 1}{+-1} are a rough guide.}
#'   \item{\code{"pmf"}}{Observed relative frequencies (bars) and fitted
#'     probabilities (points; averaged over observations for a regression).}
#'   \item{\code{"pit"}}{Non-randomised probability integral transform
#'     histogram for count data (Czado, Gneiting and Held 2009); flat at 1
#'     when the model is right, with a shaded 95\% band for a uniform
#'     histogram of the same size.}
#'   \item{\code{"qq"}}{Normal QQ plot of randomised quantile residuals (Dunn
#'     and Smyth 1996), with a pointwise 95\% envelope. These residuals use
#'     random numbers; set a seed for reproducible plots.}
#' }
#' @param style For the rootogram: \code{"hanging"} (default) or
#'   \code{"standing"}.
#' @param bins Number of bins of the PIT histogram.
#' @param envelope Number of simulated normal samples for the QQ envelope
#'   (0 for none).
#' @param ... Unused.
#'
#' @return A ggplot object.
#'
#' @references
#' Czado, C., Gneiting, T. and Held, L. (2009) Predictive model assessment for
#' count data. \emph{Biometrics} \bold{65}(4), 1254-1261.
#'
#' Dunn, P. K. and Smyth, G. K. (1996) Randomized quantile residuals.
#' \emph{Journal of Computational and Graphical Statistics} \bold{5}(3), 236-244.
#'
#' Kleiber, C. and Zeileis, A. (2016) Visualizing count data regressions using
#' rootograms. \emph{The American Statistician} \bold{70}(3), 296-303.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   fit <- skellam.fit(rpois(500, 6) - rpois(500, 4))
#'   plot(fit)                     # hanging rootogram
#'   plot(fit, type = "pmf")
#'   plot(fit, type = "pit")
#'   plot(fit, type = "qq") + ggplot2::ggtitle("Quantile residuals")
#' }
#' @name skellam-plots
NULL

#' @rdname skellam-plots
#' @export
plot.skellam_fit <- function(x, type = c("rootogram", "pmf", "pit", "qq"),
                             style = c("hanging", "standing"), bins = 10, envelope = 99, ...) {
  type <- match.arg(type); style <- match.arg(style)
  l <- x$estimate
  dd <- .skellam_diag_data(x$x, l[[1]], l[[2]], bins)
  sub <- sprintf("lambda1 = %s, lambda2 = %s, n = %d", format(l[[1]], digits = 3),
                 format(l[[2]], digits = 3), x$nobs)
  .skellam_diag_plot(dd, type, style, sub, envelope)
}

#' @rdname skellam-plots
#' @export
plot.skellam_reg <- function(x, type = c("rootogram", "pmf", "pit", "qq"),
                             style = c("hanging", "standing"), bins = 10, envelope = 99, ...) {
  type <- match.arg(type); style <- match.arg(style)
  dd <- .skellam_diag_data(x$y, x$fitted.values$lambda1, x$fitted.values$lambda2, bins)
  sub <- sprintf("Skellam regression, n = %d", x$nobs)
  .skellam_diag_plot(dd, type, style, sub, envelope)
}

#' @rdname skellam-plots
#' @exportS3Method ggplot2::autoplot
autoplot.skellam_fit <- function(object, ...) plot.skellam_fit(object, ...)

#' @rdname skellam-plots
#' @exportS3Method ggplot2::autoplot
autoplot.skellam_reg <- function(object, ...) plot.skellam_reg(object, ...)
