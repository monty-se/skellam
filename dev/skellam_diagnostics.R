## =============================================================================
## skellam_diagnostics.R  --  an annotated health check of the skellam package
## =============================================================================
##
## HOW TO RUN
##   From the root of your clone (the folder that contains DESCRIPTION):
##       Rscript dev/skellam_diagnostics.R
##   or, in RStudio with the skellam project open:
##       source("dev/skellam_diagnostics.R")
##
##   When run inside the clone, the script tests the code in R/ as it is on
##   disk right now (no install needed), so after every fix you can rerun it
##   and watch FAILs turn into PASSes. Outside the clone it tests the
##   installed package. Set SKELLAM_SRC=/some/path to test another checkout.
##
## WHAT IT PRINTS
##   For every check: PASS / FAIL / ERROR, what we expected, what we got.
##   For every FAIL or ERROR it also explains, in three short paragraphs:
##     WHY IT MATTERS  - what a user would experience, and which convention
##                       of R's own d/p/q/r functions is being broken;
##     CAUSE           - where in the code it goes wrong (file:line, current
##                       GitHub version) and the mechanism;
##     FIX             - what the correct code has to do.
##   Set EXPLAIN = "all" (below, or as an environment variable) to read the
##   explanations for passing checks too.
##
##   The IDs (N1, E3, ...) match the issue numbers in the review document.
##   Results are also saved to dev/diagnostics_results.csv.
## =============================================================================

EXPLAIN <- Sys.getenv("EXPLAIN", "fail")      # "fail" or "all"

## ---------------------------------------------------------------- load code
src <- Sys.getenv("SKELLAM_SRC")
in_clone <- file.exists("DESCRIPTION") && dir.exists("R") &&
  any(grepl("^Package:\\s*skellam\\s*$", readLines("DESCRIPTION", warn = FALSE)))
if (!nzchar(src) && in_clone) src <- "."
if (nzchar(src)) {
  for (f in list.files(file.path(src, "R"), pattern = "\\.[rR]$", full.names = TRUE)) source(f)
  cat("Testing the source files in", normalizePath(file.path(src, "R")), "\n")
} else {
  suppressPackageStartupMessages(library(skellam))
  cat("Testing the installed skellam", as.character(packageVersion("skellam")), "\n")
}
cat("R", R.version$major, ".", R.version$minor, " on ", R.version$platform, "\n", sep = "")
out_dir <- if (dir.exists("dev")) "dev" else "."
options(width = 100)

## ---------------------------------------------------------------- helpers
say <- function(label, text) {
  if (!nzchar(text)) return(invisible())
  lines <- strwrap(text, width = 88, initial = sprintf("      %-15s", label),
                   prefix = strrep(" ", 21))
  cat(lines, sep = "\n")
}
section <- function(title, intro) {
  cat("\n", strrep("=", 92), "\n", title, "\n", strrep("=", 92), "\n", sep = "")
  cat(strwrap(intro, width = 92, prefix = "  "), sep = "\n")
  cat("\n")
}
## run(): evaluate an expression and capture its value, warnings, error, time.
run <- function(expr) {
  old <- getOption("warn"); warns <- character(); t0 <- proc.time()[["elapsed"]]
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "check_error")),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  options(warn = old)
  list(value = val, warnings = unique(warns), time = proc.time()[["elapsed"]] - t0)
}
fmt <- function(v) {
  if (inherits(v, "check_error")) return(paste("ERROR:", v$msg))
  if (is.list(v)) v <- unlist(v)
  if (is.null(v)) return("NULL")
  s <- paste(format(v, digits = 6), collapse = " ")
  if (!is.null(names(v)) && length(v) <= 4) s <- paste(paste0(names(v), "=", format(v, digits = 6)), collapse = " ")
  if (nchar(s) > 70) paste0(substr(s, 1, 67), "...") else s
}
results <- list()
check <- function(id, title, r, ok, expected, why = "", cause = "", fix = "") {
  status <- if (inherits(r$value, "check_error")) "ERROR" else
    if (isTRUE(tryCatch(ok(r$value, r$warnings), error = function(e) FALSE))) "PASS" else "FAIL"
  got <- fmt(r$value)
  if (length(r$warnings))
    got <- paste0(got, "   [", length(r$warnings), " warning(s), first: ", substr(r$warnings[1], 1, 55), "]")
  cat(sprintf("%-4s %-5s %s\n", id, status, title))
  say("expected:", expected)
  say("got:", got)
  if (status != "PASS" || EXPLAIN == "all") {
    say("WHY IT MATTERS", why)
    if (status != "PASS") { say("CAUSE", cause); say("FIX", fix) }
  }
  cat("\n")
  results[[length(results) + 1]] <<- data.frame(ID = id, Status = status, Check = title,
                                                Expected = expected, Got = got, stringsAsFactors = FALSE)
}
## Exact reference: P(X = x) = sum_k P(Y1 = k + x) P(Y2 = k), summed on the log
## scale ("log-sum-exp") so that nothing underflows. Slow but simple and exact
## to about 1e-13; this is the yardstick for every accuracy check below.
lp_exact <- function(x, l1, l2) {
  if (l1 == 0) return(dpois(-x, l2, log = TRUE))
  if (l2 == 0) return(dpois(x, l1, log = TRUE))
  k <- max(0, -x):ceiling(max(l1, l2) + 60 * sqrt(max(l1, l2)) + abs(x) + 100)
  m <- dpois(k + x, l1, log = TRUE) + dpois(k, l2, log = TRUE)
  mm <- max(m); if (!is.finite(mm)) return(-Inf)
  mm + log(sum(exp(m - mm)))
}

## =============================================================================
section("0. WHAT THE PACKAGE COMPUTES, AND THE RULES IT SHOULD FOLLOW", paste(
"Let Y1 ~ Poisson(lambda1) and Y2 ~ Poisson(lambda2) be independent counts (0, 1, 2, ...).",
"The Skellam variable is their difference X = Y1 - Y2: home goals minus away goals, up-ticks",
"minus down-ticks, spikes in a stimulus window minus spikes in a baseline window. A difference",
"of two whole numbers is a whole number, so X can only take the values ..., -2, -1, 0, 1, 2, ...",
"In dskellam(x, lambda1, lambda2) the argument x is a value X might take, and the function",
"answers 'what is P(X = x)?'. For x = 2.5 the honest answer is 0: X never equals 2.5.",
"\n\nThe four functions: d = P(X = x); p = P(X <= q); q = the inverse of p, i.e. the smallest",
"integer x with P(X <= x) >= p; r = random draws.",
"\n\nEvery d/p/q/r function in base R (dpois, pbinom, qnbinom, ...) obeys the same rules, and",
"users and other packages rely on them: (1) an NA input gives an NA output, never an error;",
"(2) an impossible parameter (negative rate, probability outside [0, 1]) gives NaN with the",
"warning 'NaNs produced'; (3) outside the support the density is 0 (-Inf on the log scale);",
"(4) q undoes p exactly: q(p(x)) == x, for both tails and on the log scale; q(0) and q(1) are",
"the ends of the support; (5) log = TRUE returns the log of the same number, computed so that it",
"stays finite even when the number itself is too small for the computer; (6) vectors of",
"different lengths are recycled, with a warning when the lengths are not multiples;",
"(7) r returns the same type as the base generator (integer for counts).",
"\n\nOne fact about computers explains many checks: a double can hold numbers down to about",
"1e-308 (with reduced precision down to 5e-324); anything smaller becomes exactly 0. Its",
"logarithm, however, is an ordinary number: log(1e-1000) = -2302.6. That is why likelihoods",
"are always computed on the log scale."))

## =============================================================================
section("1. INPUT HANDLING", paste(
"These checks feed the functions the kind of input real data contain: missing values,",
"non-integers, impossible parameters, vectors of mismatched lengths. None of them is exotic;",
"a regression on a data frame with one missing row hits the first one."))

check("N1", "dskellam(NA, 3, 2) returns NA",
  run(dskellam(NA, 3, 2)), function(v, w) is.na(v), "NA  (dpois(NA, 3) is NA)",
  why = paste("Rule 1. Data have missing values, and the standard idiom",
    "sum(dskellam(x, l1, l2, log = TRUE), na.rm = TRUE) needs NA to come back as NA."),
  cause = paste("R/dskellam.r:131-132. `nonint <- x != trunc(x)` is NA when x is NA, and",
    "`if (any(nonint))` then becomes `if (NA)`, which R refuses: 'missing value where",
    "TRUE/FALSE needed'. The function stops before computing anything."),
  fix = paste("Leave NAs out of the integer test: `nonint <- !is.na(x) & abs(x - round(x)) > 1e-7 * pmax(1, abs(x))`,",
    "compute only on the non-missing entries, and put NA in the result where x was NA."))

check("N1", "dskellam(c(0, NA), 3, 2): one NA does not break the whole vector",
  run(dskellam(c(0, NA), 3, 2)), function(v, w) length(v) == 2 && is.na(v[2]) && abs(v[1] - 0.1677219) < 1e-6,
  "0.1677219 NA",
  why = paste("This is how N1 bites in practice: a log-likelihood over 10,000 observations",
    "crashes because one of them is missing, and the error message does not say why."),
  cause = "Same line as above: the if() test is on the whole vector.",
  fix = "Same fix.")

check("N11", "dskellam(2.5, 3, 2) is 0, because X is never 2.5",
  run(dskellam(2.5, 3, 2)), function(v, w) v == 0 && length(w) > 0, "0, with a warning  (dpois(2.5, 3) is 0 with a warning)",
  why = paste("Rule 3. X is an integer, so P(X = 2.5) = 0. Returning P(X = 2) instead claims a",
    "positive probability for an impossible value, and a sum of densities over a non-integer",
    "grid then adds up to more than 1. The warning is right (a non-integer x is usually a",
    "mistake), the number is wrong."),
  cause = paste("R/dskellam.r:136 `x <- trunc(x)` replaces 2.5 by 2 after warning, and the",
    "density at 2 is returned. trunc() also rounds toward zero, so -2.7 becomes -2, not -3."),
  fix = paste("Return 0 (or -Inf if log = TRUE) for non-integers. Why not an error? Because",
    "arithmetic produces values like 2.9999999999996 that are meant to be 3; base R treats",
    "anything within 1e-7 (relative) of an integer as that integer, and the fix should do the same."))

check("N12", "qskellam(1.1, 3, 2): a 'probability' above 1 gives NaN",
  run(qskellam(1.1, 3, 2)), function(v, w) is.nan(v) && length(w) > 0, "NaN, with the warning 'NaNs produced'",
  why = paste("Rule 2. There is no 110% quantile. Returning Inf silently turns a user's bug",
    "(e.g. a probability computed as 1 - p with p negative) into a plausible-looking number."),
  cause = paste("qskellam never checks p. R/qskellam.r:25 tests `p + 1.01*eps >= 1` to detect",
    "p = 1, and 1.1 passes that test, so it is treated as p = 1 and the answer is Inf."),
  fix = "Validate at the top: p outside [0, 1] (or > 0 on the log scale) -> NaN with a warning.")

check("N12", "qskellam(0.5, -1, 2): a negative rate gives NaN, not an error",
  run(qskellam(0.5, -1, 2)), function(v, w) is.nan(v), "NaN, with a warning (what dskellam and pskellam already do)",
  why = paste("Rule 2. A Poisson rate cannot be negative. dskellam and pskellam return NaN with",
    "a warning; qskellam crashes, so the three functions disagree on the same input."),
  cause = paste("qskellam has no rate check (compare R/dskellam.r:113-120). The negative rate",
    "makes sqrt(lambda1 + lambda2) NaN in the Cornish-Fisher start (line 46), and the search",
    "loop `while (any(up1))` at line 56 then meets an NA and stops with 'missing value where",
    "TRUE/FALSE needed'."),
  fix = "Use the same rate check as dskellam; set those entries to NaN and leave them out of the search.")

check("N12", "qskellam(NA, 3, 2) returns NA",
  run(qskellam(NA, 3, 2)), function(v, w) is.na(v), "NA",
  why = "Rule 1, for the quantile function.",
  cause = "Same search loop: an NA probability makes the loop condition NA.",
  fix = "Exclude NA entries before the search and return NA for them.")

check("N13", "dskellam(1:3, 1:2, 1) warns about recycling",
  run(dskellam(1:3, 1:2, 1)), function(v, w) any(grepl("multiple", w)),
  "a warning 'longer object length is not a multiple...' (as 1:3 + 1:2 gives)",
  why = paste("Rule 6. Silent recycling hides bugs: if the vector of rates is one element short,",
    "every observation after the first mismatch is paired with the wrong rate and nothing tells you."),
  cause = paste("R/dskellam.r:125 (and pskellam.r:20, qskellam.r:10) warn only if",
    "`all(len/lens != len %/% lens)`, i.e. only if EVERY argument fails to divide the longest",
    "length. The longest argument always divides itself, so the warning can never fire."),
  fix = "Replace all() by any(), and only compare the arguments that have length > 0.")

check("N17", "rskellam returns integers, like rpois",
  run({ set.seed(1); rskellam(5, 3, 2) }), function(v, w) is.integer(v), "an integer vector",
  why = paste("Rule 7. A difference of counts is a whole number; rpois returns integer, and code",
    "that tests is.integer(), or uses identical() against integer data, behaves differently on doubles."),
  cause = "R/rskellam.r:11 starts from `rep(NaN, n)`, a double vector, so the integer draws are converted to double.",
  fix = "Start from rep(NA_integer_, n) (rpois itself returns NA for invalid rates).")

check("A10", "dskellam keeps the names of x, like dpois",
  run(dskellam(c(a = 0, b = 1), 3, 2)), function(v, w) identical(names(v), c("a", "b")), "a=0.167722 b=0.183023",
  why = "Convenience: dpois(c(a = 0), 3) is named, and named results make tables self-describing.",
  cause = "The result is built in a fresh vector `ret` (R/dskellam.r:137) and names(x) are never copied.",
  fix = "`names(ret) <- names(x)` at the end (when x is the longest argument).")

## =============================================================================
section("2. THE QUANTILE FUNCTION", paste(
"For a discrete variable the quantile is defined by an inequality: qskellam(p) is the smallest",
"integer x with P(X <= x) >= p. Because of the '>=', when p is exactly one of the cdf values",
"(an 'atom'), the answer is that x itself. This is what makes q(p(x)) == x hold, and it is the",
"most sensitive test of a quantile function: being off by one at the atoms is the typical bug.",
"With lower.tail = FALSE the definition becomes: the smallest x with P(X > x) <= p."))

check("N2", "qskellam(c(0, 0.5, 1), 3, 2) gives -Inf, the median, Inf",
  run(qskellam(c(0, 0.5, 1), 3, 2)), function(v, w) identical(v[c(1, 3)], c(-Inf, Inf)),
  "-Inf 1 Inf",
  why = paste("Rule 4. With both rates positive X is unbounded in both directions, so the 0%",
    "quantile is -Inf and the 100% quantile is +Inf. Returning -Inf for p = 1 is the wrong end",
    "of the real line; any code that uses qskellam(1, ...) as an upper limit breaks."),
  cause = paste("R/qskellam.r:28-29. `p[bdry] == 0` is a logical vector with one entry per",
    "BOUNDARY value (here 2 entries: TRUE, FALSE), but it is used to index `ret`, which has one",
    "entry per INPUT (here 3). R recycles a short logical index, so TRUE, FALSE becomes",
    "TRUE, FALSE, TRUE and -Inf is written into positions 1 and 3. (Line 29 then writes into",
    "position 2, which the search later overwrites with the median.) This is also the source",
    "of the 'number of items to replace' warning."),
  fix = "Index with full-length logicals: `ret[bdry & p == 0] <- ...`, `ret[bdry & p > 0] <- ...`.")

x <- -8:10
check("N3", "upper-tail round trip: q(p(x, upper), upper) == x for x = -8..10",
  run(qskellam(pskellam(x, 3, 2, lower.tail = FALSE), 3, 2, lower.tail = FALSE) - x),
  function(v, w) all(v == 0), "all 19 differences are 0",
  why = paste("Rule 4. Every upper-tail quantile at an atom is one too large. Upper-tail",
    "quantiles are what you use for critical values ('reject if X exceeds q(0.05, upper)'), so",
    "this shifts every such test by one unit."),
  cause = paste("R/qskellam.r:53, `p <- p * (1 - 64 * eps)`. Because p(x) is computed in floating",
    "point, a p that comes from pskellam may be off in the 16th digit. The line shrinks p a tiny",
    "bit so that 'F(x) >= p' still holds at the atom. That is the right direction for the lower",
    "tail. For the upper tail the condition is 'P(X > x) <= p', and shrinking p makes it FAIL at",
    "the atom, so the search moves one step further. One line, wrong direction."),
  fix = "Shrink p for the lower tail, enlarge it for the upper tail: `p * (1 + 64 * eps)` when lower.tail = FALSE. (Tested: this alone makes all round trips pass.)")

check("N3", "lower-tail round trip on the log scale: q(p(x, log.p = TRUE), log.p = TRUE) == x",
  run(qskellam(pskellam(x, 3, 2, log.p = TRUE), 3, 2, log.p = TRUE) - x),
  function(v, w) all(v == 0), "all 19 differences are 0",
  why = "Rule 4 on the log scale, which is how extreme quantiles (p = 1e-50) have to be requested.",
  cause = paste("The same line 53. With log.p = TRUE, p is a negative number such as -1.2;",
    "multiplying it by (1 - 64 eps) makes it slightly LESS negative, i.e. a slightly LARGER",
    "probability, the opposite of what was intended."),
  fix = "On the log scale add log1p(-64 eps) (lower tail) or log1p(64 eps) (upper tail) instead of multiplying.")

check("N3", "upper and lower tails agree away from atoms: q(.3, upper) == q(.7, lower)",
  run(c(qskellam(0.3, 3, 2, lower.tail = FALSE), qskellam(0.7, 3, 2))),
  function(v, w) v[1] == v[2], "2 2  (scipy gives 2 and 2)",
  why = paste("A control check. Since P(X > x) = 1 - P(X <= x), both calls ask the same",
    "question. It passes because 0.3 is not an atom: the bug above only shows at atoms,",
    "which is why random-p tests miss it and round-trip tests catch it."))

x <- -7:0
check("N4", "lambda1 = 0 round trip: q(p(x, 0, 1), 0, 1) == x",
  run(qskellam(pskellam(x, 0, 1), 0, 1) - x), function(v, w) all(v == 0), "all 8 differences are 0",
  why = paste("With lambda1 = 0 we have X = -Y2, a Poisson turned upside down. Rule 4 again.",
    "This is the special case the package handles separately 'for speed', and it is wrong at every atom."),
  cause = paste("R/qskellam.r:21 uses `-qpois(p, lambda2, lower.tail = !lower.tail)`. Turning a",
    "distribution upside down turns 'X <= x' into 'Y >= -x'. For a continuous variable '>=' and",
    "'>' are the same; for counts they are not: P(Y >= 3) = P(Y > 2). The reflection uses the",
    "Poisson upper tail P(Y > y) where it needs P(Y >= y), which is off by one exactly at the atoms."),
  fix = paste("Either drop the shortcut and let the general search handle it, or compute",
    "y <- qpois(1 - p, lambda2) and step one further when ppois(y, lambda2) equals 1 - p (to within the tolerance)."))

r1 <- run({ options(verbose = FALSE); qskellam(c(0, 0.5, 1), 3, 0) })
r2 <- run({ options(verbose = TRUE); v <- qskellam(c(0, 0.5, 1), 3, 0); options(verbose = FALSE); v })
check("N14", "qskellam gives the same answer whatever options(verbose) is",
  list(value = list(verbose_FALSE = r1$value, verbose_TRUE = r2$value), warnings = c(r1$warnings, r2$warnings)),
  function(v, w) identical(v[[1]], v[[2]]), "0 3 Inf  0 3 Inf",
  why = paste("options(verbose) is a global R switch for printing extra messages; people turn it",
    "on while debugging something unrelated. A numerical result must never depend on it. Here",
    "qskellam(1, 3, 0) becomes 0 instead of Inf when verbose is on."),
  cause = paste("R/qskellam.r:17-23 was a developer's self-test: with verbose on, the Poisson",
    "shortcut is skipped and the general code runs, which then hits the N2 indexing bug at line 29."),
  fix = "Delete the getOption('verbose') branch; that cross-check belongs in the test suite, not in the function.")

## =============================================================================
section("3. TAILS AND THE SADDLEPOINT FALLBACKS", paste(
"Far from the centre, probabilities become smaller than the smallest double and turn into 0.",
"The package then switches to a 'saddlepoint approximation', an accurate analytic formula for",
"tail probabilities. The idea is good; the checks below show where the switch is wired wrongly.",
"These cases matter for likelihoods: one outlying observation with log-density -Inf makes the",
"whole log-likelihood -Inf and the optimizer gives up."))

check("N5", "dskellam(500, 3, 2, log = TRUE) is finite",
  run(dskellam(500, 3, 2, log = TRUE)), function(v, w) is.finite(v) && abs(v - lp_exact(500, 3, 2)) < 1e-6,
  sprintf("%.4f  (exact log-sum computed by this script)", lp_exact(500, 3, 2)),
  why = paste("Rule 5. P(X = 500) is about 10^-898: too small for a double, so dskellam(500, 3, 2)",
    "= 0 is acceptable. But its logarithm is an ordinary number, -2067, and that is what a",
    "likelihood needs. pskellam(500, 3, 2, lower.tail = FALSE, log.p = TRUE) does return -2072,",
    "so the package already knows how; dskellam does not."),
  cause = paste("R/dskellam.r:149 computes log(besselI(...)). The Bessel value itself (~1e-900)",
    "underflows to 0 before the log is taken, so log(0) = -Inf. The fallback at lines 155-171",
    "then computes the saddlepoint density on the ordinary scale, exp(K - x*s), which underflows",
    "too, and it never looks at the log argument."),
  fix = "Compute the saddlepoint (or an asymptotic Bessel formula) on the log scale: K - x*s - 0.5*log(2*pi*K2) + log1p(correction).")

check("N5", "dskellam(-500, 3, 2, log = TRUE) is finite",
  run(dskellam(-500, 3, 2, log = TRUE)), function(v, w) is.finite(v) && abs(v - lp_exact(-500, 3, 2)) < 1e-6,
  sprintf("%.4f", lp_exact(-500, 3, 2)), why = "Same as above, left tail.", cause = "Same.", fix = "Same.")

check("N5b", "dskellam(0, 1e5, 1e5, log = TRUE) is a LOG probability",
  run(dskellam(0, 1e5, 1e5, log = TRUE)), function(v, w) abs(v - lp_exact(0, 1e5, 1e5)) < 1e-8,
  sprintf("%.4f  (log of 0.000892)", lp_exact(0, 1e5, 1e5)),
  why = paste("Rule 5, and the most dangerous result in this file. For rates above about 5e4 the",
    "function returns the probability itself (0.000892) when asked for its log (-7.02). A",
    "log-likelihood built from it is silently wrong by orders of magnitude; nothing warns. This",
    "is what produces the errors of size ~39 in the accuracy table at the end."),
  cause = paste("Two things combine. (a) R's besselI() refuses arguments above 1e5 and silently",
    "returns 0 (try besselI(2e5, 0, TRUE)); the argument here is 2*sqrt(1e5*1e5) = 2e5. So the",
    "exact path gives log(0) = -Inf and the fallback is used. (b) The fallback, R/dskellam.r:170,",
    "always writes exp(...), the ordinary-scale density, into the result, even when log = TRUE."),
  fix = "Same as N5: a log-scale fallback that respects `log`. Long term: an asymptotic formula for log I_nu(z) that has no argument limit.")

check("N6", "the fallback pairs each x with its own saddlepoint when only part of a vector underflows",
  run(dskellam(c(0, 400), 3, 2)[2] - dskellam(400, 3, 2)), function(v, w) v == 0, "0",
  why = paste("A regression test. R/dskellam.r:170 uses x[usp] where the neighbouring lines use the",
    "subset xx; with mixed vectors that could mix up observations. It passes here because both",
    "values in the fallback are 0 anyway; keep the check, and fix the indexing so it cannot happen."),
  cause = "R/dskellam.r:170: `x[usp]` should be `xx[usp]`.", fix = "`xx[usp]`.")

check("N9", "pskellam(-200, 3, 2) is a probability, i.e. not negative",
  run(pskellam(-200, 3, 2)), function(v, w) v >= 0, "1.43e-317 (scipy)",
  why = paste("A probability below 0 is impossible. This value lives in the 'subnormal' range",
    "(below 2.2e-308) where doubles keep only a few digits, and the result differs between",
    "machines: it came out negative on Linux and correct on Windows. Fragile either way."),
  cause = paste("The saddlepoint cdf (R/pskellam.sp.r:15) is pnorm(-w) - dnorm(w)*(1/w - 1/u): a",
    "difference of two nearly equal tiny numbers. In the subnormal range the rounding of each",
    "term decides the sign."),
  fix = "Compute tails on the log scale and exponentiate at the end; clamp at 0 as a last guard.")

check("N9", "pskellam(-200, 3, 2, log.p = TRUE) is finite and correct",
  run(pskellam(-200, 3, 2, log.p = TRUE)), function(v, w) is.finite(v) && abs(v - (-729.5627)) < 1e-3,
  "-729.5627",
  why = "The robust version of the check above: on the log scale there is nothing subnormal about it.",
  cause = "See above.", fix = "See above.")

check("N10", "pskellam(-1, 0, 0) is 0: with both rates 0, X is always 0",
  run(pskellam(-1, 0, 0)), function(v, w) identical(v, 0), "0",
  why = "A degenerate but legal input (two Poisson(0) variables); the answer is exactly 0, not undefined.",
  cause = paste("pchisq at R/pskellam.r:30 correctly returns 0. Line 36 then treats every result",
    "below 1e-308 as 'probably underflowed' and recomputes it with the saddlepoint, which divides",
    "by lambda2 = 0 and returns NaN. The code cannot tell an exact 0 from an underflowed one."),
  fix = "Handle zero rates before anything else (X = Y1 or X = -Y2 is a Poisson), and decide on fallbacks on the log scale.")

check("N15", "dskellam.sp(0, 3, 2, log = TRUE) is a log probability",
  run(dskellam.sp(0, 3, 2, log = TRUE)), function(v, w) v < 0, "about -1.80 (the log of 0.165)",
  why = "An exported function that accepts an argument and silently ignores it.",
  cause = "R/dskellam.sp.r never uses `log`; the `terms` branch at line 10 is dead code (terms is fixed at 1).",
  fix = "`if (log) log(ret) else ret`, or better, compute on the log scale throughout.")

check("N15", "dskellam.sp with a zero rate is finite",
  run(dskellam.sp(2, 3, 0)), function(v, w) is.finite(v), "about dpois(2, 3) = 0.224",
  why = "A control check; the saddlepoint survives one zero rate here.")

## =============================================================================
section("4. LARGE RATES", paste(
"Rates of 1e4 to 1e6 are ordinary in finance (ticks per day) and neuroscience (spikes per",
"session). Two checks have exact answers we can compare against: probabilities must add up to",
"1, and when lambda1 = lambda2 the distribution is symmetric, so P(X <= 0) = (1 + P(X = 0)) / 2."))

for (L in list(c(1e4, 1e4), c(1e5, 1e5), c(1e6, 1e6), c(1e4, 1), c(1e6, 1e3))) {
  s <- sqrt(sum(L)); m <- L[1] - L[2]
  xx <- seq(floor(m - 20 * s), ceiling(m + 20 * s))
  big <- L[1] * L[2] > 2.5e9
  check("N7", sprintf("probabilities add up to 1 at lambda = (%g, %g)", L[1], L[2]),
    run(1 - sum(dskellam(xx, L[1], L[2]))), function(v, w) abs(v) < 1e-12 && length(w) == 0,
    "|1 - sum| below 1e-12, no warnings",
    why = paste("The most basic property of a distribution. A missing mass of 3e-7 sounds small,",
      "but it is not spread evenly: it means the density is systematically biased at every point,",
      "and it is 10 million times larger than the accuracy the package achieves at small rates."),
    cause = if (big) paste("besselI returns 0 for arguments above 1e5 (here 2*sqrt(lambda1*lambda2) >",
      "1e5), so EVERY point goes to the saddlepoint fallback. The fallback, R/dskellam.r:169-170,",
      "multiplies by (1 + c)/2 where c = 1 + correction: that applies only HALF of the standard",
      "second-order correction. The missing half is 1/(16 (lambda1 + lambda2)) at the centre,",
      "which is exactly the 3.1e-7 (at 1e5) and 3.1e-8 (at 1e6) we see. With the full correction",
      "the mass at 1e5 is 1 to within 2e-12 (verified).") else paste("Here lambda1 and lambda2 are very unequal, so",
      "|x| (the Bessel order) is large compared with the argument 2*sqrt(lambda1*lambda2). R's",
      "besselI loses precision in that regime, prints 'precision lost in result' thousands of",
      "times, and returns inaccurate values that are not caught by the fallback test."),
    fix = paste("Short term: apply the full correction, `* (1 + correction)`. Proper fix: a log-scale",
      "uniform asymptotic expansion of log I_nu(z) for large orders and arguments (review, N7)."))
}

for (L in c(1e2, 1e4, 1e6)) {
  r <- run(pskellam(0, L, L) - (1 + exp(lp_exact(0, L, L))) / 2)
  check("N8", sprintf("pskellam(0, %g, %g) satisfies the symmetry identity  [took %.2f s]", L, L, r$time),
    r, function(v, w) abs(v) < 1e-12 && length(w) == 0, "difference below 1e-12, no warnings",
    why = paste("At lambda = 1e6 the true value is 0.50014 and pskellam returns about 0.875: wrong",
      "by 0.37, i.e. completely wrong, with only a convergence warning from inside pchisq. Any",
      "p-value or quantile at these rates is unusable."),
    cause = paste("pskellam uses the identity P(X <= x) = pchisq(2*lambda1, 2(x+1), ncp = 2*lambda2)",
      "(R/pskellam.r:30-34). R's non-central chi-square sums a Poisson-weighted series around",
      "ncp/2 = 1e6 terms and stops at its iteration limit, returning the partial sum. The fallback",
      "on line 36 only fires for non-finite or tiny results, so the wrong 0.875 is returned as if exact."),
    fix = paste("Above lambda1 + lambda2 of about 1e4, do not use pchisq: sum the (log-scale) pmf",
      "from the nearer tail, or use the saddlepoint cdf with a documented accuracy."))
}

r <- run(qskellam(c(1e-10, 0.5, 1 - 1e-10), 1e6, 1e6))
check("N8", sprintf("qskellam at lambda = 1e6 is fast and silent  [took %.2f s]", r$time),
  r, function(v, w) length(w) == 0 && r$time < 0.05, "-8996 0 8996 in under 0.05 s, no warnings",
  why = paste("The answers happen to be right because the Cornish-Fisher starting guess is",
    "excellent and only one comparison is needed. But each comparison calls the broken pskellam,",
    "so speed and correctness both depend on N8."),
  cause = "Every search step calls pskellam, which calls the slow, non-converging pchisq.",
  fix = "Fixing pskellam (above) fixes this.")

## =============================================================================
section("5. ESTIMATION", paste(
"skellam.mle(x) estimates (lambda1, lambda2) from a sample of differences; skellam.reg(y, x)",
"lets both rates depend on covariates. Estimation code has extra duties: it must not change",
"the user's R session, must say when the data cannot support an estimate, and must return an",
"object that R's model functions (coef, vcov, logLik, AIC, summary, predict) understand."))

set.seed(42); xs <- rpois(500, 6) - rpois(500, 4)
options(warn = 0)
r <- run({ try(skellam.mle(xs), silent = TRUE); getOption("warn") })
check("E1", "skellam.mle() leaves the user's warning setting alone",
  r, function(v, w) v == 0, "0 (the value before the call)",
  why = paste("After one call to skellam.mle() every warning in the user's R session is switched",
    "off: log(-1) returns NaN silently, failed model fits go unreported, and the package's own",
    "'NaNs produced' warnings disappear. CRAN policy forbids changing global options without",
    "restoring them."),
  cause = "R/skellam.mle.R:79 `options(warn = -1)`, never undone.",
  fix = "Delete the line. If a specific warning is expected inside the objective, wrap that one call in suppressWarnings().")
options(warn = 0)
r <- run({ set.seed(1); try(skellam.reg(xs, rnorm(500)), silent = TRUE); getOption("warn") })
check("E1", "skellam.reg() leaves the user's warning setting alone (even when it fails)",
  r, function(v, w) v == 0, "0",
  why = "Same as above. Note it happens even when the fit then errors, because the option is set before optimizing.",
  cause = "R/skellam.reg.R:69 `options(warn = -1)`.", fix = "Delete the line.")
options(warn = 0)

check("E4", "the fit is an object that coef() and logLik() understand",
  run(c(coef = !is.null(coef(skellam.mle(xs))), logLik = inherits(logLik(skellam.mle(xs)), "logLik"))),
  function(v, w) all(v), "coef=TRUE logLik=TRUE",
  why = paste("Functions like coef(), vcov(), logLik(), AIC() and summary() are 'generic': they look",
    "at class(fit) and call the matching method, e.g. logLik.skellam_fit(). A plain list has no",
    "class, so AIC(skellam.mle(x)) is an error and model comparison has to be done by hand. There",
    "are also no standard errors at all."),
  cause = "R/skellam.mle.R:83 returns list(iters, loglik, param).",
  fix = "Return a classed object with estimates, vcov (from the observed information), loglik, nobs, and write the methods: this is skellam.fit() in the roadmap.")

check("E2", "an under-dispersed sample (variance < mean) gives non-negative rates",
  run(skellam.mle(c(rep(2, 45), rep(3, 5)))$param),
  function(v, w) all(is.finite(v)) && all(v >= 0), "two finite values >= 0 (here about 2.1 and 0)",
  why = paste("For a Skellam, Var(X) = lambda1 + lambda2 >= |lambda1 - lambda2| = |E(X)|. A sample",
    "with variance below |mean| cannot come from a Skellam with both rates positive; the best",
    "fit puts one rate at 0 (a Poisson). Returning lambda2 = -1 reports a negative Poisson rate,",
    "which is meaningless, and plugging it into dskellam gives NaN."),
  cause = paste("R/skellam.mle.R:71 starts at the moment estimate lambda2 = (var - mean)/2, which",
    "is negative here. The objective (line 74) then takes sqrt of a negative number, gets NaN,",
    "nlm cannot move, and the starting value is returned as the 'estimate'."),
  fix = "Optimize over log(lambda2) (or with optimize() on [0, upper]) and start from max((var - mean)/2, small positive).")

check("E3", "an all-zero sample gives the boundary estimate (0, 0)",
  run(skellam.mle(rep(0, 50))$param), function(v, w) all(is.finite(v)) && all(abs(v) < 1e-6), "0 0",
  why = paste("If every difference is 0, the likelihood is largest when both rates are 0 (then",
    "P(X = 0) = 1). It is a legitimate, if extreme, data set (e.g. a quiet sensor); the function",
    "should return the boundary estimate with a warning that standard errors do not exist."),
  cause = paste("mean = 0 and var = 0 give a start of 0; the objective contains log(1 + mean/theta)",
    "= log(1 + 0/0) = NaN at the very first evaluation, and nlm stops with 'non-finite value'."),
  fix = "Detect var(x) == 0 before optimizing and return the closed-form boundary estimate.")

check("E3", "a constant sample x = 3 gives the boundary estimate (3, 0)",
  run(skellam.mle(rep(3, 50))$param), function(v, w) all(abs(v - c(3, 0)) < 1e-4), "3 0",
  why = "Same situation shifted: every difference equals 3, best explained by lambda1 = 3, lambda2 = 0.",
  cause = "Start (0 - 3)/2 = -1.5 is negative (as in E2); the returned 'estimate' is the start, (1.5, -1.5).",
  fix = "Same as E3 above.")

check("E3", "a sample containing NA gives a clear message",
  run(skellam.mle(c(xs[1:20], NA))), function(v, w) inherits(v, "check_error") && grepl("NA", v$msg),
  "an error that mentions NA, or an na.rm argument",
  why = "The error currently comes from deep inside nlm and does not tell the user that the data contain a missing value.",
  cause = "sum(x) is NA (line 68), so the start is NA and nlm complains about its own parameter.",
  fix = "Check anyNA(x) at the top; add na.rm = FALSE as an argument.")

set.seed(7); z <- rnorm(300)
yb <- rpois(300, exp(6 + 0.2 * z)) - rpois(300, exp(5.9 + 0.2 * z))   # rates around 400
set.seed(3)
r <- run(skellam.reg(yb, z))
check("E5", "skellam.reg works when the rates are around 400",
  r, function(v, w) is.finite(v$loglik) && abs(v$param1[1, 1] - 6) < 0.1 && abs(v$param2[1, 1] - 5.9) < 0.1,
  "finite log-likelihood, intercepts near 6 and 5.9",
  why = "Rates of a few hundred are normal (ticks, spikes, page views). The regression cannot fit them at all.",
  cause = paste("Two problems in R/skellam.reg.R:66, which uses besselI(a, y) WITHOUT expon.scaled.",
    "(a) The starting values are random, rnorm(2p) at line 70, i.e. rates around 1, while the",
    "data reach |y| = 130; at rates near 1, P(X = 130) is about 1e-200 and besselI underflows to",
    "0, so log(0) = -Inf and nlm stops (the message varies by machine and seed: 'function cannot",
    "be evaluated at initial parameters' or 'non-finite finite-difference value').",
    "(b) Even from a good start, the unscaled besselI overflows to Inf once 2*sqrt(lambda1*lambda2)",
    "exceeds about 713, i.e. rates around 360 (try besselI(800, 3))."),
  fix = "Use the log-scale density (log(besselI(a, abs(y), TRUE)) + a, or dskellam(log = TRUE) once fixed) and start from two Poisson GLMs.")

wcov <- seq_along(xs) / 500
r <- run({ a <- skellam.reg(xs, wcov); b <- skellam.reg(xs, wcov); max(abs(a$param1[, 1] - b$param1[, 1])) })
check("E6", "skellam.reg gives the same answer on two identical calls",
  r, function(v, w) v < 1e-8, "0",
  why = paste("Results in a paper must be reproducible. With random starting values, two runs on",
    "the same data stop at slightly different points, or one of them fails (see E5)."),
  cause = "R/skellam.reg.R:70 starts the optimizer at stats::rnorm(2 * p), without a seed.",
  fix = "Deterministic starts from Poisson GLMs on pmax(y, 0) and pmax(-y, 0).")

## =============================================================================
section("6. ACCURACY TABLE: dskellam(log = TRUE) AGAINST THE EXACT SUM", paste(
"For each pair of rates, five points (the centre, +-3 and +-8 standard deviations) are compared",
"with the exact log-sum. max_abs_logerr is the largest |log(package) - log(exact)|; for small",
"values it equals the relative error of the probability, so 1e-14 means 14 correct digits.",
"'warnings' counts besselI precision warnings. Errors of about 39 at large rates are not",
"inaccuracy in the usual sense: they are N5b, the function returning p instead of log(p)."))

grid <- list(c(0.01, 0.01), c(3, 2), c(100, 80), c(1e3, 1e3), c(1e4, 1e4),
             c(1e5, 1e5), c(1e6, 1e6), c(1e4, 1), c(1e6, 1e3), c(50, 0.001))
acc <- do.call(rbind, lapply(grid, function(L) {
  s <- sqrt(sum(L)); m <- L[1] - L[2]
  xp <- unique(round(m + s * c(-8, -3, 0, 3, 8)))
  ex <- sapply(xp, lp_exact, l1 = L[1], l2 = L[2])
  rr <- run(dskellam(xp, L[1], L[2], log = TRUE))
  err <- if (inherits(rr$value, "check_error")) NA else abs(rr$value - ex)
  data.frame(lambda1 = L[1], lambda2 = L[2], max_abs_logerr = signif(max(err[is.finite(ex)]), 3),
             correct_digits = pmax(0, round(-log10(max(err[is.finite(ex)]) + 1e-17))),
             warnings = length(rr$warnings))
}))
print(acc, row.names = FALSE)
cat("\nTarget after the fixes: at least 12 correct digits in every row and no warnings.\n")

## =============================================================================
res <- do.call(rbind, results)
cat("\n", strrep("=", 92), "\nSUMMARY\n", strrep("=", 92), "\n", sep = "")
print(table(factor(res$Status, levels = c("PASS", "FAIL", "ERROR"))))
cat("\nNot passing:\n")
print(res[res$Status != "PASS", c("ID", "Status", "Check")], row.names = FALSE, right = FALSE)
f_out <- file.path(out_dir, "diagnostics_results.csv")
write.csv(res, f_out, row.names = FALSE)
cat("\nFull results written to", f_out, "\n")
