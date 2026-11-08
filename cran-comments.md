## Submission summary

skellam 0.2.5 is a bug-fix release with a change of maintainer team
(Montasser Ghachem remains the maintainer; Oguz Ersan joins as author).

It corrects errors in the distribution and estimation functions, among
them an inverted sign in the likelihood of `skellam.reg()` that swapped
its two parameter blocks (documented as a breaking change in NEWS.md),
wrong results of `pskellam()` and `dskellam(log = TRUE)` for rates
above about 5e4, and `skellam.mle()`/`skellam.reg()` setting
`options(warn = -1)` without restoring it. A testthat suite now covers
all exported functions.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

- local: Ubuntu 24.04, R 4.3.3
- GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)
- R-hub: linux, windows, macos (R-devel)
- win-builder: R-devel

## Downstream dependencies

There are currently no downstream dependencies for this package.
