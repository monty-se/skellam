## --------------------------------------------------------------------------
## revdep_check.R -- check the packages that depend on skellam
##
## Run from the root of the clone before every CRAN submission:
##   Rscript dev/revdep_check.R
##
## It (1) installs the development version of skellam into a temporary
## library, (2) downloads the current CRAN version of each reverse dependency
## and runs R CMD check on it against that library, and (3) compares the
## skellam functions each dependency calls between the CRAN release and the
## development version. CRAN asks submitters to report the outcome in
## cran-comments.md under "Downstream dependencies".
##
## Needs internet access to CRAN and the dependencies' own Imports installed.
##
## Authors: Montasser Ghachem and Oguz Ersan, 2026
## License: GPL (>= 2)
## --------------------------------------------------------------------------

repos <- c(CRAN = "https://cloud.r-project.org")

## Reverse dependencies and the skellam functions they call. Update the list
## with tools::package_dependencies("skellam", reverse = TRUE, which = "all").
revdeps <- list(PINstimation = "qskellam")

found <- tools::package_dependencies("skellam", db = available.packages(repos = repos),
                                     reverse = TRUE, which = "all")[[1]]
if (length(setdiff(found, names(revdeps))))
  message("New reverse dependencies on CRAN, add them to 'revdeps': ",
          paste(setdiff(found, names(revdeps)), collapse = ", "))

work <- file.path(tempdir(), "revdep")
dir.create(work, showWarnings = FALSE)
lib_dev <- file.path(work, "lib-dev"); lib_cran <- file.path(work, "lib-cran")
dir.create(lib_dev, showWarnings = FALSE); dir.create(lib_cran, showWarnings = FALSE)

## 1. Development version and the CRAN release, side by side
install.packages(".", lib = lib_dev, repos = NULL, type = "source", quiet = TRUE)
install.packages("skellam", lib = lib_cran, repos = repos, quiet = TRUE)

## 2. R CMD check of each reverse dependency against the development version
results <- character()
for (pkg in names(revdeps)) {
  tarball <- download.packages(pkg, destdir = work, repos = repos, type = "source")[, 2]
  install.packages(pkg, lib = lib_dev, repos = repos, dependencies = TRUE, quiet = TRUE)
  out <- system2(file.path(R.home("bin"), "R"),
                 c("CMD", "check", "--no-manual", "--as-cran", shQuote(tarball)),
                 env = paste0("R_LIBS=", lib_dev, .Platform$path.sep, .libPaths()[1]),
                 stdout = TRUE, stderr = TRUE)
  status <- grep("^Status:", out, value = TRUE)
  results[pkg] <- if (length(status)) status else "check did not finish"
  cat(sprintf("%-14s %s\n", pkg, results[pkg]))
}

## 3. The functions each dependency calls, CRAN release against development
grid <- expand.grid(p = c(1e-6, 0.005, 0.025, 0.05, 0.5, 0.95, 0.975, 0.995, 1 - 1e-6),
                    l1 = c(0.01, 0.5, 2, 10, 100, 1e3, 1e4),
                    l2 = c(0.01, 0.5, 2, 10, 100, 1e3, 1e4))
call_in <- function(lib, fun) {                       # load one version, call, unload
  f <- get(fun, envir = loadNamespace("skellam", lib.loc = lib))
  on.exit(unloadNamespace("skellam"))
  suppressWarnings(f(grid$p, grid$l1, grid$l2))
}
for (fun in unique(unlist(revdeps))) {
  old <- call_in(lib_cran, fun); new <- call_in(lib_dev, fun)
  same <- identical(old, new) || isTRUE(all.equal(old, new))
  cat(sprintf("%-14s %d inputs, %s\n", fun, nrow(grid),
              if (same) "identical results" else sprintf("%d differ", sum(old != new, na.rm = TRUE))))
  if (!same) print(head(cbind(grid, cran = old, dev = new)[which(old != new), ], 20))
}
