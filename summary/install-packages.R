# 2026-09-26  The R packages the release dashboard needs, pinned to a dated snapshot.
#
#   Rscript summary/install-packages.R
#
# A file rather than a line in the workflow, for two reasons. The list is reviewable here and
# nowhere else -- it is what an `renv.lock` would have given without a lockfile nobody in the
# project has `renv` installed to regenerate. And a multi-line `Rscript -e` in YAML does not work:
# the block scalar keeps the backslash continuations, which reach R as literal backslashes and
# fail with `unexpected end of line`. That is a real build failure, not a hypothetical one.
#
# DATED, not `latest`. A snapshot URL that moves is not a pin, and an unpinned upgrade is exactly
# how `guide = F`, `size =` and `..count..` broke under this document once already.
# The codename is derived, the date is pinned. `noble` was hard-coded, and Posit serves binaries
# per Ubuntu release: the day `ubuntu-latest` rolls to the next LTS, or a self-hosted runner has a
# different codename, P3M silently falls back to source tarballs - and `setup-r` with
# `use-public-rspm: true` installs no -dev headers, so tidyverse fails to compile ten minutes in.
CODENAME <- tryCatch(
  sub('"', "", sub("^VERSION_CODENAME=", "",
                   grep("^VERSION_CODENAME=",
                        suppressWarnings(readLines("/etc/os-release")), value = TRUE)[1]),
      fixed = TRUE),
  error = function(e) NA_character_)
if (is.na(CODENAME) || !nzchar(CODENAME)) CODENAME <- "noble"
SNAPSHOT <- sprintf("https://packagemanager.posit.co/cran/__linux__/%s/2026-09-01", CODENAME)

# What the RELEASE half of the Rmd needs. `vdjdb_paper_figures.Rmd` additionally needs `maps`,
# `scatterpie` and `ggrepel`; the split exists so the release path does not.
#
# `reshape2` is here because of `reshape2::dcast(...)` -- a NAMESPACED call, which appears in no
# `library()` line. Deriving the list from `library()` calls alone missed it and the build died
# fifteen minutes in with "there is no package called 'reshape2'". The check below closes that gap
# for good: it greps both documents for `pkg::` as well as `library(pkg)`.
PACKAGES <- c(
  "rmarkdown", "knitr", "ggplot2", "RColorBrewer", "data.table", "forcats",
  "ggalluvial", "circlize", "tidyverse", "stringr", "gridExtra", "cowplot", "reshape2"
)

# Packages the release document actually references, derived from the document rather than from
# memory. A package named there and missing here is a fifteen-minute failure waiting to happen.
rmd <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "vdjdb_summary.Rmd")
if (file.exists(rmd)) {
  text <- paste(readLines(rmd, warn = FALSE), collapse = "\n")
  referenced <- unique(c(
    regmatches(text, gregexpr("(?<=library\\()[A-Za-z][A-Za-z0-9.]*", text, perl = TRUE))[[1]],
    regmatches(text, gregexpr("[A-Za-z][A-Za-z0-9.]*(?=::)", text, perl = TRUE))[[1]]
  ))
  # `dplyr` arrives inside tidyverse; base namespaces need no install.
  undeclared <- setdiff(referenced, c(PACKAGES, "dplyr", "base", "stats", "utils", "grDevices"))
  if (length(undeclared)) {
    stop("vdjdb_summary.Rmd references packages this script does not install: ",
         paste(undeclared, collapse = ", "))
  }
}

options(repos = c(CRAN = SNAPSHOT), Ncpus = max(1L, parallel::detectCores()))
missing <- setdiff(PACKAGES, rownames(installed.packages()))
if (length(missing)) {
  cat("installing:", paste(missing, collapse = ", "), "\n")
  install.packages(missing)
}
still <- setdiff(PACKAGES, rownames(installed.packages()))
if (length(still)) {
  stop("failed to install: ", paste(still, collapse = ", "))
}
cat("all", length(PACKAGES), "packages present\n")

# Fail here, in seconds, rather than fifteen minutes into the render. The Rmd pins
# `dev.args = list(type = "cairo")`, which needs an R built with cairo; without it every figure
# chunk dies one at a time at the far end of the pipeline. `capabilities()` is free to ask.
if (!capabilities("cairo")) {
  stop("this R has no cairo device; the dashboard pins dev.args = list(type = \"cairo\") ",
       "for parity with the shipped PNGs. Install a cairo-enabled R.")
}
cat("cairo: available\n")
