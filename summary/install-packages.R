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
SNAPSHOT <- "https://packagemanager.posit.co/cran/__linux__/noble/2026-09-01"

# Only what the RELEASE half of the Rmd loads. `vdjdb_paper_figures.Rmd` additionally needs
# `maps`, `scatterpie` and `ggrepel`; the split exists so the release path does not.
PACKAGES <- c(
  "rmarkdown", "knitr", "ggplot2", "RColorBrewer", "data.table", "forcats",
  "ggalluvial", "circlize", "tidyverse", "stringr", "gridExtra", "cowplot"
)

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
