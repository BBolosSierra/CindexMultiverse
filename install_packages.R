### Used during Docker building image

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(warn = 2)

lib_path <- "/home/rstudio/app/r-lib"
dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(lib_path)

BiocManager::install("rhdf5", ask = FALSE, update = TRUE)
BiocManager::install(c("graph", "Rgraphviz"), ask = FALSE, update = TRUE)

install.packages("parallelly", dependencies = TRUE)

install.packages(c(
  "survival", #"MASS", "lattice", "Formula", "SparseM", 
  "data.table",
  "stringr", "dplyr", "tidyr", "gridExtra", "ggalluvial",
  "cowplot", "patchwork"
), lib = lib_path)

install.packages("Hmisc", lib = lib_path)

tryCatch({
  install.packages("rms", lib = lib_path)
}, error = function(e) {
  message("Fallback to rms 6.7-1")
  remotes::install_version("rms", version = "6.7-1", repos = "https://cloud.r-project.org", lib = lib_path)
})

install.packages(c(
  "prodlim", "pec", "riskRegression", "randomForestSRC",
  "doFuture", "future", "progressr", "foreach", "flexsurv",
  "furrr", "gtsummary", "kableExtra", "survivalmodels", "tidyverse"
), lib = lib_path)

#tinytex::install_tinytex(force = TRUE)

BiocManager::install("cBioPortalData", ask = FALSE, update = TRUE)

# Verify it loads correctly
if (!requireNamespace("cBioPortalData", quietly = TRUE)) {
  message("failed to load cBioPortalData after installation.")
  quit(status = 1)
} else {
  message("cBioPortalData installed and loaded successfully.")
}
