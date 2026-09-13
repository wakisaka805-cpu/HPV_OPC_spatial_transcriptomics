cran_packages <- c(
  "dplyr", "ggplot2", "patchwork", "readr", "readxl", "survival", "tidyr", "tidyverse"
)
cran_missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_missing)) install.packages(cran_missing)

bioc_packages <- c("BiocParallel", "GSEABase", "GSVA", "SpatialDecon")
bioc_missing <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(bioc_missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(bioc_missing, ask = FALSE, update = FALSE)
}
