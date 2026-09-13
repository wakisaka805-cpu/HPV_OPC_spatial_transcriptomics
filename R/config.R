# Repository-wide path configuration
# Run scripts from the repository root. Paths may be overridden with environment variables.

REPO_ROOT <- normalizePath(Sys.getenv("HPV_OPC_REPO_ROOT", unset = "."), mustWork = FALSE)
TCGA_ROOT <- normalizePath(Sys.getenv("TCGA_ROOT", unset = file.path(REPO_ROOT, "data", "TCGA")), mustWork = FALSE)
GEOMX_ROOT <- normalizePath(Sys.getenv("GEOMX_ROOT", unset = file.path(REPO_ROOT, "data", "GeoMx")), mustWork = FALSE)
PTT_ROOT <- normalizePath(Sys.getenv("PTT_ROOT", unset = file.path(REPO_ROOT, "data", "PTT")), mustWork = FALSE)
RESULTS_ROOT <- normalizePath(Sys.getenv("RESULTS_ROOT", unset = file.path(REPO_ROOT, "results")), mustWork = FALSE)

dir.create(RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)


# Cross-platform publication font helper.
# On Windows, register Times New Roman under a stable R family name.
if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(TimesNewRoman = grDevices::windowsFont("Times New Roman"))
  PUBLICATION_FONT <- "TimesNewRoman"
} else {
  PUBLICATION_FONT <- "Times New Roman"
}

# Common output directories
FIGURE_ROOT <- file.path(RESULTS_ROOT, "figures")
dir.create(FIGURE_ROOT, recursive = TRUE, showWarnings = FALSE)
