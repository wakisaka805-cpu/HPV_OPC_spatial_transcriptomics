source(file.path("R", "config.R"))

required_files <- c(
  file.path(TCGA_ROOT, "input", "clinical.rds"),
  file.path(TCGA_ROOT, "input", "expression_gene_symbol.rds"),
  file.path(TCGA_ROOT, "TCGA_HNSC_RNAseq_HPV_status_Nulton2017.rds"),
  file.path(REPO_ROOT, "gene_sets", "TCGA_tonsil40_gene_sets_8.rds"),
  file.path(GEOMX_ROOT, "GeoMx_ssGSEA_ROI_level_scores.csv"),
  file.path(GEOMX_ROOT, "GeoMx_ssGSEA_case_level_scores.csv"),
  file.path(GEOMX_ROOT, "GeoMx_LFR_Q3.rds"),
  file.path(GEOMX_ROOT, "GeoMx_ROI_annotation.csv"),
  file.path(PTT_ROOT, "PTT_expression_gene_symbol.rds"),
  file.path(PTT_ROOT, "PTT_ssGSEA_scores_HPV.csv"),
  file.path(REPO_ROOT, "gene_sets", "S1_Table_gene_sets.xlsx")
)

check <- data.frame(
  file = normalizePath(required_files, winslash = "/", mustWork = FALSE),
  exists = file.exists(required_files),
  stringsAsFactors = FALSE
)

print(check, row.names = FALSE)
if (all(check$exists)) {
  message("All required analysis inputs are present.")
} else {
  message("Some required analysis inputs are absent. See DATA_INPUTS.md.")
}
