# ============================================================
# Project:
# TCGA tonsillar OPC transcriptomic analysis
#
# Script:
# 05B_run_ssGSEA.R
#
# Purpose:
# Run ssGSEA for the eight validated gene sets
# in the TCGA tonsil 40-case cohort
# ============================================================


# ------------------------------------------------------------
# 0. Clear workspace
# ------------------------------------------------------------

rm(list = ls())


# ------------------------------------------------------------
# 1. Check required packages
# ------------------------------------------------------------

required_packages <- c(
  "GSVA",
  "BiocParallel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them using BiocManager::install()."
    )
  )
}

library(GSVA)
library(BiocParallel)


# ------------------------------------------------------------
# 2. Define project directories
# ------------------------------------------------------------

workdir <- file.path(
  Sys.getenv("TCGA_ROOT", unset = file.path("data", "TCGA"))
)

expression_file <- file.path(
  workdir,
  "data_expression",
  "processed",
  "TCGA_tonsil_expression_40cases.rds"
)

gene_set_file <- file.path(
  workdir,
  "results",
  "TCGA_tonsil40",
  "gene_sets",
  "TCGA_tonsil40_gene_sets_8.rds"
)

output_dir <- file.path(
  workdir,
  "results",
  "TCGA_tonsil40",
  "ssGSEA"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Confirm input files
# ------------------------------------------------------------

if (!file.exists(expression_file)) {
  stop(
    paste0(
      "Expression file was not found:\n",
      expression_file
    )
  )
}

if (!file.exists(gene_set_file)) {
  stop(
    paste0(
      "Gene-set file was not found:\n",
      gene_set_file
    )
  )
}


# ------------------------------------------------------------
# 4. Load expression and clinical data
# ------------------------------------------------------------

tonsil_40 <- readRDS(expression_file)

expr <- as.matrix(
  tonsil_40$expression
)

clinical <- as.data.frame(
  tonsil_40$clinical,
  stringsAsFactors = FALSE
)

storage.mode(expr) <- "double"


# ------------------------------------------------------------
# 5. Standardize gene symbols
# ------------------------------------------------------------

rownames(expr) <- toupper(
  trimws(
    rownames(expr)
  )
)

if (anyDuplicated(rownames(expr)) > 0) {
  stop(
    "Duplicated gene symbols were found after upper-case conversion."
  )
}


# ------------------------------------------------------------
# 6. Confirm sample alignment
# ------------------------------------------------------------

if (
  !identical(
    colnames(expr),
    clinical$patient_id
  )
) {
  stop(
    paste0(
      "Expression columns and clinical patient IDs ",
      "are not in identical order."
    )
  )
}

cat("Expression dimensions:\n")
print(dim(expr))

cat("\nHPV status:\n")
print(
  table(
    clinical$HPV_status,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 7. Load validated gene sets
# ------------------------------------------------------------

gene_set_object <- readRDS(
  gene_set_file
)

gene_sets <- gene_set_object$gene_sets_matched

if (length(gene_sets) != 8) {
  stop(
    paste0(
      "Expected 8 gene sets, but found ",
      length(gene_sets),
      "."
    )
  )
}

cat("\nGene-set sizes used for ssGSEA:\n")
print(
  vapply(
    gene_sets,
    length,
    integer(1)
  )
)


# ------------------------------------------------------------
# 8. Remove genes with zero variance
# ------------------------------------------------------------

gene_variance <- apply(
  expr,
  1,
  var,
  na.rm = TRUE
)

keep_genes <- is.finite(gene_variance) &
  gene_variance > 0

expr_filtered <- expr[
  keep_genes,
  ,
  drop = FALSE
]

cat("\nGenes before variance filtering:\n")
print(nrow(expr))

cat("\nGenes after variance filtering:\n")
print(nrow(expr_filtered))

if (nrow(expr_filtered) == 0) {
  stop("No genes remained after variance filtering.")
}


# ------------------------------------------------------------
# 9. Reconfirm gene-set overlap after filtering
# ------------------------------------------------------------

gene_sets_filtered <- lapply(
  gene_sets,
  function(genes) {
    intersect(
      genes,
      rownames(expr_filtered)
    )
  }
)

filtered_set_sizes <- vapply(
  gene_sets_filtered,
  length,
  integer(1)
)

cat("\nGene-set sizes after expression filtering:\n")
print(filtered_set_sizes)

if (any(filtered_set_sizes < 5)) {
  stop(
    paste0(
      "At least one gene set has fewer than 5 genes ",
      "after variance filtering."
    )
  )
}


# ------------------------------------------------------------
# 10. Build GSVA ssGSEA parameter object
#
# For continuous normalized RNA expression:
# kcdf = "Gaussian"
# ------------------------------------------------------------

ssgsea_parameter <- ssgseaParam(
  exprData = expr_filtered,
  geneSets = gene_sets_filtered,
  minSize = 5,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE,
  checkNA = "auto",
  use = "everything"
)


# ------------------------------------------------------------
# 11. Run ssGSEA
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("Running ssGSEA...\n")
cat("============================================\n")

ssgsea_scores <- gsva(
  ssgsea_parameter,
  verbose = TRUE,
  BPPARAM = SerialParam()
)

ssgsea_scores <- as.matrix(
  ssgsea_scores
)


# ------------------------------------------------------------
# 12. Validate ssGSEA result
# ------------------------------------------------------------

cat("\nssGSEA score dimensions:\n")
print(dim(ssgsea_scores))

if (nrow(ssgsea_scores) != 8) {
  stop(
    paste0(
      "Expected 8 ssGSEA score rows, but obtained ",
      nrow(ssgsea_scores),
      "."
    )
  )
}

if (
  !identical(
    colnames(ssgsea_scores),
    clinical$patient_id
  )
) {
  stop(
    "ssGSEA sample order does not match clinical data."
  )
}

if (anyNA(ssgsea_scores)) {
  stop("Missing values were detected in ssGSEA scores.")
}

if (any(!is.finite(ssgsea_scores))) {
  stop("Non-finite values were detected in ssGSEA scores.")
}


# ------------------------------------------------------------
# 13. Create wide-format output
# ------------------------------------------------------------

scores_wide <- data.frame(
  patient_id = colnames(ssgsea_scores),
  t(ssgsea_scores),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical_selected <- clinical[
  ,
  c(
    "patient_id",
    "HPV_status",
    "HPV_type",
    "N_group"
  ),
  drop = FALSE
]

scores_wide <- merge(
  clinical_selected,
  scores_wide,
  by = "patient_id",
  sort = FALSE
)

scores_wide <- scores_wide[
  match(
    clinical$patient_id,
    scores_wide$patient_id
  ),
  ,
  drop = FALSE
]

if (
  !identical(
    scores_wide$patient_id,
    clinical$patient_id
  )
) {
  stop(
    "Wide-format score table lost the original sample order."
  )
}


# ------------------------------------------------------------
# 14. Create long-format output
# ------------------------------------------------------------

scores_long <- do.call(
  rbind,
  lapply(
    rownames(ssgsea_scores),
    function(set_name) {
      
      data.frame(
        patient_id = colnames(ssgsea_scores),
        HPV_status = clinical$HPV_status,
        HPV_type = clinical$HPV_type,
        N_group = clinical$N_group,
        Signature = set_name,
        Score = as.numeric(
          ssgsea_scores[
            set_name,
            colnames(ssgsea_scores)
          ]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(scores_long) <- NULL


# ------------------------------------------------------------
# 15. Add publication labels
# ------------------------------------------------------------

signature_label_map <- c(
  "STEM_LIKE_T_CELL_PROGRAM" =
    "Stem-like",
  "PRE_EXHAUSTION_PROGRAM" =
    "Tpex",
  "T_CELL_ACTIVITY" =
    "T-cell activity",
  "B_CELL_ACTIVITY" =
    "B-cell activity",
  "TEX_LIKE_PROGRAM" =
    "Tex-like",
  "HALLMARK_ESTROGEN_RESPONSE_EARLY" =
    "Estrogen early",
  "HALLMARK_ESTROGEN_RESPONSE_LATE" =
    "Estrogen late",
  "HALLMARK_BILE_ACID_METABOLISM" =
    "Bile acid"
)

scores_long$Signature_label <- unname(
  signature_label_map[
    scores_long$Signature
  ]
)

if (anyNA(scores_long$Signature_label)) {
  stop(
    "At least one publication label could not be assigned."
  )
}

scores_long$HPV_group <- ifelse(
  scores_long$HPV_status == "HPV_positive",
  "HPV-positive",
  ifelse(
    scores_long$HPV_status == "HPV_negative",
    "HPV-negative",
    NA_character_
  )
)

if (anyNA(scores_long$HPV_group)) {
  stop(
    "Unexpected HPV_status values were detected."
  )
}

scores_long$HPV_group <- factor(
  scores_long$HPV_group,
  levels = c(
    "HPV-negative",
    "HPV-positive"
  )
)

scores_long$Signature_label <- factor(
  scores_long$Signature_label,
  levels = c(
    "Stem-like",
    "Tpex",
    "T-cell activity",
    "B-cell activity",
    "Tex-like",
    "Estrogen early",
    "Estrogen late",
    "Bile acid"
  )
)


# ------------------------------------------------------------
# 16. Create output object
# ------------------------------------------------------------

ssgsea_object <- list(
  scores_matrix = ssgsea_scores,
  scores_wide = scores_wide,
  scores_long = scores_long,
  clinical = clinical,
  expression_dimensions_original = dim(expr),
  expression_dimensions_filtered = dim(expr_filtered),
  gene_set_sizes_before_filtering =
    vapply(
      gene_sets,
      length,
      integer(1)
    ),
  gene_set_sizes_after_filtering =
    filtered_set_sizes,
  gsva_version =
    as.character(
      packageVersion("GSVA")
    ),
  created_at = Sys.time()
)


# ------------------------------------------------------------
# 17. Save outputs
# ------------------------------------------------------------

rds_output <- file.path(
  output_dir,
  "TCGA_tonsil40_ssGSEA_8signatures.rds"
)

matrix_output <- file.path(
  output_dir,
  "TCGA_tonsil40_ssGSEA_score_matrix.csv"
)

wide_output <- file.path(
  output_dir,
  "TCGA_tonsil40_ssGSEA_scores_wide.csv"
)

long_output <- file.path(
  output_dir,
  "TCGA_tonsil40_ssGSEA_scores_long.csv"
)

saveRDS(
  ssgsea_object,
  rds_output
)

write.csv(
  data.frame(
    Signature = rownames(ssgsea_scores),
    ssgsea_scores,
    check.names = FALSE
  ),
  matrix_output,
  row.names = FALSE
)

write.csv(
  scores_wide,
  wide_output,
  row.names = FALSE
)

write.csv(
  scores_long,
  long_output,
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Reload and validate saved object
# ------------------------------------------------------------

saved_object <- readRDS(
  rds_output
)

if (
  !identical(
    dim(saved_object$scores_matrix),
    c(8L, 40L)
  )
) {
  stop(
    "Saved ssGSEA object failed dimension validation."
  )
}


# ------------------------------------------------------------
# 19. Completion message
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("TCGA tonsil 40-case ssGSEA completed.\n")
cat("Signatures: 8\n")
cat("Samples: 40\n")
cat("\nOutput directory:\n")
cat(output_dir, "\n")
cat("\nMain RDS file:\n")
cat(rds_output, "\n")
cat("\nReady for 05C_prepare_Figure2.R\n")
cat("============================================\n")