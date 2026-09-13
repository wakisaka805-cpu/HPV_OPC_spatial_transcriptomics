source(file.path("R", "config.R"))

# ============================================================
# 15_TCGA_OPSCC_Sensitivity_Analysis.R
#
# TCGA-HNSC anatomically restricted OPSCC sensitivity analysis
# ============================================================



# Public-release convenience: if clinical/expr are not already in memory,
# load them from data/TCGA/input. These files are not distributed here.
if (!exists("clinical")) {
  f <- file.path(TCGA_ROOT, "input", "clinical.rds")
  if (file.exists(f)) clinical <- readRDS(f)
}
if (!exists("expr")) {
  f <- file.path(TCGA_ROOT, "input", "expression_gene_symbol.rds")
  if (file.exists(f)) expr <- readRDS(f)
}

# ------------------------------------------------------------
# 0. Basic objects
# ------------------------------------------------------------

# Main objects are assumed to be available in the current project:
#   clinical : TCGA-HNSC clinical data
#   expr     : gene expression matrix (genes x samples)

stopifnot(exists("clinical"))
stopifnot(exists("expr"))


# ------------------------------------------------------------
# 1. Define strict OPSCC cohort by ICD-10
# ------------------------------------------------------------

opscc_icd <- c(
  "C09.9",  # Tonsil, NOS
  "C01",    # Base of tongue, NOS
  "C10.9",  # Oropharynx, NOS
  "C10.3"   # Posterior wall of oropharynx
)

opscc_clin <- clinical[
  clinical$icd_10_code %in% opscc_icd,
  ,
  drop = FALSE
]

# Convert TCGA clinical IDs from hyphen to period format
opscc_clin$expr_id <- gsub(
  "-",
  ".",
  opscc_clin$bcr_patient_barcode
)

# RNA-seq availability
opscc_clin$RNAseq_available <-
  opscc_clin$expr_id %in% colnames(expr)

cat("\n===== OPSCC clinical candidates =====\n")
cat("Clinical candidates:", nrow(opscc_clin), "\n")

cat("\n===== RNA-seq availability =====\n")
print(table(opscc_clin$RNAseq_available))

# Keep only cases with RNA-seq
opscc_clin_expr <- opscc_clin[
  opscc_clin$RNAseq_available,
  ,
  drop = FALSE
]

cat("\nFinal OPSCC cases with RNA-seq:",
    nrow(opscc_clin_expr), "\n")

cat("\n===== Anatomical site =====\n")
print(
  table(
    opscc_clin_expr$tissue_or_organ_of_origin
  )
)


# ------------------------------------------------------------
# 2. Load Nulton2017 HPV annotation
# ------------------------------------------------------------

old_outdir <- TCGA_ROOT

hpv_file <- file.path(
  old_outdir,
  "TCGA_HNSC_RNAseq_HPV_status_Nulton2017.rds"
)

stopifnot(file.exists(hpv_file))

hpv_ann <- readRDS(hpv_file)

cat("\n===== Full TCGA-HNSC HPV status =====\n")
print(
  table(
    hpv_ann$HPV_status,
    useNA = "ifany"
  )
)

cat("\n===== Full TCGA-HNSC HPV type =====\n")
print(
  sort(
    table(
      hpv_ann$HPV_type,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)


# ------------------------------------------------------------
# 3. Merge HPV annotation
# ------------------------------------------------------------

opscc_clin_expr <- merge(
  opscc_clin_expr,
  hpv_ann[, c(
    "sample_name",
    "submitter_id",
    "HPV_status",
    "HPV_type"
  )],
  by.x = "expr_id",
  by.y = "sample_name",
  all.x = TRUE
)

cat("\n===== OPSCC HPV status =====\n")
print(
  table(
    opscc_clin_expr$HPV_status,
    useNA = "ifany"
  )
)

cat("\n===== OPSCC HPV type =====\n")
print(
  sort(
    table(
      opscc_clin_expr$HPV_type,
      useNA = "ifany"
    ),
    decreasing = TRUE
  )
)

cat("\n===== OPSCC HPV status by anatomical site =====\n")
print(
  table(
    opscc_clin_expr$tissue_or_organ_of_origin,
    opscc_clin_expr$HPV_status,
    useNA = "ifany"
  )
)

stopifnot(nrow(opscc_clin_expr) == 77)
stopifnot(length(unique(opscc_clin_expr$expr_id)) == 77)


# ------------------------------------------------------------
# 4. Prepare clinical covariates
# ------------------------------------------------------------

opscc_clin_expr$age_years <-
  opscc_clin_expr$age_at_diagnosis / 365.25

opscc_clin_expr$smoking_group <- ifelse(
  opscc_clin_expr$tobacco_smoking_status ==
    "Lifelong Non-Smoker",
  "Never",
  ifelse(
    opscc_clin_expr$tobacco_smoking_status ==
      "Unknown",
    NA,
    "Ever"
  )
)

opscc_clin_expr$clinical_stage_group <- ifelse(
  opscc_clin_expr$ajcc_clinical_stage %in%
    c("Stage I", "Stage II"),
  "I-II",
  ifelse(
    opscc_clin_expr$ajcc_clinical_stage %in%
      c(
        "Stage III",
        "Stage IVA",
        "Stage IVB",
        "Stage IVC"
      ),
    "III-IV",
    NA
  )
)

opscc_clin_expr$subsite_group <- ifelse(
  opscc_clin_expr$tissue_or_organ_of_origin ==
    "Tonsil, NOS",
  "Tonsil",
  ifelse(
    opscc_clin_expr$tissue_or_organ_of_origin ==
      "Base of tongue, NOS",
    "Base_of_tongue",
    "Other_oropharynx"
  )
)

# Factor settings
opscc_clin_expr$HPV_status <- factor(
  opscc_clin_expr$HPV_status,
  levels = c(
    "HPV_negative",
    "HPV_positive"
  )
)

opscc_clin_expr$gender <- factor(
  opscc_clin_expr$gender
)

opscc_clin_expr$smoking_group <- factor(
  opscc_clin_expr$smoking_group,
  levels = c(
    "Never",
    "Ever"
  )
)

opscc_clin_expr$clinical_stage_group <- factor(
  opscc_clin_expr$clinical_stage_group,
  levels = c(
    "I-II",
    "III-IV"
  )
)

opscc_clin_expr$subsite_group <- factor(
  opscc_clin_expr$subsite_group,
  levels = c(
    "Other_oropharynx",
    "Base_of_tongue",
    "Tonsil"
  )
)


# ------------------------------------------------------------
# 5. Clinical covariate summary
# ------------------------------------------------------------

cat("\n===== Gender =====\n")
print(
  table(
    opscc_clin_expr$gender,
    opscc_clin_expr$HPV_status
  )
)

cat("\n===== Smoking =====\n")
print(
  table(
    opscc_clin_expr$smoking_group,
    opscc_clin_expr$HPV_status,
    useNA = "ifany"
  )
)

cat("\n===== Clinical stage =====\n")
print(
  table(
    opscc_clin_expr$clinical_stage_group,
    opscc_clin_expr$HPV_status,
    useNA = "ifany"
  )
)

cat("\n===== Subsite =====\n")
print(
  table(
    opscc_clin_expr$subsite_group,
    opscc_clin_expr$HPV_status
  )
)


# ------------------------------------------------------------
# 6. Complete cases for adjusted analysis
# ------------------------------------------------------------

model_vars <- c(
  "HPV_status",
  "age_years",
  "gender",
  "smoking_group",
  "clinical_stage_group",
  "subsite_group"
)

cc <- complete.cases(
  opscc_clin_expr[, model_vars]
)

cat("\n===== Complete cases =====\n")
cat(
  "Total OPSCC:",
  nrow(opscc_clin_expr),
  "\n"
)

cat(
  "Complete cases:",
  sum(cc),
  "\n"
)

cat(
  "Excluded:",
  sum(!cc),
  "\n"
)

cat("\n===== HPV status among complete cases =====\n")
print(
  table(
    opscc_clin_expr$HPV_status[cc]
  )
)


# ------------------------------------------------------------
# 7. Create OPSCC expression matrix
# ------------------------------------------------------------

opscc_expr <- expr[
  ,
  opscc_clin_expr$expr_id,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(opscc_expr),
    opscc_clin_expr$expr_id
  )
)

cat("\n===== OPSCC expression matrix =====\n")
print(dim(opscc_expr))

cat(
  "\nAll IDs identical:",
  identical(
    colnames(opscc_expr),
    opscc_clin_expr$expr_id
  ),
  "\n"
)

cat("\n===== Final HPV distribution =====\n")
print(
  table(
    opscc_clin_expr$HPV_status
  )
)

# ------------------------------------------------------------
# 8. Load the exact 8 gene sets used in the previous analysis
# ------------------------------------------------------------

gene_set_file <- file.path(
  REPO_ROOT,
  "gene_sets",
  "TCGA_tonsil40_gene_sets_8.rds"
)

gene_set_archive <- readRDS(gene_set_file)

# Use the original/full definitions
gene_sets_full <- gene_set_archive$gene_sets_full

stopifnot(length(gene_sets_full) == 8)

cat("\n===== Eight gene sets =====\n")
print(names(gene_sets_full))

cat("\n===== Original gene-set sizes =====\n")
print(sapply(gene_sets_full, length))


# ------------------------------------------------------------
# 9. Match gene sets to OPSCC expression matrix
# ------------------------------------------------------------

gene_sets_opscc <- lapply(
  gene_sets_full,
  function(x) {
    intersect(x, rownames(opscc_expr))
  }
)

gene_sets_missing_opscc <- lapply(
  gene_sets_full,
  function(x) {
    setdiff(x, rownames(opscc_expr))
  }
)

cat("\n===== OPSCC matched gene-set sizes =====\n")
print(
  sapply(
    gene_sets_opscc,
    length
  )
)

cat("\n===== Genes missing from OPSCC expression =====\n")

for (nm in names(gene_sets_missing_opscc)) {
  
  cat("\n---", nm, "---\n")
  
  missing <- gene_sets_missing_opscc[[nm]]
  
  if (length(missing) == 0) {
    cat("None\n")
  } else {
    print(missing)
  }
}


# ------------------------------------------------------------
# 10. Compare with previous tonsil40 matching
# ------------------------------------------------------------

cat("\n===== Same matching as tonsil40? =====\n")

for (nm in names(gene_sets_opscc)) {
  
  same <- identical(
    sort(gene_sets_opscc[[nm]]),
    sort(gene_set_archive$gene_sets_matched[[nm]])
  )
  
  cat(
    nm,
    ":",
    same,
    "\n"
  )
}

# ------------------------------------------------------------
# 11. Run ssGSEA in the 77-case OPSCC cohort
# ------------------------------------------------------------

stopifnot(
  ncol(opscc_expr) == 77,
  length(gene_sets_opscc) == 8
)

# GSVA 2.x / current syntax
ssgsea_param <- GSVA::ssgseaParam(
  exprData = opscc_expr,
  geneSets = gene_sets_opscc,
  normalize = TRUE
)

opscc_ssgsea <- GSVA::gsva(
  ssgsea_param,
  verbose = FALSE
)

cat("\n===== OPSCC ssGSEA matrix =====\n")
print(dim(opscc_ssgsea))

cat("\n===== Signature names =====\n")
print(rownames(opscc_ssgsea))

cat("\n===== Sample alignment =====\n")
print(
  identical(
    colnames(opscc_ssgsea),
    opscc_clin_expr$expr_id
  )
)

cat("\n===== Score summary =====\n")
print(
  t(
    apply(
      opscc_ssgsea,
      1,
      summary
    )
  )
)

cat("\n===== Any NA? =====\n")
print(anyNA(opscc_ssgsea))

cat("\n===== Any infinite values? =====\n")
print(any(!is.finite(opscc_ssgsea)))

# ------------------------------------------------------------
# 12. Unadjusted HPV-positive vs HPV-negative comparison
# ------------------------------------------------------------

# Convert ssGSEA matrix to sample-level data frame
score_df <- as.data.frame(
  t(opscc_ssgsea),
  check.names = FALSE
)

score_df$expr_id <- rownames(score_df)

# Confirm alignment
stopifnot(
  identical(
    score_df$expr_id,
    opscc_clin_expr$expr_id
  )
)

# Add HPV status
score_df$HPV_status <- opscc_clin_expr$HPV_status


# Cliff's delta:
# positive value = higher score in HPV-positive tumors
cliffs_delta <- function(x_pos, x_neg) {
  
  x_pos <- x_pos[is.finite(x_pos)]
  x_neg <- x_neg[is.finite(x_neg)]
  
  comp <- outer(
    x_pos,
    x_neg,
    "-"
  )
  
  mean(sign(comp))
}


# Run analysis for all 8 signatures
unadjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      x_pos <- score_df[
        score_df$HPV_status == "HPV_positive",
        sig
      ]
      
      x_neg <- score_df[
        score_df$HPV_status == "HPV_negative",
        sig
      ]
      
      data.frame(
        Signature = sig,
        
        n_HPV_negative = length(x_neg),
        n_HPV_positive = length(x_pos),
        
        median_HPV_negative = median(
          x_neg,
          na.rm = TRUE
        ),
        
        median_HPV_positive = median(
          x_pos,
          na.rm = TRUE
        ),
        
        median_difference = median(
          x_pos,
          na.rm = TRUE
        ) - median(
          x_neg,
          na.rm = TRUE
        ),
        
        cliffs_delta = cliffs_delta(
          x_pos,
          x_neg
        ),
        
        raw_p = wilcox.test(
          x_pos,
          x_neg,
          exact = FALSE
        )$p.value,
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# BH correction across the 8 prespecified signatures
unadjusted_results$BH_FDR <- p.adjust(
  unadjusted_results$raw_p,
  method = "BH"
)

# Sort by signature order used in the analysis
unadjusted_results <- unadjusted_results[
  match(
    rownames(opscc_ssgsea),
    unadjusted_results$Signature
  ),
]

rownames(unadjusted_results) <- NULL


cat(
  "\n===== OPSCC HPV unadjusted results =====\n"
)

print(
  unadjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Multivariable-adjusted analysis
# ------------------------------------------------------------

# Add clinical covariates to score data
analysis_df <- cbind(
  opscc_clin_expr,
  score_df[, rownames(opscc_ssgsea), drop = FALSE]
)

# Complete cases for the prespecified clinical covariates
analysis_cc <- analysis_df[
  complete.cases(
    analysis_df[, model_vars]
  ),
  ,
  drop = FALSE
]

cat("\n===== Adjusted-analysis cohort =====\n")
cat("N =", nrow(analysis_cc), "\n")
print(table(analysis_cc$HPV_status))


# Fit the same prespecified model for each signature
adjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      
      dat$Score <- dat[[sig]]
      
      fit <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      coef_name <- "HPV_statusHPV_positive"
      
      beta <- coef(summary(fit))[
        coef_name,
        "Estimate"
      ]
      
      se <- coef(summary(fit))[
        coef_name,
        "Std. Error"
      ]
      
      p <- coef(summary(fit))[
        coef_name,
        "Pr(>|t|)"
      ]
      
      ci <- confint(
        fit,
        parm = coef_name,
        level = 0.95
      )
      
      data.frame(
        Signature = sig,
        N = nobs(fit),
        adjusted_beta = beta,
        SE = se,
        CI95_low = ci[1],
        CI95_high = ci[2],
        raw_p = p,
        stringsAsFactors = FALSE
      )
    }
  )
)


# BH correction across the 8 prespecified signatures
adjusted_results$BH_FDR <- p.adjust(
  adjusted_results$raw_p,
  method = "BH"
)

rownames(adjusted_results) <- NULL


cat(
  "\n===== OPSCC multivariable-adjusted HPV effects =====\n"
)

print(
  adjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Basic diagnostics for adjusted linear models
# ------------------------------------------------------------

diagnostic_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      dat$Score <- dat[[sig]]
      
      fit <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      cd <- cooks.distance(fit)
      
      data.frame(
        Signature = sig,
        N = nobs(fit),
        R_squared = summary(fit)$r.squared,
        adjusted_R_squared = summary(fit)$adj.r.squared,
        max_Cooks_D = max(cd),
        n_Cooks_D_gt_4_over_N =
          sum(cd > 4 / nobs(fit)),
        stringsAsFactors = FALSE
      )
    }
  )
)

cat(
  "\n===== Adjusted model diagnostics =====\n"
)

print(
  diagnostic_results,
  digits = 4,
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Influence-point sensitivity analysis
# ------------------------------------------------------------

influence_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      dat$Score <- dat[[sig]]
      
      # Primary adjusted model
      fit_full <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      # Identify observations with Cook's D > 4/N
      cd <- cooks.distance(fit_full)
      
      keep <- cd <= 4 / nobs(fit_full)
      
      # Sensitivity model after excluding influential observations
      fit_sens <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat[keep, , drop = FALSE]
      )
      
      coef_name <- "HPV_statusHPV_positive"
      
      sm <- coef(summary(fit_sens))
      
      ci <- confint(
        fit_sens,
        parm = coef_name,
        level = 0.95
      )
      
      data.frame(
        Signature = sig,
        N_full = nobs(fit_full),
        N_excluded = sum(!keep),
        N_sensitivity = nobs(fit_sens),
        
        HPV_beta_sensitivity =
          sm[coef_name, "Estimate"],
        
        CI95_low =
          ci[1],
        
        CI95_high =
          ci[2],
        
        raw_p =
          sm[coef_name, "Pr(>|t|)"],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

# BH correction across the same 8 signatures
influence_results$BH_FDR <- p.adjust(
  influence_results$raw_p,
  method = "BH"
)

cat(
  "\n===== Influence-point sensitivity analysis =====\n"
)

print(
  influence_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Save TCGA-OPSCC sensitivity-analysis results
# ------------------------------------------------------------

opscc_outdir <- file.path(RESULTS_ROOT, "TCGA_OPSCC77_revision")

dir.create(
  opscc_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Clinical cohort
write.csv(
  opscc_clin_expr,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_clinical_HPV_annotation.csv"
  ),
  row.names = FALSE
)

# ssGSEA score matrix
write.csv(
  opscc_ssgsea,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_ssGSEA_score_matrix.csv"
  )
)

saveRDS(
  opscc_ssgsea,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_ssGSEA_score_matrix.rds"
  )
)

# Unadjusted results
write.csv(
  unadjusted_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_unadjusted_HPV_results.csv"
  ),
  row.names = FALSE
)

# Primary multivariable-adjusted results
write.csv(
  adjusted_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC74_adjusted_HPV_results.csv"
  ),
  row.names = FALSE
)

# Model diagnostics
write.csv(
  diagnostic_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC74_model_diagnostics.csv"
  ),
  row.names = FALSE
)

# Influence-point sensitivity analysis
write.csv(
  influence_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC_influence_sensitivity_results.csv"
  ),
  row.names = FALSE
)

# Save the exact matched gene sets used
saveRDS(
  gene_sets_opscc,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_gene_sets_used.rds"
  )
)

# Save a compact reproducibility object
saveRDS(
  list(
    clinical = opscc_clin_expr,
    gene_sets = gene_sets_opscc,
    ssGSEA = opscc_ssgsea,
    unadjusted = unadjusted_results,
    adjusted = adjusted_results,
    diagnostics = diagnostic_results,
    influence_sensitivity = influence_results
  ),
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_revision_analysis_complete.rds"
  )
)

cat("\n===== Results saved to =====\n")
cat(opscc_outdir, "\n")

cat("\n===== Saved files =====\n")
print(list.files(opscc_outdir))

# ------------------------------------------------------------
# 17. Complementary stem-like marker-gene analysis
#     in the 77-case OPSCC cohort
# ------------------------------------------------------------

stem_marker_genes <- c(
  "TCF7",
  "IL7R",
  "CCR7",
  "LEF1",
  "BACH2",
  "SELL",
  "SLAMF6",
  "CD28"
)

cat("\n===== Stem-like marker genes present in expression matrix =====\n")

marker_presence <- data.frame(
  Gene = stem_marker_genes,
  Present = stem_marker_genes %in% rownames(opscc_expr)
)

print(marker_presence, row.names = FALSE)

stopifnot(
  all(stem_marker_genes %in% rownames(opscc_expr))
)


# Unadjusted HPV-positive vs HPV-negative comparison
marker_unadjusted_results <- do.call(
  rbind,
  lapply(
    stem_marker_genes,
    function(gene) {
      
      x_pos <- opscc_expr[
        gene,
        opscc_clin_expr$HPV_status == "HPV_positive"
      ]
      
      x_neg <- opscc_expr[
        gene,
        opscc_clin_expr$HPV_status == "HPV_negative"
      ]
      
      data.frame(
        Gene = gene,
        
        n_HPV_negative = length(x_neg),
        n_HPV_positive = length(x_pos),
        
        median_HPV_negative =
          median(x_neg, na.rm = TRUE),
        
        median_HPV_positive =
          median(x_pos, na.rm = TRUE),
        
        median_difference =
          median(x_pos, na.rm = TRUE) -
          median(x_neg, na.rm = TRUE),
        
        cliffs_delta =
          cliffs_delta(
            x_pos,
            x_neg
          ),
        
        raw_p =
          wilcox.test(
            x_pos,
            x_neg,
            exact = FALSE
          )$p.value,
        
        stringsAsFactors = FALSE
      )
    }
  )
)


# BH correction across the 8 prespecified marker genes
marker_unadjusted_results$BH_FDR <- p.adjust(
  marker_unadjusted_results$raw_p,
  method = "BH"
)


cat(
  "\n===== OPSCC complementary marker-gene analysis =====\n"
)

print(
  marker_unadjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. Multivariable-adjusted marker-gene analysis
# ------------------------------------------------------------

marker_adjusted_results <- do.call(
  rbind,
  lapply(
    stem_marker_genes,
    function(gene) {
      
      dat <- opscc_clin_expr[cc, , drop = FALSE]
      
      dat$Expression <- as.numeric(
        opscc_expr[
          gene,
          dat$expr_id
        ]
      )
      
      fit <- lm(
        Expression ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      coef_name <- "HPV_statusHPV_positive"
      
      sm <- coef(summary(fit))
      
      ci <- confint(
        fit,
        parm = coef_name,
        level = 0.95
      )
      
      data.frame(
        Gene = gene,
        N = nobs(fit),
        
        adjusted_beta =
          sm[coef_name, "Estimate"],
        
        SE =
          sm[coef_name, "Std. Error"],
        
        CI95_low =
          ci[1],
        
        CI95_high =
          ci[2],
        
        raw_p =
          sm[coef_name, "Pr(>|t|)"],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

marker_adjusted_results$BH_FDR <- p.adjust(
  marker_adjusted_results$raw_p,
  method = "BH"
)

cat(
  "\n===== OPSCC adjusted marker-gene analysis =====\n"
)

print(
  marker_adjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 19. Check overall survival variables and event counts
# ------------------------------------------------------------

cat("\n===== Candidate survival-related columns =====\n")

surv_cols <- grep(
  "surv|death|dead|vital|follow|days_to|time",
  colnames(opscc_clin_expr),
  ignore.case = TRUE,
  value = TRUE
)

print(surv_cols)


cat("\n===== Vital status =====\n")

if ("vital_status" %in% colnames(opscc_clin_expr)) {
  print(
    table(
      opscc_clin_expr$vital_status,
      useNA = "ifany"
    )
  )
}


cat("\n===== Days to death =====\n")

if ("days_to_death" %in% colnames(opscc_clin_expr)) {
  print(
    summary(
      opscc_clin_expr$days_to_death
    )
  )
}


cat("\n===== Days to last follow-up =====\n")

if ("days_to_last_follow_up" %in% colnames(opscc_clin_expr)) {
  print(
    summary(
      opscc_clin_expr$days_to_last_follow_up
    )
  )
}

# ------------------------------------------------------------
# 20. Construct overall-survival dataset
# ------------------------------------------------------------

# OS event:
# Dead = 1, Alive = 0
opscc_clin_expr$OS_event <- ifelse(
  opscc_clin_expr$vital_status == "Dead",
  1,
  ifelse(
    opscc_clin_expr$vital_status == "Alive",
    0,
    NA
  )
)

# OS time:
# Dead  -> days_to_death
# Alive -> days_to_last_follow_up
opscc_clin_expr$OS_days <- ifelse(
  opscc_clin_expr$OS_event == 1,
  opscc_clin_expr$days_to_death,
  opscc_clin_expr$days_to_last_follow_up
)

# Convert explicitly to numeric
opscc_clin_expr$OS_days <- as.numeric(
  opscc_clin_expr$OS_days
)

cat("\n===== Overall survival data =====\n")

cat(
  "Total OPSCC:",
  nrow(opscc_clin_expr),
  "\n"
)

cat(
  "Usable OS cases:",
  sum(
    !is.na(opscc_clin_expr$OS_days) &
      !is.na(opscc_clin_expr$OS_event)
  ),
  "\n"
)

cat(
  "Deaths:",
  sum(
    opscc_clin_expr$OS_event == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Censored:",
  sum(
    opscc_clin_expr$OS_event == 0,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n===== OS time summary (days) =====\n")

print(
  summary(
    opscc_clin_expr$OS_days
  )
)

cat("\n===== Death events by HPV status =====\n")

print(
  table(
    opscc_clin_expr$HPV_status,
    opscc_clin_expr$OS_event,
    useNA = "ifany"
  )
)

cat("\n===== Check invalid OS times =====\n")

print(
  opscc_clin_expr[
    is.na(opscc_clin_expr$OS_days) |
      opscc_clin_expr$OS_days <= 0,
    c(
      "expr_id",
      "HPV_status",
      "vital_status",
      "days_to_death",
      "days_to_last_follow_up",
      "OS_days",
      "OS_event"
    )
  ]
)

# ------------------------------------------------------------
# 21. Univariate Cox analysis for the 8 signatures
# ------------------------------------------------------------

stopifnot(
  requireNamespace("survival", quietly = TRUE)
)

# Add OS variables to score data
cox_df <- cbind(
  opscc_clin_expr[
    ,
    c(
      "expr_id",
      "OS_days",
      "OS_event",
      "HPV_status"
    )
  ],
  score_df[
    ,
    rownames(opscc_ssgsea),
    drop = FALSE
  ]
)

stopifnot(
  identical(
    cox_df$expr_id,
    score_df$expr_id
  )
)

cox_univ_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- cox_df[
        complete.cases(
          cox_df[
            ,
            c(
              "OS_days",
              "OS_event",
              sig
            )
          ]
        ),
        ,
        drop = FALSE
      ]
      
      dat$Score <- dat[[sig]]
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~ Score,
        data = dat
      )
      
      sm <- summary(fit)
      
      beta <- sm$coefficients[
        "Score",
        "coef"
      ]
      
      HR <- sm$coefficients[
        "Score",
        "exp(coef)"
      ]
      
      p <- sm$coefficients[
        "Score",
        "Pr(>|z|)"
      ]
      
      ci_low <- sm$conf.int[
        "Score",
        "lower .95"
      ]
      
      ci_high <- sm$conf.int[
        "Score",
        "upper .95"
      ]
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        beta = beta,
        HR = HR,
        CI95_low = ci_low,
        CI95_high = ci_high,
        raw_p = p,
        stringsAsFactors = FALSE
      )
    }
  )
)

# BH correction across the 8 prespecified signatures
cox_univ_results$BH_FDR <- p.adjust(
  cox_univ_results$raw_p,
  method = "BH"
)

cat(
  "\n===== OPSCC univariate Cox analysis =====\n"
)

print(
  cox_univ_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 22. Cox analysis per 1-SD increase in signature score
# ------------------------------------------------------------

cox_sd_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- cox_df[
        complete.cases(
          cox_df[, c("OS_days", "OS_event", sig)]
        ),
        ,
        drop = FALSE
      ]
      
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~ Score_z,
        data = dat
      )
      
      sm <- summary(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cox_sd_results$BH_FDR <- p.adjust(
  cox_sd_results$raw_p,
  method = "BH"
)

cat(
  "\n===== OPSCC Cox analysis per 1-SD increase =====\n"
)

print(
  cox_sd_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 23. Cox analysis adjusted for HPV status
#     Signature score standardized per 1 SD
# ------------------------------------------------------------

cox_hpv_adjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- cox_df[
        complete.cases(
          cox_df[
            ,
            c(
              "OS_days",
              "OS_event",
              "HPV_status",
              sig
            )
          ]
        ),
        ,
        drop = FALSE
      ]
      
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      # Ensure HPV-negative is the reference
      dat$HPV_status <- factor(
        dat$HPV_status,
        levels = c(
          "HPV_negative",
          "HPV_positive"
        )
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~
          Score_z +
          HPV_status,
        data = dat
      )
      
      sm <- summary(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        HPV_HR =
          sm$coefficients[
            "HPV_statusHPV_positive",
            "exp(coef)"
          ],
        
        HPV_p =
          sm$coefficients[
            "HPV_statusHPV_positive",
            "Pr(>|z|)"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

# BH correction across the 8 signature coefficients
cox_hpv_adjusted_results$BH_FDR <- p.adjust(
  cox_hpv_adjusted_results$raw_p,
  method = "BH"
)

cat(
  "\n===== OPSCC Cox analysis adjusted for HPV status =====\n"
)

print(
  cox_hpv_adjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 24. Save marker-gene and survival analyses
# ------------------------------------------------------------

write.csv(
  marker_unadjusted_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_marker_unadjusted_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  marker_adjusted_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC74_marker_adjusted_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_univ_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_Cox_univariate_rawscale.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_sd_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_Cox_univariate_perSD.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_hpv_adjusted_results,
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_perSD.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    marker_unadjusted = marker_unadjusted_results,
    marker_adjusted = marker_adjusted_results,
    cox_univariate_rawscale = cox_univ_results,
    cox_univariate_perSD = cox_sd_results,
    cox_HPV_adjusted_perSD = cox_hpv_adjusted_results
  ),
  file.path(
    opscc_outdir,
    "TCGA_OPSCC77_marker_survival_results.rds"
  )
)

cat("\n===== Updated OPSCC output files =====\n")
print(list.files(opscc_outdir))

# ------------------------------------------------------------
# 25. Locate the previous ssGSEA script
# ------------------------------------------------------------

list.files(
  old_outdir,
  pattern = "05B_run_ssGSEA\\.R$",
  recursive = TRUE,
  full.names = TRUE
)

# ------------------------------------------------------------
# 26. Inspect previous ssGSEA implementation
# ------------------------------------------------------------

old_ssgsea_script <-
  file.path(REPO_ROOT, "R", "reference", "05B_run_ssGSEA_original.R")

old_lines <- readLines(
  old_ssgsea_script,
  warn = FALSE
)

hit <- grep(
  "ssgsea|gsva|ssgseaParam|exprData|geneSets|normalize",
  old_lines,
  ignore.case = TRUE
)

cat(
  old_lines[
    sort(
      unique(
        unlist(
          lapply(
            hit,
            function(i) {
              max(1, i - 8):min(length(old_lines), i + 12)
            }
          )
        )
      )
    )
  ],
  sep = "\n"
)

# ------------------------------------------------------------
# 27. Inspect exact variance-filtering code
#     in the previous ssGSEA pipeline
# ------------------------------------------------------------

idx8 <- grep(
  "8\\. Remove genes with zero variance",
  old_lines
)

idx10 <- grep(
  "10\\. Build GSVA ssGSEA parameter object",
  old_lines
)

cat(
  old_lines[
    idx8:(idx10 - 1)
  ],
  sep = "\n"
)

# ------------------------------------------------------------
# 28. Re-run OPSCC ssGSEA using the exact old filtering scheme
# ------------------------------------------------------------

# Remove zero-variance genes in the 77-case OPSCC cohort
gene_variance_opscc <- apply(
  opscc_expr,
  1,
  var,
  na.rm = TRUE
)

keep_genes_opscc <- is.finite(gene_variance_opscc) &
  gene_variance_opscc > 0

opscc_expr_filtered <- opscc_expr[
  keep_genes_opscc,
  ,
  drop = FALSE
]

cat("\n===== OPSCC variance filtering =====\n")
cat("Genes before:", nrow(opscc_expr), "\n")
cat("Genes after :", nrow(opscc_expr_filtered), "\n")
cat("Removed     :", sum(!keep_genes_opscc), "\n")


# Reconfirm gene-set overlap after filtering
gene_sets_opscc_filtered <- lapply(
  gene_sets_full,
  function(genes) {
    intersect(
      genes,
      rownames(opscc_expr_filtered)
    )
  }
)

cat("\n===== Gene-set sizes after variance filtering =====\n")
print(
  vapply(
    gene_sets_opscc_filtered,
    length,
    integer(1)
  )
)

stopifnot(
  all(
    vapply(
      gene_sets_opscc_filtered,
      length,
      integer(1)
    ) >= 5
  )
)


# Run ssGSEA with old-pipeline parameters
ssgsea_param_filtered <- GSVA::ssgseaParam(
  exprData = opscc_expr_filtered,
  geneSets = gene_sets_opscc_filtered,
  minSize = 5,
  maxSize = Inf,
  alpha = 0.25,
  normalize = TRUE,
  checkNA = "auto",
  use = "everything"
)

opscc_ssgsea_filtered <- GSVA::gsva(
  ssgsea_param_filtered,
  verbose = FALSE
)

cat("\n===== Filtered ssGSEA dimensions =====\n")
print(dim(opscc_ssgsea_filtered))

cat("\n===== Compare with current OPSCC scores =====\n")
print(
  summary(
    as.vector(
      opscc_ssgsea_filtered -
        opscc_ssgsea
    )
  )
)

cat(
  "\nMaximum absolute difference:",
  max(
    abs(
      opscc_ssgsea_filtered -
        opscc_ssgsea
    )
  ),
  "\n"
)

cat(
  "\nAll.equal:",
  all.equal(
    opscc_ssgsea_filtered,
    opscc_ssgsea,
    tolerance = 1e-12
  ),
  "\n"
)

# ------------------------------------------------------------
# 29. Adopt old-pipeline-compatible ssGSEA as the formal score
# ------------------------------------------------------------

# Keep the original direct calculation for audit trail
opscc_ssgsea_direct <- opscc_ssgsea

# Adopt variance-filtered, old-pipeline-compatible scores
opscc_ssgsea <- opscc_ssgsea_filtered

cat("\n===== Formal OPSCC ssGSEA adopted =====\n")

cat(
  "Dimensions:",
  paste(dim(opscc_ssgsea), collapse = " x "),
  "\n"
)

cat(
  "Sample alignment:",
  identical(
    colnames(opscc_ssgsea),
    opscc_clin_expr$expr_id
  ),
  "\n"
)

cat(
  "Signature alignment:",
  identical(
    rownames(opscc_ssgsea),
    names(gene_sets_opscc_filtered)
  ),
  "\n"
)

cat(
  "Any NA:",
  anyNA(opscc_ssgsea),
  "\n"
)

cat(
  "Any infinite:",
  any(!is.finite(opscc_ssgsea)),
  "\n"
)

cat("\n===== Formal score summary =====\n")

print(
  t(
    apply(
      opscc_ssgsea,
      1,
      summary
    )
  )
)

# ------------------------------------------------------------
# 30. Re-run unadjusted HPV comparison
#     using the formal old-pipeline-compatible ssGSEA scores
# ------------------------------------------------------------

score_df <- as.data.frame(
  t(opscc_ssgsea),
  check.names = FALSE
)

score_df$expr_id <- rownames(score_df)

stopifnot(
  identical(
    score_df$expr_id,
    opscc_clin_expr$expr_id
  )
)

score_df$HPV_status <- opscc_clin_expr$HPV_status


unadjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      x_pos <- score_df[
        score_df$HPV_status == "HPV_positive",
        sig
      ]
      
      x_neg <- score_df[
        score_df$HPV_status == "HPV_negative",
        sig
      ]
      
      data.frame(
        Signature = sig,
        
        n_HPV_negative = length(x_neg),
        n_HPV_positive = length(x_pos),
        
        median_HPV_negative =
          median(x_neg, na.rm = TRUE),
        
        median_HPV_positive =
          median(x_pos, na.rm = TRUE),
        
        median_difference =
          median(x_pos, na.rm = TRUE) -
          median(x_neg, na.rm = TRUE),
        
        cliffs_delta =
          cliffs_delta(
            x_pos,
            x_neg
          ),
        
        raw_p =
          wilcox.test(
            x_pos,
            x_neg,
            exact = FALSE
          )$p.value,
        
        stringsAsFactors = FALSE
      )
    }
  )
)

unadjusted_results$BH_FDR <- p.adjust(
  unadjusted_results$raw_p,
  method = "BH"
)

rownames(unadjusted_results) <- NULL

cat(
  "\n===== FORMAL OPSCC HPV unadjusted results =====\n"
)

print(
  unadjusted_results,
  digits = 5,
  row.names = FALSE
)

# ------------------------------------------------------------
# 31. Re-run multivariable-adjusted analysis
#     using the formal ssGSEA scores
# ------------------------------------------------------------

analysis_df <- cbind(
  opscc_clin_expr,
  score_df[
    ,
    rownames(opscc_ssgsea),
    drop = FALSE
  ]
)

analysis_cc <- analysis_df[
  complete.cases(
    analysis_df[, model_vars]
  ),
  ,
  drop = FALSE
]

cat("\n===== FORMAL adjusted-analysis cohort =====\n")
cat("N =", nrow(analysis_cc), "\n")
print(table(analysis_cc$HPV_status))


adjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      dat$Score <- dat[[sig]]
      
      fit <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      coef_name <- "HPV_statusHPV_positive"
      
      sm <- coef(summary(fit))
      
      ci <- confint(
        fit,
        parm = coef_name,
        level = 0.95
      )
      
      data.frame(
        Signature = sig,
        N = nobs(fit),
        
        adjusted_beta =
          sm[coef_name, "Estimate"],
        
        SE =
          sm[coef_name, "Std. Error"],
        
        CI95_low =
          ci[1],
        
        CI95_high =
          ci[2],
        
        raw_p =
          sm[coef_name, "Pr(>|t|)"],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

adjusted_results$BH_FDR <- p.adjust(
  adjusted_results$raw_p,
  method = "BH"
)

rownames(adjusted_results) <- NULL

cat(
  "\n===== FORMAL OPSCC multivariable-adjusted HPV effects =====\n"
)

print(
  adjusted_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 32. Re-run model diagnostics and influence sensitivity
#     using the formal ssGSEA scores
# ------------------------------------------------------------

diagnostic_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      dat$Score <- dat[[sig]]
      
      fit <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      cd <- cooks.distance(fit)
      
      data.frame(
        Signature = sig,
        N = nobs(fit),
        R_squared = summary(fit)$r.squared,
        adjusted_R_squared = summary(fit)$adj.r.squared,
        max_Cooks_D = max(cd),
        n_Cooks_D_gt_4_over_N =
          sum(cd > 4 / nobs(fit)),
        stringsAsFactors = FALSE
      )
    }
  )
)


influence_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- analysis_cc
      dat$Score <- dat[[sig]]
      
      fit_full <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat
      )
      
      cd <- cooks.distance(fit_full)
      
      keep <- cd <= 4 / nobs(fit_full)
      
      fit_sens <- lm(
        Score ~
          HPV_status +
          age_years +
          gender +
          smoking_group +
          clinical_stage_group +
          subsite_group,
        data = dat[keep, , drop = FALSE]
      )
      
      coef_name <- "HPV_statusHPV_positive"
      
      sm <- coef(summary(fit_sens))
      
      ci <- confint(
        fit_sens,
        parm = coef_name,
        level = 0.95
      )
      
      data.frame(
        Signature = sig,
        N_full = nobs(fit_full),
        N_excluded = sum(!keep),
        N_sensitivity = nobs(fit_sens),
        
        HPV_beta_sensitivity =
          sm[coef_name, "Estimate"],
        
        CI95_low = ci[1],
        CI95_high = ci[2],
        
        raw_p =
          sm[coef_name, "Pr(>|t|)"],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

influence_results$BH_FDR <- p.adjust(
  influence_results$raw_p,
  method = "BH"
)

cat("\n===== FORMAL model diagnostics =====\n")
print(
  diagnostic_results,
  digits = 5,
  row.names = FALSE
)

cat("\n===== FORMAL influence-point sensitivity =====\n")
print(
  influence_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 33. Re-run univariable Cox models
#     using the formal ssGSEA scores
# ------------------------------------------------------------

# Add formal ssGSEA scores to survival dataset
surv_score_df <- as.data.frame(
  t(opscc_ssgsea),
  check.names = FALSE
)

surv_score_df$expr_id <- rownames(surv_score_df)

stopifnot(
  identical(
    surv_score_df$expr_id,
    opscc_clin_expr$expr_id
  )
)

# Use the existing OS variables
surv_analysis_df <- cbind(
  opscc_clin_expr,
  surv_score_df[
    ,
    rownames(opscc_ssgsea),
    drop = FALSE
  ]
)

cox_univariable_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- surv_analysis_df[
        is.finite(surv_analysis_df$OS_days) &
          !is.na(surv_analysis_df$OS_event),
        ,
        drop = FALSE
      ]
      
      # Standardize signature:
      # HR corresponds to a 1-SD increase
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~ Score_z,
        data = dat
      )
      
      sm <- summary(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cox_univariable_results$BH_FDR <- p.adjust(
  cox_univariable_results$raw_p,
  method = "BH"
)

cat(
  "\n===== FORMAL univariable Cox models =====\n"
)

print(
  cox_univariable_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 33A. Check existing OS variable names
# ------------------------------------------------------------

grep(
  "OS|surv|death|vital|follow",
  colnames(opscc_clin_expr),
  ignore.case = TRUE,
  value = TRUE
)

# ------------------------------------------------------------
# 33B. Re-run univariable Cox models
#      using OS_days and the formal ssGSEA scores
# ------------------------------------------------------------

surv_score_df <- as.data.frame(
  t(opscc_ssgsea),
  check.names = FALSE
)

surv_score_df$expr_id <- rownames(surv_score_df)

stopifnot(
  identical(
    surv_score_df$expr_id,
    opscc_clin_expr$expr_id
  )
)

surv_analysis_df <- cbind(
  opscc_clin_expr,
  surv_score_df[
    ,
    rownames(opscc_ssgsea),
    drop = FALSE
  ]
)

cox_univariable_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- surv_analysis_df[
        is.finite(surv_analysis_df$OS_days) &
          !is.na(surv_analysis_df$OS_event),
        ,
        drop = FALSE
      ]
      
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~ Score_z,
        data = dat
      )
      
      sm <- summary(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cox_univariable_results$BH_FDR <- p.adjust(
  cox_univariable_results$raw_p,
  method = "BH"
)

cat(
  "\n===== FORMAL univariable Cox models =====\n"
)

print(
  cox_univariable_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 34. Re-run HPV-adjusted Cox models
#     using the formal ssGSEA scores
# ------------------------------------------------------------

cox_HPV_adjusted_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- surv_analysis_df[
        is.finite(surv_analysis_df$OS_days) &
          !is.na(surv_analysis_df$OS_event) &
          !is.na(surv_analysis_df$HPV_status),
        ,
        drop = FALSE
      ]
      
      # Standardize signature:
      # HR corresponds to a 1-SD increase
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~
          Score_z +
          HPV_status,
        data = dat
      )
      
      sm <- summary(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_signature_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low_signature =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high_signature =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p_signature =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        HR_HPV_positive =
          sm$coefficients[
            "HPV_statusHPV_positive",
            "exp(coef)"
          ],
        
        raw_p_HPV =
          sm$coefficients[
            "HPV_statusHPV_positive",
            "Pr(>|z|)"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

# BH correction across the 8 prespecified
# signature coefficients
cox_HPV_adjusted_results$BH_FDR_signature <- p.adjust(
  cox_HPV_adjusted_results$raw_p_signature,
  method = "BH"
)

cat(
  "\n===== FORMAL HPV-adjusted Cox models =====\n"
)

print(
  cox_HPV_adjusted_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 35. Check proportional hazards assumption
#     for HPV-adjusted Cox models
# ------------------------------------------------------------

ph_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- surv_analysis_df[
        is.finite(surv_analysis_df$OS_days) &
          !is.na(surv_analysis_df$OS_event) &
          !is.na(surv_analysis_df$HPV_status),
        ,
        drop = FALSE
      ]
      
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~
          Score_z +
          HPV_status,
        data = dat,
        x = TRUE
      )
      
      zph <- survival::cox.zph(fit)
      
      data.frame(
        Signature = sig,
        
        p_signature =
          zph$table[
            "Score_z",
            "p"
          ],
        
        p_HPV =
          zph$table[
            "HPV_status",
            "p"
          ],
        
        p_global =
          zph$table[
            "GLOBAL",
            "p"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cat(
  "\n===== PH assumption: HPV-adjusted Cox models =====\n"
)

print(
  ph_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 36. HPV-stratified Cox sensitivity analysis
#     avoids PH assumption for HPV status
# ------------------------------------------------------------

cox_HPV_stratified_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      dat <- surv_analysis_df[
        is.finite(surv_analysis_df$OS_days) &
          !is.na(surv_analysis_df$OS_event) &
          !is.na(surv_analysis_df$HPV_status),
        ,
        drop = FALSE
      ]
      
      dat$Score_z <- as.numeric(
        scale(dat[[sig]])
      )
      
      fit <- survival::coxph(
        survival::Surv(
          OS_days,
          OS_event
        ) ~
          Score_z +
          strata(HPV_status),
        data = dat,
        x = TRUE
      )
      
      sm <- summary(fit)
      
      zph <- survival::cox.zph(fit)
      
      data.frame(
        Signature = sig,
        N = nrow(dat),
        Events = sum(dat$OS_event),
        
        HR_signature_per_SD =
          sm$coefficients[
            "Score_z",
            "exp(coef)"
          ],
        
        CI95_low =
          sm$conf.int[
            "Score_z",
            "lower .95"
          ],
        
        CI95_high =
          sm$conf.int[
            "Score_z",
            "upper .95"
          ],
        
        raw_p =
          sm$coefficients[
            "Score_z",
            "Pr(>|z|)"
          ],
        
        PH_p_signature =
          zph$table[
            "Score_z",
            "p"
          ],
        
        PH_p_global =
          zph$table[
            "GLOBAL",
            "p"
          ],
        
        stringsAsFactors = FALSE
      )
    }
  )
)

cox_HPV_stratified_results$BH_FDR <- p.adjust(
  cox_HPV_stratified_results$raw_p,
  method = "BH"
)

cat(
  "\n===== HPV-stratified Cox sensitivity analysis =====\n"
)

print(
  cox_HPV_stratified_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 37. Save formal survival-analysis outputs
# ------------------------------------------------------------


# ------------------------------------------------------------
# Formal output directory
# ------------------------------------------------------------
out_dir <- file.path(RESULTS_ROOT, "TCGA_OPSCC77_revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n", out_dir, "\n\n")

write.csv(
  cox_univariable_results,

  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_univariate_perSD_FORMAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_HPV_adjusted_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_perSD_FORMAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  ph_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_PH_test.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_HPV_stratified_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_stratified_sensitivity.csv"
  ),
  row.names = FALSE
)

cat("\nSaved formal survival-analysis files to:\n")
cat(out_dir, "\n")

# ------------------------------------------------------------
# 37A. Restore formal output directory
# ------------------------------------------------------------

out_dir <- file.path(RESULTS_ROOT, "TCGA_OPSCC77_revision")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n", out_dir, "\n\n")

cat(
  "Directory exists:",
  dir.exists(out_dir),
  "\n"
)

# ------------------------------------------------------------
# 37B. Save formal survival-analysis outputs
# ------------------------------------------------------------

write.csv(
  cox_univariable_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_univariate_perSD_FORMAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_HPV_adjusted_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_perSD_FORMAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  ph_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_PH_test.csv"
  ),
  row.names = FALSE
)

write.csv(
  cox_HPV_stratified_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_Cox_HPV_stratified_sensitivity.csv"
  ),
  row.names = FALSE
)

cat("\nSaved formal survival-analysis files to:\n")
cat(out_dir, "\n")

# ------------------------------------------------------------
# 38. Bootstrap 95% CI for Cliff's delta
#     unadjusted HPV comparison
# ------------------------------------------------------------

set.seed(20260911)

B <- 5000

cliff_boot_ci <- function(x_pos, x_neg, B = 5000) {
  
  obs <- cliffs_delta(x_pos, x_neg)
  
  boot_vals <- replicate(
    B,
    {
      xp <- sample(
        x_pos,
        size = length(x_pos),
        replace = TRUE
      )
      
      xn <- sample(
        x_neg,
        size = length(x_neg),
        replace = TRUE
      )
      
      cliffs_delta(xp, xn)
    }
  )
  
  ci <- quantile(
    boot_vals,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
  
  c(
    delta = obs,
    CI95_low = ci[1],
    CI95_high = ci[2]
  )
}


cliff_ci_results <- do.call(
  rbind,
  lapply(
    rownames(opscc_ssgsea),
    function(sig) {
      
      x_pos <- score_df[
        score_df$HPV_status == "HPV_positive",
        sig
      ]
      
      x_neg <- score_df[
        score_df$HPV_status == "HPV_negative",
        sig
      ]
      
      tmp <- cliff_boot_ci(
        x_pos,
        x_neg,
        B = B
      )
      
      data.frame(
        Signature = sig,
        cliffs_delta = tmp["delta"],
        CI95_low = tmp["CI95_low"],
        CI95_high = tmp["CI95_high"],
        stringsAsFactors = FALSE
      )
    }
  )
)

rownames(cliff_ci_results) <- NULL

cat(
  "\n===== Bootstrap 95% CI for Cliff's delta =====\n"
)

print(
  cliff_ci_results,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 39. Create formal unadjusted results table
#     with Cliff's delta 95% CI
# ------------------------------------------------------------

unadjusted_results_formal <- merge(
  unadjusted_results,
  cliff_ci_results[
    ,
    c(
      "Signature",
      "CI95_low",
      "CI95_high"
    )
  ],
  by = "Signature",
  all.x = TRUE,
  sort = FALSE
)

# Restore original signature order
unadjusted_results_formal <- unadjusted_results_formal[
  match(
    rownames(opscc_ssgsea),
    unadjusted_results_formal$Signature
  ),
  ,
  drop = FALSE
]

rownames(unadjusted_results_formal) <- NULL

cat(
  "\n===== FORMAL unadjusted HPV comparison with 95% CI =====\n"
)

print(
  unadjusted_results_formal,
  digits = 6,
  row.names = FALSE
)

# ------------------------------------------------------------
# 40. Save formal unadjusted HPV comparison table
# ------------------------------------------------------------

write.csv(
  unadjusted_results_formal,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_unadjusted_HPV_results_FORMAL.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSaved formal unadjusted results to:\n"
)

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC77_unadjusted_HPV_results_FORMAL.csv"
  ),
  "\n"
)

# ------------------------------------------------------------
# 41. Save formal multivariable-adjusted HPV results
# ------------------------------------------------------------

write.csv(
  adjusted_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC74_adjusted_HPV_results_FORMAL.csv"
  ),
  row.names = FALSE
)

cat(
  "\nSaved formal adjusted results to:\n"
)

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC74_adjusted_HPV_results_FORMAL.csv"
  ),
  "\n"
)

# ------------------------------------------------------------
# 42. Save formal diagnostics and influence sensitivity
# ------------------------------------------------------------

write.csv(
  diagnostic_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC74_model_diagnostics_FORMAL.csv"
  ),
  row.names = FALSE
)

write.csv(
  influence_results,
  file = file.path(
    out_dir,
    "TCGA_OPSCC74_influence_sensitivity_FORMAL.csv"
  ),
  row.names = FALSE
)

cat("\nSaved formal diagnostic files:\n")

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC74_model_diagnostics_FORMAL.csv"
  ),
  "\n"
)

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC74_influence_sensitivity_FORMAL.csv"
  ),
  "\n"
)

# ------------------------------------------------------------
# 43. Save formal old-pipeline-compatible ssGSEA scores
# ------------------------------------------------------------

write.csv(
  opscc_ssgsea,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_ssGSEA_score_matrix_FORMAL.csv"
  )
)

saveRDS(
  opscc_ssgsea,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_ssGSEA_score_matrix_FORMAL.rds"
  )
)

cat("\nSaved formal ssGSEA score matrix:\n")

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC77_ssGSEA_score_matrix_FORMAL.csv"
  ),
  "\n"
)

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC77_ssGSEA_score_matrix_FORMAL.rds"
  ),
  "\n"
)

# ------------------------------------------------------------
# 44. Freeze formal TCGA-OPSCC revision analysis
# ------------------------------------------------------------

formal_TCGA_OPSCC <- list(
  
  clinical = opscc_clin_expr,
  
  ssGSEA_scores = opscc_ssgsea,
  
  gene_sets = gene_sets_opscc,
  
  unadjusted = unadjusted_results_formal,
  
  adjusted = adjusted_results,
  
  diagnostics = diagnostic_results,
  
  influence_sensitivity = influence_results,
  
  Cox_univariable = cox_univariable_results,
  
  Cox_HPV_adjusted = cox_HPV_adjusted_results,
  
  Cox_PH_test = ph_results,
  
  Cox_HPV_stratified = cox_HPV_stratified_results
)

saveRDS(
  formal_TCGA_OPSCC,
  file = file.path(
    out_dir,
    "TCGA_OPSCC77_FORMAL_revision_analysis.rds"
  )
)

cat("\n===== FORMAL TCGA-OPSCC analysis frozen =====\n")

cat(
  file.path(
    out_dir,
    "TCGA_OPSCC77_FORMAL_revision_analysis.rds"
  ),
  "\n"
)

cat("\nObjects stored:\n")
print(names(formal_TCGA_OPSCC))
