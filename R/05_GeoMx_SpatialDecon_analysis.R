source(file.path("R", "config.R"))

GEOMX_RESULTS_DIR <- file.path(RESULTS_ROOT, "GeoMx")
dir.create(GEOMX_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 14W-1: Reload data for SpatialDecon analysis
# ============================================================

library(SpatialDecon)

base_dir <- GEOMX_ROOT

# LFR Q3-normalized expression
lfr_q3 <- readRDS(
  file.path(base_dir, "GeoMx_LFR_Q3.rds")
)

# ROI annotation
roi_anno <- read.csv(
  file.path(base_dir, "GeoMx_ROI_annotation.csv"),
  stringsAsFactors = FALSE
)

# LFR annotation
lfr_anno <- roi_anno[
  roi_anno$Region == "LFR",
]

# ImmuneCensus reference
immune_census <- download_profile_matrix(
  species = "Human",
  age_group = "Adult",
  matrixname = "ImmuneCensus_HCA"
)

# GeoMx background
lfr_bg <- SpatialDecon::derive_GeoMx_background(
  norm = lfr_q3,
  probepool = rep(1, nrow(lfr_q3)),
  negnames = "NegProbe-WTX"
)

cat("\n===== RELOAD CHECK =====\n")
print(dim(lfr_q3))
print(dim(lfr_bg))
print(dim(immune_census))

# ============================================================
# STEP 14W-2: Re-run safeTME SpatialDecon and capture warnings
# ============================================================

lfr_decon_safeTME <- SpatialDecon::spatialdecon(
  norm = lfr_q3,
  bg = lfr_bg,
  X = safeTME,
  align_genes = TRUE
)

cat("\n===== WARNINGS IMMEDIATELY AFTER RUN =====\n")
print(warnings())

cat("\n===== BASIC OUTPUT CHECK =====\n")
print(dim(lfr_decon_safeTME$beta))
print(dim(lfr_decon_safeTME$prop_of_all))

cat("\n===== NA CHECK =====\n")
cat("beta NA:",
    sum(is.na(lfr_decon_safeTME$beta)), "\n")
cat("prop_of_all NA:",
    sum(is.na(lfr_decon_safeTME$prop_of_all)), "\n")

# ============================================================
# STEP 14AH: Compare full 36-ROI results at maxit 1000 vs 10000
# ============================================================

warn_default <- character()

lfr_safeTME_default <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 1000
  ),
  warning = function(w) {
    warn_default <<- c(warn_default, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

warn_high <- character()

lfr_safeTME_high <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 10000
  ),
  warning = function(w) {
    warn_high <<- c(warn_high, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("\n===== WARNING COUNTS =====\n")
cat("default:", length(warn_default), "\n")
cat("high maxit:", length(warn_high), "\n")

cat("\n===== MAX ABS DIFFERENCE: BETA =====\n")
print(
  max(
    abs(lfr_safeTME_default$beta - lfr_safeTME_high$beta),
    na.rm = TRUE
  )
)

cat("\n===== MAX ABS DIFFERENCE: PROP_OF_ALL =====\n")
print(
  max(
    abs(lfr_safeTME_default$prop_of_all -
          lfr_safeTME_high$prop_of_all),
    na.rm = TRUE
  )
)

cat("\n===== ALL.EQUAL BETA =====\n")
print(
  all.equal(
    lfr_safeTME_default$beta,
    lfr_safeTME_high$beta,
    tolerance = 1e-6
  )
)

cat("\n===== ALL.EQUAL PROP_OF_ALL =====\n")
print(
  all.equal(
    lfr_safeTME_default$prop_of_all,
    lfr_safeTME_high$prop_of_all,
    tolerance = 1e-6
  )
)

# ============================================================
# STEP 14AH: Compare full 36-ROI results at maxit 1000 vs 10000
# ============================================================

warn_default <- character()

lfr_safeTME_default <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 1000
  ),
  warning = function(w) {
    warn_default <<- c(warn_default, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

warn_high <- character()

lfr_safeTME_high <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 10000
  ),
  warning = function(w) {
    warn_high <<- c(warn_high, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("\n===== WARNING COUNTS =====\n")
cat("default:", length(warn_default), "\n")
cat("high maxit:", length(warn_high), "\n")

cat("\n===== MAX ABS DIFFERENCE: BETA =====\n")
print(
  max(
    abs(lfr_safeTME_default$beta - lfr_safeTME_high$beta),
    na.rm = TRUE
  )
)

cat("\n===== MAX ABS DIFFERENCE: PROP_OF_ALL =====\n")
print(
  max(
    abs(lfr_safeTME_default$prop_of_all -
          lfr_safeTME_high$prop_of_all),
    na.rm = TRUE
  )
)

cat("\n===== ALL.EQUAL BETA =====\n")
print(
  all.equal(
    lfr_safeTME_default$beta,
    lfr_safeTME_high$beta,
    tolerance = 1e-6
  )
)

cat("\n===== ALL.EQUAL PROP_OF_ALL =====\n")
print(
  all.equal(
    lfr_safeTME_default$prop_of_all,
    lfr_safeTME_high$prop_of_all,
    tolerance = 1e-6
  )
)

# ============================================================
# STEP 14AH-0: Reload objects needed for SpatialDecon
# ============================================================

library(SpatialDecon)

base_dir <- GEOMX_ROOT

lfr_q3 <- readRDS(
  file.path(base_dir, "GeoMx_LFR_Q3.rds")
)

lfr_bg <- SpatialDecon::derive_GeoMx_background(
  norm = lfr_q3,
  probepool = rep(1, nrow(lfr_q3)),
  negnames = "NegProbe-WTX"
)

cat("\n===== RELOAD CHECK =====\n")
print(dim(lfr_q3))
print(dim(lfr_bg))

# ============================================================
# STEP 14AH: Compare full 36-ROI results at maxit 1000 vs 10000
# ============================================================

warn_default <- character()

lfr_safeTME_default <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 1000
  ),
  warning = function(w) {
    warn_default <<- c(warn_default, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

warn_high <- character()

lfr_safeTME_high <- withCallingHandlers(
  SpatialDecon::spatialdecon(
    norm = lfr_q3,
    bg = lfr_bg,
    X = safeTME,
    align_genes = TRUE,
    maxit = 10000
  ),
  warning = function(w) {
    warn_high <<- c(warn_high, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("\n===== WARNING COUNTS =====\n")
cat("default:", length(warn_default), "\n")
cat("high maxit:", length(warn_high), "\n")

cat("\n===== MAX ABS DIFFERENCE: BETA =====\n")
print(
  max(
    abs(lfr_safeTME_default$beta - lfr_safeTME_high$beta),
    na.rm = TRUE
  )
)

cat("\n===== MAX ABS DIFFERENCE: PROP_OF_ALL =====\n")
print(
  max(
    abs(
      lfr_safeTME_default$prop_of_all -
        lfr_safeTME_high$prop_of_all
    ),
    na.rm = TRUE
  )
)

cat("\n===== ALL.EQUAL BETA =====\n")
print(
  all.equal(
    lfr_safeTME_default$beta,
    lfr_safeTME_high$beta,
    tolerance = 1e-6
  )
)

cat("\n===== ALL.EQUAL PROP_OF_ALL =====\n")
print(
  all.equal(
    lfr_safeTME_default$prop_of_all,
    lfr_safeTME_high$prop_of_all,
    tolerance = 1e-6
  )
)

# ============================================================
# STEP 14AI: Inspect safeTME estimated proportions
# ============================================================

prop_safeTME <- lfr_safeTME_default$prop_of_all

cat("\n===== CELL TYPES =====\n")
print(rownames(prop_safeTME))

cat("\n===== COLUMN SUMS OF prop_of_all =====\n")
print(summary(colSums(prop_safeTME)))

cat("\n===== FIRST 6 ROI COLUMN SUMS =====\n")
print(colSums(prop_safeTME)[1:6])

cat("\n===== MEDIAN prop_of_all BY CELL TYPE =====\n")
print(
  sort(
    apply(prop_safeTME, 1, median, na.rm = TRUE),
    decreasing = TRUE
  )
)

# ============================================================
# STEP 14AJ: Create broad lymphoid composition estimates
# ============================================================

lfr_decon_broad <- data.frame(
  ROI_ID = colnames(prop_safeTME),
  
  B_cell_prop =
    prop_safeTME["B.naive", ] +
    prop_safeTME["B.memory", ],
  
  CD4_T_prop =
    prop_safeTME["T.CD4.naive", ] +
    prop_safeTME["T.CD4.memory", ],
  
  CD8_T_prop =
    prop_safeTME["T.CD8.naive", ] +
    prop_safeTME["T.CD8.memory", ],
  
  Treg_prop =
    prop_safeTME["Treg", ],
  
  stringsAsFactors = FALSE
)

lfr_decon_broad$Total_T_prop <-
  lfr_decon_broad$CD4_T_prop +
  lfr_decon_broad$CD8_T_prop +
  lfr_decon_broad$Treg_prop

cat("\n===== BROAD LYMPHOID ESTIMATES =====\n")
print(summary(
  lfr_decon_broad[
    c(
      "B_cell_prop",
      "CD4_T_prop",
      "CD8_T_prop",
      "Treg_prop",
      "Total_T_prop"
    )
  ]
))

cat("\n===== FIRST 6 ROIs =====\n")
print(head(lfr_decon_broad))

# ============================================================
# STEP 14AK-0: Reload LFR ROI annotation
# ============================================================

roi_anno <- read.csv(
  file.path(base_dir, "GeoMx_ROI_annotation.csv"),
  stringsAsFactors = FALSE
)

lfr_anno <- roi_anno[
  roi_anno$Region == "LFR",
  ,
  drop = FALSE
]

cat("\n===== LFR ANNOTATION CHECK =====\n")
print(dim(lfr_anno))
print(table(lfr_anno$Case, lfr_anno$HPV))
print(all(colnames(lfr_q3) %in% lfr_anno$ROI_ID))

# ============================================================
# STEP 14AK: Aggregate SpatialDecon estimates to patient level
# ============================================================

# Match ROI annotation
lfr_decon_anno <- merge(
  lfr_decon_broad,
  lfr_anno[, c("ROI_ID", "Case", "HPV")],
  by = "ROI_ID",
  all.x = TRUE,
  sort = FALSE
)

cat("\n===== ANNOTATION MATCH CHECK =====\n")
print(table(is.na(lfr_decon_anno$Case)))
print(table(lfr_decon_anno$Case, lfr_decon_anno$HPV))

# Patient-level mean of 6 LFR ROIs
lfr_decon_case <- aggregate(
  cbind(
    B_cell_prop,
    CD4_T_prop,
    CD8_T_prop,
    Treg_prop,
    Total_T_prop
  ) ~ Case + HPV,
  data = lfr_decon_anno,
  FUN = mean,
  na.rm = TRUE
)

# Number of ROIs per patient
n_roi_case <- aggregate(
  ROI_ID ~ Case + HPV,
  data = lfr_decon_anno,
  FUN = length
)

names(n_roi_case)[3] <- "n_ROI"

lfr_decon_case <- merge(
  lfr_decon_case,
  n_roi_case,
  by = c("Case", "HPV"),
  sort = FALSE
)

lfr_decon_case <- lfr_decon_case[
  order(lfr_decon_case$Case),
]

cat("\n===== PATIENT-LEVEL SpatialDecon =====\n")
print(lfr_decon_case)

# ============================================================
# STEP 14AL: Patient-level HPV comparison of SpatialDecon estimates
# ============================================================

decon_vars <- c(
  "B_cell_prop",
  "CD4_T_prop",
  "CD8_T_prop",
  "Total_T_prop"
)

# Cliff's delta: positive = higher in HPV-positive
cliff_delta_simple <- function(x_pos, x_neg) {
  
  comp <- outer(x_pos, x_neg, FUN = "-")
  
  (
    sum(comp > 0) -
      sum(comp < 0)
  ) / length(comp)
}

decon_hpv_results <- do.call(
  rbind,
  lapply(
    decon_vars,
    function(v) {
      
      x_pos <- lfr_decon_case[
        lfr_decon_case$HPV == "positive",
        v
      ]
      
      x_neg <- lfr_decon_case[
        lfr_decon_case$HPV == "negative",
        v
      ]
      
      wt <- wilcox.test(
        x_pos,
        x_neg,
        exact = TRUE
      )
      
      data.frame(
        Cell_type = v,
        Median_HPV_negative = median(x_neg),
        Median_HPV_positive = median(x_pos),
        Median_difference =
          median(x_pos) - median(x_neg),
        Cliff_delta =
          cliff_delta_simple(x_pos, x_neg),
        P_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

decon_hpv_results$FDR_BH <-
  p.adjust(
    decon_hpv_results$P_value,
    method = "BH"
  )

cat("\n===== PATIENT-LEVEL HPV COMPARISON =====\n")
print(decon_hpv_results)

# ============================================================
# STEP 14AM: Link LFR stem-like score with SpatialDecon estimates
# ============================================================

case_ssgsea <- read.csv(
  file.path(
    GEOMX_ROOT,
    "GeoMx_ssGSEA_case_level_scores.csv"
  ),
  stringsAsFactors = FALSE
)

stem_case <- case_ssgsea[
  case_ssgsea$Compartment == "LFR" &
    case_ssgsea$Signature == "STEMLIKE_CORE16_TCF7_IL7R",
  c("Case", "HPV", "Score")
]

names(stem_case)[3] <- "Stem_like_score"

stem_decon <- merge(
  stem_case,
  lfr_decon_case[
    ,
    c(
      "Case",
      "B_cell_prop",
      "CD4_T_prop",
      "CD8_T_prop",
      "Total_T_prop"
    )
  ],
  by = "Case",
  sort = FALSE
)

stem_decon <- stem_decon[
  order(stem_decon$Case),
]

cat("\n===== STEM-LIKE + SpatialDecon PATIENT DATA =====\n")
print(stem_decon)

cat("\n===== SPEARMAN CORRELATIONS =====\n")

for (v in c(
  "B_cell_prop",
  "CD4_T_prop",
  "CD8_T_prop",
  "Total_T_prop"
)) {
  
  ct <- cor.test(
    stem_decon$Stem_like_score,
    stem_decon[[v]],
    method = "spearman",
    exact = FALSE
  )
  
  cat(
    v,
    ": rho =", unname(ct$estimate),
    ", p =", ct$p.value,
    "\n"
  )
}

# ============================================================
# STEP 14AN: Exploratory composition-adjusted sensitivity models
# ============================================================

stem_decon$HPV <- factor(
  stem_decon$HPV,
  levels = c("negative", "positive")
)

# HPV + CD8 T-cell proportion
fit_stem_CD8 <- lm(
  Stem_like_score ~ HPV + CD8_T_prop,
  data = stem_decon
)

# HPV + total T-cell proportion
fit_stem_TotalT <- lm(
  Stem_like_score ~ HPV + Total_T_prop,
  data = stem_decon
)

cat("\n===== HPV + CD8 T MODEL =====\n")
print(summary(fit_stem_CD8)$coefficients)

cat("\n===== 95% CI: HPV + CD8 T MODEL =====\n")
print(confint(fit_stem_CD8))

cat("\n===== HPV + TOTAL T MODEL =====\n")
print(summary(fit_stem_TotalT)$coefficients)

cat("\n===== 95% CI: HPV + TOTAL T MODEL =====\n")
print(confint(fit_stem_TotalT))

# ============================================================
# STEP 14AO: Save formal SpatialDecon composition results
# ============================================================

write.csv(
  lfr_decon_case,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_patient_level.csv"),
  row.names = FALSE
)

write.csv(
  decon_hpv_results,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_HPV_comparison.csv"),
  row.names = FALSE
)

write.csv(
  stem_decon,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_SpatialDecon_patient_level.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(summary(fit_stem_CD8)$coefficients),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_adjusted_CD8_SpatialDecon.csv")
)

write.csv(
  as.data.frame(summary(fit_stem_TotalT)$coefficients),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_adjusted_TotalT_SpatialDecon.csv")
)

saveRDS(
  list(
    SpatialDecon_default = lfr_safeTME_default,
    SpatialDecon_high_maxit = lfr_safeTME_high,
    patient_level = lfr_decon_case,
    HPV_comparison = decon_hpv_results,
    stem_decon = stem_decon,
    fit_stem_CD8 = fit_stem_CD8,
    fit_stem_TotalT = fit_stem_TotalT,
    warnings_default = warn_default,
    warnings_high_maxit = warn_high
  ),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_STEP14.rds")
)

cat("\n===== SAVE CHECK =====\n")

files_to_check <- c(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_patient_level.csv"),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_HPV_comparison.csv"),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_SpatialDecon_patient_level.csv"),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_adjusted_CD8_SpatialDecon.csv"),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_adjusted_TotalT_SpatialDecon.csv"),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_STEP14.rds")
)

print(file.exists(files_to_check))

# ============================================================
# STEP 14AP: Save SpatialDecon convergence/stability diagnostics
# ============================================================

# Compare estimates between default maxit and high maxit
beta_max_abs_diff <- max(
  abs(
    lfr_safeTME_default$beta -
      lfr_safeTME_high$beta
  ),
  na.rm = TRUE
)

prop_max_abs_diff <- max(
  abs(
    lfr_safeTME_default$prop_of_all -
      lfr_safeTME_high$prop_of_all
  ),
  na.rm = TRUE
)

beta_identical <- isTRUE(
  all.equal(
    lfr_safeTME_default$beta,
    lfr_safeTME_high$beta
  )
)

prop_identical <- isTRUE(
  all.equal(
    lfr_safeTME_default$prop_of_all,
    lfr_safeTME_high$prop_of_all
  )
)

# Warning counts
n_warning_default <- length(warn_default)
n_warning_high <- length(warn_high)

spatialdecon_stability <- data.frame(
  Analysis = c(
    "Warning_count",
    "Beta_max_abs_difference",
    "Prop_of_all_max_abs_difference",
    "Beta_all_equal",
    "Prop_of_all_all_equal"
  ),
  maxit_1000 = c(
    n_warning_default,
    NA,
    NA,
    NA,
    NA
  ),
  maxit_10000 = c(
    n_warning_high,
    beta_max_abs_diff,
    prop_max_abs_diff,
    beta_identical,
    prop_identical
  ),
  stringsAsFactors = FALSE
)

cat("\n===== SpatialDecon STABILITY DIAGNOSTICS =====\n")
print(spatialdecon_stability)

cat("\n===== PACKAGE VERSION =====\n")
print(packageVersion("SpatialDecon"))

write.csv(
  spatialdecon_stability,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_convergence_stability.csv"),
  row.names = FALSE
)

cat("\n===== SAVE CHECK =====\n")
print(
  file.exists(
    file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_convergence_stability.csv")
  )
)

# ============================================================
# STEP 15B: Create GeoMx revision master summary
# ============================================================

# Separate HPV comparisons by compartment
hpv_lfr <- hpv_res[
  hpv_res$Compartment == "LFR",
  c(
    "Signature",
    "median_difference",
    "cliffs_delta",
    "p_value",
    "FDR_BH",
    "Cliff_CI_lower",
    "Cliff_CI_upper"
  )
]

names(hpv_lfr)[-1] <- paste0(
  "LFR_HPV_",
  names(hpv_lfr)[-1]
)

hpv_ttr <- hpv_res[
  hpv_res$Compartment == "TTR",
  c(
    "Signature",
    "median_difference",
    "cliffs_delta",
    "p_value",
    "FDR_BH",
    "Cliff_CI_lower",
    "Cliff_CI_upper"
  )
]

names(hpv_ttr)[-1] <- paste0(
  "TTR_HPV_",
  names(hpv_ttr)[-1]
)

# Paired compartment results
paired_master <- paired_res[
  ,
  c(
    "Signature",
    "median_LFR_minus_TTR",
    "n_LFR_higher",
    "n_TTR_higher",
    "p_value",
    "FDR_BH"
  )
]

names(paired_master)[-1] <- paste0(
  "Paired_",
  names(paired_master)[-1]
)

# HPV x compartment interaction
interaction_master <- interaction_res[
  ,
  c(
    "Signature",
    "interaction_estimate",
    "CI_lower",
    "CI_upper",
    "exact_permutation_p",
    "permutation_FDR_BH"
  )
]

names(interaction_master)[-1] <- paste0(
  "Interaction_",
  names(interaction_master)[-1]
)

# Merge into one master table
geomx_master <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = "Signature",
    all = TRUE,
    sort = FALSE
  ),
  list(
    hpv_lfr,
    hpv_ttr,
    paired_master,
    interaction_master
  )
)

cat("\n===== GEOMX REVISION MASTER SUMMARY =====\n")
print(geomx_master)

cat("\n===== DIMENSION =====\n")
print(dim(geomx_master))

write.csv(
  geomx_master,
  "GeoMx_revision_master_summary_FORMAL.csv",
  row.names = FALSE
)

cat("\n===== SAVE CHECK =====\n")
print(
  file.exists(
    "GeoMx_revision_master_summary_FORMAL.csv"
  )
)

# ============================================================
# STEP 15C: Create SpatialDecon revision master summary
# ============================================================

# Read formal saved outputs
decon_hpv <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_HPV_comparison.csv"),
  stringsAsFactors = FALSE
)

stem_decon_case <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_SpatialDecon_patient_level.csv"),
  stringsAsFactors = FALSE
)

stability_res <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_convergence_stability.csv"),
  stringsAsFactors = FALSE
)

# Spearman associations: Stem-like vs broad composition
spearman_vars <- c(
  "B_cell_prop",
  "CD4_T_prop",
  "CD8_T_prop",
  "Total_T_prop"
)

spearman_res <- do.call(
  rbind,
  lapply(
    spearman_vars,
    function(v) {
      
      ct <- cor.test(
        stem_decon_case$Stem_like_score,
        stem_decon_case[[v]],
        method = "spearman",
        exact = FALSE
      )
      
      data.frame(
        Cell_type = v,
        Spearman_rho = unname(ct$estimate),
        Spearman_p = ct$p.value,
        stringsAsFactors = FALSE
      )
    }
  )
)

# Exploratory HPV-adjusted models
coef_cd8 <- summary(fit_stem_CD8)$coefficients
ci_cd8 <- confint(fit_stem_CD8)

coef_totalT <- summary(fit_stem_TotalT)$coefficients
ci_totalT <- confint(fit_stem_TotalT)

adjusted_models <- data.frame(
  Model = c(
    "Stem_like ~ HPV + CD8_T_prop",
    "Stem_like ~ HPV + Total_T_prop"
  ),
  HPV_estimate = c(
    coef_cd8["HPVpositive", "Estimate"],
    coef_totalT["HPVpositive", "Estimate"]
  ),
  HPV_CI_lower = c(
    ci_cd8["HPVpositive", 1],
    ci_totalT["HPVpositive", 1]
  ),
  HPV_CI_upper = c(
    ci_cd8["HPVpositive", 2],
    ci_totalT["HPVpositive", 2]
  ),
  HPV_p_value = c(
    coef_cd8["HPVpositive", "Pr(>|t|)"],
    coef_totalT["HPVpositive", "Pr(>|t|)"]
  ),
  stringsAsFactors = FALSE
)

cat("\n===== SPATIALDECON HPV COMPARISON =====\n")
print(decon_hpv)

cat("\n===== STEM-LIKE ASSOCIATIONS =====\n")
print(spearman_res)

cat("\n===== COMPOSITION-ADJUSTED STEM-LIKE MODELS =====\n")
print(adjusted_models)

cat("\n===== CONVERGENCE/STABILITY =====\n")
print(stability_res)

# Save separate formal master components
write.csv(
  spearman_res,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_StemLike_correlations_FORMAL.csv"),
  row.names = FALSE
)

write.csv(
  adjusted_models,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_StemLike_adjusted_models_FORMAL.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    HPV_comparison = decon_hpv,
    StemLike_correlations = spearman_res,
    StemLike_adjusted_models = adjusted_models,
    Stability = stability_res
  ),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_master_summary_FORMAL.rds")
)

cat("\n===== SAVE CHECK =====\n")
print(
  file.exists(c(
    file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_StemLike_correlations_FORMAL.csv"),
    file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_StemLike_adjusted_models_FORMAL.csv"),
    file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_master_summary_FORMAL.rds")
  ))
)
