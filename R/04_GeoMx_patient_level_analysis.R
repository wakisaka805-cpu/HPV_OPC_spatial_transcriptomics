source(file.path("R", "config.R"))

GEOMX_RESULTS_DIR <- file.path(RESULTS_ROOT, "GeoMx")
dir.create(GEOMX_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: Load GeoMx revision data and check structure
# ============================================================

roi_df <- read.csv(
  file.path(GEOMX_ROOT, "GeoMx_ssGSEA_ROI_level_scores.csv"),
  stringsAsFactors = FALSE
)

case_df <- read.csv(
  file.path(GEOMX_ROOT, "GeoMx_ssGSEA_case_level_scores.csv"),
  stringsAsFactors = FALSE
)

cat("\n===== ROI-level data =====\n")
print(dim(roi_df))
print(names(roi_df))
print(head(roi_df))

cat("\n===== Case-level data =====\n")
print(dim(case_df))
print(names(case_df))
print(head(case_df))

cat("\n===== ROI-level counts =====\n")
print(table(roi_df$Compartment))
print(table(roi_df$Case))
print(table(roi_df$HPV))

cat("\n===== Case-level counts =====\n")
print(table(case_df$Compartment))
print(table(case_df$Case))
print(table(case_df$HPV))

cat("\n===== Signatures =====\n")
print(unique(case_df$Signature))

# ============================================================
# STEP 2: Verify case-level scores are ROI means
# ============================================================

library(dplyr)

case_recalc <- roi_df %>%
  group_by(Case, HPV, Compartment, Signature) %>%
  summarise(
    Score_recalc = mean(Score, na.rm = TRUE),
    n_ROI = n(),
    .groups = "drop"
  )

check_case <- case_df %>%
  left_join(
    case_recalc,
    by = c("Case", "HPV", "Compartment", "Signature")
  ) %>%
  mutate(
    diff = Score - Score_recalc
  )

cat("\n===== ROI count per case/compartment/signature =====\n")
print(table(case_recalc$n_ROI))

cat("\n===== Case-level score agreement =====\n")
cat("Number of rows:", nrow(check_case), "\n")
cat("Maximum absolute difference:",
    max(abs(check_case$diff), na.rm = TRUE), "\n")

cat("\n===== Largest differences =====\n")
print(
  check_case %>%
    arrange(desc(abs(diff))) %>%
    select(Case, HPV, Compartment, Signature, Score, Score_recalc, n_ROI, diff) %>%
    head(10)
)

# ============================================================
# STEP 3: Define the 8 prespecified signatures
# ============================================================

target_signatures <- c(
  "STEMLIKE_CORE16_TCF7_IL7R",
  "PRE_EXHAUSTION_TPEX18",
  "T_CELL_ACTIVITY_PAN",
  "B_CELL_ACTIVITY_PAN",
  "TEX_LIKE_PROGRAM_SHI2023",
  "ESTROGEN_RESPONSE_EARLY",
  "ESTROGEN_RESPONSE_LATE",
  "BILE_ACID_METABOLISM"
)

case_main <- case_df %>%
  filter(Signature %in% target_signatures)

roi_main <- roi_df %>%
  filter(Signature %in% target_signatures)

cat("\n===== Main analysis dimensions =====\n")
cat("Case-level:", dim(case_main), "\n")
cat("ROI-level :", dim(roi_main), "\n")

cat("\n===== Signatures retained =====\n")
print(sort(unique(case_main$Signature)))

cat("\n===== Expected row counts =====\n")
cat("Case-level expected = 6 cases x 2 compartments x 8 signatures = 96\n")
cat("ROI-level expected  = 6 cases x 2 compartments x 6 ROIs x 8 signatures = 576\n")

# ============================================================
# STEP 4: Patient-level primary HPV comparison
#         LFR and TTR separately
# ============================================================

# Cliff's delta:
# positive value = higher in HPV-positive
cliffs_delta <- function(x_pos, x_neg) {
  comp <- outer(x_pos, x_neg, FUN = "-")
  (sum(comp > 0) - sum(comp < 0)) / length(comp)
}

patient_primary <- case_main %>%
  group_by(Compartment, Signature) %>%
  summarise(
    n_negative = sum(HPV == "negative"),
    n_positive = sum(HPV == "positive"),
    
    median_negative = median(Score[HPV == "negative"]),
    median_positive = median(Score[HPV == "positive"]),
    
    median_difference =
      median(Score[HPV == "positive"]) -
      median(Score[HPV == "negative"]),
    
    cliffs_delta = cliffs_delta(
      Score[HPV == "positive"],
      Score[HPV == "negative"]
    ),
    
    p_value = wilcox.test(
      Score[HPV == "positive"],
      Score[HPV == "negative"],
      exact = TRUE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  group_by(Compartment) %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup()

cat("\n===== PATIENT-LEVEL PRIMARY HPV COMPARISON =====\n")
print(
  patient_primary,
)

# ============================================================
# STEP 5: Create within-patient compartment difference
#         Delta = LFR - TTR
# ============================================================

delta_df <- case_main %>%
  select(Case, HPV, Signature, Compartment, Score) %>%
  tidyr::pivot_wider(
    names_from = Compartment,
    values_from = Score
  ) %>%
  mutate(
    Delta_LFR_minus_TTR = LFR - TTR
  )

cat("\n===== WITHIN-PATIENT DELTA: LFR - TTR =====\n")
print(
  delta_df %>%
    arrange(Signature, HPV, Case),
)

cat("\n===== Structure check =====\n")
cat("Rows:", nrow(delta_df), "\n")
print(table(delta_df$HPV))
cat("Missing Delta:",
    sum(is.na(delta_df$Delta_LFR_minus_TTR)), "\n")

# ============================================================
# STEP 6: HPV effect on within-patient compartment difference
#         Delta = LFR - TTR
# ============================================================

delta_stats <- delta_df %>%
  group_by(Signature) %>%
  summarise(
    n_negative = sum(HPV == "negative"),
    n_positive = sum(HPV == "positive"),
    
    median_delta_negative =
      median(Delta_LFR_minus_TTR[HPV == "negative"]),
    
    median_delta_positive =
      median(Delta_LFR_minus_TTR[HPV == "positive"]),
    
    difference_in_median_delta =
      median(Delta_LFR_minus_TTR[HPV == "positive"]) -
      median(Delta_LFR_minus_TTR[HPV == "negative"]),
    
    cliffs_delta = cliffs_delta(
      Delta_LFR_minus_TTR[HPV == "positive"],
      Delta_LFR_minus_TTR[HPV == "negative"]
    ),
    
    p_value = wilcox.test(
      Delta_LFR_minus_TTR[HPV == "positive"],
      Delta_LFR_minus_TTR[HPV == "negative"],
      exact = TRUE
    )$p.value,
    
    .groups = "drop"
  ) %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH")
  )

cat("\n===== HPV EFFECT ON WITHIN-PATIENT DELTA =====\n")
print(
  delta_stats,
)

# ============================================================
# STEP 7: Formal patient-level HPV x compartment interaction
#         Difference-score model
#         Delta = LFR - TTR
# ============================================================

interaction_results <- lapply(
  target_signatures,
  function(sig) {
    
    d <- delta_df %>%
      filter(Signature == sig) %>%
      mutate(
        HPV = factor(HPV, levels = c("negative", "positive"))
      )
    
    fit <- lm(Delta_LFR_minus_TTR ~ HPV, data = d)
    
    sm <- summary(fit)$coefficients
    ci <- confint(fit, level = 0.95)
    
    data.frame(
      Signature = sig,
      n = nrow(d),
      interaction_estimate = sm["HPVpositive", "Estimate"],
      SE = sm["HPVpositive", "Std. Error"],
      CI_lower = ci["HPVpositive", 1],
      CI_upper = ci["HPVpositive", 2],
      p_value = sm["HPVpositive", "Pr(>|t|)"]
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH")
  )

cat("\n===== FORMAL PATIENT-LEVEL HPV x COMPARTMENT INTERACTION =====\n")
print(
  interaction_results,
)

saveRDS(
  list(
    case_main = case_main,
    roi_main = roi_main,
    delta_df = delta_df,
    patient_primary = patient_primary,
    delta_stats = delta_stats,
    interaction_results = interaction_results
  ),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_revision_STEP7.rds")
)

# ============================================================
# STEP 8: Exact patient-level permutation test
#         for HPV x compartment interaction
#         Delta = LFR - TTR
# ============================================================

exact_perm_p <- function(delta, hpv) {
  
  # observed difference in means:
  # HPV-positive minus HPV-negative
  obs <- mean(delta[hpv == "positive"]) -
    mean(delta[hpv == "negative"])
  
  n <- length(delta)
  n_pos <- sum(hpv == "positive")
  
  # all possible assignments of 4 positive patients among 6
  combs <- combn(seq_len(n), n_pos)
  
  perm_stats <- apply(
    combs,
    2,
    function(pos_idx) {
      
      neg_idx <- setdiff(seq_len(n), pos_idx)
      
      mean(delta[pos_idx]) -
        mean(delta[neg_idx])
    }
  )
  
  # two-sided exact permutation P
  p_exact <- mean(abs(perm_stats) >= abs(obs) - 1e-12)
  
  c(
    observed_mean_difference = obs,
    exact_permutation_p = p_exact
  )
}

interaction_perm <- lapply(
  target_signatures,
  function(sig) {
    
    d <- delta_df %>%
      filter(Signature == sig)
    
    res <- exact_perm_p(
      delta = d$Delta_LFR_minus_TTR,
      hpv = d$HPV
    )
    
    data.frame(
      Signature = sig,
      mean_delta_difference =
        unname(res["observed_mean_difference"]),
      permutation_p =
        unname(res["exact_permutation_p"])
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    permutation_FDR_BH =
      p.adjust(permutation_p, method = "BH")
  )

cat("\n===== EXACT PATIENT-LEVEL PERMUTATION TEST =====\n")
print(as.data.frame(interaction_perm))

# ============================================================
# STEP 9: Final patient-level interaction table
#         Estimate + 95% CI + parametric P
#         + exact permutation P + BH-FDR
# ============================================================

interaction_final <- interaction_results %>%
  select(
    Signature,
    n,
    interaction_estimate,
    SE,
    CI_lower,
    CI_upper,
    lm_p_value = p_value
  ) %>%
  left_join(
    interaction_perm %>%
      select(
        Signature,
        exact_permutation_p = permutation_p,
        permutation_FDR_BH
      ),
    by = "Signature"
  )

cat("\n===== FINAL PATIENT-LEVEL INTERACTION TABLE =====\n")
print(as.data.frame(interaction_final))

write.csv(
  interaction_final,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_HPVxCompartment_interaction_FORMAL.csv"),
  row.names = FALSE
)

# ============================================================
# STEP 10: Bootstrap 95% CI for Cliff's delta
#          Patient-level primary HPV comparison
# ============================================================

set.seed(12345)

bootstrap_cliff_ci <- function(x_pos, x_neg, B = 10000) {
  
  boot_delta <- replicate(
    B,
    {
      xb_pos <- sample(
        x_pos,
        size = length(x_pos),
        replace = TRUE
      )
      
      xb_neg <- sample(
        x_neg,
        size = length(x_neg),
        replace = TRUE
      )
      
      cliffs_delta(xb_pos, xb_neg)
    }
  )
  
  quantile(
    boot_delta,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )
}

cliff_ci_results <- lapply(
  seq_len(nrow(patient_primary)),
  function(i) {
    
    comp <- patient_primary$Compartment[i]
    sig  <- patient_primary$Signature[i]
    
    d <- case_main %>%
      filter(
        Compartment == comp,
        Signature == sig
      )
    
    ci <- bootstrap_cliff_ci(
      x_pos = d$Score[d$HPV == "positive"],
      x_neg = d$Score[d$HPV == "negative"],
      B = 10000
    )
    
    data.frame(
      Compartment = comp,
      Signature = sig,
      Cliff_CI_lower = ci[1],
      Cliff_CI_upper = ci[2]
    )
  }
) %>%
  bind_rows()

patient_primary_formal <- patient_primary %>%
  left_join(
    cliff_ci_results,
    by = c("Compartment", "Signature")
  )

cat("\n===== PATIENT-LEVEL PRIMARY HPV COMPARISON WITH 95% CI =====\n")
print(as.data.frame(patient_primary_formal))

write.csv(
  patient_primary_formal,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_HPV_comparison_FORMAL.csv"),
  row.names = FALSE
)

# ============================================================
# STEP 11: Exploratory ROI-level direction check
#          Compare direction with patient-level primary results
# ============================================================

roi_exploratory <- lapply(
  unique(roi_main$Compartment),
  function(comp) {
    
    lapply(
      target_signatures,
      function(sig) {
        
        d <- roi_main %>%
          filter(
            Compartment == comp,
            Signature == sig
          )
        
        med_neg <- median(
          d$Score[d$HPV == "negative"],
          na.rm = TRUE
        )
        
        med_pos <- median(
          d$Score[d$HPV == "positive"],
          na.rm = TRUE
        )
        
        delta_roi <- cliffs_delta(
          d$Score[d$HPV == "positive"],
          d$Score[d$HPV == "negative"]
        )
        
        data.frame(
          Compartment = comp,
          Signature = sig,
          ROI_median_negative = med_neg,
          ROI_median_positive = med_pos,
          ROI_median_difference = med_pos - med_neg,
          ROI_cliffs_delta = delta_roi
        )
      }
    ) %>%
      bind_rows()
  }
) %>%
  bind_rows()

roi_vs_patient <- patient_primary_formal %>%
  select(
    Compartment,
    Signature,
    patient_median_difference = median_difference,
    patient_cliffs_delta = cliffs_delta
  ) %>%
  left_join(
    roi_exploratory,
    by = c("Compartment", "Signature")
  ) %>%
  mutate(
    same_direction_median =
      sign(patient_median_difference) ==
      sign(ROI_median_difference),
    
    same_direction_cliff =
      sign(patient_cliffs_delta) ==
      sign(ROI_cliffs_delta)
  )

cat("\n===== ROI vs PATIENT DIRECTION CHECK =====\n")
print(as.data.frame(roi_vs_patient))

# ============================================================
# STEP 12: Patient-level paired LFR vs TTR comparison
#          All 6 patients
# ============================================================

paired_compartment <- lapply(
  target_signatures,
  function(sig) {
    
    d <- delta_df %>%
      filter(Signature == sig)
    
    wt <- wilcox.test(
      d$Delta_LFR_minus_TTR,
      mu = 0,
      alternative = "two.sided",
      exact = TRUE
    )
    
    data.frame(
      Signature = sig,
      n = nrow(d),
      median_LFR_minus_TTR =
        median(d$Delta_LFR_minus_TTR),
      min_delta =
        min(d$Delta_LFR_minus_TTR),
      max_delta =
        max(d$Delta_LFR_minus_TTR),
      n_LFR_higher =
        sum(d$Delta_LFR_minus_TTR > 0),
      n_TTR_higher =
        sum(d$Delta_LFR_minus_TTR < 0),
      p_value = wt$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR_BH = p.adjust(p_value, method = "BH")
  )

cat("\n===== PATIENT-LEVEL PAIRED LFR vs TTR =====\n")
print(as.data.frame(paired_compartment))

write.csv(
  paired_compartment,
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_paired_LFR_vs_TTR_FORMAL.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    patient_primary = patient_primary_formal,
    delta_df = delta_df,
    interaction = interaction_final,
    interaction_permutation = interaction_perm,
    paired_compartment = paired_compartment,
    roi_vs_patient = roi_vs_patient
  ),
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_revision_STEP12.rds")
)

file.exists(file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_revision_STEP12.rds"))

# ============================================================
