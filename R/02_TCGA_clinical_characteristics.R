source(file.path("R", "config.R"))

# ============================================================
# Reviewer #1 Comment 9
# TCGA-OPSCC clinicopathological characteristics by HPV status
# n = 77 (HPV-positive 51, HPV-negative 26)
# ============================================================

# ------------------------------------------------------------
# 0. Check input data
# ------------------------------------------------------------

required_vars <- c(
  "HPV_status",
  "age_years",
  "gender",
  "smoking_group",
  "clinical_stage_group",
  "subsite_group"
)

stopifnot(exists("opscc_clin_expr"))

missing_vars <- setdiff(required_vars, colnames(opscc_clin_expr))

if (length(missing_vars) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_vars, collapse = ", ")
  )
}

dat <- opscc_clin_expr[, required_vars]

cat("============================================================\n")
cat("TCGA-OPSCC clinicopathological comparison by HPV status\n")
cat("============================================================\n\n")

cat("Total N =", nrow(dat), "\n\n")
print(table(dat$HPV_status, useNA = "ifany"))


# ------------------------------------------------------------
# 1. Ensure HPV order
# ------------------------------------------------------------

dat$HPV_status <- factor(
  dat$HPV_status,
  levels = c("HPV_negative", "HPV_positive")
)


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

fmt_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p < 0.001) {
    return(formatC(p, format = "e", digits = 2))
  }
  sprintf("%.3f", p)
}

fmt_n_pct <- function(n, denom) {
  sprintf("%d (%.1f%%)", n, 100 * n / denom)
}

fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  
  q <- quantile(
    x,
    probs = c(0.25, 0.50, 0.75),
    names = FALSE
  )
  
  sprintf(
    "%.1f [%.1f–%.1f]",
    q[2], q[1], q[3]
  )
}


# ------------------------------------------------------------
# 3. AGE
# median [IQR], Wilcoxon rank-sum test
# ------------------------------------------------------------

age_neg <- dat$age_years[
  dat$HPV_status == "HPV_negative"
]

age_pos <- dat$age_years[
  dat$HPV_status == "HPV_positive"
]

age_test <- wilcox.test(
  age_years ~ HPV_status,
  data = dat,
  exact = FALSE
)

cat("\n============================================================\n")
cat("AGE\n")
cat("============================================================\n")

cat(
  "HPV-negative:",
  fmt_median_iqr(age_neg),
  " n =", sum(!is.na(age_neg)), "\n"
)

cat(
  "HPV-positive:",
  fmt_median_iqr(age_pos),
  " n =", sum(!is.na(age_pos)), "\n"
)

cat(
  "Wilcoxon P =",
  format(age_test$p.value, digits = 10),
  "\n"
)


# ------------------------------------------------------------
# 4. Function for categorical variables
# Fisher's exact test, excluding missing values from test
# Percentages use non-missing denominator within each HPV group
# ------------------------------------------------------------

analyze_categorical <- function(data, variable, label) {
  
  x <- data[[variable]]
  
  keep <- !is.na(x) & !is.na(data$HPV_status)
  
  tab <- table(
    x[keep],
    data$HPV_status[keep]
  )
  
  test <- fisher.test(tab)
  
  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  
  print(tab)
  
  cat("\nNon-missing denominators:\n")
  print(colSums(tab))
  
  cat("\nPercentages within HPV group:\n")
  print(round(prop.table(tab, margin = 2) * 100, 1))
  
  cat(
    "\nFisher exact P =",
    format(test$p.value, digits = 10),
    "\n"
  )
  
  missing_tab <- table(
    data$HPV_status,
    is.na(data[[variable]])
  )
  
  cat("\nMissingness:\n")
  print(missing_tab)
  
  invisible(
    list(
      table = tab,
      p = test$p.value,
      missing = missing_tab
    )
  )
}


# ------------------------------------------------------------
# 5. SEX / GENDER
# ------------------------------------------------------------

res_gender <- analyze_categorical(
  dat,
  "gender",
  "SEX / GENDER"
)


# ------------------------------------------------------------
# 6. SMOKING
# ------------------------------------------------------------

res_smoking <- analyze_categorical(
  dat,
  "smoking_group",
  "SMOKING STATUS"
)


# ------------------------------------------------------------
# 7. CLINICAL STAGE
# ------------------------------------------------------------

res_stage <- analyze_categorical(
  dat,
  "clinical_stage_group",
  "CLINICAL STAGE"
)


# ------------------------------------------------------------
# 8. ANATOMICAL SUBSITE
# ------------------------------------------------------------

res_subsite <- analyze_categorical(
  dat,
  "subsite_group",
  "ANATOMICAL SUBSITE"
)


# ------------------------------------------------------------
# 9. Construct publication-ready S7 Table
# ------------------------------------------------------------

get_cat_rows <- function(data, variable, display_variable) {
  
  x <- data[[variable]]
  
  levs <- levels(factor(x))
  
  neg_denom <- sum(
    data$HPV_status == "HPV_negative" & !is.na(x)
  )
  
  pos_denom <- sum(
    data$HPV_status == "HPV_positive" & !is.na(x)
  )
  
  tab <- table(
    factor(x, levels = levs),
    data$HPV_status
  )
  
  p <- fisher.test(
    table(
      x[!is.na(x)],
      data$HPV_status[!is.na(x)]
    )
  )$p.value
  
  out <- data.frame(
    Characteristic = c(
      display_variable,
      paste0("  ", levs)
    ),
    HPV_negative = c(
      "",
      sapply(
        levs,
        function(z) {
          fmt_n_pct(
            sum(
              data$HPV_status == "HPV_negative" &
                x == z,
              na.rm = TRUE
            ),
            neg_denom
          )
        }
      )
    ),
    HPV_positive = c(
      "",
      sapply(
        levs,
        function(z) {
          fmt_n_pct(
            sum(
              data$HPV_status == "HPV_positive" &
                x == z,
              na.rm = TRUE
            ),
            pos_denom
          )
        }
      )
    ),
    P_value = c(
      fmt_p(p),
      rep("", length(levs))
    ),
    stringsAsFactors = FALSE
  )
  
  out
}


# Age row
age_row <- data.frame(
  Characteristic = "Age, years, median [IQR]",
  HPV_negative = fmt_median_iqr(age_neg),
  HPV_positive = fmt_median_iqr(age_pos),
  P_value = fmt_p(age_test$p.value),
  stringsAsFactors = FALSE
)

# Categorical rows
gender_rows <- get_cat_rows(
  dat,
  "gender",
  "Sex/gender, n (%)"
)

smoking_rows <- get_cat_rows(
  dat,
  "smoking_group",
  "Smoking status, n (%)"
)

stage_rows <- get_cat_rows(
  dat,
  "clinical_stage_group",
  "Clinical stage, n (%)"
)

subsite_rows <- get_cat_rows(
  dat,
  "subsite_group",
  "Anatomical subsite, n (%)"
)


# ------------------------------------------------------------
# 10. Add missing-value rows where applicable
# ------------------------------------------------------------

make_missing_row <- function(data, variable) {
  
  n_neg <- sum(
    data$HPV_status == "HPV_negative" &
      is.na(data[[variable]])
  )
  
  n_pos <- sum(
    data$HPV_status == "HPV_positive" &
      is.na(data[[variable]])
  )
  
  data.frame(
    Characteristic = "  Missing",
    HPV_negative = as.character(n_neg),
    HPV_positive = as.character(n_pos),
    P_value = "",
    stringsAsFactors = FALSE
  )
}

if (any(is.na(dat$smoking_group))) {
  smoking_rows <- rbind(
    smoking_rows,
    make_missing_row(dat, "smoking_group")
  )
}

if (any(is.na(dat$clinical_stage_group))) {
  stage_rows <- rbind(
    stage_rows,
    make_missing_row(dat, "clinical_stage_group")
  )
}

if (any(is.na(dat$age_years))) {
  age_missing <- make_missing_row(dat, "age_years")
} else {
  age_missing <- NULL
}


# ------------------------------------------------------------
# 11. Final S7 table
# ------------------------------------------------------------

S7_table <- rbind(
  age_row,
  age_missing,
  gender_rows,
  smoking_rows,
  stage_rows,
  subsite_rows
)

colnames(S7_table) <- c(
  "Characteristic",
  "HPV-negative (n=26)",
  "HPV-positive (n=51)",
  "P value"
)

cat("\n\n")
cat("============================================================\n")
cat("FINAL S7 TABLE\n")
cat("============================================================\n\n")

print(
  S7_table,
  row.names = FALSE,
  right = FALSE
)


# ------------------------------------------------------------
# 12. Save CSV
# ------------------------------------------------------------

write.csv(
  S7_table,
  file = "S7_TCGA_OPSCC_clinicopathological_characteristics.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nSaved:\n",
  "S7_TCGA_OPSCC_clinicopathological_characteristics.csv\n"
)


# ------------------------------------------------------------
# 13. Raw exact P values for manuscript / Response letter
# ------------------------------------------------------------

p_summary <- data.frame(
  Variable = c(
    "Age",
    "Sex/gender",
    "Smoking status",
    "Clinical stage",
    "Anatomical subsite"
  ),
  Test = c(
    "Wilcoxon rank-sum",
    "Fisher's exact",
    "Fisher's exact",
    "Fisher's exact",
    "Fisher's exact"
  ),
  P_value = c(
    age_test$p.value,
    res_gender$p,
    res_smoking$p,
    res_stage$p,
    res_subsite$p
  )
)

cat("\n============================================================\n")
cat("EXACT P VALUES\n")
cat("============================================================\n\n")

print(
  p_summary,
  row.names = FALSE,
  digits = 10
)