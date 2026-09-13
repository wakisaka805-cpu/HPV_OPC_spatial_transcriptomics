source(file.path("R", "config.R"))

# ============================================================
# Supplementary Figure
# Exploratory overall-survival analysis in TCGA-OPSCC
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Set formal TCGA revision directory
# ------------------------------------------------------------

tcga_dir <- file.path(RESULTS_ROOT, "TCGA_OPSCC77_revision")

# ------------------------------------------------------------
# 2. Load formal Cox results
# ------------------------------------------------------------

cox_uni <- read.csv(
  file.path(
    tcga_dir,
    "TCGA_OPSCC77_Cox_univariate_perSD_FORMAL.csv"
  ),
  stringsAsFactors = FALSE
)

cox_hpv <- read.csv(
  file.path(
    tcga_dir,
    "TCGA_OPSCC77_Cox_HPV_adjusted_perSD_FORMAL.csv"
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 3. Check
# ------------------------------------------------------------

dim(cox_uni)
colnames(cox_uni)
cox_uni

dim(cox_hpv)
colnames(cox_hpv)
cox_hpv

# ------------------------------------------------------------
# 4. Load HPV-stratified Cox sensitivity analysis
# ------------------------------------------------------------

cox_strat <- read.csv(
  file.path(
    tcga_dir,
    "TCGA_OPSCC77_Cox_HPV_stratified_sensitivity.csv"
  ),
  stringsAsFactors = FALSE
)

dim(cox_strat)
colnames(cox_strat)

cox_strat

# ------------------------------------------------------------
# 5. Prepare data for survival forest plots
# ------------------------------------------------------------

signature_levels <- rev(c(
  "Stem-like",
  "Tpex",
  "B-cell activity",
  "T-cell activity",
  "Tex-like",
  "Estrogen response early",
  "Estrogen response late",
  "Bile acid metabolism"
))

recode_signature <- function(x) {
  recode(
    x,
    "STEM_LIKE_T_CELL_PROGRAM" = "Stem-like",
    "PRE_EXHAUSTION_PROGRAM" = "Tpex",
    "T_CELL_ACTIVITY" = "T-cell activity",
    "B_CELL_ACTIVITY" = "B-cell activity",
    "TEX_LIKE_PROGRAM" = "Tex-like",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY" = "Estrogen response early",
    "HALLMARK_ESTROGEN_RESPONSE_LATE" = "Estrogen response late",
    "HALLMARK_BILE_ACID_METABOLISM" = "Bile acid metabolism"
  )
}

surv_a <- cox_uni %>%
  transmute(
    Signature_label = recode_signature(Signature),
    HR = HR_per_SD,
    CI_low = CI95_low,
    CI_high = CI95_high,
    P = raw_p,
    FDR = BH_FDR
  )

surv_b <- cox_hpv %>%
  transmute(
    Signature_label = recode_signature(Signature),
    HR = HR_signature_per_SD,
    CI_low = CI95_low_signature,
    CI_high = CI95_high_signature,
    P = raw_p_signature,
    FDR = BH_FDR_signature
  )

surv_c <- cox_strat %>%
  transmute(
    Signature_label = recode_signature(Signature),
    HR = HR_signature_per_SD,
    CI_low = CI95_low,
    CI_high = CI95_high,
    P = raw_p,
    FDR = BH_FDR
  )

surv_a$Signature_label <- factor(
  surv_a$Signature_label,
  levels = signature_levels
)

surv_b$Signature_label <- factor(
  surv_b$Signature_label,
  levels = signature_levels
)

surv_c$Signature_label <- factor(
  surv_c$Signature_label,
  levels = signature_levels
)

# Check
surv_a
surv_b
surv_c

# ------------------------------------------------------------
# 6. Create survival forest plots
# ------------------------------------------------------------

library(patchwork)

make_forest <- function(dat, panel_title) {
  
  ggplot(
    dat,
    aes(
      x = HR,
      y = Signature_label
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_errorbarh(
      aes(
        xmin = CI_low,
        xmax = CI_high
      ),
      height = 0.18,
      linewidth = 0.6
    ) +
    geom_point(
      size = 2.5
    ) +
    scale_x_log10(
      limits = c(0.25, 2.0),
      breaks = c(0.25, 0.5, 1, 2)
    ) +
    labs(
      title = panel_title,
      x = "Hazard ratio per 1-SD increase",
      y = NULL
    ) +
    theme_classic(
      base_family = PUBLICATION_FONT,
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        family = PUBLICATION_FONT,
        face = "bold",
        size = 12,
        hjust = 0
      ),
      axis.text = element_text(
        family = PUBLICATION_FONT,
        color = "black"
      ),
      axis.title = element_text(
        family = PUBLICATION_FONT,
        color = "black"
      )
    )
}

p_surv_a <- make_forest(
  surv_a,
  "A  Univariable Cox model"
)

p_surv_b <- make_forest(
  surv_b,
  "B  HPV-adjusted Cox model"
)

p_surv_c <- make_forest(
  surv_c,
  "C  HPV-stratified Cox sensitivity analysis"
)

fig_survival <- (
  p_surv_a /
    p_surv_b /
    p_surv_c
) +
  plot_layout(
    heights = c(1, 1, 1)
  )

fig_survival

# ------------------------------------------------------------
# 7. Save final survival Supplementary Figure
# ------------------------------------------------------------

revision_fig_dir <- file.path(RESULTS_ROOT, "PLOS_ONE_Revision_Figures")

dir.create(
  revision_fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    revision_fig_dir,
    "FigureS_TCGA_OS_FINAL.tiff"
  ),
  plot = fig_survival,
  width = 7.5,
  height = 6.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    revision_fig_dir,
    "FigureS_TCGA_OS_FINAL.pdf"
  ),
  plot = fig_survival,
  width = 7.5,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)

file.exists(
  file.path(revision_fig_dir, "FigureS_TCGA_OS_FINAL.tiff")
)

file.exists(
  file.path(revision_fig_dir, "FigureS_TCGA_OS_FINAL.pdf")
)