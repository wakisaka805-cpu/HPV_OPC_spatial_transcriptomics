source(file.path("R", "config.R"))

# ============================================================
# Figure 2
# TCGA-OPSCC HPV-associated transcriptional programs
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Locate formal revision-analysis files
# ------------------------------------------------------------

list.files(
  path = ".",
  pattern = "TCGA_OPSCC.*FORMAL.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# ------------------------------------------------------------
# 2. Search configured results folder for formal TCGA files
# ------------------------------------------------------------

list.files(
  path = RESULTS_ROOT,
  pattern = "TCGA_OPSCC.*FORMAL.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# ------------------------------------------------------------
# 3. Load formal unadjusted and adjusted OPSCC results
# ------------------------------------------------------------

tcga_dir <- file.path(RESULTS_ROOT, "TCGA_OPSCC77_revision")

unadj_res <- read.csv(
  file.path(
    tcga_dir,
    "TCGA_OPSCC77_unadjusted_HPV_results_FORMAL.csv"
  ),
  stringsAsFactors = FALSE
)

adj_res <- read.csv(
  file.path(
    tcga_dir,
    "TCGA_OPSCC74_adjusted_HPV_results_FORMAL.csv"
  ),
  stringsAsFactors = FALSE
)

# Check
dim(unadj_res)
colnames(unadj_res)
unadj_res

dim(adj_res)
colnames(adj_res)
adj_res

# ------------------------------------------------------------
# 4. Prepare Fig 2A
#    Unadjusted HPV effect in strict TCGA-OPSCC cohort
#    N = 77 (HPV- 26, HPV+ 51)
# ------------------------------------------------------------

fig2a_df <- unadj_res %>%
  mutate(
    Signature_label = recode(
      Signature,
      "STEM_LIKE_T_CELL_PROGRAM" = "Stem-like",
      "PRE_EXHAUSTION_PROGRAM" = "Tpex",
      "T_CELL_ACTIVITY" = "T-cell activity",
      "B_CELL_ACTIVITY" = "B-cell activity",
      "TEX_LIKE_PROGRAM" = "Tex-like",
      "HALLMARK_ESTROGEN_RESPONSE_EARLY" = "Estrogen response early",
      "HALLMARK_ESTROGEN_RESPONSE_LATE" = "Estrogen response late",
      "HALLMARK_BILE_ACID_METABOLISM" = "Bile acid metabolism"
    ),
    Signature_label = factor(
      Signature_label,
      levels = rev(c(
        "Stem-like",
        "Tpex",
        "T-cell activity",
        "B-cell activity",
        "Tex-like",
        "Estrogen response early",
        "Estrogen response late",
        "Bile acid metabolism"
      ))
    )
  )

fig2a_df %>%
  select(
    Signature_label,
    cliffs_delta,
    CI95_low,
    CI95_high,
    raw_p,
    BH_FDR
  )

# ------------------------------------------------------------
# 5. Draw Fig 2A
# ------------------------------------------------------------

fig2a <- ggplot(
  fig2a_df,
  aes(
    x = cliffs_delta,
    y = Signature_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = CI95_low,
      xmax = CI95_high
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(-0.30, 0.90),
    breaks = c(-0.25, 0, 0.25, 0.50, 0.75)
  ) +
  labs(
    x = "Cliff's delta (HPV+ vs HPV−)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig2a

# ------------------------------------------------------------
# 6. Prepare Fig 2B
#    Multivariable-adjusted HPV association
#    Complete-case N = 74
# ------------------------------------------------------------

fig2b_df <- adj_res %>%
  mutate(
    Signature_label = recode(
      Signature,
      "STEM_LIKE_T_CELL_PROGRAM" = "Stem-like",
      "PRE_EXHAUSTION_PROGRAM" = "Tpex",
      "T_CELL_ACTIVITY" = "T-cell activity",
      "B_CELL_ACTIVITY" = "B-cell activity",
      "TEX_LIKE_PROGRAM" = "Tex-like",
      "HALLMARK_ESTROGEN_RESPONSE_EARLY" = "Estrogen response early",
      "HALLMARK_ESTROGEN_RESPONSE_LATE" = "Estrogen response late",
      "HALLMARK_BILE_ACID_METABOLISM" = "Bile acid metabolism"
    ),
    Signature_label = factor(
      Signature_label,
      levels = rev(c(
        "Stem-like",
        "Tpex",
        "T-cell activity",
        "B-cell activity",
        "Tex-like",
        "Estrogen response early",
        "Estrogen response late",
        "Bile acid metabolism"
      ))
    )
  )

fig2b_df %>%
  select(
    Signature_label,
    adjusted_beta,
    CI95_low,
    CI95_high,
    raw_p,
    BH_FDR
  )

# ------------------------------------------------------------
# 7. Draw Fig 2B
# ------------------------------------------------------------

fig2b <- ggplot(
  fig2b_df,
  aes(
    x = adjusted_beta,
    y = Signature_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = CI95_low,
      xmax = CI95_high
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(-0.03, 0.24),
    breaks = c(0, 0.05, 0.10, 0.15, 0.20)
  ) +
  labs(
    x = "Adjusted HPV coefficient",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig2b

# ------------------------------------------------------------
# 8. Combine Fig 2A and Fig 2B
# ------------------------------------------------------------

library(patchwork)

fig2_main <- 
  fig2a /
  fig2b +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig2_main

# ------------------------------------------------------------
# 9. Finalize Main Figure 2
#    Times New Roman
# ------------------------------------------------------------

fig2_main_final <- fig2_main &
  theme(
    text = element_text(
      family = PUBLICATION_FONT,
      color = "black"
    ),
    axis.text = element_text(
      family = PUBLICATION_FONT,
      color = "black"
    ),
    axis.title = element_text(
      family = PUBLICATION_FONT,
      color = "black"
    ),
    plot.tag = element_text(
      family = PUBLICATION_FONT,
      face = "bold",
      size = 14
    )
  )

fig2_main_final

# ------------------------------------------------------------
# 10. Export Main Figure 2
# ------------------------------------------------------------

dir.create(
  "Revision_Figures",
  showWarnings = FALSE
)

ggsave(
  filename = "Revision_Figures/Figure2_TCGA_OPSCC_FINAL.tiff",
  plot = fig2_main_final,
  width = 7.5,
  height = 6.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = "Revision_Figures/Figure2_TCGA_OPSCC_FINAL.pdf",
  plot = fig2_main_final,
  width = 7.5,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)

file.exists(
  "Revision_Figures/Figure2_TCGA_OPSCC_FINAL.tiff"
)

file.exists(
  "Revision_Figures/Figure2_TCGA_OPSCC_FINAL.pdf"
)

# ------------------------------------------------------------
# 11. Save Figure 2 to configured revision figure folder
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
    "Figure2_TCGA_OPSCC_FINAL.tiff"
  ),
  plot = fig2_main_final,
  width = 7.5,
  height = 6.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    revision_fig_dir,
    "Figure2_TCGA_OPSCC_FINAL.pdf"
  ),
  plot = fig2_main_final,
  width = 7.5,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)

normalizePath(revision_fig_dir)