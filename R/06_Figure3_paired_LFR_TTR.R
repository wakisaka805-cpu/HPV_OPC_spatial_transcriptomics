source(file.path("R", "config.R"))

GEOMX_RESULTS_DIR <- file.path(RESULTS_ROOT, "GeoMx")
REVISION_FIGURE_DIR <- file.path(RESULTS_ROOT, "PLOS_ONE_Revision_Figures")
dir.create(REVISION_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Fig 3A
# Within-patient paired comparison:
# immune transcriptional programs in LFR vs TTR
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Load case-level ssGSEA scores
# ------------------------------------------------------------

case_df <- read.csv(
  file.path(GEOMX_ROOT, "GeoMx_ssGSEA_case_level_scores.csv"),
  stringsAsFactors = FALSE
)

# Check
dim(case_df)
colnames(case_df)
unique(case_df$Signature)

# ------------------------------------------------------------
# 2. Select the five immune programs for Fig 3A
# ------------------------------------------------------------

immune_signatures <- c(
  "STEMLIKE_CORE16_TCF7_IL7R",
  "PRE_EXHAUSTION_TPEX18",
  "T_CELL_ACTIVITY_PAN",
  "B_CELL_ACTIVITY_PAN",
  "TEX_LIKE_PROGRAM_SHI2023"
)

fig3a_df <- case_df %>%
  filter(Signature %in% immune_signatures) %>%
  mutate(
    Signature_label = recode(
      Signature,
      "STEMLIKE_CORE16_TCF7_IL7R" = "Stem-like",
      "PRE_EXHAUSTION_TPEX18" = "Tpex",
      "T_CELL_ACTIVITY_PAN" = "T-cell activity",
      "B_CELL_ACTIVITY_PAN" = "B-cell activity",
      "TEX_LIKE_PROGRAM_SHI2023" = "Tex-like"
    ),
    Signature_label = factor(
      Signature_label,
      levels = c(
        "Stem-like",
        "Tpex",
        "T-cell activity",
        "B-cell activity",
        "Tex-like"
      )
    ),
    Compartment = factor(
      Compartment,
      levels = c("LFR", "TTR")
    )
  )

# Check
dim(fig3a_df)

fig3a_df %>%
  count(Signature_label, Compartment)

fig3a_df %>%
  arrange(Signature_label, Case, Compartment) %>%
  select(Case, HPV, Signature_label, Compartment, Score)

# ------------------------------------------------------------
# 3. Draw Fig 3A
#    Paired LFR vs TTR comparison for 5 immune programs
# ------------------------------------------------------------

fig3a <- ggplot(
  fig3a_df,
  aes(
    x = Compartment,
    y = Score,
    group = Case
  )
) +
  geom_line(
    aes(linetype = HPV),
    linewidth = 0.7,
    alpha = 0.8
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank(),
    legend.position = "right"
  )

fig3a

# ------------------------------------------------------------
# 4. Final styling of Fig 3A
# ------------------------------------------------------------

fig3a_final <- ggplot(
  fig3a_df,
  aes(
    x = Compartment,
    y = Score,
    group = Case
  )
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.70
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      "negative" = 16,
      "positive" = 17
    ),
    labels = c(
      "negative" = "HPV−",
      "positive" = "HPV+"
    )
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

fig3a_final

# ------------------------------------------------------------
# 5. Add FDR annotation to Fig 3A
# ------------------------------------------------------------

fig3a_ann <- fig3a_df %>%
  group_by(Signature_label) %>%
  summarise(
    y = max(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    label = "FDR = 0.036"
  )

fig3a_final <- ggplot(
  fig3a_df,
  aes(
    x = Compartment,
    y = Score,
    group = Case
  )
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.70
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  geom_text(
    data = fig3a_ann,
    aes(
      x = 1.5,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 3.2
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      "negative" = 16,
      "positive" = 17
    ),
    labels = c(
      "negative" = "HPV−",
      "positive" = "HPV+"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

fig3a_final

# ------------------------------------------------------------
# 6. Prepare Fig 3B
#    Estrogen-response and bile-acid programs
# ------------------------------------------------------------

metabolic_signatures <- c(
  "ESTROGEN_RESPONSE_EARLY",
  "ESTROGEN_RESPONSE_LATE",
  "BILE_ACID_METABOLISM"
)

fig3b_df <- case_df %>%
  filter(Signature %in% metabolic_signatures) %>%
  mutate(
    Signature_label = recode(
      Signature,
      "ESTROGEN_RESPONSE_EARLY" = "Estrogen response early",
      "ESTROGEN_RESPONSE_LATE"  = "Estrogen response late",
      "BILE_ACID_METABOLISM"    = "Bile acid metabolism"
    ),
    Signature_label = factor(
      Signature_label,
      levels = c(
        "Estrogen response early",
        "Estrogen response late",
        "Bile acid metabolism"
      )
    ),
    Compartment = factor(
      Compartment,
      levels = c("LFR", "TTR")
    )
  )

# Check
dim(fig3b_df)

fig3b_df %>%
  count(Signature_label, Compartment)

fig3b_df %>%
  arrange(Signature_label, Case, Compartment) %>%
  select(Case, HPV, Signature_label, Compartment, Score)

# ------------------------------------------------------------
# 7. Draw Fig 3B
#    Estrogen-response and bile-acid programs
# ------------------------------------------------------------

fig3b_ann <- tibble(
  Signature_label = factor(
    c(
      "Estrogen response early",
      "Estrogen response late",
      "Bile acid metabolism"
    ),
    levels = c(
      "Estrogen response early",
      "Estrogen response late",
      "Bile acid metabolism"
    )
  ),
  label = c(
    "FDR = 0.036",
    "FDR = 0.036",
    "FDR = 0.438"
  )
)

fig3b_ann <- fig3b_ann %>%
  left_join(
    fig3b_df %>%
      group_by(Signature_label) %>%
      summarise(
        y = max(Score, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Signature_label"
  )

fig3b_final <- ggplot(
  fig3b_df,
  aes(
    x = Compartment,
    y = Score,
    group = Case
  )
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.70
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  geom_text(
    data = fig3b_ann,
    aes(
      x = 1.5,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = -1.0,
    size = 3.2
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      "negative" = 16,
      "positive" = 17
    ),
    labels = c(
      "negative" = "HPV−",
      "positive" = "HPV+"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

fig3b_final

# ------------------------------------------------------------
# 8. Load formal paired-analysis results for Fig 3C
# ------------------------------------------------------------

paired_res <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_paired_LFR_vs_TTR_FORMAL.csv"),
  stringsAsFactors = FALSE
)

# Check structure
dim(paired_res)
colnames(paired_res)

# Display full table
paired_res

# ------------------------------------------------------------
# 9. Prepare Fig 3C
#    Summary of within-patient LFR - TTR differences
# ------------------------------------------------------------

fig3c_df <- paired_res %>%
  mutate(
    Signature_label = recode(
      Signature,
      "STEMLIKE_CORE16_TCF7_IL7R" = "Stem-like",
      "PRE_EXHAUSTION_TPEX18" = "Tpex",
      "T_CELL_ACTIVITY_PAN" = "T-cell activity",
      "B_CELL_ACTIVITY_PAN" = "B-cell activity",
      "TEX_LIKE_PROGRAM_SHI2023" = "Tex-like",
      "ESTROGEN_RESPONSE_EARLY" = "Estrogen response early",
      "ESTROGEN_RESPONSE_LATE" = "Estrogen response late",
      "BILE_ACID_METABOLISM" = "Bile acid metabolism"
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
    ),
    FDR_label = ifelse(
      FDR_BH < 0.05,
      "FDR = 0.036",
      "FDR = 0.438"
    )
  )

fig3c_df %>%
  select(
    Signature_label,
    median_LFR_minus_TTR,
    min_delta,
    max_delta,
    FDR_BH
  )

# ------------------------------------------------------------
# 10. Draw Fig 3C
#     Summary of within-patient LFR - TTR differences
# ------------------------------------------------------------

fig3c <- ggplot(
  fig3c_df,
  aes(
    x = median_LFR_minus_TTR,
    y = Signature_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_segment(
    aes(
      x = min_delta,
      xend = max_delta,
      yend = Signature_label
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      x = max_delta + 0.035,
      label = FDR_label
    ),
    hjust = 0,
    size = 3.2
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.20))
  ) +
  labs(
    x = "Within-patient difference in ssGSEA score (LFR − TTR)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig3c

# ------------------------------------------------------------
# 11. Combine Fig 3A, 3B, and 3C
# ------------------------------------------------------------

library(patchwork)

fig3_combined <-
  fig3a_final /
  fig3b_final /
  fig3c +
  plot_layout(
    heights = c(1.0, 1.0, 1.15)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig3_combined

# ------------------------------------------------------------
# 12. Main Figure 3 candidate
#     A: immune programs
#     B: estrogen / bile-acid programs
#     No statistical text inside panels
# ------------------------------------------------------------

fig3a_clean <- fig3a_final +
  theme(
    legend.position = "none"
  )

# Remove annotation layer by rebuilding from the pre-annotation plot
fig3a_clean <- ggplot(
  fig3a_df,
  aes(x = Compartment, y = Score, group = Case)
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.70
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c("negative" = 16, "positive" = 17)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )

fig3b_clean <- ggplot(
  fig3b_df,
  aes(x = Compartment, y = Score, group = Case)
) +
  geom_line(
    linewidth = 0.65,
    alpha = 0.70
  ) +
  geom_point(
    aes(shape = HPV),
    size = 2.8
  ) +
  facet_wrap(
    ~ Signature_label,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c("negative" = 16, "positive" = 17),
    labels = c("negative" = "HPV−", "positive" = "HPV+")
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    x = NULL,
    y = "ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

fig3_main <-
  fig3a_clean /
  fig3b_clean +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig3_main

# ------------------------------------------------------------
# 13. Finalize Figure 3
#     Times New Roman + publication-quality export
# ------------------------------------------------------------

# Register Times New Roman on Windows
# Apply Times New Roman to all panels
fig3_main_final <- fig3_main &
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
    strip.text = element_text(
      family = PUBLICATION_FONT,
      face = "bold",
      color = "black"
    ),
    legend.text = element_text(
      family = PUBLICATION_FONT,
      color = "black"
    ),
    plot.tag = element_text(
      family = PUBLICATION_FONT,
      face = "bold",
      size = 14
    )
  )

# Preview
fig3_main_final

# ------------------------------------------------------------
# 14. Save final Figure 3
# ------------------------------------------------------------

dir.create(REVISION_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# High-resolution TIFF
ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "Figure3_patient_level_LFR_TTR_FINAL.tiff"),
  plot = fig3_main_final,
  width = 7.5,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

# PDF for checking / vector archive
ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "Figure3_patient_level_LFR_TTR_FINAL.pdf"),
  plot = fig3_main_final,
  width = 7.5,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

# Also save the summary panel C separately
ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "FigureS_summary_LFR_TTR_candidate.pdf"),
  plot = fig3c,
  width = 7.5,
  height = 4.5,
  units = "in",
  device = cairo_pdf
)
