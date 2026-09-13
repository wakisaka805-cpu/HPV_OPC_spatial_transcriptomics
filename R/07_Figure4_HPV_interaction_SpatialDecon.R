source(file.path("R", "config.R"))

GEOMX_RESULTS_DIR <- file.path(RESULTS_ROOT, "GeoMx")
REVISION_FIGURE_DIR <- file.path(RESULTS_ROOT, "PLOS_ONE_Revision_Figures")
dir.create(REVISION_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Figure 4
# HPV-associated differences, HPV-by-compartment interaction,
# and SpatialDecon sensitivity analysis
# ============================================================

library(tidyverse)

# ------------------------------------------------------------
# 1. Load formal patient-level HPV comparison results
# ------------------------------------------------------------

hpv_res <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_HPV_comparison_FORMAL.csv"),
  stringsAsFactors = FALSE
)

# Check
dim(hpv_res)
colnames(hpv_res)

hpv_res

# ------------------------------------------------------------
# 2. Prepare Fig 4A
#    Patient-level HPV effect sizes within LFR and TTR
# ------------------------------------------------------------

fig4a_df <- hpv_res %>%
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
    Compartment = factor(
      Compartment,
      levels = c("LFR", "TTR")
    )
  )

fig4a_df %>%
  select(
    Compartment,
    Signature_label,
    cliffs_delta,
    Cliff_CI_lower,
    Cliff_CI_upper,
    p_value,
    FDR_BH
  )

# ------------------------------------------------------------
# 3. Draw Fig 4A
#    Cliff's delta for HPV+ vs HPV- within LFR and TTR
# ------------------------------------------------------------

fig4a <- ggplot(
  fig4a_df,
  aes(
    x = cliffs_delta,
    y = Signature_label,
    shape = Compartment
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_errorbarh(
    aes(
      xmin = Cliff_CI_lower,
      xmax = Cliff_CI_upper
    ),
    height = 0.18,
    position = position_dodge(width = 0.45),
    linewidth = 0.7
  ) +
  geom_point(
    position = position_dodge(width = 0.45),
    size = 2.8
  ) +
  scale_shape_manual(
    values = c(
      "LFR" = 16,
      "TTR" = 17
    )
  ) +
  scale_x_continuous(
    limits = c(-1.05, 1.05),
    breaks = c(-1, -0.5, 0, 0.5, 1)
  ) +
  labs(
    x = "Cliff's delta (HPV+ vs HPV−)",
    y = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

fig4a

# ------------------------------------------------------------
# 4. Load formal HPV-by-compartment interaction results
# ------------------------------------------------------------

interaction_res <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_patient_level_HPVxCompartment_interaction_FORMAL.csv"),
  stringsAsFactors = FALSE
)

dim(interaction_res)
colnames(interaction_res)

interaction_res

# ------------------------------------------------------------
# 5. Prepare and draw Fig 4B
#    HPV-by-compartment interaction
# ------------------------------------------------------------

fig4b_df <- interaction_res %>%
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
    P_label = paste0(
      "Pperm = ",
      formatC(
        exact_permutation_p,
        format = "f",
        digits = 3
      )
    )
  )

fig4b <- ggplot(
  fig4b_df,
  aes(
    x = interaction_estimate,
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
      xmin = CI_lower,
      xmax = CI_upper
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_text(
    aes(
      x = CI_upper + 0.045,
      label = P_label
    ),
    hjust = 0,
    size = 3.1
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  labs(
    x = "HPV-by-compartment interaction estimate",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig4b

# ------------------------------------------------------------
# 6. Load formal SpatialDecon HPV comparison results
# ------------------------------------------------------------

decon_res <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_HPV_comparison.csv"),
  stringsAsFactors = FALSE
)

dim(decon_res)
colnames(decon_res)

decon_res

# ------------------------------------------------------------
# 7. Load patient-level SpatialDecon estimates for Fig 4C
# ------------------------------------------------------------

decon_patient <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_patient_level.csv"),
  stringsAsFactors = FALSE
)

dim(decon_patient)
colnames(decon_patient)

decon_patient

# ------------------------------------------------------------
# 8. Prepare patient-level SpatialDecon data for Fig 4C
# ------------------------------------------------------------

fig4c_df <- decon_patient %>%
  select(
    Case,
    HPV,
    B_cell_prop,
    CD4_T_prop,
    CD8_T_prop,
    Total_T_prop
  ) %>%
  pivot_longer(
    cols = c(
      B_cell_prop,
      CD4_T_prop,
      CD8_T_prop,
      Total_T_prop
    ),
    names_to = "Cell_type",
    values_to = "Proportion"
  ) %>%
  mutate(
    Cell_type = recode(
      Cell_type,
      "B_cell_prop" = "B cells",
      "CD4_T_prop" = "CD4 T cells",
      "CD8_T_prop" = "CD8 T cells",
      "Total_T_prop" = "Total T cells"
    ),
    Cell_type = factor(
      Cell_type,
      levels = c(
        "B cells",
        "CD4 T cells",
        "CD8 T cells",
        "Total T cells"
      )
    ),
    HPV = factor(
      HPV,
      levels = c("negative", "positive"),
      labels = c("HPV−", "HPV+")
    )
  )

# Check
fig4c_df

# ------------------------------------------------------------
# 9. Draw Fig 4C
#    Patient-level broad cell-composition estimates in LFRs
# ------------------------------------------------------------

fig4c <- ggplot(
  fig4c_df,
  aes(
    x = HPV,
    y = Proportion,
    shape = HPV
  )
) +
  geom_point(
    size = 3,
    position = position_jitter(
      width = 0.06,
      height = 0
    )
  ) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    linewidth = 0.6
  ) +
  facet_wrap(
    ~ Cell_type,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_shape_manual(
    values = c(
      "HPV−" = 16,
      "HPV+" = 17
    )
  ) +
  labs(
    x = NULL,
    y = "Estimated proportion"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )

fig4c

# ------------------------------------------------------------
# 10. Combine Fig 4A, Fig 4B, and Fig 4C
# ------------------------------------------------------------

library(patchwork)

fig4_combined <-
  (fig4a | fig4b) /
  fig4c +
  plot_layout(
    heights = c(1.25, 0.85),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig4_combined

# ------------------------------------------------------------
# 11. Clean Fig 4B for the final multi-panel figure
# ------------------------------------------------------------

fig4b_clean <- ggplot(
  fig4b_df,
  aes(
    x = interaction_estimate,
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
      xmin = CI_lower,
      xmax = CI_upper
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_text(
    data = fig4b_df %>%
      filter(Signature == "ESTROGEN_RESPONSE_EARLY"),
    aes(
      x = CI_upper + 0.035,
      label = "Permutation P = 0.067"
    ),
    hjust = 0,
    size = 3
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.08, 0.35))
  ) +
  labs(
    x = "HPV-by-compartment interaction estimate",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig4b_clean

# ------------------------------------------------------------
# 12. Recombine Figure 4 with improved layout
# ------------------------------------------------------------

fig4_main <- (
  fig4a + fig4b_clean +
    plot_layout(widths = c(1, 1))
) /
  fig4c +
  plot_layout(
    heights = c(1.25, 0.90),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig4_main

# ------------------------------------------------------------
# 13. Main Figure 4
#     A: HPV effect within compartments
#     B: HPV-by-compartment interaction
# ------------------------------------------------------------

fig4_main_clean <-
  fig4a /
  fig4b_clean +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig4_main_clean

# ------------------------------------------------------------
# 14. Final refinement of Fig 4B
# ------------------------------------------------------------

fig4b_final <- ggplot(
  fig4b_df,
  aes(
    x = interaction_estimate,
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
      xmin = CI_lower,
      xmax = CI_upper
    ),
    height = 0.18,
    linewidth = 0.7
  ) +
  geom_point(
    size = 2.8
  ) +
  geom_text(
    data = fig4b_df %>%
      filter(Signature == "ESTROGEN_RESPONSE_EARLY"),
    aes(
      x = 0.20,
      label = "Permutation P = 0.067"
    ),
    hjust = 0,
    size = 3
  ) +
  coord_cartesian(
    xlim = c(-0.52, 0.62),
    clip = "off"
  ) +
  labs(
    x = "HPV-by-compartment interaction estimate",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black")
  )

fig4_main_final <-
  fig4a /
  fig4b_final +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

fig4_main_final

# ------------------------------------------------------------
# 15. Finalize and save Main Figure 4
# ------------------------------------------------------------

fig4_main_final_font <- fig4_main_final &
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
fig4_main_final_font

# ------------------------------------------------------------
# 16. Export Main Figure 4
# ------------------------------------------------------------

dir.create(REVISION_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "Figure4_HPV_effect_interaction_FINAL.tiff"),
  plot = fig4_main_final_font,
  width = 7.5,
  height = 6.0,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "Figure4_HPV_effect_interaction_FINAL.pdf"),
  plot = fig4_main_final_font,
  width = 7.5,
  height = 6.0,
  units = "in",
  device = cairo_pdf
)

file.exists(
  file.path(REVISION_FIGURE_DIR, "Figure4_HPV_effect_interaction_FINAL.tiff")
)

file.exists(
  file.path(REVISION_FIGURE_DIR, "Figure4_HPV_effect_interaction_FINAL.pdf")
)

# ------------------------------------------------------------
# 17. Load formal Stem-like vs SpatialDecon correlation results
# ------------------------------------------------------------

decon_cor <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_SpatialDecon_StemLike_correlations_FORMAL.csv"),
  stringsAsFactors = FALSE
)

dim(decon_cor)
colnames(decon_cor)

decon_cor

# ------------------------------------------------------------
# 18. Load patient-level Stem-like + SpatialDecon data
# ------------------------------------------------------------

stem_decon_patient <- read.csv(
  file.path(GEOMX_RESULTS_DIR, "GeoMx_LFR_StemLike_SpatialDecon_patient_level.csv"),
  stringsAsFactors = FALSE
)

dim(stem_decon_patient)
colnames(stem_decon_patient)

stem_decon_patient

# ------------------------------------------------------------
# 19. Supplementary Fig Sx-B
#     Stem-like score vs CD8 T-cell proportion
# ------------------------------------------------------------

sup_b <- ggplot(
  stem_decon_patient,
  aes(
    x = CD8_T_prop,
    y = Stem_like_score,
    shape = HPV
  )
) +
  geom_point(size = 3.2) +
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
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "\u03c1 = -0.371, P = 0.468",
    hjust = 1.1,
    vjust = 1.4,
    size = 3.3
  ) +
  labs(
    x = "Estimated CD8 T-cell proportion",
    y = "Stem-like ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

sup_b

# ------------------------------------------------------------
# 20. Supplementary Fig Sx-C
#     Stem-like score vs total T-cell proportion
# ------------------------------------------------------------

sup_c <- ggplot(
  stem_decon_patient,
  aes(
    x = Total_T_prop,
    y = Stem_like_score,
    shape = HPV
  )
) +
  geom_point(size = 3.2) +
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
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "\u03c1 = 0.829, P = 0.042",
    hjust = 1.1,
    vjust = 1.4,
    size = 3.3
  ) +
  labs(
    x = "Estimated total T-cell proportion",
    y = "Stem-like ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

sup_c

# ------------------------------------------------------------
# 21. Combine SpatialDecon Supplementary Figure
# ------------------------------------------------------------

sup_spatialdecon <- 
  fig4c /
  (sup_b | sup_c) +
  plot_layout(
    heights = c(1, 1)
  ) +
  plot_annotation(
    tag_levels = "A"
  )

sup_spatialdecon

# ------------------------------------------------------------
# 22. Final refinement of Supplementary panels B and C
# ------------------------------------------------------------

sup_b_final <- ggplot(
  stem_decon_patient,
  aes(
    x = CD8_T_prop,
    y = Stem_like_score,
    shape = HPV
  )
) +
  geom_point(size = 3.2) +
  scale_shape_manual(
    values = c("negative" = 16, "positive" = 17),
    labels = c("negative" = "HPV−", "positive" = "HPV+")
  ) +
  annotate(
    "text",
    x = 0.043,
    y = 0.451,
    label = "\u03c1 = -0.371, P = 0.468",
    hjust = 0,
    size = 3.2
  ) +
  labs(
    x = "Estimated CD8 T-cell proportion",
    y = "Stem-like ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )

sup_c_final <- ggplot(
  stem_decon_patient,
  aes(
    x = Total_T_prop,
    y = Stem_like_score,
    shape = HPV
  )
) +
  geom_point(size = 3.2) +
  scale_shape_manual(
    values = c("negative" = 16, "positive" = 17),
    labels = c("negative" = "HPV−", "positive" = "HPV+")
  ) +
  annotate(
    "text",
    x = 0.30,
    y = 0.451,
    label = "\u03c1 = 0.829, P = 0.042",
    hjust = 0,
    size = 3.2
  ) +
  labs(
    x = "Estimated total T-cell proportion",
    y = "Stem-like ssGSEA score",
    shape = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right"
  )

sup_spatialdecon_final <-
  fig4c /
  (sup_b_final | sup_c_final) +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A"
  )

sup_spatialdecon_final

# ------------------------------------------------------------
# 23. Finalize SpatialDecon Supplementary Figure
#     Times New Roman
# ------------------------------------------------------------

sup_spatialdecon_final_font <- sup_spatialdecon_final &
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

sup_spatialdecon_final_font

# ------------------------------------------------------------
# 24. Export SpatialDecon Supplementary Figure
# ------------------------------------------------------------

dir.create(REVISION_FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "FigureS_SpatialDecon_FINAL.tiff"),
  plot = sup_spatialdecon_final_font,
  width = 7.5,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)

ggsave(
  filename = file.path(REVISION_FIGURE_DIR, "FigureS_SpatialDecon_FINAL.pdf"),
  plot = sup_spatialdecon_final_font,
  width = 7.5,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

file.exists(
  file.path(REVISION_FIGURE_DIR, "FigureS_SpatialDecon_FINAL.tiff")
)

file.exists(
  file.path(REVISION_FIGURE_DIR, "FigureS_SpatialDecon_FINAL.pdf")
)