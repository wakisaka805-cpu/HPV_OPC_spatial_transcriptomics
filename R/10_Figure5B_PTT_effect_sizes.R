source(file.path("R", "config.R"))

# =====================================================
# Figure 5B Revision
# PTT transcriptional programs
# Case-level Cliff's delta
# =====================================================

library(dplyr)
library(ggplot2)

PTT_RESULTS_DIR <- file.path(RESULTS_ROOT, "PTT")
dir.create(PTT_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------
# Load data
# Tex-like included version
# -----------------------------------------------------

ptt <- read.csv(
  file.path(PTT_ROOT, "PTT_ssGSEA_scores_HPV_v2_TexLike.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(Signature = recode(Signature,
    "TEX_LIKE_CORE18_SHI2023" = "TEX_LIKE_PROGRAM_SHI2023"
  ))

# -----------------------------------------------------
# Final 8 signatures used in the revised manuscript
# -----------------------------------------------------

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

label_map <- c(
  "STEMLIKE_CORE16_TCF7_IL7R" = "Stem-like",
  "PRE_EXHAUSTION_TPEX18" = "Tpex",
  "T_CELL_ACTIVITY_PAN" = "T-cell activity",
  "B_CELL_ACTIVITY_PAN" = "B-cell activity",
  "TEX_LIKE_PROGRAM_SHI2023" = "Tex-like",
  "ESTROGEN_RESPONSE_EARLY" = "Estrogen early",
  "ESTROGEN_RESPONSE_LATE" = "Estrogen late",
  "BILE_ACID_METABOLISM" = "Bile acid metabolism"
)

# -----------------------------------------------------
# Cliff's delta
# Positive = higher in HPV-positive
# Negative = higher in HPV-negative
# -----------------------------------------------------

cliffs_delta <- function(x_pos, x_neg) {
  comp <- outer(x_pos, x_neg, FUN = "-")
  (sum(comp > 0) - sum(comp < 0)) / length(comp)
}

# -----------------------------------------------------
# Calculate statistics
# -----------------------------------------------------

effect_df <- ptt %>%
  filter(Signature %in% target_signatures) %>%
  mutate(Score = as.numeric(Score)) %>%
  group_by(Signature) %>%
  summarise(
    CliffsDelta = cliffs_delta(
      Score[HPV == "HPV-positive"],
      Score[HPV == "HPV-negative"]
    ),
    p_value = wilcox.test(
      Score[HPV == "HPV-positive"],
      Score[HPV == "HPV-negative"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    Signature_label = label_map[Signature]
  )

# -----------------------------------------------------
# Biological order
# -----------------------------------------------------

effect_df$Signature_label <- factor(
  effect_df$Signature_label,
  levels = rev(c(
    "Stem-like",
    "Tpex",
    "T-cell activity",
    "B-cell activity",
    "Tex-like",
    "Estrogen early",
    "Estrogen late",
    "Bile acid metabolism"
  ))
)

# -----------------------------------------------------
# Check and save values
# -----------------------------------------------------

print(effect_df)

write.csv(
  effect_df,
  file.path(PTT_RESULTS_DIR, "Fig5B_PTT_Revision_values.csv"),
  row.names = FALSE
)

# -----------------------------------------------------
# Times New Roman
# -----------------------------------------------------

# -----------------------------------------------------
# Plot
# -----------------------------------------------------

p5b <- ggplot(
  effect_df,
  aes(
    x = CliffsDelta,
    y = Signature_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_col(
    width = 0.68,
    fill = "grey45",
    color = "black",
    linewidth = 0.3
  ) +
  scale_x_continuous(
    limits = c(-0.60, 0.60),
    breaks = seq(-0.6, 0.6, by = 0.2),
    labels = c(
      "-0.6",
      "-0.4",
      "-0.2",
      "0",
      "0.2",
      "0.4",
      "0.6"
    )
  ) +
  labs(
    x = "Higher in HPV-negative        \u2190   Cliff's delta   \u2192        Higher in HPV-positive",
    y = NULL
  ) +
  theme_classic(
    base_size = 14,
    base_family = PUBLICATION_FONT
  ) +
  theme(
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    axis.title.x = element_text(
      size = 12,
      margin = margin(t = 8)
    ),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

print(p5b)