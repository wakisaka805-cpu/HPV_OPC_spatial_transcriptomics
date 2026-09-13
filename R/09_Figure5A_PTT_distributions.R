source(file.path("R", "config.R"))

# =====================================================
# Figure 5A Revision
# PTT: Stem-like and Bile acid metabolism
# HPV-positive vs HPV-negative
# =====================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------
# Load final PTT ssGSEA data
# -----------------------------------------------------

ptt <- read.csv(
  file.path(PTT_ROOT, "PTT_ssGSEA_scores_HPV_v2_TexLike.csv"),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------
# Prepare data
# -----------------------------------------------------

plot_df <- ptt %>%
  filter(
    Signature %in% c(
      "STEMLIKE_CORE16_TCF7_IL7R",
      "BILE_ACID_METABOLISM"
    )
  ) %>%
  mutate(
    Score = as.numeric(Score),
    HPV = factor(
      HPV,
      levels = c("HPV-negative", "HPV-positive")
    ),
    Program = case_when(
      Signature == "STEMLIKE_CORE16_TCF7_IL7R" ~ "Stem-like",
      Signature == "BILE_ACID_METABOLISM" ~ "Bile acid metabolism"
    )
  ) %>%
  filter(!is.na(Score))

# -----------------------------------------------------
# Statistics
# -----------------------------------------------------

stat_df <- plot_df %>%
  group_by(Program) %>%
  summarise(
    n_neg = sum(HPV == "HPV-negative"),
    n_pos = sum(HPV == "HPV-positive"),
    p_value = wilcox.test(
      Score[HPV == "HPV-positive"],
      Score[HPV == "HPV-negative"],
      exact = FALSE
    )$p.value,
    .groups = "drop"
  )

print(stat_df)

# -----------------------------------------------------
# Font
# -----------------------------------------------------

# -----------------------------------------------------
# Function to create each panel
# -----------------------------------------------------

make_panel <- function(data, title_text) {
  
  ymax <- max(data$Score, na.rm = TRUE)
  ymin <- min(data$Score, na.rm = TRUE)
  yrange <- ymax - ymin
  
  ggplot(
    data,
    aes(x = HPV, y = Score, fill = HPV)
  ) +
    geom_violin(
      trim = FALSE,
      width = 0.75,
      alpha = 0.55,
      color = "black",
      linewidth = 0.4
    ) +
    geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      alpha = 0.75,
      linewidth = 0.4
    ) +
    geom_jitter(
      width = 0.08,
      size = 2.2,
      shape = 21,
      color = "black",
      alpha = 0.9
    ) +
    annotate(
      "segment",
      x = 1,
      xend = 2,
      y = ymax + 0.08 * yrange,
      yend = ymax + 0.08 * yrange,
      linewidth = 0.5
    ) +
    annotate(
      "text",
      x = 1.5,
      y = ymax + 0.18 * yrange,
      label = "n.s.",
      size = 4.5,
      family = PUBLICATION_FONT
    ) +
    scale_fill_manual(
      values = c(
        "HPV-negative" = "#3B6FB6",
        "HPV-positive" = "#D84B4B"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.22))
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "ssGSEA score"
    ) +
    theme_classic(
      base_size = 14,
      base_family = PUBLICATION_FONT
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 15
      ),
      axis.text = element_text(
        color = "black",
        size = 11
      ),
      axis.title.y = element_text(
        size = 12
      )
    )
}

# -----------------------------------------------------
# Create panels
# -----------------------------------------------------

p_stem <- make_panel(
  filter(plot_df, Program == "Stem-like"),
  "Stem-like"
)

p_bile <- make_panel(
  filter(plot_df, Program == "Bile acid metabolism"),
  "Bile acid metabolism"
)

# -----------------------------------------------------
# Combine as Figure 5A
# -----------------------------------------------------

p5a <- p_stem + p_bile

print(p5a)

# =====================================================
# Combine final Figure 5
# A: Representative PTT distributions
# B: Cliff's delta summary
# =====================================================

library(patchwork)

