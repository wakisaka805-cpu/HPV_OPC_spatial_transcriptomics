source(file.path("R", "config.R"))

# =====================================================
# Rebuild final PTT ssGSEA scores from the processed
# GSE228432 expression matrix using the exact archived
# eight-gene-set definitions used in the manuscript.
# =====================================================

library(dplyr)
library(tidyr)
library(GSVA)

expr_file <- file.path(PTT_ROOT, "PTT_expression_gene_symbol.rds")
anno_file <- file.path(PTT_ROOT, "PTT_ssGSEA_scores_HPV.csv")
gs_file <- file.path(REPO_ROOT, "gene_sets", "TCGA_tonsil40_gene_sets_8.rds")
out_file <- file.path(PTT_ROOT, "PTT_ssGSEA_scores_HPV_v2_TexLike.csv")

for (f in c(expr_file, anno_file, gs_file)) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

expr <- readRDS(expr_file)
expr <- as.matrix(expr)
mode(expr) <- "numeric"

archive <- readRDS(gs_file)
if (is.null(archive$gene_sets_full)) {
  stop("Archived gene-set object does not contain gene_sets_full.")
}
gene_sets <- archive$gene_sets_full

# Harmonize the historical Tex-like internal name with the manuscript name.
names(gene_sets)[names(gene_sets) == "TEX_LIKE_CORE18_SHI2023"] <-
  "TEX_LIKE_PROGRAM_SHI2023"

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

missing_sets <- setdiff(target_signatures, names(gene_sets))
if (length(missing_sets) > 0) {
  stop("Missing archived gene sets: ", paste(missing_sets, collapse = ", "))
}

gene_sets <- gene_sets[target_signatures]
gene_sets <- lapply(gene_sets, function(x) intersect(x, rownames(expr)))

ssgsea_param <- GSVA::ssgseaParam(
  exprData = expr,
  geneSets = gene_sets,
  normalize = TRUE
)
ssgsea <- GSVA::gsva(ssgsea_param, verbose = FALSE)

scores <- as.data.frame(t(ssgsea), check.names = FALSE)
scores$GSM <- rownames(scores)

long_df <- scores %>%
  pivot_longer(
    cols = -GSM,
    names_to = "Signature",
    values_to = "Score"
  )

anno_raw <- read.csv(anno_file, stringsAsFactors = FALSE)
anno <- anno_raw %>%
  select(GSM, HPV, Metastasis) %>%
  distinct()

final_df <- long_df %>%
  left_join(anno, by = "GSM") %>%
  select(Signature, GSM, Score, HPV, Metastasis)

if (any(is.na(final_df$HPV))) {
  stop("HPV annotation is missing for one or more PTT samples.")
}

write.csv(final_df, out_file, row.names = FALSE)
cat("Wrote: ", out_file, "\n", sep = "")
print(table(final_df$Signature, final_df$HPV))
