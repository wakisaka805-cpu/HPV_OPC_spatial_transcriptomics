# HPV-associated transcriptional programs across OPSCC tissue compartments

Analysis code accompanying the manuscript **“Spatial Organization of HPV-Associated Immune Programs Across Tissue Compartments in Oropharyngeal Cancer.”**

## Purpose

This repository contains the author-generated R code underlying the revised analyses of TCGA-OPSCC, GeoMx DSP, SpatialDecon, surrounding palatine tonsillar tissue (PTT), and exploratory overall survival. Superseded prototypes, local working files, machine-specific paths, and non-public study data are excluded from the release.

## Analysis scripts

1. `R/01_TCGA_OPSCC_primary_marker_survival.R` — strict TCGA-OPSCC cohort, eight-signature ssGSEA, HPV analyses, marker-gene analyses, and exploratory survival analyses.
2. `R/02_TCGA_clinical_characteristics.R` — clinicopathological summary by HPV status.
3. `R/03_Figure2_TCGA_OPSCC.R` — TCGA Figure 2.
4. `R/04_GeoMx_patient_level_analysis.R` — patient-level HPV comparisons, paired LFR–TTR analysis, exact HPV-by-compartment permutation, and exploratory ROI checks.
5. `R/05_GeoMx_SpatialDecon_analysis.R` — SpatialDecon/safeTME analysis and exploratory composition-adjusted models.
6. `R/06_Figure3_paired_LFR_TTR.R` — paired GeoMx Figure 3.
7. `R/07_Figure4_HPV_interaction_SpatialDecon.R` — Figure 4 and the SpatialDecon supplementary figure.
8. `R/08_PTT_ssGSEA_rebuild.R` — rebuilds final PTT ssGSEA scores from the processed GSE228432 expression matrix using the archived eight-gene-set definitions.
9. `R/09_Figure5A_PTT_distributions.R` and `R/10_Figure5B_PTT_effect_sizes.R` — PTT Figure 5 analyses.
10. `R/11_FigureS1_TCGA_survival.R` — exploratory TCGA-OPSCC survival figure.

`R/reference/05B_run_ssGSEA_original.R` is retained as an archival reference for the original ssGSEA implementation.

## Quick start

Run scripts from the repository root. First place the required inputs as described in `DATA_INPUTS.md`, then run:

```r
source("R/check_inputs.R")
```

Paths may be overridden without editing code:

```r
Sys.setenv(TCGA_ROOT = "D:/my_tcga_inputs")
source("R/01_TCGA_OPSCC_primary_marker_survival.R")
```

Recommended analysis order is `01` through `11`; Figure scripts depend on results produced by preceding analysis scripts.

## Gene sets

The repository contains `gene_sets/S1_Table_gene_sets.xlsx` and the exact archived eight-gene-set R object used by the computational workflow, `gene_sets/TCGA_tonsil40_gene_sets_8.rds`. The historical internal Tex-like name is harmonized in the release code to the manuscript name `TEX_LIKE_PROGRAM_SHI2023`.

## Data availability

Public third-party datasets should be obtained from their original repositories. Study-derived GeoMx inputs are not redistributed in this code archive and must be handled according to the manuscript Data Availability Statement and institutional requirements. See `DATA_INPUTS.md`.

## Statistical framework

- Patient is the primary statistical unit for revised GeoMx inference.
- ROI-level analyses are exploratory.
- Benjamini–Hochberg correction is applied within prespecified testing families.
- LFR–TTR comparisons use paired within-patient inference.
- HPV-by-compartment analysis uses exact patient-label permutation.
- TCGA overall-survival analyses are exploratory and include HPV-adjusted and HPV-stratified sensitivity analyses.

## Reproducibility validation

The curated release was numerically re-run on the original Windows analysis environment on 2026-09-13. The TCGA primary/adjusted analyses, GeoMx patient-level analysis, SpatialDecon analysis, Figures 3–5, and Supplementary Figure S1 were checked against the revised manuscript workflow. Minor differences in the final decimal digits of some exploratory Cox P values reflect the current R/package environment and do not alter manuscript conclusions. `sessionInfo.txt` records the validated software environment.

## Outputs

Generated tables are written below `results/`; manuscript figure files are written below `results/PLOS_ONE_Revision_Figures/`. The `results/` directory is ignored by Git except for its placeholder.

## License

Released under the MIT License. See `LICENSE`.

## Citation

Citation metadata are provided in `CITATION.cff`. After Zenodo archives the GitHub release, add the assigned DOI to `CITATION.cff` and the manuscript Data Availability/Code Availability statement as appropriate.
