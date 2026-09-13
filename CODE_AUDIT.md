# Code curation and validation audit

## Curated release

The original working archive contained 123 R files, including superseded prototypes, duplicate versions, empty files, and scripts from earlier analytical frameworks. The public release retains only the scripts supporting the revised manuscript's final analytical framework.

## Safety and portability review

- Machine-specific absolute Windows paths were removed from retained executable code.
- No passwords, GitHub tokens, or API keys were detected during curation.
- Private/local inputs are excluded by `.gitignore` (`data/*`).
- Generated results are excluded by `.gitignore` (`results/*`).
- Repository-wide relative paths and environment-variable overrides are defined in `R/config.R`.
- The exact archived eight-gene-set definition object is included under `gene_sets/`.
- The uncurated 123-script source archive must not be uploaded to GitHub.

## Numerical validation

Local validation was performed on the original Windows analysis environment on 2026-09-13. The following components were run successfully for the revised workflow:

- TCGA-OPSCC primary, adjusted, marker, and exploratory survival analyses.
- GeoMx patient-level HPV and paired LFR–TTR analyses.
- Exact HPV-by-compartment permutation analysis.
- SpatialDecon/safeTME analyses and composition-adjusted exploratory models.
- Figure 3.
- Figure 4 and SpatialDecon supplementary visualization.
- PTT Figure 5A and Figure 5B analyses.
- Supplementary TCGA survival figure.
- `sessionInfo.txt` generation.

The SpatialDecon working script originally contained an obsolete post-save tail; the retained release copy ends with the validated formal save/check block. Figure 4 no longer depends on a manually created root-level CSV; it reads the formal SpatialDecon output from `results/GeoMx/`.

Minor last-decimal differences in some exploratory Cox P values were observed under the current package environment; they do not change inference or manuscript conclusions. The manuscript/Supplementary tables remain the reporting reference.

## Release readiness

This code package is suitable for the versioned GitHub release after the repository metadata are reviewed. Keep the GitHub repository private until the final repository diff confirms that no restricted `data/` files are staged. After that, create release `v1.0.0`, archive it with Zenodo, and add the Zenodo DOI to `CITATION.cff`.
