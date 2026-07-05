# ccRCC-disulfidptosis-signature

Analysis code for a disulfidptosis-related prognostic signature in clear cell
renal cell carcinoma.

This repository contains the R and Python scripts used for data processing,
model construction, validation, visualization, and reproducibility checks in
the manuscript. Public source datasets are not redistributed here; see
[Data sources](#data-sources) below.

## Repository layout

```text
code/kirc/              Analysis, validation, and figure-generation scripts
FIGURE_CODE_MAP.md      Mapping between manuscript figures and source scripts
```

## Main workflow

Run scripts from the repository root. The expected working structure is:

```text
.
├── code/kirc/
├── data/raw/
├── results/
└── submission_BMC/
```

The consolidated figure builder is:

```bash
Rscript code/kirc/assemble_v3.R
```

This script sources the figure-panel scripts, writes outputs to `results/`, and
copies manuscript-ready Figures 1-9 into `submission_BMC/figures`,
`submission_BMC/figures_TIFF_300dpi`, and `submission_BMC/supplementary`.

## Pipeline order

| Step | Script | Main output |
|---|---|---|
| Prognosis core | `code/kirc/01_disulfidptosis_prognosis.R` | Cox/LASSO-Cox model, risk score, time ROC |
| Subtypes and immune context | `code/kirc/02_subtype_immune.R` | Consensus subtypes and immune contrasts |
| DEG and enrichment | `code/kirc/03_DEG_enrichment.R` | Tumour-vs-normal DEGs and GO enrichment |
| Internal validation | `code/kirc/internal_validation.R` | Split-cohort validation |
| External validation | `code/kirc/extval_cptac.R`, `code/kirc/external_validation_v2.R` | CPTAC and GSE29609 validation |
| GSEA | `code/kirc/supp_gsea.R`, `code/kirc/supp_gsea_v2.R` | Supplementary Fig. S1 |
| Drug sensitivity | `code/kirc/run_gdsc.R`, `code/kirc/fig9_drug_v2.R` | GDSC2 oncoPredict outputs |
| Figure assembly | `code/kirc/assemble_v3.R` | Main Figures 1-9 and supplementary figures |

See `FIGURE_CODE_MAP.md` for the detailed figure-to-script mapping.

## Data sources

The scripts use public datasets that should be downloaded by users from the
original sources:

- TCGA-KIRC expression, survival, and clinical data: UCSC Xena
- ccRCC single-cell RNA-seq: GEO GSE159115
- CPTAC ccRCC: cBioPortal study `rcc_cptac_gdc` and GDC CPTAC-3 follow-up
- GSE29609 external validation cohort: GEO
- GDSC2 drug response data: Genomics of Drug Sensitivity in Cancer

Large input data, intermediate RDS files, and generated figures are intentionally
excluded from this repository.

## Software

The analysis was written primarily in R, with small Python helper scripts for
public data retrieval and reference checks. Package versions should be reported
from the analysis environment used for manuscript submission.

## License

This code is released under the MIT License.
