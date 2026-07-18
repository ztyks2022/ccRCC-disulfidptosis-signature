# ccRCC disulfidptosis signature

Analysis code accompanying the manuscript, *An externally evaluated disulfidptosis-related signature stratifies prognosis and reflects immune context in clear cell renal cell carcinoma*.

This repository contains the analysis and figure-generation source code for a public-data study of a 12-gene disulfidptosis-related prognostic signature in clear cell renal cell carcinoma (ccRCC). It includes model development, internal and external evaluation, ClearCode34 benchmarking, immune and molecular-subtype analyses, single-cell analyses, drug-sensitivity prediction, machine-learning robustness checks, and CellChat analysis.

## Repository layout

```text
code/
  kirc/
    README.md               workflow and figure-to-script map
    write_session_info.R    regenerates the software-environment snapshot
    R and Python analysis scripts
sessionInfo.txt             R version and package versions used for this release
```

## Reproducibility

The public datasets are not redistributed here. The data sources, accessions, and script order are described in `code/kirc/README.md`. Some scripts expect the local `data/` and `results/` directories used during the analysis; paths and source-specific access requirements should be configured before rerunning.

`sessionInfo.txt` records the R version and package versions used for the release. To update it after rerunning the pipeline from the repository root:

```bash
Rscript code/kirc/write_session_info.R
```

## Public data sources

- TCGA-KIRC expression, clinical, and survival data: UCSC Xena and the Genomic Data Commons.
- CPTAC ccRCC RNA-seq data: cBioPortal study `rcc_cptac_gdc`, with survival follow-up from GDC project CPTAC-3.
- ccRCC single-cell RNA-seq: GEO accession GSE159115.
- External two-colour array cohort: GEO accession GSE29609.
- Drug-response training data: Genomics of Drug Sensitivity in Cancer (GDSC2).

## License

The repository is distributed under the MIT License.
