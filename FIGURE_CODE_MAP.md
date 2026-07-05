# Figure ↔ code map (9 main figures + 3 supplementary)

> **Every figure in this submission has real, runnable source code, listed below. No figure is
> produced by any step outside these scripts.**

**Working directory = the repository/project root, where `code/`, `data/`, `results/`, and
`submission_BMC/` live. Edit scripts under `code/kirc/`.**

**Single command regenerates ALL 9 main + 3 supplementary figures AND refreshes the submission
package (PNG/SVG where available + 300-dpi TIFF; PDF is optional and depends on the local R font device):**

```bash
Rscript code/kirc/assemble_v3.R
```

`assemble_v3.R` sources most panel scripts in an isolated environment (`grab()`, with
`ASSEMBLE_MODE=TRUE` so the panel scripts build their panels but skip their own save),
recomposes the panels with continuous A.. lettering, writes `results/Fig*_v3.{png,svg,pdf}`,
then copies them to `submission_BMC/figures/Figure{N}.*` and writes the TIFFs. So **editing a
panel script and re-running the one command updates both `results/` and the package** — one-to-one.

| Main figure | Panels | Controlled by |
|---|---|---|
| **Figure 1** landscape | A score · B DEG volcano · C key genes · D gene-correlation heatmap | `fig1_landscape_v2.R` (pA/pB/pC/pD) |
| **Figure 2** signature | A univariate-Cox forest · B risk KM · C time-ROC | `fig3_prognosis_v2.R` (pA/pB/pC) |
| | D–F gene-class subset (cyto vs metab) | `subset_cyto_vs_metab.R` (pA/pB/pC) |
| **Figure 3** prognosis + subtype | A nomogram | `fig4_nomogram_v2.R` (object `nom`) |
| | B calibration · C multivariate-Cox forest · D decision-curve analysis | **coded inside `assemble_v3.R`** (`calP`, `forP`) + `fig_dca.R` (`pDCA`) |
| | E subtype heatmap · F subtype KM · G subtype immune · H subtype activity | `fig2_subtype_v2.R` (pA/pB/pC/pD) |
| **Figure 4** validation | A internal KM · B CPTAC KM · C GSE29609 KM | `fig11_validation.R` (pA/pB/pC) |
| | D ClearCode34 C-index | `supp_cindex.R` (p) |
| **Figure 5** immune | A infiltration · B checkpoints · C ESTIMATE | `fig5_immune_v2.R` (pA/pB/pC) |
| **Figure 6** context | A scUMAP · B activity · C violin · D dotplot | `fig6_singlecell_v2.R` (pA/pB/pC/pD) |
| | E WGCNA module–trait | `fig7_wgcna_v2.R` (pHt) |
| | F candidate Venn · G SYCE1L KM | `fig8_candidate_v2.R` (pV/pK) |
| **Figure 7** therapy | A drug lollipop · B IC50 | `fig9_drug_v2.R` (pA/pB) |
| | C ACTB-KO perturbed genes · D KO GO | `fig10_virtualKO_v2.R` (pA/pB) |
| **Figure 8** machine learning | A cross-validated AUC · B GBM importance · C SHAP importance · D SHAP summary | `fig_ml_shap.R` |
| **Figure 9** cell-cell communication | A interaction count · B interaction strength · C sender/receiver roles · D MHC-II heatmap · E tumour-to-immune signalling | `fig_cellchat.R` |
| **Suppl. Fig. S1** GSEA (high vs low risk) | — | `supp_gsea_v2.R` (whole figure) |
| **Suppl. Fig. S2** co-expression network | — | `fig12_network.R` (whole figure) |
| **Suppl. Fig. S3** somatic mutation landscape | — | `fig_tmb.R` (`results/Fig_oncoplot.png`, copied by `assemble_v3.R`) |

Supplementary Figures S1 and S2 are produced by standalone scripts; S3 is the oncoplot rendered by
`fig_tmb.R` and copied by `assemble_v3.R`. They are published to `submission_BMC/supplementary/` +
TIFF in the same single command.

> **Two panels live in `assemble_v3.R` itself, not in a `fig*_v2.R` script:** Figure 4B (calibration)
> and 4C (forest) were rebuilt as clean ggplot inside the assembler. To change them, edit
> `assemble_v3.R` (the `calP` / `forP` blocks), not `fig4_nomogram_v2.R`.
>
> The old standalone `fig*_v2.R` scripts still also write their own `results/Fig*_v2.*` files when run
> directly (without `ASSEMBLE_MODE`); those v2 files are the pre-consolidation versions and are **not**
> used in the submission package.
