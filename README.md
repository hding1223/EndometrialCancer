# TCGA-UCEC Immune-Related Prognostic Analysis Scripts

This directory contains the R scripts used for the computational analysis of immune-related genes, microsatellite status, prognosis, tumor microenvironment, single-gene associations, drug sensitivity, and model visualization in uterine corpus endometrial carcinoma (TCGA-UCEC).

The scripts are numbered according to the current analysis workflow. The numbering is not completely continuous because scripts `19`, `26`, and `27` are not present in this directory.

## Analysis Workflow

The recommended high-level workflow is:

```text
Input data
  ├─ WGCNA and immune-gene preprocessing
  ├─ Differential-expression analysis
  ├─ Univariate Cox screening
  ├─ Stability selection
  ├─ Machine-learning prognostic modeling
  ├─ Model interpretation and nomogram construction
  ├─ Functional enrichment and immune-microenvironment analysis
  ├─ TMB, immune-checkpoint, IPS, and drug-sensitivity analysis
  └─ Independent-cohort assessment and visualization
```

The main analysis branch is:

```text
01.WGCNA.R
  → 02.IRG-Differ-Analysis-wilcoxon.R
  → 03.IRG-Cox.R
  → 04.stability_pipeline.R
  → 05.IRG-MSS-Machine-Learning.R / 07.IRG-MSS-Machine-Learning.R
  → 16–18 nomogram and model-evaluation scripts
  → 20–44 downstream analyses
```

Some scripts are alternative implementations or visualization-only scripts. They should not all be executed blindly as one linear pipeline.

## Script Inventory

### 1. Gene Preprocessing, WGCNA, and Cox Screening

| File | Main purpose | Typical outputs or downstream inputs |
|---|---|---|
| [`01.WGCNA.R`](01.WGCNA.R) | Filters TCGA-UCEC expression data, extracts immune-related genes, performs a 70:30 OS-stratified split, constructs a WGCNA network, identifies modules, and exports module genes. | `trainSamples.txt`, WGCNA figures, module gene lists, `module_all.txt`, `hubGenes*.txt` |
| [`02.IRG-Differ-Analysis-wilcoxon.R`](02.IRG-Differ-Analysis-wilcoxon.R) | Performs tumor-versus-normal differential-expression analysis for immune-related genes using the Wilcoxon rank-sum test. | Differential-expression table, heatmap, volcano plot, DEG lists |
| [`03.IRG-Cox.R`](03.IRG-Cox.R) | Intersects WGCNA genes and DEGs, restricts the analysis to the WGCNA training samples, and performs univariate Cox regression. | Univariate Cox results and forest plot |
| [`04.stability_pipeline.R`](04.stability_pipeline.R) | Repeats the patient-level 7:3 split and the WGCNA → DEG → Cox workflow across multiple seeds; estimates gene-selection frequency, Jaccard stability, proportional-hazards diagnostics, and optional bootstrap intervals. | `summary_long.tsv`, `run_summary.tsv`, `stage_sets.rds`, `gene_frequency.tsv`, `robust_genes.tsv`, stability figures |
| [`09.MMS-WGCNA.R`](09.MMS-WGCNA.R) | Performs a WGCNA workflow for the MSS/MSI-related analysis branch, including sample clustering, soft-threshold selection, module detection, and module-trait relationships. | WGCNA module files and figures |
| [`10.MMS-WGCNA.R`](10.MMS-WGCNA.R) | Alternative or extended WGCNA implementation for the MSS/MSI-related branch. | WGCNA module files and figures |
| [`11.DEGs-MMS-WGCNA.R`](11.DEGs-MMS-WGCNA.R) | Performs differential-expression analysis associated with the MSS/MSI-WGCNA branch. | DEG tables and visualization files |
| [`12.MMS-WGCNA.R`](12.MMS-WGCNA.R) | Additional MSS/MSI WGCNA analysis and module visualization. | Module-trait plots, module gene lists, and network files |
| [`13.Genes-Relation.R`](13.Genes-Relation.R) | Examines the relationship or overlap between selected gene sets. | Gene-overlap tables or plots |
| [`14.Genes-Prognosis.R`](14.Genes-Prognosis.R) | Performs prognosis-related analysis for selected genes or gene sets. | Prognostic association results and plots |

### 2. Machine-Learning Modeling and Model Visualization

| File | Main purpose | Typical outputs or downstream inputs |
|---|---|---|
| [`05.IRG-MSS-Machine-Learning.R`](05.IRG-MSS-Machine-Learning.R) | Builds and evaluates prognostic models using the MIME/ML.Dev.Prog.Sig framework and multiple machine-learning combinations, including StepCox and GBM. Includes validation-set C-index bootstrap analysis and risk-score reproducibility code. | Model comparison tables, risk scores, model objects, C-index results, ROC/KM-related outputs |
| [`06.MMS-Machine-Learning-drawmap.R`](06.MMS-Machine-Learning-drawmap.R) | Generates figures from the MSS/MSI machine-learning results. | Model-performance plots |
| [`07.IRG-MSS-Machine-Learning.R`](07.IRG-MSS-Machine-Learning.R) | Fits or inspects the final immune-related MSS prognostic model, including StepCox[forward] + GBM gene parameters, relative influence, partial dependence, and risk-score direction. | Model parameter tables, PDP plots, gene-effect summaries |
| [`08.IRG-MSS-Machine-Learning-drawmap.R`](08.IRG-MSS-Machine-Learning-drawmap.R) | Generates visualization figures for the immune-related MSS machine-learning results. | Model visualization figures |
| [`15.MSS-Machine-Learning-drawmap.R`](15.MSS-Machine-Learning-drawmap.R) | Produces additional plots for MSS machine-learning model results. | Model-performance and risk-score figures |
| [`16.MSS-nomogram-1.R`](16.MSS-nomogram-1.R) | Prepares clinical variables and risk scores for nomogram construction. | Nomogram input data |
| [`17.MSS-nomogram-2.R`](17.MSS-nomogram-2.R) | Constructs a Cox-based nomogram using clinical variables and the prognostic risk score. | Nomogram figure and model object |
| [`18.MSS-nomogram-3.R`](18.MSS-nomogram-3.R) | Alternative or extended nomogram preparation and analysis. | Nomogram-related data or figures |
| [`32.MSS-calibration.R`](32.MSS-calibration.R) | Evaluates calibration of predicted 1-, 3-, and 5-year survival probabilities. | Calibration plots |
| [`33.MSS-DCA.R`](33.MSS-DCA.R) | Performs decision-curve analysis for the prognostic model. | DCA curves and net-benefit summaries |

### 3. Functional Enrichment and Tumor-Microenvironment Analysis

| File | Main purpose | Typical outputs or downstream inputs |
|---|---|---|
| [`20.MSS-Gene-GO.R`](20.MSS-Gene-GO.R) | Performs Gene Ontology enrichment analysis for the selected gene set. | `IRG_GO.txt`, GO bar plot, bubble plot, and circular plot |
| [`21.MSS-Gene21.KEGG.R`](21.MSS-Gene21.KEGG.R) | Performs KEGG pathway enrichment analysis. | `IRG_KEGG.txt`, KEGG bar plot, and bubble plot |
| [`22.MSS-Gene22.GSEA.R`](22.MSS-Gene22.GSEA.R) | Performs gene-set enrichment analysis using ranked expression or gene statistics. | GSEA result table and enrichment plot |
| [`23.MSS-Gene23.estimate.R`](23.MSS-Gene23.estimate.R) | Calculates ESTIMATE-derived stromal, immune, and composite scores. | ESTIMATE score files |
| [`24.MSS-Gene24.TMEvioplot.R`](24.MSS-Gene24.TMEvioplot.R) | Visualizes tumor-microenvironment scores between risk or molecular subgroups. | TME violin plots |
| [`25.MSS-CiberSort.R`](25.MSS-CiberSort.R) | Estimates immune-cell proportions using the CIBERSORT LM22 signature matrix. | Immune-cell proportion tables and plots |
| [`28.MSS-Gene28.checkpoint.R`](28.MSS-Gene28.checkpoint.R) | Correlates the prognostic risk score with immune-checkpoint gene expression. | Checkpoint correlation tables and heatmaps |
| [`29.MSS-Gene29.TMBcor.R`](29.MSS-Gene29.TMBcor.R) | Calculates tumor mutational burden and evaluates its association with risk scores. | TMB table, correlation statistics, and plots |
| [`31.MSS-Gene31.IPS.R`](31.MSS-Gene31.IPS.R) | Compares TCIA immune phenotype scores between prognostic groups. | IPS comparison plots and statistics |
| [`39.MSS_SubGroup.IPS.R`](39.MSS_SubGroup.IPS.R) | Performs IPS analysis within MSS-related subgroups. | Subgroup IPS plots |
| [`40.MSI-H_SubGroup.IPS.R`](40.MSI-H_SubGroup.IPS.R) | Performs IPS analysis within MSI-H subgroups. | MSI-H subgroup IPS plots |

### 4. Drug Sensitivity and Drug-Relationship Analysis

| File | Main purpose | Typical outputs or downstream inputs |
|---|---|---|
| [`30.MSS-Gene30.pRRophetic.R`](30.MSS-Gene30.pRRophetic.R) | Predicts drug sensitivity using pRRophetic and GDSC/cgp2016 expression data. | Predicted IC50 values and drug-sensitivity plots |
| [`41.MSS-MSI-H.pRRophetic.R`](41.MSS-MSI-H.pRRophetic.R) | Compares predicted drug sensitivity between MSS and MSI-H or related risk groups. | Differential drug-sensitivity results |
| [`42.Drugs-Relation.R`](42.Drugs-Relation.R) | Examines overlaps among drug candidate sets. | Drug-overlap tables or Venn diagrams |

### 5. External Assessment and Single-Gene Visualization

| File | Main purpose | Typical outputs or downstream inputs |
|---|---|---|
| [`34.MSS-Gene16.indep.R`](34.MSS-Gene16.indep.R) | Assesses selected genes in independent GEO datasets or external expression cohorts. | External differential-expression results and plots |
| [`35.MSS-KS-Analysis.R`](35.MSS-KS-Analysis.R) | Performs Kolmogorov-Smirnov or related distributional comparisons for selected genes. | KS-test statistics and plots |
| [`36.MSS_TYPE_.ScoreVioplot.R`](36.MSS_TYPE_.ScoreVioplot.R) | Visualizes risk-score or gene-expression differences across molecular types. | Molecular-type violin plots |
| [`37.MSS_Group_.ScoreVioplot.R`](37.MSS_Group_.ScoreVioplot.R) | Visualizes score differences between MSS-related groups. | Group comparison violin plots |
| [`38.MSI-H_SubGroup_.ScoreVioplot.R`](38.MSI-H_SubGroup_.ScoreVioplot.R) | Visualizes score differences among MSI-H subgroups. | MSI-H subgroup violin plots |
| [`43.MSS_Group_TPMVioplot.R`](43.MSS_Group_TPMVioplot.R) | Visualizes TPM expression differences among MSS-related groups. | TPM violin plots |
| [`44.Normal_MSS_Group_TPMVioplot.R`](44.Normal_MSS_Group_TPMVioplot.R) | Compares gene TPM expression between normal tissue and MSS-related groups. | Normal/MSS expression violin plots |

## Main Input Data

The scripts refer to project-specific files under a hard-coded workspace structure. Common inputs include:

| Input | Role |
|---|---|
| `TCGA_UCEC_TPM.txt` | TCGA-UCEC TPM expression matrix; genes are rows and samples are columns. |
| `TCGA_UCEC_clinical_selected_above_32_MMS.txt` | Clinical and survival data used in the WGCNA/Cox branch. |
| `clinical_above30_ml.txt` | Clinical data used for machine-learning and nomogram analyses. |
| `clinical_above30_ml_MMS.txt` | MSS/MSI-related clinical data for model analyses. |
| `IRG_2499.txt` | Immune-related gene list used as the WGCNA input background. |
| `IRG_UCEC_351.txt` | Alternative immune-related gene list used by some older or alternative scripts. |
| `trainSamples.txt` | Patient/sample IDs generated by the WGCNA split and used by downstream screening scripts. |
| `module_grey.txt` | Gene list exported from WGCNA for the selected module. |
| `DEGs_362.txt` or related DEG lists | Differentially expressed gene sets used for intersection and Cox screening. |
| `final_mode_riskscore.txt` | Model-derived risk scores used by nomogram, survival, subgroup, and downstream analyses. |
| `GDSC/cgp2016` data | Reference expression and drug-response data for pRRophetic analyses. |
| `CIBERSORT LM22` | Immune-cell signature matrix for immune-infiltration estimation. |
| `TCIA IPS` data | Immune phenotype scores used for immunotherapy-related exploratory analysis. |

The exact input filename and directory should be checked in each script before execution.

## Recommended Execution Order

### Core prognostic-model branch

1. Run [`01.WGCNA.R`](01.WGCNA.R) to generate training samples and WGCNA module files.
2. Run [`02.IRG-Differ-Analysis-wilcoxon.R`](02.IRG-Differ-Analysis-wilcoxon.R) for immune-related differential-expression analysis.
3. Run [`03.IRG-Cox.R`](03.IRG-Cox.R) for univariate Cox screening.
4. Run [`04.stability_pipeline.R`](04.stability_pipeline.R) to evaluate selection stability and proportional-hazards diagnostics.
5. Run [`05.IRG-MSS-Machine-Learning.R`](05.IRG-MSS-Machine-Learning.R) or [`07.IRG-MSS-Machine-Learning.R`](07.IRG-MSS-Machine-Learning.R) for model construction and interpretation.
6. Run [`06.MMS-Machine-Learning-drawmap.R`](06.MMS-Machine-Learning-drawmap.R), [`08.IRG-MSS-Machine-Learning-drawmap.R`](08.IRG-MSS-Machine-Learning-drawmap.R), or [`15.MSS-Machine-Learning-drawmap.R`](15.MSS-Machine-Learning-drawmap.R) for model plots.
7. Run [`16.MSS-nomogram-1.R`](16.MSS-nomogram-1.R), [`17.MSS-nomogram-2.R`](17.MSS-nomogram-2.R), and [`18.MSS-nomogram-3.R`](18.MSS-nomogram-3.R) for nomogram analyses.
8. Run [`32.MSS-calibration.R`](32.MSS-calibration.R) and [`33.MSS-DCA.R`](33.MSS-DCA.R) for calibration and decision-curve analysis.

### Downstream biological-analysis branch

After the relevant gene lists and risk scores have been generated, the following scripts can be run as independent downstream analyses:

1. Functional enrichment: `20`–`22`.
2. Tumor microenvironment and immune infiltration: `23`–`25` and `28`–`29`.
3. Drug sensitivity: `30`, `41`, and `42`.
4. Immune phenotype score: `31`, `39`, and `40`.
5. Independent-cohort and subgroup analyses: `34`–`38`, `43`, and `44`.

## R Packages

The exact package requirements vary by script. Frequently used packages include:

- `WGCNA`, `survival`, `survminer`, `caret`, `tidyverse`, and `data.table`;
- `limma`, `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, and `enrichplot`;
- `ggplot2`, `ggpubr`, `pheatmap`, `VennDiagram`, `ggrepel`, `cowplot`, and `ggthemes`;
- `glmnet`, `gbm`, `CoxBoost`, `randomForestSRC`, `survivalSVM`, and `Mime`;
- `regplot`, `rms`, `ggDCA`, and `timeROC`;
- `Seurat` and related single-cell packages where applicable;
- `estimate`, `e1071`, `pRRophetic`, `parallel`, and `RColorBrewer`.

Package versions should be recorded with `sessionInfo()` before final analysis. The scripts were developed around R 4.3.x, but package-version differences may affect model fitting, random-number generation, plots, or output formatting.

## Configuration Before Running

1. Replace hard-coded paths such as `D:/workspace/TCGA/UCEC/...` with paths valid on the local machine.
2. Confirm that every input file exists and that row names, sample IDs, and gene symbols use consistent formats.
3. Confirm whether sample IDs are full TCGA barcodes or truncated 12-character patient IDs.
4. Confirm whether expression values are TPM, log2(TPM + 1), counts, or another normalized scale.
5. Check the training/validation split and random seed before rerunning model development.
6. Confirm that the required supplementary data files, GDSC reference data, CIBERSORT LM22 matrix, and TCIA IPS data are available.
7. Run scripts in a clean R session when possible, because several scripts reuse common object names such as `data`, `clinical`, `train`, `test`, `MEs`, and `risk score` objects.

## Reproducibility Notes

- [`04.stability_pipeline.R`](04.stability_pipeline.R) uses repeated patient-level splits and fixed seeds to assess feature-selection stability.
- The final machine-learning model should be saved as an R object together with the optimal tree number, model parameters, input-gene list, preprocessing rules, and risk-group cutoff.
- When applying a frozen model to an external cohort, use the same gene identifiers, expression transformation, probe-to-gene mapping, and risk-score direction as in the training cohort.
- Keep a copy of the generated `trainSamples.txt`, model-comparison table, final risk-score table, and session information.
- Do not interpret pRRophetic drug sensitivity as therapeutic efficacy without experimental or independent pharmacogenomic validation.
- Do not interpret single-cell cell-level comparisons as patient-level validation when only a small number of patients are available.

## Known Caveats

- Several scripts contain hard-coded input and output directories and may require path editing before use.
- Some scripts are alternative versions of similar analyses, particularly the WGCNA and machine-learning branches. Compare their parameters before deciding which version is the primary analysis.
- The numbering contains gaps: scripts `19`, `26`, and `27` are not present in the current directory.
- Running a script may overwrite files with the same name in the configured output directory.
- The scripts are analysis code rather than an installable R package; functions and objects may depend on files generated by earlier scripts.
- A successful syntax check does not guarantee that all external data files, package versions, or working directories are correctly configured.

## Validation

All R scripts in this directory should be syntax-checked before execution:

```r
parse(file = "src/01.WGCNA.R")
```

For a batch syntax check from PowerShell:

```powershell
Get-ChildItem .\src\*.R | ForEach-Object {
    Rscript -e "parse(file='$($_.FullName)')"
}
```

