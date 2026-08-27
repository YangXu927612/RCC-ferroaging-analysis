# RCC-ferroaging-analysis
R code for bioinformatics and statistical analyses of ferro-aging-related molecular features in renal cell carcinoma.
# Ferro-aging-related molecular characteristics and immunological heterogeneity in renal cell carcinoma — Code Repository

This repository contains the **R scripts** (in the `R/` sub-folder) that support the analyses reported in the manuscript *"Integrated single-cell and bulk transcriptomic analyses reveal ferro-aging-related molecular characteristics and immunological heterogeneity in renal cell carcinoma."*

The code covers data integration, batch-effect correction, differential expression, machine-learning-based biomarker discovery, immune infiltration, consensus clustering, WGCNA, GO/KEGG enrichment, single-cell RNA-seq analysis, and in-silico gene knockout (scTenifoldKnk).

> **No absolute paths** are referenced in any script. All inputs/outputs are read from / written to the **current working directory** of the R session.

---

## Repository contents

| #    | File                                    | Purpose                                                      | Manuscript section                   |
| ---- | --------------------------------------- | ------------------------------------------------------------ | ------------------------------------ |
| 01   | `01_data_merge_batch_correction.R`      | Merge the **training** GEO expression matrices (`GSE40435.csv`, `GSE66272.csv`); remove batch effects via `limma::removeBatchEffect`; PCA + boxplot diagnostics. **Explicitly excludes `GSE53757.csv`** (external validation) even if it sits in the same working directory. | Methods §1 (Data collection)         |
| 02   | `02_gene_set_scoring.R`                 | **Optional.** Scores pathway-level mean expression between Control / Disease for user-supplied gene sets. Not required to reproduce any figure / table in the manuscript. | (optional)                           |
| 03   | `03_differential_expression_analysis.R` | LIMMA differential expression; volcano, heatmap, PCA, Venn diagrams with ferro-aging-related gene set. | Methods §1                           |
| 04   | `04_svm_rfe_feature_selection.R`        | SVM-RFE with 10-fold cross-validation; ranks candidate genes. | Methods §2 (Diagnostic model)        |
| 05   | `05_random_forest_gene_selection.R`     | Random-forest importance ranking (MeanDecreaseGini).         | Methods §2                           |
| 06   | `06_candidate_gene_intersection.R`      | Identifies the overlapping candidate genes selected by **SVM-RFE and Random Forest** (the "candidate core genes" reported in the manuscript) and generates the corresponding Venn diagram. | Methods §2                           |
| 07   | `07_diagnostic_model_construction.R`    | Builds a multi-gene logistic-regression diagnostic model; produces nomogram, calibration curve, decision curve, ROC/AUC. | Methods §2                           |
| 08   | `08_external_validation.R`              | Rebuilds the diagnostic model on the training cohort and validates it on the external GEO cohort **GSE53757** (ROC + 95 % CI). | Methods §2                           |
| 09   | `09_cibersort_immune_infiltration.R`    | CIBERSORT-based estimation of 22 immune-cell fractions; gene–immune Spearman correlations. | Methods §3 (Immune microenvironment) |
| 10   | `10_mcp_counter_immune_abundance.R`     | MCP-counter estimation of immune/stromal cell abundance and gene–immune correlations. | Methods §3                           |
| 11   | `11_hub_genes_inflammatory_cytokines.R` | Spearman correlation of candidate core genes with inflammatory cytokines. | Methods §3                           |
| 12   | `12_consensus_clustering.R`             | Consensus clustering on hub candidate gene expression; outputs sample-to-cluster labels (`geneCluster.txt`). | Methods §4 (Consensus clustering)    |
| 13   | `13_clinical_distribution_heatmap.R`    | Heatmap of clinical features across consensus clusters.      | Methods §4                           |
| 14   | `14_clinical_distribution_boxplot.R`    | Box-plot / stacked-bar comparison of clinical variables across clusters. | Methods §4                           |
| 15   | `15_mcp_counter_cluster_comparison.R`   | MCP-counter immune-cell abundance differences among clusters. | Methods §4                           |
| 16   | `16_cytokine_cluster_comparison.R`      | Inflammatory-cytokine differences among clusters.            | Methods §4                           |
| 17   | `17_wgcna_analysis.R`                   | WGCNA on the integrated cohort; module–trait relationships vs cluster subtype. | Methods §5 (WGCNA)                   |
| 18   | `18_go_kegg_enrichment.R`               | GO and KEGG enrichment analyses for hub candidate genes / module hub genes. | Methods §5                           |
| 19   | `19_scrnaseq_comprehensive_analysis.R`  | Single-cell QC, normalization, PCA, UMAP, cell-type annotation (manual), marker genes, Ro/e, module-score, scRNAtoolVis panels, in-silico knockout via **scTenifoldKnk**. | Methods §6 (scRNA-seq & virtual KO)  |

---

## Recommended analysis order

```text
01 → 03 → 04 → 05 → 06 → 07 → 08
                       ↓
                       09 / 10 / 11 (immune-correlation)
                       ↓
                       12 → 13 / 14 (clinical)
                              ↓
                              15 / 16 (immune in clusters)
                              ↓
                              17 → 18 (WGCNA + enrichment)
                              ↓
                              19 (scRNA-seq + virtual KO)
```

`R/02_gene_set_scoring.R` is **not** part of the main pipeline (it is an optional utility for exploratory pathway scoring) and can be run independently from the same working directory if you want to inspect additional pathway-level comparisons.

---

## Working-directory convention (important!)

Every script sets `workDir <- getwd()` (or uses the current working directory). All input files are read **relative to that directory**, so:

1. **Download / clone** this repository to a local folder.

2. Place **all input files** (the GEO expression CSVs, `gene set/` folder, `clinical.xlsx`, `refer.txt`, `scRNA_matrix.rds`, …) in that same folder.

3. In R / RStudio, set the working directory to that folder:

   ```r
   setwd("/path/to/this/repository")
   ```

   or, in RStudio: `Session ▸ Set Working Directory ▸ To Source File Location` after opening any script.

4. Run the scripts in the order indicated above.

All outputs (CSVs, PDFs, RDS files, …) are written into the working directory or a timestamped sub-folder inside it.

---

## Required input files (place in the working directory)

| Input file                                                   | Used by              | Source                                                       |
| ------------------------------------------------------------ | -------------------- | ------------------------------------------------------------ |
| `GSE40435.csv` (geneSymbol × sample)                         | `01`                 | GEO series matrix + probe annotation. Sample names should keep `_con` / `_tre` suffix. |
| `GSE66272.csv` (geneSymbol × sample)                         | `01`                 | GEO series matrix + probe annotation. Sample names should keep `_con` / `_tre` suffix. |
| `merged_after_batch_removal.csv`                             | 03 – 17              | Output of `01_data_merge_batch_correction.R`. Sample names must end with `_con` (Control) or `_tre` (Disease). |
| `gene set/` (folder of `.txt` files, one gene per line, header optional) | 02, 03               | Ferro-aging-related gene lists; e.g. the 95-gene ferro-aging set reported by Liu *et al.* (Cell Metab, 2026). |
| `Intersect_Genes.csv`                                        | 06 – 11, 15 – 16, 19 | Output of `06_candidate_gene_intersection.R` (intersection of SVM-RFE & random-forest gene sets). |
| `geneCluster.txt` (sample × cluster)                         | 12 – 16              | Output of `12_consensus_clustering.R`.                       |
| `clinical.xlsx`                                              | 13, 14               | User-provided clinical annotation.                           |
| `refer.txt` (LM22 reference)                                 | 09                   | CIBERSORT LM22 signature matrix. Download from the CIBERSORT portal. |
| `Inflammatory Cytokines gene.txt`                            | 11, 16               | User-provided list of inflammatory cytokines.                |
| `GSE53757.csv` (geneSymbol × sample)                         | 08                   | Independent external validation cohort.                      |
| `scRNA_matrix.rds`                                           | 19                   | Single-cell expression matrix (`dgCMatrix` / `matrix`) with cell-barcode names prefixed `control.` or `Disease.`. |

---

## Public data sources

| Accession     | Type                                  | Used in     |
| ------------- | ------------------------------------- | ----------- |
| **GSE40435**  | Bulk microarray (training)            | 01, 03 – 17 |
| **GSE66272**  | Bulk microarray (training)            | 01, 03 – 17 |
| **GSE53757**  | Bulk microarray (external validation) | 08          |
| **GSE159115** | Single-cell RNA-seq                   | 19          |

All datasets can be retrieved from the NCBI **Gene Expression Omnibus (GEO)**.

---

## Data preparation chain (GEO → per-cohort CSV → merged CSV)

The analysis scripts start from already-formatted matrices. The full chain expected by the scripts is:

```
GEO accession       (NCBI)
   ↓ download series matrix + annotate probes to gene symbols
   ├── GSE40435.csv ─┐
   │                 ├──→ R/01_data_merge_batch_correction.R
   ├── GSE66272.csv ─┘        (only these two are read;
   │                          GSE53757.csv is ignored here)
   │
   └── GSE53757.csv ──→ R/08_external_validation.R
                          (used independently, never merged with training)
                              │
   ↓                              ↓
merged_after_batch_removal.csv  external-validation ROC
   ↓                              ↓
scripts in R/02 ... R/17         Figures (Fig S1)
   ↓
figures / tables (Fig 1 ... Fig S5)
```

### Inputs expected by R/01_data_merge_batch_correction.R

`01_data_merge_batch_correction.R` reads **only the two training cohorts listed in its `input_files` vector** (`GSE40435.csv`, `GSE66272.csv`). It does NOT use `list.files()`, so any other CSVs in the working directory (e.g. the external validation `GSE53757.csv`, or output CSVs produced by earlier analysis steps) are explicitly ignored. If `GSE53757.csv` is present in the working directory the script emits a warning, but the file is still never merged into the training matrix.

| Place in working directory           | Produced by                                  | Used by                              |
| ------------------------------------ | -------------------------------------------- | ------------------------------------ |
| `GSE40435.csv` (geneSymbol × sample) | GEO series matrix + probe → gene annotation  | `R/01_data_merge_batch_correction.R` |
| `GSE66272.csv` (geneSymbol × sample) | GEO series matrix + probe → gene annotation  | `R/01_data_merge_batch_correction.R` |
| → `merged_after_batch_removal.csv`   | Output of `01_data_merge_batch_correction.R` | Scripts `R/02 – R/17`                |

Other inputs listed below in [Required input files](#required-input-files-place-in-the-working-directory) are read by the rest of the pipeline.

### Building the per-cohort CSVs from GEO

The scripts do not download from GEO or annotate probes — that preprocessing is expected to be done by the user (or by separate preprocessing scripts) before reaching this repository. A reasonable workflow:

1. Download the series matrices for **GSE40435**, **GSE66272**, and **GSE53757** from NCBI GEO. (`GSE159115`, used by `R/19`, is single-cell RNA-seq and arrives in a different format — see the `scRNA_matrix.rds` step below.)
2. For each microarray cohort, map probes to gene symbols using its platform annotation package (e.g. `hgu133plus2.db`, `hgu133a.db`).
3. Save each cohort as a `geneSymbol × sample` CSV. Sample names should keep a `_con` / `_tre` suffix when applicable (this is how the downstream scripts infer group labels).
4. **Place `GSE40435.csv` and `GSE66272.csv` in the working directory and run `R/01_data_merge_batch_correction.R`.** `GSE53757.csv` is consumed separately by `R/08_external_validation.R` and must NEVER be merged with the training data.

---

## Required R packages

The pipeline depends on packages from three sources. The lists below are installation hints only — install once, before running the pipeline. Note: this repository does not ship a captured `sessionInfo()`; the package versions you record should match the R installation that produced the manuscript figures.

R version used for the manuscript: **4.2.2** (per Methods §1).

### CRAN packages

```r
install.packages(c(
  # Core data manipulation and plotting
  "data.table", "dplyr", "tidyr", "magrittr", "stringr", "broom",
  "ggplot2", "ggpubr", "ggrepel", "ggsignif", "ggsci",
  "RColorBrewer", "viridis", "scales", "patchwork",
  "pheatmap", "ggvenn", "circlize",
  "reshape2", "readxl", "openxlsx", "DT", "progress",
  "svglite", "corrplot", "cowplot",
  "tools", "parallel",

  # Statistical / ML
  "e1071", "randomForest", "glmnet",
  "rms", "regplot", "rmda", "pROC",

  # Single-cell — CRAN-only parts
  "Seurat", "DoubletFinder",
  "Matrix", "clustree", "cluster", "ggunchull",
  "scRNAtoolVis", "plot1cell", "rstatix",

  # QC normalisation
  "preprocessCore"
))
```

### Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "limma",                  # differential expression (script 03)
  "ComplexHeatmap",         # complex heatmaps         (script 03, 19)
  "clusterProfiler",        # GO/KEGG enrichment       (script 18)
  "org.Hs.eg.db",           # human gene IDs           (script 18)
  "MCPcounter",             # immune cell abundance    (scripts 10, 15)
  "WGCNA",                  # co-expression network    (script 17)
  "SingleR",                # cell-type annotation     (script 19)
  "celldex",                # reference datasets       (script 19)
  "celda"                   # decontX ambient-RNA removal (script 19)
))
```

### GitHub-only packages (not on CRAN / Bioconductor)

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# In-silico gene knockout (script 19 — scTenifoldKnk virtual KO pipeline)
remotes::install_github("cailab-tamu/scTenifoldKnk")

# GOplot: chord / circular plots (script 18).
# GOplot is on CRAN at older releases; if it has been archived, install from GitHub:
#   remotes::install_github("maxplanck-ie/Goplot")
```

> **Why no version pins here?** Manually pinning dozens of package versions here drifts out of date. Pin versions in your local `renv.lock` / `sessionInfo()` only if strict reproducibility is required.

---

## Key findings reported by these scripts

The figure / panel labels below are taken from the locked-in numbering in the manuscript (with track changes resolved). The supplementary figures **Fig S1, Fig S2, Figs S3 & S4, and Fig S5** below are final; the main-figure numbers (Fig 2 / Fig 3 / Fig 5 / Fig 6 / Fig 7 / Fig 8 / Fig 9) follow the layout described in the manuscript `手稿.docx` and should be re-verified if the manuscript is re-numbered.

| Output                                                       | Script(s)  | Manuscript Figure / Table                          |
| ------------------------------------------------------------ | ---------- | -------------------------------------------------- |
| Batch-corrected expression matrix                            | `01`       | — (upstream input)                                 |
| Differential expression (volcano / heatmap / PCA / Venn with ferro-aging gene set) | `03`       | Fig 2 (panels 2A–2F)                               |
| Candidate-gene selection via SVM-RFE & Random Forest         | `04`, `05` | Fig 3A                                             |
| Five candidate core genes (`Intersect_Genes.csv`)            | `06`       | Abstract; Methods §2                               |
| Diagnostic model — training-cohort ROC, nomogram, calibration, DCA | `07`       | Fig 3B–3E                                          |
| **External validation ROC on GSE53757 (AUC = 0.978)**        | `08`       | **Fig S1**                                         |
| CIBERSORT 22 immune-cell fractions + cross-correlations      | `09`       | main-text immune figure; full matrix in **Fig S2** |
| MCP-counter immune / stromal abundance + cross-correlations  | `10`       | main-text immune figure; full matrix in **Fig S2** |
| Inflammatory-cytokine correlation (full matrix + coefficients) | `11`       | main-text immune figure; full matrix in **Fig S2** |
| Consensus clustering (C1 / C2)                               | `12`       | Fig 5A–5B                                          |
| Clinical-feature heatmap & stacked plots                     | `13`, `14` | Fig 5C                                             |
| Cluster-wise immune / cytokine differences                   | `15`, `16` | Fig 5H–5J                                          |
| WGCNA modules & module–trait heatmap                         | `17`       | Fig 6                                              |
| GO / KEGG enrichment                                         | `18`       | Fig 6                                              |
| Single-cell QC / PCA / UMAP / cell-type annotation           | `19`       | Fig 7                                              |
| **Single-cell cellular origin of the five candidate core genes (cell-type distribution & module score)** | `19`       | **Figs S3, S4**                                    |
| scTenifoldKnk in-silico knockout                             | `19`       | Fig 8                                              |
| Experimental validation (qRT-PCR, CCK-8, Transwell)          | (wet-lab)  | Fig 9                                              |
| **BODIPY-C11 lipid-peroxidation flow cytometry (LOX knockdown)** | (wet-lab)  | **Fig S5**                                         |

---

## Notes on portability and reproducibility

- **No absolute paths.** All paths (input CSVs, output folders) are relative to the working directory you set before running each script.

- **Random seeds** are fixed where applicable (e.g. `set.seed(12345)` in SVM-RFE, `set.seed(2025)` / `set.seed(2026)` in random-forest and consensus clustering).

- **Package versions** may influence figure appearance; the analysis was developed with R 4.2.2 on Windows.

- If you run from the **command line**, `cd` into the repository folder first. All R scripts live in the `R/` sub-folder:

  ```bash
  cd /path/to/this/repository
  Rscript R/01_data_merge_batch_correction.R
  Rscript R/03_differential_expression_analysis.R
  # … and so on
  ```

---

## Citation

If you use this code, please cite the accompanying manuscript.

## License

This project is licensed under the MIT License. See the [`LICENSE`](./LICENSE) file for details.
