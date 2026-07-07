# gene_level_differential_expression_analysis_with_DESeq2

A self-directed bioinformatics project implementing the full DESeq2 differential
gene expression (DGE) workflow, adapted from the Harvard Chan Bioinformatics Core (HBC) "[Introduction to DGE](https://hbctraining.github.io/Intro-to-DGE/schedule/links-to-lessons.html)" self-learning workshop.

This project uses the Mov10 RNA-seq dataset (control vs. MOV10 knockdown vs.
MOV10 overexpression) to demonstrate an end-to-end bulk RNA-seq DE analysis in R:
from raw Salmon quantification files through statistical testing, visualization,
and functional (GO/KEGG) enrichment analysis.

## Experimental design (Mov10 dataset)

- control: baseline condition (3 replicates)
- MOV10_knockdown: MOV10 gene expression knocked down (2 replicates)
- MOV10_overexpression: MOV10 gene over-expressed (3 replicates)

Two pairwise contrast are tested against the control: overexpression vs. control and knockdown vs. control.

## Repository files and folders
1. README.md (this file)
2. DE_analysis_script.R (fully-annotated DESeq2 analysis script)
4. output (figures)
5. results (results tables)
6. kegg_pathways 

## Workflow summary:

The `DE_analysis_script.R` script is organized into the following sections:
 
1. **Package installation & setup** — install/load all required CRAN and
   Bioconductor packages; record session info for reproducibility.
2. **Project setup** — recommended RStudio project structure (`data/`, `meta/`, `results/`).
3. **Data import** — import Salmon transcript-level quantifications with
   `tximport`, summarized to gene level using the GRCh38 (Ensembl v94)
   tx2gene mapping table.
4. **Exploratory QC on raw counts** — inspect count distribution and the
   mean-variance relationship to confirm over-dispersion (justifying
   DESeq2's negative binomial model).
5. **Count normalization** — build the `DESeqDataSet`, estimate size
   factors, and generate normalized counts and a variance-stabilized
   (rlog) matrix for visualization.
6. **Sample-level QC** — Principal Component Analysis (PCA) and
   hierarchical clustering (sample correlation heatmap) to check for
   outliers, batch effects, and expected grouping by condition.
7. **Differential expression testing** — run `DESeq()` (Wald test) for two
   pairwise contrasts (overexpression vs. control, knockdown vs. control);
   inspect result tables, apply `apeglm` log2 fold-change shrinkage, and
   generate MA plots.
8. **Significant gene extraction** — filter each result table to genes
   passing the adjusted p-value (FDR) cutoff (`padj < 0.05`).
9. **Result visualization** — single-gene expression plot (MOV10), heatmap
   of all significant genes, and volcano plots (with top-10 gene labels).
10. **Alternate hypothesis testing** — Likelihood Ratio Test (LRT) as a
    complementary approach to the Wald test, useful for detecting genes
    affected by the experimental design across any group.
11. **Expression pattern clustering** — group significant LRT genes into
    clusters of shared expression trajectories using `DEGreport::degPatterns()`.
12. **Genomic annotation** — retrieve gene symbol/Entrez ID annotations from
    `AnnotationHub`/`ensembldb` and resolve duplicate/missing mappings.
13. **Functional analysis (ORA)** — GO Biological Process over-representation
    analysis with `clusterProfiler::enrichGO()`; visualized via dot plot,
    enrichment map, and category netplots.
14. **Gene Set Enrichment Analysis (GSEA)** — rank-based KEGG pathway
    enrichment with `gseKEGG()`, visualized with pathway diagrams via `pathview`.

## How to Run
 
1. Source "DE_analysis_script.R" into R or RStudio.
2. Create a new project and `meta/`, and `results/` subfolders.
3. Download the Salmon pseudocounts (one per sample, each containing
   `quant.sf`) and the `tx2gene_grch38_ens94.txt` annotation file inside project directory from [HBC training](https://hbctraining.github.io/Intro-to-DGE/lessons/01b_DGE_setup_and_overview.html).
4. Source `DE_analysis_script.R` and run it.
5. Export generated plots to `outputs/` and results tables to `results/`.
