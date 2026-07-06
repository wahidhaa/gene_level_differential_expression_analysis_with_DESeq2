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

# Set-up
1. Create new project in RStudio with 

CRAN packages:

    install.packages("BiocManager")
    install.packages("tidyverse")
    install.packages("RColorBrewer")
    install.packages("pheatmap")
    install.packages("ggrepel")
    install.packages("cowplot")

Bioconductor packages:

    library(BiocManager)
    install("DESeq2")
    install("clusterProfiler")
    install("DOSE")
    install("org.Hs.eg.db")
    install("pathview")
    install("DEGreport")
    install("tximport")
    install("AnnotationHub")
    install("ensembldb")
    install("apeglm")

Load packages

    library(DESeq2)
    library(tidyverse)
    library(RColorBrewer)
    library(pheatmap)
    library(ggrepel)
    library(cowplot)
    library(clusterProfiler)
    library(DEGreport)
    library(org.Hs.eg.db)
    library(DOSE)
    library(pathview)
    library(tximport)
    library(AnnotationHub)
    library(ensembldb)
    library(apeglm)




