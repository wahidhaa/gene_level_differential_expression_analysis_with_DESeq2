#######################################################################
# Gene-level Differential Expression Analysis using DESeq2

# Adapted from the Harvard Chan Bioinformatics Core (HBC) self-learning
# workshop: "Introduction to DGE" 
# https://hbctraining.github.io/Intro-to-DGE/schedule/links-to-lessons.html

# Dataset: Mov10 RNA-seq (control vs. MOV10 knockdown vs. MOV10
# overexpression), quantified with Salmon and imported with tximport.
#######################################################################

# ============================= 1. PACKAGES ==========================

## ---- Load packages ----
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

# Record R version and package versions for reproducibility
sessionInfo()

# ============================= 2. PROJECT SET-UP =====================

# File menu > New Project > New Directory > New Project > name "DE_analysis"
# This keeps all scripts, data, and outputs self-contained in one folder.

getwd()  # confirm the working directory is the project root

# Create "meta", "output", "results"  and "kegg_pathways" subfolders in the Files panel (bottom-right).
# Download the Salmon pseudocount folders (data) and the tx2gene annotation file
# into project directory before continuing.
# Note: Make sure that the unzipped "data" folder does not have a child "data"
# folder. If so, move all the contents to the parent "data" folder.

# ============================= 3. LOAD DATA ===========================

## List every sample subdirectory produced by Salmon (each ends in "salmon")
samples <- list.files(path = "./data", full.names = T, pattern = "salmon$")

## Build full file paths to each sample's quant.sf (Salmon quantification file)
files <- file.path(samples, "quant.sf")

## Strip path/suffix so each file is labelled with just its sample name
names(files) <- str_replace(samples, "./data/", "") %>%
  str_replace(".salmon", "")

## Load the transcript-to-gene mapping table (GRCh38, Ensembl release 94)
## This lets tximport collapse transcript-level counts to gene-level counts
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")

## Sanity-check the annotation table
tx2gene %>% View()

## Import and summarize transcript abundances to the gene level
## countsFromAbundance = "lengthScaledTPM" corrects for library size AND
## average transcript length, recommended before DESeq2 input
txi <- tximport(files, type = "salmon", tx2gene = tx2gene[, c("tx_id", "ensgene")],
                countsFromAbundance = "lengthScaledTPM")

## Inspect the structure of the tximport object (counts, abundance, length)
attributes(txi)

## Check the imported (un-rounded) gene-level counts
txi$counts %>% View()

## DESeq2 requires integer counts -> round the scaled counts
data <- txi$counts %>%
  round() %>%
  data.frame()

## Build the sample metadata table describing experimental groups
## Order must exactly match the column order of txi$counts
sampletype <- factor(c(rep("control", 3), rep("MOV10_knockdown", 2), rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))

# ==================== 4. EXPLORATORY DATA QC (RAW COUNTS) =============

## Inspect the count distribution of a single sample - RNA-seq counts are
## typically highly skewed with many low counts and a long right tail
ggplot(data) +
  geom_histogram(aes(x = Mov10_oe_1), stat = "bin", bins = 200) +
  xlab("Raw expression counts") +
  ylab("Number of genes")
# Save figure to "outputs/"

## Assess the mean-variance relationship within a group (replicates 6-8)
## RNA-seq count data is over-dispersed: variance > mean (departs from
## the red y=x line), which is why DESeq2 models counts with a
## negative binomial distribution rather than a Poisson distribution
mean_counts <- apply(data[, 6:8], 1, mean)
variance_counts <- apply(data[, 6:8], 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x = mean_counts, y = variance_counts)) +
  scale_y_log10(limits = c(1, 1e9)) +
  scale_x_log10(limits = c(1, 1e9)) +
  geom_abline(intercept = 0, slope = 1, color = "red")
# Save figure to "outputs/"

# ==================== 5. COUNT NORMALIZATION (DESeq2) ==================

## Confirm sample names in counts and metadata match, and are in the same order
## (DESeq2 will NOT reorder columns automatically)
all(colnames(txi$counts) %in% rownames(meta))
all(colnames(txi$counts) == rownames(meta))

## Build the DESeqDataSet object directly from the tximport output
## design = ~ sampletype tells DESeq2 which variable defines the comparison groups
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

## Inspect the raw (un-normalized) count matrix stored inside the object
View(counts(dds))

## Estimate size factors to correct for differences in sequencing depth/
## library composition across samples (median-of-ratios method)
dds <- estimateSizeFactors(dds)

## Inspect the size factor estimated per sample
sizeFactors(dds)

## Extract normalized counts (raw counts divided by size factors)
normalized_counts <- counts(dds, normalized = TRUE)

## Save the normalized count matrix for downstream use/sharing
write.table(normalized_counts, file = "results/normalized_counts.txt", sep = "\t", quote = F, col.names = NA)

## Apply a regularized-log (rlog) transformation for visualization only
## (stabilizes variance across the range of mean values; blind=TRUE ignores
## the experimental design so it does not bias exploratory QC plots)
rld <- rlog(dds, blind = TRUE)

# ============================= 6. PCA =================================

## Plot PC1 vs PC2 of the rlog-transformed data, colored by sample group
## Used to visually assess how samples cluster by experimental condition
## and to flag potential outliers or batch effects
plotPCA(rld, intgroup = "sampletype")
## Save figure to "outputs/"

# ==================== 7. HIERARCHICAL CLUSTERING =======================

## Extract the rlog-transformed expression matrix
rld_mat <- assay(rld)

## Compute pairwise sample-to-sample correlations
rld_cor <- cor(rld_mat)

## Confirm the correlation matrix's row/column names line up with metadata
head(rld_cor)
head(meta)

## Visualize the correlation matrix as a heatmap annotated by sample group
pheatmap(rld_cor, annotation = meta)
## Save figure to "outputs/".

## High correlation across all samples (no isolated/outlying samples) supports
## proceeding with the full dataset in the DE analysis

# ================ 8. DIFFERENTIAL EXPRESSION ANALYSIS ===================

## Re-create the DESeqDataSet object (fresh start using raw counts + design)
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sampletype)

## Run the complete DESeq2 pipeline: size factor estimation, dispersion
## estimation, and negative binomial GLM fitting + the Wald significance test
dds <- DESeq(dds)

# ---------------------- Results: Overexpression vs Control -------------

## Define the contrast: (variable, level of interest, reference level)
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")

## Extract the Wald test results table with an FDR (padj) threshold of 0.05
res_tableOE <- results(dds, contrast = contrast_oe, alpha = 0.05)

## Confirm object class (DESeqResults)
class(res_tableOE)

res_tableOE %>%
  data.frame() %>%
  View()

## Inspect column descriptions (baseMean, log2FoldChange, pvalue, padj, etc.)
mcols(res_tableOE, use.names = T)

## ---- Manual gene-level filtering (understanding NA results) ----

## Genes with zero counts across all samples -> filtered out automatically
res_tableOE[which(res_tableOE$baseMean == 0), ] %>%
  data.frame() %>%
  View()

## Genes flagged as extreme count outliers by Cook's distance (NA p/padj
## despite non-zero expression)
res_tableOE[which(is.na(res_tableOE$pvalue) &
                    is.na(res_tableOE$padj) &
                    res_tableOE$baseMean > 0), ] %>%
  data.frame() %>%
  View()

## Genes excluded by independent filtering for having a low mean count
## (padj NA but pvalue present)
res_tableOE[which(!is.na(res_tableOE$pvalue) &
                    is.na(res_tableOE$padj) &
                    res_tableOE$baseMean > 0), ] %>%
  data.frame() %>%
  View()

## ---- Shrink log2 fold changes for more reliable ranking/visualization ----

## Keep a copy of the unshrunken results for comparison
res_tableOE_unshrunken <- res_tableOE

## Apply apeglm shrinkage: pulls noisy, low-count-driven fold changes toward
## zero while preserving true large effects - improves gene ranking/plots
res_tableOE <- lfcShrink(dds, coef = "sampletype_MOV10_overexpression_vs_control", type = "apeglm")

## ---- MA plots: mean normalized expression vs. log2 fold change ----

## Before shrinkage: note the noisy fold changes at low mean counts
plotMA(res_tableOE_unshrunken, ylim = c(-2, 2))
## Save figure to "outputs/".

## After shrinkage: fold-change noise at low counts is reduced
plotMA(res_tableOE, ylim = c(-2, 2))
## Save figure to "outputs/".

## Summarize number of up/down-regulated genes at the chosen alpha
summary(res_tableOE, alpha = 0.05)

## ---- Extract significant genes (overexpression) ----

padj.cutoff <- 0.05  # FDR threshold for calling significance

## Convert results to a tidy tibble with gene IDs as a column
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

## Keep only genes passing the significance threshold
sigOE <- res_tableOE_tb %>%
  dplyr::filter(padj < padj.cutoff)

sigOE

# Save the full and significant-only OE results tables to "results/"
write.csv(res_tableOE_tb, "results/deseq2_results_OEvsControl_all.csv", row.names = FALSE)
write.csv(sigOE, "results/deseq2_results_OEvsControl_sig.csv", row.names = FALSE)

# ---------------------- Results: Knockdown vs Control -------------------

## Define the contrast for the knockdown comparison
contrast_kd <- c("sampletype", "MOV10_knockdown", "control")

## Extract Wald test results
res_tableKD <- results(dds, contrast = contrast_kd, alpha = 0.05)

## Keep unshrunken copy for comparison
res_tableKD_unshrunken <- res_tableKD

## Apply apeglm shrinkage for the knockdown coefficient
res_tableKD <- lfcShrink(dds, coef = "sampletype_MOV10_knockdown_vs_control", type = "apeglm")

## MA plots before/after shrinkage
plotMA(res_tableKD_unshrunken, ylim = c(-2, 2))
plotMA(res_tableKD, ylim = c(-2, 2))
## Save both figures to "outputs/"

## Summarize results
summary(res_tableKD, alpha = 0.05)

## ---- Extract significant genes (knockdown) ----
padj.cutoff <- 0.05

res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

sigKD <- res_tableKD_tb %>%
  dplyr::filter(padj < padj.cutoff)

sigKD

# Save the full and significant-only KD results tables to "results/"
write.csv(res_tableKD_tb, "results/deseq2_results_KDvsControl_all.csv", row.names = FALSE)
write.csv(sigKD, "results/deseq2_results_KDvsControl_sig.csv", row.names = FALSE)

# ============================= 9. VISUALIZING RESULTS ===================

## ---- Metadata as a tibble (for downstream joins with plotting data) ----
mov10_meta <- meta %>%
  rownames_to_column(var = "samplename") %>%
  as_tibble()

## ---- Prepare normalized counts with gene symbols attached ----

## Convert normalized counts matrix to a data frame with gene IDs as a column
normalized_counts <- counts(dds, normalized = T) %>%
  data.frame() %>%
  rownames_to_column(var = "gene")

## Build a lookup table of Ensembl ID -> gene symbol from tx2gene
grch38annot <- tx2gene %>%
  dplyr::select(ensgene, symbol) %>%
  dplyr::distinct()

## Attach gene symbols to the normalized counts table
normalized_counts <- merge(normalized_counts, grch38annot, by.x = "gene", by.y = "ensgene")

normalized_counts <- normalized_counts %>%
  as_tibble()

## ---- Plot expression of a single gene of interest (MOV10) ----

## Look up the Ensembl ID corresponding to the MOV10 gene symbol
grch38annot[grch38annot$symbol == "MOV10", "ensgene"]

## Quick built-in DESeq2 plot of normalized counts by group
plotCounts(dds, gene = "ENSG00000155363", intgroup = "sampletype")

## Same plot, but exported as a data frame so it can be customized with ggplot2
d <- plotCounts(dds, gene = "ENSG00000155363", intgroup = "sampletype", returnData = TRUE)

d %>% View()

ggplot(d, aes(x = sampletype, y = count, color = sampletype)) +
  geom_point(position = position_jitter(w = 0.1, h = 0)) +
  geom_text_repel(aes(label = rownames(d))) +
  theme_bw() +
  ggtitle("MOV10") +
  theme(plot.title = element_text(hjust = 0.5))
# Save figure to "outputs/".

## ---- Heatmap of all significant DE genes (overexpression vs control) ----

## Subset normalized counts to the OE + control samples and significant genes only
norm_OEsig <- normalized_counts[, c(1:4, 7:9)] %>%
  dplyr::filter(gene %in% sigOE$gene)

## Choose a sequential color palette for expression intensity
heat_colors <- brewer.pal(6, "YlOrRd")

## Plot heatmap; scale="row" z-scores each gene across samples so patterns
## of relative up/down-regulation are visible regardless of absolute expression
pheatmap(norm_OEsig[2:7],
         color = heat_colors,
         cluster_rows = T,
         show_rownames = F,
         annotation = meta,
         border_color = NA,
         fontsize = 10,
         scale = "row",
         fontsize_row = 10,
         height = 20)
# Save figure to "outputs/"

## ---- Volcano plot (overexpression vs control) ----

## Flag genes passing both a significance and effect-size threshold
## (padj < 0.05 and at least a ~1.5-fold change, log2(1.5) ~ 0.58)
res_tableOE_tb <- res_tableOE_tb %>%
  dplyr::mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

ggplot(res_tableOE_tb) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
# Save figure to "outputs".

## ---- Label the top 10 most significant genes on the volcano plot ----

## Attach gene symbols to the results tibble
res_tableOE_tb <- bind_cols(res_tableOE_tb, symbol = grch38annot$symbol[match(res_tableOE_tb$gene, grch38annot$ensgene)])

## Create an empty label column (only top genes will be filled in)
res_tableOE_tb <- res_tableOE_tb %>% dplyr::mutate(genelabels = "")

## Sort by adjusted p-value so the most significant genes come first
res_tableOE_tb <- res_tableOE_tb %>% dplyr::arrange(padj)

## Populate labels for just the top 10 genes
res_tableOE_tb$genelabels[1:10] <- as.character(res_tableOE_tb$symbol[1:10])

View(res_tableOE_tb)

## Volcano plot with gene symbol labels for the top hits
ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold_OE)) +
  geom_text_repel(aes(label = genelabels)) +
  ggtitle("Mov10 overexpression") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
# Save figure to "outputs/".

# ==================== 10. ALTERNATE HYPOTHESIS TEST: LRT ================

## Likelihood Ratio Test (LRT): compares the full model (~ sampletype) to a
## reduced model (~ 1, intercept only) to find genes where sampletype
## explains significant variation across ANY group (not just one pairwise
## contrast at a time, unlike the Wald test used above)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)

## Extract LRT results
res_LRT <- results(dds_lrt)

res_LRT

## ---- Identify significant genes from the LRT ----

res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()

sigLRT_genes <- res_LRT_tb %>%
  dplyr::filter(padj < padj.cutoff)

## Number of genes significant by the LRT
nrow(sigLRT_genes)

## Compare against the Wald test significant gene counts
nrow(sigOE)
nrow(sigKD)

# ============= 11. GENE CLUSTERING WITH SHARED EXPRESSION PATTERNS =======

## Subset to the top 1000 most significant LRT genes for faster clustering
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n = 1000)

## Pull the rlog expression values for just these genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

## Use DEGreport to group genes into clusters of similar expression
## trajectories across sample groups
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col = NULL)
## Save figure to "outputs/"

## e.g. one useful cluster pattern to look for: genes with decreased
## expression in knockdown samples and increased expression in
## overexpression samples (consistent with direct MOV10 regulation)

class(clusters)
names(clusters)
head(clusters$df)

## Extract genes belonging to cluster 1 specifically
group1 <- clusters$df %>%
  dplyr::filter(cluster == 1)

# ========================= 12. GENOMIC ANNOTATIONS =======================

## Connect to AnnotationHub (large repository of curated genomic resources)
ah <- AnnotationHub()

## Search AnnotationHub for a human Ensembl annotation database
human_ens <- query(ah, c("Homo sapiens", "EnsDb"))

## Select a specific annotation version/record by its AnnotationHub ID
human_ens <- human_ens[["AH75011"]]

## Extract gene-level annotation info (swap "genes" for "transcripts"/"exons"
## to retrieve those levels instead)
genes(human_ens, return.type = "data.frame") %>% View()

## Build a clean gene-level annotation data frame limited to genes present
## in the OE results table
annotations_ahb <- genes(human_ens, return.type = "data.frame") %>%
  dplyr::select(gene_id, gene_name, entrezid, gene_biotype) %>%
  dplyr::filter(gene_id %in% res_tableOE_tb$gene)

## entrezid is a list-column because some genes map to multiple Entrez IDs -
## keep only the first ID per gene for a 1:1 mapping
class(annotations_ahb$entrezid)
which(map(annotations_ahb$entrezid, length) > 1)

annotations_ahb$entrezid <- map(annotations_ahb$entrezid, 1) %>% unlist()

## Check annotation completeness/uniqueness
which(is.na(annotations_ahb$gene_name)) %>% length()      # genes missing a symbol
which(duplicated(annotations_ahb$gene_name)) %>% length()  # duplicated symbols

## Identify rows to keep (first occurrence of each gene symbol)
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)

annotations_ahb %>% nrow()

## Remove duplicate gene symbol rows
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

annotations_ahb %>% nrow()

## Check remaining missing Entrez IDs (these genes can't be used in
## Entrez-ID-based analyses like KEGG GSEA later on)
which(is.na(annotations_ahb$entrezid)) %>% length()

# ========================= 13. FUNCTIONAL ANALYSIS ========================

## ---- Over-representation analysis (ORA) of GO Biological Process terms ----

## Keep only genes that were actually tested (i.e. padj is not NA)
res_tableOE_tb_noNAs <- dplyr::filter(res_tableOE_tb, padj != "NA")

## Merge in gene symbol/Entrez annotations
res_ids <- left_join(res_tableOE_tb_noNAs, annotations_ahb, by = c("gene" = "gene_id"))

## Background gene universe = all genes tested in the DE analysis
allOE_genes <- as.character(res_ids$gene)

## Gene list of interest = significant DE genes only
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE_genes <- as.character(sigOE$gene)

## Run GO enrichment (hypergeometric test) restricted to Biological Process (BP)
ego <- enrichGO(gene = sigOE_genes,
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

## Export GO enrichment results to a table
cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")

## ---- Visualize ORA results ----

## Dot plot of the top enriched GO terms (size = gene count, color = padj)
dotplot(ego, showCategory = 50)
## Save figure to "outputs/".

## Compute term-term similarity needed for the enrichment map plot
ego <- enrichplot::pairwise_termsim(ego)

## Enrichment map: clusters GO terms by shared genes/similarity
emapplot(ego, showCategory = 50)
## Save figure to "outputs/"

## ---- Category netplot: links genes to the GO terms they belong to ----

## Extract log2FC for the significant genes, named by Ensembl ID
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene

cnetplot(ego,
         showCategory = 5,
         foldChange = OE_foldchanges, # new ver.: foldChange nests inside color.params
         )                            # new ver.: vertex.label.font deprecated
# Save figure to "outputs/".

## Cap fold-change color scale at +/-2 so extreme outliers don't wash out
## the color contrast for the rest of the genes
OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)

cnetplot(ego,
         showCategory = 5,
         foldChange = OE_foldchanges,
         )
## Save figure to "outputs/"

## Subset to specific GO terms of interest (not necessarily the top 5 by rank)
ego2 <- ego
ego2@result <- ego@result[c(1, 3, 4, 8, 9), ]

cnetplot(ego2,
         foldChange = OE_foldchanges,
         categorySizeBy = ~pvalue,
         showCategory = 5)

# ================ 14. GENE SET ENRICHMENT ANALYSIS (GSEA) =================
# Unlike ORA above (which needs a hard significance cutoff), GSEA uses the
# full ranked gene list, so it can detect coordinated but individually
# sub-threshold shifts in a pathway.

## ---- Prepare the ranked gene list (by Entrez ID) ----

## Remove genes without an Entrez ID (required for KEGG)
res_entrez <- dplyr::filter(res_ids, entrezid != "NA")

## Remove duplicate Entrez IDs
res_entrez <- res_entrez[which(duplicated(res_entrez$entrezid) == F), ]

## Extract fold changes and name them by Entrez ID
foldchanges <- res_entrez$log2FoldChange
names(foldchanges) <- res_entrez$entrezid

## GSEA requires the gene list sorted in decreasing order of fold change
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

## ---- Run GSEA against KEGG pathways ----

set.seed(123456)  # GSEA uses permutation testing - fix seed for reproducibility

gseaKEGG <- gseKEGG(geneList = foldchanges,  # ranked, named fold-change vector
                    organism = "hsa",        # Homo sapiens KEGG code
                    minGSSize = 20,          # minimum pathway gene-set size to test
                    pvalueCutoff = 0.05,     # significance threshold (padj)
                    verbose = FALSE)

## Extract and save GSEA results table
gseaKEGG_results <- gseaKEGG@result

write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote = F)

View(gseaKEGG_results)

## Enrichment plot for a specific KEGG pathway of interest
gseaplot(gseaKEGG, geneSetID = 'hsa03008')
## Save figure to "outputs/"

## ---- Overlay fold changes onto KEGG pathway diagrams with pathview ----

## dplyr's select() masks pathview's internal use of the same name -
## detach dplyr before running pathview to avoid conflicts
detach("package:dplyr", unload = TRUE)

## Create a destination folder for pathway diagrams
dir.create("outputs/kegg_pathways", recursive = TRUE)

## Save the original working directory so we can return to it afterward
original_wd <- getwd()

## Render a single KEGG pathway diagram colored by fold change
setwd("outputs/kegg_pathways")
pathview(gene.data = foldchanges,
         pathway.id = "hsa03008",
         species = "hsa",
         limit = list(gene = 2,  # caps the color scale at +/-2 log2FC
                      cpd = 1))
setwd(original_wd)
## pathview() writes its .png/.xml/.png files straight to the
## current working directory and does not take an output-path argument,
## so set the working directory to the target folder before each call

## Batch-render pathway diagrams for every significant KEGG pathway from GSEA
## while skipping large pathways or non-standard image layout.
get_kegg_plots <- function(x) {
  pid <- gseaKEGG_results$ID[x]
  tryCatch({
    setwd("outputs/kegg_pathways")
    pathview(gene.data = foldchanges,
             pathway.id = gseaKEGG_results$ID[x],
             species = "hsa",
             limit = list(gene = 2, cpd = 1))
    setwd(original_wd)
  }, error = function(e) {
    setwd(original_wd)
    message("Skipping pathway ", gseaKEGG_results$ID[x], ": ", conditionMessage(e))
    NULL
  })
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)
# ============================== END OF SCRIPT ==============================
