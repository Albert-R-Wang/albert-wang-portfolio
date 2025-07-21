# ------------------------------------------------------------------------------
# Title:       DESeq2 analysis pipeline with volcano plot and Gene Ontology (GO) analysis
# Author:      Albert Wang
# Last updated:        2025-06-15
# ------------------------------------------------------------------------------


### Summary
# This R script implements a complete RNA-seq differential expression analysis pipeline 
# using the DESeq2 package (https://doi.org/10.1186/s13059-014-0550-8), 
# including data preprocessing, statistical testing, result visualization, and biological interpretation.
#
# It begins by importing raw count data and sample metadata, 
# and constructing a DESeqDataSet object using experimental group annotations. 
# The pipeline supports optional pre-filtering of low-count genes 
# to improve computational efficiency and visual clarity.
# 
# The script then performs differential expression analysis using DESeq2, 
# applying log2 fold change shrinkage via the apeglm method 
# to stabilize estimates for low-count genes while preserving statistical robustness.
# Significant differentially expressed genes (DEGs) are identified
# based on user-defined thresholds for adjusted p-value and fold change, 
# and the results are exported for downstream use.
#
# For visualization, the script generates a simple volcano plot 
# that highlights significantly upregulated and downregulated genes. 
# It also applies variance-stabilizing transformation (VST) to normalize 
# gene expression variance across samples, facilitating more reliable visualization and clustering.
#
# Finally, Gene Ontology (GO) enrichment analysis is performed on the set of DEGs
# using the clusterProfiler package (https://doi.org/10.1089/omi.2011.0118), 
# and the results are visualized via dot plots, enrichment trees, and category network plots 
# to reveal functional themes and gene-pathway relationships.



###=============================================================================
### Required packages

library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

## DESeq2 package and references
# https://bioconductor.org/packages/release/bioc/html/DESeq2.html
# https://github.com/thelovelab/DESeq2
# https://doi.org/10.1186/s13059-014-0550-8

## clusterProfiler package and references
# https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# https://github.com/YuLab-SMU/clusterProfiler
# https://doi.org/10.1089/omi.2011.0118


### ============================================================================
### DESeq2 Analysis Pipeline

# Sections:
# 1. Data Import
# 2. DESeq2 Dataset Setup and Preprocessing
# 3. Differential Expression Analysis
# 4. Identify Differentially Expressed Genes (DEGs)
# 5. Variance Stabilizing Transformation (VST)
# 6. Visualization
# 7. GO Enrichment Analysis



### ============================================================================
### 1. Data Import

## Loads raw count data and sample information files.

# Import count matrix (e.g., DESeq2_example_counts.txt)
countFile <- file.choose()
cts <- read.table(countFile,header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)

# Import sample metadata (e.g., DESeq2_sample_info.txt)
annFile <- file.choose()
coldata <- read.table(annFile,header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)


## Ensure sample names match between count data and metadata
# Optionally clean sample names if needed
# rownames(coldata) <- sub("fb", "", rownames(coldata))

# Check name consistency
all(rownames(coldata) %in% colnames(cts)) 

# Reorder count matrix to match sample annotation
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) 


### ============================================================================
### 2. DESeq2 Dataset Setup and Preprocessing

## Create DESeq2 dataset and perform optional pre-filtering.

# Construct DESeq2 dataset using the "Group" column from sample metadata
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Group) #specify a column called "Group" in the metadata file. For example, "treated" and "untreated"
dds

# Optional: attach additional gene-level metadata
# featureData <- data.frame(gene = rownames(cts))
# mcols(dds) <- DataFrame(mcols(dds), featureData)

### Pre-filtering ###
# it's not necessary to pre-filter low count genes. Two useful reasons: 1) reduce memory size of the dds data object, thus increase speed within DESeq2, and 2) improve visualizations as low counts gene are not plotted

## Option 1: filter by sum of counts
keep <- rowSums(counts(dds)) >= 10 # adjust as needed
dds <- dds[keep,]
## Option 2: filter by ensuring at least X samples with a count of 10 or more
x = ncol(cts)*0.2 # 20% of the number of samples, adjust as needed
keep <- rowSums(counts(dds) >= 10) >= X
dds <- dds[keep,]

### Set reference level ###
dds$condition <- relevel(dds$Group, ref = "Lung Primary Tumor") #specify reference group (ex. "untreated")

### ============================================================================
### 3. Differential Expression Analysis

# Runs the DESeq function to perform normalization, dispersion estimation, and model fitting.
# Applies log2 fold change (LFC) shrinkage using lfcShrink with the apeglm method.
# Exports results and normalized counts to CSV files.

dds <- DESeq(dds, betaPrior = F)
# Starting in version 1.16, betaPrior is F by default, this means that DESeq2 will not perform any shrinkage of log2FC. This is moved to the lfcShrink function so that newer shrinkage methods can be more easily apply as needed (https://support.bioconductor.org/p/95695/). Turn betaPrior on if want to match the result of earlier versions of DESEq2.
# The purpose of shrinkage of log2FC is to address the problem that lowly expressed genes tend to have high variability (https://www.biostars.org/p/340269/)
# It looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution. So the large fold changes from genes with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk. This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior. (https://support.bioconductor.org/p/77461/)

# Extract DE results using Bonferroni-adjusted p-values
res <- results(dds, pAdjustMethod = "bonferroni") #default is Benjamini and Hochberg method (BH or fdr)
res

## LFC (log fold change) shrinkage
resultsNames(dds) # Show all possible comparisons that can be used in lfcShrink as coef. Can input name directly, or use the order number

resLFC = lfcShrink(dds, coef=2, type="apeglm", res = res) # The newer shrinkage methods/estimators (apeglm and ashr) outperform the Normal prior (or betaprior) in most cases (https://support.bioconductor.org/p/115900/). The default is "apeglm"
# Unlike betaPrior, all estimators used in lfcShrink do not modify the adjusted p-value (padj). When betaPrior = TRUE, DESeq2 provides p-values associated with the shrunken log2 fold changes. In contrast, lfcShrink returns only the shrunken LFC values while preserving the original p-values from the initial results.
resLFC



#### Export results ####
# Convert results to data frame and fix gene name format (e.g., HLA.E → HLA-E)
res.DE = as.data.frame(resLFC)
rownames(res.DE) <- gsub("\\.", "-", rownames(res.DE)) # DESeq2 sets the gene name with "." instead of "-" (ex. HLA.E instead of HLA-E). Here change gene names with "." to "-" (i.e. HLA.E becomes HLA-E)

# Save DEG results and normalized counts to CSV
fileName = basename(countFile) #get file name
fileName = tools::file_path_sans_ext(fileName) #remove extension
write.csv(res.DE, 
          file=paste("DESeq2Result_", fileName, ".csv", sep=""))


nc <- counts(dds, normalized = T)
write.csv(as.data.frame(nc), 
          file=paste("DESeq2NormCount_", fileName, ".csv", sep=""))


### ============================================================================
### 4. Identify Differentially Expressed Genes (DEGs)

## Filter results to identify significantly differentially expressed genes using defined LFC and adjusted p-value thresholds.

# Set thresholds
FC <- 1       # Minimum |log2FoldChange|
adjp <- 0.05  # Maximum adjusted p-value

# Identify significant genes
sigGenes <- rownames(subset(res.DE, (abs(log2FoldChange)>=FC & padj<=adjp )))

# Extract subset of significant results
sig.res.DE = subset(res.DE, rownames(res.DE) %in% sigGenes)


### ============================================================================
### 5. Variance Stabilizing Transformation (VST)

## Apply VST to stabilize variance across samples for downstream visualization.

# Adjusts for sequencing depth and stabilizes variance across the range of mean values.
# This transformation is important for visualization and clustering, as it reduces the influence of highly variable genes and makes expression values more comparable across samples.
vsd = vst(dds, blind = F) #transform counts into log2 scale (modeled to stablize variance) for visualization

# Extract VST-transformed expression matrix
vst_mat = assay(vsd)
sig.vst_mat = subset(vst_mat, rownames(vst_mat) %in% sigGenes)


### ============================================================================
### 6. Visualization

## Generates a volcano plot highlighting upregulated, downregulated, and non-significant genes.

#### Volcano plot ####
vol = res.DE %>% filter(!is.na(padj)) # Exclude genes with NA padj

# Determine which genes are upregulated, downregulated, or not significant (NS) based on the DEG cutoff
vol$group = with(vol, ifelse(padj < adjp & log2FoldChange > FC, "Upregulated",
                             ifelse(padj < adjp & log2FoldChange < -FC, "Downregulated", "NS")))
numDEG = table(vol$group) # Count the number of genes that are upregulated, downregulated, or NS

# Plot volcano plot
ggplot(vol, aes(x=log2FoldChange, y=-log10(padj), color = group)) + geom_point(size = 2, alpha = 0.75) +
  scale_color_manual(values = c("Upregulated" = "red3", "Downregulated" = "royalblue3", "NS" = "gray20")) +
  annotate("text", x = max(vol$log2FoldChange), y = max(-log10(vol$padj), na.rm = TRUE), 
           label = paste("Upregulated:", numDEG["Upregulated"]), color = "red3", hjust = 1) +
  annotate("text", x = min(vol$log2FoldChange), y = max(-log10(vol$padj), na.rm = TRUE), 
           label = paste("Downregulated:", numDEG["Downregulated"]), color = "royalblue3", hjust = 0) +
  theme_minimal()

### ============================================================================
### 7. GO Enrichment Analysis

# Performs Gene Ontology enrichment on DEGs (Biological Process terms).
# Visualizes enrichment results using dot plots, enrichment trees, and category network plots.

### Gene Ontology (GO) analysis ###

go_enrich = enrichGO(
  gene = sigGenes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05)

# Convert result to data frame and preview
go_results = as.data.frame(go_enrich)
head(go_results)


## Dot plot
# The dot plot provides an overview of enriched terms along with their gene ratios and adjusted p-values.
dotplot(go_enrich, showCategory = 10)

## Enrichment Map
# The tree plot displays hierarchical relationships between GO terms to highlight functional clusters.
enrich_result = pairwise_termsim(go_enrich)
treeplot(enrich_result, showCategory = 15, label_format = 20)

## Network plot
# The cnetplot (category network plot) shows how individual genes are connected to multiple enriched terms, offering insights into shared functional roles.
cnetplot(go_enrich, foldChange = NULL, showCategory = 3, node_label = "all")
