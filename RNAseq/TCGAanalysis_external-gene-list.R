# ------------------------------------------------------------------------------
# Title:       TCGA RNA-seq Analysis Pipeline Using External Gene Signatures
# Author:      Albert Wang
# Last updated:        2025-05-22
# ------------------------------------------------------------------------------


### Summary
# This R script implements an RNA-seq analysis pipeline that integrates 
# an external lung–brain metastasis (LBM) gene signature with
# transcriptomic and clinical data from The Cancer Genome Atlas 
# Lung Adenocarcinoma (TCGA-LUAD) dataset.
#
# The workflow begins by downloading and preprocessing RNA-seq data from TCGA 
# using the TCGAbiolinks package (https://doi.org/10.1093/nar/gkv1507). 
# Formalin-fixed paraffin-embedded (FFPE) and replicate samples are removed, 
# followed by optional filtering based on clinical data and low gene expression 
# to retain informative samples and genes. Differentially expressed genes (DEGs) 
# or genes of interest (e.g. from the LBM project) are then imported 
# and intersected with the TCGA dataset for downstream analysis. 
#
# The script generates annotated heatmaps with flexible customization options, 
# including log2 fold-change labeling, clustering distance metrics, 
# and sample annotations (e.g., tissue type or tumor stage). 
# Gene and sample clusters can be extracted and correlated with clinical variables. 
# Based on sample groupings, survival analysis is performed 
# using Kaplan–Meier curves and Cox proportional hazards models. 
# All results are exported to Excel spreadsheets as documentation.
#
# Additionally, the pipeline supports single-gene analysis in TCGA 
# using various stratification methods, including median, quartile, or custom cutoffs.
# Users may also export TPM expression data stratified by tissue type and tumor stage.
#
# This pipeline enables integrative analysis of TCGA RNA-seq data and 
# metastasis-related signatures to identify clinically relevant patterns.



###=============================================================================
### Required packages

library(TCGAbiolinks)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(dplyr)
library(limma) # need org.Hs.eg.db package to convert alias name to offical gene symbols (alias2SymbolTable)
library(org.Hs.eg.db)
library(openxlsx)
library(survival)
library(survminer)
library(gtsummary)
library(rstatix)


## TCGAbiolinks package and references
# http://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
# https://github.com/BioinformaticsFMRP/TCGAbiolinks
# https://doi.org/10.1093/nar/gkv1507

# For manually checking clinical data
# https://www.cbioportal.org/study/clinicalData?id=luad_tcga_pan_can_atlas_2018

### ============================================================================
### TCGA LUAD Analysis Pipeline (using external gene list)

# Sections:
# 1. Download TCGA-LUAD data
# 2. Filter FFPE and replicate samples
# 3. Additional filtering based on clinical information (optional)
# 4. Process TPM expression data
# 5. Import and match external LBM gene list
# 6. Correlate gene expression with clinical outcomes (optional)
# 7. Generate annotated heatmaps
# 8. Cluster analysis and group extraction
# 9. Survival analysis and statistical tests
# 10. Export results to Excel
# 11. Analyze single gene of interest in TCGA (optional)
# 12. Export TPM values by tissue type and tumor stage (optional)



### ============================================================================
### 1. Download TCGA-LUAD data (Lung Adenocarcinoma)


# Query and download TCGA-LUAD RNA-seq gene expression data using STAR - Counts workflow
query.luad <- GDCquery(project = "TCGA-LUAD",
                       #legacy = F,
                       #sample.type = c("Primary Tumor", "Recurrent Tumor"), #Primary Tumor, Recurrent Tumor, Metastatic, Solid Tissue Normal
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       experimental.strategy = "RNA-Seq",
                       workflow.type = "STAR - Counts")

GDCdownload(query.luad)
data.luad <- GDCprepare(query.luad)



### ============================================================================
### 2. Filter FFPE and replicate samples

# TCGA barcode
#https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/

# Filtering of FFPE samples due to potential RNA degradation in old FFPE samples
#data.luad.orig = data.luad #save the original data (before filtering)
data.luad = data.luad[, !data.luad$is_ffpe] #remove FFPE samples

barcode.all = data.luad$barcode # all barcode ID

# Based on Broad Institute GDAC (https://gdac.broadinstitute.org/ and https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334036/FAQ), 
# When multiple barcodes exist for a given sample, the following rules are apply:
# 1) Prefer H analytes over R, and R analytes over T for RNA (H > R > T). In the case of LUAD, all samples are R
# 2) Otherwise prefer the barcode with the highest lexicographical sort value, to ensure that higher portion and/or plate numbers are selected when all other barcode fields are identical


# Input the list of barcodes that should be filtered, return the filtered barcodes (one barcode for each patient, no replicates)
tcga_replicateFilter = function(tsb, analyte_position=20, plate=c(22,25), portion=c(18,19), decreasing = TRUE){ 
  sampleID = substr(tsb, start = 1, stop = 15) # sample ID (identified at position 1 to 15 in barcode)
  dp_samples = unique(sampleID[duplicated(sampleID)]) # find replicate samples
  
  if(length(dp_samples)==0){
    message("No duplicated barcodes, return original input.")
    tsb
  }else{
    uniq_tsb = tsb[! sampleID %in% dp_samples] # unique barcodes
    dp_tsb = setdiff(tsb, uniq_tsb) # duplicate barcodes list
    
    add_tsb = c() # after filtering the replicates, to add to this temp list to combine with the unique barcodes (uniq_tsb) later
    
    
    for(x in dp_samples){
      mulaliquots = dp_tsb[substr(dp_tsb, 1, 15) == x] # for each loop, get all barcodes that associate with each patient
      analytes = substr(mulaliquots, 
                        start = analyte_position, # at position 20 in barcode
                        stop = analyte_position) # obtain all analytes for each patient
      # First filter, prefer H analytes over R, and R over T. If all sample have the same analyte, move to the next filter
      if(any(analytes == "H") & !(all(analytes == "H"))){ # if there's any analyte that is H, and not all analytes are H
        aliquot = mulaliquots[which(analytes == "H")] # find which analytes is H
        add_tsb = c(add_tsb, aliquot) # add to the temp list
        dp_tsb = setdiff(dp_tsb, mulaliquots) # remove mulaliquots from the duplicate list
      }else if(any(analytes == "R") & !(all(analytes == "R"))){
        aliquot = mulaliquots[which(analytes == "R")]
        add_tsb = c(add_tsb, aliquot)
        dp_tsb = setdiff(dp_tsb, mulaliquots)
      }else if(any(analytes == "T") & !(all(analytes == "T"))){
        aliquot = mulaliquots[which(analytes == "T")]
        add_tsb = c(add_tsb, aliquot)
        dp_tsb = setdiff(dp_tsb, mulaliquots)
      }
    }
  }
  
  
  if(length(dp_tsb) == 0){
    message("Filter barcodes successfully based on analyte")
    c(uniq_tsb, add_tsb)
  } else{ # Second filter, either with portion or plate number
    # filter according to portion number (at position 18 to 19 in barcode)
    sampleID_res = substr(dp_tsb, start=1, stop=15) #remaining sample ID
    dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
    
    for(x in dp_samples_res){
      mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
      portion_codes = substr(mulaliquots,
                             start = portion[1],
                             stop = portion[2]) # obtain all portion codes for each patient
      portion_keep = sort(portion_codes, decreasing = decreasing)[1] # sort the portion code in decreasing order (larger value first), and get the largest value (first value)
      if(!all(portion_codes == portion_keep)){ # if not all portion codes are the same
        if(length(which(portion_codes == portion_keep)) == 1){ # if the largest portion code is unique to one sample, then keep that sample barcode, remove all the other replicates
          add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
          dp_tsb = setdiff(dp_tsb, mulaliquots)
        }else{
          dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)]) # if the largest portion code is not unique to one sample, then only remove those replicates that have smaller portion codes
        }
        
      }
    }
    
    if(length(dp_tsb)==0){
      message("Filter barcodes successfully based on portion number")
      c(uniq_tsb, add_tsb)
    }else{
      # filter according to plate number (at position between 22 to 25 in barcode)
      sampleID_res = substr(dp_tsb, start=1, stop=15)
      dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
      for(x in dp_samples_res){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        plate_codes = substr(mulaliquots,
                             start = plate[1],
                             stop = plate[2]) #obtain all plate numbers for each patient
        plate_keep = sort(plate_codes, decreasing = decreasing)[1] # sort the plate number in decreasing order, and get the first value
        add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)]) # add that barcode to the temp list
        dp_tsb = setdiff(dp_tsb, mulaliquots) # remove all barcodes from this patient from the duplicate list
      }
      
      if(length(dp_tsb)==0){
        message("Filter barcodes successfully based on plate number")
        c(uniq_tsb, add_tsb)
      }else{
        message("Barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
        c(uniq_tsb, add_tsb)
      }
    }
  }
}


barcode.filt = tcga_replicateFilter(barcode.all) # use the function tcga_replicateFilter defined earlier to get the one barcode for each patient (no replicate)
data.luad = data.luad[, data.luad$barcode %in% barcode.filt] # for LUAD, there are 59 normal tissue and 518 tumor samples, total of 577 samples after filtering (as of 04/12/24)

# Number of samples after filtering out low-quality or duplicate entries
length(data.luad$barcode)

### ============================================================================
### 3. Additional filtering based on clinical information (optional, run if needed)


# Construct clinical data table and filter out tumor/normal samples without clinical info

# Extract clinical info
clinical.info = data.frame(
  ID = data.luad$barcode,
  vital_status = data.luad$vital_status, #indicate whether the patients were dead or alive at the end of study
  days_to_death = data.luad$days_to_death, #only for patients that were dead
  days_to_last_follow_up = data.luad$days_to_last_follow_up, #only for patients that were alive at the end of study
  patientID = data.luad$patient,
  FFPE = data.luad$is_ffpe,
  stage = data.luad$ajcc_pathologic_stage, # Stage I, etc.
  tissue_type = data.luad$tissue_type, # normal or tumor
  tumorType = data.luad$definition # ex. Primary solid Tumor, Recurrent Solid Tumor, etc.
)

# Create a unified "days" column by prioritizing "days_to_death" and using "days_to_last_follow_up" if "days_to_death" is missing
#https://groups.google.com/g/ucsc-cancer-genomics-browser/c/YvKnWZSsw1Q?pli=1
clinical.info = clinical.info %>% mutate(days = coalesce(days_to_death, days_to_last_follow_up))

# Merge sub-stages (e.g., IA, IB, IIA, IIB, IIIA, IIIB) into their corresponding main stages (I, II, III, IV)
clinical.info = clinical.info %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ "Stage I",
                                   stage == "Stage IA" ~ "Stage I",
                                   stage == "Stage IB" ~ "Stage I",
                                   stage == "Stage II" ~ "Stage II",
                                   stage == "Stage IIA" ~ "Stage II",
                                   stage == "Stage IIB" ~ "Stage II",
                                   stage == "Stage IIIA" ~ "Stage III",
                                   stage == "Stage IIIB" ~ "Stage III",
                                   stage == "Stage IV" ~ "Stage IV"))

clinical.info.subset = clinical.info %>%
  #filter(overall.Stage != "Stage IV") %>%     # (optional) remove certain stages
  filter(tissue_type == "Tumor") %>%           # (optional) remove Normal or Tumor samples
  filter(!is.na(days)) %>%                     # (optional) remove samples without survival data
  #filter(days != 0) %>%                       # (optional) remove samples with survival day = 0
  filter(!is.na(overall.Stage)) %>%            # remove samples without tumor stage
  filter(tumorType != "Recurrent Solid Tumor") # remove recurrent tumor

data.luad = data.luad[, data.luad$barcode %in% clinical.info.subset$ID]

# Number of total samples after filtering based on clinical info
length(data.luad$barcode)

# Number of tumor samples
length(data.luad[, data.luad$tissue_type == "Tumor"]$barcode)
# Number of normal samples
length(data.luad[, data.luad$tissue_type == "Normal"]$barcode)

extName <- rowData(data.luad)$gene_name # get the official gene name

### ============================================================================
### 4. Process TPM expression data

### TPM
# Available assays in SummarizedExperiment: 1) unstranded, 2) stranded_first, 3) stranded_second, 4) tpm_unstrand, 5) fpkm_unstrand, 6) fpkm_uq_unstrand
dataPrep.TPM.luad = TCGAanalyze_Preprocessing(object = data.luad, cor.cut = 0.6, datatype = "tpm_unstrand")
rownames(dataPrep.TPM.luad) = extName
dataFilt.TPM.luad <- TCGAanalyze_Filtering(tabDF = dataPrep.TPM.luad, method = "quantile", qnt.cut = 0.25) #quantile filtering to remove genes with low count




### ============================================================================
### 5. Import and match external LBM gene list


## Import LBM data (example file: LBM_RNAseq_DEresult_all.txt)
LBMfile <- file.choose()
LBMdata <- read.table(LBMfile,header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE) #import LBM gene list

## Find differentially expressed genes (DEGs)

DEG_filter <- function(LBM_data, comparison = "BrMvsLT", FCcut = 0.6, PVcut = 0.05) {
  # Construct dynamic column names
  logFC_col <- paste0("log2FC_", comparison)
  padj_col <- paste0("padj_", comparison)
  
  # Check that the required columns exist
  if (!all(c(logFC_col, padj_col) %in% colnames(LBM_data))) {
    stop("Required columns not found in LBM_data for comparison: ", comparison)
  }
  
  ## Get DEGs based on the log2 fold change (FCcut) and adjusted p-value (PVcut) cutoff
  deg_genes <- rownames(LBM_data)[abs(LBM_data[[logFC_col]]) > FCcut & LBM_data[[padj_col]] < PVcut]
  deg_genes <- deg_genes[!is.na(deg_genes)]
  
  # filter the data to keep only genes that are differentially expressed
  DEtable <- LBM_data[rownames(LBM_data) %in% deg_genes,]
  
  return(DEtable)
}

# Specify the main comparison group of interest. Options: "BrMvsLT" and "LTvsLN"
comp = "BrMvsLT"
LBMdata = DEG_filter(LBMdata, comparison = comp, FCcut = 0.6, PVcut = 0.05) # get only DEGs
origData = LBMdata
fileName = basename(LBMfile) # Get file name
fileName = tools::file_path_sans_ext(fileName) # Remove extension

## Convert alias gene names to official gene symbols (requires org.Hs.eg.db package)
rownames(LBMdata) = alias2SymbolTable(rownames(LBMdata), species = "Hs")

# Subset TCGA expression data using the official gene names from the LBM list
TCGA.LBMlist <- dataFilt.TPM.luad[rownames(dataFilt.TPM.luad) %in% rownames(LBMdata), ]
LBM.TCGAlist = LBMdata[rownames(LBMdata) %in% rownames(TCGA.LBMlist), ]

# Record gene name conversions and identify genes excluded from the TCGA analysis due to absence in the TCGA dataset
nameConvert = setdiff(rownames(origData), rownames(dataFilt.TPM.luad)) # find out which genes are different from TCGA list
officialName = alias2SymbolTable(nameConvert, species = "Hs") # convert to 
notInTCGA = setdiff(rownames(LBMdata), rownames(TCGA.LBMlist)) #find if any genes are not included in TCGA
nameConvertTable = data.frame(matrix(ncol = 3, nrow = length(nameConvert)))
colnames(nameConvertTable) = c("Original Name", "Official Name", "TCGA")
nameConvertTable$`Original Name` = nameConvert
nameConvertTable$`Official Name` = officialName

nameConvertTable$`TCGA`[nameConvertTable$`Official Name` %in% notInTCGA] = "Not found" # note the genes that are not found in TCGA data
nameConvertTable[is.na(nameConvertTable)] = "" #replace NA value




## Filter out low-expressing genes in the TCGA dataset (adjust thresholds if needed)
idx = rowSums(TCGA.LBMlist >= 5) >= (ncol(TCGA.LBMlist) * 0.1) # filter out genes where less than 10% of samples with expressions greater than or equal to 5
TCGA.LBMlist.orig = TCGA.LBMlist
TCGA.LBMlist = TCGA.LBMlist[idx,]
TCGA.LBMlist.filterOut = TCGA.LBMlist.orig[!idx,]
LBM.TCGAlist.filterOut = LBM.TCGAlist[rownames(LBM.TCGAlist) %in% rownames(TCGA.LBMlist.filterOut),]




## Define heatmap cluster distance ##
clusterDist = "euclidean" #options: "euclidean", "correlation", and more
#clusterDist = "correlation"

# set up annotations
AnnData <- data.frame(stage = data.luad$ajcc_pathologic_stage,
                      tissue = data.luad$sample_type) # annotation for TCGA samples
rownames(AnnData) <- data.luad$barcode


### ============================================================================
### 6. Correlate gene expression with clinical outcomes (optional, run if needed)

clinical.cor = clinical.info.subset

clinical.cor = clinical.cor %>% mutate(status = case_when(vital_status == "Alive" ~ 1,
                                                    vital_status == "Dead" ~ 2)) # add another column "status" to indicate that the patient is alive (1) or dead (2)
clinical.cor = clinical.cor %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ 1,
                                   stage == "Stage IA" ~ 1,
                                   stage == "Stage IB" ~ 1,
                                   stage == "Stage II" ~ 2,
                                   stage == "Stage IIA" ~ 2,
                                   stage == "Stage IIB" ~ 2,
                                   stage == "Stage IIIA" ~ 3,
                                   stage == "Stage IIIB" ~ 3,
                                   stage == "Stage IV" ~ 4))
clinical.cor = subset(clinical.cor, tissue_type != "Normal") # remove Normal, only interested in correlation in tumor samples



testExp = as.data.frame(t(TCGA.LBMlist))
testExp$ID = rownames(testExp)
clinical.test = merge(testExp, clinical.cor, by = "ID")
corMatrix = subset(clinical.test, select = -c(ID, vital_status, days_to_death, days_to_last_follow_up, patientID, FFPE, stage, tissue_type, tumorType, status)) # remove unwanted columns


corResult = cor(corMatrix, use = "complete.obs") # ignore NA
corResult.gene = as.data.frame(corResult[,c("overall.Stage", "days")]) # select correlation of all genes with overall.Stage and days (column data)
corResult.gene = corResult.gene[!(rownames(corResult.gene) %in% c("overall.Stage", "days")),] # remove the correlation of clinical info with itself (cor=1, row data)
corResult.gene.abs = abs(corResult.gene) # absolute value of correlation result
corResult.gene.abs = corResult.gene.abs[order(corResult.gene.abs$overall.Stage, decreasing = T),] # sort based on correlation of gene expression and stage
geneNum = nrow(corResult.gene.abs)/2 #define the number of genes to use
topCorGene = rownames(corResult.gene.abs[1:geneNum,]) # take the top x number of genes

TCGA.LBMlist.preCorFilt = TCGA.LBMlist
TCGA.LBMlist = subset(TCGA.LBMlist, rownames(TCGA.LBMlist) %in% topCorGene)

# Export result as excel file
#list_of_datasets <- list("corResult" = corResult.gene, "corResult abs" = corResult.gene.abs, "heatmap genes" = TCGA.corList)
#write.xlsx(list_of_datasets, file = paste("TCGAbiolinks_LBM_CorFilter_", fileName, "_", clusterDist, ".xlsx", sep=""), rowNames = T)



AnnData_gene = data.frame(BrMvsLT.logFC = LBM.TCGAlist$log2FC_BrMvsLT) # annotation for genes
AnnData_gene$UW.logFC.sign[LBM.TCGAlist$log2FC_BrMvsLT<0] = "Downregulated"
AnnData_gene$UW.logFC.sign[LBM.TCGAlist$log2FC_BrMvsLT>0] = "Upregulated"
rownames(AnnData_gene) = rownames(LBM.TCGAlist)
Ann_color = list("stage" = c("Stage I" = "yellow1", "Stage IA" = "yellow2", "Stage IB" = "yellow3",
                             "Stage II" = "green1", "Stage IIA" = "green2", "Stage IIB" = "green3",
                             "Stage IIIA" = "blue2", "Stage IIIB" = "blue3",
                             "Stage IV" = "red3", "NA" = "azure1"), 
                 "tissue" = c(
                   "Primary Tumor" = "cyan",
                   "Solid Tissue Normal" = "tan",
                   "Recurrent Tumor" = "cyan3"),
                 "UW.logFC.sign" = c("Downregulated" = "purple3", "Upregulated" = "orange2"))


breaksList = seq(-2, 2, by = 0.1)


TCGAheatmap = pheatmap(log2(TCGA.LBMlist+1), scale="row",cluster_row=T, cluster_col=T, 
                       clustering_distance_rows=clusterDist, clustering_distance_cols=clusterDist, 
                       clustering_method="ward.D2",color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(length(breaksList)),
                       border_color="grey10", show_rownames = F, show_colnames = F, fontsize = 11, #cellwidth = 10, cellheight = 10,
                       annotation_col = AnnData, annotation_row = AnnData_gene,breaks = breaksList, treeheight_row = 100,annotation_colors = Ann_color)

tiff(paste("TCGAbiolinks_LBM_corHeatmap_", fileName, "_", clusterDist,".tiff", sep = ""), width = 7000, height = 5000, units = "px", res = 600, compression = 'lzw')
TCGAheatmap
dev.off()


### Skip to Section 8 to continue with cluster analysis and group assignment








### ============================================================================
### 7. Generate annotated heatmaps


### Heatmap (without labeling log2FC sign from LBM data) ###

## Detailed annotation (legacy option — use only if detailed staging is required) ##
Ann_color = list("stage" = c("Stage I" = "yellow1", "Stage IA" = "yellow2", "Stage IB" = "yellow3",
                             "Stage II" = "green1", "Stage IIA" = "green2", "Stage IIB" = "green3",
                             "Stage IIIA" = "blue2", "Stage IIIB" = "blue3",
                             "Stage IV" = "red3", "NA" = "azure1"), 
                 "tissue" = c(
                   "Primary Tumor" = "cyan",
                   "Solid Tissue Normal" = "tan",
                   "Recurrent Tumor" = "cyan3"))



## Simplified annotation ##
AnnData <- data.frame(stage = data.luad$ajcc_pathologic_stage,
                      tissue = data.luad$sample_type) # annotation for TCGA samples
rownames(AnnData) <- data.luad$barcode
newAnn = AnnData %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ "Stage I",
                                   stage == "Stage IA" ~ "Stage I",
                                   stage == "Stage IB" ~ "Stage I",
                                   stage == "Stage II" ~ "Stage II",
                                   stage == "Stage IIA" ~ "Stage II",
                                   stage == "Stage IIB" ~ "Stage II",
                                   stage == "Stage IIIA" ~ "Stage III",
                                   stage == "Stage IIIB" ~ "Stage III",
                                   stage == "Stage IV" ~ "Stage IV")) #simplified labeling of stages

# Select only the stage annotation for the heatmap (tissue type is excluded for simplicity, as only tumor samples are included in the final analysis).
# Comment out this line if you want to include tissue type as well.
AnnData = newAnn %>% select("overall.Stage") 

Ann_color = list("overall.Stage" = c("Stage I" = "grey95", "Stage II" = "#d4c0d9",
                             "Stage III" = "#a982b4", "Stage IV" = "#643d6e")) # purple gradient
# Alternative color option
Ann_color = list("overall.Stage" = c("Stage I" = "grey97", "Stage II" = "#9bf09d",
                                     "Stage III" = "#19a337", "Stage IV" = "#337a2c")) # green gradient



## Generate heatmap wihtout labeling log2FC ##
breaksList = seq(-2, 2, by = 0.1)

TCGAheatmap = pheatmap(log2(TCGA.LBMlist+1), scale="row",cluster_row=T, cluster_col=T, 
                       clustering_distance_rows=clusterDist, clustering_distance_cols=clusterDist, 
                       clustering_method="ward.D2",color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(length(breaksList)),
                       border_color="grey10", show_rownames = F, show_colnames = F, fontsize = 11, #cellwidth = 10, cellheight = 10,
                       annotation_col = AnnData, breaks = breaksList, treeheight_row = 100,annotation_colors = Ann_color)


tiff(paste("TCGAbiolinks_LBM_heatmap_", fileName, "_", clusterDist,".tiff", sep = ""), width = 7000, height = 5000, units = "px", res = 600, compression = 'lzw')
TCGAheatmap
dev.off()



### Heatmap (with labeling log2FC sign from LBM data) ###
# Can only be run with data file with columns labeled with log2FC_BrMvsLT, padj_BrMvsLT, log2FC_LTvsLN, and padj_LTvsLN

## Single comparison (ex. BrMvsLT or LTvsLN) ##
create_gene_annotation <- function(LBM_data, comparison = "BrMvsLT", detail.Ann = F) {
  # Construct dynamic column names
  logFC_col <- paste0("log2FC_", comparison)
  padj_col <- paste0("padj_", comparison)
  
  # Check that the required columns exist
  if (!all(c(logFC_col, padj_col) %in% colnames(LBM_data))) {
    stop("Required columns not found in LBM_data for comparison: ", comparison)
  }
  
  # Create gene annotation dataframe
  AnnData_gene <- data.frame(logFC = LBM_data[[logFC_col]])
  AnnData_gene$UW.logFC.sign <- NA
  AnnData_gene$UW.logFC.sign[LBM_data[[logFC_col]] < 0 & LBM_data[[padj_col]] < 0.05] <- "Downregulated"
  AnnData_gene$UW.logFC.sign[LBM_data[[logFC_col]] > 0 & LBM_data[[padj_col]] < 0.05] <- "Upregulated"
  rownames(AnnData_gene) <- rownames(LBM_data)
  colnames(AnnData_gene) <- c(paste0(comparison, ".logFC"), "UW.logFC.sign")
  
  # Define annotation colors
  if (detail.Ann){
    Ann_color <- list(
      "stage" = c("Stage I" = "yellow1", "Stage IA" = "yellow2", "Stage IB" = "yellow3",
                  "Stage II" = "green1", "Stage IIA" = "green2", "Stage IIB" = "green3",
                  "Stage IIIA" = "blue2", "Stage IIIB" = "blue3",
                  "Stage IV" = "red3", "NA" = "azure1"),
      "tissue" = c(
        "Primary Tumor" = "cyan",
        "Solid Tissue Normal" = "tan",
        "Recurrent Tumor" = "cyan3"),
      "UW.logFC.sign" = c("Downregulated" = "purple3", "Upregulated" = "orange2")
    )
  } else {
    Ann_color = list("overall.Stage" = c("Stage I" = "grey95", "Stage II" = "#d4c0d9",
                                         "Stage III" = "#a982b4", "Stage IV" = "#643d6e"),
                     "tissue" = c(#"Solid Tissue Normal" = "tan",
                       #"Recurrent Tumor" = "cyan3",
                       "Primary Tumor" = "cyan"),
                     "UW.logFC.sign" = c("Downregulated" = "gray50", "Upregulated" = "goldenrod3"))
    
  }
  
  
  return(list(AnnData_gene = AnnData_gene, Ann_color = Ann_color))
}

## Set up annotation ##
result_Ann <- create_gene_annotation(LBM.TCGAlist, comparison = comp, detail.Ann = F)
AnnData_gene <- result_Ann$AnnData_gene
Ann_color <- result_Ann$Ann_color



### Alternative option: Two comparisons (ex. BrMvsLT AND LTvsLN) (skip to "Heatmap with log2FC" if only using single comparison) ###

# Select a second comparison group of interest. For example, if "BrMvsLT" was selected as the main comparison (variable `comp`), then can select "LTvsLN" as a second comparison. This is only used for heatmap annotation.
if (comp == "BrMvsLT"){
  comp2 = "LTvsLN"
} else if (comp == "LTvsLN"){
  comp2 = "BrMvsLT"
} else {
  warning("The main comparison is neither BrMvsLT nor LTvsLN; defaulting to BrMvsLT as the secondary comparison (comp2).")
  comp2 = "BrMvsLT"
}

Ann1 = create_gene_annotation(LBM.TCGAlist, comparison = comp, detail.Ann = F)
colnames(Ann1$AnnData_gene) = c(paste0(comp, ".log2FC"), paste0(comp, ".log2FC.sign"))
Ann2 = create_gene_annotation(LBM.TCGAlist, comparison = comp2, detail.Ann = F)
colnames(Ann2$AnnData_gene) = c(paste0(comp2 ,".log2FC"), paste0(comp2, ".log2FC.sign"))


AnnData_gene = cbind(Ann1$AnnData_gene, Ann2$AnnData_gene) # annotation for genes


# Set value limits for visualization (truncate extreme log2FC values for plotting purpose)
upLimit = 3
lowLimit = -3
AnnData_gene <- AnnData_gene %>%
  mutate(across(where(is.numeric), ~ pmin(pmax(., lowLimit), upLimit)))
#where(is.numeric): selects only numeric columns
#pmax(..., lowLimit): raises any value below lowLimit to lowLimit
#pmin(..., upLimit): lowers any value above upLimit to upLimit


# Remove raw log2FC columns, keeping only direction annotations to simplify the plot.
# Comment out this line if you wish to plot the actual log2FC values instead.
AnnData_gene = subset(AnnData_gene, select = -c(1,3))

# Check if AnnData_gene and LBM.TCGAlist has the same rownames
identical(rownames(AnnData_gene), rownames(LBM.TCGAlist))


Ann_color = list("overall.Stage" = c("Stage I" = "grey97", "Stage II" = "#9bf09d",
                                     "Stage III" = "#19a337", "Stage IV" = "#337a2c"),
                 "tissue" = c(#"Solid Tissue Normal" = "tan",
                   #"Recurrent Tumor" = "cyan3",
                   "Primary Tumor" = "cyan"),
                 "LTvsLN.log2FC.sign" = c("Downregulated" = "purple3", "Upregulated" = "orange2"),
                 "BrMvsLT.log2FC.sign" = c("Downregulated" = "purple3", "Upregulated" = "orange2"))



### Heatmap with log2FC ###
breaksList = seq(-2, 2, by = 0.1)


TCGAheatmap = pheatmap(log2(TCGA.LBMlist+1), scale="row",cluster_row=T, cluster_col=T, 
                       clustering_distance_rows=clusterDist, clustering_distance_cols=clusterDist, 
                       clustering_method="ward.D2",color=colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(length(breaksList)),
                       border_color="grey10", show_rownames = F, show_colnames = F, fontsize = 11, #cellwidth = 10, cellheight = 10,
                       annotation_col = AnnData, annotation_row = AnnData_gene,breaks = breaksList, treeheight_row = 100,annotation_colors = Ann_color)


## Export heatmap ##
tiff(paste("TCGAbiolinks_LBM_heatmap_", fileName,"_geneAnn_", clusterDist, ".tiff", sep = ""), width = 7000, height = 5000, units = "px", res = 600, compression = 'lzw')
TCGAheatmap
dev.off()


### ============================================================================
### 8. Cluster analysis and group extraction

### Export gene cluster order and group assignments from heatmap ###
cutHeight_g = 150 # Adjust based on the cluster

geneOrder <- rownames(TCGA.LBMlist[TCGAheatmap$tree_row[["order"]],]) # Reorder original gene data to match clustering order in the heatmap

tiff(paste("TCGAbiolinks_LBM_heatmapGeneGroup_", fileName, "_", clusterDist, "_H", cutHeight_g, ".tiff", sep = ""), width = 7000, height = 4000, units = "px", res = 600, compression = 'lzw')
plot(TCGAheatmap$tree_row) # Used to visually determine the appropriate cut height (h) for clustering
abline(h=cutHeight_g, col="red", lty=2, lwd=2)
dev.off()

geneGroup <- as.data.frame(cutree(TCGAheatmap$tree_row, h=cutHeight_g))
geneGroup_order <- geneGroup[order(match(rownames(geneGroup), geneOrder)), , drop = F]
OrigGroup_gene = unique(geneGroup_order[,1]) # find the group numbers (from column 2) in original order
for(i in 1:length(OrigGroup_gene)){
  geneGroup_order$heatmapGroup[geneGroup_order[,1]==OrigGroup_gene[i]] = i # assign group number based on the order from heatmap (order from left to right of heatmap, group 1 starts on the left)
}
geneGroup_order$geneName = rownames(geneGroup_order)
#writeClipboard(geneOrder)



### Export sample cluster order and group assignments from heatmap ###
AnnData <- data.frame(stage = data.luad$ajcc_pathologic_stage,
                      tissue = data.luad$sample_type) # annotation for TCGA samples, reset for code afterwards (need stage and tissue to work)
rownames(AnnData) <- data.luad$barcode

cutHeight_s = 170 # Adjust based on the cluster

tiff(paste("TCGAbiolinks_LBM_heatmapSampleGroup_", fileName, "_", clusterDist, "_H", cutHeight_s, ".tiff", sep = ""), width = 7000, height = 4000, units = "px", res = 600, compression = 'lzw')
plot(TCGAheatmap$tree_col)
abline(h=cutHeight_s, col="red", lty=2, lwd=2)
dev.off()

sampleGroup <- as.data.frame(cutree(TCGAheatmap$tree_col, h=cutHeight_s))
sampleOrder <- colnames(TCGA.LBMlist[,TCGAheatmap$tree_col[["order"]]])
#writeClipboard(sampleOrder)
sampleGroup_order <- sampleGroup[order(match(rownames(sampleGroup), sampleOrder)), ,drop = F] #match and reorder the sampleGroup with the order of sampleOrder
#sampleGroup_Ann <- AnnData[match(rownames(sampleGroup_order),rownames(AnnData)),]
sampleGroup_all <- merge(sampleGroup_order, AnnData, by = "row.names") #merge AnnData and sampleGroup based on row name
sampleGroup_all <- sampleGroup_all[match(rownames(sampleGroup_order),sampleGroup_all$Row.names),] #reorder based on the sampleOrder
OrigGroup_sample = unique(sampleGroup_all[,2]) # find the group numbers (from column 2) in original order
for(i in 1:length(OrigGroup_sample)){
  sampleGroup_all$heatmapGroup[sampleGroup_all[,2]==OrigGroup_sample[i]] = i # assign group number based on the order from heatmap (order from left to right of heatmap, group 1 starts on the left)
}





### ============================================================================
### 9. Survival analysis and statistical tests



### Survival test and Fisher's exact test ###

# https://docs.gdc.cancer.gov/Data_Dictionary/viewer/ (description of the clinical information)
clinical = data.frame(
  ID = data.luad$barcode,
  vital_status = data.luad$vital_status, #indicate whether the patients were dead or alive at the end of study
  days_to_death = data.luad$days_to_death, #only for patients that were dead
  days_to_last_follow_up = data.luad$days_to_last_follow_up, #only for patients that were alive at the end of study
  patientID = data.luad$patient,
  FFPE = data.luad$is_ffpe
  )


colnames(sampleGroup_all)[1] = "ID" # Rename the first column

clinical_all = merge(clinical, sampleGroup_all, by="ID") # Merge clinical information with sample group annotations

clinical_all_surv = clinical_all %>% mutate(days = coalesce(days_to_death, days_to_last_follow_up)) # Combine "days_to_death" and "days_to_last_follow_up" to one column named "days"

clinical_all_surv = clinical_all_surv %>% mutate(status = case_when(vital_status == "Alive" ~ 1,
                                        vital_status == "Dead" ~ 2)) # Add another column "status" to indicate that the patient is alive (1) or dead (2)




clinical_all_surv = clinical_all_surv %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ "Stage I or II",
                                   stage == "Stage IA" ~ "Stage I or II",
                                   stage == "Stage IB" ~ "Stage I or II",
                                   stage == "Stage II" ~ "Stage I or II",
                                   stage == "Stage IIA" ~ "Stage I or II",
                                   stage == "Stage IIB" ~ "Stage I or II",
                                   stage == "Stage IIIA" ~ "Stage III or IV",
                                   stage == "Stage IIIB" ~ "Stage III or IV",
                                   stage == "Stage IV" ~ "Stage III or IV")) # consider stage 1 and 2 as early stage, while stage 3 and 4 as late stage


tumorGroup = data.frame(matrix(nrow = length(unique(clinical_all_surv$heatmapGroup))-1, ncol = 1))
colnames(tumorGroup) = c("heatmapGroup")
m=1
for (n in sort(unique(clinical_all_surv$heatmapGroup))) { # record the heatmap groups that doesn't have normal tissue as the majority
  temp = subset(clinical_all_surv, heatmapGroup %in% n)
  if (tail(names(sort(table(temp$tissue))),1) != "Solid Tissue Normal"){
    tumorGroup[m,] = n
    m = m+1
  }
}
clinical_all_surv_noNormal = subset(clinical_all_surv, heatmapGroup %in% tumorGroup$heatmapGroup) # if no group with normal tissue as the majority, then it's not a good cluster because it cannot distinguish normal and tumor tissues

df = clinical_all_surv_noNormal[,c("heatmapGroup", "overall.Stage")]


## Fisher's exact test ##
padjMethod = "holm"
fisherResult = pairwise_fisher_test(table(df), detailed = T, p.adjust.method = padjMethod)
fisherResult$padj.method = padjMethod




### Pairwise comparisons ###
nLoop = 1
fisherResult[,c("Total # sample (1)", "# stage I/II (1)", "# stage III/IV (1)", "% stage III/IV (1)", "Total # sample (2)", "# stage I/II (2)", "# stage III/IV (2)", "% stage III/IV (2)", "logrank pval", "HazardRatio", "95 upper", "95 lower", "HR pvalue")] = NA # add empty columns to store values in the for loop
for (first in sort(unique(clinical_all_surv_noNormal$heatmapGroup))){ # sort the group number alphabetically
  for (second in sort(unique(clinical_all_surv_noNormal$heatmapGroup))){
    if (first < second) { # To rule out comparisons that have the same group number (ex. 1 vs 1) and smaller group number in the second comparison because those are repeats (ex. rule out 2 vs 1 because 1 vs 2 has already been done in previous loop)
      
      
      selectGroups = c(as.character(first), as.character(second)) # choose two groups to compare (based on the numbering in heatmapGroup). Smaller number has to go first, otherwise the labeling in the survival plot would be switched.
      selectGroups_export = paste(selectGroups, collapse = "-") # ex. "3" "4" becomes "3-4". This is only for naming the export tiff file.
      clinical_all_surv_subset = subset(clinical_all_surv_noNormal, heatmapGroup %in% as.numeric(selectGroups)) # select only certain groups
      
      
      # Fit survival data using the Kaplan-Meier method
      surv_patient <- Surv(time = clinical_all_surv_subset$days, event = clinical_all_surv_subset$status) # A '+' behind survival time indicates censored data points (i.e. alive vs dead)
      
      fit <- survfit(surv_patient ~ clinical_all_surv_subset$heatmapGroup, data = clinical_all_surv_subset) # calculate survival fit based on groups
      
      #summary(fit)
      surv_table <- surv_summary(fit, data = clinical_all_surv_subset) # create a summary table
      surv_med_result <- surv_median(fit,combine=T) # calculate survival median
      
      #get log-rank p-value
      fisherResult$`logrank pval`[nLoop] = surv_pvalue(fit)$pval
      
      ## Cox Proportional Hazards Model ##
      fit.coxph = coxph(surv_patient ~ heatmapGroup, data = clinical_all_surv_subset)
      sum.HR = summary(fit.coxph)
      fisherResult$HazardRatio[nLoop] = sum.HR$conf.int[1]
      fisherResult$`95 upper`[nLoop] = sum.HR$conf.int[4]
      fisherResult$`95 lower`[nLoop] = sum.HR$conf.int[3]
      fisherResult$`HR pvalue`[nLoop] = sum.HR$waldtest[3]
      table.HR = fit.coxph%>% tbl_regression(exp = T, label = heatmapGroup~paste("Group", as.character(second), " vs. Group", as.character(first), sep="")) # https://www.danieldsjoberg.com/gtsummary/
      gt::gtsave(as_gt(table.HR), file = paste("TCGAbiolinks_LBM_HazardRatio_", fileName, "_group", selectGroups_export, ".png", sep = "")) # save table as screenshot
      #ggforest(fit.coxph, data = clinical_all_surv_subset) # Good for comparing how different factors (ex. sex, age, etc.)
      
      
      
      
      ### Record number of samples in each stage ###
      
      
      totalCount = table(clinical_all_surv_subset$heatmapGroup) # find out the total number of samples in each group
      clinical_all_surv_subset1 = subset(clinical_all_surv_subset, heatmapGroup %in% as.numeric(selectGroups[1])) # get only group 1 samples
      stageCount1 = table(clinical_all_surv_subset1$overall.Stage) # find the number of early stage and late stage samples in group 1
      percentLate1 = stageCount1[2]/totalCount[1] # find the percentage of late stage samples (stage III and IV) in group 1. NOTE: samples with no stage label will be counted towards total number of samples, but not in either early or late stage samples
      clinical_all_surv_subset2 = subset(clinical_all_surv_subset, heatmapGroup %in% as.numeric(selectGroups[2])) # get only group 2 samples
      stageCount2 = table(clinical_all_surv_subset2$overall.Stage) # find the number of early stage and late stage samples in group 2
      percentLate2 = stageCount2[2]/totalCount[2] # find the percentage of late stage samples (stage III and IV) in group 2
      
      
      fisherResult$`Total # sample (1)`[nLoop] = totalCount[1]
      fisherResult$`# stage I/II (1)`[nLoop] = stageCount1[1]
      fisherResult$`# stage III/IV (1)`[nLoop] = stageCount1[2]
      fisherResult$`% stage III/IV (1)`[nLoop] = percentLate1
      fisherResult$`Total # sample (2)`[nLoop] = totalCount[2]
      fisherResult$`# stage I/II (2)`[nLoop] = stageCount2[1]
      fisherResult$`# stage III/IV (2)`[nLoop] = stageCount2[2]
      fisherResult$`% stage III/IV (2)`[nLoop] = percentLate2
      
      
      
      
      surv_plot <- ggsurvplot(fit, 
                              pval = T,                 # show p-value
                              conf.int = F,             # show confidence intervals
                              #conf.int.style = "step",  # customize style of confidence intervals
                              xlab = "Time since diagnosis (year)",   # X axis label
                              xscale = "d_y", # transform time labels from one unit to another. For example, "d_y" transform labels from days to years. 
                              break.time.by = 730.5,       # specify X axis in time intervals. (note: 1 year = 365.25 days, 2 years = 730.5 days, etc.)
                              ggtheme = theme_classic(), # customize plot and risk table with a theme.
                              risk.table = "abs_pct",  # absolute number and percentage at risk.
                              risk.table.y.text.col = T,# colour risk table text annotations.
                              risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table.
                              #ncensor.plot = TRUE,      # plot the number of censored subjects at time t
                              #pval.method = T,        # display the test name used for calculating the p-value 
                              surv.median.line = "none",  # add the median survival pointer.
                              legend.labs = paste("Group ", selectGroups, " (n=", totalCount, ")", sep = ""),    # change legend labels.
                              palette = 
                                c("grey35", "cyan3"), # custom color palettes.
                              legend.title="",          #remove legend title (Strata)
                              pval.size = 9,           # p-value font size
                              size =3,                 # survival line width
                              #font.tickslab = 15,     # change font size of tick labels. Can also change color and format with font.tickslab = c(14, "bold.italic", "darkred")
                              #font.x = 15,
                              #font.y = 15,
                              fontsize = 5.5,      # change font size in the risk table
                              surv.scale="percent"    # change y-axis scale to percentage
      )
      # adjust text font size on the plot
      surv_plot$plot <- surv_plot$plot + 
        theme(legend.title = element_text(size = 30, color = "black"), #can add face = "bold" to change the font format
              legend.text = element_text(size = 27, color = "black"),
              axis.text.x = element_text(size = 30, color = "black"),
              axis.text.y = element_text(size = 30, color = "black"),
              axis.title.x = element_text(size = 30, color = "black"),
              axis.title.y = element_text(size = 30, color = "black"))
      surv_plot$plot$theme$line$size <- 1.5 # change line width in the plot
      
      
      surv_plot$table$theme$axis.title.x$size <- 30 # risk table title font size
      surv_plot$table$theme$line$size <-1.5  # risk table line width
      surv_plot$table$theme$axis.text.x <- element_text(size = 30, color = "black")  # risk table x tick labels font size
      surv_plot$table$theme$text <- element_text(size = 30, color = "black")  # font size of "Number at risk"
      
      
      tiff(paste("TCGAbiolinks_LBM_survival_", fileName, "_group", selectGroups_export, ".tiff", sep = ""), width = 6000, height = 6000, units = "px", res = 600, compression = 'lzw')
      print(surv_plot)
      dev.off()
      
      
      
      
      nLoop = nLoop + 1
    }
    
  }
  
}


### ============================================================================
### 10. Export results to Excel

# Move row names to the first column named "gene" (write.xlsx doesn't export row names)
TCGA.LBMlist.export = as.data.frame(TCGA.LBMlist)
TCGA.LBMlist.export = tibble::rownames_to_column(TCGA.LBMlist.export, "gene") 
TCGA.LBMlist.filterOut.export = as.data.frame(TCGA.LBMlist.filterOut)
TCGA.LBMlist.filterOut.export = tibble::rownames_to_column(TCGA.LBMlist.filterOut.export, "gene")
LBM.TCGAlist.export = as.data.frame(LBM.TCGAlist)
LBM.TCGAlist.export = tibble::rownames_to_column(LBM.TCGAlist.export, "gene") 

# Export renamed genes, filtered-out genes, TCGA expression data for heatmap, gene/sample clusters, and Fisher test results
list_of_datasets <- list("SampleGroup" = sampleGroup_all, "GeneGroup" = geneGroup_order, "FisherTest & HR" = fisherResult, "Clinical" = clinical_all_surv,
                         "TCGA LBMlist" = TCGA.LBMlist.export, "low count removed" = TCGA.LBMlist.filterOut.export, "LBM expression" = LBM.TCGAlist.export, "GeneNameConverted" = nameConvertTable)
write.xlsx(list_of_datasets, file = paste("TCGAbiolinks_LBM_compiledResult_", fileName, ".xlsx", sep=""))



## run the following if correlation filter was used (skip the export above)
list_of_datasets <- list("SampleGroup" = sampleGroup_all, "GeneGroup" = geneGroup_order, "FisherTest & HR" = fisherResult, "Clinical" = clinical_all_surv,
                         "TCGA LBMlist" = TCGA.LBMlist, "low count removed" = TCGA.LBMlist.filterOut, "LBM expression" = LBM.TCGAlist, "GeneNameConverted" = nameConvertTable,
                         "corResult" = corResult.gene, "corResult abs" = corResult.gene.abs, "preCorFilter" = TCGA.LBMlist.preCorFilt)
write.xlsx(list_of_datasets, file = paste("TCGAbiolinks_LBM_compiledResult_", fileName, "_corFilter.xlsx", sep=""), rowNames = T)



### ============================================================================
# 11. Analyze single gene of interest in TCGA (optional)


### Select gene of interest for expression and survival (Run only if needed) ###

# For gene of interest, only interest in tumor samples' expression
data.luad.tumor = data.luad[, data.luad$tissue_type == "Tumor"] # get only tumor samples
dataPrep.TPM.luad.tumor = TCGAanalyze_Preprocessing(object = data.luad.tumor, cor.cut = 0.6, datatype = "tpm_unstrand")
rownames(dataPrep.TPM.luad.tumor) = extName
dataFilt.TPM.luad.tumor <- TCGAanalyze_Filtering(tabDF = dataPrep.TPM.luad.tumor, method = "quantile", qnt.cut = 0.25) #quantile filtering to remove genes with low count
TCGA.tumor = as.data.frame(dataFilt.TPM.luad.tumor)

# Specify gene of interest
selectGene = "GNG2"
geneInterest = t(subset(TCGA.tumor, rownames(TCGA.tumor) %in% selectGene))


medianExp = median(geneInterest) # median expression of the gene of interest
boxStat = boxplot.stats(geneInterest)  # calculate stats for this gene's expressions (the extreme of the lower whisker, the lower ‘hinge’, the median, the upper ‘hinge’ and the extreme of the upper whisker, etc.)
geneInterest = as.data.frame(geneInterest)

clinical.tumor = data.frame(
  ID = data.luad.tumor$barcode,
  vital_status = data.luad.tumor$vital_status, #indicate whether the patients were dead or alive at the end of study
  days_to_death = data.luad.tumor$days_to_death, #only for patients that were dead
  days_to_last_follow_up = data.luad.tumor$days_to_last_follow_up, #only for patients that were alive at the end of study
  type = data.luad.tumor$sample_type,
  stage = data.luad.tumor$ajcc_pathologic_stage
)


## Choose one of the following methods to define high vs. low expression groups ##
# Option 1: Use the median expression value as the cutoff
geneInterest = geneInterest %>% mutate(group = case_when(geneInterest[1] > medianExp ~ paste("High ", selectGene, " expression", sep = ""),
                                                         geneInterest[1] <= medianExp ~ paste("Low ", selectGene, " expression", sep = "")))

# Option 2: Use quartiles — assign "High" to values in the upper quartile, "Low" to those in the lower quartile
geneInterest = geneInterest %>% mutate(group = case_when(geneInterest[1] >= boxStat[4] ~ paste("High ", selectGene, " expression", sep = ""),
                                                         geneInterest[1] <= boxStat[2] ~ paste("Low ", selectGene, " expression", sep = "")))

# Option 3: Use fixed sample count cutoffs (e.g., top and bottom N samples)
cutoff.high = geneInterest[order(geneInterest[1], decreasing = T)][126]
cutoff.low = geneInterest[order(geneInterest[1], decreasing = F)][376]
geneInterest = geneInterst %>% 
  mutate(group = case_when(geneInterest[1] >= cutoff.high ~ paste("High ", selectGene, " expression", sep = ""),
                           geneInterest[1] <= cutoff.low ~ paste("Low ", selectGene, " expression", sep = "")))





geneInterest$ID = rownames(geneInterest)

clinical_all.tumor = merge(clinical.tumor, geneInterest, by="ID")


clinical_all_surv.tumor = clinical_all.tumor %>% mutate(days = coalesce(days_to_death, days_to_last_follow_up)) # combine "days_to_death" and "days_to_last_follow_up" to one column named "days"

clinical_all_surv.tumor = clinical_all_surv.tumor %>% mutate(status = case_when(vital_status == "Alive" ~ 1,
                                                              vital_status == "Dead" ~ 2)) # add another column "status" to indicate that the patient is alive (1) or dead (2)


# Fit survival data using the Kaplan-Meier method
surv_patient.tumor <- Surv(time = clinical_all_surv.tumor$days, event = clinical_all_surv.tumor$status) # A '+' behind survival time indicates censored data points (i.e. alive vs dead)

fit.tumor <- survfit(surv_patient.tumor ~ clinical_all_surv.tumor$group, data = clinical_all_surv.tumor) # calculate survival fit based on groups

#summary(fit)
surv_table.tumor <- surv_summary(fit.tumor, data = clinical_all_surv.tumor) # create a summary table

surv_med_result.tumor <- surv_median(fit.tumor,combine=T) # calculate survival median

surv_plot.tumor <- ggsurvplot(fit.tumor, 
                        pval = T,                 # show p-value
                        conf.int = F,             # show confidence intervals
                        #conf.int.style = "step",  # customize style of confidence intervals
                        xlab = "Time since diagnosis (year)",   # X axis label
                        xscale = "d_y", # transform time labels from one unit to another. For example, "d_y" transform labels from days to years. 
                        break.time.by = 365.25,       # specify X axis in time intervals. (note: 1 year = 365.25 days, 2 years = 730.5 days, etc.)
                        #xlim = c(0, 1826.25), # 5 years = 365.25*5 = 1826.25 limit. Use only if needed
                        ggtheme = theme_classic(), # customize plot and risk table with a theme.
                        risk.table = "absolute",  # absolute number and percentage at risk.
                        risk.table.y.text.col = T,# colour risk table text annotations.
                        risk.table.y.text = FALSE,# show bars instead of names in text annotations in legend of risk table.
                        #ncensor.plot = TRUE,      # plot the number of censored subjects at time t
                        #pval.method = T,        # display the test name used for calculating the p-value 
                        surv.median.line = "none",  # add the median survival pointer.
                        legend.labs = c(paste("High ", selectGene, sep = ""), paste("Low ", selectGene, sep = "")),    # change legend labels.
                        palette = 
                          c("grey35", "cyan3"), # custom color palettes.
                        legend.title="",          #remove legend title (Strata)
                        pval.size = 9,           # p-value font size
                        size =3,                 # survival line width
                        #font.tickslab = 15,     # change font size of tick labels. Can also change color and format with font.tickslab = c(14, "bold.italic", "darkred")
                        #font.x = 15,
                        #font.y = 15,
                        fontsize = 5.5,      # change font size in the risk table
                        surv.scale="percent"    # change y-axis scale to percentage
)
# adjust text font size on the plot
surv_plot.tumor$plot <- surv_plot.tumor$plot + 
  theme(legend.title = element_text(size = 30, color = "black"), #can add face = "bold" to change the font format
        legend.text = element_text(size = 27, color = "black"),
        axis.text.x = element_text(size = 30, color = "black"),
        axis.text.y = element_text(size = 30, color = "black"),
        axis.title.x = element_text(size = 30, color = "black"),
        axis.title.y = element_text(size = 30, color = "black"))
surv_plot.tumor$plot$theme$line$size <- 2 # change line width in the plot


surv_plot.tumor$table$theme$axis.title.x$size <- 30 # risk table title font size
surv_plot.tumor$table$theme$line$size <-2  # risk table line width
surv_plot.tumor$table$theme$axis.text.x <- element_text(size = 30, color = "black")  # risk table x tick labels font size
surv_plot.tumor$table$theme$text <- element_text(size = 30, color = "black")  # font size of "Number at risk"


tiff(paste("TCGAbiolinks_LBM_selectGene_survival_", selectGene, ".tiff", sep = ""), width = 5000, height = 6000, units = "px", res = 600, compression = 'lzw')
surv_plot.tumor
dev.off()


# To calculate fisher test
#clinical_all_surv.tumor = subset(clinical_all_surv.tumor, !stage %in% "Stage IV") # if want to filter out certain group
clinical_all_surv.tumor = clinical_all_surv.tumor %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ "Stage I or II",
                                   stage == "Stage IA" ~ "Stage I or II",
                                   stage == "Stage IB" ~ "Stage I or II",
                                   stage == "Stage II" ~ "Stage I or II",
                                   stage == "Stage IIA" ~ "Stage I or II",
                                   stage == "Stage IIB" ~ "Stage I or II",
                                   stage == "Stage IIIA" ~ "Stage III or IV",
                                   stage == "Stage IIIB" ~ "Stage III or IV",
                                   stage == "Stage IV" ~ "Stage III or IV")) # consider stage 1 and 2 as early stage, while stage 3 and 4 as late stage



df.gene = clinical_all_surv.tumor[,c("group", "overall.Stage")]


# Fisher's exact test
padjMethod = "holm"
fisherResult.gene = pairwise_fisher_test(table(df.gene), detailed = T, p.adjust.method = padjMethod)
fisherResult.gene$padj.method = padjMethod

stageTable.gene = as.data.frame(table(df.gene))
stageTable.gene.1 = subset(stageTable.gene, stageTable.gene$group %in% paste("High ", selectGene, " expression", sep = ""))
totalCount.gene.1 = sum(stageTable.gene.1$Freq)
stageTable.gene.2 = subset(stageTable.gene, stageTable.gene$group %in% paste("Low ", selectGene, " expression", sep = ""))
totalCount.gene.2 = sum(stageTable.gene.2$Freq)

stageTable.gene$Percentage[1] = stageTable.gene$Freq[1]/totalCount.gene.1
stageTable.gene$Percentage[2] = stageTable.gene$Freq[2]/totalCount.gene.2
stageTable.gene$Percentage[3] = stageTable.gene$Freq[3]/totalCount.gene.1
stageTable.gene$Percentage[4] = stageTable.gene$Freq[4]/totalCount.gene.2

list_of_datasets <- list("FisherResult" = fisherResult.gene, "stageTable" = stageTable.gene, "Data" = clinical_all_surv.tumor)
write.xlsx(list_of_datasets, file = paste("TCGAbiolinks_LBM_selectGene_data_", selectGene, ".xlsx", sep=""))





### ============================================================================
### 12. Export TPM values by tissue type and tumor stage (optional)


### Export TPM for all samples and different stages (Run only if needed after section 4) ###
# Note: Some samples were filtered out in Sections 2 and 3. To retain specific samples, consider modifying or skipping these sections as needed.

LUAD.TPM = as.data.frame(dataFilt.TPM.luad) # convert to dataframe so that it can be exported as xlsx

# Export all samples
write.xlsx(LUAD.TPM, file = paste("TCGAbiolinks_LBM_", query.luad$project, "_TPM_all.xlsx"), rowNames = T)

# Export by tissue type and tumor stage
clinical = data.frame(
  ID = data.luad$barcode,
  vital_status = data.luad$vital_status, #indicate whether the patients were dead or alive at the end of study
  days_to_death = data.luad$days_to_death, #only for patients that were dead
  days_to_last_follow_up = data.luad$days_to_last_follow_up, #only for patients that were alive at the end of study
  type = data.luad$sample_type,
  stage = data.luad$ajcc_pathologic_stage
)
clinical = clinical %>% 
  mutate(overall.Stage = case_when(stage == "Stage I" ~ "Stage I",
                                   stage == "Stage IA" ~ "Stage I",
                                   stage == "Stage IB" ~ "Stage I",
                                   stage == "Stage II" ~ "Stage II",
                                   stage == "Stage IIA" ~ "Stage II",
                                   stage == "Stage IIB" ~ "Stage II",
                                   stage == "Stage IIIA" ~ "Stage III",
                                   stage == "Stage IIIB" ~ "Stage III",
                                   stage == "Stage IV" ~ "Stage IV"))
clinical_tumor = subset(clinical, clinical$type != "Solid Tissue Normal") #remove normal sample to get only tumor samples
clinical_normal = subset(clinical, clinical$type == "Solid Tissue Normal") #get only normal samples
clinical_stage1 = subset(clinical_tumor, clinical_tumor$overall.Stage == "Stage I") #get only samples of specific stages (only considering tumor samples)
clinical_stage2 = subset(clinical_tumor, clinical_tumor$overall.Stage == "Stage II")
clinical_stage3 = subset(clinical_tumor, clinical_tumor$overall.Stage == "Stage III")
clinical_stage4 = subset(clinical_tumor, clinical_tumor$overall.Stage == "Stage IV")

LUAD.tumor = LUAD.TPM %>% select(clinical_tumor$ID) #select columns by column name
LUAD.normal = LUAD.TPM %>% select(clinical_normal$ID)
LUAD.stage1 = LUAD.TPM %>% select(clinical_stage1$ID)
LUAD.stage2 = LUAD.TPM %>% select(clinical_stage2$ID)
LUAD.stage3 = LUAD.TPM %>% select(clinical_stage3$ID)
LUAD.stage4 = LUAD.TPM %>% select(clinical_stage4$ID)

# Export multiple dataframes into separate sheets of an excel file
list_export_TN = list("Tumor" = LUAD.tumor, "Normal" = LUAD.normal) # "Name of sheet1" = dataframe1, "Name of sheet2" = dataframe2, etc.
list_export_stage = list("Normal" = LUAD.normal, "Stage1" = LUAD.stage1, "Stage2" = LUAD.stage2, "Stage3" = LUAD.stage3, "Stage4" = LUAD.stage4)
write.xlsx(list_export_TN, file = paste("TCGAbiolink_LBM_", query.luad$project, "_TPM_Tumor_Normal.xlsx"), rowNames = T)
write.xlsx(list_export_stage, file = paste("TCGAbiolinks_LBM_", query.luad$project, "_TPM_individualStages.xlsx"), rowNames = T)

