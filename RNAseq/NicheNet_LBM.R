# ------------------------------------------------------------------------------
# Title:       NicheNet analysis for Lung-Brain Metastasis project
# Author:      Albert Wang
# Last updated:        2025-06-28
# ------------------------------------------------------------------------------


### Summary
# This script performs a comprehensive NicheNet analysis (https://doi.org/10.1038/s41592-019-0667-5) tailored for my lung-brain metastasis (LBM) project.
# NicheNet is a computational method that predicts how ligands from one cell population influence gene expression in another, helping to uncover cell–cell communication mechanisms. 
# By integrating prior knowledge of ligand–receptor interactions and intracellular signaling pathways, NicheNet identifies ligands most likely responsible for observed transcriptional changes in target (receiver) cells. 
# 
# The analysis begins by loading the NicheNet model (v2.1.5) and pre-processed bulk RNA-seq expression data from the LBM project, and specifies gene sets of interest as target genes for downstream analysis.
# The pipeline defines sender and receiver cell populations, identifies their expressed genes, and filters potential ligands using the NicheNet model alongside criteria such as expression thresholds, differential expression, and optionally provided custom ligand lists. 
# Using NicheNet’s model, the script prioritizes ligands based on their ability to predict target gene expression, infers ligand–receptor and ligand–target regulatory networks, and visualizes these relationships via heatmaps and Circos plots.
# Additional visualizations include ligand–target activity networks and signaling pathway inference, with optional export for Cytoscape integration. 

### Version iteration
# v3: Changed the criteria for selecting expressed genes (v2 mean -> v3 median) and group specific ligands (v2 used counts -> v3 use log2FC and p-values)
# v4: Add a filter for getting tumor-specific and DE ligands. These ligands are common in tumors (LT and BrM), but differentially expressed compared to LN. Also add the option to add receiver ligands to the "expressed_ligands". 
#     Add Matrisome labeling to target genes if applicable. Update NicheNet model from version 1.1.0 to version 2.1.5. Criteria for expressed genes is changed back to using mean to include more genes as suggested for bulk sequencing.
# v5: Add ability to specify the ligand to consider (by providing a list of genes) and identify ligands as up or down regulated in the circo plot
# LBM: Final version for the lung-brain metastasis (LBM) project. Fine tune code and comment.


###=============================================================================
### Required packages

library(nichenetr)
library(tidyverse)
library(openxlsx)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(ggpubr)
library(circlize)

## NicheNet package and references
# https://github.com/saeyslab/nichenetr
# https://doi.org/10.1038/s41592-019-0667-5


### ============================================================================
### NicheNet Analysis Pipeline

# Sections:
# 1. Load NicheNet Model and Expression Dataset of Interest
# 2. Define Expressed Genes in Sender and Receiver Cell Populations
# 3. Define Gene Set of Interest and Background Genes
# 4. Determine the Potential Ligands
# 5. Perform NicheNet's Ligand Activity Analysis on the Gene Set of Interest
# 6. Infer Target Genes of Top-Ranked Ligands
# 7. Ligand-Receptor Network Inference for Top-Ranked Ligands
# 8. Visualize Expression of Top-Predicted Ligands and Their Target Genes
# 9. Visualize Top Ligand–Target Interactions using Circos Plots
# 10. Visualize Top Ligand–Receptor Interactions using Circos Plots
# 11. Infer Ligand-to-Target Signaling Paths


### ============================================================================
### 1. Load NicheNet Model and Expression Dataset of Interest

## Legacy code from version 1.1.0 (skip)
#ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))


## Using NicheNet model (v2.1.5), downloaded on 2024-08-14 for faster access. Replaces previous version (v1.1.0).
# https://arxiv.org/abs/2404.16358
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/model_evaluation.md

## Load NicheNet model from online repository (Human/Mouse)
#organism <- "human"

#if(organism == "human"){
#  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
#  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
#  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
#} else if(organism == "mouse"){
#  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#  }



## Load NicheNet model from local directory (pre-downloaded files)
ligand_target_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")

ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns


## Select the expression dataset of interest
LBMfile <- file.choose() # ex. LBM_readcounts_NicheNetFormat_expression_v2.txt (Originally from LBM_readcounts_NicheNetFormat_v2.xlsx)
LBMexpression <- read.table(LBMfile, header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)
expression <- t(LBMexpression)

# Sample annotation file
LBMinfoFile <- file.choose()# ex. LBM_readcounts_NicheNetFormat_sample_info_v2.txt
sample_info <- read.table(LBMinfoFile, header=TRUE ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)

# Differential expression results used to label Circos plot
DEresultFile <- file.choose() # ex. LBM_readcounts_NicheNetFormat_DEresult.txt
DEresult <- read.table(DEresultFile, header=TRUE ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)



### ============================================================================
### 2. Define Expressed Genes in Sender and Receiver Cell Populations

# Select group for comparison
compSelect = "BrMvsLT" # <--------------------------------------------------------------------------------------------------------------- option: "BrMvsLT" or "LTvsLN"

# Consider including receiver ligands in expressed_ligands? (This includes ligands that are highly expressed in receiver but not in sender)
AddRecLigand = "Yes" # <------------------------------------------------------------------------------------------------------------- option: "Yes" or "No"

# Keep only tumor-specific or DE ligands?
# For BrMvsLT, this option retains tumor-specific ligands (i.e. remove ligands not DE in both LTvsLN and BrMvsLT). For LTvsLN, this option keeps the DE ligands in LTvsLN.
tumorLigandOnly = "Yes" # <-------------------------------------------------------------------------------------------------------------option: "Yes" or "No"

# Use a specified list of ligands? 
# Provide a custom list of ligand (genes). Only expressed ligands from this list will be included in the analysis.
specLigand = "Yes" # <-------------------------------------------------------------------------------------------------------------------option: "Yes" or "No"

# Label the ligands in the circo plot as up- or down-regulated? 
# By default, ligands are labeled as sender- or receiver-specific. This option changes labeling to indicate up- or down-regulation.
ligandLabel = "DE" # <-------------------------------------------------------------------------------------------------------------option: "DE" or "default"



# Add file name suffixes based on selected options
if (AddRecLigand == "Yes" & tumorLigandOnly != "Yes"){
  AddFileName = "_AddRecLig"
} else if (AddRecLigand != "Yes" & tumorLigandOnly == "Yes"){
  AddFileName = "_tumorLigOnly"
} else if (AddRecLigand == "Yes" & tumorLigandOnly == "Yes"){
  AddFileName = "_AddRecLig-tumorLigOnly"
} else { # if none of the options were used
  AddFileName = ""
}

# Extract sample IDs for sender and receiver from annotation data
if (compSelect == "BrMvsLT"){
  sender_ids = sample_info %>% filter(`Tissue` == "Lung tumor") %>% pull(Sample) # possible options: "Lung normal", "Lung tumor", or "Brain Met"
  receiver_ids = sample_info %>% filter(`Tissue` == "Brain Met") %>% pull(Sample)
} else if (compSelect == "LTvsLN"){
  sender_ids = sample_info %>% filter(`Tissue` == "Lung normal") %>% pull(Sample)
  receiver_ids = sample_info %>% filter(`Tissue` == "Lung tumor") %>% pull(Sample)
} else {stop("No comparison selected or incorrect comparison")}



### Determine expressed genes in sender and receiver

## Legacy code (skip)
#expressed_genes_sender = expression[sender_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 50] %>% names()
#expressed_genes_receiver = expression[receiver_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 50] %>% names()


## If using DEGs as the expressed gene set, load the DEG file first
#DEGfile <- file.choose()
#DEG <- read.table(DEGfile, header=TRUE, sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE) # DEG is abs(log2FC) > 0.6 and padj < 0.05

#DEG_ids = DEG %>% filter(`BTvsLT` == "DEG") %>% pull(Gene) #'LTvsLN', 'BTvsLT', or 'BTvsLN'
#expressed_genes_sender = expression[sender_ids,DEG_ids] %>% colnames() 
#expressed_genes_receiver = expression[receiver_ids,DEG_ids] %>% colnames()


## Use mean count values to identify expressed genes
expressed_genes_sender = expression[sender_ids,] %>% apply(2,function(x){mean(x)}) %>% .[. >= 5] %>% names() # mean counts greater than 5
expressed_genes_receiver = expression[receiver_ids,] %>% apply(2,function(x){mean(x)}) %>% .[. >= 5] %>% names() 


# Check the number of expressed genes: should be a 'reasonable' number of total expressed genes in a cell type, e.g. between 5000-10000 (and not 500 or 20000) for single cell, and 10000-15000 for bulk
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/faq.md
length(expressed_genes_sender)
## [1] 6706
length(expressed_genes_receiver)
## [1] 6351


### ============================================================================
### 3. Define Gene Set of Interest and Background Genes

## Select gene set of interest (target)
genesetFile <- file.choose() # ex. iPathway_GOBP_BrMvsLT_ECM organization.txt
geneset_oi = readr::read_tsv(genesetFile, col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] %>% .[. %in% colnames(expression)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(geneset_oi)
## [1] "SERPINE1" "TGFBI"    "MMP10"    "LAMC2"    "P4HA2"    "PDPN"

fileName = basename(genesetFile) #get file name
fileName = tools::file_path_sans_ext(fileName) #remove extension
compName <- str_extract(fileName, "BrMvsLT|LTvsLN") # determine if the comparison is BrMvsLT or LTvsLN based on gene set of interest

if (compSelect != compName){
  stop("Grouping of the gene set of interest differs from the expression data grouping")
}

# If the gene set is from the ECM matrisome (http://matrisomeproject.mit.edu/), provide a gene label file specifying categories.
# For example, ECM matrisome_compiled_BTvsLT_geneLabel.txt
if (grepl("matrisome", fileName)){
  geneCatFile = file.choose()
  geneCat = read.table(geneCatFile, header=TRUE ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)
  colnames(geneCat) = c("target", "Matrisome category")
  
  labelFileName = basename(geneCatFile)
  if (str_extract(labelFileName, "BrMvsLT|LTvsLN") != compName){
    stop("Grouping of the gene set of interest differs from the expression data grouping")
  } else if (!grepl("geneLabel", basename(geneCatFile))){
    stop("This file is not a gene label file")
  }
}

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)
## [1] "RPS11"   "ELMO2"   "PNMA1"   "MMP2"    "TMEM216" "ERCC5"



### ============================================================================
### 4. Determine the Potential Ligands

## Legacy code from version 1.1.0 (skip)
#lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#lr_network = readRDS("lr_network.rds")

## version 2.1.5
lr_network = readRDS("lr_network_human_21122021.rds")


# Optional: users can exclude ligand-receptor interactions predicted based on protein-protein interactions, retaining only those described in curated databases.
# To do this: uncomment following line of code:
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

# Add the receiver ligands to the expressed_ligands if this option is selected. These ligands are considered as autocrine signaling
if (AddRecLigand == "Yes"){
  expressed_ligands_receiver = intersect(ligands, expressed_genes_receiver)
  expressed_ligands = union(expressed_ligands_receiver, expressed_ligands)
}


# Only consider specific ligands if a list is provided
if (specLigand == "Yes"){
  specLigandFile = file.choose()
  specLigandList = read.table(specLigandFile, header=F ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)
  
  expressed_ligands = intersect(specLigandList$V1, expressed_ligands)
  
}



# Filter out non-DE ligands to yield tumor-specific ligands only
DEresult_notDE = DEresult %>% select(c("gene", "padj_LTvsLN", "log2FC_LTvsLN", "padj_BrMvsLT", "log2FC_BrMvsLT"))
if (tumorLigandOnly == "Yes" & compSelect == "BrMvsLT"){
  
  # find ligands that are not differentially expressed across all LN, LT, and BrM (i.e. not DE in both LTvsLN and BrMvsLT)
  DEresult_notDE = DEresult_notDE %>% filter((padj_LTvsLN>0.05 | abs(log2FC_LTvsLN)<0.6) & (padj_BrMvsLT>0.05 | abs(log2FC_BrMvsLT)<0.6))
  
  ligand_nonDE = intersect(DEresult_notDE$gene, expressed_ligands)
  ligand_tumorOnly = setdiff(expressed_ligands, DEresult_notDE$gene)
  
  expressed_ligands = ligand_tumorOnly
  
} else if (tumorLigandOnly == "Yes" & compSelect == "LTvsLN"){
  # find ligands that are not differentially expressed across LN and LT (i.e. not DE in LTvsLN)
  DEresult_notDE = DEresult_notDE %>% filter((padj_LTvsLN>0.05 | abs(log2FC_LTvsLN)<0.6))
  
  ligand_nonDE = intersect(DEresult_notDE$gene, expressed_ligands)
  ligand_tumorOnly = setdiff(expressed_ligands, DEresult_notDE$gene)
  
  expressed_ligands = ligand_tumorOnly
}

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)
## # A tibble: 6 x 4
##   from    to        source         database
##   <chr>   <chr>     <chr>          <chr>   
## 1 HGF     MET       kegg_cytokines kegg    
## 2 TNFSF10 TNFRSF10A kegg_cytokines kegg    
## 3 TNFSF10 TNFRSF10B kegg_cytokines kegg    
## 4 TGFB2   TGFBR1    kegg_cytokines kegg    
## 5 TGFB3   TGFBR1    kegg_cytokines kegg    
## 6 INHBA   ACVR2A    kegg_cytokines kegg

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)
## [1] "HGF"     "TNFSF10" "TGFB2"   "TGFB3"   "INHBA"   "CD99"



### ============================================================================
### 5. Perform NicheNet's Ligand Activity Analysis on the Gene Set of Interest

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-pearson) 
## # A tibble: 131 x 4
##    test_ligand auroc   aupr pearson
##    <chr>       <dbl>  <dbl>   <dbl>
##  1 PTHLH       0.667 0.0720   0.128
##  2 CXCL12      0.680 0.0507   0.123
##  3 AGT         0.676 0.0581   0.120
##  4 TGFB3       0.689 0.0454   0.117
##  5 IL6         0.693 0.0510   0.115
##  6 INHBA       0.695 0.0502   0.113
##  7 ADAM17      0.672 0.0526   0.113
##  8 TNC         0.700 0.0444   0.109
##  9 CTGF        0.680 0.0473   0.108
## 10 FN1         0.679 0.0505   0.108
## # ... with 121 more rows
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
head(best_upstream_ligands)
## [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") + #pearson correlation coefficient (PCC)
  theme_classic()
p_hist_lig_activity




### ============================================================================
### 6. Infer Target Genes of Top-Ranked Ligands

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)
## # A tibble: 6 x 3
##   ligand target  weight
##   <chr>  <chr>    <dbl>
## 1 PTHLH  COL1A1 0.00399
## 2 PTHLH  MMP1   0.00425
## 3 PTHLH  MMP2   0.00210
## 4 PTHLH  MYH9   0.00116
## 5 PTHLH  P4HA2  0.00190
## 6 PTHLH  PLAU   0.00401

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)

nrow(active_ligand_target_links_df)
## [1] 143
head(active_ligand_target_links_df)
## # A tibble: 6 x 3
##   ligand target  weight
##   <chr>  <chr>    <dbl>
## 1 PTHLH  COL1A1 0.00399
## 2 PTHLH  MMP1   0.00425
## 3 PTHLH  MMP2   0.00210
## 4 PTHLH  MYH9   0.00116
## 5 PTHLH  P4HA2  0.00190
## 6 PTHLH  PLAU   0.00401

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets = intersect(order_targets, row.names(active_ligand_target_links))
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized sender-ligands","Genes in cancer cells", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

#tiff("LBM_NicheNet_ligand_targetGenes_BTvsLT.tiff", width = 5000, height = 3000, units = "px", res = 600, compression = 'lzw')
p_ligand_target_network
#dev.off()


### ============================================================================
### 7. Ligand-Receptor Network Inference for Top-Ranked Ligands

# Get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()


## Extract the interaction weights for ligand-receptor pairs as defined in the NicheNet model

## Legacy code from version 1.1.0 (skip)
#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
#weighted_networks = readRDS("weighted_networks.rds")

## Version 2.1.5
weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")

lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# Convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# Perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized sender-ligands","Receptors expressed by cancer cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")


#tiff("LBM_NicheNet_ligand_receptor_interactionPotential_BTvsLT.tiff", width = 12000, height = 3000, units = "px", res = 600, compression = 'lzw')
p_ligand_receptor_network

#dev.off()


### ============================================================================
### 8. Visualize Expression of Top-Predicted Ligands and Their Target Genes


ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized sender-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\n(target gene prediction ability)")
p_ligand_pearson

# Expression of ligands in sender cells per tumor
expression_df_sender = expression[sender_ids,order_ligands] %>% data.frame() %>% rownames_to_column("Sample") %>% as_tibble() %>% inner_join(sample_info %>% select(Sample), by =  "Sample")

aggregated_expression_sender = expression_df_sender %>% select(-Sample) %>% summarise_all(median)

aggregated_expression_df_sender = aggregated_expression_sender %>% t() %>% magrittr::set_colnames("Median") %>% data.frame() %>% rownames_to_column("ligand") %>% as_tibble() 


aggregated_expression_matrix_sender = aggregated_expression_df_sender %>% select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_sender$ligand)


vis_ligand_tumor_expression = aggregated_expression_matrix_sender

color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
p_ligand_tumor_expression = vis_ligand_tumor_expression %>% make_heatmap_ggplot("Prioritized sender-ligands","Sender", color = color[100],legend_position = "top", x_axis_position = "top", legend_title = "Ligand expression") + theme(axis.text.y = element_text(face = "italic"))
#tiff("LBM_NicheNet_Median_expression_ligand_BTvsLT.tiff", width = 3000, height = 3000, units = "px", res = 600, compression = 'lzw')
p_ligand_tumor_expression
#dev.off()

# Expression of target genes in malignant cells per tumor
expression_df_target = expression[receiver_ids,geneset_oi] %>% data.frame() %>% rownames_to_column("Sample") %>% as_tibble() %>% inner_join(sample_info %>% select(Sample), by =  "Sample") 

aggregated_expression_target = expression_df_target %>% select(-Sample) %>% summarise_all(median)

aggregated_expression_df_target = aggregated_expression_target %>% t() %>% magrittr::set_colnames("Median") %>% data.frame() %>% rownames_to_column("target") %>% as_tibble() 

aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)

vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% .[,order_targets] %>% as.matrix() %>% magrittr::set_colnames("Median expression") %>% t()

p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Tumor","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Tumor expression") + theme(axis.text.x = element_text(face = "italic"))
#tiff("LBM_NicheNet_Median_expression_target_BTvsLT.tiff", width = 3500, height = 3000, units = "px", res = 600, compression = 'lzw')
p_target_tumor_scaled_expression
#dev.off()

# Combine different heatmaps in one overview figure
tiff(paste("LBM_NicheNet_overview_", fileName, AddFileName, ".tiff", sep = ""), width = 12000, height = 5000, units = "px", res = 600, compression = 'lzw')
figures_without_legend = plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  p_ligand_tumor_expression + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""), 
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none", axis.ticks = element_blank()) + xlab(""), 
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_pearson)+ 5, ncol(vis_ligand_tumor_expression)+5, ncol(vis_ligand_target)) -2,
  rel_heights = c(nrow(vis_ligand_pearson), nrow(vis_target_tumor_expression_scaled) + 4)) 

legends = plot_grid(
  as_ggplot(get_legend(p_ligand_pearson)),
  as_ggplot(get_legend(p_ligand_tumor_expression)),
  as_ggplot(get_legend(p_ligand_target_network)),
  as_ggplot(get_legend(p_target_tumor_scaled_expression)),
  nrow = 2,
  align = "h")

plot_grid(figures_without_legend, 
          legends, 
          rel_heights = c(10,2), nrow = 2, align = "hv")
dev.off()



### ============================================================================
### 9. Visualize Top Ligand–Target Interactions using Circos Plots


## Cutoff threshold for circo plot: Remove 66% (default from vignette, 75% seems to suit better for most LBM data) of links with lowest score
cutoff_circos = 0.75 # <---------------------------------------------------------------------------------------------------------------------- cutoff_include_all_ligands

expressed_ligands_sender = intersect(ligands, expressed_genes_sender)
expressed_ligands_receiver = intersect(ligands, expressed_genes_receiver)

best_upstream_ligands %>% intersect(expressed_ligands_sender) 
##  [1] "PTHLH"  "CXCL12" "AGT"    "TGFB3"  "IL6"    "INHBA"  "ADAM17" "TNC"    "CTGF"   "FN1"    "BMP5"   "IL24"  
## [13] "CXCL11" "MMP9"   "COL4A1" "PSEN1"  "CXCL9"
best_upstream_ligands %>% intersect(expressed_ligands_receiver)
##  [1] "EDN1"   "CXCL12" "IL6"    "ADAM17" "VWF"    "CTGF"   "FN1"    "SPP1"   "CXCL11" "COL4A1" "PSEN1"  "CXCL9"

### Legacy code – not applicable to LBM project (skip)
## Many ligands are potentially expressed in both cell types. Identify which ligands show stronger expression in each cell type.

#ligand_expression_tbl = tibble(
#  ligand = best_upstream_ligands, 
#  sender = expression[sender_ids,best_upstream_ligands] %>% 
#    #apply(2,function(x){10*(2**x - 1)}) %>% 
#    apply(2,function(x){log2(mean(x) + 1)}), # 10*(2**x - 1) transform linear scale (2^x). Next, calculate mean, then apply log2 transformation
#  receiver = expression[receiver_ids,best_upstream_ligands] %>% 
#    #apply(2,function(x){10*(2**x - 1)}) %>% 
#    apply(2,function(x){log2(mean(x) + 1)})) # can treat as autocrine?

#sender_specific_ligands = ligand_expression_tbl %>% filter(sender > receiver + 2) %>% pull(ligand) #filter expression level of a ligand in the sender cells is greater than in receiver cells by more than 2 units
#receiver_specific_ligands = ligand_expression_tbl %>% filter(receiver > sender + 2) %>% pull(ligand)
#common_ligands = setdiff(best_upstream_ligands,c(sender_specific_ligands,receiver_specific_ligands))


## Use log2FC and adjusted p-value (padj) to identify ligands more strongly expressed in either the sender or receiver.
# Ligands that are not differentially expressed are considered common ligands to both cell types.
DEresult_comp = DEresult %>% select(c("gene", paste("padj_", compName, sep = ""), paste("log2FC_", compName, sep = ""))) # extract the padj and log2FC for the comparison
DEresult_comp = DEresult_comp %>% filter(gene %in% best_upstream_ligands)
colnames(DEresult_comp) = c("ligand", "padj", "log2FC")

sender_specific_ligands = DEresult_comp %>% filter(padj < 0.05) %>% filter(log2FC < 0) %>% pull(ligand) # In BrMvsLT, downregulated genes suggest LT (sender) has higher expression than BrM (receiver). Adjust if needed.
receiver_specific_ligands = DEresult_comp %>% filter(padj < 0.05) %>% filter(log2FC > 0) %>% pull(ligand)
common_ligands = setdiff(best_upstream_ligands,c(sender_specific_ligands,receiver_specific_ligands))



ligand_type_indication_df = tibble(
  ligand_type = c(rep("Sender-specific", times = sender_specific_ligands %>% length()),
                  rep("Common", times = common_ligands %>% length()),
                  rep("Receiver-specific", times = receiver_specific_ligands %>% length())),
  ligand = c(sender_specific_ligands, common_ligands, receiver_specific_ligands))

## Get the active ligand-target links
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

if (grepl("matrisome", fileName)){ # If the gene set of interest is from the ECM matrisome, a gene category file should be provided (stored in geneCat)
  active_ligand_target_links_df = active_ligand_target_links_df %>% inner_join(geneCat) %>% # classify target genes based on Matrisome category
    rename(target_type = `Matrisome category`) %>% # rename column name from "Matrisome category" to "target_type
    inner_join(ligand_type_indication_df) 
} else { # default if no gene label is supplied. All target genes will be labeled as "target" in target_type
  active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "target") %>% inner_join(ligand_type_indication_df) # if you want to make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
  
}



## remove 66% (default from vignette; adjustable via cutoff_circos) of links with lowest scores
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(cutoff_circos)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
circos_target = circos_links

## Prepare circos plot visualization (specify colors for ligands and targets)
# Ligand colors
if (ligandLabel == "DE"){
  grid_col_ligand = c("Sender-specific" = "royalblue", # If it's sender-specific, then it's a downregulated gene (ex. BrMvsLT, sender is LT and reciever is BrM). Sender-specific means LT has a higher expression than BrM
                      "Receiver-specific" = "red3", # If it's receiver-specific, then it's a up-regulated gene
                      "Common" = "grey40" # If it's common, then it's a non-DE gene in this comparison
                      ) 
} else { # default
  if (compName == "BrMvsLT"){
    # For BrMvsLT
    grid_col_ligand =c("Common" = "royalblue",
                       "Sender-specific" = "grey40",
                       "Receiver-specific" = "forestgreen")
  } else if (compName == "LTvsLN"){
    # For LTvsLN
    grid_col_ligand =c("Common" = "royalblue",
                       "Sender-specific" = "burlywood3",
                       "Receiver-specific" = "grey40")
  } else {stop("Color not set up for this comparison")}
}



# Target colors
if (grepl("matrisome", fileName)){# If the gene set of interest is matrisome, assign colors for matrisome categories
  grid_col_target = c(
    "Collagen" = "black",
    "Secreted factors" = "orange",
    "Proteoglycan" = "skyblue",
    "Glycoprotein" = "purple3",
    "ECM regulator" = "tomato",
    "ECM affiliated" = "lawngreen"
  )
  
} else { 
  if (compName == "BrMvsLT" & ligandLabel == "DE"){
    # For BrMvsLT
    grid_col_target =c(
      "target" = "forestgreen")
  } else if (compName == "LTvsLN" & ligandLabel == "DE"){
    # For LTvsLN
    grid_col_target =c(
      "target" = "grey40")
  } else {# default
    grid_col_target =c(
      "target" = "tomato")}
}



grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# Provide option to set link transparency in the Circos plot based on ligand–target potential scores
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Order ligands and targets
if (grepl("matrisome", fileName)){ # If using Matrisome as gene set, sort target order based on matrisome category
  circos_links_reorder= circos_links %>% arrange(desc(target_type))
  target_order = circos_links_reorder$target %>% unique()
} else { # default, no sorting
  target_order = circos_links$target %>% unique()
}

ligand_order = c(sender_specific_ligands,common_ligands,receiver_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

# define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5


SenderNum = circos_links %>% filter(ligand_type == "Sender-specific") %>% distinct(ligand) %>% nrow()
CommonNum = circos_links %>% filter(ligand_type == "Common") %>% distinct(ligand) %>% nrow()
ReceiverNum = circos_links %>% filter(ligand_type == "Receiver-specific") %>% distinct(ligand) %>% nrow()
if (grepl("matrisome", fileName)){
  TargetNum = circos_links %>%  distinct(target) %>% nrow()
} else { # default
  TargetNum = circos_links %>% filter(target_type == "target") %>% distinct(target) %>% nrow()
}



# Initialize the gaps vector
gaps <- c()

# Add gaps for Sender-specific ligands if they exist
if (SenderNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (SenderNum - 1)))
  gaps <- c(gaps, width_different_cell)
}

# Add gaps for Common ligands if they exist
if (CommonNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (CommonNum - 1)))
  gaps <- c(gaps, width_different_cell)
}

# Add gaps for Receiver-specific ligands if they exist
if (ReceiverNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (ReceiverNum - 1)))
  gaps = c(gaps, width_ligand_target)
}

# Add gaps for receptors
if (TargetNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_target_type, times = (TargetNum - 1)))
  gaps <- c(gaps, width_ligand_target)
}





## Render the circos plot (degree of transparency determined by the regulatory potential value of a ligand-target interaction)
tiff(paste("LBM_NicheNet_circo_", fileName, "_ligand-target", AddFileName, ".tiff", sep = ""), width = 5500, height = 5500, units = "px", res = 600, compression = 'lzw')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# Go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()
dev.off()

## Optional: Alternate style
tiff(paste("LBM_NicheNet_circo_", fileName, "_ligand-target", AddFileName, "_alt.tiff", sep = ""), width = 5500, height = 5500, units = "px", res = 600, compression = 'lzw')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "triangle", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()
dev.off()



### ============================================================================
### 10. Visualize Top Ligand–Receptor Interactions using Circos Plots

# Visualize ligand-receptor interactions of the prioritized ligands in a circos plot

## Get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()


## Retrieve the weights of ligand–receptor interactions as defined in the NicheNet model
## Legacy code from version 1.1.0 (skip)
#weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

weighted_networks = readRDS("weighted_networks.rds")
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "target_receptor") %>% inner_join(ligand_type_indication_df)

# Ligand colors
if (ligandLabel == "DE"){
  grid_col_ligand = c("Sender-specific" = "royalblue", # If it's sender-specific, then it's a downregulated gene (ex. BrMvsLT, sender is LT and reciever is BrM). Sender-specific means LT has a higher expression than BrM
                      "Receiver-specific" = "red3", # If it's receiver-specific, then it's a up-regulated gene
                      "Common" = "grey40" # If it's common, then it's a non-DE gene in this comparison
                      )
} else { # default
  if (compName == "BrMvsLT"){
    # For BrMvsLT
    grid_col_ligand =c("Common" = "royalblue",
                       "Sender-specific" = "grey40",
                       "Receiver-specific" = "forestgreen")
  } else if (compName == "LTvsLN"){
    # For LTvsLN
    grid_col_ligand =c("Common" = "royalblue",
                       "Sender-specific" = "burlywood3",
                       "Receiver-specific" = "grey40")
  } else {stop("Color not set up for this comparison")}
}



# Receptor colors
if (compName == "BrMvsLT" & ligandLabel == "DE"){
  # For BrMvsLT
  grid_col_receptor =c(
    "target_receptor" = "darkgreen")
} else if (compName == "LTvsLN" & ligandLabel == "DE"){
  # For LTvsLN
  grid_col_receptor =c(
    "target_receptor" = "grey20")
} else {# default
  grid_col_receptor =c(
    "target_receptor" = "darkred")}



grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor) %>% filter(weight > 0.25) #filter to reduce the number of receptors shown in figure (default weight > 0.25)
links_circle = circos_links %>% select(ligand,receptor, weight)

circos_receptor = circos_links

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# Provide option to set link transparency in the Circos plot based on ligand–target potential scores
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Order ligands and receptors
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(sender_specific_ligands,common_ligands,receiver_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

# Define gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

SenderNum = circos_links %>% filter(ligand_type == "Sender-specific") %>% distinct(ligand) %>% nrow()
CommonNum = circos_links %>% filter(ligand_type == "Common") %>% distinct(ligand) %>% nrow()
ReceiverNum = circos_links %>% filter(ligand_type == "Receiver-specific") %>% distinct(ligand) %>% nrow()
ReceptorNum = circos_links %>% filter(receptor_type == "target_receptor") %>% distinct(receptor) %>% nrow()


## Initialize the gaps vector
gaps <- c()

# Add gaps for Sender-specific ligands if they exist
if (SenderNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (SenderNum - 1)))
  gaps <- c(gaps, width_different_cell) # gap for the next section
}

# Add gaps for Common ligands if they exist
if (CommonNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (CommonNum - 1)))
  gaps <- c(gaps, width_different_cell)
}

# Add gaps for Receiver-specific ligands if they exist
if (ReceiverNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_ligand_type, times = (ReceiverNum - 1)))
  gaps <- c(gaps, width_ligand_receptor)
}

# Add gaps for receptors
if (ReceptorNum > 0) {
  gaps <- c(gaps, rep(width_same_cell_same_receptor_type, times = (ReceptorNum - 1)))
  gaps <- c(gaps, width_ligand_receptor)
}





## Render circos plot (degree of transparency determined by the prior interaction weight of the ligand-receptor interaction)
tiff(paste("LBM_NicheNet_circo_", fileName, "_ligand-receptor", AddFileName, ".tiff", sep = ""), width = 5500, height = 5500, units = "px", res = 600, compression = 'lzw')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA)
circos.clear()
dev.off()

## Optional: Alternate style
tiff(paste("LBM_NicheNet_circo_", fileName, "_ligand-receptor", AddFileName, "_alt.tiff", sep = ""), width = 5500, height = 5500, units = "px", res = 600, compression = 'lzw')
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "triangle", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA)
circos.clear()
dev.off()



## Export data
list_export = list("Geneset" = geneset_oi, "Background genes" = background_expressed_genes, "Sender ligand Expression" = expression_df_sender, "Receiver target expression" = expression_df_target,
                   "Ligand activity" = ligand_activities, "Circo_target" = circos_target, "Circo_receptor" = circos_receptor, "Ligand-Receptor top" = lr_network_top_df, "Ligand-Receptor expressed" = lr_network_expressed)
write.xlsx(list_export, file = paste("LBM_NicheNet_data_", fileName, AddFileName, ".xlsx", sep = ""), rowNames = F)




### ============================================================================
### 11. Infer Ligand-to-Target Signaling Paths

# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_target_signaling_path.md

library(DiagrammeR)
library(DiagrammeRsvg)
library(magick)
library(cowplot)


weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_tf_matrix_nsga2r_final.rds"))

lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
sig_network <- readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))
gr_network <- readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))


ligands_oi <- "COL1A1" # this can be a list of multiple ligands if required
targets_oi <- c("COL5A1","COL10A1", "MMP9", "GDF15", "MMP2", "CCL5")
receptors_oi = c("ITGB1", "MAG", "CD44", "DDR1", "ITGAV")

# Get network with the specified ligand and target
active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 4,
                                                      minmax_scaling = TRUE) 

# Optional: Also specify receptor if needed (Run this or the previous one)
active_signaling_network = get_ligand_signaling_path_with_receptor(ligand_tf_matrix = ligand_tf_matrix,
                                                                   weighted_networks = weighted_networks,
                                                                   ligands_all = ligands_oi,
                                                                   targets_all = targets_oi,
                                                                   receptors_all = receptors_oi,
                                                                   top_n_regulators = 4)





graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

# To render the graph in RStudio Viewer, uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")

# To export/draw the svg, you need to install DiagrammeRsvg
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))

# Write the SVG content to a temporary file
tmp_svg <- tempfile(fileext = ".svg")
writeLines(graph_svg, con = tmp_svg)

# Read the SVG file as an image with magick
graph_image <- magick::image_read(tmp_svg)

# Set the resolution to 600 DPI and save as a png file
#magick::image_write(graph_image, path = "graph_image.png", format = "png", density = 600)
magick::image_write(graph_image, path = "NichNet_network_graph.svg")

# Display it using cowplot
cowplot::ggdraw() + cowplot::draw_image(graph_image)
#cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))



## Data sources support the interactions in this network
data_source_network <- infer_supporting_datasources(signaling_graph_list = active_signaling_network,
                                                    lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network)


## Export to Cytoscape
output_path <- ""
write_output <- T # change to TRUE for writing output

# weighted networks ('import network' in Cytoscape)
if(write_output){
  bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"),
            active_signaling_network$gr %>% mutate(layer = "regulatory")) %>%
    write_tsv(paste0(output_path,"NicheNet_Cytoscape_weighted_signaling_network.txt")) 
}

# networks with information of supporting data sources ('import network' in Cytoscape)
if(write_output){
  data_source_network %>% write_tsv(paste0(output_path,"NicheNet_Cytoscape_data_source_network.txt"))
}

# Node annotation table ('import table' in Cytoscape)
specific_annotation_tbl <- bind_rows(
  tibble(gene = ligands_oi, annotation = "ligand"),
  tibble(gene = targets_oi, annotation = "target"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_oi,ligands_oi)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
  tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_oi,ligands_oi)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
)
non_specific_annotation_tbl <- tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")

if(write_output){
  bind_rows(specific_annotation_tbl, non_specific_annotation_tbl) %>%
    write_tsv(paste0(output_path,"NicheNet_Cytoscape_annotation_table.txt"))
}