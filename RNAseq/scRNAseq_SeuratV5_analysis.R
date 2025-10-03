# ------------------------------------------------------------------------------
# Title:       Single-Cell RNA-seq Analysis using Seurat (v5)
# Author:      Albert Wang
# Last updated:        2025-10-02
# ------------------------------------------------------------------------------


### Summary
# This script performs single-cell RNA-seq analysis using the Seurat package (v5).
# It adapts content from the Harvard Chan Bioinformatics Core tutorial 
# (https://github.com/hbctraining/Intro-to-scRNAseq), 
# the NIH Center for Cancer Research (CCR) scRNA-seq seminar series 
# (https://bioinformatics.ccr.cancer.gov/docs/getting-started-with-scrna-seq/),
# and the official Seurat vignettes (https://satijalab.org/seurat/), 
# with additional modifications and enhancements by the script author.

###=============================================================================
### Installation of Seurat package
# https://satijalab.org/seurat/articles/install_v5

## Seurat package (CRAN)
install.packages('Seurat')
# From Github
#remotes::install_github("satijalab/seurat", "seurat5")
# Developer version
#remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')


## Optional packages to enhance performance
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 
                                       'https://bnprks.r-universe.dev/'))
# BPCells: enables memory-efficient storage and manipulation of large scRNA-seq.
#          It utilizes bit-packing compression to store counts matrices on disk 
#          and C++ code to cache operations. 

install.packages(c("BPCells", 
                   "presto", 
                   "glmGamPoi"))


## Optional packages to enhance the functionality of Seurat
# Signac: analysis of single-cell chromatin data
# SeuratData: automatically load datasets pre-packaged as Seurat objects
# Azimuth: local annotation of scRNA-seq and scATAC-seq queries across multiple organs and tissues
# SeuratWrappers: enables use of additional integration and differential expression methods

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages('Signac')
remotes::install_github("satijalab/seurat-data")
remotes::install_github("satijalab/azimuth")
remotes::install_github("satijalab/seurat-wrappers")

###=============================================================================
### Required packages

library(Seurat)
library(dplyr)
library(tidyr)
library(patchwork) # for plotting
library(ggplot2)
library(sctransform)
library(presto) # for more efficient differential expression
library(glmGamPoi) # for more efficient sctransform
library(BPCells)
library(tools)
library(SeuratWrappers) # collection of community-provided methods for Seurat



### Optional packages ###

### Enable parallel execution for the `PrepSCTIntegration` functions
library(future) # manages parallel plans
library(future.apply) # provides parallel versions of base apply functions

# set global parallel plan: 2 workers = process 2 Seurat objects at once
#plan(multisession, workers = 2)

# NOTE:
# - Each worker gets its own copy of the data it needs, so parallel execution
#   can use significantly more RAM.
# - After you’re done with heavy parallel steps, you can reset to serial mode with:
#       plan(sequential)
#   This releases the workers and frees their memory.
# - Manage plan locally (e.g., define in `apply_seurat_list`) to have better 
#   control of memory usage


### scVIIntegration
library(reticulate)

# Needs to install Miniconda (conda+Python) or Anaconda (conda+python+packages)
# Install scvi-tools, R, and reticulate in the conda environment
# https://docs.scvi-tools.org/en/stable/installation.html
# https://github.com/satijalab/seurat/issues/7164

## Installation
# 1) Ensure a clean Miniconda
install_miniconda(force = TRUE, update = TRUE)

# 2) Tell reticulate exactly where conda is
miniconda <- miniconda_path()
options(reticulate.conda_binary = file.path(miniconda, "condabin", "conda.bat"))
Sys.setenv(RETICULATE_MINICONDA_PATH = miniconda)

# 3) Confirm R now sees conda
reticulate::conda_binary()   # should print .../r-miniconda/condabin/conda.bat

# 4) Create env + install
conda_create("r-scvi", python_version = "3.11.13")
use_condaenv("r-scvi", required = TRUE)

py_install("pip", envname = "r-scvi")
py_install("torch", envname = "r-scvi", pip = TRUE)     # CPU torch
py_install("scvi-tools", envname = "r-scvi", pip = TRUE)

py_module_available("scvi")  # should be TRUE


## If installation is complete, start here
# 1) Use the env you created
use_condaenv("r-scvi", required = TRUE)

# 2) Import scvi and set workers BEFORE training/integration
scvi <- import("scvi")

workers <- max(1L, min(8L, parallel::detectCores(logical = TRUE) - 1L))
scvi$settings$dl_num_workers <- as.integer(workers)
scvi$settings$dl_persistent_workers <- TRUE   # set FALSE if Windows gives issues

# (Optional) sanity check
scvi$settings$dl_num_workers
scvi$settings$dl_persistent_workers



### Tutorial dataset
library(SeuratData)
# For tutorial datasets, define the apply_seurat_list() function first
# before proceeding directly to the normalization section.

## ifnb
InstallData("ifnb")
ifnb <- LoadData("ifnb")
ifnb[["RNA"]] <- split(ifnb[["RNA"]], f = ifnb$stim)
ifnb
filtered_seurat <- ifnb

## pbmcsca
# example dataset for integration 
# (https://satijalab.org/seurat/articles/seurat5_integration)
#
# Install Azimuth:
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("TFBSTools")
# devtools::install_github("satijalab/azimuth", "seurat5")
library(Azimuth)
InstallData("pbmcsca")
pbmcsca <- LoadData("pbmcsca")
pbmcsca <- subset(pbmcsca, nFeature_RNA > 1000)
pbmcsca <- RunAzimuth(pbmcsca, reference = "pbmcref")
# This data contains 9 different batches (stored in Method)
# To integrate the batches, split them into different layers in the Seurat object
# (previous versions require to be split into 9 different Seurat objects)
# After splitting, there are 18 layers (a counts and data layer for each batch)
# Note that normalization and variable feature identifications is performed for
# each layer/batch independently.
# After integrative analysis, the layers can be rejoin (`JoinLayers') and recreates
# counts and data layers. This is required to perform differential expression.
pbmcsca[["RNA"]] <- split(pbmcsca[["RNA"]], f = pbmcsca$Method)
pbmcsca
filtered_seurat <- pbmcsca


## Seurat package references
# https://satijalab.org/seurat/
# https://github.com/satijalab/seurat


## Seurat package index
# https://satijalab.org/seurat/reference/


### ============================================================================
### Script outlines

# Sections:
# 1. Setup the Seurat Object
# 2. Cell Filtering Based on QC
# 3. Data Normalization
# 4. Dimensional Reduction with Optional Integration
# (Additional sections to be added)


### ============================================================================
### 1. Setup the Seurat Object

## For loading 10X scRNA-Seq data, the following functions will be most relevant:  
# -`Read10X()`: primary argument is a directory from CellRanger containing the 
#               matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv files.  
# -`Read10X_h5()`: reads the hdf5 file from 10X CellRanger (scATAC-Seq or scRNA-Seq).  
# -`ReadMtx()`: read from local or remote (arguments for the .mtx file, 
#               cells/barcodes file, and features/genes file). 
#               This option is not just for 10X data, 
#               but is useful for most platforms and pipelines. 

## Importing data in a .tsv or .csv  
# If the available count data was saved as a `.txt`, `.tsv`, or `.csv` file, 
# you can import these data using regular base R import functions 
# (e.g., `read.table()`, `read.csv()`). 
# The data does not need to be in a sparse format to create a Seurat object. 
# However, if you would like to save memory upon import by loading as a sparse data matrix, 
# see `readSparseCounts()`(https://rdrr.io/bioc/scuttle/man/readSparseCounts.html)
# from the `scuttle` package (https://www.bioconductor.org/packages/release/bioc/html/scuttle.html). 

## Alternative method of loading data
# `readMM()` function from the `Matrix` package: convert standard matrix into a sparse matrix
# https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/lessons/readMM_loadData.md



## Using filtered_feature_bc_matrix or raw_feature_bc_matrix?
# The filtering performed by Cell Ranger when generating 
# the filtered_feature_bc_matrix is often good, but should consider a few things:
# - If the data has very good quality, then Cell Ranger can remove high quality cells
# - Might want to explore your own data while taking into account the biology
#   of the experiment for applying thresholds during filtering.
# - For example, if you expect a particular cell type in the dataset to be smaller
#   and/or not as transcriptionally active as other cell types, these cells could
#   be filtered out. However, Cell Ranger v3 have now tried to account for cells
#   of different sizes (e.g. tumor vs. infiltrating lymphocytes).

## Example data used in this tutorial
# The data originate from Kang et al., 2017 (https://www.nature.com/articles/nbt.4042)
# and were obtained via the Harvard Biostatistic Core GitHub repository 
# (https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/lessons/03_SC_quality_control-setup.md).
# The count matrix available from GEO (GSE96583) did not include mitochondrial reads, 
# so the BAM files were downloaded from the SRA (SRP102802), converted back to FASTQ files, 
# and processed with Cell Ranger to generate the count data used in this tutorial.


# Define data directories
data_dir <- "data"

# Use BPCells to reduce memory usage (enable only if needed)
use_bp <- FALSE
if (use_bp) {
  bp_dir <- "bpcounts"
  dir.create(bp_dir, showWarnings = FALSE)
}

# List .h5 files and Matrix Market folders
h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)
mtx_dirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

# Remove folders with matching .h5 files to prevent duplicate processing
mtx_dirs <- mtx_dirs[!basename(mtx_dirs) %in% file_path_sans_ext(basename(h5_files))]

# Initialize empty list to store Seurat objects
list_seurat <- list()


# Unified handler function
handle_sample <- function(sample, mat, 
                          min.features = 200, 
                          min.cells = 3,
                          replace.text = "_raw_feature_bc_matrix") {
  if (use_bp) {
    bp_path <- file.path(bp_dir, sample)
    write_matrix_dir(mat = mat, dir = bp_path)
    mat <- open_matrix_dir(dir = bp_path)
  }
  
  # Simplified sample name from file name
  clean_sample <- gsub(replace.text, "", sample)
  
  # Initialize the Seurat object with the raw (non-normalized data).
  seurat_obj <- CreateSeuratObject(counts = mat,
                                   min.features = min.features,
                                   min.cells = min.cells,
                                   project = clean_sample)
  # `min.cells = 3` will filter genes / features that are not present across a minimum of 3 cells
  # `min.feature=200` will filter cells that do not contain a minimum of 200 genes / features.
  # Harvard bioinformatic core also uses min.features=100
  
  
  
  list_seurat[[clean_sample]] <<- seurat_obj
}

# Process .h5 files
for (file in h5_files) {
  sample <- file_path_sans_ext(basename(file))
  mat <- if (use_bp) open_matrix_10x_hdf5(file) else Read10X_h5(file)
  handle_sample(sample, mat)
}

# Process Matrix Market directories
for (dir in mtx_dirs) {
  sample <- basename(dir)
  mat <- Read10X(data.dir = dir)
  handle_sample(sample, mat)
}


# Quick look at the list
list_seurat


## Apply a function to each Seurat object in a list
#'
#' Iterates over a named list of Seurat objects and applies either:
#' - A single (global) function to all objects, or
#' - Sample-specific functions provided as a named list.
#'
#' Two execution modes are available:
#' - `"inspect"`: Apply the function and print a preview (via `head()`) of the 
#'   result for each sample (e.g., for checking metadata or expression values).
#' - `"apply"`: Apply the function and return the full results as a list.
#'
#' @param seurat_list A named list of Seurat objects. 
#'                    Names must be unique and non-NULL.
#' @param fn A function to apply to each Seurat object, or a named list of functions 
#'           corresponding exactly to `names(seurat_list)`.  
#'           - If a single function is provided, it will be applied to all objects.  
#'           - If a list is provided, each function will be applied to its matching object.  
#'           Each function should return a printable or storable result.
#' @param n Integer. Number of rows to display per sample in `"inspect"` mode.
#' @param mode Character. Either `"inspect"` (print previews) or `"apply"` (return list).
#' @param show_progress Logical. If `TRUE`, show a progress bar 
#'                      (requires the **progressr** package).
#'                      
#' @return  
#' - If `mode = "inspect"`: prints output previews for each object.
#' - If `mode = "apply"`: Returns a named list of results for all samples.
#'
#' @details  
#' - Names in `seurat_list` and `fn` (if provided as a list) must match exactly.  
#'
#' @examples
#' # Apply the same function to all Seurat objects and preview results:
#' apply_seurat_list(seurat_list, function(obj) obj@meta.data, mode = "inspect")
#'
#' # Apply sample-specific functions and return a list:
#' fn_list <- list(sample1 = my_func1, sample2 = my_func2)
#' apply_seurat_list(seurat_list, fn_list, mode = "apply")

apply_seurat_list <- function(seurat_list, fn, n = 5,
                              mode = c("inspect", "apply"),
                              show_progress = FALSE) {
  ### Validate mode argument and seurat_list input ###
  ## Match and validate mode argument ("inspect" or "apply")
  # If not specified, "inspect" (the first value) is used by default.
  mode <- match.arg(mode)
  
  ## Check seurat_list
  # - Must be a list
  # - Must contain at least one object
  # - Must be named
  # - Names must be unique
  stopifnot(is.list(seurat_list), 
            length(seurat_list) > 0, 
            !is.null(names(seurat_list)),
            anyDuplicated(names(seurat_list)) == 0)
  
  
  ### Prepare function mapper ###
  # The `fn` argument can be:
  #  - a single function to apply to all samples
  #  - a named list of functions, one per sample
  
  if (!is.list(fn)) {
    # Case 1: Single function for all samples
    # Allows both anonymous and named functions
    fn_global <- rlang::as_function(fn)
    fn_map <- function(sample, obj) fn_global(obj)
    
  } else {
    # Case 2: Per-sample functions (list): apply matching function to each sample
    
    # Check names in fn and seurat_list
    if (is.null(names(fn))) {
      # If function list has no names, match by position and rename automatically
      
      if (length(fn) != length(seurat_list))
        stop("When `fn` is unnamed, its length must equal length(seurat_list).")
      
      names(fn) <- names(seurat_list) # explicit mapping by position
      warning("Functions in `fn` were automatically renamed by position to match names(seurat_list).")
      
    } else {
      
      # Named function list: validate names
      if (anyNA(names(fn)) || anyNA(names(seurat_list)))
        stop("Names in `fn` and seurat_list must not be NA.")
      if (anyDuplicated(names(fn)) > 0)
        stop("Duplicate names detected in `fn`.")
      
      # Ensure that the set of names matches exactly
      if (!setequal(names(fn), names(seurat_list)))
        stop("Names in `fn` must match names(seurat_list).")
      
      # Reorder functions to match the order of seurat_list
      fn <- fn[names(seurat_list)]
    }
    
    # Resolve all provided functions into callable form
    fn_resolved <- lapply(fn, rlang::as_function)
    fn_map <- function(sample, obj) fn_resolved[[sample]](obj)
  }
  
  
  ### Progress handling ###
  # Prepare sample names
  samples <- names(seurat_list)
  
  ## Helper function: Process a single Seurat sample
  # Applies the mapped function (`fn_map`) to one sample.
  # In "inspect" mode: prints a preview and returns NULL.
  # In "apply" mode: returns the full result for storage in the output list. 
  process_sample <- function(sample, obj) {
    tryCatch({
      result <- fn_map(sample, obj)
      if (mode == "inspect") {
        # Display a preview of each result for quick inspection
        cat("\n", strrep("-", 15), sample, strrep("-", 15), "\n", sep = "")
        if (is.data.frame(result) || is.matrix(result)) print(utils::head(result, n))
        else if (inherits(result, "Seurat")) print(utils::head(result@meta.data, n))
        else print(utils::head(result, n))
        return(NULL) # inspect mode
      } else {
        return(result) # apply mode
      }
    }, error = function(e) {
      warning(sprintf("'%s' failed: %s", sample, e$message))
      NULL
    })
  }
  
  ## Helper function: Run mapply over all samples
  # Wraps mapply with consistent arguments.
  run_apply <- function(fun) {
    mapply(FUN = fun,
          sample = samples,
          obj = seurat_list,
          SIMPLIFY = FALSE,
          USE.NAMES = TRUE)
    # For future.apply you can add: future.seed = TRUE, future.globals = FALSE
  }
  
  ### Execute with optional progress reporting ###
  if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
    # Show a progress bar if {progressr} is installed
    # Increment progress bar after each sample
    out <- progressr::with_progress({
      p <- progressr::progressor(steps = length(seurat_list))
      run_apply(function(sample, obj) { on.exit(p(), add = TRUE);
        process_sample(sample, obj) })
    })
    cat("\n")  # Ensure next console output starts on a new line
  } else {
    # No progress bar
    out <- run_apply(process_sample)
  }
  
  ### Output handling ###
  if (mode == "inspect") {
    # Already printed inside process_sample
    # Return the NULL list invisibly (no console spam)
    return(invisible(out))
  } else {
    # mode == "apply": return the list of results
    return(out)
  }
}






## Explore the metadata
# Seurat automatically creates some metadata for each of the cells when using
# `Read10X()` function to read in data. This is stored in the `meta.data` slot.
apply_seurat_list(list_seurat, fn = function(x) x@meta.data)

# Legacy code
#for (sample in names(list_seurat)){
#  cat("\n---", sample, "---\n")
#  print(head(list_seurat[[sample]]@meta.data))
#}
#head(list_seurat[[1]]@meta.data)

# It contains:
# -`orig.ident`: this contains the sample identity if known (or "SeuratProject" by default)
# -`nCount_RNA`: number of RNA molecules (UMIs) per cell (i.e. count depth)
#    - Each unique RNA molecule (non-PCR duplicates) will have its own UMI
#    - High total count: potential doublets or multiplets
#    - Low total count: potential ambient mRNA (not real cells)
#    - Cell Ranger threshold set at 500 UMIs
# -`nFeature_RNA`: number of genes detected per cell
#    - High number of detected genes: potential doublets or multiplets
#    - Low number of detected genes: potential ambient mRNA (not real cells)


## Examine a few genes in the first thirty cells
apply_seurat_list(
  list_seurat, 
  fn = function(x) GetAssayData(x, assay = "RNA", layer = "counts")
  [c("CD3D", "TCL1A", "MS4A1"), 1:30]
)

## Example output
# 3 x 30 sparse Matrix of class "dgCMatrix"
#                                                                    
# CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
# TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
# MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .

# The . values in the matrix represent 0s (no molecules detected). 
# Since most values in an scRNA-seq matrix are 0, 
# Seurat uses a sparse-matrix representation whenever possible. 
# This results in significant memory and speed savings for Drop-seq/inDrop/10x data.



## Check the number of cells or features
# Return number of cells across all layers
#num_cells <- ncol(list_seurat[[1]])
# Return number of features across all layers
#num_features <- nrow(list_seurat[[1]])
# Return row by column information
#dim(list_seurat[[1]])




### Merge all objects together into a single Seurat object ###
# This make it easier to run the QC steps for all samples together and 
# compare the data quality.
#
# NOTE: if samples are fairly different, QC thresholds should be applied 
# on a per sample basis (using list_seurat)

if (length(list_seurat) > 1) {
  merged_seurat <- merge(
    x = list_seurat[[1]],
    y = list_seurat[-1], # `-1` excludes the first element, returning everything else
    add.cell.ids = names(list_seurat)
  )
} else {
  # If the list only has one sample
  merged_seurat <- list_seurat[[1]]
}

# Check that the merged object has the appropriate sample-specific prefixes
#head(merged_seurat@meta.data)
#tail(merged_seurat@meta.data)
lapply(split(merged_seurat@meta.data, merged_seurat$orig.ident), head)


# Use `JoinLayers()` to combine counts from multiple samples (layers) into a single unified matrix.
# This is useful for downstream tools that require a flat structure, but increases memory usage.
#merged_seurat <- JoinLayers(merged_seurat)   


## Check memory usage of all objects in the environment
sapply(ls(), function(x) format(object.size(get(x)), units = "auto"))


### ============================================================================
### 2. Cell Filtering Based on QC

# Explore QC metrics and filter cells based on user-defined criteria.

## A few commonly used QC metrics:
# 1) The number of unique genes detected in each cell.
#       - Low-quality cells or empty droplets often have very few genes
#       - Cell doublets or multiplets may exhibit an aberrantly high gene count
# 2) The total number of molecules detected within a cell 
#    (correlates strongly with unique genes)
# 3) The percentage of reads that map to the mitochondrial genome
#       - Low-quality/dying cells often exhibit extensive mitochondrial contamination
#       - In Seurat, mitochondial QC is calculated with the `PercentageFeatureSet()`
#         function, which calculates the percentage of counts originating from a set of features
#       - All genes starting with "MT-" is set as mitochondrial genes
#
# NOTE: Doublets and multiplets are also not as easy to predict. 
# There are dedicated tools for this purpose (e.g., DoubletFinder, scDblFinder,
# scrublet (python package)).
# Eliminating ambient RNA is also not as straight forward. 
# Ambient RNA can be included in GEMs along with intact cells. 
# Tools that can be used to predict and remove ambient RNA (e.g., SoupX and DecontX).


### Mitochondrial Percentage ###

# Store percentage of mitochondrial genes
#
# The `PercentageFeatureSet()` function calculates the percentage of counts 
# from genes starting with the pattern "MT-" by dividing their total count 
# by the total counts for all genes, then multiplying by 100 to obtain a percentage.
#
# NOTE: The pattern ("^MT-") is specific to human mitochondrial gene names. 
# You may need to adjust the pattern argument for other organisms. 
# Additionally, this function will not work if gene names were not used as gene IDs. 
# To compute this metric manually, refer to:
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md

### For individual objects in the list
list_seurat <- apply_seurat_list(list_seurat,
                       fn = function(x){PercentageFeatureSet(x, pattern = "^MT-", 
                                                             col.name = "percent.mt")},
                       mode = "apply")
# Quick view after storing percentage of mitochondial genes
apply_seurat_list(list_seurat, fn = function(x) x@meta.data)

# Alternative:
#list_seurat <- apply_seurat_list(list_seurat, 
#                                 fn = function(x){x[["percent.mt"]] <- 
#                                   PercentageFeatureSet(x, pattern = "^MT-")
#                                 return(x)}, mode = "apply")

#apply_seurat_list(list_seurat, fn = function(x) x@meta.data)




### For merged object
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
# Alternative 1:
# merged_seurat <- PercentageFeatureSet(merged_seurat, pattern = "^MT-", col.name = "percent.mt")
#
# Alternative 2 (without using the `PercentageFeatureSet` function:
# For example, when not working with human gene or not using gene names as ID
# https://github.com/hbctraining/scRNA-seq/blob/master/lessons/mitoRatio.md

# The number of unique genes and total molecules are calculated during `CreateSeuratObject()`
# Can find them stored in the object mneta data
# Show QC metrics for the first 5 cells
head(merged_seurat@meta.data, 5)



### Novelty score ###
# https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/lessons/04_SC_quality_control.md
# Divide the log10 number of genes by the log10 number of UMIs
# 
# The novelty score reflects the complexity of RNA species in a cell.
# If a cell has a high number of UMIs (nCount_RNA) but 
# a low number of detected genes (nFeature_RNA),
# it suggests that only a few genes were captured and repeatedly sequenced,
# indicating low transcriptomic complexity.
# Such low-novelty cells may correspond to specific cell types (e.g., red blood cells)
# that lack a typical transcriptome, or they may result from artifacts or contamination.
# As a general guideline, high-quality cells typically have a novelty score above 0.8.

### For individual objects in the list
list_seurat <- apply_seurat_list(list_seurat, 
                                 fn = function(x){x$log10GenesPerUMI <- 
                                   log10(x$nFeature_RNA)/ log10(x$nCount_RNA)
                                   return(x)},
                                 mode = "apply")
# Quick view after calculating Novelty score
apply_seurat_list(list_seurat, fn = function(x) x@meta.data)

### For merged object
merged_seurat$log10GenesPerUMI <- 
  log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)

head(merged_seurat)
#head(merged_seurat@meta.data[merged_seurat$orig.ident == "stim", ])



### Metadata can be manipulated as a standard data frame ####################

### DO NOT RUN: Only applicable to the example dataset used in this tutorial.

## Create metadata dataframe
metadata <- merged_seurat@meta.data

## Add cell IDs to metadata
metadata$cells <- rownames(metadata)

## Create sample column
# Each cell ID has a prefix such as `ctrl_` or `stim_`, 
# which can be used to create a new column indicating 
# the condition each cell belongs to.
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"

# Rename columns to be more intuitive
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file="data/merged_filtered_seurat.RData")

### More ways to add metadata (from NIH seminar) ############################

### DO NOT RUN

# Add condition to metadata (e.g. wildtype or double knockout)
merged_seurat$condition <- ifelse(str_detect(merged_seurat@meta.data$orig.ident, "^W"),
                        "WT","DKO")

# Add time information to metadata (e.g. Day 0 or Day 6)
merged_seurat$time_point <- ifelse(str_detect(merged_seurat@meta.data$orig.ident, "0"),
                         "Day 0","Day 6")

# Add Condition + Time to the metadata
merged_seurat$cond_tp <- paste(merged_seurat$condition, merged_seurat$time_point)

############################################################################


### Visualize QC metrics and filter cells ###

### Cell counts

# Cell counts are determined by the number of unique cellular barcodes 
# (e.g., the column names of the Seurat object or the row names of its metadata).
# 
# In this example dataset, approximately 12,000–13,000 cells were loaded. 
# However, only a fraction of cells are typically captured during library prep 
# (e.g., 50–60% for 10X Genomics and 70–80% for inDrops).
# 
# Cell numbers can vary by protocol and may appear higher than expected.
#  - In inDrops, some hydrogels may contain multiple barcodes.
#  - In 10X, GEMs can contain barcodes without cells.
# These artifacts, along with dying cells, can lead to inflated barcode counts.
# 
# In this example dataset, each sample contains approximately 15,000 cells, 
# which exceeds the expected range of 12,000 to 13,000 cells. 
# Low-quality cells or artifacts may have contributed to the inflated cell count.
# 
# NOTE: Cell concentration for library preparation should be measured using 
# a hemocytometer or an automated cell counter, rather than a FACS machine 
# or Bioanalyzer, as these are not accurate for determining cell concentration.


## Visualize the number of cell counts per sample
merged_seurat@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar(color="black") +
  stat_count(geom = "text", colour = "black", size = 3.5, 
             aes(label = after_stat(count)),
             position=position_stack(vjust=0.5))+
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of Cells per Sample")



### UMI counts (transcripts) per cell (nCount_RNA)

# The UMI counts should generally be at least above 500.
# Counts between 500 and 1,000 are acceptable but suggest under-sequencing.

## Visualize the number UMIs/transcripts per cell
merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("nUMI") +
  geom_vline(xintercept = 500, color = "grey20", linetype = "dotted")
# Most samples in the example data have 1000 UMIs or greater



### Genes detected per cell (nFeature_RNA)

# A high-quality dataset typically shows a single major peak in the histogram, 
# representing the encapsulated cells. A shoulder to the left of the major peak 
# or a bimodal distribution may indicate failed cells, 
# distinct biological populations (i.e. quiescent cell), 
# or differences in cell size or complexity.

## Visualize the distribution of genes detected per cell via histogram
merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("nGene") +
  geom_vline(xintercept = 300, color = "grey20", linetype = "dotted")



### Complexity (novelty score)

# The novelty score is computed by taking the ratio of nGenes over nUMI.

## Visualize the overall complexity of the gene expression 
# by visualizing the genes detected per UMI (novelty score)
merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=log10GenesPerUMI, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("log10GenesPerUMI") +
  geom_vline(xintercept = 0.8, color = "grey20", linetype = "dotted")


### Mitochondrial percentage

# Check for high mitochondrial content, which can signal dead or dying cells

## Visualize the distribution of mitochondrial gene expression detected per cell
merged_seurat@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  xlab("Mitochondrial Percentage") +
  geom_vline(xintercept = 20, color = "grey20", linetype = "dotted")


### Joint Visualization ###

## Visualize QC metrics as a violin plot
VlnPlot(merged_seurat, 
        features = c("nFeature_RNA", "nCount_RNA", 
                     "percent.mt", "log10GenesPerUMI"),
        layer = "counts", ncol = 2, alpha = 0.2, group.by = "orig.ident") & 
  theme(axis.title.x = element_blank())


## Visualize the correlation between genes detected and number of UMIs 
# and determine whether strong presence of cells with low numbers of genes/UMIs
merged_seurat@meta.data %>% 
  ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt),shape=21,alpha=0.4) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm, color = "blue") +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500, color = "grey20", linetype = "dashed") +
  geom_hline(yintercept = 250, color = "grey20", linetype = "dashed") +
  facet_grid(.~orig.ident)
# Cells that are poor quality are likely to have low genes and UMIs per cell, 
# and correspond to the data points in the bottom left quadrant of the plot.
#
# Data points in the bottom right hand quadrant of the plot have a high number 
# of UMIs but only a few number of genes. These could be dying cells, but also 
# could represent a population of a low complexity celltype (i.e red blood cells).
#
# Cells with high mitochondrial percentages are shown as darker-colored data points. 
# If these cells also have low UMI counts and few detected genes, 
# it may indicate damaged or dying cells whose cytoplasmic mRNA has leaked out 
# through compromised membranes, leaving mostly mitochondrial mRNA preserved.



## FeatureScatter is commonly used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, 
# i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2




### Filtering ###

### Filter threshold examples:
# Filtering thresholds vary depending on the dataset. Here are a few examples:
# - In the Seurat vignette, cells are filtered if they have fewer than 200 or 
#   more than 2,500 detected features, or if mitochondrial content exceeds 5%.
# - In the Harvard Bioinformatics Core tutorial, the following criteria are used:
#   (nCount_RNA >= 500) & (nFeature_RNA >= 250) & 
#   (log10GenesPerUMI > 0.80) & (percent.mt < 20)
# - In the NIH CCR seminar series, thresholds are tailored per sample. For example:
#     - Day 0 samples:     nFeature_RNA > 350, nCount_RNA > 650, percent.mt < 10
#     - Day 6 samples:     nFeature_RNA > 350, nCount_RNA > 650, percent.mt < 25



## Apply quality filtering with universal thresholds or on a per-sample basis
# It is ideal to use universal thresholds to apply quality filtering on all samples.
# However, filtering should be considered to be adjusted at a per-sample level

filter.method <- "individual" # options: "universal" or "individual"

if (filter.method == "individual"){
  message("Applying individual filtering thresholds to each sample.")
  
  # Define individual filtering thresholds for each sample
  filters <- list(
    "ctrl" = function(x) subset(x, subset = (nCount_RNA >= 500) & 
                                     (nFeature_RNA >= 250) &
                                     (log10GenesPerUMI > 0.80) &
                                     (percent.mt < 5)),
    "stim" = function(x) subset(x, subset = (nCount_RNA >= 200) & 
                                     (nFeature_RNA >= 250) &
                                     (log10GenesPerUMI > 0.80) &
                                     (percent.mt < 10))
  )
  
  
  
  ## Apply sample-specific filtering to each Seurat object
  list_filtered <- apply_seurat_list(seurat_list = list_seurat, 
                                     fn = filters, mode = "apply")
  
  
  ## Merge all filtered objects into one
  filtered_seurat <- merge(
    x = list_filtered[[1]],
    y = list_filtered[-1], # `-1` excludes the first element, returning everything else
    add.cell.ids = names(list_filtered)
  )
  
} else {
  if (filter.method != "universal"){
    warning("Invalid filtering method specified. 
            Using the default: universal threshold.")
  }
  message("Applying universal filtering thresholds to all samples.")
  
  ## Using subset.Seurat to filter the merged Seurat object
  filtered_seurat <- subset(x = merged_seurat, 
                            subset= (nCount_RNA >= 500) & 
                              (nFeature_RNA >= 250) & 
                              (log10GenesPerUMI > 0.80) & 
                              (percent.mt < 20))
  
  ## Split the merged object by orig.ident
  list_filtered <- SplitObject(filtered_seurat, split.by = "orig.ident")
  # Rename project.name for each split object
  for(name in names(list_filtered)){
    list_filtered[[name]]@project.name <- name
  }
}


### Combine metadata from original and filtered Seurat objects
meta_orig <- merged_seurat@meta.data %>%
  mutate(dataset = "Original")

meta_filtered <- filtered_seurat@meta.data %>%
  mutate(dataset = "Filtered")

meta_combined <- bind_rows(meta_orig, meta_filtered)


# Bin percent.mt into categorical intervals
meta_combined$mt_bin <- cut(meta_combined$percent.mt, 
                            breaks = c(-Inf, 5, 10, 20, Inf), 
                            labels = c("<5%", "5–10%", "10–20%", ">20%"))

# Quick summary of cell counts in each mitochondrial percentage bin
table(meta_combined$mt_bin)


### Plot QC metrics to compare original and filtered Seurat objects.
QCplot.style <- "gradient" 
# Options: "gradient" (continuous color scale for % mitochondrial content)
#          "bin"      (categorical color bins for % mitochondrial content)

if (QCplot.style == "gradient"){
  
  message("Plotting QC metrics using a continuous color gradient for percent.mt")
  
  # Use a continuous gradient to represent percent.mt
  ggplot(meta_combined, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt), 
         shape = 21, alpha = 0.2) +
    geom_point() +
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm, color = "blue") +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500, color = "grey20", linetype = "dashed") +
    geom_hline(yintercept = 250, color = "grey20", linetype = "dashed") +
    facet_grid(dataset ~ orig.ident) +
    labs(title = "QC metrics before and after filtering")
  
} else if (QCplot.style == "bin"){
  
  message("Plotting QC metrics using binned categories for percent.mt")
  
  # Use categorical color scheme for binned mitochondrial percentages
  ggplot(meta_combined, aes(x = nCount_RNA, y = nFeature_RNA, color = mt_bin)) +
    geom_point(shape = 21, stroke = 1.5, alpha = 0.5, size = 2) +
    scale_color_manual(
      values = c("<5%"      = "grey50",  # medium gray
                 "5–10%"    = "grey10",  # dark gray
                 "10–20%"   = "#fb6a4a",  # light red
                 ">20%"     = "#a50f15"   # dark red
      ),
      name = "% MT"
    ) +
    stat_smooth(method = lm, color = "blue") +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500, color = "grey20", linetype = "dashed") +
    geom_hline(yintercept = 250, color = "grey20", linetype = "dashed") +
    facet_grid(dataset ~ orig.ident) +
    labs(title = "QC metrics before and after filtering")
  
} else {
  stop("Invalid plot style specified. Please use either 'gradient' or 'bin'.")
}




### ============================================================================
### 3. Data Normalization

# 10X genomics' normalization guide of existing methods: 
# https://www.10xgenomics.com/analysis-guides/single-cell-rna-seq-data-normalization

## Normalization is the process of adjusting raw count values to account for
# technical factors that influence variation in mRNA counts, while preserving
# true biological differences. This ensures that expression levels are
# comparable across genes and/or samples.
#
# The main source of technical variation in scRNA-seq is sequencing depth.
# Cells with greater depth may appear to have higher expression due to more
# sampling. Potential causes of variation in sequencing depth include:
#
# - Library size: differences in cell lysis efficiency, mRNA capture, 
#   amplification, ambient RNA contamination, and batch effects
# - Cell-specific RNA content: individual cells naturally contain different total
#   amounts of RNA, which can vary due to cell size, cell cycle stage, or 
#   metabolic activity
# - Gene length: longer genes may appear to have higher expression than shorter
#   ones. However, in scRNA-seq, especially with UMI-based technologies, gene
#   length has less influence because each transcript is ideally counted once via
#   its UMI, regardless of length (typically only the 5' or 3' end is sequenced).
#   In contrast, full-length protocols may require normalization for transcript
#   length.


### Normalize a Seurat Object Using LogNormalize or SCTransform Workflow
#'
#' This function normalizes a single Seurat object using either the traditional
#' log normalization workflow or the SCTransform workflow. The choice of method
#' depends on analysis goals and dataset characteristics.
#' 
#' @param seurat.object A Seurat object containing raw counts in the RNA assay.
#' @param norm.method   A character string specifying the normalization method:
#'                      - "LogNormalize": The traditional workflow in Seurat, 
#'                               which uses NormalizeData -> FindVariableFeatures 
#'                               -> ScaleData to rescale counts, identify variable 
#'                               genes, and standardize expression values.
#'                      - "SCT": This workflow uses `SCTransform` function, 
#'                                       which replaces the traditional steps by
#'                                       fitting a regularized negative binomial 
#'                                       model per gene to perform normalization, 
#'                                       variance stabilization, and regression 
#'                                       of technical effects in a single step.
#'                      Default = "SCT".
#'
#' @return A Seurat object with normalized data stored in either:
#'           - RNA assay (log method)
#'           - SCT assay (sctransform method)
#'           
#' @note 
#' - This function accepts a single Seurat object as input
#' - Modify parameters, such as `vars.to.regress`, within the function to 
#'   account for unwanted data-specific variation (e.g., mitochondrial content, 
#'   cell cycle effects)
norm.seurat <- function(seurat.object, norm.method = "SCT"){
  
  ## Validate normalization method
  if (!norm.method %in% c("SCT", "LogNormalize")) {
    warning("Invalid normalization method; using default 'sctransform'.")
    norm.method <- "SCT"
  }
  
  
  
  if (norm.method=="LogNormalize"){
    
    ### Traditional normalization method ###
    
    message("Applying log normalization method...\n")
    
    ## Normalization
    # By default, Seurat uses a global-scaling normalization method "LogNormalize".
    # It normalizes the feature expression measurements for each cell by the
    # total expression, multiples this by a scale factor (10,000 by default), 
    # and then log-transforms the result.
    # In Seurat v5, normalized values are stored in seurat.object[["RNA"]]$data.
    seurat.object <- NormalizeData(seurat.object, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)
    
    ## Identify highly variable features
    # Identify a subset of genes/features that exhibit high cell-to-cell variation.
    # (i.e., they are highly expressed in some cells, and lowly expressed in others).
    # These are often more informative for downstream analyses because they highlight 
    # biological differences between cells.
    # By default, 2000 features are returned per dataset.
    seurat.object <- FindVariableFeatures(seurat.object, 
                                          selection.method = "vst", 
                                          nfeatures = 2000)
    
    # Identify the 10 most variable genes
    top10 <- head(VariableFeatures(seurat.object), 10)
    
    # Plot variable features with and without labels
    #plot1 <- VariableFeaturePlot(seurat.object)
    #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    #print(plot1 + plot2)
    
    ## Scaling the data
    # Apply a linear transformation ("scaling") using `ScaleData()`
    # - Shifts the expression of each gene, so that:
    #     - The mean expression across cells is 0
    #     - The variance across cells is 1. This gives equal weight in downstream
    #       analysis, so that highly-expressed genes do not dominate
    # - The results of this are stored in `seurat.object[["RNA"]]$scale.data`
    # - By default, only variable features are scaled. Can scale additional features
    #   by specifying the `features` argument
    # - Can also use the `vars.to.regress` argument to remove unwanted sources
    #   of variation. For example, we could regress out heterogeneity associated
    #   with cell cycle stage or mitochondrial contamination.
    all.genes <- rownames(seurat.object)
    seurat.object <- ScaleData(seurat.object, 
                               features = all.genes,
                               #vars.to.regress = "percent.mt"
                               )
    
    
    # NOTE: While this workflow of normalization is standard and widely used,
    # global-scaling relies on an assumption that each cell originally contains
    # the same number of RNA molecules.
    
    
  } else if (norm.method=="SCT"){
    
    ### SCTransform ###
    # https://satijalab.org/seurat/articles/sctransform_vignette
    
    message("Applying SCTransform workflow...\n")
    
    # The use of `SCTransform` replaces the need to run `NormalizeData`, 
    # `FindVariableFeatures`, or `ScaleData`
    #
    # Compared to traditional log normalization:
    # - Does not assume equal RNA content per cell as the traditional workflow
    # - Better accounts for technical factors such as sequencing depth
    # - Uses Pearson residuals for transformation, giving more weight to lowly expressed
    #   but biologically informative genes, and reducing the dominance of broadly
    #   expressed high-abundance genes.
    #
    # Additional notes:
    # - Transformed data will be available in the SCT assay,
    #   which is set as the default after running sctransform
    # - During normalization, can also remove confounding sources of variation,
    #   for example, mitochondrial mapping percentage
    # - Can apply different SCT version by setting vst.flavor = 'v1'. (default is v2)
    # - The `glmGamPoi` package substantially improves speed.
    
    
    # Run sctransform
    seurat.object <- SCTransform(seurat.object, 
                                 vars.to.regress = "percent.mt", 
                                 verbose = TRUE)
    
    
    # The results of sctransform are stored in the "SCT" assay.
    # - `seurat.object[["SCT"]]$scale.data` contains the residuals (normalized values), 
    #   and is used directly as input to PCA
    #   - Note: this matrix is non-sparse, thus can take up lot of memory if stored for all genes
    #   - To save memory, only store values of variable genes (return.only.var.genes = TRUE)
    #     by default in `SCTransform()`
    # - The 'corrected' UMI counts are stored in `seurat.object[["SCT"]]$counts` and
    #   log-normalized versions of these corrected counts in `seurat.object[["SCT"]]$data`
  }
  
  
  
  ## Return normalized object
  return(seurat.object)
}



### Specify input format: merged object or individual samples
input.format <- "merged"
# options: "merged" (normalize the combined Seurat object containing all samples)
#          "individual" (normalize each sample/Seurat object separately)

## Note: In the Harvard Bioinformatics Core and NIH CCR seminar tutorials,
#        normalization was performed on a merged Seurat object containing all samples.
#        In practice, whether to normalize after merging or individually before 
#        merging depends on the analysis goals and the presence of batch effects
#        (an ongoing discussion). 
#        If batch effects are minimal, merging first and normalizing 
#        together is acceptable. However, if batch effects are expected or if
#        using Seurat’s integration workflow, each sample should be normalized
#        separately before running integration:
#        https://satijalab.org/seurat/articles/integration_introduction
#
# References:
# https://github.com/satijalab/seurat/issues/7407
# https://github.com/satijalab/sctransform/issues/182
# https://github.com/satijalab/sctransform/issues/55
# https://github.com/satijalab/seurat/issues/6116
# https://satijalab.org/seurat/archive/v4.3/sctransform_v2_vignette


### Specify normalization method (option: "LogNormalize" or "SCT")
specify_norm_method = "SCT"

if (!(specify_norm_method %in% c("LogNormalize", "SCT"))) {
  stop("Invalid normalization method. Please choose 'LogNormalize' or 'SCT'.")
}



### Apply normalization
if (input.format=="merged"){
  
  message("Normalizing the merged Seurat object.\n")
  
  normalized_seurat <- norm.seurat(filtered_seurat, 
                                   norm.method = specify_norm_method)
  
} else if (input.format=="individual"){
  
  message("Normalizing individual Seurat objects separately.\n")
  
  list_normalized <- apply_seurat_list(seurat_list = list_filtered, 
                    fn = function(x){
                      norm.seurat(x, norm.method = specify_norm_method)},
                    mode = "apply")
  
  ## Merge all normalized objects into one
  # NOTE: This merge is for visualization purposes only. Please run the 
  #       integration workflow before proceeding with downstream analysis.
  normalized_seurat <- merge(
    x = list_normalized[[1]],
    y = list_normalized[-1],
    add.cell.ids = names(list_normalized))
  
} else{
  
  stop("Invalid input format specified. Please use either 'merged' or 'individual'.")
  
}


### Compare nCount_RNA before and after normalization

# nCount_RNA: Total number of transcripts (UMIs) detected per cell.
#             This metric reflects sequencing depth or cellular RNA content.
# 
# nCount_SCT: Total corrected UMI counts in the SCT assay
# 
# After normalization: To confirm that technical variability, such as
#   differences in sequencing depth, has been mitigated.
#   A more uniform distribution across samples or groups indicates that
#   normalization was effective in reducing technical bias.
#
# NOTE: Both normalization methods store normalized data separately from the raw data.
# - For the traditional workflow (NormalizeData), the log-normalized expression
#   is stored in the "data" layer of the RNA assay.
# - For SCTransform, the normalized data is stored in a separate assay called "SCT".
#   - "counts": corrected UMI counts
#   - "data": log-transformed corrected counts
#   - "scale.data": Pearson residuals used for PCA
# - The raw (untransformed) counts remain unchanged in the "counts" layer of the RNA assay.


if (specify_norm_method == "LogNormalize"){
  
  ## Specify gene of interest
  specify_gene <- "GAPDH"
  
  ## Before normalization
  qc_before <- VlnPlot(normalized_seurat, features = specify_gene, 
                       layer = "counts") & 
    theme(axis.title.x = element_blank()) &
    ggtitle("Before normalization") & 
    ylab(paste0(specify_gene, " Expression Level"))
  
  ## After normalization
  qc_after <- VlnPlot(normalized_seurat, features = specify_gene, 
                      layer = "data") & 
    theme(axis.title.x = element_blank()) &
    ggtitle("After normalization") & 
    ylab(paste0(specify_gene, " Expression Level"))
  
} else if (specify_norm_method == "SCT"){
  
  ## QC plot before SCTransform
  qc_before <- VlnPlot(
    filtered_seurat,
    features = c("nCount_RNA"),
    ncol = 3, pt.size = 0.1, layer = "counts", alpha = 0.2
  ) & theme(axis.title.x = element_blank()) &
    ggtitle("Before normalization") & ylab("nCount_RNA")
  
  ## QC plot after SCTransform
  qc_after <- VlnPlot(
    normalized_seurat,
    features = c("nCount_SCT"),
    ncol = 3, pt.size = 0.1, layer = "counts", alpha = 0.2) & 
    theme(axis.title.x = element_blank()) &
    ggtitle("After normalization")  & ylab("nCount_SCT")
}


## Combine before & after QC plots
qc_before / qc_after
# Patchwork uses "|" for side-by-side, "/" for vertical layout



# Add Cell cycle scoring:
# https://satijalab.org/seurat/articles/cell_cycle_vignette
# https://github.com/satijalab/seurat/issues/1679
# https://github.com/satijalab/seurat/issues/3692
# https://github.com/hbctraining/Intro-to-scRNAseq/blob/master/lessons/06_SC_SCT_normalization.md

# Integration
# https://satijalab.org/seurat/articles/integration_introduction
# https://satijalab.org/seurat/articles/seurat5_integration



### ============================================================================
### 4. Dimensional Reduction with Optional Integration

## Dimensional reduction is the process of projecting this data into a smaller 
# set of variables that capture the most important patterns of variation. 
# In Seurat, this typically means running PCA to summarize global variation, 
# followed by UMAP or t-SNE to visualize cells in two dimensions.

## When working with multiple samples or batches, each dataset may contain both 
# biological signal and technical variation (e.g., sequencing depth, batch 
# effects). Integration refers to methods that align or “harmonize” multiple 
# datasets so that cells with similar expression profiles cluster together 
# regardless of their origin.
#
# Integration is useful when your scientific question requires joint analysis 
# across samples—for example, comparing disease vs. control, or pooling 
# patients to identify shared cell types. Without integration, technical 
# effects can dominate dimensional reduction, leading to clustering by sample 
# rather than by biology. Integration helps remove these unwanted sources of 
# variation so that downstream analyses reflect meaningful biological structure.
#
# It’s good practice to run dimensional reduction before integration to 
# understand what each dataset looks like on its own. This helps answer key 
# questions: Are major cell types separable? Is the clustering reasonable? 
# Do you see strong batch effects? This diagnostic step helps decide whether 
# integration is necessary. 
#
# Importantly, integration is not always needed. If datasets already overlap 
# well, integration may not add value and can even risk over-correcting and 
# removing real biological differences. Conversely, if strong batch effects 
# are evident, integration can be justified.


## SCTransform with v5 integration
# https://jasonbiology.tokyo/2024/05/06/seurat-v5-sctransform-data-integration/
# https://github.com/satijalab/seurat/issues/7251


### Perform dimensional reduction
#'
#' This function performs dimensional reduction on a single Seurat object 
#' using either a linear or nonlinear workflow. The linear option runs PCA 
#' to identify major axes of variation, while the nonlinear option performs 
#' PCA followed by UMAP or t-SNE for visualization. Clustering is applied 
#' after PCA to identify transcriptionally similar groups of cells. Wrapping 
#' these steps into a single function keeps the dimensional reduction workflow 
#' organized and makes it easy to rerun with consistent parameters if 
#' integration or other downstream analyses require repeating the reduction.
#' 
#' @param seurat.object A Seurat object as input.
#' @param projection.method A character string specifying the dimensional 
#' reduction workflow:
#'     - "linear": Perform PCA only (Default)
#'     - "nonlinear": Perform PCA followed by UMAP or t-SNE
#'
#' @param reduction.use A character string specifying the reduction layer 
#' or integration result to use:
#'    - "pca": Use PCA without integration (Default)
#'    - "integrated.cca": Use PCA results from CCA-based integration
#'    - "integrated.rpca": Use PCA results from RPCA-based integration
#'    - "harmony": Use Harmony-corrected embeddings
#'    - "integrated.mnn": Use MNN (mutual nearest neighbor) integration
#'    - "integrated.scvi": Use scVI latent space
#' 
#' @param umap.name A character string specifying the reduction.name in 
#' `RunUMAP`:
#'    - "umap": Default setting for RunUMAP.
#'
#' @param cluster.name A character string specifying the cluster.name in 
#' `FindClusters`:
#'    - "seurat_clusters": Default setting for `FindClusters`. This will be
#'                         overwritten everytime `FindClusters` is run.
#' 
#' @param dims.use Integer vector specifying which principal components to
#' use for downstream analysis. These PCs are used for neighbor graph
#' construction, clustering, and UMAP. Typically, the first 20–50 PCs capture
#' the most meaningful biological variation. Default: 1:30.
#' 
#' @param resolution Numeric value that controls the granularity of the
#' clustering result. Higher values lead to a larger number of smaller
#' clusters, while lower values result in fewer, broader clusters. For datasets 
#' with ~3,000 cells, values between 0.4 and 1.2 yield reasonable results. 
#' Consider increasing this value for larger datasets. Default: 0.8.
#'
#' @return A Seurat object with dimensional reduction data including:
#'    - PCA (linear method)
#'    - UMAP or t-SNE (nonlinear method)
#'           
#' @note 
#' - This function accepts a single Seurat object as input
#' - If `reduction.use != "pca"`, that reduction must already exist
#'   (e.g., created via `IntegrateLayers()`).
#' - Adjust PCs (`dims.use`), clustering resolution, and UMAP parameters as needed.

DimRed.seurat <- function(seurat.object, 
                          projection.method = "nonlinear",
                          reduction.use = "pca",
                          umap.name = "umap",
                          cluster.name = "seurat_clusters",
                          dims.use = 1:30,
                          resolution = 0.8){
  
  ### PCA (linear dimensional reduction) ###
  if (reduction.use == "pca"){
    obj <- RunPCA(seurat.object)
  } else {
    obj <- seurat.object
    
    if (!(reduction.use %in% Reductions(obj))) {
      stop(sprintf(
        "Requested reduction '%s' not found. Available reductions: %s. ",
        reduction.use, paste(Reductions(obj), collapse = ", ")
      ), call. = FALSE)
    }
  }
  
  
  
  # Visualize PCA
  #print(DimPlot(obj, label = F))
  
  
  ## Return PCA result if only asking for linear reduction
  if (projection.method == "linear"){
    message("Return linear result")
    return(obj)
  } else if (projection.method != "nonlinear") {
    warning("Invalid projection method; defaulting to 'nonlinear'.")
    projection.method <- "nonlinear"
  }
  
  
  ### Clustering ###
  
  # Seurat applies a graph-based clustering approach that embed cells in a graph structure,
  # for example a K-nearest neighbor (KNN) graph, with edges drawn between cells 
  # with similar feature expression patterns, and attempt to partition this graph 
  # into highly interconnected 'quasi-cliques' or 'communities'
  
  # First, construct a KNN graph based on the Euclidean distance in PCA space, 
  # and refine the edge weights between any two cells based on the shared overlap 
  # in their local neighborhoods (Jaccard similarity).
  obj <- FindNeighbors(obj, dims = dims.use, verbose = TRUE, reduction = reduction.use)
  
  # Next, apply modularity optimization such as Louvain algorithm (default) or SLM,
  # to iteratively group cells together, with the goal of optimizing 
  # the standard modularity function.
  # The `FindClusters()` function implements this procedure, and contains 
  # a resolution parameter that sets the 'granularity' of the downstream clustering,
  # with increased values leading to a greater number of clusters.
  # Setting this parameter between 0.4-1.2 typically returns good results for around 3K cells.
  # Optimal resolution often increases for larger datasets
  obj <- FindClusters(obj, resolution = resolution, verbose = TRUE,
                      cluster.name = cluster.name)
  
  # Look at cluster IDs of the first 5 cells
  head(Idents(obj), 5)
  
  
  ### Non-linear dimensional reduction (UMAP/tSNE)
  # To learn the underlying structure in the dataset,
  # in order to place similar cells together in low-dimensional space
  
  # These methods aim to preserve local distances in the dataset 
  # (i.e. ensuring that cells with very similar gene expression profiles co-localize),
  # but often do not preserve more global relationships
  # Note: avoid drawing biological conclusions solely based on visualization techniques
  
  # Run UMAP
  obj <- RunUMAP(obj, dims = dims.use, verbose = TRUE, 
                 reduction = reduction.use, reduction.name = umap.name)
  
  #print(DimPlot(obj, 
  #        reduction = "umap", 
  #        label = TRUE,
  #        group.by = c("stim", "seurat_clusters")
  #        #group.by = "orig.ident"
  #        ))
  
  message("Return nonlinear result")
  return(obj)
  
}

### Check clustering results before deciding on integration
obj <- DimRed.seurat(normalized_seurat, projection.method = "nonlinear",
                     dims.use = 1:30, resolution = 2,
                     cluster.name = "unintegrated_clusters",
                     umap.name = "umap.unintegrated")


### Examine and visualize PCA results a few different ways
print(obj[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(obj, dims = 1:2, reduction = "pca")

DimPlot(obj, label = F, reduction = "pca")

## `DimHeatmap()` allows for easy exploration of the primary sources of heterogeneity,
# and can be useful when trying to decide which PCs to include for further analyses.
# Both cells and features are ordered according to their PCA scores.
DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(obj, dims = 1:15, cells = 500, balanced = TRUE)

## `ElbowPlot`: a ranking of priciple components based on 
# the percentage of variance explained by each one
ElbowPlot(obj, ndims = 20)

## Things to consider when choosing number of PCs
# - Rare cell types may be associated with higher PCs, but they are difficult
#   to distinguish from background noise for a small dataset without prior knowledge.
# - Encourage to repeat downstream analyses with a different number of PCs (10, 15, etc.)
# - Advise to err on the higher side when choosing this parameter
# The results are similar with more PCs used.


### Visualize UMAP results
p1 <- DimPlot(obj, 
              reduction = "umap.unintegrated", 
              label = FALSE,
              #group.by = c("Method", "predicted.celltype.l2", 
              #             "unintegrated_clusters") # for pbmcsca dataset
              #group.by = c("stim", "seurat_clusters") # for ifnb dataset
              group.by = c("unintegrated_clusters")
              ) + ggtitle("Unintegrated Clusters")
p1


### Perform integrative analysis ###

### Supported integration method methods (v5)
# - CCAIntegration
#   - Canonical Correlation Analysis
#   - Classic Seurat v3/v4 approach (FindIntegrationAnchors + IntegrateData)
#   - Suitable for smaller dataset
#
# - RPCAIntegration
#   - Reciprocal PCA
#   - Faster and more scalable than CCA. Best if share many common cell states
#
# - HarmonyIntegration
#   - Uses the Harmony algorithm
#   - Iteratively adjusts PCs to remove batch effects while preserving biological variation
#
# - FastMNNIntegration
#   - Uses mutual nearest neighbors (MNN) batch-correction algorithm
#   - Can over-correct if biology differs strongly across batches (i.e. different tissue/lines)
#
# - JointPCAIntegration
#   - Performs PCA on the combined datasets simultaneuosly (instead of finding anchors)
#   - Very lightweight compared to other methods, good for well-matched datasets
#   - Weaker batch correction and not suitable for strong batch effect
#
# - SCVIIntegration
#   - scVI (single-cell Variational Inference)
#   - Uses scvi-tools Python library (so requires reticulate + scvi-tools installed)
#   - Accounts for both biological variation and batch effects
#   - Often outperforms traditional methods on large or complex datasets, but slow

## Benchmarking integration methods:
# https://www.nature.com/articles/s41592-021-01336-8

## Specify integration method (options: CCA, RPCA, Harmony, MNN, scVI)
# Enter any value other than the listed options to skip integration (i.e., use PCA)
specify_int_method <- "CCA"
int_method <- tolower(specify_int_method)

if (int_method %in% c("cca", "rpca", "mnn","harmony", "scvi")){
  
  # Map the specified method to its corresponding integration function
  integration_fun <- list(
    cca     = CCAIntegration,
    rpca    = RPCAIntegration,
    mnn     = FastMNNIntegration,
    harmony = HarmonyIntegration,
    scvi    = scVIIntegration
  )[[int_method]]
  
  # Set standardized names for the integrated reduction, clusters, and UMAP
  integration.use <-  paste0("integrated.", int_method)
  int.cluster.name <- paste0(int_method, "_clusters")
  int.umap.name <- paste0("umap.", int_method)
  if (int_method == "harmony"){
    # Using the standard naming convention for other methods causes Harmony to issue a warning.
    # This warning does not affect the results, but setting the reduction name to "harmony" (default)
    # prevents the message and keeps the output consistent.
    integration.use = "harmony"
  }
  
  # Run integration
  if (int_method == "mnn"){
    # FastMNNIntegration can't take additional arguments
    # https://github.com/satijalab/seurat/issues/8631
    
    obj <- IntegrateLayers(
      object = obj,
      method = integration_fun,
      new.reduction = integration.use,
      verbose = TRUE)
    
  } else if(int_method == "scvi"){
    # scVI needs to specify conda environment
    
    obj <- IntegrateLayers(
      object = obj,
      method = integration_fun,
      #orig.reduction = "pca", 
      new.reduction = integration.use,
      #normalization.method = specify_norm_method,
      conda_env = file.path(miniconda_path(), "envs", "r-scvi"),
      max_epochs = 400,
      verbose = TRUE)
    
  } else{
    # For CCA, RPCA, and Harmony
    
    obj <- IntegrateLayers(
      object = obj,
      method = integration_fun,
      orig.reduction = "pca", 
      new.reduction = integration.use,
      normalization.method = specify_norm_method,
      verbose = TRUE)
    
  }
  
  
  # Run clustering and umap with the new integration
  obj <- DimRed.seurat(obj, projection.method = "nonlinear",
                       reduction.use = integration.use,
                       dims.use = 1:30, resolution = 2,
                       cluster.name = int.cluster.name,
                       umap.name = int.umap.name)
  
  # Visualize
  p2 <- DimPlot(obj, 
                reduction = int.umap.name, 
                label = FALSE,
                #group.by = c("Method", "predicted.celltype.l2", int.cluster.name) # for pbmcsca dataset
                #group.by = c("stim", "seurat_clusters") # for ifnb dataset
                group.by = c(int.cluster.name)
                ) + ggtitle(paste(specify_int_method, "Clusters"))
  
} else{
  warning("No valid integration method specified. Proceeding without integration.")
}

wrap_plots(c(p1, p2), nrow = 2)



### Save the object for easy reload later
# Save a checkpoint so previous steps do not need to be rerun

#saveRDS(obj, file = "../output/seurat_tutorial.rds")


### Annotate clusters using canonical markers

## Some canonical markers for immune cells
# CD8 T cell populations (naive, memory, effector): CD8A, GZMK, CCL5, CCR7
# CD4 T cell populations (naive, memory, IFN-activated): S100A4, CCR7, IL32, ISG15
# B cell (developmental sub-structure): TCL1A, FCER2
# NK cells (CD56dim vs. bright clusters): XCL1 and FCGR3A

## Visualize canonical marker genes as violin plots.
VlnPlot(obj, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
        pt.size = 0.2, ncol = 4, group.by = "unintegrated_clusters")

VlnPlot(obj, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7", "ISG15", "CD3D"),
        pt.size = 0.2, ncol = 4, group.by = int.cluster.name)

# Compare specific gene
p3 <- VlnPlot(obj, features = "CD8A", group.by = "unintegrated_clusters") + 
  NoLegend() + ggtitle("Unintegrated Clusters")

p4 <- VlnPlot(obj, features = "S100A4", group.by = int.cluster.name) + 
  NoLegend() + ggtitle(paste(specify_int_method, "Clusters"))

p3 | p4

p2 | p4

## Visualize canonical marker genes on the sctransform embedding.
p5 <- FeaturePlot(obj, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), 
                  pt.size = 0.2, ncol = 3, reduction = "umap.unintegrated")

p6 <- FeaturePlot(obj, features = c("CD8A", "GZMK", "CCL5", "S100A4", "ANXA1", "CCR7"), 
                  pt.size = 0.2, ncol = 3, reduction = int.umap.name)

p2 | p6


### Rejoin the layers after integrative analysis is complete
# This collapses the individual datasets together and recreates the original
# `counts` and `data` layers.
# This is required before performing any differential expression analysis
obj <- JoinLayers(obj)
obj



### ============================================================================
### 5. Differentially Expression Analysis

## Additional sections to be added
