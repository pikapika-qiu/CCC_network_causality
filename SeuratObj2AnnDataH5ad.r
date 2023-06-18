# This is an R script that reading a Seurat object and convert it to an AnnData file in h5ad format.

# Usage: Rscript SeuratObj2AnnDataH5ad.r my_Seuratobj

library(Seurat)
library(SeuratDisk)

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  Seuratobj_file <- args[1]
} else {
  stop("No Seurat object provided")
}

Seuratobj_file <- "/data/ICI_exprs/EGAS00001004809/SCTransform_tcells_cohort1.rds"
# Load the Seurat object from disk
my_Seuratobj <- readRDS(Seuratobj_file)
# Convert the Seurat object to an AnnData file 
SaveH5Seurat(my_Seuratobj, filename = "/data/ICI_exprs/EGAS00001004809/SCTransform_tcells_cohort1.h5Seurat")
Convert("/data/ICI_exprs/EGAS00001004809/SCTransform_tcells_cohort1.h5Seurat", dest = "h5ad")
file.remove("my_Seuratobj.h5Seurat") # Remove the intermediate file to free space
