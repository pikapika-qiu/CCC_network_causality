##
# RDS to H5ad file conversion

# Path: rds2h5ad.r

library(HDF5Array)
library(anndata)
library(Seurat)

# parse arguments
args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
        rds_in <- args[1]
    }
    else {
        stop("No RDS file provided")
        return(NULL)
    }
#rds_in <- "/data/ICI_exprs/EGAS00001004809/1863-counts_cells_cohort1.rds"
print(paste("Reading RDS file:", rds_in))
rds.obj <- readRDS(rds_in)
showClass(class(rds.obj))

X <- t(as(rds.obj, "Matrix"))

adata <- AnnData(X = X)
adata$obs_names <- rds.obj@Dimnames[[2]]
adata$var_names <- rds.obj@Dimnames[[1]]


# Write the AnnData object to disk as h5ad file
h5ad_out <- gsub(x = rds_in, pattern = ".rds", replacement = ".h5ad")
print(paste("Writing h5ad file:", h5ad_out))
write_h5ad(adata, h5ad_out)
print(paste("Done!"))