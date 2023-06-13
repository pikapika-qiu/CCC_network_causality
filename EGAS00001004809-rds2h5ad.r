##
# RDS to H5ad file conversion

# Path: rds2h5ad.r

library(HDF5Array)
library(anndata)
library(Seurat)



rdsfiles  = c('1863-counts_cells_cohort1.rds',   '1864-counts_tcell_cohort1.rds',  '1865-counts_myeloid_cohort1.rds',  '1866-counts_DC_cohort1.rds',  '1867-counts_cells_cohort2.rds')
for (i in 1:seq_along(rdsfiles)) {
    rds_in = paste("/data/ICI_exprs/EGAS00001004809/", rdsfiles[i], sep="")
    print(paste("Reading RDS file:", rds_in))
    rds.obj <- readRDS(rds_in)
    #showClass(class(rds.obj))

    X <- t(as(rds.obj, "Matrix"))

    adata <- AnnData(X = X)
    adata$obs_names <- rds.obj@Dimnames[[2]]
    adata$var_names <- rds.obj@Dimnames[[1]]


    # Write the AnnData object to disk as h5ad file
    h5ad_out <- gsub(x = rds_in, pattern = ".rds", replacement = ".h5ad")
    print(paste("Writing h5ad file:", h5ad_out))
    write_h5ad(adata, h5ad_out)
    print(paste("Done with", rds_in))
}