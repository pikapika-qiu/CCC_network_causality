# Install libraries if not avaliable
# install.packages("Seurat")
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)

# Convert the Seurat object to an AnnData file 
SaveH5Seurat(my_Seuratobj, filename = "my_Seuratobj.h5Seurat")
Convert("my_Seuratobj.h5Seurat", dest = "h5ad")
file.remove("my_Seuratobj.h5Seurat") # Remove the intermediate file to free space
