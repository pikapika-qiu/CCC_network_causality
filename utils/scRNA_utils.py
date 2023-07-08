''' 
load_10X_matrices.py

This function module will load a fold of 10X matrix files into a single sparse matrix
and return a single AnnData object by concatenating the matrices.
'''

import os as os
import sys as sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scRNA_utils import *
from scipy import stats

def clustering_adata(adata):
    '''
    This function will cluster an AnnData object 

    Parameters:
        adata: AnnData object

    Returns:
        adata: AnnData object with a new column in adata.obs called 'leiden' that contains the cluster label for each cell
    '''

    # check if adata is AnnData object
    if not isinstance(adata, ad.AnnData):
        print ("Input adata is not an AnnData object")
        return None
    
    # check if adata has raw data
    if adata.raw is None:
        print ("Input adata does not have raw data")
        return None
    
    # check if adata has more than 1000 cells
    if adata.shape[0] < 1000:
        print ("Input adata has less than 1000 cells")
        return None
    
    # check if adata has more than 1000 genes
    if adata.shape[1] < 1000:
        print ("Input adata has less than 1000 genes")
        return None 
    
    # check if adata has more than 10000 genes
    if adata.shape[1] > 10000:
        # select high veriable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        # filter adata
        adata = adata[:, adata.var['highly_variable']]

        # check if X is log transformed
        if not 'log1p' in adata.layers:
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata, base = 2)

    # run PCA
    sc.tl.pca(adata, svd_solver='arpack', n_comps=50)   
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    sc.tl.leiden(adata, resolution=0.5)

    #plot UMAP
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=['leiden'], legend_loc='on data', title='leiden')

    return adata



def load_10X_matrices(matrix_dir):
    '''
    load_10X_matrices(matrix_dir)
    
    This function will load a fold of 10X matrix files into a single sparse matrix.

    Parameters:
        matrix_dir (str): The directory containing the 10X matrix files.
    
    Returns:
        A single AnnData object containing the concatenated matrices.

    '''

    #check if matrix_dir is a directory
    if not os.path.isdir(matrix_dir):
        print ("Input " + matrix_dir + " is not a directory")
        return None
    
    # Open the matrix directory and get the list of files
    matrix_files = os.listdir(matrix_dir)

    # Loop through the files and concatenate the matrices
    mtx_files = [x for x in matrix_files if '.mtx' in x]
    if len(mtx_files) == 0:
        print ("Input directory " + matrix_dir + " has no mtx files.\n")
        return None
        
    # find out whether files has prefix, use set to remove duplicates
    prefixes = set([x.split('matrix.mtx')[0] for x in mtx_files if 'matrix.mtx' in x]) 

    if len(prefixes) > 0:
        # create a list to hold adatas
        adata_list = []   

        # interate through prefixes and index mtx_files        
        for index, prefix in enumerate(prefixes): 
            print("Loading " + prefix)
            tmp = sc.read_10x_mtx(matrix_dir, prefix = prefix, cache = True)
            tmp.obs['sample_ID'] = prefix 
            adata_list.append(tmp)

        # concatenate adata_list
        overall_adata = ad.concat(adata_list, join='outer')
        return overall_adata
    

def labelClusterWithCellType(adata, cell_type_markers, cluster_column='leiden'):
    '''
    This function will label each cluster with the cell type that is most abundant in that cluster.

    Parameters:
        adata: AnnData object
        cell_type_markers: a dictionary where the key is the cell type and the value is a list of markers for that cell type
        cluster_column: the column in adata.obs that contains the cluster labels

    Returns:
        adata: AnnData object with a new column in adata.obs called 'cell_type' that contains the cell type label for each cell
    
    '''

    # check if adata is AnnData object
    if not isinstance(adata, ad.AnnData):
        print ("Input adata is not an AnnData object")
        return None
    
    # check if cell_type_markers is a dictionary
    if not isinstance(cell_type_markers, dict):
        print ("Input cell_type_markers is not a dictionary")
        return None
    
    # Check if cell_type_markers is empty
    if len(cell_type_markers) == 0:
        print ("Input cell_type_markers is empty")
        return None
    
    # Check if 'cluster_column' is in adata.obs
    if cluster_column not in adata.obs.columns:
        print ("Input cluster_column " + cluster_column + " is not in adata.obs")
        return None
    

    # find total number of clusters   
    cls_ids = adata.obs[cluster_column].unique()
    
    # iterate through all cluster
    for i in cls_ids:
        # find cells in cluster i        
        cell_in_cls_i = adata.obs[cluster_column] == i  
        # this will return a vector of True/False where True means the cell is in cluster i
        # print('processing cluster: ' + str(i) + ' with ' + str(sum(cell_in_cls_i)) + ' cells')

        # keep track of which cell type is most abundant in cluster i
        cell_type_cluster_overlapp_pct = dict()

        #iterate through key and value of cell_type_markers
        for cell_type, marker_genes in cell_type_markers.items():   
            # Extract the expression of all marker genes for cells in cluster i
            # this will return a sparse matrix of cells x markers
            cell_w_marker_genes = adata.raw.X[:, adata.raw.var_names.isin(marker_genes)] > 0  
            
            # change cell_in_cls_i to numpy array and repeat it to match the shape of cell_w_marker_genes
            cell_in_cls_i_m = np.tile(cell_in_cls_i.to_numpy(), (cell_w_marker_genes.shape[1], 1)).T

            # find cells in cluster i that express the marker
            # this create a matrix of cells x markers where True means the cell express the marker and in cluster i
            cell_w_marker_genes = cell_w_marker_genes.toarray() & cell_in_cls_i_m

            # caclualte average markers expressed in each cell in Marker_genes_i
            nmarker_per_cell = np.sum(cell_w_marker_genes, axis=0) / cell_w_marker_genes.shape[1]
            #print(nmarker_per_cell)

            # keep track of which cell type is most abundant in cluster i
            # assuming the cell type with the highest average marker present is the most abundant         
            cell_type_cluster_overlapp_pct[cell_type] = np.sum(nmarker_per_cell) / sum(cell_in_cls_i)

        # check with cell type is most abundant in cluster i
        max_type = max(cell_type_cluster_overlapp_pct, key=cell_type_cluster_overlapp_pct.get)
        print('Cluster ' + str(i) + ' is most likely ' + max_type + ' with ' + str(cell_type_cluster_overlapp_pct[max_type]) + ' overlap')
        adata.obs.loc[cell_in_cls_i, 'cell_type'] = max_type    

    return adata

def find_pseudoBolk_RNA(adata, sample_id_col = None): 
    '''
        This function checks if adata and sample id column exists, then group cells from individual samples and produce tpm RNA to find pseudoBolk RNA
         
        Parameters:
            adata: anndata object
        
        Returns:
        
    '''
    if not isinstance(adata, ad.AnnData):
        print ("Input adata is not an AnnData object")
        return None
    if not sample_id_col:
        print ("sample id column not provided")
        return None
    # check if adata have sample id col
    if sample_id_col not in adata.obs.columns:
        print ("sample id", sample_id_col, "column not available in adata.obs")
        return None
    
    X = np.zeros((len(adata.obs['sample_id'].unique()), len(adata.var_names)))
    df_obs = pd.DataFrame(index=adata.obs['sample_id'].unique(), columns=['patient_id', 'timepoint', 'batch'])
    
    for sample in adata.obs['sample_id'].unique():
        tpm = np.sum(adata.X[adata.obs['sample_id'] == sample, :], axis=0)
        tpm = tpm / np.sum(tpm) * 1e6  # Normalize to TPM per cell(check for division by 0)
        df_tpm.loc[sample, :] = np.log2(tpm + 1) # check

        # Populate df_obs
        df_obs.loc[sample, 'patient_id'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'patient_id'].unique()[0]
        df_obs.loc[sample, 'timepoint'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'timepoint'].unique()[0]
        df_obs.loc[sample, 'batch'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'batch'].unique()[0]

    # Create an AnnData object for the pseudo-bulk RNA data
    adata_sample_tpm = ad.AnnData(X=X, obs=df_obs, var=adata.var)
    

    # Perform t-test
    cluster_data = adata_sample_tpm[adata_sample_tpm.obs[condition_key].isin(['pre', 'on'])]
    cluster_data.uns["pseudoBulk"] = "log_tpm"
    
    return adata_sample_tpm

def find_cluster_DEGs_pairwise(adata, cluster_label, condition_key):
    '''
    This function will find differentially expressed genes between two conditions for a given cluster.
    Steps in the process:
        1. Identify cells from a sample that belong to a specific cluster.
        2. Create pseudo-bulk RNA data for each sample.
        3. Match samples from the same patient.
        4. Perform pairwise t-test between two conditions for each gene.
    '''
    # assume data is already pseudo bulk, check
    # 
    
    # Filter cells based on the cluster
    cluster_mask = adata.obs['cluster'] == cluster_label
    adata_cluster = adata[cluster_mask].copy()
    # Create pseudo-bulk RNA data for each sample
    bulk_data = {}
    for sample in adata.obs['sample_id'].unique():
        # Find cells that belong to the specific cluster in this sample
        # Produce pseudo-bulk RNA data
        sample_mask = adata_cluster.obs['sample_id'] == sample
        bulk_data[sample] = np.array(adata_cluster.X[sample_mask].sum(axis=0)).flatten()

    # A dictionary to match samples from the same patient under two conditions.
    # Produce a matrix with the following axes: pre/on, N-patients, N-Genes.

    # create list for storing data
    DEGs = []

    # looping through var names
    for gene in cluster_data.var_names:
        gene_data = cluster_data[:, gene]

        # split data into two conditions
        pre_data = gene_data[gene_data.obs[condition_key] == 'pre']
        on_data = gene_data[gene_data.obs[condition_key] == 'on']

        # perform t-test using scipy
        t_stat, p_value = stats.ttest_ind(pre_data.X, on_data.X)

        # store statistics in a dict
        gene_stats[gene] = {'t_stat': t_stat, 'p_value': p_value}

        # check if differentially expressed
        if np.abs(t_stat) > 0:
            DEGs.append(gene)

    return DEGs



