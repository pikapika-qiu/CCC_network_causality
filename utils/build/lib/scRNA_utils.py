''' 
scRNA_utils.py

This file contains utility functions for scRNA-seq analysis commonly used in the lab.
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

def clustering_adata(adata, resolution = 0.5):
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
    sc.tl.leiden(adata, resolution = resolution)

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
            tmp.obs[sample_id_col] = prefix 
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

def scRNA2PseudoBulkAnnData(adata, sample_id_col = None): 
    '''        
        This function convert a scRNA AnnData oboject to an AnnData object,
           where gene expression from the same sample is merged and normalized as 
           transcript per million (TPM) format.  
         
        Parameters:
            adata: anndata object
            sample_id_col: the column in adata.obs that contains the sample id
        
        Returns:
            adata: AnnData object with adata.X in TPM format.  The annData object 
            is annoted with uns["pseudoBulk"] = "log_2_tpm"
        
    '''
    # check if input adata is AnnData object
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
    
    # check if adata have raw data
    if not adata.raw:
        print ("adata.raw is not available")
        return None

    col_to_remove = ['ncount_rna', 'nfeature_rna', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'n_genes_by_counts', 'log1p_n_genes_by_counts']
    col_to_keep_in_obs = [x for x in adata.obs.columns.str.lower() if x not in col_to_remove]

    nSamples = len(adata.obs['sample_id'].unique()) 
    nGenes = len(adata.var_names)
    X = np.zeros((nSamples, nGenes), dtype=np.float32)
    df_tpm = pd.DataFrame(X, index=adata.obs['sample_id'].unique(), columns = adata.var_names)

    # remove obs columns that are added by sc.pp functions
    col_to_remove = ['ncount_rna', 'nfeature_rna', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt', 'n_genes_by_counts', 'log1p_n_genes_by_counts']
    col_to_keep_in_obs = [x for x in adata.obs.columns.str.lower() if x not in col_to_remove]
    df_obs = pd.DataFrame(index=adata.obs['sample_id'].unique(), columns = col_to_keep_in_obs)

    for sample in adata.obs['sample_id'].unique():
        tpm = np.sum(adata.X[adata.obs['sample_id'] == sample, :], axis = 0)
        tpm = np.array(tpm / np.sum(tpm) * 1e6, dtype=np.float32) # normalize to TPM/per cell and force to float32
        df_tpm.loc[sample,:] = tpm

        # Populate df_obs
        for col in adata.obs.columns:
            df_obs.loc[sample, col] = adata.obs.loc[adata.obs[sample_id_col] == sample, col].unique()[0]
 

    # Create an AnnData object for the pseudo-bulk RNA data
    adata_sample_tpm = ad.AnnData(df_tpm.values, obs=df_obs, var=adata.var)
    adata_sample_tpm.uns["pseudoBulk"] = "tpm"
    adata_sample_tpm.raw = adata_sample_tpm

    return adata_sample_tpm

def paird_ttest(adata, condition_key = None, sample_id_col = None, patient_id_col = None, pval_cutoff = 0.05, log2fc_cutoff = 1):
    '''
    This function is to find the genes or gene modules that are differentially expression
    between two conditions collected from a same subject, e.g., tumor-vs-normal or before or after a 
    specific treatment. The function will perform pairwise t-test between two conditions for each gene.

    Steps in the process:
        1. Create pseudo-bulk RNA data for each sample 
        2. Identify cells from a sample that belong to a specific sample.
        3. Match samples from the same patient.
        4. Perform pairwise t-test between two conditions for each gene.


    Parameters:
        adata: AnnData object with adata.X in TPM format.  The annData object
            If annoted with uns["pseudoBulk"] = "log_2_tpm", the data is pseudo-bulk RNA in log2(TPM+1) format.
        Condition_key: the column in adata.obs that contains the condition information based on which pairwise t-test will be performed.
        sample_id_col: the column in adata.obs that contains the sample id
        patient_id_col: the column in adata.obs that contains the patient id
    
    return:
        A dataframe consisting of a list of genes and statistics of pair-wise t-test between two conditions.
    
    '''

    # check inputs
    if not isinstance(adata, ad.AnnData):
        print ("Input adata is not an AnnData object")
        return None 
    if not condition_key:
        print ("Condition key not provided")
        return None
    # check if condition to compare is binary
    if len(adata.obs[condition_key].unique()) != 2:
        print ("Condition to compare is not binary")
        return None
    if not sample_id_col:
        print ("sample id column not provided")
        return None
    if not patient_id_col:
        print ("patient id column not provided")
        return None
    # check if adata have raw data
    if not adata.raw:
        print ("adata.raw is not available")
        return None
    
    # assume data is already pseudo bulk, check
    if not adata.uns["pseudoBulk"] :
        print ("Input adata is not pseudo-bulk RNA data. Convert to pseudo-bulk RNA data.")
        adata = scRNA2PseudoBulkAnnData(adata, sample_id_col=sample_id_col)
    
    # Create a 3-d matrix, one dimension is the patient, the other is the gene, the third is the condition
    nPatients = len(adata.obs[patient_id_col].unique())
    nGenes = len(adata.var_names)
    nConditions = len(adata.obs[condition_key].unique())
    X = np.zeros((nConditions, nPatients, nGenes), dtype=np.float32)

    res_df = pd.DataFrame(index=adata.var_names, columns = ['pval', 'log2fc', 'mean_condition1', 'mean_condition2'])
    patients = adata.obs[patient_id_col].unique()  # this is a numpy array
    for index, patient in np.ndenumerate(patients):
        indx_p = index[0]
        # print ("Processing patient %s" % patient)
        # check if the patient has two conditions
        if len(adata.obs[condition_key][adata.obs[patient_id_col] == patient].unique()) < 2:
            # print ("Patient %s does not have two conditions" % patient)
            continue
        # extract data from the patient under condition 1 and condition 2
        condition1 = adata.obs[condition_key].unique()[0]
        condition2 = adata.obs[condition_key].unique()[1]
        # print ("Extract data from patient %s under condition %s & %s" % (patient, condition1, condition2))
        X[0, indx_p, :] = adata.raw.X[(adata.obs[patient_id_col] == patient) & (adata.obs[condition_key] == condition1), :]
        X[1, indx_p, :] = adata.raw.X[(adata.obs[patient_id_col] == patient) & (adata.obs[condition_key] == condition2), :]
        
    # perform paired t-test 
    # for each gene, perform t-test between two conditions of the same patient
    for i in range(nGenes):  # need check how to parallelize this loop, maybe use cupy
        x_1 = X[0, :, i]
        x_2 = X[1, :, i]
        pval = stats.ttest_rel(x_1, x_2)[1]
        gene_name = adata.var_names[i]        
        mean_condition1 = np.mean(x_1)
        mean_condition2 = np.mean(x_2)
        if mean_condition1 == 0 or mean_condition2 == 0:
            log2fc = np.nan
        else:
            log2fc = np.log2(np.mean(x_1) / np.mean(x_2))
        res_df.loc[gene_name, 'pval'] = pval
        res_df.loc[gene_name, 'log2fc'] = log2fc
        res_df.loc[gene_name, 'mean_condition1'] = mean_condition1
        res_df.loc[gene_name, 'mean_condition2'] = mean_condition2

    # estimate q-value based on p-value        
    qvalue = importr('qvalue')
    r_p_values = robjects.FloatVector(res_df['pval'])
    r_q_values = qvalue.qvalue(r_p_values)
    res_df['qval'] = np.array(r_q_values.rx2('qvalues'))

    return res_df


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
    for sample in adata.obs[sample_id_col].unique():
        # Find cells that belong to the specific cluster in this sample
        # Produce pseudo-bulk RNA data
        sample_mask = adata_cluster.obs[sample_id_col] == sample
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

def findDEGs(adata, cluster_id, condition_col,  method = 'Wilcoxon'):
    pass
'''
   We want to find genes that are differentially expressed between cells from common cluster but collected
    under different conditions from a study. 
    
    The input is adata, which is a scRNA-seq data with the following structure:
    adata.obs:
      - sample_id: sample id
      - cluster_id: cluster label for cell to be examined.  For example cells belonging to a specific T cell cluster
      _ condition: The conditions under which cells are collected and to be compared. 

    Parameters:
       adata
       cluster_id: a list of two cluster labels
       condition_col: the column name in adata.obs that contains the condition information
       method: 'Wilcoxon' or 't-test'
    
    Return:
     a list of DEGs
    
'''
    # Steps:
    # 1. Filter cells based on the cluster --> adata object containing cells from the two clusters
    # 2. call rank_genes_groups with the two clusters as groups with method = 'wilcoxon'
    # 3. extract the list of DEGs