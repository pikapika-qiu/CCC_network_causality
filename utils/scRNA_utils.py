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

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

def plotGEMs(adata, cluster_id_col, cluster_id_to_plot, ncols=4):
    '''
    Plot the expression of genes in a cluster
    Parameters:
        adata: AnnData object
        cluster_id_col: column name in adata.obs that contains the cluster id
        cluster_id_to_plot: A list of cluster ids to plot
        ncols: number of columns in the plot
    '''

    # check input values
    if cluster_id_col not in adata.obs.columns:
        print("Error: cluster_id_col not found in adata.obs.columns")
        return
    elif not cluster_id_col:
        print("Error: cluster_id_col is empty")
        return
    
    # check if cluster_id_to_plot is in adata.obs[cluster_id_col]
    if cluster_id_to_plot not in adata.obs[cluster_id_col].unique():
        print("Error: cluster_id_to_plot not found in adata.obs[cluster_id_col]")
        return
    
    # identify the cells assigned in the cluster_id_to_plot
    adata_tmp = adata[adata.obs[cluster_id_col] == cluster_id_to_plot, :].copy()  
    nCells, nGenes = adata_tmp.shape
    # search for GEMs expressed in the cells of this cluster
    GEMs_exprs_in_cls = adata_tmp.var_names[(np.sum(adata_tmp.X > 25, axis= 0) / nCells > .05)].tolist()
    GEMs_exprs_in_cls = ['timepoint'] + GEMs_exprs_in_cls
    sc.pl.umap(adata_tmp, color = GEMs_exprs_in_cls )

    
def parseNHDP_RData(file_name):
    '''
    This function parses the RData file generated by NHDP R package
    and returns two adata objects: cell_by_GEM and GEM_by_gene

    Input: file_name: the name of the RData file returned by Han's NHDP R package
            The content of the RData file is a ListVector with the following elements:
            "nHDP_init", "nHDP_trained_mb", "X", which are initial tree, model trained with minibatch, 
            and the count matrix, respectively.
    Output: a dictionary with two keys: 'cell_by_GEM' and 'GEM_by_gene', in the format of adata object
    '''
    
    # Load the RData file
    robjects.r['load'](file_name)
    nHDP_trained_mb = robjects.r['nHDP_trained_mb']

    # extract $cell from nHDP_trained_mb, which is a ListVector
    cell_name_vector = nHDP_trained_mb.rx('cell')
    cell_names = list(cell_name_vector[0])  

    # extract $count_matrix from nHDP_trained_mb, which is a Matrix    
    count_matrix_obj = nHDP_trained_mb.rx('count_matrix')
    count_matrix_mat = np.array(count_matrix_obj[0]).T

    # extract $gene from nHDP_trained_mb, which is a ListVector
    gene_obj = nHDP_trained_mb.rx('gene')
    genes = list(gene_obj[0])

    # extract $centroids from nHDP_trained_mb, which is a matrix
    centroids_obj = nHDP_trained_mb.rx('centroids')
    centroids_mat = np.array(centroids_obj[0])
    
    # create two adata objects
    cell_by_GEM_adata = ad.AnnData(X=count_matrix_mat, obs=pd.DataFrame(index=cell_names), var=pd.DataFrame(index=range(count_matrix_mat.shape[1])))
    GEM_by_gene_adata = ad.AnnData(X=centroids_mat, obs=pd.DataFrame(index=genes), var=pd.DataFrame(index=range(centroids_mat.shape[1])))
    return {'cell_by_GEM': cell_by_GEM_adata, 'GEM_by_gene': GEM_by_gene_adata}

def clustering_adata(adata, resolution = 0.5, n_top_genes=5000):
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
    
    # check if adata has already been through selection of high variance genes
    if not 'highly_variable' in adata.var.columns:
        print ("Select ", n_top_genes, " high variance genes")
        # select high veriable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
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
    # sc.pl.umap(adata, color=['leiden'],  title='leiden')

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

def analyze_cell_type(adata, cell_type, markers, adata_name):
    '''
    this function should automatically extract a desired cell type for the user to save to a .h5ad file.  

    cell_type: the type of cell you wish to extract ("T cell" or "Myeloid", etc)
    markers: a dictionary of cell markers for the desired cell type, provided by the user 
    adata_name: as this makes a copy of adata, provide what you would like the copy to be named (ex: "adata_T")
    '''

    # set Scanpy plotting parameters
    sc.set_figure_params()

    # make a copy of adata
    adata = adata.copy()

    # extract cells and create a new AnnData object
    adata_type = adata[adata.obs['cell_type'] == cell_type].copy()

    # restore the X to the original raw.X for re-processing
    adata_type = ad.AnnData(X=adata_type.raw.X, obs=adata_type.obs, var=adata_type.raw.var, obsm=adata_type.obsm, uns=adata_type.uns)
    adata_type.raw = adata_type
    print(str(adata_type.shape))

    # re-calculate highly variable genes
    if 'highly_variable' in adata_type.var.columns:
        adata_type.var.drop(columns="highly_variable", inplace=True)          
    sc.pp.highly_variable_genes(adata_type, n_top_genes=5000)
    sc.pl.highly_variable_genes(adata_type)
    len(adata_type.var_names)

    # re-cluster the specified cell type
    sc.tl.pca(adata_type, svd_solver='arpack', n_comps=40)
    sc.pp.neighbors(adata_type, n_neighbors=80, n_pcs=40)
    sc.tl.leiden(adata_type, resolution=.5)
    
    sc.tl.umap(adata_type)
    sc.pl.umap(adata_type, color='leiden', legend_loc='on data')    

    # apply cell type labels using the marker dictionary
    adata_type.obs.drop(columns="cell_type", inplace=True)
    labelClusterWithCellType(adata_type, markers, cluster_column='leiden')

    # UMAP
    sc.pl.umap(adata_type, color='cell_type')

    # more UMAPs
    sc.pl.umap(adata_type, color=['timepoint', 'cell_type'])

    # save a copy of adata_type under custom name
    globals()[adata_name] = adata_type


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
    if not sample_id_col and not 'sample_id' in adata.obs.columns:
        print ("sample id column not provided")
        return None
    
    if not sample_id_col and 'sample_id' in adata.obs.columns:
        sample_id_col = 'sample_id'
    
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
    if 'pseudoBulk' not in adata.uns.keys():
        print ("Input adata is not pseudo-bulk RNA data. Convert to pseudo-bulk RNA data.")
        adata = scRNA2PseudoBulkAnnData(adata, sample_id_col=sample_id_col)
    
    # Create a 3-d matrix, one dimension is the patient, the other is the gene, the third is the condition
    nPatients = len(adata.obs[patient_id_col].unique())
    nGenes = len(adata.var_names)
    nConditions = len(adata.obs[condition_key].unique())
    if nConditions != 2:
        print ("Number of conditions is not 2")
        return None
    
    X = np.zeros((nConditions, nPatients, nGenes), dtype=np.float32)

    condition1 = adata.obs[condition_key].unique()[0]
    condition2 = adata.obs[condition_key].unique()[1]
    condition1_mean_name = condition1 + '_mean'
    condition2_mean_name = condition2 + '_mean'

    # create a dataframe to store the results
    colNames = ['pval', 'log2fc', condition1_mean_name, condition2_mean_name]
    res_df = pd.DataFrame(index=adata.var_names, columns = colNames)
    patients = adata.obs[patient_id_col].unique()  # this is a numpy array
    
    for index, patient in np.ndenumerate(patients):
        indx_p = index[0]
        # print ("Processing patient %s" % patient)
        # check if the patient has two conditions
        if len(adata.obs[condition_key][adata.obs[patient_id_col] == patient].unique()) < 2:
            # print ("Patient %s does not have two conditions" % patient)
            continue
        # extract data from the patient under condition 1 and condition 2

        # print ("Extract data from patient %s under condition %s & %s" % (patient, condition1, condition2))
        X[0, indx_p, :] = adata.raw.X[(adata.obs[patient_id_col] == patient) & (adata.obs[condition_key] == condition1), :]
        X[1, indx_p, :] = adata.raw.X[(adata.obs[patient_id_col] == patient) & (adata.obs[condition_key] == condition2), :]
        
    # perform paired t-test 
    # for each gene, perform t-test between two conditions of the same patient
    for i in range(nGenes):  # need check how to parallelize this loop, maybe use cupy
        x_1 = X[0, :, i]
        x_2 = X[1, :, i]
        
        # check if x_1 and x_2 are all zeros
        if np.sum(x_1) == 0 or np.sum(x_2) == 0:
            continue
        
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
        res_df.loc[gene_name, condition1_mean_name] = mean_condition1
        res_df.loc[gene_name, condition2_mean_name] = mean_condition2

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
        pre_data = gene_data[gene_data.obs[condition_key] == 'Pre']
        on_data = gene_data[gene_data.obs[condition_key] == 'On']

        # perform t-test using scipy
        t_stat, p_value = stats.ttest_ind(pre_data.X, on_data.X)

        # store statistics in a dict
        gene_stats[gene] = {'t_stat': t_stat, 'p_value': p_value}

        # check if differentially expressed
        if np.abs(t_stat) > 0:
            DEGs.append(gene)

    return DEGs

def findDEGs(adata, cluster_id, condition_col, cluster_id_col = 'leiden', method = 'wilcoxon'):
    pass
    # '''
    #     We want to find genes that are differentially expressed between cells from a cluster but collected
    #         under different conditions from a study. Assume cluster_id col is contained in the leiden col.
            
    #         The input is adata, which is a scRNA-seq data with the following structure:
    #         adata.obs:
    #         - sample_id: sample id
    #         - cluster_id: cluster label for cell to be examined.  For example cells belonging to a specific T cell cluster
    #         _ condition: The conditions under which cells are collected and to be compared. 

    #         Parameters:
    #         adata
    #         cluster_id: a scalar(leiden: string) indicating the interested cluster labels
    #         condition_col: the column name in adata.obs that contains the conditions information
    #         cluster_id_col: the column is contianed in the leiden column(cluster ids)
    #         method: 'Wilcoxon' or 't-test'
            
    #         Return:
    #         a list of DEGs
            
    # '''
    # # Steps/pseudocode:
    # # 1. Filter cells based on the cluster --> adata object containing cells from the two clusters
    # # subset = cluster_id

    # adata_subset = adata[adata.obs[cluster_id_col]==cluster_id,:].copy()
    
    # # determine the rank of each gene in a given cell.  That is, we need find a function in python which determines rank of of a gene
    # #  in a row and put the rank-order value in place.  # 2. call rank_genes_groups with the two clusters as groups with method = 'wilcoxon'
    # # condition/treatment col
    # gene_ranks = np.argsort(adata_subset.X, axis=1) # need to check if this is correct

    # # loop through genes
    # # for each gene, find the rank of the gene in each cell
    # # find the cells (rows) that belong a specific condition
    # # collect rank values of the gene in the cells that belong to a specific condition
    # # perform wilcoxon test between the two sets of rank values
    # for i in range(adata_subset.shape[1]):

    # unique_conditions = adata_subset.obs[condition_col].unique()
    # sc.tl.rank_genes_groups(adata_subset, groupby= condition_col , groups=unique_conditions, method=method)
    
    # # 3. extract the list of DEGs
    # result = adata_subset.uns['rank_genes_groups']
    # de_genes = result['names'][condition_col][cluster_2]
    
    # return de_genes.tolist()
    
# def performDEG(adata, groupby, group1 = 'pre', group2 = 'on'):
#     """
#     Perform differential expression analysis between two clusters in an AnnData object

#     Parameters:
#         adata (scanpy.AnnData)
#         groupby (str): Key in `adata.obs` grouping
#         group1 (str): Name of the first group to contrast, default 'pre'
#         group2 (str): Name of the second group to contrast, default 'on'

#     Returns:
#         pandas.DataFrame
#     """
#     # Check if 'louvain' column exists in adata.obs
#     if 'louvain' not in adata.obs.columns:
#         print("Column 'louvain' not found in adata.obs. Run clustering first.")
        
        
#     # Identify clusters
#     sc.tl.louvain(adata)
#     clusters = adata.obs['louvain'].unique()

#     for cluster_id in clusters:
#         # Extract the data only for the current cluster
#         sub_adata = adata[adata.obs['louvain'] == cluster_id]

#         # Perform differential expression analysis between two groups
#         sc.tl.rank_sum_test(sub_adata, groupby=groupby, groups=[group1, group2])

#         # Access the results
#         result = sub_adata.uns['rank_sum_test']
#         if 'gene_names' in result and 'reject' in result:
#             cluster_diff_genes = result['gene_names'][result['reject']]  # Differentially expressed genes for this cluster
#             cluster_diff_genes = pd.DataFrame(cluster_diff_genes, columns=[f'Cluster_{cluster_id}_DEGs'])
#             groups_df = pd.concat([groups_df, cluster_diff_genes], axis=1)
            
#     # for cluster_id in adata.obs['']
#     # # Perform differential expression analysis between two clusters
#     # sc.tl.rank_sum_test(adata, groupby=groupby, groups=[group1, group2])

#     # # Access the results
#     # result = adata.uns['rank_sum_test']
#     # groups_df = result['gene_names'][result['reject']]  # Differentially expressed genes
    
#     return groups_df

# # Example usage:
# # Assuming you have already loaded and preprocessed your AnnData object as 'adata'
# # diff_genes = performDEG(adata, groupby='louvain', group1='cluster_1', group2='cluster_2')


def findDEGsFromClusters(adata, condition_col=None, condition_1=None, condition_2=None, top_n_degs=100):
    '''
    This function searches for clusters and then finds DEGs with each cluster conditioning on specified conditions.

    Parameters
    --------
    adata: AnnData object
        Annotated data matrix with rows for cells and columns for genes.
    condition_col: the column name of the condition in the adata.obs
    condition_1: the condition_1    
    condition_2: the condition_2
    top_n_degs: Number of top DEGs to consider for plotting

    Returns:
    --------
    DEGs: A dataframe with DEGs and their logFC, pval, pval_adj, etc.
    significant_genes_dict: A dictionary containing significant genes for each cluster.

    pseudocode:
    1. find clusters by call leiden or louvian by clustering_adata function
    2. loop through each cluster:
        2.1. extract cells belonging to the cluster (adata.copy())
        2.2. Call paird_ttest funciton using the adata_cluster find DEGs conditioning on the condition_1 and condition_2
        2.3. return the dataframe of DEGs
        2.4. find significant genes using sc.tl.rank_genes_groups
        2.5. plot UMAP for significant genes
        2.6. save a dataframe of significant DEGs
        2.7. plot volcano plot of DEGs

    '''

    # 1: find clusters using leiden or louvain by calling clustering_adata function
    if condition_col is None or condition_1 is None or condition_2 is None:
        print("Error: Missing condition information.")
        return None

    adata_clusters = clustering_adata(adata) # call clustering_adata function

    # 2: loop through each cluster, extract cells belonging to the cluster, and find DEGs
    clusters = adata_clusters.obs['leiden'].unique()
    result_dfs = []  # store DEG dataframes for each cluster
    significant_genes_df = {} # store dataframe for significant DEGs

    for cluster in clusters:
        print(f"Finding DEGs for cluster {cluster}")

        # 2.1. extrac cells belonging to the cluster (adata.copy())
        adata_cluster = adata_clusters[adata_clusters.obs['leiden'] == cluster].copy()

        # 2.2. Call paired_ttest function using the adata_cluster to find DEGs conditioning on condition_1 and condition_2
        DEGs_cluster = paird_ttest(adata_cluster, condition_key=condition_col, sample_id_col='sample_id', patient_id_col='patient_id', pval_cutoff=0.05, log2fc_cutoff=1)

        # 2.3. return the dataframe of DEGs
        if DEGs_cluster is not None:
            result_dfs.append(DEGs_cluster)
        
        # Create a copy of the original adata and apply log1p transformation to the adata_cluster
        adata_copy = adata.copy()
        sc.pp.log1p(adata_copy[adata_copy.obs.index.isin(adata_cluster.obs.index)])
                
        # 2.4. find significant genes using sc.tl.rank_genes_groups
        print(f"plotting significant genes for cluster {cluster}")

        # Set the 'base' value in adata.uns['log1p']
        adata_cluster.uns['log1p'] = {'base': 2} # check if already logged
        
        sc.tl.rank_genes_groups(adata_cluster, groupby=condition_col, method='wilcoxon')
        sc.pl.rank_genes_groups(adata_cluster, n_genes=25, sharey=False)


        print(f"DEGs: \n{DEGs_cluster}")

        # 2.5 some UMAPs
        sc.pp.neighbors(adata_cluster, n_neighbors=30, n_pcs=50)
        sc.tl.umap(adata_cluster)
        sc.pl.umap(adata_cluster, color=['cell_type', 'timepoint'], legend_loc='on data', title = f"cluster {cluster}")
        
        # UMAP for DEGs
        if not DEGs_cluster.empty:
            # Convert 'pval' column to numeric type
            DEGs_cluster['pval'] = pd.to_numeric(DEGs_cluster['pval'])
            
            top_n_degs_cluster = DEGs_cluster.nsmallest(top_n_degs, 'pval')
            sc.pl.umap(adata_cluster, color=top_n_degs_cluster.index.tolist(), use_raw=False, cmap='viridis', legend_loc='on data')


        # 2.6. save a dataframe of significant DEGs
        significant_genes_df = pd.DataFrame(columns=['pval', 'log2FC', 'mean_1', 'mean_2', 'qval'])
        
        if not DEGs_cluster.empty:
            significant_genes = DEGs_cluster[(DEGs_cluster['pval'] < 0.05) & (DEGs_cluster['qval'] < 0.1)]
            significant_genes_df = significant_genes_df.append(significant_genes, ignore_index=True)
            # seperate clusters and save as csv
            
        # 2.7 volcano plot: (still a WIP)

        # if not DEGs_cluster.empty:
        #     print('creating volcano plot:')
        #     # Create a volcano plot using hvplot
        #     vc_df = significant_genes_df
        #     vc_df['names'] = vc_df.index
        #     vc_df.hvplot.scatter(
        #         "log2FC", "pval", 
        #         flip_yaxis=True, logy=True, 
        #         hover_cols=["names"]
        #     ).opts(title='Volcano Plot of DEGs')

        # # Set the significance thresholds
        # pval_threshold = 0.05
        # qval_threshold = 0.1

        # # Select significant DEGs
        # significant_DEGs = significant_genes_df[(significant_genes_df['pval'] < pval_threshold) & (significant_genes_df['qval'] < qval_threshold)]

        # # Create a volcano plot
        # plt.figure(figsize=(10, 6))
        # plt.scatter(significant_DEGs['log2FC'], -1 * np.log10(significant_DEGs['pval']), c='blue', label='Significant DEGs')

        # # Add vertical lines for fold change cutoffs
        # plt.axvline(x=1, color='red', linestyle='--', label='Log2FC = 1')
        # plt.axvline(x=-1, color='red', linestyle='--', label='Log2FC = -1')

        # # Add labels and legend
        # plt.xlabel('log2FC')
        # plt.ylabel('-log10(p-value)')
        # plt.title('Volcano Plot of DEGs')
        # plt.legend()
        # plt.grid(True)
        # plt.show()


    # Combine all the DEG dataframes into a single DataFrame
    DEGs = pd.concat(result_dfs)

    return DEGs, significant_genes_df