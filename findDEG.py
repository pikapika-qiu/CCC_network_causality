import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

def find_cluster_DEGs_pairwise(adata, cluster_label, condition_key):
    '''
    This function will find differentially expressed genes between two conditions for a given cluster
    Steps in the process:
        1. Identify cells from a sample belong to a specific cluster
        2. Create pseudo-bulk RNA data for each sample
        3. Match samples from same patient
        4. Perform pair-wise t-test between two conditions for each gene
    
    
    '''
    # filter cells based on cluster
    cluster_mask = adata.obs['cluster'] == cluster_label
    adata_cluster = adata[cluster_mask].copy()
    
    # create pseudo-bulk RNA data for each sample
    bulk_data = {}
    for sample in adata.obs['sample_id'].unique():
        # find cells belong to the specific cluster in this sample
        # produce pseudo-bulk RNA data
        sample_mask = adata_cluster.obs['sampe_id'] == sample
        bulk_data[sample] = np.array(adata_cluster.X[sample_mask].sum(axis=0)).flatten()
    
    # A dictionary match samples from same patient under two conditions
    # Produce a matrix with the following axis: pre/on, N-patients, N-Genes

    # Loop throught all genes
    # extract match pseudo-bulk RNA data for the gene in all patients
    # perform pairwised t-test between two conditions for the gene; Check scipy for pair-wise t-test
    # 
    # 
    #    
    # create data frames for 'pre' and 'on' matrices
    pre_df = pd.DataFrame(bulk_data['pre'], columns=['pre'])
    on_df = pd.DataFrame(bulk_data['on'], columns=['on'])
    
    # merge 'pre' and 'on' matrices based on patient ID
    merged_df = pd.concat([pre_df, on_df], axis=1)
    
    # perform t test
    ttest_results = {}
    for gene in merged_df.index:
        pre_values = merged_df.loc[gene, 'pre']
        on_values = merged_df.loc[gene, 'on']
        t_statistic, p_value = stats.ttest_rel(pre_values, on_values)
        ttest_results[gene] = {'t_statistic': t_statistic, 'p_value': p_value}
    
    # convert t-test results to a data frame
    ttest_df = pd.DataFrame.from_dict(ttest_results, orient='index')
    
    # filter?
    significant_DEGs = ttest_df[ttest_df['p_value'] < 0.05]
    
    return significant_DEGs

