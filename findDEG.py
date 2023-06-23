import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

def find_cluster_DEGs(adata, cluster_label, condition_key):
    # filter cells based on cluster
    cluster_mask = adata.obs['cluster'] == cluster_label
    adata_cluster = adata[cluster_mask].copy()
    
    # create pseudo-bulk RNA data for each sample
    samples = ['pre', 'on']
    bulk_data = {}
    for sample in samples:
        sample_mask = adata_cluster.obs[condition_key] == sample
        bulk_data[sample] = np.array(adata_cluster.X[sample_mask].sum(axis=0)).flatten()
    
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

