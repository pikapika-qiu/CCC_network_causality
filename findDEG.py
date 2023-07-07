import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

def find_cluster_DEGs_pairwise(adata, cluster_label, condition_key):
    '''
    This function will find differentially expressed genes between two conditions for a given cluster.
    Steps in the process:
        1. Identify cells from a sample that belong to a specific cluster.
        2. Create pseudo-bulk RNA data for each sample.
        3. Match samples from the same patient.
        4. Perform pairwise t-test between two conditions for each gene.
    '''

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

    # Loop through all genes
    # Extract matching pseudo-bulk RNA data for the gene in all patients
    # Perform pairwise t-test between two conditions for the gene; Check scipy for pair-wise t-test
    df_tpm = pd.DataFrame(index=adata.obs['sample_id'].unique(), columns=adata.var_names)
    df_obs = pd.DataFrame(index=adata.obs['sample_id'].unique(), columns=['patient_id', 'timepoint', 'batch'])

    for sample in adata.obs['sample_id'].unique():
        tpm = np.sum(adata.X[adata.obs['sample_id'] == sample, :], axis=0)
        tpm = tpm / np.sum(tpm) * 1e6  # Normalize to TPM per cell
        df_tpm.loc[sample, :] = np.log2(tpm + 1)

        # Populate df_obs
        df_obs.loc[sample, 'patient_id'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'patient_id'].unique()[0]
        df_obs.loc[sample, 'timepoint'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'timepoint'].unique()[0]
        df_obs.loc[sample, 'batch'] = adata.obs.loc[adata.obs['sample_id'] == sample, 'batch'].unique()[0]

    # Create an AnnData object for the pseudo-bulk RNA data
    adata_sample_tpm = sc.AnnData(X=df_tpm.values, obs=df_obs, var=adata.var)

    # Perform t-test
    cluster_data = adata_sample_tpm[adata_sample_tpm.obs[condition_key].isin(['pre', 'on'])]

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
