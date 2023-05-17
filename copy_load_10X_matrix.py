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


'''
load_10X_matrices(matrix_dir)

This function will load a fold of 10X matrix files into a single sparse matrix.

Parameters:
    matrix_dir (str): The directory containing the 10X matrix files.

Returns:
    A single AnnData object containing the concatenated matrices.

'''

def load_10X_matrices(matrix_dir):
    # Open the matrix directory and get the list of files
    matrix_dir = "/home/data/PanCanSC/CRC/GEO/GSE161277/GSE161277_RAW"
    matrix_files = os.listdir(matrix_dir)

    # Initialize the AnnData object
    adata = sc.AnnData()

    # Loop through the files and concatenate the matrices
    mtx_files = [x for x in matrix_files if '.mtx' in x]
        
    # find out whehter files has prefix, use set to remove duplicates
    prefixes = set([x.split('matrix.mtx.gz')[0] for x in mtx_files if '_matrix.mtx.gz' in x]) 

    if len(prefixes) > 0:
        # create a list to hold adatas
        adata_list = []
        matrix_dir_p = os.path.join(matrix_dir)

        # interate through prefixes and index mtx_files        
        for index, prefix in enumerate(prefixes): 
            tmp = sc.read_10x_mtx(matrix_dir_p, prefix=prefix)
            adata_list.append(tmp)



    # concatenate adata_list
    overall_adata = ad.concat(adata_list, join='outer')
    overall_adata.obs_names_make_unique()



    # TESTS
    print(adata_list)
    print(overall_adata.X)
    print(f"OVERALL OBS: \n {overall_adata.obs}")
    print(f"OVERALL VAR: \n {overall_adata.var}")