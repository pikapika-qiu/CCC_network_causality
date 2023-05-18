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
import AnnData as ad

def load_10X_matrices(matrix_dir):
    '''
    load_10X_matrices(matrix_dir)
    
    This function will load a fold of 10X matrix files into a single sparse matrix.

    Parameters:
        matrix_dir (str): The directory containing the 10X matrix files.
    
    Returns:
        A single AnnData object containing the concatenated matrices.

    '''
    # Open the matrix directory and get the list of files
    matrix_files = os.listdir(matrix_dir)

    # Initialize the AnnData object
    adata = sc.AnnData()
    
    # Loop through the files and concatenate the matrices
    mtx_files = [x for x in matrix_files if x.find('.mtx') > -1]
    
    # find out whehter files has prefix, use set to remove duplicates
    prefixes = set([x for x in mtx_files if x.split('matrix.mtx')[0] != ''])  

    if len(prefixes) > 0:
        # create a list to hold adatas
        adata_list = []
        matrix_dir_p = os.path(matrix_dir)

        # interate through prefixes and index mtx_files        
        for index, prefix in enumerate(prefixes): 
            tmp = sc.read_10x_mtx(matrix_dir_p, prefix=prefix)
            adata_list.append(tmp)

    # concatenate adata_list
    overall_adata = ad.concatenate(adata_list, join='outer')
    return overall_adata
