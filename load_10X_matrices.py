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