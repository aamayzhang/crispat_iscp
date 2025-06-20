import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import glob
from scipy.sparse import csr_matrix

def create_anndata_from_csv(csv_file, save_dir):
    '''
    Creates an AnnData object from a csv file (rows: gRNAs, columns: cells; name of first column: gRNA)
    
    Args:
        csv_file (str): path of the csv file
        save_dir (str): directory in which to save the AnnData object that is used as input for guide assignment
    
    Returns:
        None
    '''
    # Load csv files with the gRNA counts
    print('Load data')
    data = pd.read_csv(csv_file).set_index('gRNA')

    # Create anndata object
    print('Create anndata object')
    data = ad.AnnData(data.transpose())
    data.obs['batch'] = 1
    data.X = csr_matrix(data.X)
    
    # Save as h5ad objects
    data.write(save_dir + 'gRNA_counts.h5ad')
    print('Done: AnnData object is saved in ' + save_dir)
    
    
def create_anndata_from_dragen_or_cellranger(input_dir, save_dir=''):
    '''
    Creates an AnnData object from cellranger output (one folder per batch named batch1, batch2,...)
    
    Args:
        batch_list (list): list of batch numbers 
        input_dir (list): directory containing the cellranger output of each batch
        save_dir (str): directory in which to save the AnnData object that is used as input for guide assignment 
    
    Returns:
        None
    '''
    print('Load data')
    # Load data of one batch
    prefix = ''
    save_dir = save_dir if save_dir != '' else input_dir

    if not os.path.exists(input_dir):
        raise FileNotFoundError('Could not find input dir. Ensure path is correct and try again.')

    # First search for dragen filtered fastqs in input_dir
    matrix_file_glob = os.path.join(input_dir, '*scRNA.filtered.matrix.mtx.gz')
    features_file_glob = os.path.join(input_dir, '*scRNA.filtered.features.tsv.gz')
    barcodes_file_glob = os.path.join(input_dir, '*scRNA.filtered.barcodes.tsv.gz')

    if matrix_file_glob:

        matrix_file = glob.glob(matrix_file_glob)[0]
        features_file = glob.glob(features_file_glob)[0]
        barcodes_file = glob.glob(barcodes_file_glob)[0]

    if matrix_file:
        matrix_filename = os.path.basename(matrix_file)
        prefix = matrix_filename.split('matrix.mtx')[0]
        if not (features_file or barcodes_file):
            raise FileNotFoundError('Found matrix, but unable to find a corresponding DRAGEN feature/barcode tsv.gz')
    else:
        print('Unable to find a DRAGEN filtered matrix (*filtered.matrix.mtx.gz) file in the input dir. \n'
              'Falling back to CellRanger matrix output format.')

    adata = sc.read_10x_mtx(input_dir,
                            var_names='gene_symbols',
                            gex_only=False,
                            prefix=prefix)

    adata.obs_names = [name.split('-')[0] for name in adata.obs_names]

    # Concatenate batches into one AnnData object
    print('Create concatenated anndata object')

    # Subset to CRISPR gRNA features
    mask = adata.var["feature_types"].str.contains(r"crispr|grna", case=False, na=False)

    # subset your AnnData accordingly
    adata_crispr = adata[:, mask].copy()
    
    # Save as h5ad object
    adata_crispr.write(os.path.join(save_dir, 'gRNA_counts.h5ad'))
    print('Done: AnnData object is saved in ' + save_dir)
    
