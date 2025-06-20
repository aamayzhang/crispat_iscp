<p align="center">
    <img src="https://github.com/velten-group/crispat/blob/main/crispat_logo.png" alt="logo" width="200"/>
</p>

# crispat: CRISPR guide assignment tool
This fork of the [crispat repo]([url](https://github.com/velten-group/crispat)) was prepared to perform guide calling on Illumina Single-Cell Prep (ISCP) samples processed by DRAGEN. 

It does the following: 
- Read in DRAGEN filtered matrix files from an DRAGEN output directory
- Creates and reads in a .h5ad file (via scanpy)
- Performs 2-component Gaussian-Gaussian mixture modeling on the gRNA counts
- Outputs gRNA assignments for each cell barcode along with plots showing the fit of the model to the raw counts

# Installation
Clone this repository and then run `pip install .` inside of the crispat directory.

## Gaussian Mixture Model Implementation Details
This fork is customized specifically for processing ISCP data and uses the across-cells implementation of the gaussian mixture model (GMM) with tied covariance, only including cells with nonzero gRNA counts `ga_gauss`. 
Only one GMM implementation is provided here, which uses the scikit-learn implementation with expectation maximization (EM) on nonzero gRNA count data. 
We have found that this model tends to reliably identify the positive distribution under most circumstances. 
The original crispat implementation provided a pyro-based Bayesian variational inference GMM. We found that this approach not perform as well as with EM for ISPC samples, so it was removed here for simplicity.

For details on other guide calling approaches provided in the original repo, please refer to [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392). 
Note that updates have been made here to properly read in DRAGEN-processed ISCP samples, so crispat will not work for these other approaches without modification.

## Getting started
An example use case of crispat is shown in the [`guide_assignment.ipynb`](tutorial/guide_assignment.ipynb) script. 

Tutorials on how to process an ISCP sample can be found here: [`tutorials/iscp_demo.py`](tutorial/downstream_analyses). 

## Documentation

For details on individual functions in the package refer to the [crispat documentation](https://crispat.readthedocs.io/en/latest/index.html).
