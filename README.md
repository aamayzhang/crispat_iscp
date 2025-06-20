<p align="center">
    <img src="https://github.com/velten-group/crispat/blob/main/crispat_logo.png" alt="logo" width="200"/>
</p>

## Installation
Clone this repository and then run `pip install .` inside of the crispat directory.

# crispat: CRISPR guide assignment tool
Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is guide assignment, where cells are assigned to specific guides and corresponding genetic targets. For this, thresholds on the measured gRNA counts per cell are used to distinguish between background contamination and the actual guide presence in a cell. 

This fork of the crispat repo has been prepared specifically to perform guide calling on Illumina Single-Cell Prep (ISCP) samples processed by DRAGEN. 
It does the following

## Gaussian Mixture Model Implementation Details
This fork is customized specifically for processing ISCP data and uses the across-cells implementation of the gaussian mixture model (GMM) with tied covariance, only including cells with nonzero gRNA counts `ga_gauss`. 
Only one GMM implementation is provided here, which uses the scikit-learn implementation with expectation maximization (EM) on nonzero gRNA count data. 
We have found that this model tends to reliably identify the positive distribution under most circumstances. 
The original crispat implementation provided a pyro-based Bayesian variational inference GMM, which did not perform as well as with the EM approach for ISPC samples.

For details on other guide calling approaches provided in the original repo, please refer to [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392). Note that updates have been made here to properly read in DRAGEN-processed ISCP samples, so crispat will not work for these other approaches without modification.

## Getting started
An example use case of crispat is shown in the [`guide_assignment.ipynb`](tutorial/guide_assignment.ipynb) script. 

Tutorials on how to process an ISCP sample can be found here: [`tutorials/iscp_demo.py`](tutorial/downstream_analyses). 

## Documentation

For details on individual functions in the package refer to the [crispat documentation](https://crispat.readthedocs.io/en/latest/index.html).
