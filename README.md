<p align="center">
    <img src="https://github.com/velten-group/crispat/blob/main/crispat_logo.png" alt="logo" width="200"/>
</p>

# crispat: CRISPR guide assignment tool
Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is guide assignment, where cells are assigned to specific guides and corresponding genetic targets. For this, thresholds on the measured gRNA counts per cell are used to distinguish between background contamination and the actual guide presence in a cell. However, lots of different guide assignment strategies and thresholds are used by different labs without any guidance on what model or threshold to choose when. 

As demonstrated on *low MOI CRISPRi* screens in [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392) the choice of guide assignment strategy strongly influences the results, highlighting the need to choose a suitable strategy for the data at hand for a reliable and powerful analysis of the data. To help with this choice the **crispat** package implements 11 different assignment methods and facilitates their comparison. 

## Installation
To install the Illumina Single-Cell Prep (ISCP) version of crispat, clone this repository and then run `pip install .` inside of the crispat directory.

## Recommended Guide assignment methods
This fork is customized specifically for processing ISCP data and uses the across-cells implementation of the gaussian mixture model with tied covariance, only including cells with nonzero gRNA counts `ga_gauss`

For details on the individual methods please refer to [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392).

## Getting started
An example use case of crispat is shown in the [`guide_assignment.ipynb`](tutorial/guide_assignment.ipynb) script. 

Tutorials on how to process an ISCP sample can be found here: [`tutorials/iscp_demo.py`](tutorial/downstream_analyses). 

## Documentation

For details on individual functions in the package refer to the [crispat documentation](https://crispat.readthedocs.io/en/latest/index.html).
