"""
Runs an optimized version of the crispat gaussian-gaussian mixture model that works best for Illumina Single-Cell CRISPR Prep data.
"""
# -----------------------------------------------------------------------------
# Imports
# -----------------------------------------------------------------------------
import os
import crispat
import anndata as ad

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
data_dir = '../example_data/iscp_dragen'

# -----------------------------------------------------------------------------
# Section 2: Load AnnData
# -----------------------------------------------------------------------------
crispat.create_anndata_from_dragen_or_cellranger(data_dir)
h5ad_path = os.path.join(data_dir, 'gRNA_counts.h5ad')
adata = ad.read_h5ad(h5ad_path)
print(adata)

guide_assignments_dir = os.path.join(data_dir, 'results')
os.makedirs(guide_assignments_dir, exist_ok=True)


# -----------------------------------------------------------------------------
# Section 3: Run Gaussian Mixture Modeling
# -----------------------------------------------------------------------------
gauss_output_dir = os.path.join(guide_assignments_dir, 'gauss')
os.makedirs(gauss_output_dir, exist_ok=True)
crispat.ga_gauss(h5ad_path, gauss_output_dir, nonzero=True)


