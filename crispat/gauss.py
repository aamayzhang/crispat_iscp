import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse, stats
from scipy.sparse import csr_matrix
from sklearn import mixture
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time


def plot_loss(losses, gRNA, output_dir, inference_tag):
    '''
    Saves a plot of the loss over the SVI steps
    Args:
        losses: (list) loss over the SVI steps
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    plt.figure(figsize=(8, 3), dpi=300).set_facecolor("white")
    plt.plot(losses)
    plt.xlabel("iters")
    plt.ylabel("loss")
    plt.title("Convergence of SVI for " + gRNA)
    plt.savefig(os.path.join(output_dir, "loss_plots_" + inference_tag, "loss_" + gRNA + ".png"), bbox_inches="tight")
    plt.close()


def plot_fitted_model(data, weights, locs, scales, threshold, gRNA, output_dir):
    '''
    Saves a plot of the data histogram and the fitted mixture model
    Args:
        data: (tensor or array) observed transformed gRNA counts
        weights: (np array) estimated proportion for the two Normal components
        locs: (np array) MAP estimate for the means of the Normal distributions
        scales: (np array) MAP estimate for the scales of the Normal distributions
        threshold: (float) threshold for assigning to the higher-mean component
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    X = np.arange(0, max(data) + 1, 0.01)
    Y1 = weights[0] * stats.norm.pdf(X, locs[0], scales)
    Y2 = weights[1] * stats.norm.pdf(X, locs[1], scales)

    fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
    sns.histplot(data, binwidth=0.1, color='grey', stat="proportion")
    ax.plot(X, Y1, "r-", label="Normal 1")
    ax.plot(X, Y2, "b-", label="Normal 2")
    ax.plot(X, Y1 + Y2, "k--", label="Mixture model")
    ax.set_ylim(0, 1)
    ax.axvline(threshold, c="green", label="Threshold")
    plt.legend()
    plt.title("Gaussian mixture model for " + gRNA)
    plt.ylabel("Probability Density")
    plt.xlabel("Log10 " + gRNA + " UMI counts")
    plt.savefig(os.path.join(output_dir, "fitted_model_plots", "fitted_model_" + gRNA + ".png"), bbox_inches="tight")
    plt.close()


def fit_em(gRNA, adata_crispr, output_dir, nonzero=True):
    '''
    Fits Gaussian mixture model for log10 of UMI counts of one gRNA using EM
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts
        output_dir: (str) directory in which plots will be saved
        nonzero: (bool) if True, only nonzero counts are fit
    Returns:
        DataFrame with perturbed cells and gRNA
    '''
    # 1) pull out raw counts and log-transform
    selected_guide = adata_crispr[:, [gRNA]].X
    counts = selected_guide.toarray().reshape(-1)
    data = np.log10(counts + 1)
    if nonzero:
        mask = data != 0
        data = data[mask]
        counts = counts[mask]
        perturbed_cell_barcodes = adata_crispr.obs_names[mask]

    # 2) fit the mixture
    gmm = mixture.GaussianMixture(
        n_components=2,
        n_init=10,
        covariance_type="tied",
        random_state=0
    )
    gmm.fit(data.reshape(-1, 1))

    # 3) compute posterior & threshold
    posterior = gmm.predict_proba(data.reshape(-1, 1))
    high_component = np.argmax(gmm.means_.flatten())
    threshold_count = counts[posterior[:, high_component] > 0.5].min()
    threshold_log = np.log10(threshold_count + 1)

    # 4) pull out perturbed cells
    perturbed_idx = counts >= threshold_count
    if nonzero:
        perturbed_cells = perturbed_cell_barcodes[perturbed_idx].tolist()
    else:
        perturbed_cells = adata_crispr.obs_names[perturbed_idx].tolist()

    perturbations = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})

    # 5) extract parameters for plotting
    weights = gmm.weights_
    locs = gmm.means_.flatten()
    scales = np.array([np.sqrt(gmm.covariances_.item())])

    # 6) plot the fitted model
    plot_fitted_model(
        data,
        weights,
        locs,
        scales,
        threshold_log,
        gRNA,
        output_dir,
    )

    return perturbations


def ga_gauss(input_file, output_dir, start_gRNA=0, step=None,
             UMI_threshold=0, nonzero=True):
    '''
    Guide assignment using a Gaussian mixture model per batch
    '''
    print('Guide assignment using a Gaussian mixture model per batch')
    adata_crispr = sc.read_h5ad(input_file)

    # create output dirs
    if os.path.exists(output_dir):
        os.makedirs(os.path.join(output_dir, "fitted_model_plots"), exist_ok=True)
    else:
        raise(f'output dir {output_dir} does not exist')

    adata_crispr_batch = adata_crispr
    gRNA_list = adata_crispr_batch.var_names.tolist()
    if step is not None:
        end = min(start_gRNA + step, len(gRNA_list))
        gRNA_list = gRNA_list[start_gRNA:end]

    perturbations = pd.DataFrame(columns=['cell', 'gRNA'])

    for gRNA in gRNA_list:
        df = fit_em(gRNA, adata_crispr_batch, output_dir, nonzero)

        # add raw UMI counts
        read_counts = adata_crispr_batch[df['cell'], [gRNA]].X.toarray().reshape(-1)
        df['read_counts'] = read_counts
        perturbations = pd.concat([perturbations, df], ignore_index=True)

    # filter by UMI_threshold
    if not perturbations.empty:
        perturbations = perturbations[perturbations['read_counts'] >= UMI_threshold]

    perturbations.to_csv(os.path.join(output_dir, 'assignments.csv'), index=False)

    print('Done: outputs are saved in ' + output_dir)


if __name__ == "__main__":
    # example usage:
    # ga_gauss('data/crispr_counts.h5ad', 'results/')
    pass
