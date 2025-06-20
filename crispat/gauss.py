import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from scipy import sparse, stats
from scipy.sparse import csr_matrix
from sklearn import mixture
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate


@config_enumerate
def model(data):
    '''
    Gaussian-Gaussian Mixture Model 
    '''
    # Global variables
    weights = pyro.sample("weights", dist.Dirichlet(torch.tensor([0.99, 0.01])))
    with pyro.plate("components", 2):
        locs = pyro.sample("locs", dist.Normal(1.0, 2.0))

    scale = pyro.sample("scales", dist.LogNormal(-3.0, 2.0))

    with pyro.plate("data", len(data)):
        # Local variables
        assignment = pyro.sample("assignment", dist.Categorical(weights))
        pyro.sample("obs", dist.Normal(locs[assignment], scale), obs=data)


def init_loc_fn(site):
    '''
    Define initial parameter values
    '''
    if site["name"] == "weights":
        return torch.tensor([0.99, 0.01])
    if site["name"] == "locs":
        return torch.tensor([0.0, 1.0])
    if site["name"] == "scales":
        return torch.tensor([0.01])
    raise ValueError(site["name"])


def initialize(seed, optim, elbo, data):
    '''
    Initialization for SVI
    Args:
        seed: (str) seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data: (tensor) observed transformed gRNA counts
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["weights", "locs", "scales"]),
        init_loc_fn=init_loc_fn,
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data)


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


def prob_normal_component(X, weights, locs, scales):
    '''
    Calculates the probability for belonging to the Gaussian component given observations
    Args:
        X: (list) values for which the probability is calculated
        weights: (np array) estimated proportions for the two components
        locs: (np array) MAP estimates for the means
        scales: (np array) MAP estimates for the scales
    Returns:
        Boolean mask of whether each X favors the higher-mean component
    '''
    # determine which component has the higher mean
    if locs[0] > locs[1]:
        h, l = 0, 1
    else:
        h, l = 1, 0

    high = dist.Normal(torch.tensor(locs[h]), torch.tensor(scales)).log_prob(torch.tensor(X)) + np.log(weights[h])
    low = dist.Normal(torch.tensor(locs[l]), torch.tensor(scales)).log_prob(torch.tensor(X)) + np.log(weights[l])

    return (high > low).numpy()


def fit_em(gRNA, adata_crispr, output_dir, nonzero):
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


def ga_gauss(input_file, output_dir, start_gRNA=0, step=None, batch_list=None,
             UMI_threshold=0, n_iter=250, nonzero=False):
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
        UMI_counts = adata_crispr_batch[df['cell'], [gRNA]].X.toarray().reshape(-1)
        df['UMI_counts'] = UMI_counts
        perturbations = pd.concat([perturbations, df], ignore_index=True)

    # filter by UMI_threshold
    if not perturbations.empty:
        perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]

    perturbations.to_csv(os.path.join(output_dir, 'assignments.csv'), index=False)

    print('Done: outputs are saved in ' + output_dir)


if __name__ == "__main__":
    # example usage:
    # ga_gauss('data/crispr_counts.h5ad', 'results/', inference='em')
    pass
