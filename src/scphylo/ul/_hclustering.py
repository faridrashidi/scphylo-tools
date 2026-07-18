import numba
import numpy as np
import pandas as pd
import scipy as sp
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import pairwise_distances

import scphylo as scp
from scphylo.external._betabinom import pmf_BetaBinomial


@numba.jit(nopython=True)
def _l1_ignore_na(a, b):
    distance = 0.0
    observed = 0
    for i in range(a.shape[0]):
        if a[i] != 3 and b[i] != 3 and not np.isnan(a[i]) and not np.isnan(b[i]):
            distance += np.abs(a[i] - b[i])
            observed += 1
    if observed == 0:
        return np.nan
    return distance / observed


def _l1_ignore_na_wrapper(x, y, **kwargs):
    return _l1_ignore_na(x, y)


def dist_l1_ignore_na(I_mtr, n_jobs=1):
    """Compute pairwise L1 distances while ignoring missing values."""
    dist = pairwise_distances(
        I_mtr, metric=_l1_ignore_na_wrapper, force_all_finite="allow-nan", n_jobs=n_jobs
    )
    np.fill_diagonal(dist, 0)
    return dist


# https://gist.github.com/FedericoV/0e7d6d8c8794a99a7a42
@numba.jit(nopython=True)
def _cosine_ignore_na(u, v):
    m = u.shape[0]
    udotv = 0
    u_norm = 0
    v_norm = 0
    for i in range(m):
        if (np.isnan(u[i])) or (np.isnan(v[i])):
            continue
        udotv += u[i] * v[i]
        u_norm += u[i] * u[i]
        v_norm += v[i] * v[i]
    u_norm = np.sqrt(u_norm)
    v_norm = np.sqrt(v_norm)
    if (u_norm == 0) or (v_norm == 0):
        ratio = 1.0
    else:
        ratio = 1 - udotv / (u_norm * v_norm)
    if ratio < 0:
        return 0
    return ratio


def _cosine_ignore_na_wrapper(x, y, **kwargs):
    return _cosine_ignore_na(x, y)


def dist_cosine_ignore_na(I_mtr, n_jobs=1):
    """Compute pairwise cosine distances while ignoring missing values."""
    dist = pairwise_distances(
        I_mtr,
        metric=_cosine_ignore_na_wrapper,
        force_all_finite="allow-nan",
        n_jobs=n_jobs,
    )
    np.fill_diagonal(dist, 0)
    return dist


def _dist_dendro(T, V, I_mtr):
    PROB_SEQ_ERROR = 0.001

    def logSum_1(x, y):
        return np.logaddexp(x, y)

    D = np.full(V.shape, np.nan, dtype=np.float64)
    np.divide(V, T, out=D, where=T != 0)

    observed = np.isfinite(D)
    observed_count = observed.sum(axis=0)
    Mu = np.full(D.shape[1], np.nan, dtype=np.float64)
    np.divide(
        np.where(observed, D, 0).sum(axis=0),
        observed_count,
        out=Mu,
        where=observed_count > 0,
    )
    centered = np.where(observed, D - Mu, 0)
    Var = np.full(D.shape[1], np.nan, dtype=np.float64)
    np.divide(
        np.square(centered).sum(axis=0),
        observed_count - 1,
        out=Var,
        where=observed_count > 1,
    )
    ratio = np.full(D.shape[1], np.nan, dtype=np.float64)
    np.divide(
        (1 - Mu) * Mu,
        Var,
        out=ratio,
        where=np.isfinite(Var) & (Var != 0),
    )
    a = (ratio - 1) * Mu
    b = (ratio - 1) * (1 - Mu)
    bad_muts = (
        (a <= 0) | (b <= 0) | np.isnan(a) | np.isnan(b) | np.isinf(a) | np.isinf(b)
    )
    V = V[:, ~bad_muts]
    T = T[:, ~bad_muts]
    # D = D[:, ~bad_muts]
    I_mtr = I_mtr[:, ~bad_muts]
    a = a[~bad_muts]
    b = b[~bad_muts]

    lPz0 = np.zeros(T.shape, dtype=np.float64)
    lPz1 = np.zeros(T.shape, dtype=np.float64)
    for i in range(T.shape[0]):
        for j in range(T.shape[1]):
            if T[i, j] != 0:
                probability0 = sp.stats.binom.pmf(V[i, j], T[i, j], PROB_SEQ_ERROR)
                probability1 = pmf_BetaBinomial(V[i, j], T[i, j], a[j], b[j])
                lPz0[i, j] = -np.inf if probability0 == 0 else np.log(probability0)
                lPz1[i, j] = -np.inf if probability1 == 0 else np.log(probability1)

    Pg = np.sum(I_mtr == 1, axis=0) / I_mtr.shape[0]
    lPg = np.full(Pg.shape, -np.inf, dtype=np.float64)
    np.log(Pg, out=lPg, where=Pg > 0)
    one_minus_Pg = 1 - Pg
    l1Pg = np.full(Pg.shape, -np.inf, dtype=np.float64)
    np.log(one_minus_Pg, out=l1Pg, where=one_minus_Pg > 0)
    lupiall = logSum_1(lPz0 + l1Pg, lPz1 + lPg)

    dist = np.zeros((T.shape[0], T.shape[0]), dtype=np.float64)
    for i in range(T.shape[0]):
        ldowni = logSum_1(lPz0[i, :] + lPz0 + l1Pg, lPz1[i, :] + lPz1 + lPg)
        lupi = logSum_1(lupiall[i, :] + lupiall, ldowni)
        contribution = np.zeros_like(lupi)
        same_infinity = np.isinf(lupi) & (lupi == ldowni)
        np.subtract(lupi, ldowni, out=contribution, where=~same_infinity)
        dist[i, :] = np.sum(contribution, axis=1)

    dist = dist - np.min(dist) + 1
    return dist, bad_muts


def dist_dendro(adata):
    """Compute DENDRO distances and remove unsupported mutations."""
    T = adata.layers["total"]
    V = adata.layers["mutant"]
    G = adata.layers["genotype"]
    G[(G == 1) | (G == 3)] = 1
    G[G == 2] = 0
    dist, bad_muts = _dist_dendro(T, V, G)
    scp.pp.remove_mut_by_list(adata, bad_muts)
    scp.logg.info(f"{sum(bad_muts)} mutations filtered")
    return dist


def hclustering(df, metric="l1", method="ward", return_dist=False):
    """Hierarchical clustering.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        The genotype matrix.
    metric: :obj:`str`, optional
        The metric option. Can be:

            - `l1`
            - `cosine`
    method : :obj:`str`, optional
        The method for the hierarchical clustering, by default "ward"

    Returns
    -------
    :obj:`dict`
        A dictionary in which keys are the number of clusters and
        values are the cluster labels for each item.
    """
    if metric == "l1":
        dist = dist_l1_ignore_na(df.values)
    elif metric == "cosine":
        dist = dist_cosine_ignore_na(df.values)
    else:
        scp.logg.error("Wroing `metric` choice!")
    clust = linkage(dist[np.triu_indices(dist.shape[0], 1)], method=method)
    clusters = {}

    for i in range(2, dist.shape[0]):
        fc = fcluster(clust, i, criterion="maxclust")
        clusters[i] = pd.Series(fc, index=df.index)

    if return_dist:
        return dist

    return clusters
