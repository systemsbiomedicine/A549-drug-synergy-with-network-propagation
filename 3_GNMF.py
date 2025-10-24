"""
Graph‑regularised Non‑Negative Matrix Factorisation (GNMF) workflow
------------------------------------------------------------------
This annotated script accompanies the methods described in the first
part of the manuscript.  It reproduces the complete
pipeline that:
  1.  Loads −log‑transformed network‑propagated profiles for 607
      drug combinations across 2100 genes.
  2.  Builds a protein‑protein‑interaction (PPI) adjacency matrix for
      the A549‑specific cancer sub‑network (2100 nodes).
  3.  Performs Graph NMF (GNMF) across a range of component counts *k*
      while enforcing smoothness on the PPI graph (λ = 0.9).
  4.  Uses the KneeLocator elbow method to pick the optimal k.
  5.  Runs GNMF one final time with the chosen k and extracts the W & H
      matrices.
"""

# %% ------------------------------------------------------------------
# 1. Imports
# ---------------------------------------------------------------------
# WHY:  Standard scientific Python libraries cover numerical operations
#       (NumPy, SciPy), data frames (Pandas), sparse matrices, and
#       plotting.  Additional tools (kneed, NetworkX) support elbow‑
#       point detection and graph handling respectively.
import numpy as np              # base numerical package
import pandas as pd             # tabular data handling
import scipy.sparse as sp       # sparse‑matrix support
from numpy.linalg import norm   # fast Frobenius norms
from sklearn.utils import (
    check_random_state,
    extmath
)  # misc. sklearn utilities
from sklearn.decomposition._nmf import _initialize_nmf as sk_initialize_nmf
import warnings                 # surface convergence warnings
!pip install kneed
from kneed import KneeLocator   # detect an ‘elbow’ in error curve
import networkx as nx           # graph utilities
import matplotlib.pyplot as plt # basic plotting
# NOTE: seaborn is only used for pretty cluster heatmaps later on
import seaborn as sns

# %% ------------------------------------------------------------------
# 2. Load the −log‑transformed propagation matrix (X)
# ---------------------------------------------------------------------
# WHY:  X has shape (n_samples × n_features) = (607 × 2 100) after
#       transpose.  We store samples (drug combinations) in rows and
#       genes in columns so factorisation splits genes into metagenes
#       and samples into weighted mixtures of those metagenes.
file_path = '/filesData/Network_Popagation_Result_alpha0.5(-log).csv'
X = pd.read_csv(file_path, index_col=0).T  # transpose so rows = samples
X_np = X.to_numpy()                        # convert to ndarray for speed
print('Shape of input matrix X (samples × genes):', X_np.shape)

# %% ------------------------------------------------------------------
# 3. Helper utility functions
# ---------------------------------------------------------------------
# WHY:  We re‑implement several private Scikit‑learn helpers so that the
#       whole pipeline is self‑contained and reproducible without
#       relying on sklearn internals that may change in future.

def check_non_negative(mat, context):
    """Raise if *mat* contains negative entries (GNMF assumes ≥0)."""
    arr = mat.data if sp.issparse(mat) else mat
    if (arr < 0).any():
        raise ValueError(f"Negative values passed to {context}")

def _sparseness(x):
    """Hoyer’s sparsity metric (0 = dense, 1 = maximally sparse)."""
    if sp.issparse(x):
        x = x.data
    sqrt_n = np.sqrt(len(x))
    return (sqrt_n - np.linalg.norm(x, 1) / np.linalg.norm(x)) / (sqrt_n - 1)

def _initialize_nmf_custom(X, n_components, eps=1e-6, random_state=None):
    """NNDSVD initialisation (copied to remove sklearn private import)."""
    U, S, Vt = extmath.randomized_svd(X, n_components)
    W = np.zeros_like(U)
    H = np.zeros_like(Vt)
    # first component is taken as absolute of leading singular vectors
    W[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
    H[0, :] = np.sqrt(S[0]) * np.abs(Vt[0, :])
    # remaining components follow Boutsidis et al.
    for j in range(1, n_components):
        u, v = U[:, j], Vt[j, :]
        u_p, v_p = np.maximum(u, 0), np.maximum(v, 0)
        u_n, v_n = np.maximum(-u, 0), np.maximum(-v, 0)
        m_p = np.linalg.norm(u_p) * np.linalg.norm(v_p)
        m_n = np.linalg.norm(u_n) * np.linalg.norm(v_n)
        if m_p > m_n:
            uu, vv, sigma = u_p, v_p, m_p
        else:
            uu, vv, sigma = u_n, v_n, m_n
        W[:, j] = np.sqrt(S[j] * sigma) * uu / (np.linalg.norm(uu) + eps)
        H[j, :] = np.sqrt(S[j] * sigma) * vv / (np.linalg.norm(vv) + eps)
    W[W < eps] = eps
    H[H < eps] = eps
    return W, H

def NBS_init(X, n_components, random_state=None):
    """Wrapper choosing NNDSVD or random non‑negative seed."""
    rng = check_random_state(random_state)
    if n_components < min(X.shape):
        return _initialize_nmf_custom(X, n_components, random_state=rng)
    # fallback: fully random non‑negative
    return np.abs(rng.randn(X.shape[0], n_components)), np.abs(rng.randn(n_components, X.shape[1]))

def as_float_array(mat):
    """Ensure *mat* is float64 dense or CSR for numerical stability."""
    return np.asarray(mat, dtype=np.float64) if not sp.issparse(mat) else mat.astype(np.float64)

# %% ------------------------------------------------------------------
# 4. Core Graph‑regularised NMF implementation
# ---------------------------------------------------------------------
# WHY:  GNMF constrains classic NMF with two Laplacian terms (Lp/Lm)
#       encouraging smooth H along graph edges controlled by λ (lambd).

def GNMF(X, L, lambd=0.0, n_components=10, tol=1e-4, max_iter=200, verbose=False):
    """Factorise X ≈ WH subject to graph Laplacian smoothness on H."""
    X = as_float_array(X)
    check_non_negative(X, 'GNMF.fit')
    n_samples, n_features = X.shape
    W, H = NBS_init(X, n_components)
    eps = 1e-8  # small constant to avoid /0
    Lp = (np.abs(L) + L) / 2.0  # positive part of Laplacian
    Lm = (np.abs(L) - L) / 2.0  # negative part
    err_prev = norm(X - W @ H)
    rng = check_random_state(None)
    for i in range(1, max_iter + 1):
        # multiplicative H update with graph regularisation
        H *= (lambd * H @ Lm + W.T @ (X / (W @ H + eps))) / (
            lambd * H @ Lp + W.T @ np.ones_like(X) + eps)
        H[H < eps] = eps
        # multiplicative W update (standard)
        W *= ((X / (W @ H + eps)) @ H.T) / (np.ones_like(X) @ H.T + eps)
        W[W < eps] = eps
        err = norm(X - W @ H)
        if verbose and i % 50 == 0:
            print(f'Iter {i:4d}/{max_iter} – reconstruction err = {err:.4e}')
        # early stopping
        if abs(err_prev - err) < tol:
            if verbose:
                print('Tolerance reached – stopping.')
            break
        err_prev = err
    else:
        warnings.warn('Maximum iterations reached before convergence')
    return W, H, err_prev

# %% ------------------------------------------------------------------
# 5. Build the gene‑interaction Laplacian (L)
# ---------------------------------------------------------------------
# WHY:  GNMF enforces metagenes to respect biological proximity.  We
#       construct L from the cancer‑specific PPI edge list.
ppi_path = '/filesData/gene_interaction_cancer_subnetwork.csv'
ppi_edges = pd.read_csv(ppi_path, header=None, names=["Source", "Target"])
G = nx.from_pandas_edgelist(ppi_edges, 'Source', 'Target')
L = nx.laplacian_matrix(G).todense()  # shape (2100 × 2100)
print('Graph built with', G.number_of_nodes(), 'nodes and', G.number_of_edges(), 'edges')

# %% ------------------------------------------------------------------
# 6. Determine optimal number of components *k* (elbow on RE curve)
# ---------------------------------------------------------------------
# WHY:  A too‑small k underfits (high error); a too‑large k overfits +
#       worsens interpretability.  We scan k=1…60 and pick elbow.
max_k = 60
reconstruction_errors = np.empty(max_k)
for k in range(1, max_k + 1):
    _, _, err = GNMF(X_np, L, lambd=0.9, n_components=k, tol=1e-4, max_iter=500)
    reconstruction_errors[k - 1] = err
    if k % 5 == 0:
        print(f'GNMF finished for k = {k}')
# find elbow (convex decreasing curve)
kn = KneeLocator(range(1, max_k + 1), reconstruction_errors,
                 curve='convex', direction='decreasing')
optimal_k = kn.elbow
print('Optimal number of components (k):', optimal_k)

# Visualise the elbow (optional)
plt.figure(figsize=(8, 5))
plt.plot(range(1, max_k + 1), reconstruction_errors, marker='o')
plt.axvline(optimal_k, linestyle='--')
plt.xlabel('k (components)')
plt.ylabel('Reconstruction error')
plt.title('Elbow method for GNMF')
plt.tight_layout()
plt.show()

# %% ------------------------------------------------------------------
# 7. Run GNMF once with optimal k and λ = 0.9
# ---------------------------------------------------------------------
H, W, final_err = GNMF(X_np, L, lambd=0.9, n_components=optimal_k,
                       tol=1e-4, max_iter=3000, verbose=True)
print('Final GNMF reconstruction error:', final_err)

df_W = pd.DataFrame(W.T, index=range(2100), columns=range(26))
df_H = pd.DataFrame(H.T, index=range(26), columns=range(607))
df_W.to_csv('W.csv')
df_H.to_csv('H.csv')
