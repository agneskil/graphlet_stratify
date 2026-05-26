"""
Shared helpers used by both per-gene and per-orbit-type clustering.

"""

import logging

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

log = logging.getLogger(__name__)


def normalise_id(s: str, id_length: int | None = None) -> str:
    s = str(s).strip()
    return s[:id_length] if id_length else s


def remove_outliers(
    X_log: np.ndarray,
    names: list[str],
    outlier_z: float,
) -> tuple[np.ndarray, list[str], int]:
    """
    Remove samples whose row L2-norm is more than ``outlier_z`` standard
    deviations from the mean.

    Returns: (filtered_matrix, filtered_names, n_removed).
    """
    norms = np.linalg.norm(X_log, axis=1)
    z     = (norms - norms.mean()) / (norms.std() + 1e-9)
    keep  = np.abs(z) < outlier_z
    return (
        X_log[keep],
        [n for n, k in zip(names, keep) if k],
        int((~keep).sum()),
    )


def ward_cluster(
    X_log: np.ndarray,
    names: list[str],
    n_clusters: int,
    outlier_z: float,
) -> dict | None:
    """
    Full Ward hierarchical clustering pipeline for one feature matrix.

    Steps:  outlier removal → StandardScaler → Ward linkage → cut tree
            → silhouette score.

    Returns:
    dict with keys: names, labels (0-based), sil, X_scaled, n_removed
    or None if too few samples remain after outlier removal.
    """
    X_log, names, n_removed = remove_outliers(X_log, names, outlier_z)

    if len(names) < n_clusters:
        return None

    X_scaled = StandardScaler().fit_transform(X_log)
    Z        = linkage(X_scaled, method="ward")
    labels   = fcluster(Z, t=n_clusters, criterion="maxclust") - 1  # 0-based

    try:
        sil = silhouette_score(X_scaled, labels)
    except Exception:
        sil = float("nan")

    return {
        "names"    : names,
        "labels"   : labels,
        "sil"      : sil,
        "X_scaled" : X_scaled,
        "n_removed": n_removed,
    }


PALETTE = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B2",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
]


def plot_pca(
    X_scaled: np.ndarray,
    labels: np.ndarray,
    names: list[str],
    title: str,
    out_path: str,
    subtype_map: dict[str, str] | None = None,
) -> None:
    """
    Two-panel PCA: left coloured by subtype, right by cluster assignment.
    Subtype colours are assigned automatically from PALETTE in sorted order.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA

    if X_scaled.shape[0] < 2 or X_scaled.var(axis=0).sum() == 0:
        log.warning(f"Skipping PCA plot '{out_path}' — data has no variance.")
        return

    pca   = PCA(n_components=2, random_state=42)
    X_pca = pca.fit_transform(X_scaled)
    var   = pca.explained_variance_ratio_ * 100

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(title, fontsize=13, fontweight="bold")

    # Left: subtype labels (dynamically coloured from PALETTE)
    ax = axes[0]
    if subtype_map:
        subtype_labels = [subtype_map.get(n, "Unknown") for n in names]
        unique_subtypes = sorted(set(subtype_labels))
        subtype_colors  = {s: PALETTE[i % len(PALETTE)]
                           for i, s in enumerate(unique_subtypes)}
        for subtype in unique_subtypes:
            mask  = [s == subtype for s in subtype_labels]
            ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                       color=subtype_colors[subtype], s=30, alpha=0.7,
                       edgecolors="white", linewidths=0.3,
                       label=f"{subtype}  (n={sum(mask)})")
        ax.legend(fontsize=8)
        ax.set_title("Subtype labels")
    else:
        ax.scatter(X_pca[:, 0], X_pca[:, 1], s=30, alpha=0.7)
        ax.set_title("PCA (no subtype labels provided)")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% var)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% var)")
    ax.grid(True, alpha=0.3)

    # Right: cluster assignments
    ax = axes[1]
    for c in sorted(set(labels)):
        mask = labels == c
        ax.scatter(X_pca[mask, 0], X_pca[mask, 1],
                   color=PALETTE[c % len(PALETTE)], s=30, alpha=0.7,
                   edgecolors="white", linewidths=0.3,
                   label=f"Cluster {c}  (n={mask.sum()})")
    ax.set_xlabel(f"PC1 ({var[0]:.1f}% var)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}% var)")
    ax.set_title("Assigned clusters")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close()
