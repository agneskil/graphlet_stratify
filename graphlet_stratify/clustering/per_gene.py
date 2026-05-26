"""
Per-gene Ward hierarchical clustering of graphlet orbit vectors.

Produces one 'gene<ID>_clusters.csv' per retained gene.

"""

import glob
import logging
import os

import numpy as np
import pandas as pd

from .cluster import normalise_id, ward_cluster, plot_pca

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_sample_files(input_dir: str, id_length: int | None = None) -> tuple[dict, list[str]]:
    """
    Read all .txt Jesse output files in 'input_dir'.

    Returns:
    gene_data   : { gene_id (int) -> [(sample_name, orbit_vector), ...] }
    sample_list : list of normalised sample names
    """
    file_paths  = sorted(glob.glob(os.path.join(input_dir, "*.txt")))
    gene_data   = {}
    sample_list = []

    for fp in file_paths:
        sample_name = normalise_id(os.path.splitext(os.path.basename(fp))[0], id_length)
        sample_list.append(sample_name)
        try:
            df = pd.read_csv(fp, sep=r"\s+", header=None, dtype=float)
            for _, row in df.iterrows():
                gid = int(row.iloc[0])
                vec = row.iloc[1:].values
                gene_data.setdefault(gid, []).append((sample_name, vec))
        except Exception as e:
            log.warning(f"Could not read {fp}: {e}")

    return gene_data, sample_list


def filter_genes(
    gene_data: dict,
    sample_list: list[str],
    min_presence: float,
) -> tuple[dict, int]:
    """Keep genes present in at least ``min_presence`` fraction of samples."""
    min_count = int(np.ceil(min_presence * len(sample_list)))
    retained  = {gid: recs for gid, recs in gene_data.items()
                 if len(recs) >= min_count}
    return retained, min_count


# ─────────────────────────────────────────────────────────────────────────────
# Main per-gene clustering
# ─────────────────────────────────────────────────────────────────────────────

def run_per_gene_clustering(
    input_dir: str,
    output_dir: str,
    n_clusters: int,
    min_presence: float = 0.90,
    outlier_z: float = 5.0,
    save_pca: bool = False,
    subtype_map: dict[str, str] | None = None,
    id_length: int | None = None,
) -> None:
    """
    Cluster samples per gene using their graphlet orbit vectors.

    Parameters:
    input_dir    : folder with Jesse output .txt files (one per sample)
    output_dir   : folder for result CSVs (and PCA PNGs if save_pca)
    n_clusters   : number of Ward clusters to cut into
    min_presence : minimum fraction of samples a gene must appear in
    outlier_z    : z-score threshold for outlier removal
    save_pca     : whether to save PCA plots
    subtype_map  : optional { sample_id: subtype } for PCA colouring
    """
    os.makedirs(output_dir, exist_ok=True)

    log.info(f"Loading sample files from '{input_dir}' ...")
    gene_data, sample_list = load_sample_files(input_dir, id_length)
    log.info(f"  {len(sample_list)} samples, {len(gene_data)} unique genes")

    retained, min_count = filter_genes(gene_data, sample_list, min_presence)
    log.info(
        f"Gene presence filter (>={min_presence*100:.0f}% = {min_count} samples): "
        f"{len(gene_data)} → {len(retained)} genes retained"
    )

    log.info(f"Clustering {len(retained)} genes (k={n_clusters}) ...")
    n_saved = 0

    for gene_id, records in sorted(retained.items()):
        names = [r[0] for r in records]
        X     = np.vstack([r[1] for r in records])
        X_log = np.log1p(X)

        result = ward_cluster(X_log, names, n_clusters, outlier_z)
        if result is None:
            log.debug(f"Gene {gene_id}: skipped — too few samples after outlier removal.")
            continue

        out_csv = os.path.join(output_dir, f"gene{gene_id}_clusters.csv")
        pd.DataFrame({"sample": result["names"], "cluster": result["labels"]})\
          .to_csv(out_csv, index=False)

        if save_pca:
            plot_pca(
                result["X_scaled"], result["labels"], result["names"],
                title=f"Gene {gene_id}  |  k={n_clusters}",
                out_path=os.path.join(output_dir, f"gene{gene_id}_pca.png"),
                subtype_map=subtype_map,
            )

        n_saved += 1
        log.info(
            f"  Gene {gene_id:>6}: n={len(result['names']):>4}  "
            f"outliers_removed={result['n_removed']}  sil={result['sil']:.4f}"
        )

    log.info(f"Done. {n_saved} genes clustered → '{output_dir}'")
