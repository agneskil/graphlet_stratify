"""
Per-orbit-type Ward hierarchical clustering of samples.

For each orbit index, builds a (samples × genes) matrix of counts for that
orbit type, then clusters samples.  Produces one 'orbit<N>_clusters.csv'
per orbit type.
"""

import glob
import logging
import os

import numpy as np
import pandas as pd

from .cluster import normalise_id, ward_cluster, plot_pca

log = logging.getLogger(__name__)


def load_sample_files(input_dir: str, id_length: int | None = None) -> tuple[dict, list[str], int | None]:
    """
    Read all .txt Jesse output files in 'input_dir'.

    Returns:
    gene_data   : { gene_id -> [(sample_name, orbit_vector), ...] }
    sample_list : list of normalised sample names
    n_orbits    : number of orbit dimensions
    """
    file_paths  = sorted(glob.glob(os.path.join(input_dir, "*.txt")))
    gene_data   = {}
    sample_list = []
    n_orbits    = None

    for fp in file_paths:
        sample_name = normalise_id(os.path.splitext(os.path.basename(fp))[0], id_length)
        sample_list.append(sample_name)
        try:
            df = pd.read_csv(fp, sep=r"\s+", header=None, dtype=float)
            for _, row in df.iterrows():
                gid = int(row.iloc[0])
                vec = row.iloc[1:].values
                if n_orbits is None:
                    n_orbits = len(vec)
                gene_data.setdefault(gid, []).append((sample_name, vec))
        except Exception as e:
            log.warning(f"Could not read {fp}: {e}")

    return gene_data, sample_list, n_orbits


def filter_genes(
    gene_data: dict,
    sample_list: list[str],
    min_presence: float,
) -> tuple[dict, int]:
    min_count = int(np.ceil(min_presence * len(sample_list)))
    retained  = {gid: recs for gid, recs in gene_data.items()
                 if len(recs) >= min_count}
    return retained, min_count


def build_orbit_matrices(
    retained: dict,
    sample_list: list[str],
    n_orbits: int,
) -> list[pd.DataFrame]:
    """
    Build one (samples × genes) DataFrame per orbit index.
    Missing values (gene absent in a sample) are filled with 0.
    """
    gene_ids = sorted(retained.keys())
    matrices = []

    for orbit_idx in range(n_orbits):
        col_data = {
            gid: {sn: vec[orbit_idx] for sn, vec in retained[gid]}
            for gid in gene_ids
        }
        df = pd.DataFrame(col_data, index=sample_list).fillna(0.0)
        matrices.append(df)

    return matrices

# ─────────────────────────────────────────────────────────────────────────────
# Main per-orbit clustering
# ─────────────────────────────────────────────────────────────────────────────

def run_per_orbit_clustering(
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
    Cluster samples per orbit type using gene-level orbit counts.

    Parameters:
    input_dir    : folder with Jesse output .txt files
    output_dir   : folder for result CSVs (and PCA PNGs if save_pca)
    n_clusters   : number of Ward clusters
    min_presence : minimum fraction of samples a gene must appear in
    outlier_z    : z-score threshold for outlier removal
    save_pca     : whether to save PCA plots
    subtype_map    : optional { sample_id: subtype } for PCA colouring
    """
    os.makedirs(output_dir, exist_ok=True)

    log.info(f"Loading sample files from '{input_dir}' ...")
    gene_data, sample_list, n_orbits = load_sample_files(input_dir, id_length)
    log.info(f"  {len(sample_list)} samples, {len(gene_data)} genes, {n_orbits} orbit types")

    retained, min_count = filter_genes(gene_data, sample_list, min_presence)
    log.info(
        f"Gene presence filter (>={min_presence*100:.0f}% = {min_count} samples): "
        f"{len(gene_data)} → {len(retained)} genes retained"
    )
    if not retained:
        raise RuntimeError("No genes passed the presence filter. Lower --min_presence.")

    log.info(f"Building {n_orbits} orbit matrices ({len(sample_list)} × {len(retained)}) ...")
    orbit_matrices = build_orbit_matrices(retained, sample_list, n_orbits)

    log.info(f"Clustering {n_orbits} orbit types (k={n_clusters}) ...")
    n_saved = 0

    for orbit_idx, df_orbit in enumerate(orbit_matrices):
        names = list(df_orbit.index)
        X_log = np.log1p(df_orbit.values)

        result = ward_cluster(X_log, names, n_clusters, outlier_z)
        if result is None:
            log.debug(f"Orbit {orbit_idx}: skipped — too few samples.")
            continue

        out_csv = os.path.join(output_dir, f"orbit{orbit_idx}_clusters.csv")
        pd.DataFrame({"sample": result["names"], "cluster": result["labels"]})\
          .to_csv(out_csv, index=False)

        if save_pca:
            plot_pca(
                result["X_scaled"], result["labels"], result["names"],
                title=f"Orbit {orbit_idx}  |  k={n_clusters}",
                out_path=os.path.join(output_dir, f"orbit{orbit_idx}_pca.png"),
                subtype_map=subtype_map,
            )

        n_saved += 1
        log.info(
            f"  Orbit {orbit_idx:>3}: n={len(result['names']):>4}  "
            f"outliers_removed={result['n_removed']}  sil={result['sil']:.4f}"
        )

    log.info(f"Done. {n_saved} orbit types clustered → '{output_dir}'")
