"""
Validate cluster assignments against known subtype labels.

Works for both per-gene (gene*_clusters.csv) and per-orbit-type
(orbit*_clusters.csv) result sets.

Metrics used:
- ARI  (Adjusted Rand Index)      — measures cluster/subtype agreement,
                                    corrected for chance [-1, 1]
- NMI  (Normalised Mutual Info)   — information-theoretic overlap [0, 1]
- Jaccard similarity              — per cluster-subtype pair overlap [0, 1]
- % correctly assigned            — after Hungarian matching of clusters
                                    to subtypes
"""

import glob
import logging
import os
import re

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def normalise_id(s: str, id_length: int | None = None) -> str:
    s = str(s).strip()
    return s[:id_length] if id_length else s


def load_label_file(subtype_file: str, id_length: int | None = None) -> pd.Series:
    """
    Load a two-column label file (TSV or CSV).
    Column 0 = sample ID, Column 1 = subtype label.
    Returns a Series indexed by normalised sample ID.
    """
    sep = "\t" if subtype_file.endswith(".tsv") else ","
    df  = pd.read_csv(subtype_file, sep=sep, header=0, usecols=[0, 1])
    df.columns = ["sample_id", "subtype"]
    df = df.dropna(subset=["subtype"])
    df["sample_id"] = df["sample_id"].astype(str).apply(lambda s: normalise_id(s, id_length))
    df = df.drop_duplicates(subset="sample_id", keep="first")
    return df.set_index("sample_id")["subtype"]


def _jaccard(set_a: set, set_b: set) -> float:
    union = len(set_a | set_b)
    return len(set_a & set_b) / union if union else 0.0


def _build_jaccard_matrix(
    cluster_labels: np.ndarray,
    subtype_labels: np.ndarray,
    cluster_ids: list,
    subtype_ids: list,
) -> pd.DataFrame:
    c_sets = {c: set(np.where(cluster_labels == c)[0]) for c in cluster_ids}
    s_sets = {s: set(np.where(subtype_labels == s)[0]) for s in subtype_ids}
    mat = np.zeros((len(cluster_ids), len(subtype_ids)))
    for i, c in enumerate(cluster_ids):
        for j, s in enumerate(subtype_ids):
            mat[i, j] = _jaccard(c_sets[c], s_sets[s])
    return pd.DataFrame(mat, index=cluster_ids, columns=subtype_ids)


def _hungarian_match(jac_df: pd.DataFrame) -> list[tuple]:
    row_ind, col_ind = linear_sum_assignment(1 - jac_df.values)
    return [(jac_df.index[r], jac_df.columns[c], jac_df.values[r, c])
            for r, c in zip(row_ind, col_ind)]


# ─────────────────────────────────────────────────────────────────────────────
# Core validation
# ─────────────────────────────────────────────────────────────────────────────

def validate_cluster_file(
    cluster_csv: str,
    label_series: pd.Series,
    id_length: int | None = None,
) -> dict:
    """
    Compute validation metrics for one cluster assignment CSV.

    Returns:
    dict with keys: status, n_overlap, ari, nmi, mean_jaccard,
                    max_jaccard, pct_correct, matched, jac_df
    """
    df = pd.read_csv(cluster_csv)
    cluster_col = "cluster" if "cluster" in df.columns else "assignment"
    df = df[["sample", cluster_col]].dropna()
    df["sample"]     = df["sample"].astype(str).apply(lambda s: normalise_id(s, id_length))
    df[cluster_col]  = df[cluster_col].astype(int)
    df = df.drop_duplicates(subset="sample", keep="first")

    label_norm = label_series[~label_series.index.duplicated(keep="first")]
    cluster_series = df.set_index("sample")[cluster_col]
    common = cluster_series.index.intersection(label_norm.index)

    if len(common) < 2:
        return {"status": "too_few_overlap", "n_overlap": len(common)}

    cl = cluster_series.loc[common].values.astype(int)
    pa = label_norm.loc[common].values

    subtypes    = sorted(pd.unique(pa))
    cluster_ids = sorted(pd.unique(cl))

    if len(cluster_ids) < 2:
        return {"status": "one_cluster_only", "n_overlap": len(common)}

    pa_enc = np.array(pd.Categorical(pa, categories=subtypes).codes)

    jac_df  = _build_jaccard_matrix(cl, pa_enc, cluster_ids,
                                    list(range(len(subtypes))))
    jac_df.columns = subtypes

    matched    = _hungarian_match(jac_df)
    jac_scores = [j for _, _, j in matched]
    mapping    = {c: s for c, s, _ in matched}
    predicted  = np.array([mapping.get(c) for c in cl])

    return {
        "status"      : "ok",
        "n_overlap"   : len(common),
        "ari"         : adjusted_rand_score(pa_enc, cl),
        "nmi"         : normalized_mutual_info_score(pa_enc, cl,
                                                     average_method="arithmetic"),
        "mean_jaccard": float(np.mean(jac_scores)),
        "max_jaccard" : float(np.max(jac_scores)),
        "pct_correct" : int((predicted == pa).sum()) / len(common) * 100,
        "matched"     : matched,
        "jac_df"      : jac_df,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Batch runners
# ─────────────────────────────────────────────────────────────────────────────

def _extract_id(filepath: str, prefix: str) -> int | str:
    m = re.search(rf"{prefix}(\d+)_clusters\.csv$", os.path.basename(filepath))
    return int(m.group(1)) if m else os.path.basename(filepath)


def run_validation(
    input_dir: str,
    output_dir: str,
    subtype_file: str,
    mode: str = "gene",
    save_jaccard: bool = False,
    id_length: int | None = None,
) -> pd.DataFrame:
    """
    Validate all cluster CSVs in 'input_dir' against 'subtype_file'.

    Returns:
    summary DataFrame ranked by ARI (also written to ranking_summary.csv)
    """
    os.makedirs(output_dir, exist_ok=True)

    log.info(f"Loading subtype labels from '{subtype_file}' ...")
    label_series = load_label_file(subtype_file, id_length)
    log.info(f"  {len(label_series)} labelled samples, "
             f"subtypes: {sorted(label_series.unique())}")

    prefix     = "gene" if mode == "gene" else "orbit"
    id_col     = "gene_id" if mode == "gene" else "orbit_type"
    pattern    = os.path.join(input_dir, f"{prefix}*_clusters.csv")
    csv_files  = sorted(glob.glob(pattern),
                        key=lambda p: _extract_id(p, prefix))

    log.info(f"Found {len(csv_files)} cluster files.")
    if not csv_files:
        raise FileNotFoundError(
            f"No {prefix}*_clusters.csv files found in {input_dir}."
        )

    summary_rows = []

    for i, fp in enumerate(csv_files, 1):
        item_id = _extract_id(fp, prefix)
        try:
            result = validate_cluster_file(fp, label_series, id_length)
        except Exception as e:
            log.warning(f"  [{i}/{len(csv_files)}] {prefix} {item_id}: ERROR — {e}")
            continue

        if result["status"] != "ok":
            log.debug(f"  {prefix} {item_id}: skipped ({result['status']}, "
                      f"overlap={result['n_overlap']})")
            continue

        if save_jaccard:
            jac_path = os.path.join(output_dir, f"{prefix}{item_id}_jaccard.csv")
            result["jac_df"].to_csv(jac_path)

        summary_rows.append({
            id_col        : item_id,
            "n_overlap"   : result["n_overlap"],
            "ari"         : round(result["ari"], 4),
            "nmi"         : round(result["nmi"], 4),
            "mean_jaccard": round(result["mean_jaccard"], 4),
            "max_jaccard" : round(result["max_jaccard"], 4),
            "pct_correct" : round(result["pct_correct"], 1),
            "best_match"  : "; ".join(
                f"C{c}->{s}({j:.3f})" for c, s, j in result["matched"]
            ),
        })

        log.info(
            f"  [{i:>4}/{len(csv_files)}] {prefix} {item_id}: "
            f"overlap={result['n_overlap']:>4}  ARI={result['ari']:.4f}  "
            f"NMI={result['nmi']:.4f}  mJacc={result['mean_jaccard']:.4f}  "
            f"correct={result['pct_correct']:.1f}%"
        )

    if not summary_rows:
        raise RuntimeError(
            "No items validated successfully.  "
            "Check that sample IDs match between cluster files and the label file."
        )

    summary = (
        pd.DataFrame(summary_rows)
        .sort_values("ari", ascending=False)
        .reset_index(drop=True)
    )
    summary.index += 1

    summary_path = os.path.join(output_dir, "ranking_summary.csv")
    summary.to_csv(summary_path, index_label="rank")
    log.info(f"Ranking summary → {summary_path}")

    log.info(f"\nTop 10 by ARI:\n"
             + summary[[id_col, "n_overlap", "ari", "nmi",
                         "mean_jaccard", "pct_correct"]].head(10).to_string())
    return summary
