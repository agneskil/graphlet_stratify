"""
Build one isoform-aware, DDI-filtered PPI network per sample.


This module reads three files from the DIGGER data archive (Zenodo):

    <digger_dir>/container/domain/data/<species>/PPI.pkl
    <digger_dir>/container/domain/data/<species>/DDI.pkl
    <digger_dir>/container/domain/data/<species>/all_Proteins.csv

You must have cloned / installed DIGGER. See README.md "External dependencies" 
for setup instructions.

"""

import itertools
import logging
import pickle
from collections import defaultdict
from pathlib import Path

import networkx as nx
import pandas as pd

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# DIGGER data loading
# ─────────────────────────────────────────────────────────────────────────────

def _digger_paths(digger_dir: str | Path, species: str) -> tuple[Path, Path, Path]:
    """
    Return (PPI.pkl, DDI.pkl, all_Proteins.csv) paths for the given species.
    """
    data = Path(digger_dir) / "container" / "domain" / "data" / species
    return data / "PPI.pkl", data / "DDI.pkl", data / "all_Proteins.csv"


def load_digger_data(
    digger_dir: str | Path,
    species: str = "Homo sapiens[human]",
) -> tuple[dict, nx.Graph, dict, nx.Graph]:
    """
    Load DIGGER reference data needed to build isoform-aware networks.

    """
    ppi_pkl, ddi_pkl, all_proteins = _digger_paths(digger_dir, species)

    for p in (ppi_pkl, ddi_pkl, all_proteins):
        if not p.exists():
            raise FileNotFoundError(
                f"DIGGER data file not found: {p}\n"
                f"Make sure --digger_dir points to the DIGGER root folder and "
                f"--species matches the sub-folder name under data/."
            )

    proteins_df = pd.read_csv(all_proteins, low_memory=False)
    proteins_df = proteins_df[
        ["Transcript stable ID", "NCBI gene ID", "Pfam ID"]
    ].dropna(subset=["Transcript stable ID", "NCBI gene ID"])
    proteins_df["NCBI gene ID"] = (
        proteins_df["NCBI gene ID"].astype(float).astype(int).astype(str)
    )
    proteins_df["Transcript stable ID"] = proteins_df["Transcript stable ID"].str.strip()

    transcript_to_entrez: dict[str, str] = (
        proteins_df.drop_duplicates("Transcript stable ID")
        .set_index("Transcript stable ID")["NCBI gene ID"]
        .to_dict()
    )

    t2d: dict[str, set] = defaultdict(set)
    for _, row in proteins_df.dropna(subset=["Pfam ID"]).iterrows():
        t2d[row["Transcript stable ID"]].add(row["Pfam ID"])
    transcript_to_domains = dict(t2d)

    with open(ppi_pkl, "rb") as f:
        ppi_graph: nx.Graph = pickle.load(f)

    with open(ddi_pkl, "rb") as f:
        ddi_graph: nx.Graph = pickle.load(f)

    return transcript_to_entrez, ppi_graph, transcript_to_domains, ddi_graph


# ─────────────────────────────────────────────────────────────────────────────
# Per-sample helpers
# ─────────────────────────────────────────────────────────────────────────────

def read_expression_matrix(path: str | Path) -> pd.DataFrame:
    """
    Read a transcript × sample expression matrix.

    """
    path = Path(path)
    with open(path, encoding="utf-8") as f:
        first = f.readline()
    sep = "\t" if "\t" in first else ","

    df = pd.read_csv(path, sep=sep, index_col=0, comment="#")
    df.index = df.index.astype(str).str.strip().str.split(".").str[0]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def select_dominant_isoforms(
    transcript_tpm: dict[str, float],
    transcript_to_entrez: dict[str, str],
) -> dict[str, str]:
    """
    For each gene, pick the single transcript with the highest TPM.

    Returns:
    dict { entrez_id: dominant_transcript_id }
    """
    gene_to_transcripts: dict[str, list[tuple[str, float]]] = defaultdict(list)
    for tx, tpm in transcript_tpm.items():
        gene = transcript_to_entrez.get(tx)
        if gene:
            gene_to_transcripts[gene].append((tx, tpm))

    return {
        gene: max(txs, key=lambda x: x[1])[0]
        for gene, txs in gene_to_transcripts.items()
    }


def _dominant_isoforms_have_ddi(
    tx_a: str,
    tx_b: str,
    transcript_to_domains: dict,
    ddi_graph: nx.Graph,
) -> bool:
    """Return True if any domain pair between tx_a and tx_b has DDI support."""
    domains_a = transcript_to_domains.get(tx_a, set())
    domains_b = transcript_to_domains.get(tx_b, set())
    if not domains_a or not domains_b:
        return False
    return any(ddi_graph.has_edge(da, db) for da in domains_a for db in domains_b)


def build_sample_network(
    transcript_tpm: dict[str, float],
    transcript_to_entrez: dict[str, str],
    ppi_graph: nx.Graph,
    transcript_to_domains: dict,
    ddi_graph: nx.Graph,
) -> pd.DataFrame:
    """
    Build an isoform-aware, DDI-filtered PPI network for one sample.

    An edge is retained only if both genes are expressed AND their dominant
    isoforms carry at least one domain pair with known DDI support.

    Returns:
    DataFrame with columns: Protein 1, Protein 2, weight,
                            dominant_isoform_1, dominant_isoform_2
    """
    dominant = select_dominant_isoforms(transcript_tpm, transcript_to_entrez)
    expressed_genes = list(dominant.keys())
    rows = []

    for gene_a, gene_b in itertools.combinations(expressed_genes, 2):
        if not ppi_graph.has_edge(gene_a, gene_b):
            continue
        tx_a, tx_b = dominant[gene_a], dominant[gene_b]
        if _dominant_isoforms_have_ddi(tx_a, tx_b, transcript_to_domains, ddi_graph):
            rows.append({
                "Protein 1": min(gene_a, gene_b),
                "Protein 2": max(gene_a, gene_b),
                "weight": 1.0,
                "dominant_isoform_1": tx_a,
                "dominant_isoform_2": tx_b,
            })

    if not rows:
        return pd.DataFrame(columns=["Protein 1", "Protein 2", "weight",
                                     "dominant_isoform_1", "dominant_isoform_2"])
    df = pd.DataFrame(rows).drop_duplicates(subset=["Protein 1", "Protein 2"])
    df = df.iloc[df[["Protein 1", "Protein 2"]]
                 .apply(lambda r: (int(r["Protein 1"]), int(r["Protein 2"])), axis=1)
                 .argsort()]
    return df.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# Batch runner
# ─────────────────────────────────────────────────────────────────────────────

def run_batch(
    expression_matrix: str | Path,
    output_dir: str | Path,
    digger_dir: str | Path,
    species: str = "Homo sapiens[human]",
    threshold: float = 0.0,
) -> None:
    """
    Build one network TSV per sample column in the expression matrix.

    """
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    transcript_to_entrez, ppi_graph, transcript_to_domains, ddi_graph = \
        load_digger_data(digger_dir, species)

    log.info(f"Reading expression matrix: {expression_matrix}")
    matrix = read_expression_matrix(expression_matrix)
    n_samples = len(matrix.columns)
    log.info(f"  {len(matrix):,} transcripts × {n_samples:,} samples | log2 threshold: {threshold}")


    total_edges = 0
    for i, sample_id in enumerate(matrix.columns, 1):
        log.info(f"[{i}/{n_samples}] {sample_id}")
        col = matrix[sample_id].dropna()
        transcript_tpm = {
            str(tx): float(tpm)
            for tx, tpm in col.items()
            if float(tpm) >= threshold
        }

        network_df = build_sample_network(
            transcript_tpm, transcript_to_entrez,
            ppi_graph, transcript_to_domains, ddi_graph,
        )
        out_file = out_path / (str(sample_id) + "_network.txt")
        network_df[["Protein 1", "Protein 2"]].to_csv(
            out_file, sep="\t", index=False, header=False
        )
        log.info(f"  {len(network_df):,} edges → {out_file.name}")

    log.info(f"Done. Networks for {n_samples} samples written to {out_path}.")
