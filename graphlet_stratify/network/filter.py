"""

Identify genes that show isoform dominance switching across samples,
and write a whitelist of transcript IDs.

Input format - a single expression matrix file:
  - Rows    = transcripts (Ensembl IDs such as ENST00000366570 or ENST00000366570.3)
  - Columns = patient samples (column names = sample IDs)
  - First column = transcript identifiers; remaining columns = TPM values
  - TSV or CSV (auto-detected from first line)
"""

import logging
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Matrix reader
# ─────────────────────────────────────────────────────────────────────────────

def read_expression_matrix(path: str | Path) -> pd.DataFrame:
    """
    Read a transcript × sample expression matrix.

    Returns a DataFrame with transcript IDs as index and sample IDs as columns.
    Values are expected to be log2(TPM+0.001).
    """
    path = Path(path)
    with open(path, encoding="utf-8") as f:
        first = f.readline()
    sep = "\t" if "\t" in first else ","

    df = pd.read_csv(path, sep=sep, index_col=0, comment="#")
    df.index = df.index.astype(str).str.strip().str.split(".").str[0]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Transcript → gene mapping
# ─────────────────────────────────────────────────────────────────────────────

def load_transcript_to_gene(all_proteins_csv: str | Path) -> dict[str, str]:
    """
    Build a transcript-to-Entrez-gene mapping from DIGGER's all_Proteins.csv.

    Returns dict { "ENST00000366570": "57829" }
    """
    path = Path(all_proteins_csv)
    if not path.exists():
        raise FileNotFoundError(f"all_Proteins.csv not found: {path}")

    log.info(f"Loading transcript→gene mapping from {path} ...")
    df = pd.read_csv(path, low_memory=False)
    df = df[["Transcript stable ID", "NCBI gene ID"]].dropna()
    df["Transcript stable ID"] = df["Transcript stable ID"].str.strip()
    df["NCBI gene ID"] = df["NCBI gene ID"].astype(float).astype(int).astype(str)
    mapping = (
        df.drop_duplicates("Transcript stable ID")
        .set_index("Transcript stable ID")["NCBI gene ID"]
        .to_dict()
    )
    log.info(f"  {len(mapping):,} transcripts mapped")
    return mapping


# ─────────────────────────────────────────────────────────────────────────────
# Dominance tracking
# ─────────────────────────────────────────────────────────────────────────────

def build_dominance_map(
    matrix: pd.DataFrame,
    tx_to_gene: dict[str, str],
    threshold: float,
) -> dict[str, dict[str, str]]:
    """
    For each gene and each sample column, record which transcript was dominant
    (highest TPM at or above threshold).

    Returns
    -------
    dominance[gene_id][sample_id] = transcript_id
    """
    dominance: dict[str, dict[str, str]] = defaultdict(dict)
    n = len(matrix.columns)

    for i, sample_id in enumerate(matrix.columns, 1):
        if i % 100 == 0 or i == n:
            log.info(f"  Processing sample {i}/{n} ...")

        col = matrix[sample_id].dropna()
        gene_best: dict[str, tuple[str, float]] = {}
        for tx, tpm in col.items():
            tpm = float(tpm)
            if tpm < threshold:
                continue
            gene = tx_to_gene.get(str(tx))
            if gene and (gene not in gene_best or tpm > gene_best[gene][1]):
                gene_best[gene] = (str(tx), tpm)

        for gene, (tx, _) in gene_best.items():
            dominance[gene][str(sample_id)] = tx

    return dominance


def compute_switching(
    dominance: dict[str, dict[str, str]],
    min_samples: int,
) -> pd.DataFrame:
    """
    Compute per-gene isoform switching statistics.

    Returns DataFrame sorted by switching_rate descending, with columns:
        entrez_gene_id, n_samples, modal_isoform, modal_count,
        switching_rate, n_competing_isoforms
    """
    records = []
    for gene_id, sample_dom in dominance.items():
        n = len(sample_dom)
        if n < min_samples:
            continue
        counter = Counter(sample_dom.values())
        modal_tx, modal_count = counter.most_common(1)[0]
        records.append({
            "entrez_gene_id":       gene_id,
            "n_samples":            n,
            "modal_isoform":        modal_tx,
            "modal_count":          modal_count,
            "switching_rate":       round(1.0 - modal_count / n, 4),
            "n_competing_isoforms": len(counter),
        })
    df = pd.DataFrame(records).sort_values("switching_rate", ascending=False)
    return df.reset_index(drop=True)


def build_whitelist(
    dominance: dict[str, dict[str, str]],
    retained_genes: set[str],
) -> set[str]:
    """Collect every transcript that ever dominated in a retained gene."""
    whitelist: set[str] = set()
    for gene in retained_genes:
        whitelist.update(dominance[gene].values())
    return whitelist


# ─────────────────────────────────────────────────────────────────────────────
# Optional: write filtered expression matrix
# ─────────────────────────────────────────────────────────────────────────────

def filter_matrix(
    matrix: pd.DataFrame,
    whitelist: set[str],
    output_path: Path,
) -> None:
    """Write a copy of the expression matrix keeping only whitelisted transcript rows."""
    filtered = matrix[matrix.index.isin(whitelist)]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sep = "\t" if str(output_path).endswith(".tsv") else ","
    filtered.to_csv(output_path, sep=sep)
    log.info(f"Filtered matrix ({len(filtered):,} transcripts) → {output_path}")


# ─────────────────────────────────────────────────────────────────────────────
# Run filter
# ─────────────────────────────────────────────────────────────────────────────

def run_filter(
    expression_matrix: str | Path,
    all_proteins_csv: str | Path,
    report_path: str | Path,
    whitelist_path: str | Path,
    filtered_matrix_path: str | Path,
    switch_rate: float = 0.10,
    min_samples: int = 50,
    threshold: float = 0.0,
) -> None:
    """Full isoform switching filter run."""
    tx_to_gene = load_transcript_to_gene(all_proteins_csv)

    log.info(f"Reading expression matrix: {expression_matrix}")
    matrix = read_expression_matrix(expression_matrix)
    log.info(f"  {len(matrix):,} transcripts × {len(matrix.columns):,} samples")

    log.info("Building dominance map ...")
    dominance = build_dominance_map(matrix, tx_to_gene, threshold)
    log.info(f"  Genes observed: {len(dominance):,}")

    log.info("Computing switching rates ...")
    stats_df = compute_switching(dominance, min_samples)
    log.info(f"  Genes passing min_samples >= {min_samples}: {len(stats_df):,}")

    retained = stats_df[stats_df["switching_rate"] > switch_rate]
    log.info(f"  Genes retained (switch_rate > {switch_rate}): {len(retained):,}")

    report_path = Path(report_path)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(report_path, sep="\t", index=False)
    log.info(f"Switching report → {report_path}")

    retained_genes = set(retained["entrez_gene_id"])
    whitelist = build_whitelist(dominance, retained_genes)

    whitelist_path = Path(whitelist_path)
    whitelist_path.parent.mkdir(parents=True, exist_ok=True)
    whitelist_path.write_text("\n".join(sorted(whitelist)))
    log.info(f"Whitelist ({len(whitelist):,} transcripts) → {whitelist_path}")

    log.info("Writing filtered expression matrix ...")
    filter_matrix(matrix, whitelist, Path(filtered_matrix_path))
