#!/usr/bin/env python3

"""
─────────────────────────────────────────────────────────────────────────────

The script runs all or selected stages of the patient-stratification pipeline:

Stage 1 (filter):   Identify genes where the dominant isoform varies across samples, 
                    write whitelist.

Stage 2 (network):  Build one isoform-aware DDI-filtered PPI network per sample.

Stage 3 (orbits):   Run Jesse to compute graphlet orbit frequencies for every network.

Stage 4 (cluster):  Hierarchical clustering of samples, either per-gene or 
                    per-orbit-type (or both).

Stage 5 (validate): Score clusters against known subtype labels. 
                    Requires: a two-column label file.

─────────────────────────────────────────────────────────────────────────────

RUN THE FULL PIPELINE WITH A CONFIG FILE

  cp demo/config_template.yaml my_config.yaml
  # edit my_config.yaml to match your paths
  python run_pipeline.py --config my_config.yaml

─────────────────────────────────────────────────────────────────────────────
    
RUN INDIVIDUAL STAGES

  python run_pipeline.py --config my_config.yaml --stages filter network
  python run_pipeline.py --config my_config.yaml --stages orbits
  python run_pipeline.py --config my_config.yaml --stages cluster validate

─────────────────────────────────────────────────────────────────────────────

"""

import argparse
import logging
import sys
from pathlib import Path

import yaml

# Logging for tracking and debugging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("graphlet_stratify")

ALL_STAGES = ["filter", "network", "orbits", "cluster", "validate"]


# Loading configuration
def load_config(config_path: str) -> dict:
    with open(config_path, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    return cfg


def cfg_get(cfg: dict, *keys, default=None):
    """Safely navigate nested config keys."""
    node = cfg
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            return default
        node = node[k]
    return node


# ─────────────────────────────────────────────────────────────────────────────
# Stage runners
# ─────────────────────────────────────────────────────────────────────────────

def stage_filter(cfg: dict) -> None:
    """Stage 1: isoform switching filter."""
    log.info("-" * 60)
    log.info("STAGE 1 — Isoform switching filter")
    log.info("-" * 60)

    from graphlet_stratify.network.filter import run_filter

    digger_dir   = cfg_get(cfg, "digger", "dir")
    species      = cfg_get(cfg, "digger", "species", default="Homo sapiens[human]")
    all_proteins = (
        Path(digger_dir) / "container" / "domain" / "data" / species / "all_Proteins.csv"
    )

    run_filter(
        expression_matrix    = cfg_get(cfg, "input", "expression_matrix"),
        all_proteins_csv     = all_proteins,
        report_path          = cfg_get(cfg, "filter", "report_path",
                                       default="results/switching_report.tsv"),
        whitelist_path       = cfg_get(cfg, "filter", "whitelist_path",
                                       default="results/whitelist.txt"),
        switch_rate          = cfg_get(cfg, "filter", "switch_rate", default=0.10),
        min_samples          = cfg_get(cfg, "filter", "min_samples", default=50),
        threshold            = cfg_get(cfg, "input", "tpm_threshold", default=0.0),
        filtered_matrix_path = cfg_get(cfg, "filter", "filtered_matrix_path",
                                       default="results/expression_filtered.tsv"),
    )


def stage_network(cfg: dict) -> None:
    """Stage 2: build per-sample isoform-aware networks."""
    log.info("-" * 60)
    log.info("STAGE 2 — Network building")
    log.info("-" * 60)

    from graphlet_stratify.network.build import run_batch

    filtered_path = cfg_get(cfg, "filter", "filtered_matrix_path")
    if filtered_path and Path(filtered_path).exists():
        expression_matrix = filtered_path
        log.info(f"Using filtered matrix: {filtered_path}")
    else:
        expression_matrix = cfg_get(cfg, "input", "expression_matrix")
        log.warning(
            "filtered_matrix_path not set or file not found — "
            "using full unfiltered matrix."
            "Run stage 'filter' first, or set filter.filtered_matrix_path in your config."
        )

    run_batch(
        expression_matrix = expression_matrix,
        output_dir        = cfg_get(cfg, "network", "output_dir",
                                    default="results/networks"),
        digger_dir        = cfg_get(cfg, "digger", "dir"),
        species           = cfg_get(cfg, "digger", "species",
                                    default="Homo sapiens[human]"),
        threshold         = cfg_get(cfg, "input", "tpm_threshold", default=0.0),
    )


def stage_orbits(cfg: dict) -> None:
    """Stage 3: run Jesse on all network files."""
    log.info("-" * 60)
    log.info("STAGE 3 — Graphlet orbit counting (Jesse)")
    log.info("-" * 60)

    from graphlet_stratify.orbits.jesse import run_jesse

    run_jesse(
        input_dir     = cfg_get(cfg, "network", "output_dir",
                                default="results/networks"),
        output_dir    = cfg_get(cfg, "orbits", "output_dir",
                                default="results/orbits"),
        jesse_jar     = cfg_get(cfg, "jesse", "jar"),
        graphlet_size = cfg_get(cfg, "jesse", "graphlet_size", default=5),
        file_pattern  = cfg_get(cfg, "orbits", "network_file_pattern",
                                default="*.txt"),
    )


def stage_cluster(cfg: dict) -> None:
    """Stage 4: hierarchical clustering."""
    log.info("-" * 60)
    log.info("STAGE 4 — Clustering")
    log.info("-" * 60)

    mode = cfg_get(cfg, "clustering", "mode", default="gene")
    if mode not in ("gene", "orbit", "both"):
        raise ValueError(f"clustering.mode must be 'gene', 'orbit', or 'both', got '{mode}'")

    n_clusters   = cfg_get(cfg, "clustering", "n_clusters", default=5)
    min_presence = cfg_get(cfg, "clustering", "min_presence", default=0.90)
    outlier_z    = cfg_get(cfg, "clustering", "outlier_z", default=5.0)
    save_pca     = cfg_get(cfg, "clustering", "save_pca", default=False)
    id_length    = cfg_get(cfg, "clustering", "id_length", default=None)
    orbits_dir   = cfg_get(cfg, "orbits", "output_dir", default="results/orbits")

    # Load optional subtype map for PCA colouring
    subtype_map = None
    label_file = cfg_get(cfg, "validation", "subtype_file")
    if save_pca and label_file:
        from graphlet_stratify.validation.validate import load_label_file
        label_series = load_label_file(label_file, id_length)
        subtype_map  = label_series.to_dict()

    if mode in ("gene", "both"):
        from graphlet_stratify.clustering.per_gene import run_per_gene_clustering
        run_per_gene_clustering(
            input_dir    = orbits_dir,
            output_dir   = cfg_get(cfg, "clustering", "gene_output_dir",
                                   default="results/clusters/per_gene"),
            n_clusters   = n_clusters,
            min_presence = min_presence,
            outlier_z    = outlier_z,
            save_pca     = save_pca,
            subtype_map  = subtype_map,
            id_length    = id_length,
        )

    if mode in ("orbit", "both"):
        from graphlet_stratify.clustering.per_orbit import run_per_orbit_clustering
        run_per_orbit_clustering(
            input_dir    = orbits_dir,
            output_dir   = cfg_get(cfg, "clustering", "orbit_output_dir",
                                   default="results/clusters/per_orbit"),
            n_clusters   = n_clusters,
            min_presence = min_presence,
            outlier_z    = outlier_z,
            save_pca     = save_pca,
            subtype_map  = subtype_map,
            id_length    = id_length,
        )


def stage_validate(cfg: dict) -> None:
    """Stage 5: validate clusters against subtype labels."""
    log.info("-" * 60)
    log.info("STAGE 5 — Validation")
    log.info("-" * 60)

    from graphlet_stratify.validation.validate import run_validation

    subtype_file = cfg_get(cfg, "validation", "subtype_file")
    if not subtype_file:
        log.warning("validation.subtype_file not set in config — skipping validation.")
        return

    mode      = cfg_get(cfg, "clustering", "mode", default="gene")
    id_length = cfg_get(cfg, "clustering", "id_length", default=None)

    if mode in ("gene", "both"):
        run_validation(
            input_dir    = cfg_get(cfg, "clustering", "gene_output_dir",
                                   default="results/clusters/per_gene"),
            output_dir   = cfg_get(cfg, "validation", "gene_output_dir",
                                   default="results/validation/per_gene"),
            subtype_file = subtype_file,
            mode         = "gene",
            save_jaccard = cfg_get(cfg, "validation", "save_jaccard", default=False),
            id_length    = id_length,
        )

    if mode in ("orbit", "both"):
        run_validation(
            input_dir    = cfg_get(cfg, "clustering", "orbit_output_dir",
                                   default="results/clusters/per_orbit"),
            output_dir   = cfg_get(cfg, "validation", "orbit_output_dir",
                                   default="results/validation/per_orbit"),
            subtype_file = subtype_file,
            mode         = "orbit",
            save_jaccard = cfg_get(cfg, "validation", "save_jaccard", default=False),
            id_length    = id_length,
        )


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline manager
# ─────────────────────────────────────────────────────────────────────────────

STAGE_FN = {
    "filter"  : stage_filter,
    "network" : stage_network,
    "orbits"  : stage_orbits,
    "cluster" : stage_cluster,
    "validate": stage_validate,
}


_TEMPLATE = Path(__file__).parent / "demo" / "config_template.yaml"


def parse_args():
    parser = argparse.ArgumentParser(
        description="graphlet_stratify — unified patient stratification pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Simple mode (no configuration file)
    parser.add_argument(
        "--matrix", metavar="FILE",
        help=(
            "Expression matrix: rows = transcripts, columns = samples, values = TPM.  "
            "TSV or CSV.  When given, --config is not needed."
        ),
    )
    parser.add_argument(
        "--subtypes", metavar="FILE",
        help="Two-column label file (sample_id, subtype).  TSV or CSV.",
    )
    parser.add_argument(
        "--output", metavar="DIR", default="results",
        help="Root folder for all output files (default: results/).",
    )

    # Using configuration file for parameter tuning
    parser.add_argument(
        "--config", metavar="FILE",
        help=(
            "YAML config file for full control over all settings.  "
            "Not needed when --matrix is given."
        ),
    )
    parser.add_argument(
        "--stages", nargs="+", choices=ALL_STAGES, default=ALL_STAGES,
        metavar="STAGE",
        help=f"Which stages to run (default: all).  Choices: {', '.join(ALL_STAGES)}",
    )
    return parser.parse_args()


def _build_config_from_flags(args) -> dict:
    """Build a config dict from CLI flags, using the template for DIGGER/Jesse paths."""
    if _TEMPLATE.exists():
        with open(_TEMPLATE, encoding="utf-8") as f:
            cfg = yaml.safe_load(f) or {}
    else:
        cfg = {}

    out = Path(args.output)

    cfg.setdefault("input", {})["expression_matrix"] = str(Path(args.matrix).resolve())
    filtered_matrix = str(out / "expression_filtered.tsv")
    cfg.setdefault("filter", {}).update({
        "report_path":         str(out / "switching_report.tsv"),
        "whitelist_path":      str(out / "whitelist.txt"),
        "filtered_matrix_path": filtered_matrix,
    })
    cfg.setdefault("network",  {})["output_dir"] = str(out / "networks")
    cfg.setdefault("orbits",   {})["output_dir"] = str(out / "orbits")
    cfg.setdefault("clustering", {}).update({
        "gene_output_dir":  str(out / "clusters" / "per_gene"),
        "orbit_output_dir": str(out / "clusters" / "per_orbit"),
    })
    cfg.setdefault("validation", {}).update({
        "subtype_file":     str(Path(args.subtypes).resolve()) if args.subtypes else None,
        "gene_output_dir":  str(out / "validation" / "per_gene"),
        "orbit_output_dir": str(out / "validation" / "per_orbit"),
    })
    return cfg


def preflight_check(cfg: dict, stages: list) -> None:
    """
    Path checks that run before any stage starts.

    Catches missing paths, unset config keys, and Java not on PATH immediately.
    """
    import shutil
    errors = []
    stages_set = set(stages)

    # Expression matrix — always required
    matrix_path = cfg_get(cfg, "input", "expression_matrix")
    if not matrix_path:
        errors.append(
            "input.expression_matrix is not set in your config.\n"
            "  Set it to your transcript × sample TPM matrix file."
        )
    elif not Path(matrix_path).is_file():
        errors.append(
            f"input.expression_matrix does not exist or is not a file:\n"
            f"  {matrix_path}"
        )

    # DIGGER — needed for filter and network stages
    if stages_set & {"filter", "network"}:
        digger_dir = cfg_get(cfg, "digger", "dir")
        species    = cfg_get(cfg, "digger", "species", default="Homo sapiens[human]")
        if not digger_dir:
            errors.append(
                "digger.dir is not set in your config.\n"
            )
        else:
            data_dir = Path(digger_dir) / "container" / "domain" / "data" / species
            for fname in ("PPI.pkl", "DDI.pkl", "all_Proteins.csv"):
                fpath = data_dir / fname
                if not fpath.exists():
                    errors.append(
                        f"DIGGER data file not found: {fpath}\n"
                        f"  Run 'python setup_dependencies.py' to download DIGGER data."
                    )

    # Jesse JAR and Java — needed for orbits stage
    if "orbits" in stages_set:
        jesse_jar = cfg_get(cfg, "jesse", "jar")
        if not jesse_jar:
            errors.append(
                "jesse.jar is not set in your config.\n"
            )
        elif not Path(jesse_jar).exists():
            errors.append(
                f"Jesse JAR not found: {jesse_jar}\n"
            )
        if not shutil.which("java"):
            errors.append(
                "'java' is not on your PATH — required to run Jesse.\n"
            )

    # Subtype file — needed for validate stage
    if "validate" in stages_set:
        subtype_file = cfg_get(cfg, "validation", "subtype_file")
        if subtype_file and not Path(subtype_file).exists():
            errors.append(
                f"validation.subtype_file not found: {subtype_file}"
            )

    if errors:
        log.error("Pre-processing check failed — fix these issues before running analysis:")
        for i, err in enumerate(errors, 1):
            log.error(f"  [{i}] {err}")
        sys.exit(1)

    log.info("Pre-processing check passed.")


def main():
    args = parse_args()

    if args.matrix:
        cfg = _build_config_from_flags(args)
    elif args.config:
        cfg = load_config(args.config)
    else:
        log.error(
            "Provide either --matrix <file> or --config <file>.\n"
        )
        sys.exit(1)

    stages = args.stages

    preflight_check(cfg, stages)

    for stage in ALL_STAGES:
        if stage not in stages:
            continue
        try:
            STAGE_FN[stage](cfg)
        except Exception as e:
            log.error(f"Stage '{stage}' failed: {e}")
            sys.exit(1)

    log.info("Pipeline complete.")


if __name__ == "__main__":
    main()
