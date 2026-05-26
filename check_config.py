#!/usr/bin/env python3
"""
check_config.py  —  Validate your config file before running the pipeline.

Run this after editing your config, before running run_pipeline.py.

    python check_config.py --config my_analysis.yaml

Checks every path, file, and external binary referenced by the config.
Prints a clear PASS / FAIL / WARN for each item with specific fix instructions.
Exits with code 0 only if all required checks pass.

Nothing is modified. Nothing is downloaded. Safe to run repeatedly.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

# ─────────────────────────────────────────────────────────────────────────────
# Terminal colours
# ─────────────────────────────────────────────────────────────────────────────

try:
    import ctypes
    ctypes.windll.kernel32.SetConsoleMode(
        ctypes.windll.kernel32.GetStdHandle(-11), 7)
except Exception:
    pass

GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
BOLD   = "\033[1m"
DIM    = "\033[2m"
RESET  = "\033[0m"

PASS_LABEL = f"{GREEN}PASS{RESET}"
FAIL_LABEL = f"{RED}FAIL{RESET}"
WARN_LABEL = f"{YELLOW}WARN{RESET}"
SKIP_LABEL = f"{DIM}SKIP{RESET}"


# ─────────────────────────────────────────────────────────────────────────────
# Result accumulator
# ─────────────────────────────────────────────────────────────────────────────

class CheckResult:
    def __init__(self):
        self.failures  = []   # hard failures — pipeline cannot run
        self.warnings  = []   # soft warnings — pipeline may still work
        self.passes    = []

    def passed(self, label: str, detail: str = ""):
        tag = f"{PASS_LABEL}  {label}"
        if detail:
            tag += f"  {DIM}({detail}){RESET}"
        print(f"  {tag}")
        self.passes.append(label)

    def failed(self, label: str, problem: str, fix: str):
        print(f"  {FAIL_LABEL}  {label}")
        print(f"          {RED}Problem:{RESET} {problem}")
        for line in fix.strip().splitlines():
            print(f"          {YELLOW}Fix:{RESET}     {line}")
        self.failures.append(label)

    def warned(self, label: str, problem: str, suggestion: str = ""):
        print(f"  {WARN_LABEL}  {label}")
        print(f"          {YELLOW}Note:{RESET} {problem}")
        if suggestion:
            print(f"          Suggestion: {suggestion}")
        self.warnings.append(label)

    def skipped(self, label: str, reason: str):
        print(f"  {SKIP_LABEL}  {label}  {DIM}({reason}){RESET}")

    @property
    def ok(self) -> bool:
        return len(self.failures) == 0


# ─────────────────────────────────────────────────────────────────────────────
# Config helpers
# ─────────────────────────────────────────────────────────────────────────────

def load_config(path: str) -> dict:
    with open(path) as f:
        return yaml.safe_load(f) or {}


def cfg_get(cfg: dict, *keys, default=None):
    node = cfg
    for k in keys:
        if not isinstance(node, dict) or k not in node:
            return default
        node = node[k]
    return node


# ─────────────────────────────────────────────────────────────────────────────
# Individual check functions
# ─────────────────────────────────────────────────────────────────────────────

def check_python_version(r: CheckResult) -> None:
    v = sys.version_info
    label = f"Python version  ({v.major}.{v.minor}.{v.micro})"
    if v >= (3, 10):
        r.passed(label)
    else:
        r.failed(
            label,
            f"Python {v.major}.{v.minor} found; graphlet_stratify requires >= 3.10",
            "Install Python 3.10+ from https://www.python.org/downloads/\n"
            "        or via conda: conda create -n gs python=3.11"
        )


def check_python_packages(r: CheckResult) -> None:
    required = {
        "numpy"          : "numpy",
        "pandas"         : "pandas",
        "scipy"          : "scipy",
        "sklearn"        : "scikit-learn",
        "networkx"       : "networkx",
        "matplotlib"     : "matplotlib",
        "yaml"           : "PyYAML",
        "PyComplexHeatmap": "PyComplexHeatmap",
    }
    missing = []
    for module, pkg in required.items():
        try:
            __import__(module)
        except ImportError:
            missing.append(pkg)

    if not missing:
        r.passed("Python packages", f"{len(required)} required packages installed")
    else:
        r.failed(
            "Python packages",
            f"Missing: {', '.join(missing)}",
            f"pip install {' '.join(missing)}\n"
            "        or:  pip install -r requirements.txt"
        )


def check_java(r: CheckResult) -> None:
    if not shutil.which("java"):
        r.failed(
            "Java (runtime for Jesse)",
            "'java' not found on PATH",
            "Linux:    sudo apt install default-jdk\n"
            "        macOS:    brew install openjdk  (then follow the symlink instruction)\n"
            "        Windows:  https://adoptium.net/\n"
            "        After installing, open a NEW terminal and re-run this check."
        )
        return

    result = subprocess.run(["java", "-version"], capture_output=True, text=True)
    version_line = (result.stderr or result.stdout).split("\n")[0].strip()

    # Extract version number — handles both "1.8.0_xxx" and "11.0.xx" styles
    import re
    m = re.search(r'"(\d+)[\._](\d+)', version_line)
    if m:
        major = int(m.group(1))
        # old-style: "1.8.0" means Java 8
        if major == 1:
            major = int(m.group(2))
        if major >= 8:
            r.passed("Java", version_line)
        else:
            r.failed(
                "Java",
                f"Java {major} found; Jesse requires Java >= 8",
                "Install a newer Java version (11+ recommended)."
            )
    else:
        r.warned("Java", f"Could not parse version from: {version_line}",
                 "Run 'java -version' manually to confirm it is >= 8.")


def check_digger(cfg: dict, r: CheckResult) -> None:
    digger_dir = cfg_get(cfg, "digger", "dir")
    species    = cfg_get(cfg, "digger", "species", default="Homo sapiens[human]")

    if not digger_dir:
        r.failed(
            "DIGGER directory",
            "'digger.dir' is not set in your config",
            "Set it to the root of your DIGGER clone:\n"
            "        digger:\n"
            "          dir: /full/path/to/DIGGER"
        )
        return

    digger_path = Path(digger_dir)
    if not digger_path.exists():
        r.failed(
            "DIGGER directory",
            f"Path does not exist: {digger_path}",
            f"Clone DIGGER:  git clone https://github.com/daisybio/DIGGER\n"
            f"        Then set digger.dir to the cloned folder's absolute path.\n"
            f"        Or run:  python setup_dependencies.py"
        )
        return

    if not digger_path.is_dir():
        r.failed(
            "DIGGER directory",
            f"Path exists but is a file, not a directory: {digger_path}",
            "Set digger.dir to the folder that contains the 'data/' sub-directory."
        )
        return

    # Check the three required data files (downloaded from Zenodo by setup_dependencies.py)
    data_dir = digger_path / "container" / "domain" / "data" / species
    required_files = {
        "PPI.pkl"         : "protein–protein interaction graph",
        "DDI.pkl"         : "domain–domain interaction graph",
        "all_Proteins.csv": "transcript → gene + domain mapping",
    }

    if not data_dir.exists():
        data_parent = digger_path / "container" / "domain" / "data"
        if data_parent.exists():
            available = [d.name for d in data_parent.iterdir() if d.is_dir()]
            hint = (
                f"Available species folders: {available}\n"
                f"        Update 'digger.species' in your config to match one of these."
            )
        else:
            hint = (
                f"Data folder not found: {data_parent}\n"
                f"        Download DIGGER data by running:  python setup_dependencies.py"
            )
        r.failed(
            f"DIGGER data  ({species})",
            f"Species directory not found: {data_dir}",
            hint
        )
        return

    for fname, desc in required_files.items():
        fpath = data_dir / fname
        if fpath.exists():
            size_mb = fpath.stat().st_size / 1_048_576
            r.passed(f"DIGGER data — {fname}", f"{size_mb:.1f} MB  ({desc})")
        else:
            r.failed(
                f"DIGGER data — {fname}",
                f"File not found: {fpath}",
                f"Download DIGGER data files by running:  python setup_dependencies.py\n"
                f"        (downloads the Zenodo archive automatically)"
            )


def check_jesse_jar(cfg: dict, r: CheckResult) -> None:
    jesse_jar = cfg_get(cfg, "jesse", "jar")

    if not jesse_jar:
        r.failed(
            "Jesse JAR",
            "'jesse.jar' is not set in your config",
            "Set it to the full path of the compiled Jesse JAR:\n"
            "        jesse:\n"
            "          jar: /full/path/to/jesse/target/jesse-1.1.0.jar\n"
            "        Or run:  python setup_dependencies.py"
        )
        return

    jar_path = Path(jesse_jar)
    if not jar_path.exists():
        # Try to help: look for any Jesse JAR nearby
        parent = jar_path.parent
        if parent.exists():
            nearby = list(parent.glob("jesse-*.jar"))
            if nearby:
                r.failed(
                    "Jesse JAR",
                    f"JAR not found at: {jar_path}",
                    f"A Jesse JAR was found nearby: {nearby[0]}\n"
                    f"        Update jesse.jar in your config to: {nearby[0]}"
                )
                return

        r.failed(
            "Jesse JAR",
            f"File not found: {jar_path}",
            "Build Jesse with Maven:\n"
            "        cd /path/to/jesse\n"
            "        mvn package -DskipTests\n"
            "        The JAR will appear at target/jesse-*.jar\n"
            "        Or run:  python setup_dependencies.py"
        )
        return

    if not jar_path.suffix == ".jar":
        r.warned("Jesse JAR", f"File does not have .jar extension: {jar_path}",
                 "Make sure jesse.jar points to the compiled JAR file.")
        return

    size_kb = jar_path.stat().st_size / 1024
    r.passed("Jesse JAR", f"{jar_path.name}  ({size_kb:.0f} KB)")

    # Quick smoke-test: try running the JAR with no input (expect it to start, not crash)
    result = subprocess.run(
        ["java", "-jar", str(jar_path)],
        input="",
        capture_output=True, text=True,
        timeout=10,
    )
    # Jesse expects interactive input — it either exits cleanly or waits.
    # A failed JAR load would produce a specific Java error.
    stderr = result.stderr.lower()
    if "error: unable to access jarfile" in stderr:
        r.failed(
            "Jesse JAR (smoke test)",
            "Java cannot read the JAR file",
            "The file may be corrupted. Rebuild Jesse:\n"
            "        cd /path/to/jesse && mvn clean package -DskipTests"
        )
    elif "exception in thread" in stderr and "classnotfound" in stderr:
        r.failed(
            "Jesse JAR (smoke test)",
            "JAR loaded but main class not found — possibly wrong JAR file",
            "Make sure jesse.jar points to the Jesse application JAR,\n"
            "        not a dependency or sources JAR."
        )
    else:
        r.passed("Jesse JAR (smoke test)", "JAR loads correctly")


def check_graphlet_size(cfg: dict, r: CheckResult) -> None:
    gs = cfg_get(cfg, "jesse", "graphlet_size", default=5)
    label = f"jesse.graphlet_size  (={gs})"
    if gs in (3, 4, 5):
        r.passed(label)
    else:
        r.failed(
            label,
            f"graphlet_size must be 3, 4, or 5, got {gs}",
            "Set jesse.graphlet_size to 5 (73 orbits, default), 4 (15 orbits), or 3."
        )


def check_input_matrix(cfg: dict, r: CheckResult) -> None:
    matrix_path = cfg_get(cfg, "input", "expression_matrix")

    if not matrix_path:
        r.failed(
            "input.expression_matrix",
            "Not set in config",
            "Set input.expression_matrix to your transcript × sample TPM matrix file.\n"
            "        Rows = transcripts, columns = patient samples."
        )
        return

    path = Path(matrix_path)
    if not path.exists():
        r.failed(
            "input.expression_matrix",
            f"File does not exist: {path}",
            "Check the path and make sure the file is accessible."
        )
        return

    if not path.is_file():
        r.failed(
            "input.expression_matrix",
            f"Path is a directory, not a file: {path}",
            "input.expression_matrix must point to a FILE (TSV or CSV), not a folder."
        )
        return

    # Spot-check format
    try:
        with open(path, encoding="utf-8") as f:
            header     = f.readline()
            first_data = f.readline()

        sep         = "\t" if "\t" in header else ","
        header_cols = header.strip().split(sep)
        data_cols   = first_data.strip().split(sep) if first_data else []
        n_samples   = len(header_cols) - 1  # first column = transcript IDs

        if n_samples < 1:
            r.failed(
                "input.expression_matrix format",
                "Only 1 column detected — expected transcript_id column + ≥1 sample columns",
                "Format: first column = transcript IDs, remaining columns = samples."
            )
            return

        if data_cols and len(data_cols) >= 2:
            try:
                float(data_cols[1])
                r.passed(
                    "input.expression_matrix",
                    f"{path.name}  |  {n_samples} sample(s)  |  "
                    f"sep={'TAB' if sep == chr(9) else 'COMMA'}"
                )
            except ValueError:
                r.warned(
                    "input.expression_matrix format",
                    f"Second column doesn't look numeric in first data row: '{data_cols[1]}'",
                    "TPM values should be numbers. Check your file format."
                )
        else:
            r.warned("input.expression_matrix format",
                     "File appears to have only a header row — no data rows found.")

    except UnicodeDecodeError:
        r.failed("input.expression_matrix", "File is not valid UTF-8 text",
                 "Expression matrix must be plain text (TSV or CSV).")
    except Exception as e:
        r.warned("input.expression_matrix", f"Could not read file for format check: {e}")


def check_subtype_file(cfg: dict, r: CheckResult) -> None:
    subtype_file = cfg_get(cfg, "validation", "subtype_file")

    if not subtype_file:
        r.skipped("validation.subtype_file",
                  "not set — validation stage will be skipped")
        return

    path = Path(subtype_file)
    if not path.exists():
        r.failed(
            "validation.subtype_file",
            f"File not found: {path}",
            "Set validation.subtype_file to a two-column TSV or CSV file:\n"
            "        Column 1: sample ID   Column 2: subtype label"
        )
        return

    # Spot-check format
    try:
        sep = "\t" if subtype_file.endswith(".tsv") else ","
        import pandas as pd
        df = pd.read_csv(path, sep=sep, header=0, nrows=5)
        if df.shape[1] < 2:
            r.failed(
                "validation.subtype_file",
                f"Only {df.shape[1]} column(s) found — expected at least 2",
                "The subtype file needs two columns: sample_id and subtype_label."
            )
        else:
            r.passed("validation.subtype_file",
                     f"{path.name}  |  columns: {list(df.columns[:2])}")
    except Exception as e:
        r.warned("validation.subtype_file", f"Could not read file: {e}")


def check_clustering_mode(cfg: dict, r: CheckResult) -> None:
    mode = cfg_get(cfg, "clustering", "mode", default="gene")
    label = f"clustering.mode  (={mode})"
    if mode in ("gene", "orbit", "both"):
        r.passed(label)
    else:
        r.failed(
            label,
            f"Invalid value '{mode}'",
            "clustering.mode must be 'gene', 'orbit', or 'both'."
        )

    n_clusters = cfg_get(cfg, "clustering", "n_clusters", default=5)
    if isinstance(n_clusters, int) and n_clusters >= 2:
        r.passed(f"clustering.n_clusters  (={n_clusters})")
    else:
        r.failed(
            f"clustering.n_clusters  (={n_clusters})",
            "Must be an integer >= 2",
            "Set clustering.n_clusters to the number of groups you expect\n"
            "        (e.g. 5 for the 5 PAM50 subtypes)."
        )


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Validate a graphlet_stratify config file before running the pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--config", required=True,
                   help="Path to your YAML config file.")
    p.add_argument("--stages", nargs="+",
                   choices=["filter", "network", "orbits", "cluster", "validate"],
                   default=["filter", "network", "orbits", "cluster", "validate"],
                   help="Which stages to check (default: all).")
    return p.parse_args()


def main():
    args = parse_args()

    # ── Load config ──────────────────────────────────────────────────────────
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"{RED}Config file not found: {config_path}{RESET}")
        print(f"  Copy the template first:  cp demo/config_template.yaml {args.config}")
        sys.exit(1)

    try:
        cfg = load_config(args.config)
    except Exception as e:
        print(f"{RED}Could not parse config YAML: {e}{RESET}")
        print("  Check for indentation errors or missing quotes in your config file.")
        sys.exit(1)

    stages = set(args.stages)
    r = CheckResult()

    print(f"\n{BOLD}graphlet_stratify — config check{RESET}")
    print(f"Config: {config_path.resolve()}")
    print("=" * 60)

    # ── Environment ──────────────────────────────────────────────────────────
    print(f"\n{BOLD}Environment{RESET}")
    check_python_version(r)
    check_python_packages(r)
    if "orbits" in stages:
        check_java(r)

    # ── Input data ───────────────────────────────────────────────────────────
    print(f"\n{BOLD}Input data{RESET}")
    check_input_matrix(cfg, r)

    # ── DIGGER (stages 1 & 2) ────────────────────────────────────────────────
    if stages & {"filter", "network"}:
        print(f"\n{BOLD}DIGGER  (stages: filter, network){RESET}")
        check_digger(cfg, r)

    # ── Jesse (stage 3) ──────────────────────────────────────────────────────
    if "orbits" in stages:
        print(f"\n{BOLD}Jesse  (stage: orbits){RESET}")
        check_jesse_jar(cfg, r)
        check_graphlet_size(cfg, r)

    # ── Clustering ───────────────────────────────────────────────────────────
    if "cluster" in stages:
        print(f"\n{BOLD}Clustering{RESET}")
        check_clustering_mode(cfg, r)

    # ── Validation ───────────────────────────────────────────────────────────
    if "validate" in stages:
        print(f"\n{BOLD}Validation{RESET}")
        check_subtype_file(cfg, r)

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    n_pass = len(r.passes)
    n_fail = len(r.failures)
    n_warn = len(r.warnings)

    print(f"  {GREEN}{n_pass} passed{RESET}   "
          f"{RED}{n_fail} failed{RESET}   "
          f"{YELLOW}{n_warn} warnings{RESET}")

    if r.ok and n_warn == 0:
        print(f"\n{GREEN}{BOLD}All checks passed.{RESET}")
        print("  You can now run the pipeline:")
        print(f"  python run_pipeline.py --config {args.config}\n")
    elif r.ok:
        print(f"\n{YELLOW}{BOLD}Checks passed with warnings.{RESET}")
        print("  Review the warnings above, then run the pipeline:")
        print(f"  python run_pipeline.py --config {args.config}\n")
    else:
        print(f"\n{RED}{BOLD}Config check failed.{RESET}")
        print("  Fix the issues above before running the pipeline.\n")
        print(f"  Failed checks:")
        for f in r.failures:
            print(f"    - {f}")
        print()
        sys.exit(1)


if __name__ == "__main__":
    main()
