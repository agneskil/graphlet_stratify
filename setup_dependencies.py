#!/usr/bin/env python3
"""
setup_dependencies.py  —  Download and set up external dependencies.

Run this ONCE before using graphlet_stratify for the first time.

    python setup_dependencies.py

What this script does automatically
-------------------------------------
  1. Checks that Java (>= 8) is on your PATH.
  2. Checks that Git is on your PATH.
  3. Clones DIGGER into  ./dependencies/DIGGER/
  4. Clones Jesse into   ./dependencies/jesse/
  5. Builds Jesse (runs `mvn package` inside the Jesse directory).
     Maven must be installed — the script checks for it and gives
     installation instructions if it is missing.
  6. Finds the Jesse JAR and writes the correct paths into
     demo/config_template.yaml so you don't have to fill them in manually.

What you still need to do yourself
-------------------------------------
  - Install Java (>= 8)  if not already installed.
  - Install Maven         if not already installed.
  - Install Git           if not already installed.

The script prints precise installation instructions for any missing tool.

After this script succeeds
---------------------------
  cp demo/config_template.yaml my_analysis.yaml
  # Edit my_analysis.yaml: set input.sample_dir and validation.subtype_file
  python run_pipeline.py --config my_analysis.yaml
"""

import os
import shutil
import subprocess
import sys
import glob
import urllib.request
import zipfile
from pathlib import Path


def _run(cmd: list, **kwargs) -> subprocess.CompletedProcess:
    """Run a command, resolving the executable via shutil.which first.

    On Windows, commands like 'mvn' are actually 'mvn.cmd' batch files.
    subprocess.run won't find them without the full resolved path.
    """
    resolved = shutil.which(cmd[0])
    if resolved:
        cmd = [resolved] + cmd[1:]
    return subprocess.run(cmd, **kwargs)

DIGGER_URL = "https://github.com/daisybio/DIGGER"
JESSE_URL  = "https://github.com/biointec/jesse"

DEPS_DIR   = Path(__file__).parent / "dependencies"
DIGGER_DIR = DEPS_DIR / "DIGGER"
JESSE_DIR  = DEPS_DIR / "jesse"


# ─────────────────────────────────────────────────────────────────────────────
# Terminal colours (graceful fallback on Windows)
# ─────────────────────────────────────────────────────────────────────────────

try:
    import ctypes
    ctypes.windll.kernel32.SetConsoleMode(ctypes.windll.kernel32.GetStdHandle(-11), 7)
except Exception:
    pass

GREEN  = "\033[32m"
YELLOW = "\033[33m"
RED    = "\033[31m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

def ok(msg):    print(f"  {GREEN}✓{RESET}  {msg}")
def warn(msg):  print(f"  {YELLOW}!{RESET}  {msg}")
def fail(msg):  print(f"  {RED}✗{RESET}  {msg}")
def section(msg): print(f"\n{BOLD}{msg}{RESET}")


# ─────────────────────────────────────────────────────────────────────────────
# Prerequisite checks
# ─────────────────────────────────────────────────────────────────────────────

def check_git() -> bool:
    if shutil.which("git"):
        ok("Git found.")
        return True
    fail("Git not found on PATH.")
    print("""
  Git is required to clone DIGGER and Jesse.

  Install it:
    Linux (Debian/Ubuntu):  sudo apt install git
    Linux (Fedora/RHEL):    sudo dnf install git
    macOS:                  brew install git   (or: xcode-select --install)
    Windows:                https://git-scm.com/download/win
""")
    return False


def check_java() -> bool:
    """Check that Java >= 8 is available."""
    if not shutil.which("java"):
        fail("Java not found on PATH.")
        print("""
  Java (>= 8) is required to run Jesse.

  Install it:
    Linux (Debian/Ubuntu):  sudo apt install default-jdk
    Linux (Fedora/RHEL):    sudo dnf install java-latest-openjdk
    macOS:                  brew install openjdk
                            (follow the 'For the system Java wrappers' symlink
                             instruction that Homebrew prints)
    Windows:                https://adoptium.net/
                            Download and run the MSI installer.

  After installing, open a NEW terminal and re-run this script.
""")
        return False

    result = _run(["java", "-version"], capture_output=True, text=True)
    version_line = (result.stderr or result.stdout).split("\n")[0]
    ok(f"Java found: {version_line.strip()}")
    return True


def check_maven() -> bool:
    """Check that Maven is available (needed to build Jesse)."""
    if shutil.which("mvn"):
        result = _run(["mvn", "--version"], capture_output=True, text=True)
        ok(f"Maven found: {result.stdout.split(chr(10))[0].strip()}")
        return True

    fail("Maven not found on PATH.")
    print("""
  Maven is required to build Jesse from source.

  Install it:
    Linux (Debian/Ubuntu):  sudo apt install maven
    Linux (Fedora/RHEL):    sudo dnf install maven
    macOS:                  brew install maven
    Windows:
      1. Download the binary ZIP from https://maven.apache.org/download.cgi
      2. Extract it (e.g. to C:\\Program Files\\apache-maven-3.x.x)
      3. Add  C:\\Program Files\\apache-maven-3.x.x\\bin  to your PATH:
           System Properties → Environment Variables → Path → New

  After installing, open a NEW terminal and re-run this script.
""")
    return False


# ─────────────────────────────────────────────────────────────────────────────
# Clone helpers
# ─────────────────────────────────────────────────────────────────────────────

def git_clone(url: str, dest: Path, name: str) -> bool:
    if dest.exists():
        warn(f"{name} already cloned at {dest}  (skipping clone, will update).")
        result = subprocess.run(["git", "-C", str(dest), "pull", "--ff-only"],
                                capture_output=True, text=True)
        if result.returncode == 0:
            ok(f"{name} updated.")
        else:
            warn(f"Could not update {name} (this is fine if you have local changes).")
        return True

    print(f"  Cloning {name} from {url} ...")
    result = subprocess.run(
        ["git", "clone", "--depth", "1", url, str(dest)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        fail(f"Git clone failed for {name}:")
        print(f"    {result.stderr.strip()}")
        return False
    ok(f"{name} cloned to {dest}")
    return True


# ─────────────────────────────────────────────────────────────────────────────
# DIGGER data download
# ─────────────────────────────────────────────────────────────────────────────

DIGGER_DATA_URL = "https://zenodo.org/records/14770005/files/data.zip"
DIGGER_DATA_DIR = DIGGER_DIR / "container" / "domain" / "data"
DIGGER_SPECIES  = "Homo sapiens[human]"


def _show_progress(block_num: int, block_size: int, total_size: int) -> None:
    downloaded = block_num * block_size
    if total_size > 0:
        pct = min(downloaded / total_size * 100, 100)
        mb_done  = downloaded   / 1024 / 1024
        mb_total = total_size   / 1024 / 1024
        print(f"\r    {pct:5.1f}%  ({mb_done:.0f} / {mb_total:.0f} MB)", end="", flush=True)


def download_digger_data(digger_dir: Path) -> bool:
    """
    Download DIGGER data files from Zenodo and extract into
    <digger_dir>/container/domain/data/.
    Skips the download if the required files are already present.
    """
    species_dir   = digger_dir / "container" / "domain" / "data" / DIGGER_SPECIES
    required      = ["PPI.pkl", "DDI.pkl", "all_Proteins.csv"]
    already_there = [f for f in required if (species_dir / f).exists()]

    if len(already_there) == len(required):
        ok("DIGGER data files already present — skipping download.")
        return True

    print(f"  Downloading DIGGER data from Zenodo:")
    print(f"    {DIGGER_DATA_URL}")
    print(f"  This file is ~200 MB — may take a few minutes on a slow connection.")

    zip_path = digger_dir / "digger_data.zip"
    try:
        urllib.request.urlretrieve(DIGGER_DATA_URL, zip_path, _show_progress)
        print()  # newline after progress bar
    except Exception as exc:
        fail(f"Download failed: {exc}")
        return False

    print(f"  Extracting {zip_path.name} ...")
    extract_to = digger_dir / "container" / "domain"
    extract_to.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            zf.extractall(extract_to)
        zip_path.unlink()
    except Exception as exc:
        fail(f"Extraction failed: {exc}")
        return False

    missing = [f for f in required if not (species_dir / f).exists()]
    if missing:
        fail(f"Expected files still missing after extraction: {missing}")
        print(f"  Extracted to: {extract_to}")
        print(f"  Check that the zip contained a 'data/{DIGGER_SPECIES}/' folder.")
        return False

    ok("DIGGER data files downloaded and verified.")
    return True


# ─────────────────────────────────────────────────────────────────────────────
# Jesse build
# ─────────────────────────────────────────────────────────────────────────────

def build_jesse(jesse_dir: Path) -> Path | None:
    """
    Run `mvn package -DskipTests` inside the Jesse directory.
    Returns the path to the built JAR, or None on failure.
    """
    # Check whether a JAR already exists (skip rebuild)
    existing_jars = glob.glob(str(jesse_dir / "target" / "jesse-*.jar"))
    existing_jars = [j for j in existing_jars if "sources" not in j]
    if existing_jars:
        jar_path = Path(existing_jars[0])
        ok(f"Jesse JAR already built: {jar_path.name}  (skipping rebuild)")
        return jar_path

    print("  Building Jesse with Maven (this may take 1-2 minutes) ...")
    result = _run(
        ["mvn", "package", "-DskipTests", "--no-transfer-progress"],
        cwd=jesse_dir,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        fail("Maven build failed.")
        print(f"\n  Maven output (last 30 lines):\n")
        for line in result.stdout.split("\n")[-30:]:
            print(f"    {line}")
        print(f"\n  If you see a network error, check your internet connection.")
        print(f"  If you see a Java version error, install Java 11 or higher.")
        return None

    jars = glob.glob(str(jesse_dir / "target" / "jesse-*.jar"))
    jars = [j for j in jars if "sources" not in j]
    if not jars:
        fail("Maven succeeded but no JAR found in target/.")
        return None

    jar_path = Path(jars[0])
    ok(f"Jesse built: {jar_path}")
    return jar_path


# ─────────────────────────────────────────────────────────────────────────────
# Config patching
# ─────────────────────────────────────────────────────────────────────────────

def patch_config(digger_dir: Path, jesse_jar: Path) -> None:
    """
    Write the resolved DIGGER and Jesse paths into the config template
    so users don't have to find and paste them manually.
    """
    config_path = Path(__file__).parent / "demo" / "config_template.yaml"
    if not config_path.exists():
        warn("demo/config_template.yaml not found — skipping config patch.")
        return

    text = config_path.read_text()

    # Replace placeholder values
    text = text.replace(
        "  dir: /path/to/DIGGER",
        f"  dir: {digger_dir.resolve()}"
    )
    text = text.replace(
        "  jar: /path/to/jesse/target/jesse-1.1.0.jar",
        f"  jar: {jesse_jar.resolve()}"
    )

    config_path.write_text(text)
    ok(f"Config template updated with correct paths.")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print(f"\n{BOLD}graphlet_stratify — dependency setup{RESET}")
    print("=" * 50)

    # ── Step 1: check prerequisites ──────────────────────────────────────────
    section("Step 1/4  Checking prerequisites")
    git_ok  = check_git()
    java_ok = check_java()
    mvn_ok  = check_maven()

    missing = []
    if not git_ok:  missing.append("git")
    if not java_ok: missing.append("java")
    if not mvn_ok:  missing.append("maven")

    if missing:
        print(f"\n{RED}{BOLD}Setup cannot continue.{RESET}")
        print(f"  Please install the following and re-run this script:")
        for m in missing:
            print(f"    - {m}")
        sys.exit(1)

    DEPS_DIR.mkdir(parents=True, exist_ok=True)

    # ── Step 2: clone DIGGER ─────────────────────────────────────────────────
    section("Step 2/5  Cloning DIGGER (web app + preprocessing code)")
    print(f"  Target: {DIGGER_DIR}")
    if not git_clone(DIGGER_URL, DIGGER_DIR, "DIGGER"):
        print(f"\n{RED}DIGGER clone failed.  Check your internet connection and try again.{RESET}")
        sys.exit(1)

    # ── Step 3: download DIGGER data files from Zenodo ───────────────────────
    section("Step 3/5  Downloading DIGGER reference data (PPI/DDI/proteins)")
    if not download_digger_data(DIGGER_DIR):
        print(f"\n{RED}DIGGER data download failed.  See error output above.{RESET}")
        sys.exit(1)

    # ── Step 4: clone and build Jesse ────────────────────────────────────────
    section("Step 4/5  Cloning and building Jesse (graphlet orbit counter)")
    print(f"  Target: {JESSE_DIR}")
    if not git_clone(JESSE_URL, JESSE_DIR, "Jesse"):
        print(f"\n{RED}Jesse clone failed.  Check your internet connection and try again.{RESET}")
        sys.exit(1)

    jesse_jar = build_jesse(JESSE_DIR)
    if jesse_jar is None:
        print(f"\n{RED}Jesse build failed.  See error output above.{RESET}")
        print(f"""
  Common fixes:
    - Make sure Java 11+ is installed (not just Java 8):  java -version
    - If you have Java 25 and the build fails, install Java 17 LTS, set
      JAVA_HOME to the Java 17 folder, then re-run this script.
    - Make sure Maven can reach the internet (it downloads dependencies).
    - Try running manually:
        cd {JESSE_DIR}
        mvn package -DskipTests
      and check the full error message.
""")
        sys.exit(1)

    # ── Step 5: patch config template ────────────────────────────────────────
    section("Step 5/5  Updating config template with resolved paths")
    patch_config(DIGGER_DIR, jesse_jar)

    # ── Done ─────────────────────────────────────────────────────────────────
    print(f"\n{'=' * 50}")
    print(f"{GREEN}{BOLD}Setup complete!{RESET}")
    print(f"""
Next steps:

  1. Copy the config template and edit it:
       cp demo/config_template.yaml my_analysis.yaml

     You only need to fill in TWO things:
       input.sample_dir           — folder with your per-sample TPM files
       validation.subtype_file    — your subtype label file (optional)

     The DIGGER and Jesse paths have already been filled in for you.

  2. Run the pipeline:
       python run_pipeline.py --config my_analysis.yaml

  3. Or run individual stages:
       python run_pipeline.py --config my_analysis.yaml --stages filter network
       python run_pipeline.py --config my_analysis.yaml --stages orbits
       python run_pipeline.py --config my_analysis.yaml --stages cluster validate
""")


if __name__ == "__main__":
    main()
