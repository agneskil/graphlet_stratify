"""
Python wrapper around the Jesse graphlet orbit counter.

Jesse is an external Java tool.  This module does NOT modify Jesse and
only automates calling it across many network files.

Jesse supports graphlet sizes 3–5 (and higher with generated files).

You must have cloned / installed Jesse. See README.md "External dependencies" 
for setup instructions.

"""

import logging
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

_JESSE_PROMPTS = {
    3: "3\nn\ny\nn\ny\n{input}\n{output}\n",
    4: "4\nn\ny\nn\ny\n{input}\n{output}\n",
    5: "5\nn\ny\nn\ny\n{input}\n{output}\n",
}


def run_jesse(
    input_dir: str | Path,
    output_dir: str | Path,
    jesse_jar: str | Path,
    graphlet_size: int = 5,
    file_pattern: str = "*.txt",
) -> None:
    """
    Run Jesse on every network file in 'input_dir'.

    """
    input_dir  = Path(input_dir)
    output_dir = Path(output_dir)
    jesse_jar  = Path(jesse_jar)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not jesse_jar.exists():
        raise FileNotFoundError(
            f"Jesse JAR not found: {jesse_jar}\n"
            "Build Jesse with: cd <jesse_dir> && mvn package\n"
            "Then pass the resulting JAR path as --jesse_jar."
        )

    if graphlet_size not in _JESSE_PROMPTS:
        raise ValueError(f"graphlet_size must be 3, 4, or 5, got {graphlet_size}")

    network_files = sorted(input_dir.glob(file_pattern))
    if not network_files:
        raise FileNotFoundError(
            f"No files matching '{file_pattern}' in {input_dir}"
        )

    log.info(f"Running Jesse on {len(network_files)} network files "
             f"(graphlet_size={graphlet_size}) ...")

    n_ok = 0
    n_fail = 0

    for i, net_file in enumerate(network_files, 1):
        out_file = output_dir / (net_file.stem + "_results.txt")
        stdin_str = _JESSE_PROMPTS[graphlet_size].format(
            input=net_file,
            output=out_file,
        )

        log.info(f"  [{i}/{len(network_files)}] {net_file.name} ...")

        try:
            result = subprocess.run(
                ["java", "-jar", str(jesse_jar)],
                input=stdin_str,
                capture_output=True,
                text=True,
                timeout=600,  # 10-minute per-file timeout
            )
            if result.returncode != 0:
                log.warning(
                    f"  Jesse returned non-zero exit code for {net_file.name}:\n"
                    f"  {result.stderr.strip()[:300]}"
                )
                n_fail += 1
            else:
                log.info(f"  → {out_file.name}")
                n_ok += 1
        except subprocess.TimeoutExpired:
            log.warning(f"  Timeout on {net_file.name} — skipping.")
            n_fail += 1
        except FileNotFoundError:
            raise RuntimeError(
                "Java not found on PATH.  Install Java (>= 8) and make sure "
                "'java' is accessible from your terminal."
            )

    log.info(
        f"Jesse done: {n_ok} succeeded, {n_fail} failed.  "
        f"Results in {output_dir}."
    )
