# graphlet_stratify

**Patient stratification from isoform-aware graphlet orbit signatures.**

A generalised pipeline for unsupervised patient stratification using
isoform-aware protein–protein interaction (PPI) networks and graphlet
orbit frequency vectors.

For each sample, only the dominant isoform (highest-TPM transcript) per
gene is used, and a PPI edge is retained only if the two dominant isoforms
have domain–domain interaction (DDI) support. Graphlet orbit features of
the resulting networks are clustered to stratify patients, with optional
validation against known subtype labels.

---

## Installation and setup

### Step 1 — Prerequisites (install once)

Install these four tools first (skip any you already have):

- **Python ≥ 3.10** — [python.org](https://www.python.org/downloads/)
- **Git** — [git-scm.com](https://git-scm.com/downloads)
- **Java ≥ 8** — [adoptium.net](https://adoptium.net/)
- **Maven** — [maven.apache.org](https://maven.apache.org/download.cgi)

Verify each is installed:

​```bash
python --version && git --version && java -version && mvn --version
​```

### Step 2 — Clone this repository and install Python dependencies

```bash
git clone https://github.com/YOUR_USERNAME/graphlet_stratify
cd graphlet_stratify
pip install -r requirements.txt
```

**Important:** always run all scripts from the `graphlet_stratify/` project root directory. The pipeline imports its own modules using relative paths, so running from elsewhere will fail.

### Step 3 — Run the setup script

```bash
python setup_dependencies.py
```

The script will:
- Clone DIGGER into `dependencies/DIGGER/` and download data files from Zenodo
- Clone Jesse into `dependencies/jesse/` and build the JAR with Maven
- Verify all expected data files are present
- Write the correct paths into `demo/config_template.yaml`

---

## Running the analysis

You need two input files: a transcript expression matrix (rows = transcripts, columns = samples, values = log2(TPM+1)) and a subtype label file (sample ID, subtype).

**Option 1 — straight from the command line:**

​```bash
python run_pipeline.py \
    --matrix   path/to/tpm_matrix.tsv \
    --subtypes path/to/subtypes.tsv \
    --output   results/
​```

This runs all five stages with defaults.

**Option 2 — config file (for tuning parameters):**

See `demo/config_template.yaml` for the full list of options with explanation comments.

​```bash
cp demo/config_template.yaml my_analysis.yaml
python run_pipeline.py --config my_analysis.yaml
​```

### Running analysis stages seperately 

Pass `--stages` to skip work you don't need to repeat (useful when iterating on clustering after the orbits stage has already finished):

​```bash
python run_pipeline.py --matrix tpm_matrix.tsv --stages filter network
python run_pipeline.py --matrix tpm_matrix.tsv --stages orbits
python run_pipeline.py --matrix tpm_matrix.tsv --stages cluster validate
​```

## Input format

**Expression matrix** — a single transcript × sample file:

```
transcript_id     TCGA-A2-A04P  TCGA-A8-A06R  TCGA-A8-A08H
ENST00000000233   12.5          8.3            0.0
ENST00000000412   0.0           134.2          45.9
```

- Rows = transcripts (Ensembl IDs, e.g. `ENST00000366570` or `ENST00000366570.3`)
- Columns = patient samples
- Values = log2(TPM+0.001)
- TSV or CSV (auto-detected)

**Subtype label file** (for validation and PCA colouring, optional):

```
sample_id          PAM50
TCGA-A1-A0SB-01A   LumA
TCGA-A1-A0SD-01A   LumB
```

- Two columns, any header names
- TSV or CSV

---

## Output structure

After a full run, results are organised as:

```
results/
  switching_report.tsv          gene-level isoform switching statistics
  whitelist.txt                 retained transcript IDs
  expression_filtered.tsv       expression matrix, whitelisted transcripts only
  networks/                     one .txt network per sample (Jesse-compatible)
  orbits/                       one orbit result file per sample (Jesse output)
  clusters/
    per_gene/
      gene<ID>_clusters.csv     cluster assignment per sample, per gene
      gene<ID>_pca.png          PCA plot (if clustering.save_pca: true)
    per_orbit/
      orbit<N>_clusters.csv     cluster assignment per sample, per orbit type
      orbit<N>_pca.png          PCA plot (if clustering.save_pca: true)
  validation/
    per_gene/
      ranking_summary.csv       genes ranked by ARI
    per_orbit/
      ranking_summary.csv       orbit types ranked by ARI
```

---
## Reproducing the BRCA analysis

The original analysis used:
- TCGA-BRCA RNA-seq data (1,212 samples, log2(TPM+0.001) from UCSC Xena)
- PAM50 subtype labels from `TCGA-BRCA_Xena_Subtypes.tsv`
- `tpm_threshold: 0`, `switch_rate: 0.20`, `n_clusters: 5`, `graphlet_size: 5`
