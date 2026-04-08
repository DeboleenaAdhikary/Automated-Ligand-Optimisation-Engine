# ALOE

ALOE is a Flask-based web application for enzyme-ligand design, scaffold-guided analog generation, and docking-driven hit ranking. Starting from an enzyme `.pdb` file and a ligand `.sdf` file, the app lets you visualize a scaffold, mark substitution sites interactively, generate fragment-substituted candidates, run docking with AutoDock Vina, and review ranked outputs directly in the browser.

## Key capabilities

- Upload an enzyme structure (`.pdb`) and ligand structure (`.sdf`)
- Choose between using the original ligand as the scaffold or extracting a BRICS-derived scaffold
- Interactively mark up to 3 carbon substitution positions in a 3D viewer
- Generate substitution jobs from the bundled fragment library with molecular-weight filtering
- Prepare ligands and receptors automatically for docking
- Detect pockets with P2Rank and dock candidates with AutoDock Vina
- Download ranked docking outputs as CSV and ZIP bundles
- Persist user login and run history in the web interface

## Golden validation benchmark

A reference validation case has been established using **1TCA + ethyl butyrate**.

### Validation summary

| Case | System | Best predicted affinity (kcal/mol) | Evidence |
| --- | --- | ---: | --- |
| Baseline ligand | 1TCA + ethyl butyrate | `-3.672` | `validation/Validation_exp.png` |
| Best ALOE-ranked redesigned hit | 1TCA + `ligand_402` | `-8.841` | `validation/docking_results_ranked-10.csv` |
| Improvement over baseline | `ligand_402` vs ethyl butyrate | `5.169` kcal/mol more favorable | Calculated from the two values above |

### Top redesigned hits from the ranked docking CSV

| Rank | Candidate | Conformation 1 affinity (kcal/mol) |
| --- | --- | ---: |
| 1 | `ligand_402` | `-8.841` |
| 2 | `ligand_275` | `-8.571` |
| 3 | `ligand_274` | `-8.555` |

> Note: these values are AutoDock Vina docking scores (predicted affinities), not experimental binding constants.

## Repository layout

```text
.
├── app.py
├── index.html
├── 350_frag_sorted.csv
├── 350_frag_cleaned.txt
├── static/
│   └── ALOE_logo.png
├── validation/
│   ├── Validation_exp.png
│   └── docking_results_ranked-10.csv
└── enzyme_processor_results/
```

## Prerequisites

Make sure the following are installed before running the app:

- Python 3.9 or newer
- Flask
- RDKit
- Biopython
- Pillow
- Open Babel
  - Python bindings for `openbabel`
  - Command-line executable `obabel`
- AutoDock Vina
  - executable `vina` available on `PATH`, or exposed via `VINA_PATH`
- P2Rank
  - executable `prank` available for automatic pocket prediction

The repository also expects the bundled fragment library files to be present:

- `350_frag_sorted.csv`
- `350_frag_cleaned.txt`

## Installation

### Recommended: Conda environment

```bash
conda create -n aloe python=3.10 -y
conda activate aloe
conda install -c conda-forge flask rdkit biopython pillow openbabel vina
```

### Pip-based alternative

If you are not using Conda, the Python packages can also be installed with pip:

```bash
pip install flask biopython pillow rdkit openbabel-wheel
```

### External tools

#### AutoDock Vina

Make sure `vina` is on your shell `PATH`. If it is installed in a non-standard location, set:

```bash
export VINA_PATH=/full/path/to/vina
```

#### P2Rank

Install P2Rank so that the `prank` executable is discoverable. The app can auto-locate common installs such as:

- `~/Downloads/p2rank_*`
- a nearby extracted `p2rank_*` directory

## Running the application

From the project directory:

```bash
python app.py
```

The server prints the local URL in the terminal, for example:

```text
http://127.0.0.1:5001
```

The app automatically selects a free port starting at `5001`, so use the exact URL shown in the terminal.

## Web workflow

1. Open the local server URL in your browser.
2. Create an account or sign in if you want persistent run history.
3. Upload:
   - an enzyme file in `.pdb` format
   - a ligand file in `.sdf` format
4. Choose scaffold mode:
   - use the original ligand as-is, or
   - generate a BRICS-derived scaffold
5. Mark carbon atoms in the 3D viewer to define substitution sites.
6. Review the generated substitution jobs.
7. Click **Process** to run:
   - scaffold preparation
   - apo enzyme preparation
   - pocket prediction
   - fragment substitution
   - PDBQT conversion
   - AutoDock Vina docking
   - result ranking and packaging
8. Download the ranked CSV and ZIP outputs from the results panel.

## Outputs

The pipeline writes results under:

```text
enzyme_processor_results/jobs/<job_id>/
```

Common output artifacts include:

- `docking_results_ranked.csv`
- `top_ligands.zip`
- intermediate uploaded files and docking outputs per job

The web app also creates small persistence files for authentication and run history:

- `aloe_users.json`
- `aloe_history.json`
- `aloe_secret.key`

## Validation artifacts

- Screenshot of baseline validation run: [validation/Validation_exp.png](validation/Validation_exp.png)
- Ranked docking CSV used for the redesigned-hit comparison: [validation/docking_results_ranked-10.csv](validation/docking_results_ranked-10.csv)

## Notes

- ALOE is intended for computational screening and prioritization.
- Docking scores should be interpreted as relative ranking signals, not direct experimental binding measurements.
- If `vina`, `obabel`, or `prank` are missing, the server logs will indicate which executable could not be found.
