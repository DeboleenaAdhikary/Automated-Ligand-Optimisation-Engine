# ALOE — Automated Ligand Optimisation Engine

> **Fragment-based ligand design webtool for enzyme-focused drug discovery**  
> University of Hyderabad · School of Life Sciences · Institution of Eminence

---

## Overview

ALOE is a browser-based, full-stack cheminformatics pipeline that automates **fragment-based ligand design** for any enzyme target. Given an enzyme structure (`.pdb`) and a known ligand (`.sdf`/`.mol`/`.pdb`), ALOE systematically substitutes drug-relevant chemical fragments at user-selected positions on the ligand scaffold, docks every generated analogue into the enzyme's binding pocket, and returns a ranked list of candidates with their predicted binding affinities — all without requiring any specialist software on the user's machine.

---

## Features

- **Interactive 3D ligand viewer** — 3Dmol.js renders the scaffold directly in the browser; users click carbon atoms to select substitution sites (up to 3)
- **Automatic SMILES generation** — clicking an atom instantly generates a valence-safe SMILES with a `[*]` attachment point via RDKit
- **BRICS scaffold fragmentation** — optionally decomposes the ligand at retrosynthetically meaningful bonds and uses the largest fragment as the scaffold
- **86,000-fragment library** — sourced from 1.4 lakh+ BRENDA enzyme substrates; filtered to ≤ 350 Da following Lipinski-inspired rules
- **Binary-search fragment screening** — a pre-sorted CSV indexed by average MW allows O(log n + k) lookup instead of an O(n) full scan, dramatically reducing screening time
- **SMILES deduplication** — canonicalised RDKit SMILES prevent chemically identical products from being docked twice
- **Automated pocket detection** — P2Rank ML model locates binding pockets from the apo enzyme surface
- **AutoDock Vina docking** — every generated ligand PDBQT is docked into the detected pocket; blind docking fallback when pocket data is insufficient
- **Ranked results CSV** — top 3 ligands with binding affinities per conformation (number of columns matches substitution sites chosen)
- **Downloadable ZIP** — all docked ligand PDBQTs + the prepared `apo_enzyme.pdbqt` ready for PyMOL or further studies
- **Run history** — previous runs stored per user via localStorage
- **Interactive tutorial** — step-by-step guide with annotated SVG diagrams built into the interface

---

## Pipeline Architecture

```
Upload enzyme.pdb + ligand.sdf
        │
        ▼
Scaffold Preparation
  RDKit parses ligand → optional BRICS fragmentation
  → canonical SMILES round-trip (strips artefacts)
  → ETKDGv3(seed=42) + UFF 3D optimisation → SDF for viewer
        │
        ▼
Apo Enzyme Preparation
  BioPython strips heteroatoms & non-chain-A residues
  → clean apo PDB → OpenBabel → PDBQT
        │
        ▼
Pocket Detection
  P2Rank ML model scores protein surface
  → ranked pockets with 3D centre + box dimensions
        │
        ▼
Fragment Substitution
  Binary search on sorted fragment CSV (Avg_MW column)
  → only MW-window rows loaded
  → RDKit adds [*] branch at clicked carbon (GetTotalNumHs check)
  → _combine_scaffold_fragment() grafts each fragment
  → canonical SMILES deduplication
        │
        ▼
PDBQT Conversion
  RDKit ETKDGv3 → MMFF/UFF geometry optimisation
  → OpenBabel SDF → PDBQT (Gasteiger charges)
        │
        ▼
AutoDock Vina Docking
  Per-ligand config.txt (centre, box, num_modes, energy_range)
  → REMARK VINA RESULT parsed from docked PDBQT
  → padded conformations marked with * for transparency
        │
        ▼
Ranked Results CSV + ZIP
```

---

## Fragment Library

| Stage | Count |
|---|---|
| Source substrate `.mol` files (BRENDA) | 1,40,000+ |
| Fragmentation methods | Simple bond cuts · BRICS · RECAP |
| Raw fragments generated | 4,63,000+ |
| MW range before filtering | ~1 – 20,000 Da |
| Cutoff applied | ≤ 350 Da (Lipinski-inspired) |
| **Final library size** | **~86,000 fragments** |

The library is distributed as a pre-sorted CSV (`350_frag_sorted.csv`) with columns `Fragment_SMILES`, `Avg_MW`, `Exact_MW`, sorted ascending by `Avg_MW` to enable binary-search screening.

---

## Tech Stack

| Layer | Technology |
|---|---|
| Backend | Python 3.10 · Flask |
| Cheminformatics | RDKit · OpenBabel · BioPython |
| 3D viewer | 3Dmol.js |
| Pocket detection | P2Rank |
| Docking | AutoDock Vina 1.2.5 |
| Frontend | Vanilla HTML/CSS/JS |
| Concurrency | Python `threading` — non-blocking pipeline with real-time progress polling |

---

## Installation

### Prerequisites

```bash
conda create -n rdkit_env python=3.10
conda activate rdkit_env
pip install flask rdkit openbabel-wheel biopython pillow werkzeug
conda install -c conda-forge vina
```

P2Rank must be downloaded separately:  
[https://github.com/rdk/p2rank/releases](https://github.com/rdk/p2rank/releases)

### Fragment Library

Place the sorted fragment CSV in your chosen directory and update the default path in `Fullcode_updated_fixed.py`:

```python
_base = "/path/to/Fragment_lib/350_frag_sorted.csv"
```

### Run

```bash
conda activate rdkit_env
cd /path/to/Run/
python Fullcode_updated_fixed.py
# Open http://127.0.0.1:5001 in your browser
```

---

## Usage

1. **Upload** enzyme `.pdb` and ligand `.sdf` / `.mol` / `.pdb`
2. **Choose scaffold mode** — use molecule as-is, or apply BRICS fragmentation
3. **Click carbon atoms** in the 3D viewer to select substitution sites (up to 3)
4. **Set molecular weight range** per job (default 0–350 Da)
5. **Click Process** — the pipeline runs automatically
6. **Download** the ranked CSV and PDBQT ZIP when complete

---

## Example Results

Using **1BFK** (mandelate racemase, *Pseudomonas putida*) with **methyl mandelate** as the baseline ligand:

| | Baseline | Best Fragment Variant |
|---|---|---|
| Ligand | Methyl mandelate | ligand_65 |
| Best Vina Affinity | −4.117 kcal/mol | −4.657 kcal/mol |
| Improvement | — | **+0.54 kcal/mol** |

All fragment-substituted variants outperformed the baseline, confirming that the binding pocket contained unfilled space exploitable by fragment growing.

---

## Tested Enzyme Systems

| Enzyme | Organism | PDB | Application |
|---|---|---|---|
| Mandelate racemase | *Pseudomonas putida* | 1BFK | Baseline validation |
| α-Amylase | *Bacillus subtilis* | 1BAG | Starch processing / food industry |
| TEM-1 β-Lactamase | *E. coli* | 1BTL | Antibiotic resistance |
| Lipase B (CALB) | *Candida antarctica* | 1TCA | Pharmaceutical / biodiesel synthesis |
| AHAS | *Candida albicans* | 6DEQ | Antifungal drug target |

---

## Project Structure

```
Run/
├── Fullcode_updated_fixed.py     # Flask backend — full pipeline
├── enzyme_processor.html         # Frontend — single-page app
├── 350_frag_sorted.csv           # Pre-sorted fragment library
├── enzyme_processor_results/     # Output directory (auto-created)
│   ├── uploads/                  # Uploaded enzyme & ligand files
│   ├── job{N}_ligands/           # Generated ligand SMI + SDF + PDBQT
│   ├── docking_pocket_1/         # Vina docking output per pocket
│   ├── docking_results_ranked.csv
│   └── top_ligands.zip
└── README.md
```

---

## Limitations & Future Work

- Currently processes **chain A** of the enzyme only; multi-chain support planned
- Fragment library covers ≤ 350 Da; a higher-MW tier is under consideration
- Docking scoring uses Vina's empirical force field — MD-based rescoring is a planned extension
- Authentication system uses browser localStorage; a proper database backend is planned for multi-user deployment

---

## Citation & Acknowledgements

- AutoDock Vina: Eberhardt et al., *J. Chem. Inf. Model.* 2021; Trott & Olson, *J. Comput. Chem.* 2010
- P2Rank: Krivák & Hoksza, *J. Cheminformatics* 2018
- Fragment library sourced from the **BRENDA Enzyme Database** (brenda-enzymes.org)
- RDKit: [rdkit.org](https://www.rdkit.org)

---

## Author

**Deboleena Adhikary**  
School of Life Sciences, University of Hyderabad  
Institution of Eminence — *National Needs, Global Standards*

---

*ALOE is developed as a capstone project for academic and research purposes.*
