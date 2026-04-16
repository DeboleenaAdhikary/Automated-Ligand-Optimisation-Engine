# ALOE — Automated Ligand Optimisation Engine

---

## Overview

**ALOE** is a browser-based, zero-installation web application that automates fragment-based ligand design for enzyme targets. It integrates a biologically curated chemical fragment library (~86,000 entries derived from BRENDA enzyme substrates), automated binding-pocket detection, and molecular docking into a single end-to-end pipeline — all accessible from a standard web browser without any local software installation by the end user.

The project addresses a core bottleneck in computational Fragment-Based Drug Discovery (FBDD): the need for domain expertise to orchestrate fragmentation, 3D coordinate generation, pocket detection, and docking into a coherent workflow. ALOE solves this by wrapping all stages in a REST API backend and a reactive single-page frontend with real-time progress feedback.

---

## Key Results

| Enzyme | PDB | Baseline Ligand | Baseline (kcal/mol) | Best ALOE Ligand (kcal/mol) | Improvement |
|---|---|---|---|---|---|
| Acetolactate Synthase | 1BFK | Methyl mandelate | −4.128 | −7.628 | +3.500 (+84.7%) |
| *Thermomyces lanuginosus* Lipase | 1TCA | Ethyl butyrate | −3.672 | −8.841 | +5.169 (+140.8%) |
| Alpha-Glucosidase | 1BAG | DNJ | −3.436 | −10.312 | +6.876 (+200.1%) |

---

## Architecture

ALOE uses a three-tier design:

```
┌─────────────────────────────────────────┐
│   Presentation Layer                    │
│   index.html  (HTML / CSS / JavaScript) │
│   3Dmol.js v2.0.4  ·  SSE progress     │
└────────────────┬────────────────────────┘
                 │ REST API (HTTP/JSON)
┌────────────────▼────────────────────────┐
│   Application Layer                     │
│   app.py  (Flask · Python 3.10)         │
│   RDKit · OpenBabel · BioPython         │
│   AutoDock Vina · P2Rank               │
└────────────────┬────────────────────────┘
                 │ File I/O
┌────────────────▼────────────────────────┐
│   Data Layer                            │
│   Sorted-CSV fragment library           │
│   UUID-based per-session temp dirs      │
└─────────────────────────────────────────┘
```

The frontend communicates with the backend exclusively via predefined HTTP endpoints, making the two layers independently deployable (e.g., frontend on GitHub Pages / Netlify; backend on AWS EC2, Google Cloud Run, or a university HPC node).

---

## Features

- **Interactive 3D Molecular Viewer** — 3Dmol.js surface/cartoon enzyme rendering with ball-and-stick ligand display; click-to-select substitution points directly on the 3D structure.
- **BRICS Fragmentation** — RDKit `BRICSDecompose` generates scaffold + attachment points from any uploaded ligand SMOL/SDF.
- **Fragment Library (86,000 entries)** — Built from >1.36 lakh BRENDA substrate molecules using three complementary fragmentation algorithms: BRICS, RECAP, and rotatable single-bond cleavage. Stored as a sorted CSV enabling binary-search retrieval by molecular weight range.
- **8-Stage Automated Docking Pipeline**:
  1. Scaffold SMILES generation & attachment-point assignment
  2. Apo-enzyme preparation (chain A selection via BioPython)
  3. Binding-pocket detection via **P2Rank** (fpocket fallback / blind docking fallback)
  4. Fragment substitution across user-defined MW windows
  5. Multi-fallback 3D coordinate generation (ETKDGv3 → UFF → OpenBabel)
  6. PDBQT conversion via OpenBabel
  7. AutoDock Vina docking (parallelised, threaded)
  8. Result ranking, ZIP packaging, and download
- **Real-Time Progress** — Server-Sent Events (SSE) stream pipeline stage, percentage, and log messages live to the browser.
- **User Authentication & Job History** — Session-based login, run history, and legacy-import support.
- **7-Step Interactive Tutorial Overlay** — Embedded SVG diagrams walk new users through the full workflow without leaving the page.
- **Cross-Platform Tool Discovery** — Auto-locates Vina and P2Rank executables across conda envs, Homebrew, system paths, and `VINA_PATH` / user-set overrides.

---

## Fragment Library Construction

| Stage | Detail |
|---|---|
| Source | BRENDA Enzyme Database (.mol files, all 6 EC classes) |
| Raw substrates | ~140,000 |
| After sanitisation | ~136,000 |
| Fragmentation methods | BRICS · RECAP · Single-bond cleavage |
| Post-processing | Canonical SMILES · deduplication · valence/ring/radical filter · Lipinski-RO3 filter |
| Final library size | **~86,000 unique fragments** |
| Storage format | Sorted CSV (binary-search by Avg_MW) |

---

## Tech Stack

| Component | Technology |
|---|---|
| Backend framework | Flask (Python 3.10) |
| Cheminformatics | RDKit |
| Format conversion / 3D gen | OpenBabel 3.1.0 |
| Structural biology | BioPython |
| Molecular docking | AutoDock Vina ≥ 1.2.5 |
| Pocket detection | P2Rank |
| 3D visualisation | 3Dmol.js v2.0.4 |
| Frontend | Vanilla HTML5 / CSS3 / JavaScript (single-page) |
| Concurrency | Python `ThreadPoolExecutor` |
| Progress streaming | Server-Sent Events (SSE) |

---

## Requirements

### Python packages

```
flask
rdkit
openbabel-wheel
biopython
pillow
werkzeug
```

Install with:

```bash
pip install flask rdkit openbabel-wheel biopython pillow --break-system-packages
```

### External binaries (must be installed separately)

| Tool | Purpose | Notes |
|---|---|---|
| **AutoDock Vina** ≥ 1.2.5 | Molecular docking | Set `VINA_PATH` env var if not on `PATH` |
| **P2Rank** | Binding-pocket detection | Unzip anywhere; ALOE auto-discovers `p2rank_*` directories |

### Python version

Python **3.10** or later is recommended.

### Java installation

Please install Java since it is required to run P2Rank.

---

## Installation & Usage

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/aloe.git
cd aloe
```

### 2. Install Python dependencies

```bash
pip install flask rdkit openbabel-wheel biopython pillow --break-system-packages
```

### 3. Place the fragment library

Ensure `fragment_library.csv` (the sorted-CSV file) is located in the same directory as `app.py`. The backend resolves paths relative to the script location automatically.

### 4. Configure external tools (optional)

```bash
# If Vina is not on your PATH:
export VINA_PATH=/path/to/vina

# P2Rank: just unzip the distribution anywhere under ~/ or next to app.py
# ALOE will auto-scan common locations (~/Downloads/p2rank_*, etc.)
```

### 5. Run the server

```bash
python app.py
```

ALOE automatically finds a free port starting from 5001 (port 5000 is skipped to avoid macOS AirPlay conflicts) and prints the URL:

```
============================================================
Enzyme Drug Discovery Pipeline - Server Starting
============================================================

Server URL: http://127.0.0.1:5001
```

### 6. Open in your browser

Navigate to the printed URL. The single-page app is served directly by Flask.

---

## API Endpoints

| Method | Endpoint | Description |
|---|---|---|
| `GET` | `/` | Serve the frontend SPA |
| `POST` | `/generate-scaffold` | Parse uploaded PDB/SDF, run BRICS, return scaffold SMILES |
| `POST` | `/generate-job-smiles` | Build substitution-job SMILES with `[*]` attachment points |
| `POST` | `/start-pipeline` | Launch the 8-stage docking pipeline; returns `job_id` |
| `GET` | `/pipeline-status` | Poll job status (stage, percent, done flag) |
| `GET` | `/pipeline-log` | Cursor-based incremental log fetch for live console output |
| `GET` | `/download` | Download result ZIP for a completed job |
| `GET/POST/DELETE` | `/auth/history` | Job history CRUD |
| `POST` | `/auth/register` · `/auth/login` · `/auth/logout` | User auth |

All endpoints return JSON. Binary downloads use `send_file` with `as_attachment=True`.

---

## Validation Enzymes

- **1BFK** — Acetolactate synthase (ALS); herbicide target; baseline ligand: methyl mandelate
- **1TCA** — *Thermomyces lanuginosus* Lipase (TLL); industrial ester synthesis; baseline: ethyl butyrate
- **1BAG** — Alpha-glucosidase (GH31 family); anti-diabetic target; baseline: deoxynojirimycin (DNJ)

---

## Repository Structure

```
aloe/
├── app.py                         # Flask backend — full pipeline logic
├── index.html                     # Single-page frontend
├── fragment_library.csv           # Sorted-CSV fragment library (~86k entries)
├── 22BBT0006_Project2Final_Draft.pdf  # Full project report
└── README.md
```

---

## Abbreviations

| Abbreviation | Meaning |
|---|---|
| ALOE | Automated Ligand Optimisation Engine |
| BRICS | Breaking of Retrosynthetically Interesting Chemical Substructures |
| BRENDA | Braunschweig Enzyme Database |
| FBDD | Fragment-Based Drug Discovery |
| FBLD | Fragment-Based Ligand Design |
| HTS | High-Throughput Screening |
| P2Rank | Protein Binding Site Prediction using Random Forests |
| PDBQT | PDB with partial charges and torsion tree (AutoDock format) |
| RECAP | Retrosynthetic Combinatorial Analysis Procedure |
| SBDD | Structure-Based Drug Design |
| SSE | Server-Sent Events |
| TLL | *Thermomyces lanuginosus* Lipase |

---

## Author

**Deboleena Adhikary**
B.Tech. Biotechnology  

---

## License

This project was submitted as an academic capstone. Please contact the author before reuse or redistribution.
