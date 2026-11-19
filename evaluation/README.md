# OpenBench Evaluation Guide

This directory provides the reproducible Docker environment and end-to-end scoring scripts for all six OpenBench tasks:

```text
evaluation/
├── Dockerfile               # OpenStructure 2.12.0 (`ost`) + RDKit + PoseBusters + Gemmi
├── README.md                # Detailed evaluation instructions & protocol definitions
├── eval_polymer.py          # Polymer tasks (PP, AbAg, ProtRNA, ProtDNA, Monomer)
└── eval_protein_ligand.py   # Protein-ligand task (Runs N' Poses protocol)
```

---

## 1. Build Evaluation Container

Run from the `OpenBench/` root directory:

```bash
docker build -t openbench-eval -f evaluation/Dockerfile evaluation/
```

Verify the environment against ground-truth structures (`selftest`):

```bash
docker run --rm -v "$(pwd):/workspace" -w /workspace openbench-eval bash -c "
  python3 evaluation/eval_polymer.py selftest --benchmark-root . &&
  python3 evaluation/eval_protein_ligand.py selftest --benchmark-root .
"
```

---

## 2. Prediction Input Format

Both `eval_polymer.py` and `eval_protein_ligand.py` accept either:

1. **`--predictions-dir <dir>`** (automatic discovery for OpenDDE / Protenix output layouts):
   - `<dir>/<case>/<case>/seed_<seed>/predictions/<case>_sample_<0..4>.cif`
   - `<dir>/<case>/<case>/seed_<seed>/predictions/<case>_summary_confidence_sample_<0..4>.json`
   - (Also supports `<dir>/<case>/seed_<seed>/predictions/...` and `<dir>/<task>/<case>/...`.)
2. **`--manifest <predictions.csv>`** (any model output layout):
   - CSV columns: `pdb_id,seed,sample,prediction_path,confidence_path`
   - Each case is evaluated over 5 seeds (`101, 102, 103, 104, 105`) × 5 diffusion samples = **25 predictions**.

---

## 3. Polymer Evaluation (`eval_polymer.py`)

Covers `protein-protein` (310 interfaces / 297 cases), `antibody-antigen` (223 interfaces / 117 cases), `protein-rna` (115 interfaces / 93 cases), `protein-dna` (208 interfaces / 148 cases), and `protein-monomer` (504 chains / 504 cases).

```bash
docker run --rm -v "$(pwd):/workspace" -v "/path/to/preds:/preds:ro" -w /workspace openbench-eval \
  python3 evaluation/eval_polymer.py run \
    --benchmark-root . \
    --predictions-dir /preds \
    --output-dir results/polymer \
    --workers 16
```

### Protocol Details
- **Interfaces (`protein-protein`, `antibody-antigen`, `protein-rna`, `protein-dna`)**:
  - Scored with `ost compare-structures --fault-tolerant --min-pep-length 4 --min-nuc-length 4 --lddt --rigid-scores --tm-score --dockq`.
  - Reference chain pairs `(interface_chain_id_1, interface_chain_id_2)` (`label_asym_id` in `structures/<task>/<case>.cif`) are matched against `dockq_interfaces`. Primary metric: **`DockQ >= 0.23`** (medium `>= 0.49`, high `>= 0.80` also reported).
  - Three reference-mapping edge cases are handled automatically inside `eval_polymer.py`:
    - `9r2z-assembly1`: reference converted to PDB without `SEQRES` (`A1 -> A`, `B1 -> B`) so chains differing only by an N-terminal Met cluster into one chemical group, then mapped back to `A1`/`B1`.
    - `10yz-assembly1`: evaluated with `--min-nuc-length 3` (4 resolved nucleotides).
    - `24br-assembly1`: evaluated with `--chem-map-seqid-thresh 0` (incomplete 5' nucleotide).
- **Monomers (`protein-monomer`)**:
  - Uses `chain_maps/protein-monomer/<case>.json` (`reference_to_input`) and Gemmi to extract the single target chain from both reference and prediction (removing waters, ligands, and hydrogens and reassigning `label_seq_id`), then runs `ost compare-structures` to obtain chain-level **`lDDT`** (primary) and **`TM-score`**.
- **Aggregation**:
  - `rank` (**Top-ranked**): single prediction with highest case-level `ranking_score` among the 25 predictions (ties broken by lower `seed`, `sample`).
  - `oracle` (**Best-of-25**): best score across the 25 predictions per evaluated unit.
  - Missing or unscored units count as `0.0`.
- **Outputs**: `units_rank.csv`, `units_oracle.csv`, `summary.csv`.

---

## 4. Protein–Ligand Evaluation (`eval_protein_ligand.py`, Runs N' Poses Protocol)

Covers `protein-ligand` (**952 systems / 1,129 proper ligands / 514 SuCOS-pocket clusters**), strictly following **Runs N' Poses** (Škrinjar et al., 2025, Methods 5.2–5.4).

```bash
docker run --rm -v "$(pwd):/workspace" -v "/path/to/preds:/preds:ro" -w /workspace openbench-eval \
  python3 evaluation/eval_protein_ligand.py run \
    --benchmark-root . \
    --predictions-dir /preds \
    --output-dir results/protein-ligand \
    --workers 16
```

### Protocol Details
1. **SMILES-Templated Complex Splitting**:
   - Because `inputs/protein-ligand/<system_id>.json` encodes all bioassembly-1 ligands as SMILES (including ions and crystallization additives), predicted mmCIFs do not carry wwPDB CCD residue/atom names.
   - `eval_protein_ligand.py` splits each predicted complex into `receptor.cif` (polymer chains) plus one `<chain_id>.sdf` per ligand chain with bond orders restored from the input SMILES and predicted 3D heavy-atom coordinates attached.
   - Dative bonds (SDF bond type `9`) in metallo-ligand SDFs are normalized to single bonds (`1`) in temporary copies.
2. **OpenStructure & PoseBusters Scoring**:
   - Runs `ost compare-ligand-structures --fault-tolerant --lddt-pli --lddt-pli-add-mdl-contacts --rmsd` against `structures/protein-ligand/<system_id>/receptor.cif` and `ligand_files/*.sdf` using the same unified OpenStructure runtime (`2.12.0`) as the polymer tasks.
   - Evaluates only the **1,129 proper ligands** listed in `targets/protein-ligand.csv` (ions/additives are modeled in the input but excluded from scoring).
   - Runs **PoseBusters** (`PoseBusters(config="dock")`) on each assigned proper ligand SDF conditioned on the predicted receptor (`pb_valid` = all default checks pass).
3. **Metrics, Ranking, & Stratified Reporting**:
   - **Primary success criterion**: **`BiSyRMSD < 2.0 Å` and `LDDT-PLI > 0.8`** (Runs N' Poses Methods 5.3 / Fig. 1A).
   - **Physical-validity subset**: **`BiSyRMSD < 2.0 Å`, `LDDT-PLI > 0.8`, and `PB-valid`** (Runs N' Poses Fig. 1E), alongside standalone `BiSyRMSD < 2.0 Å` and `PB-valid` rates.
   - **Primary Top-1 selection (`top_iptm`)**: per proper ligand, selects the prediction with the highest mean **ligand-to-protein chain-pair iPTM** (`lig_prot_iptm`) across the 25 predictions (Methods 5.2); also reports `top_ranking_score`, `best_of_25`, and `mean_of_25` with 1,000-sample bootstrap 95% confidence intervals. Missing/unscored proper ligands count as failure (`0`).
   - **Stratified reporting**: stratified across the 8 SuCOS-pocket `similarity_bin` buckets (`0-20` through `80-100`, relative to the 2025-06-30 cutoff) over both **all proper ligands (`n = 1,129`)** and **SuCOS-pocket cluster representatives (`is_cluster_representative = true`, `n = 514`)**.
- **Outputs**: `per_prediction_ligand.csv`, `per_ligand_selection.csv`, `summary_overall.csv`, `summary_by_similarity_bin.csv`.
