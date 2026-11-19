# OpenBench

All-atom biomolecular structure-prediction benchmark built from PDB entries released **2025-07-01 to 2026-10-07**.

## Benchmark Summary

| Task | Cases | Evaluated Units | Metric |
|---|---:|---:|---|
| `protein-protein` | 297 | 310 interfaces | DockQ ≥ 0.23 |
| `antibody-antigen` | 117 | 223 interfaces | DockQ ≥ 0.23 |
| `protein-rna` | 93 | 115 interfaces | DockQ ≥ 0.23 |
| `protein-dna` | 148 | 208 interfaces | DockQ ≥ 0.23 |
| `protein-monomer` | 504 | 504 chains | Chain lDDT |
| `protein-ligand` | 952 | 1,129 ligands | BiSyRMSD < 2 Å & LDDT-PLI > 0.8 (+ PB-valid) |

## Directory Layout

All data directories share the same six task names (`<task>`):

```text
OpenBench/
├── targets/<task>.csv                     # One row per evaluated unit (2,489 total units)
├── metadata/
│   ├── cases.csv                          # Global index of all 2,085 distinct cases
│   └── <task>.csv                         # Per-task metadata (protein-ligand.csv includes PLINDER/similarity/cluster fields)
├── inputs/<task>/<case>.json              # OpenDDE evaluation JSON inputs (modelSeeds: [101..105], relative MSA paths)
├── structures/
│   ├── <polymer-task>/<case>.cif          # Biological assembly 1 reference mmCIF
│   └── protein-ligand/<system_id>/        # PLINDER ground truth: receptor.cif, system.cif, sequences.fasta, ligand_files/*.sdf
├── chain_maps/<polymer-task>/<case>.json  # Reference-to-input chain ID mapping (5 polymer tasks)
├── msa/<task>/<md5>/                      # LocalColabFold protein MSAs (no templates, no RNA MSA)
│   ├── pairing.a3m                        # TaxID-annotated UniRef100 hits (>UniRef100_<acc>_<taxid>/)
│   └── non_pairing.a3m                    # Environmental (colabfold_envdb) and unmapped hits
└── evaluation/                            # Docker environment & scoring scripts (see evaluation/README.md)
    ├── Dockerfile
    ├── README.md
    ├── eval_polymer.py
    └── eval_protein_ligand.py
```

## Data Download

This repository contains `targets/`, `metadata/`, `inputs/`, `chain_maps/`, and `evaluation/`. Reference structures (`structures/`) and precomputed protein MSAs (`msa/`) are distributed via [GitHub Releases (`v1.0.0`)](https://github.com/aurekaresearch/OpenBench/releases/tag/v1.0.0).

Download, verify, and extract directly in the repository root:

```bash
wget https://github.com/aurekaresearch/OpenBench/releases/download/v1.0.0/structures.tar.gz
wget https://github.com/aurekaresearch/OpenBench/releases/download/v1.0.0/msa.tar.gz
wget https://github.com/aurekaresearch/OpenBench/releases/download/v1.0.0/SHA256SUMS
sha256sum -c SHA256SUMS
tar -xzf structures.tar.gz
tar -xzf msa.tar.gz
```

## Inference

Run from the `OpenBench/` root so relative `msa/<task>/<md5>/{pairing,non_pairing}.a3m` paths resolve directly:

- **Inputs**: `inputs/<task>/<case>.json` are formatted for OpenDDE evaluation with 5 default seeds (`"modelSeeds": [101, 102, 103, 104, 105]`) × 5 diffusion samples = **25 predictions** per case.
- **Composition**: Polymer tasks include all bioassembly-1 polymer chains (`proteinChain`, `rnaSequence`, `dnaSequence`); `protein-ligand` follows Runs N' Poses and includes all bioassembly-1 protein chains and ligands encoded as SMILES (including ions and additives; selenocysteine `U` encoded as `C` + `CCD_SEC`).

## Evaluation

Build the evaluation container via Docker and run the automated scoring scripts in [`evaluation/`](evaluation/README.md):

```bash
# 1. Build evaluation image
docker build -t openbench-eval -f evaluation/Dockerfile evaluation/

# 2. Score Polymer tasks (PP, AbAg, ProtRNA, ProtDNA, Monomer)
docker run --rm -v "$(pwd):/workspace" -v "/path/to/preds:/preds:ro" -w /workspace openbench-eval \
  python3 evaluation/eval_polymer.py run --benchmark-root . --predictions-dir /preds --output-dir results/polymer

# 3. Score Protein-Ligand task (Runs N' Poses protocol: OpenStructure + PoseBusters `dock`)
docker run --rm -v "$(pwd):/workspace" -v "/path/to/preds:/preds:ro" -w /workspace openbench-eval \
  python3 evaluation/eval_protein_ligand.py run --benchmark-root . --predictions-dir /preds --output-dir results/protein-ligand
```

See [`evaluation/README.md`](evaluation/README.md) for full protocol specifications, ranking rules (`ranking_score` for polymer tasks; per-ligand `lig_prot_iptm` for `protein-ligand`), special-case handling, and similarity-stratified reporting.

## References

- FoldBench: [github.com/BEAM-Labs/FoldBench](https://github.com/BEAM-Labs/FoldBench)
- Runs N' Poses: Škrinjar et al. (2025), [github.com/plinder-org/runs-n-poses](https://github.com/plinder-org/runs-n-poses)
- PLINDER: [github.com/plinder-org/plinder](https://github.com/plinder-org/plinder) · OpenStructure: [openstructure.org](https://openstructure.org/) · PoseBusters: [github.com/maabuu/posebusters](https://github.com/maabuu/posebusters)

## License

- **Code**: Evaluation scripts and container definitions (`evaluation/`) are released under the [MIT License](LICENSE).
- **Dataset & Upstream Assets**: Benchmark annotations, inputs, and PLINDER / Runs N' Poses-derived `protein-ligand` ground-truth assets are distributed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/); experimental macromolecular structures originate from the [wwPDB](https://www.wwpdb.org/) ([CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)).
