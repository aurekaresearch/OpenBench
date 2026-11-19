#!/usr/bin/env python3
"""OpenBench Polymer Evaluation Suite.

Scores and aggregates predictions across the five polymer tasks:
  - protein-protein   (310 interfaces across 297 cases; primary metric: DockQ >= 0.23)
  - antibody-antigen  (223 interfaces across 117 cases; primary metric: DockQ >= 0.23)
  - protein-rna       (115 interfaces across 93 cases;  primary metric: DockQ >= 0.23)
  - protein-dna       (208 interfaces across 148 cases; primary metric: DockQ >= 0.23)
  - protein-monomer   (504 chains across 504 cases;     primary metric: Chain lDDT)

Protocol details:
  - Sampling: 5 seeds (101..105) x 5 diffusion samples = 25 predictions per case.
  - Interface tasks: OpenStructure `compare-structures --fault-tolerant --min-pep-length 4
    --min-nuc-length 4 --lddt --rigid-scores --tm-score --dockq`, matching reference chain
    pairs `(interface_chain_id_1, interface_chain_id_2)` against `dockq_interfaces`.
    Includes automatic handling of the 3 special reference-mapping cases:
      * `9r2z-assembly1`: reference converted to PDB without SEQRES (`A1->A`, `B1->B`) so
        chains differing only by an N-terminal Met form one chemical group, then mapped back.
      * `10yz-assembly1`: `--min-nuc-length 3` (4 resolved nucleotides).
      * `24br-assembly1`: `--chem-map-seqid-thresh 0` (incomplete 5' nucleotide).
  - Monomer task: single target chain extracted from reference (`chain_id`) and prediction
    (`chain_maps/protein-monomer/<case>.json`) via Gemmi (waters, ligands, and hydrogens
    removed; `label_seq_id` reassigned) before running `ost compare-structures`.
  - Aggregation:
      * `rank` (Top-ranked): single prediction with highest case-level `ranking_score`
        among the 25 predictions (ties broken by lower `seed`, `sample`).
      * `oracle` (Best-of-25): best score across the 25 predictions for each evaluated unit.
      * Missing or unscored units count as failure (`0.0`).
"""
import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

POLYMER_TASKS = [
    "protein-protein",
    "antibody-antigen",
    "protein-rna",
    "protein-dna",
    "protein-monomer",
]
INTERFACE_TASKS = {
    "protein-protein",
    "antibody-antigen",
    "protein-rna",
    "protein-dna",
}
SPECIAL_INTERFACE_FLAGS = {
    "9r2z-assembly1": [],
    "10yz-assembly1": ["--min-nuc-length", "3"],
    "24br-assembly1": ["--chem-map-seqid-thresh", "0"],
}
DEFAULT_SEEDS = [101, 102, 103, 104, 105]


def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not fieldnames:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def resolve_ost_bin(ost_bin: str | None = None) -> str:
    if ost_bin:
        return ost_bin
    found = shutil.which("ost")
    if found:
        return found
    for candidate in ("/opt/conda/envs/ost/bin/ost", "/usr/local/bin/ost"):
        if Path(candidate).is_file():
            return candidate
    raise RuntimeError("Could not find `ost` executable on PATH.")


def extract_single_chain(src_cif: Path, chain_id: str, dst_cif: Path) -> None:
    """Extract a single chain (matched by subchain/label_asym_id, fallback to chain name) into chain A."""
    import gemmi

    st = gemmi.read_structure(str(src_cif))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.remove_hydrogens()

    out = gemmi.Structure()
    out.spacegroup_hm = st.spacegroup_hm
    out.cell = st.cell
    ch = gemmi.Chain("A")
    for c in st[0]:
        for r in c:
            if r.subchain == chain_id:
                ch.add_residue(r)
    if len(ch) == 0:
        for c in st[0]:
            if c.name == chain_id:
                for r in c:
                    ch.add_residue(r)
    if len(ch) == 0:
        avail = sorted({r.subchain for c in st[0] for r in c})[:10]
        raise KeyError(f"Chain {chain_id!r} not found in {src_cif} (subchains: {avail})")

    model = gemmi.Model("1")
    model.add_chain(ch)
    out.add_model(model)
    out.setup_entities()
    out.assign_label_seq_id()
    out.make_mmcif_document().write_file(str(dst_cif))


def prepare_9r2z_ref_pdb(ref_cif: Path, dst_pdb: Path) -> dict[str, str]:
    """Write 9r2z-assembly1 as PDB without SEQRES (A1->A, B1->B) and return reverse map."""
    import gemmi

    st = gemmi.read_structure(str(ref_cif))
    st.remove_ligands_and_waters()
    rev_map = {}
    for ch in st[0]:
        sub = ch[0].subchain
        new_name = {"A1": "A", "B1": "B"}[sub]
        rev_map[new_name] = sub
        ch.name = new_name
    st.write_pdb(str(dst_pdb))
    return rev_map


def run_ost(ost_bin: str, args: list[str], out_json: Path) -> dict:
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    proc = subprocess.run(
        [ost_bin] + args + ["-o", str(out_json)],
        capture_output=True,
        text=True,
        env=env,
        timeout=1800,
    )
    data = json.loads(out_json.read_text(encoding="utf-8")) if out_json.is_file() else {}
    if data.get("status") != "SUCCESS":
        raise RuntimeError(
            f"OST {args[0]} failed (rc={proc.returncode}, status={data.get('status')}): {proc.stderr[-400:]}"
        )
    return data


def score_one_case_prediction(job: dict) -> dict:
    """Score one prediction CIF for a single (task, case, seed, sample)."""
    ost_bin = job["ost_bin"]
    task = job["task"]
    case = job["case"]
    pred_cif = Path(job["prediction_path"])
    ref_cif = Path(job["reference_path"])
    chain_map = job["chain_map"]
    units = job["units"]
    cache_path = Path(job["cache_path"]) if job.get("cache_path") else None

    if cache_path and cache_path.is_file():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("status") == "SUCCESS":
                return cached
        except json.JSONDecodeError:
            pass

    result = {"status": "SUCCESS", "task": task, "case": case, "unit_scores": []}
    with tempfile.TemporaryDirectory(prefix="openbench_poly_") as tmp_s:
        tmp = Path(tmp_s)
        if task in INTERFACE_TASKS:
            extra_flags = SPECIAL_INTERFACE_FLAGS.get(case, [])
            rev_map = None
            ref_arg = str(ref_cif)
            if case == "9r2z-assembly1":
                ref_pdb = tmp / "ref_9r2z.pdb"
                rev_map = prepare_9r2z_ref_pdb(ref_cif, ref_pdb)
                ref_arg = str(ref_pdb)
            ost_data = run_ost(
                ost_bin,
                [
                    "compare-structures",
                    "-m", str(pred_cif),
                    "-r", ref_arg,
                    "--fault-tolerant",
                    "--min-pep-length", "4",
                    "--min-nuc-length", "4",
                    "--lddt",
                    "--rigid-scores",
                    "--tm-score",
                    "--dockq",
                ] + extra_flags,
                tmp / "iface_ost.json",
            )
            if rev_map:
                ost_data["dockq_interfaces"] = [
                    [rev_map[i[0]], rev_map[i[1]], i[2], i[3]]
                    for i in ost_data.get("dockq_interfaces", [])
                ]
            for u in units:
                target_pair = {u["chain_1"], u["chain_2"]}
                matched = None
                for idx, ifc in enumerate(ost_data.get("dockq_interfaces", [])):
                    if {ifc[0], ifc[1]} == target_pair:
                        matched = {
                            "dockq": ost_data["dockq"][idx],
                            "irmsd": ost_data["irmsd"][idx],
                            "lrmsd": ost_data["lrmsd"][idx],
                            "lddt": ost_data.get("lddt"),
                            "tm_score": ost_data.get("tm_score"),
                        }
                        break
                result["unit_scores"].append(matched)
        else:
            for idx, u in enumerate(units):
                ref_chain = u["chain_1"]
                pred_chain = chain_map[ref_chain]
                ref_c_path = tmp / f"ref_{idx}.cif"
                pred_c_path = tmp / f"pred_{idx}.cif"
                extract_single_chain(ref_cif, ref_chain, ref_c_path)
                extract_single_chain(pred_cif, pred_chain, pred_c_path)
                ost_data = run_ost(
                    ost_bin,
                    [
                        "compare-structures",
                        "-m", str(pred_c_path),
                        "-r", str(ref_c_path),
                        "--fault-tolerant",
                        "--min-pep-length", "4",
                        "--min-nuc-length", "4",
                        "--lddt",
                        "--rigid-scores",
                        "--tm-score",
                    ],
                    tmp / f"mono_ost_{idx}.json",
                )
                result["unit_scores"].append({
                    "dockq": None,
                    "irmsd": None,
                    "lrmsd": None,
                    "lddt": ost_data.get("lddt"),
                    "tm_score": ost_data.get("tm_score"),
                    "gdt_ts": ost_data.get("oligo_gdtts"),
                    "rmsd": ost_data.get("rmsd"),
                })

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_out = cache_path.with_suffix(".tmp")
        tmp_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        tmp_out.replace(cache_path)
    return result


def load_units(benchmark_root: Path, tasks: list[str]) -> dict[tuple[str, str], list[dict]]:
    units_by_tc: dict[tuple[str, str], list[dict]] = {}
    for task in tasks:
        for r in read_csv(benchmark_root / "targets" / f"{task}.csv"):
            case = r["pdb_id"]
            if task == "protein-monomer":
                u = {"task": task, "pdb_id": case, "chain_1": r["chain_id"], "chain_2": ""}
            else:
                u = {
                    "task": task,
                    "pdb_id": case,
                    "chain_1": r["interface_chain_id_1"],
                    "chain_2": r["interface_chain_id_2"],
                }
            units_by_tc.setdefault((task, case), []).append(u)
    return units_by_tc


def read_confidence_score(conf_path: Path) -> float:
    if conf_path.suffix == ".npz":
        import numpy as np
        with np.load(conf_path, allow_pickle=False) as data:
            val = float(np.asarray(data["aggregate_score"]).reshape(-1)[0])
    else:
        obj = json.loads(conf_path.read_text(encoding="utf-8"))
        raw = obj.get("ranking_score", obj.get("confidence_score"))
        if isinstance(raw, list):
            raw = raw[0]
        val = float(raw)
    if not math.isfinite(val):
        raise ValueError(f"Non-finite ranking score in {conf_path}")
    return val


def discover_predictions(
    predictions_dir: Path | None,
    manifest_path: Path | None,
    needed_cases: set[str],
) -> dict[str, list[dict]]:
    """Return case -> list of {seed, sample, prediction_path, ranking_score}."""
    by_case: dict[str, list[dict]] = {}
    if manifest_path:
        for r in read_csv(manifest_path):
            case = r["pdb_id"]
            if case not in needed_cases:
                continue
            pred_p = Path(r["prediction_path"])
            if not pred_p.is_file():
                continue
            if r.get("ranking_score") not in (None, ""):
                rank = float(r["ranking_score"])
            else:
                rank = read_confidence_score(Path(r["confidence_path"]))
            by_case.setdefault(case, []).append({
                "seed": int(r["seed"]),
                "sample": int(r["sample"]),
                "prediction_path": str(pred_p),
                "ranking_score": rank,
            })
        return by_case

    assert predictions_dir is not None
    for case in sorted(needed_cases):
        candidates = [
            predictions_dir / case / case,
            predictions_dir / case,
        ]
        for task in POLYMER_TASKS:
            candidates.append(predictions_dir / task / case / case)
            candidates.append(predictions_dir / task / case)
        base = next((d for d in candidates if d.is_dir() and any(d.glob("seed_*"))), None)
        if base is None:
            continue
        for seed in DEFAULT_SEEDS:
            pdir = base / f"seed_{seed}" / "predictions"
            if not pdir.is_dir():
                continue
            for sample in range(5):
                cif = pdir / f"{case}_sample_{sample}.cif"
                if not cif.is_file():
                    cif = pdir / f"{case}_seed_{seed}_sample_{sample}.cif"
                conf = pdir / f"{case}_summary_confidence_sample_{sample}.json"
                if not conf.is_file():
                    conf = pdir / f"{case}_seed_{seed}_summary_confidence_sample_{sample}.json"
                if cif.is_file() and conf.is_file():
                    by_case.setdefault(case, []).append({
                        "seed": seed,
                        "sample": sample,
                        "prediction_path": str(cif),
                        "ranking_score": read_confidence_score(conf),
                    })
    return by_case


def cmd_selftest(args: argparse.Namespace) -> int:
    root = Path(args.benchmark_root).resolve()
    ost_bin = resolve_ost_bin(args.ost_bin)
    test_items = []
    if args.case:
        test_items.append((args.task, args.case))
    else:
        test_items = [
            ("protein-protein", "9r2z-assembly1"),
            ("antibody-antigen", "10gh-assembly1"),
            ("protein-rna", "10yz-assembly1"),
            ("protein-rna", "9w2k-assembly1"),
            ("protein-dna", "9d3p-assembly1"),
            ("protein-monomer", "9ynm-assembly1"),
        ]
    units_map = load_units(root, POLYMER_TASKS)
    for task, case in test_items:
        ref = root / "structures" / task / f"{case}.cif"
        cm_raw = json.loads((root / "chain_maps" / task / f"{case}.json").read_text())["reference_to_input"]
        ident_map = {k: k for k in cm_raw}
        res = score_one_case_prediction({
            "ost_bin": ost_bin,
            "task": task,
            "case": case,
            "prediction_path": str(ref),
            "reference_path": str(ref),
            "chain_map": ident_map,
            "units": units_map[(task, case)],
            "cache_path": None,
        })
        print(f"[SELFTEST OK] {task}/{case}: {res['unit_scores']}")
        for us in res["unit_scores"]:
            if us is None:
                raise RuntimeError(f"Selftest unit returned None for {task}/{case}")
            val = us["lddt"] if task == "protein-monomer" else us["dockq"]
            if val is None or val < 0.99:
                raise RuntimeError(f"Unexpected selftest score < 0.99 for {task}/{case}: {us}")
    return 0


def _safe_score(j: dict):
    try:
        return j, score_one_case_prediction(j), None
    except Exception as exc:
        return j, None, str(exc)


def cmd_run(args: argparse.Namespace) -> int:
    root = Path(args.benchmark_root).resolve()
    out_dir = Path(args.output_dir).resolve()
    ost_bin = resolve_ost_bin(args.ost_bin)
    tasks = POLYMER_TASKS if "all" in args.tasks else args.tasks

    units_map = load_units(root, tasks)
    needed_cases = {case for _, case in units_map}
    preds_by_case = discover_predictions(
        Path(args.predictions_dir).resolve() if args.predictions_dir else None,
        Path(args.manifest).resolve() if args.manifest else None,
        needed_cases,
    )

    jobs = []
    for (task, case), units in sorted(units_map.items()):
        ref_cif = root / "structures" / task / f"{case}.cif"
        cmap = json.loads((root / "chain_maps" / task / f"{case}.json").read_text())["reference_to_input"]
        for p in preds_by_case.get(case, []):
            s, k = p["seed"], p["sample"]
            cache_f = out_dir / "ost_cache" / task / f"{case}__s{s}__n{k}.json"
            jobs.append({
                "ost_bin": ost_bin,
                "task": task,
                "case": case,
                "seed": s,
                "sample": k,
                "ranking_score": p["ranking_score"],
                "prediction_path": p["prediction_path"],
                "reference_path": str(ref_cif),
                "chain_map": cmap,
                "units": units,
                "cache_path": str(cache_f),
            })

    print(f"Discovered {len(jobs)} prediction scoring jobs across {len(units_map)} (task, case) pairs.", flush=True)

    scored_map: dict[tuple[str, str, int, int], dict] = {}
    errors = []

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for j, res, err in pool.map(_safe_score, jobs):
            if err is not None:
                errors.append({"task": j["task"], "case": j["case"], "seed": j["seed"], "sample": j["sample"], "error": err})
            else:
                scored_map[(j["task"], j["case"], j["seed"], j["sample"])] = res

    rank_rows, oracle_rows = [], []
    for (task, case), units in sorted(units_map.items()):
        cand_preds = [
            p for p in preds_by_case.get(case, [])
            if (task, case, p["seed"], p["sample"]) in scored_map
        ]
        top_pred = (
            max(cand_preds, key=lambda x: (x["ranking_score"], -x["seed"], -x["sample"]))
            if cand_preds else None
        )
        for idx, u in enumerate(units):
            base = {
                "task": task,
                "pdb_id": case,
                "chain_id_1": u["chain_1"],
                "chain_id_2": u["chain_2"],
                "n_scored_predictions": len(cand_preds),
            }
            if top_pred is not None:
                m_top = scored_map[(task, case, top_pred["seed"], top_pred["sample"])]["unit_scores"][idx] or {}
                rank_rows.append({
                    **base,
                    "selected_seed": top_pred["seed"],
                    "selected_sample": top_pred["sample"],
                    "ranking_score": top_pred["ranking_score"],
                    "dockq": m_top.get("dockq", ""),
                    "irmsd": m_top.get("irmsd", ""),
                    "lrmsd": m_top.get("lrmsd", ""),
                    "lddt": m_top.get("lddt", ""),
                    "tm_score": m_top.get("tm_score", ""),
                })
            else:
                rank_rows.append({
                    **base,
                    "selected_seed": "",
                    "selected_sample": "",
                    "ranking_score": "",
                    "dockq": "",
                    "irmsd": "",
                    "lrmsd": "",
                    "lddt": "",
                    "tm_score": "",
                })

            unit_cands = [
                scored_map[(task, case, p["seed"], p["sample"])]["unit_scores"][idx]
                for p in cand_preds
            ]
            unit_cands = [c for c in unit_cands if c is not None]

            def best_val(key: str, fn):
                vals = [c[key] for c in unit_cands if c.get(key) is not None]
                return fn(vals) if vals else ""

            oracle_rows.append({
                **base,
                "dockq": best_val("dockq", max),
                "irmsd": best_val("irmsd", min),
                "lrmsd": best_val("lrmsd", min),
                "lddt": best_val("lddt", max),
                "tm_score": best_val("tm_score", max),
            })

    write_csv(out_dir / "units_rank.csv", rank_rows)
    write_csv(out_dir / "units_oracle.csv", oracle_rows)
    if errors:
        (out_dir / "scoring_errors.json").write_text(json.dumps(errors, indent=2) + "\n", encoding="utf-8")

    summary_rows = []
    for mode, rows in (("rank", rank_rows), ("oracle", oracle_rows)):
        for task in tasks:
            sub = [r for r in rows if r["task"] == task]
            n_units = len(sub)
            n_cases = len({r["pdb_id"] for r in sub})
            if task == "protein-monomer":
                lddts = [float(r["lddt"]) if r["lddt"] != "" else 0.0 for r in sub]
                tms = [float(r["tm_score"]) if r["tm_score"] != "" else 0.0 for r in sub]
                summary_rows.append({
                    "mode": mode,
                    "task": task,
                    "cases": n_cases,
                    "units": n_units,
                    "primary_metric": "mean_chain_lddt",
                    "primary_value": round(sum(lddts) / n_units, 4),
                    "dockq_ge_023_pct": "",
                    "dockq_ge_049_pct": "",
                    "dockq_ge_080_pct": "",
                    "mean_tm_score": round(sum(tms) / n_units, 4),
                })
            else:
                dqs = [float(r["dockq"]) if r["dockq"] != "" else 0.0 for r in sub]
                acc = round(100.0 * sum(x >= 0.23 for x in dqs) / n_units, 2)
                med = round(100.0 * sum(x >= 0.49 for x in dqs) / n_units, 2)
                high = round(100.0 * sum(x >= 0.80 for x in dqs) / n_units, 2)
                summary_rows.append({
                    "mode": mode,
                    "task": task,
                    "cases": n_cases,
                    "units": n_units,
                    "primary_metric": "dockq_ge_023_pct",
                    "primary_value": acc,
                    "dockq_ge_023_pct": acc,
                    "dockq_ge_049_pct": med,
                    "dockq_ge_080_pct": high,
                    "mean_tm_score": "",
                })
    write_csv(out_dir / "summary.csv", summary_rows)
    print(f"Wrote {out_dir / 'summary.csv'} ({len(errors)} scoring errors).", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="OpenBench Polymer Evaluation Suite")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_self = sub.add_parser("selftest", help="Verify OST polymer scoring against ground-truth structures")
    p_self.add_argument("--benchmark-root", default=".", help="Path to OpenBench/ root directory")
    p_self.add_argument("--task", choices=POLYMER_TASKS, default="protein-protein")
    p_self.add_argument("--case", default=None, help="Optional single case ID to test")
    p_self.add_argument("--ost-bin", default=None, help="Optional path to `ost` executable")

    p_run = sub.add_parser("run", help="Score and aggregate polymer predictions")
    p_run.add_argument("--benchmark-root", default=".", help="Path to OpenBench/ root directory")
    grp = p_run.add_mutually_exclusive_group(required=True)
    grp.add_argument("--predictions-dir", help="Root directory of model predictions")
    grp.add_argument("--manifest", help="CSV manifest with columns: pdb_id,seed,sample,prediction_path,confidence_path")
    p_run.add_argument("--output-dir", required=True, help="Output directory for CSV results and OST cache")
    p_run.add_argument("--tasks", nargs="+", default=["all"], choices=["all"] + POLYMER_TASKS)
    p_run.add_argument("--workers", type=int, default=16)
    p_run.add_argument("--ost-bin", default=None, help="Optional path to `ost` executable")

    args = parser.parse_args()
    if args.cmd == "selftest":
        return cmd_selftest(args)
    return cmd_run(args)


if __name__ == "__main__":
    sys.exit(main())
