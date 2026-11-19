#!/usr/bin/env python3
"""OpenBench Protein-Ligand Evaluation Suite (Runs N' Poses Protocol).

Implements the exact Runs N' Poses evaluation protocol (Skrinjar et al., 2025, Methods 5.2-5.4)
over the 952 protein-ligand systems (1,129 proper ligands, 514 SuCOS-pocket clusters):

  1. Prediction Preparation (SMILES -> Receptor CIF + Ligand SDFs):
     - Splits each predicted complex mmCIF into `receptor.cif` (polymer chains) plus one
       `<chain_id>.sdf` per ligand chain reconstructed from the input SMILES in
       `inputs/protein-ligand/<system_id>.json` and the predicted heavy-atom coordinates.
     - Normalizes SDF dative bonds (bond type 9 -> single bond 1) in temporary copies of
       both predicted and reference ligand SDFs so OpenStructure 2.8.0 parses metallo-ligands.
     - Extracts `ranking_score` and per-ligand-chain `lig_prot_iptm` (mean ligand->protein
       chain-pair iPTM across protein chains, Runs N' Poses Methods 5.2).

  2. Accuracy & Physical Plausibility Scoring (Methods 5.3):
     - Runs OpenStructure `compare-ligand-structures` with flags:
         `--fault-tolerant --lddt-pli --lddt-pli-add-mdl-contacts --rmsd`
       against `structures/protein-ligand/<system_id>/receptor.cif` and `ligand_files/*.sdf`.
     - Runs PoseBusters (`PoseBusters(config="dock")`) on each assigned proper ligand SDF
       conditioned on the predicted receptor structure (`pb_valid` = all default checks pass).

  3. Aggregation & Stratified Reporting (Methods 5.3 & 5.4):
     - Primary success criterion: `rmsd < 2.0` A AND `lddt_pli > 0.8`
     - PoseBusters-valid success: `rmsd < 2.0` A AND `lddt_pli > 0.8` AND `pb_valid`
     - Primary Top-1 selection (`top_iptm`): highest `lig_prot_iptm` per proper ligand across
       the 25 predictions (5 seeds x 5 samples); also reports `top_ranking_score`, `best_of_25`,
       and `mean_of_25` with 1,000-sample bootstrap 95% confidence intervals.
     - Stratified across the 8 SuCOS-pocket `similarity_bin` buckets (`0-20` .. `80-100`)
       over both All proper ligands (1,129) and Cluster representatives (514).
"""
import argparse
import csv
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

SIMILARITY_BINS = ["0-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-100"]
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
    for candidate in ("/usr/local/bin/ost", "/opt/conda/envs/ost/bin/ost"):
        if Path(candidate).is_file():
            return candidate
    raise RuntimeError("Could not find `ost` executable on PATH.")


def rewrite_sdf_dative_bonds(src_sdf: Path, dst_sdf: Path) -> bool:
    """Rewrite V2000/V3000 SDF bond type 9 (dative) to 1 (single) for OST 2.8.0 compatibility."""
    lines = src_sdf.read_text(encoding="utf-8", errors="replace").splitlines()
    changed = False
    if len(lines) > 3 and "V2000" in lines[3]:
        na, nb = int(lines[3][:3]), int(lines[3][3:6])
        for j in range(nb):
            idx = 4 + na + j
            bl = lines[idx]
            if len(bl) >= 9 and int(bl[6:9]) == 9:
                lines[idx] = bl[:6] + "  1" + bl[9:]
                changed = True
    elif len(lines) > 3 and "V3000" in lines[3]:
        in_bond = False
        for idx, l in enumerate(lines):
            if "M  V30 BEGIN BOND" in l:
                in_bond = True
                continue
            if "M  V30 END BOND" in l:
                in_bond = False
                continue
            if in_bond and l.startswith("M  V30 "):
                parts = l[7:].split()
                if len(parts) >= 4 and parts[1] == "9":
                    parts[1] = "1"
                    lines[idx] = "M  V30 " + " ".join(parts)
                    changed = True
    if changed:
        dst_sdf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    else:
        shutil.copy2(src_sdf, dst_sdf)
    return changed


def load_input_ligand_smiles(input_json: Path) -> dict[str, str]:
    raw = json.loads(input_json.read_text(encoding="utf-8"))
    item = raw[0] if isinstance(raw, list) else raw
    mapping = {}
    for seq in item["sequences"]:
        if "ligand" in seq:
            lig = seq["ligand"]
            for cid in lig["id"]:
                mapping[cid] = lig["ligand"]
    return mapping


def parse_confidence(conf_path: Path, chains: list[str], prot_chains: list[str], lig_chains: list[str]) -> tuple[float, dict[str, float]]:
    import numpy as np

    if conf_path.suffix == ".npz":
        with np.load(conf_path, allow_pickle=False) as data:
            rank = float(np.asarray(data["aggregate_score"]).reshape(-1)[0])
            mat = np.asarray(data["per_chain_pair_iptm"], dtype=float)
            if mat.ndim == 3:
                mat = mat[0]
    else:
        c = json.loads(conf_path.read_text(encoding="utf-8"))
        if "pair_chains_iptm" in c:
            pm = c["pair_chains_iptm"]
            n = len(pm)
            mat = np.array([[pm[str(i)][str(j)] for j in range(n)] for i in range(n)], dtype=float)
            raw_rank = c["confidence_score"]
        else:
            mat = np.array(c["chain_pair_iptm"], dtype=float)
            raw_rank = c["ranking_score"]
        if isinstance(raw_rank, list):
            raw_rank = raw_rank[0]
        rank = float(raw_rank)

    idx = {ch: i for i, ch in enumerate(chains)}
    prot_idx = [idx[p] for p in prot_chains if p in idx]
    lig_iptm = {}
    for lc in lig_chains:
        li = idx.get(lc)
        if li is not None and prot_idx and mat.shape[0] > max([li] + prot_idx):
            lig_iptm[lc] = float(np.mean(mat[li, prot_idx]))
    return rank, lig_iptm


def split_predicted_complex(
    pred_cif: Path,
    smiles_by_chain: dict[str, str],
    out_dir: Path,
) -> tuple[Path, Path, list[Path], list[str], list[str], list[str], list[str]]:
    """Split predicted mmCIF into receptor.cif, receptor.pdb, and SMILES-templated <chain>.sdf files."""
    import gemmi
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    st = gemmi.read_structure(str(pred_cif))
    st.setup_entities()
    st.remove_hydrogens()

    all_chains = [ch.name for ch in st[0]]
    lig_chains = [ch.name for ch in st[0] if ch.name in smiles_by_chain]
    prot_chains = [ch.name for ch in st[0] if ch.name not in smiles_by_chain]

    pred_sdfs = []
    problems = []
    for cname in lig_chains:
        ch = st[0][cname]
        atoms = [a for res in ch for a in res]
        mol = Chem.MolFromSmiles(smiles_by_chain[cname])
        if mol is None:
            problems.append(f"{cname}: RDKit failed to parse SMILES")
            continue
        el_m = [a.GetSymbol().upper() for a in mol.GetAtoms()]
        el_p = [a.element.name.upper() for a in atoms]
        elem_mismatch = len(el_m) != len(el_p) or any(
            m != p and not (len(el_m) == 1 or (len(m) == 2 and a.name.upper().startswith(m)))
            for m, p, a in zip(el_m, el_p, atoms)
        )
        if elem_mismatch:
            problems.append(f"{cname}: heavy-atom order mismatch ({len(el_m)} vs {len(el_p)})")
            continue

        conf3d = Chem.Conformer(mol.GetNumAtoms())
        for i, a in enumerate(atoms):
            conf3d.SetAtomPosition(i, Point3D(a.pos.x, a.pos.y, a.pos.z))
        mol.AddConformer(conf3d, assignId=True)
        mol.SetProp("_Name", cname)
        for b in mol.GetBonds():
            if b.GetBondType() == Chem.BondType.DATIVE:
                b.SetBondType(Chem.BondType.SINGLE)

        sdf_path = out_dir / f"{cname}.sdf"
        w = Chem.SDWriter(str(sdf_path))
        w.write(mol)
        w.close()
        pred_sdfs.append(sdf_path)

    for cname in lig_chains:
        st[0].remove_chain(cname)
    st.setup_entities()

    rec_cif = out_dir / "receptor.cif"
    rec_pdb = out_dir / "receptor.pdb"
    st.make_mmcif_document().write_file(str(rec_cif))
    st_pdb = st.clone()
    st_pdb.shorten_chain_names()
    st_pdb.write_pdb(str(rec_pdb))
    return rec_cif, rec_pdb, pred_sdfs, all_chains, prot_chains, lig_chains, problems


def run_posebusters_dock(pred_sdf: Path, receptor_pdb: Path) -> bool | None:
    """Run PoseBusters `dock` suite and return True iff all default checks pass."""
    try:
        from posebusters import PoseBusters
    except ImportError:
        return None
    buster = PoseBusters(config="dock")
    df = buster.bust([str(pred_sdf)], None, str(receptor_pdb))
    bool_df = df.select_dtypes(include=["bool"])
    if bool_df.empty:
        return False
    return bool(bool_df.iloc[0].all())


def score_one_pl_prediction(job: dict) -> dict:
    ost_bin = job["ost_bin"]
    sid = job["system_id"]
    pred_cif = Path(job["prediction_path"])
    conf_path = Path(job["confidence_path"]) if job.get("confidence_path") else None
    input_json = Path(job["input_json_path"])
    gt_dir = Path(job["gt_dir"])
    proper_ligands = set(job["proper_ligands"])
    check_pb = job["run_posebusters"]
    cache_path = Path(job["cache_path"]) if job.get("cache_path") else None

    if cache_path and cache_path.is_file():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("status") == "SUCCESS":
                return cached
        except json.JSONDecodeError:
            pass

    with tempfile.TemporaryDirectory(prefix="openbench_pl_") as tmp_s:
        tmp = Path(tmp_s)
        prep_dir = tmp / "prep"
        ref_lig_dir = tmp / "ref_ligs"
        prep_dir.mkdir()
        ref_lig_dir.mkdir()

        ref_sdfs = []
        for s in sorted((gt_dir / "ligand_files").glob("*.sdf")):
            dst = ref_lig_dir / s.name
            rewrite_sdf_dative_bonds(s, dst)
            ref_sdfs.append(str(dst))

        if job.get("selftest"):
            rec_cif = gt_dir / "receptor.cif"
            import gemmi
            st = gemmi.read_structure(str(rec_cif))
            st.shorten_chain_names()
            rec_pdb = prep_dir / "receptor.pdb"
            st.write_pdb(str(rec_pdb))
            pred_sdfs = [Path(p) for p in ref_sdfs]
            ranking_score = 1.0
            lig_iptm = {Path(p).stem: 1.0 for p in ref_sdfs}
            problems = []
        else:
            smiles_map = load_input_ligand_smiles(input_json)
            rec_cif, rec_pdb, pred_sdfs, all_ch, prot_ch, lig_ch, problems = split_predicted_complex(
                pred_cif, smiles_map, prep_dir
            )
            if conf_path and conf_path.is_file():
                ranking_score, lig_iptm = parse_confidence(conf_path, all_ch, prot_ch, lig_ch)
            else:
                ranking_score = float(job.get("ranking_score") or 0.0)
                lig_iptm = {}

        ligand_rows = []
        if pred_sdfs:
            ost_out = tmp / "ost.json"
            cmd = [
                ost_bin,
                "compare-ligand-structures",
                "--fault-tolerant",
                "-m", str(rec_cif),
                "-ml", *[str(p) for p in pred_sdfs],
                "-r", str(gt_dir / "receptor.cif"),
                "-rl", *ref_sdfs,
                "-o", str(ost_out),
                "--lddt-pli",
                "--lddt-pli-add-mdl-contacts",
                "--rmsd",
            ]
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
            proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1800)
            ost_data = json.loads(ost_out.read_text(encoding="utf-8")) if ost_out.is_file() else {}
            if ost_data.get("status") != "SUCCESS":
                raise RuntimeError(
                    f"OST compare-ligand-structures failed on {sid} (rc={proc.returncode}): {proc.stderr[-400:]}"
                )

            lp_by_ref = {
                Path(x["reference_ligand"]).stem: x
                for x in ost_data.get("lddt_pli", {}).get("assigned_scores", [])
            }
            rm_by_ref = {
                Path(x["reference_ligand"]).stem: x
                for x in ost_data.get("rmsd", {}).get("assigned_scores", [])
            }
            pred_sdf_by_name = {p.stem: p for p in pred_sdfs}

            for lig in sorted(proper_ligands):
                if lig not in lp_by_ref and lig not in rm_by_ref:
                    continue
                assigned = rm_by_ref.get(lig) or lp_by_ref.get(lig)
                mdl_name = Path(assigned["model_ligand"]).stem
                rmsd_val = rm_by_ref[lig]["score"] if lig in rm_by_ref else None
                lddt_lp_val = rm_by_ref[lig].get("lddt_lp") if lig in rm_by_ref else None
                lddt_pli_val = lp_by_ref[lig]["score"] if lig in lp_by_ref else None

                pb_valid = None
                if check_pb and mdl_name in pred_sdf_by_name:
                    pb_valid = run_posebusters_dock(pred_sdf_by_name[mdl_name], rec_pdb)

                ligand_rows.append({
                    "ligand_chain_id": lig,
                    "model_ligand_chain": mdl_name,
                    "rmsd": rmsd_val,
                    "lddt_pli": lddt_pli_val,
                    "lddt_lp": lddt_lp_val,
                    "pb_valid": pb_valid,
                    "lig_prot_iptm": lig_iptm.get(mdl_name),
                    "ranking_score": ranking_score,
                })

    result = {
        "status": "SUCCESS",
        "system_id": sid,
        "ranking_score": ranking_score,
        "problems": problems,
        "ligands": ligand_rows,
    }
    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        tmp_out = cache_path.with_suffix(".tmp")
        tmp_out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        tmp_out.replace(cache_path)
    return result


def bootstrap_ci(vals: list[float], n_boot: int = 1000) -> tuple[float, float]:
    if not vals:
        return 0.0, 0.0
    rng = random.Random(0)
    n = len(vals)
    means = sorted(sum(rng.choice(vals) for _ in range(n)) / n for _ in range(n_boot))
    return means[int(0.025 * n_boot)], means[int(0.975 * n_boot)]


def discover_pl_predictions(
    predictions_dir: Path | None,
    manifest_path: Path | None,
    needed_systems: set[str],
) -> dict[str, list[dict]]:
    by_sys: dict[str, list[dict]] = {}
    if manifest_path:
        for r in read_csv(manifest_path):
            sid = r["pdb_id"]
            if sid not in needed_systems:
                continue
            pred_p = Path(r["prediction_path"])
            if not pred_p.is_file():
                continue
            by_sys.setdefault(sid, []).append({
                "seed": int(r["seed"]),
                "sample": int(r["sample"]),
                "prediction_path": str(pred_p),
                "confidence_path": r.get("confidence_path", ""),
                "ranking_score": float(r["ranking_score"]) if r.get("ranking_score") not in (None, "") else None,
            })
        return by_sys

    assert predictions_dir is not None
    for sid in sorted(needed_systems):
        candidates = [
            predictions_dir / sid / sid,
            predictions_dir / sid,
            predictions_dir / "protein-ligand" / sid / sid,
            predictions_dir / "protein-ligand" / sid,
        ]
        base = next((d for d in candidates if d.is_dir() and any(d.glob("seed_*"))), None)
        if base is None:
            continue
        for seed in DEFAULT_SEEDS:
            pdir = base / f"seed_{seed}" / "predictions"
            if not pdir.is_dir():
                continue
            for sample in range(5):
                cif = pdir / f"{sid}_sample_{sample}.cif"
                if not cif.is_file():
                    cif = pdir / f"{sid}_seed_{seed}_sample_{sample}.cif"
                conf = pdir / f"{sid}_summary_confidence_sample_{sample}.json"
                if not conf.is_file():
                    conf = pdir / f"{sid}_seed_{seed}_summary_confidence_sample_{sample}.json"
                if cif.is_file() and conf.is_file():
                    by_sys.setdefault(sid, []).append({
                        "seed": seed,
                        "sample": sample,
                        "prediction_path": str(cif),
                        "confidence_path": str(conf),
                        "ranking_score": None,
                    })
    return by_sys


def cmd_selftest(args: argparse.Namespace) -> int:
    root = Path(args.benchmark_root).resolve()
    ost_bin = resolve_ost_bin(args.ost_bin)
    targets = read_csv(root / "targets" / "protein-ligand.csv")
    proper_by_sys: dict[str, list[str]] = defaultdict(list)
    for r in targets:
        proper_by_sys[r["pdb_id"]].append(r["ligand_chain_id"])

    test_sids = [args.system] if args.system else ["10aj__1__1.A__1.C"]
    for sid in test_sids:
        res = score_one_pl_prediction({
            "ost_bin": ost_bin,
            "system_id": sid,
            "prediction_path": "",
            "confidence_path": "",
            "input_json_path": str(root / "inputs" / "protein-ligand" / f"{sid}.json"),
            "gt_dir": str(root / "structures" / "protein-ligand" / sid),
            "proper_ligands": proper_by_sys[sid],
            "run_posebusters": not args.skip_posebusters,
            "cache_path": None,
            "selftest": True,
        })
        print(f"[SELFTEST OK] protein-ligand/{sid}: {res['ligands']}")
        if len(res["ligands"]) != len(proper_by_sys[sid]):
            raise RuntimeError(f"Expected {len(proper_by_sys[sid])} ligands, got {len(res['ligands'])}")
        for lr in res["ligands"]:
            if lr["rmsd"] is None or lr["rmsd"] > 1e-3 or lr["lddt_pli"] is None or lr["lddt_pli"] < 0.99:
                raise RuntimeError(f"Unexpected selftest ligand metrics for {sid}: {lr}")
    return 0


def _safe_score(j: dict):
    try:
        return j, score_one_pl_prediction(j), None
    except Exception as exc:
        return j, None, str(exc)


def cmd_run(args: argparse.Namespace) -> int:
    root = Path(args.benchmark_root).resolve()
    out_dir = Path(args.output_dir).resolve()
    ost_bin = resolve_ost_bin(args.ost_bin)

    targets = read_csv(root / "targets" / "protein-ligand.csv")
    proper_by_sys: dict[str, list[str]] = defaultdict(list)
    for r in targets:
        proper_by_sys[r["pdb_id"]].append(r["ligand_chain_id"])

    needed_systems = set(proper_by_sys.keys())
    preds_by_sys = discover_pl_predictions(
        Path(args.predictions_dir).resolve() if args.predictions_dir else None,
        Path(args.manifest).resolve() if args.manifest else None,
        needed_systems,
    )

    jobs = []
    for sid in sorted(needed_systems):
        for p in preds_by_sys.get(sid, []):
            s, k = p["seed"], p["sample"]
            cache_f = out_dir / "ost_cache" / f"{sid}__s{s}__n{k}.json"
            jobs.append({
                "ost_bin": ost_bin,
                "system_id": sid,
                "seed": s,
                "sample": k,
                "prediction_path": p["prediction_path"],
                "confidence_path": p["confidence_path"],
                "ranking_score": p["ranking_score"],
                "input_json_path": str(root / "inputs" / "protein-ligand" / f"{sid}.json"),
                "gt_dir": str(root / "structures" / "protein-ligand" / sid),
                "proper_ligands": proper_by_sys[sid],
                "run_posebusters": not args.skip_posebusters,
                "cache_path": str(cache_f),
                "selftest": False,
            })

    print(f"Discovered {len(jobs)} protein-ligand prediction scoring jobs across {len(needed_systems)} systems.", flush=True)

    per_pred_rows = []
    errors = []

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for j, res, err in pool.map(_safe_score, jobs):
            if err is not None:
                errors.append({"system_id": j["system_id"], "seed": j["seed"], "sample": j["sample"], "error": err})
                continue
            for lr in res["ligands"]:
                rmsd = lr["rmsd"]
                lpli = lr["lddt_pli"]
                pb = lr["pb_valid"]
                succ = bool(rmsd is not None and lpli is not None and rmsd < 2.0 and lpli > 0.8)
                succ_rmsd = bool(rmsd is not None and rmsd < 2.0)
                succ_pb = bool(succ and pb is True) if pb is not None else ""
                per_pred_rows.append({
                    "pdb_id": j["system_id"],
                    "ligand_chain_id": lr["ligand_chain_id"],
                    "seed": j["seed"],
                    "sample": j["sample"],
                    "model_ligand_chain": lr["model_ligand_chain"],
                    "rmsd": rmsd if rmsd is not None else "",
                    "lddt_pli": lpli if lpli is not None else "",
                    "lddt_lp": lr["lddt_lp"] if lr["lddt_lp"] is not None else "",
                    "pb_valid": pb if pb is not None else "",
                    "lig_prot_iptm": lr["lig_prot_iptm"] if lr["lig_prot_iptm"] is not None else "",
                    "ranking_score": lr["ranking_score"],
                    "success_rnp": succ,
                    "success_rmsd_lt2": succ_rmsd,
                    "success_rnp_pb": succ_pb,
                })

    if per_pred_rows:
        write_csv(out_dir / "per_prediction_ligand.csv", per_pred_rows)
    if errors:
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "scoring_errors.json").write_text(json.dumps(errors, indent=2) + "\n", encoding="utf-8")

    by_unit: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for r in per_pred_rows:
        by_unit[(r["pdb_id"], r["ligand_chain_id"])].append(r)

    has_pb = any(r["pb_valid"] != "" for r in per_pred_rows)
    selection_rows = []
    for t in targets:
        sid = t["pdb_id"]
        lch = t["ligand_chain_id"]
        cands = by_unit.get((sid, lch), [])
        if cands:
            top_iptm = max(
                cands,
                key=lambda x: (
                    float(x["lig_prot_iptm"]) if x["lig_prot_iptm"] != "" else float(x["ranking_score"]),
                    -int(x["seed"]),
                    -int(x["sample"]),
                ),
            )
            top_rank = max(cands, key=lambda x: (float(x["ranking_score"]), -int(x["seed"]), -int(x["sample"])))
            selection_rows.append({
                "pdb_id": sid,
                "ligand_chain_id": lch,
                "ligand_ccd_code": t["ligand_ccd_code"],
                "similarity_bin": t["similarity_bin"],
                "is_cluster_representative": t["is_cluster_representative"],
                "n_scored_predictions": len(cands),
                "top_iptm_success": int(bool(top_iptm["success_rnp"])),
                "top_iptm_rmsd_lt2": int(bool(top_iptm["success_rmsd_lt2"])),
                "top_iptm_rnp_pb": int(bool(top_iptm["success_rnp_pb"])) if has_pb else "",
                "top_iptm_pb_valid": int(top_iptm["pb_valid"] is True) if has_pb else "",
                "top_ranking_success": int(bool(top_rank["success_rnp"])),
                "best_of_25_success": int(any(bool(c["success_rnp"]) for c in cands)),
                "mean_of_25_success": round(sum(1.0 if c["success_rnp"] else 0.0 for c in cands) / 25.0, 4),
            })
        else:
            selection_rows.append({
                "pdb_id": sid,
                "ligand_chain_id": lch,
                "ligand_ccd_code": t["ligand_ccd_code"],
                "similarity_bin": t["similarity_bin"],
                "is_cluster_representative": t["is_cluster_representative"],
                "n_scored_predictions": 0,
                "top_iptm_success": 0,
                "top_iptm_rmsd_lt2": 0,
                "top_iptm_rnp_pb": 0 if has_pb else "",
                "top_iptm_pb_valid": 0 if has_pb else "",
                "top_ranking_success": 0,
                "best_of_25_success": 0,
                "mean_of_25_success": 0.0,
            })

    write_csv(out_dir / "per_ligand_selection.csv", selection_rows)

    subsets = [
        ("all_proper_ligands", selection_rows),
        ("cluster_representatives", [r for r in selection_rows if r["is_cluster_representative"].lower() == "true"]),
    ]
    overall_rows = []
    bin_rows = []
    for subset_name, sub in subsets:
        n = len(sub)
        vals = [float(r["top_iptm_success"]) for r in sub]
        lo, hi = bootstrap_ci(vals)
        overall_rows.append({
            "subset": subset_name,
            "n_ligands": n,
            "top_iptm_rnp_success_pct": round(100.0 * sum(vals) / n, 2),
            "top_iptm_ci95_low": round(100.0 * lo, 2),
            "top_iptm_ci95_high": round(100.0 * hi, 2),
            "top_iptm_rmsd_lt2_pct": round(100.0 * sum(float(r["top_iptm_rmsd_lt2"]) for r in sub) / n, 2),
            "top_iptm_rnp_and_pb_pct": round(100.0 * sum(float(r["top_iptm_rnp_pb"]) for r in sub) / n, 2) if has_pb else "",
            "top_iptm_pb_valid_pct": round(100.0 * sum(float(r["top_iptm_pb_valid"]) for r in sub) / n, 2) if has_pb else "",
            "top_ranking_rnp_success_pct": round(100.0 * sum(float(r["top_ranking_success"]) for r in sub) / n, 2),
            "best_of_25_rnp_success_pct": round(100.0 * sum(float(r["best_of_25_success"]) for r in sub) / n, 2),
            "mean_of_25_rnp_success_pct": round(100.0 * sum(float(r["mean_of_25_success"]) for r in sub) / n, 2),
            "n_unscored_ligands": sum(1 for r in sub if int(r["n_scored_predictions"]) == 0),
        })
        for b in SIMILARITY_BINS:
            bsub = [r for r in sub if r["similarity_bin"] == b]
            bn = len(bsub)
            bvals = [float(r["top_iptm_success"]) for r in bsub]
            blo, bhi = bootstrap_ci(bvals)
            bin_rows.append({
                "subset": subset_name,
                "similarity_bin": b,
                "n_ligands": bn,
                "top_iptm_rnp_success_pct": round(100.0 * sum(bvals) / bn, 2) if bn else "",
                "top_iptm_ci95_low": round(100.0 * blo, 2) if bn else "",
                "top_iptm_ci95_high": round(100.0 * bhi, 2) if bn else "",
                "top_iptm_rnp_and_pb_pct": round(100.0 * sum(float(r["top_iptm_rnp_pb"]) for r in bsub) / bn, 2) if (bn and has_pb) else "",
                "best_of_25_rnp_success_pct": round(100.0 * sum(float(r["best_of_25_success"]) for r in bsub) / bn, 2) if bn else "",
            })

    write_csv(out_dir / "summary_overall.csv", overall_rows)
    write_csv(out_dir / "summary_by_similarity_bin.csv", bin_rows)
    print(f"Wrote {out_dir / 'summary_overall.csv'} and {out_dir / 'summary_by_similarity_bin.csv'} ({len(errors)} scoring errors).", flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="OpenBench Protein-Ligand Evaluation Suite (Runs N' Poses Protocol)")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_self = sub.add_parser("selftest", help="Verify OpenStructure + PoseBusters scoring on ground truth")
    p_self.add_argument("--benchmark-root", default=".", help="Path to OpenBench/ root directory")
    p_self.add_argument("--system", default=None, help="Optional PLINDER system_id to test")
    p_self.add_argument("--skip-posebusters", action="store_true", help="Skip PoseBusters check")
    p_self.add_argument("--ost-bin", default=None, help="Optional path to `ost` executable")

    p_run = sub.add_parser("run", help="Score and aggregate protein-ligand predictions (Runs N' Poses protocol)")
    p_run.add_argument("--benchmark-root", default=".", help="Path to OpenBench/ root directory")
    grp = p_run.add_mutually_exclusive_group(required=True)
    grp.add_argument("--predictions-dir", help="Root directory of model predictions")
    grp.add_argument("--manifest", help="CSV manifest with columns: pdb_id,seed,sample,prediction_path,confidence_path")
    p_run.add_argument("--output-dir", required=True, help="Output directory for CSV results and OST cache")
    p_run.add_argument("--workers", type=int, default=16)
    p_run.add_argument("--skip-posebusters", action="store_true", help="Skip PoseBusters `dock` check")
    p_run.add_argument("--ost-bin", default=None, help="Optional path to `ost` executable")

    args = parser.parse_args()
    if args.cmd == "selftest":
        return cmd_selftest(args)
    return cmd_run(args)


if __name__ == "__main__":
    sys.exit(main())
