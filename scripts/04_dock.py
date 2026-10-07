#!/usr/bin/env python3
"""04 — Docking with a positive control, or no docking at all.

Project rule (AGENTS.md #2): a docking table without a **positive control** — the
crystal ligand re-docked into its own site — is rejected. This script therefore
re-docks the co-crystal ligand first and refuses to publish scores for a target whose
control pose is worse than CONTROL_RMSD_LIMIT.

Targets (both human, both chosen because they are what the *measured* literature points
at, not because they score well):

  4GQR  human pancreatic alpha-amylase + myricetin, 1.2 A
        -> the plant's crude extract inhibits alpha-amylase (IC50 34 ug/mL, aljoufi2022)
           and myricetin is reported in the plant (makhlouf2024, LC-MS, tentative).
  2P54  human PPARalpha LBD + agonist GW735, 1.79 A
        -> two pure saponins of this plant activate PPARalpha 2.25- and 1.86-fold in
           HepG2 cells (salaheldine2019). That is the only pure-compound molecular
           target datum that exists for this species.

What this script cannot do: turn a Vina score into an affinity. Vina's error is about
+-2 kcal/mol and it is a *ranking* tool. Glycosides with many rotatable bonds are
docked only to show how unreliable they are, and are flagged.

Outputs: data/04_docking.csv, data/04_docking_summary.md, work/dock/** (ignored by git)
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

import numpy as np
from openbabel import pybel
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolDescriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")
pybel.ob.obErrorLog.SetOutputLevel(0)

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work" / "dock"
CACHE = ROOT / "data" / "cache" / "pdb"
OUT_CSV = ROOT / "data" / "04_docking.csv"
OUT_MD = ROOT / "data" / "04_docking_summary.md"
VINA = ROOT / "bin" / "vina"

EXHAUSTIVENESS = 24
N_MODES = 9
SEED = 20261007
CPUS = 16
CONTROL_RMSD_LIMIT = 2.5  # A; above this the target's scores are not reported

TARGETS = {
    "4GQR": {
        "name": "human pancreatic alpha-amylase",
        "ligand_code": "MYC",
        "chain": "A",
        "keep_het": ("CA", "CL"),
        "why": "crude extract inhibits alpha-amylase (aljoufi2022); myricetin co-crystal",
    },
    "2P54": {
        "name": "human PPARalpha ligand-binding domain",
        "ligand_code": "735",
        "chain": "A",
        "keep_het": (),
        "why": "pure saponins of this plant activate PPARalpha (salaheldine2019)",
    },
}

# Ligands: plant compounds, reference actives, and a deliberate negative control.
LIGANDS = {
    # --- plant triterpenoid aglycones (from 01b) -----------------------------
    "AGL-1 oleanolic acid": "01b",
    "AGL-2 23-al-30-nor-oleanadienoic acid": "01b",
    "AGL-3 30-nor-oleanadienoic acid": "01b",
    "AGL-4 20-OH-30-nor-oleanene-23,28-dioic acid": "01b",
    "AGL-5 23-al-lupenoic acid": "01b",
    "AGL-6 lupene-23,28-dioic acid": "01b",
    # --- plant saponins (flagged: too flexible for reliable docking) ---------
    "SAP-3 nor-saponin (PPARa 2.25-fold)": "01b",
    "SAP-6 lupane 23-O-Glc ester": "01b",
    # --- plant phenolics (PubChem) ------------------------------------------
    "quercetin": "Oc1cc(O)c2c(=O)c(O)c(-c3ccc(O)c(O)c3)oc2c1",
    "myricetin": "Oc1cc(O)c2c(=O)c(O)c(-c3cc(O)c(O)c(O)c3)oc2c1",
    "kaempferol": "Oc1ccc(-c2oc3cc(O)cc(O)c3c(=O)c2O)cc1",
    "catechin": "Oc1cc(O)c2c(c1)O[C@@H](c1ccc(O)c(O)c1)[C@H](O)C2",
    "chlorogenic acid": "O[C@@H]1C[C@](O)(C(=O)O)C[C@H](OC(=O)/C=C/c2ccc(O)c(O)c2)[C@@H]1O",
    "gallic acid": "OC(=O)c1cc(O)c(O)c(O)c1",
    # --- reference compounds -------------------------------------------------
    "REF acarbose (amylase drug)": "C[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O[C@H]4[C@H](O)[C@@H](O)C(O)O[C@@H]4CO)O[C@@H]3CO)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1N[C@H]1C=C(CO)[C@@H](O)[C@H](O)[C@H]1O",
    "REF fenofibric acid (PPARa drug)": "CC(C)(Oc1ccc(C(=O)c2ccc(Cl)cc2)cc1)C(=O)O",
    "REF GW7647 (PPARa agonist)": "CCCCCC(CCCCC)NC(=O)Nc1ccc(SC(C)(C)C(=O)O)c(C)c1",
    # --- negative control ----------------------------------------------------
    "NEG betaine (osmolyte, expect no fit)": "C[N+](C)(C)CC(=O)[O-]",
}


def http_get(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": "bagel-dry-lab/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    return dest


def component_smiles(code: str) -> str:
    """Canonical SMILES of a PDB chemical component, straight from RCSB."""
    path = CACHE / f"comp_{code}.json"
    http_get(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}", path)
    data = json.loads(path.read_text())
    for desc in data.get("pdbx_chem_comp_descriptor", []):
        if desc.get("type") == "SMILES_CANONICAL" and desc.get("program") == "OpenEye OEToolkits":
            return desc["descriptor"]
    for desc in data.get("pdbx_chem_comp_descriptor", []):
        if desc.get("type", "").startswith("SMILES"):
            return desc["descriptor"]
    raise RuntimeError(f"no SMILES for component {code}")


def split_pdb(pdb_id: str, spec: dict, outdir: Path) -> tuple[Path, Path]:
    src = http_get(f"https://files.rcsb.org/download/{pdb_id}.pdb", CACHE / f"{pdb_id}.pdb")
    rec_lines, lig_lines = [], []
    for line in src.read_text().splitlines(keepends=True):
        if line.startswith("ATOM") and line[21] == spec["chain"]:
            rec_lines.append(line)
        elif line.startswith("HETATM"):
            res = line[17:20].strip()
            if res == spec["ligand_code"] and line[21] == spec["chain"]:
                lig_lines.append(line)
            elif res in spec["keep_het"] and line[21] == spec["chain"]:
                rec_lines.append(line)
    rec = outdir / f"{pdb_id}_receptor.pdb"
    lig = outdir / f"{pdb_id}_ligand_crystal.pdb"
    rec.write_text("".join(rec_lines) + "END\n")
    lig.write_text("".join(lig_lines) + "END\n")
    return rec, lig


def receptor_pdbqt(rec_pdb: Path) -> Path:
    out = rec_pdb.with_suffix(".pdbqt")
    if out.exists():
        return out
    mol = next(pybel.readfile("pdb", str(rec_pdb)))
    mol.OBMol.CorrectForPH(7.4)
    mol.addh()
    mol.write("pdbqt", str(out), overwrite=True, opt={"r": None})
    return out


def ligand_pdbqt(name: str, smiles: str, outdir: Path) -> tuple[Path, Chem.Mol] | tuple[None, None]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[04] SMILES parse failed: {name}")
        return None, None
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        print(f"[04] 3D embedding failed: {name}")
        return None, None
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    setups = MoleculePreparation().prepare(mol)
    written = PDBQTWriterLegacy.write_string(setups[0])
    text = written[0] if isinstance(written, tuple) else written
    path = outdir / (re.sub(r"[^A-Za-z0-9]+", "_", name)[:50] + ".pdbqt")
    path.write_text(text)
    return path, mol


def box_from_ligand(lig_pdb: Path) -> tuple[np.ndarray, np.ndarray]:
    coords = []
    for line in lig_pdb.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    arr = np.array(coords)
    center = arr.mean(axis=0)
    size = np.clip(arr.max(axis=0) - arr.min(axis=0) + 10.0, 22.0, 30.0)
    return center, size


def run_vina(rec: Path, lig: Path, center, size, out: Path) -> tuple[float | None, Path | None]:
    cmd = [str(VINA), "--receptor", str(rec), "--ligand", str(lig),
           "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}",
           "--center_z", f"{center[2]:.3f}",
           "--size_x", f"{size[0]:.1f}", "--size_y", f"{size[1]:.1f}",
           "--size_z", f"{size[2]:.1f}",
           "--exhaustiveness", str(EXHAUSTIVENESS), "--num_modes", str(N_MODES),
           "--seed", str(SEED), "--cpu", str(CPUS), "--out", str(out)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"[04] vina failed for {lig.name}: {proc.stderr.strip()[:300]}")
        return None, None
    best = None
    for line in proc.stdout.splitlines():
        m = re.match(r"\s*1\s+(-?\d+\.\d+)", line)
        if m:
            best = float(m.group(1))
            break
    return best, out


def control_rmsd(template_smiles: str, crystal_pdb: Path, docked_pdbqt: Path) -> float | None:
    """Symmetry-aware RMSD between the re-docked best pose and the crystal pose."""
    try:
        template = Chem.MolFromSmiles(template_smiles)
        cryst = Chem.MolFromPDBBlock(crystal_pdb.read_text(), removeHs=True,
                                     proximityBonding=True)
        cryst = AllChem.AssignBondOrdersFromTemplate(template, cryst)

        pose_pdb = docked_pdbqt.with_suffix(".model1.pdb")
        model: list[str] = []
        for line in docked_pdbqt.read_text().splitlines():
            if line.startswith("ENDMDL"):
                break
            if line.startswith(("ATOM", "HETATM")):
                model.append(line[:66] + "\n")
        pose_pdb.write_text("".join(model) + "END\n")
        probe = Chem.MolFromPDBBlock(pose_pdb.read_text(), removeHs=True,
                                     proximityBonding=True)
        probe = AllChem.AssignBondOrdersFromTemplate(template, probe)
        return float(AllChem.GetBestRMS(probe, cryst))
    except Exception as err:  # noqa: BLE001
        print(f"[04] control RMSD not computable: {err}")
        return None


def main() -> int:
    if not VINA.exists():
        print("[04] bin/vina missing — run scripts/setup_lab.sh first")
        return 1
    WORK.mkdir(parents=True, exist_ok=True)

    # ligand set: SMILES from this file + the constructed triterpenoids
    triterp = {}
    tri_csv = ROOT / "data" / "01b_triterpenoids.csv"
    if tri_csv.exists():
        with tri_csv.open() as fh:
            triterp = {r["id"]: r for r in csv.DictReader(fh)}

    rows: list[dict[str, object]] = []
    summary: dict[str, dict] = {}

    for pdb_id, spec in TARGETS.items():
        tdir = WORK / pdb_id
        tdir.mkdir(parents=True, exist_ok=True)
        rec_pdb, lig_pdb = split_pdb(pdb_id, spec, tdir)
        rec = receptor_pdbqt(rec_pdb)
        center, size = box_from_ligand(lig_pdb)
        ref_smiles = component_smiles(spec["ligand_code"])
        print(f"[04] {pdb_id} {spec['name']}: box centre "
              f"{np.round(center, 2).tolist()} size {np.round(size, 1).tolist()}")
        print(f"[04] {pdb_id} control ligand {spec['ligand_code']} = {ref_smiles}")

        # ---- positive control: re-dock the crystal ligand --------------------
        ctrl_lig, _ = ligand_pdbqt(f"CONTROL_{spec['ligand_code']}", ref_smiles, tdir)
        ctrl_score, ctrl_out = (None, None)
        if ctrl_lig:
            ctrl_score, ctrl_out = run_vina(rec, ctrl_lig, center, size,
                                            tdir / "control_out.pdbqt")
        rmsd = control_rmsd(ref_smiles, lig_pdb, ctrl_out) if ctrl_out else None
        valid = rmsd is not None and rmsd <= CONTROL_RMSD_LIMIT
        summary[pdb_id] = {"name": spec["name"], "why": spec["why"],
                           "control_code": spec["ligand_code"],
                           "control_score": ctrl_score, "control_rmsd": rmsd,
                           "valid": valid, "center": center.tolist(),
                           "size": size.tolist()}
        print(f"[04] {pdb_id} CONTROL score {ctrl_score} kcal/mol, RMSD to crystal "
              f"{rmsd if rmsd is None else round(rmsd, 2)} A -> "
              f"{'VALID' if valid else 'INVALID (scores will be withheld)'}")
        rows.append({
            "target": pdb_id, "target_name": spec["name"],
            "ligand": f"CONTROL {spec['ligand_code']} (crystal ligand, re-docked)",
            "ligand_class": "positive control", "vina_kcal_mol": ctrl_score,
            "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(
                Chem.MolFromSmiles(ref_smiles)),
            "control_rmsd_A": None if rmsd is None else round(rmsd, 2),
            "target_run_valid": valid, "flag": "",
        })

        # ---- the real ligands ------------------------------------------------
        for name, source in LIGANDS.items():
            if source == "01b":
                key = name.split()[0]
                if key not in triterp:
                    print(f"[04] missing {key} in 01b output, skipped")
                    continue
                smiles = triterp[key]["smiles"]
                lclass = ("plant triterpenoid aglycone" if key.startswith("AGL")
                          else "plant saponin")
            else:
                smiles = source
                lclass = ("reference drug" if name.startswith("REF")
                          else "negative control" if name.startswith("NEG")
                          else "plant phenolic")
            lig, mol = ligand_pdbqt(name, smiles, tdir)
            if lig is None:
                continue
            rotb = rdMolDescriptors.CalcNumRotatableBonds(Chem.MolFromSmiles(smiles))
            score, _ = run_vina(rec, lig, center, size,
                                tdir / (lig.stem + "_out.pdbqt"))
            flag = ""
            if rotb > 10:
                flag = f"{rotb} rotatable bonds: Vina score not interpretable"
            print(f"[04] {pdb_id} {name[:46]:46s} {score if score is not None else 'FAIL':>8} "
                  f"kcal/mol  {flag}")
            rows.append({
                "target": pdb_id, "target_name": spec["name"], "ligand": name,
                "ligand_class": lclass, "vina_kcal_mol": score, "rot_bonds": rotb,
                "control_rmsd_A": None if rmsd is None else round(rmsd, 2),
                "target_run_valid": valid, "flag": flag,
            })

    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    write_summary(rows, summary)
    print(f"[04] wrote {OUT_CSV.relative_to(ROOT)} and {OUT_MD.relative_to(ROOT)}")
    return 0


def write_summary(rows, summary) -> None:
    lines = [
        "# 04 — docking with positive controls",
        "",
        f"AutoDock Vina 1.2.5, exhaustiveness {EXHAUSTIVENESS}, {N_MODES} modes, "
        f"seed {SEED}, {CPUS} CPUs. Receptors prepared with Open Babel at pH 7.4, "
        "ligands with Meeko from RDKit ETKDGv3 + MMFF94 geometries.",
        "",
        "## Control check (this decides whether the numbers below mean anything)",
        "",
        "| target | protein | control ligand | control score | RMSD to crystal pose | verdict |",
        "|---|---|---|---|---|---|",
    ]
    for pdb_id, s in summary.items():
        rmsd = "not computable" if s["control_rmsd"] is None else f"{s['control_rmsd']:.2f} A"
        lines.append(f"| {pdb_id} | {s['name']} | {s['control_code']} | "
                     f"{s['control_score']} kcal/mol | {rmsd} | "
                     f"{'**valid**' if s['valid'] else '**INVALID — scores withheld**'} |")
    lines += [
        "",
        f"Acceptance limit: the re-docked crystal ligand must land within "
        f"{CONTROL_RMSD_LIMIT} A of its crystallographic pose.",
        "",
    ]
    for pdb_id, s in summary.items():
        lines += [f"## {pdb_id} — {s['name']}", "", f"Why this target: {s['why']}.", ""]
        if not s["valid"]:
            lines += [
                "**The positive control failed, so no score for this target is reported "
                "here.** Per AGENTS.md rule 2 the table is withheld rather than "
                "published with a caveat. The raw numbers remain in "
                "`data/04_docking.csv` with `target_run_valid=False` for anyone who "
                "wants to see why the run was discarded.",
                "",
            ]
            continue
        lines += ["| ligand | class | Vina (kcal/mol) | rot. bonds | note |",
                  "|---|---|---|---|---|"]
        sel = [r for r in rows if r["target"] == pdb_id and r["vina_kcal_mol"] is not None]
        for r in sorted(sel, key=lambda r: float(r["vina_kcal_mol"])):
            lines.append(f"| {r['ligand']} | {r['ligand_class']} | "
                         f"{r['vina_kcal_mol']} | {r['rot_bonds']} | {r['flag']} |")
        lines.append("")
    lines += [
        "## How to read this",
        "",
        "- A Vina score is a **ranking**, with roughly +-2 kcal/mol of error. It is not "
        "an affinity and it is not evidence of activity.",
        "- Ligands with more than 10 rotatable bonds (the saponins, acarbose) are listed "
        "for completeness only; their scores are dominated by conformational search "
        "noise.",
        "- The negative control (betaine) should score clearly worse than the real "
        "ligands. If it does not, the pocket definition is wrong, not the chemistry.",
        "- A good score for a plant compound means exactly one thing: it is worth "
        "measuring. The measurement is an enzyme assay, and it costs little.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines))


if __name__ == "__main__":
    sys.exit(main())
