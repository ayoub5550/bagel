#!/usr/bin/env python3
"""01b — Build the triterpenoid structures this plant actually makes.

PubChem does not know most of them, so every structure here is *constructed* from a
PubChem parent (oleanolic acid CID 10494, betulinic acid CID 64971) by explicit,
auditable edits on labelled atoms, then checked against an independently computed
molecular formula (see EXPECTED below).

Why this script exists
----------------------
The published chemistry of *Anabasis articulata* is NOT plain oleanolic-acid
glycosides. Two isolation papers with NMR report 23-oxygenated **30-nor-oleanane**
and **lupane** glycosides:

  * salaheldine2019 — Planta Med 85:274-281, doi:10.1055/a-0762-0885
  * gamal2022      — Nat Prod Res 36:4076-4084, doi:10.1080/14786419.2021.1961769

Compound names are copied verbatim from those abstracts (read 2026-10-07).

Honest limits (kept in the output column `stereo_note`)
-------------------------------------------------------
1. C-4 carries two methyls (C-23, C-24). The edit oxidises one of them; which one the
   SMILES stereochemistry makes "C-23" is NOT resolved here. Constitution, formula and
   exact mass are unaffected; 3D geometry marginally so.
2. For 3beta,20alpha-dihydroxy-30-nor-olean-12-ene-23,28-dioic acid the C-20
   configuration is left unspecified (reported as 20alpha).
3. Sugar stereochemistry comes from the PubChem pyranose parents (beta-D-Glc, -Gal,
   -Xyl, -GlcA); the papers' anomeric/linkage assignments are accepted as written.

Output: data/01b_triterpenoids.csv , QC picture: data/qc/triterpenoids.png
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdMolDescriptors

ROOT = Path(__file__).resolve().parents[1]
OUT_CSV = ROOT / "data" / "01b_triterpenoids.csv"
QC_PNG = ROOT / "data" / "qc" / "triterpenoids.png"

SITE = "site"

# ---------------------------------------------------------------- parents (PubChem)
OLEANOLIC = "C[C@]12CC[C@@H](C([C@@H]1CC[C@@]3([C@@H]2CC=C4[C@]3(CC[C@@]5([C@H]4CC(CC5)(C)C)C(=O)O)C)C)(C)C)O"  # CID 10494
BETULINIC = "CC(=C)[C@@H]1CC[C@]2([C@H]1[C@H]3CC[C@@H]4[C@]5(CC[C@@H](C([C@@H]5CC[C@]4([C@@]3(CC2)C)C)(C)C)O)C)C(=O)O"  # CID 64971

# Glycosyl donors, written from the anomeric carbon (atom 0), taken from the PubChem
# pyranose parents: beta-D-glucopyranose CID 64689, beta-D-galactopyranose CID 439353,
# beta-D-xylopyranose CID 125409, beta-D-glucopyranuronic acid CID 94715.
SUGARS = {  # `*` marks where the glycosidic bond is formed (anomeric carbon)
    "Glc": "*[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O",
    "Gal": "*[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O",
    "Xyl": "*[C@@H]1OC[C@@H](O)[C@H](O)[C@H]1O",
    "GlcA": "*[C@@H]1O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]1O",
}


# ---------------------------------------------------------------- labelling
def _atom(mol: Chem.Mol, site: str) -> int:
    hits = [a.GetIdx() for a in mol.GetAtoms()
            if a.HasProp(SITE) and a.GetProp(SITE) == site]
    if len(hits) != 1:
        raise RuntimeError(f"site {site!r}: expected 1 atom, found {len(hits)}")
    return hits[0]


def _has(mol: Chem.Mol, site: str) -> bool:
    return any(a.HasProp(SITE) and a.GetProp(SITE) == site for a in mol.GetAtoms())


def label_parent(smiles: str, skeleton: str) -> Chem.Mol:
    """Label C-3/O-3, C-4 and its two methyls, the C-28 acid, and (oleananes) C-20."""
    mol = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts("[OX2H]-[CX4H1]-[CX4]([CH3])([CH3])")
    hits = mol.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise RuntimeError(f"expected one 3-OH/4,4-dimethyl motif, found {len(hits)}")
    o3, c3, c4, me_a, me_b = hits[0]
    mol.GetAtomWithIdx(o3).SetProp(SITE, "O3")
    mol.GetAtomWithIdx(c3).SetProp(SITE, "C3")
    mol.GetAtomWithIdx(c4).SetProp(SITE, "C4")
    mol.GetAtomWithIdx(me_a).SetProp(SITE, "C23")  # the methyl we oxidise
    mol.GetAtomWithIdx(me_b).SetProp(SITE, "C24")

    acids = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
    if len(acids) != 1:
        raise RuntimeError(f"expected one COOH in the parent, found {len(acids)}")
    c28, _o, o28 = acids[0]
    mol.GetAtomWithIdx(c28).SetProp(SITE, "C28")
    mol.GetAtomWithIdx(o28).SetProp(SITE, "O28H")

    if skeleton == "oleanane":
        found = []
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "C" or not atom.IsInRing():
                continue
            if atom.GetIdx() == c4:
                continue
            methyls = [a.GetIdx() for a in atom.GetNeighbors()
                       if a.GetSymbol() == "C" and a.GetTotalNumHs() == 3]
            if len(methyls) == 2:
                found.append((atom.GetIdx(), methyls))
        if len(found) != 1:
            raise RuntimeError(f"expected one oleanane C-20, found {len(found)}")
        c20, (me_c, me_d) = found[0]
        mol.GetAtomWithIdx(c20).SetProp(SITE, "C20")
        mol.GetAtomWithIdx(me_c).SetProp(SITE, "C29")
        mol.GetAtomWithIdx(me_d).SetProp(SITE, "C30")  # the methyl removed by 30-nor
    return mol


# ---------------------------------------------------------------- edits
def oxidise_c23(mol: Chem.Mol, to: str) -> Chem.Mol:
    """C-23 methyl -> aldehyde ('al') or carboxylic acid ('oic')."""
    rw = Chem.RWMol(mol)
    c23 = _atom(rw, "C23")
    o_dbl = rw.AddAtom(Chem.Atom(8))
    rw.GetAtomWithIdx(o_dbl).SetProp(SITE, "O23")
    rw.AddBond(c23, o_dbl, Chem.BondType.DOUBLE)
    if to == "oic":
        o_h = rw.AddAtom(Chem.Atom(8))
        rw.GetAtomWithIdx(o_h).SetProp(SITE, "O23H")
        rw.AddBond(c23, o_h, Chem.BondType.SINGLE)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def nor30_ene20(mol: Chem.Mol) -> Chem.Mol:
    """30-nor + 20(29)-ene: drop C-30, double-bond C-20=C-29 (exocyclic =CH2)."""
    rw = Chem.RWMol(mol)
    c30 = _atom(rw, "C30")
    rw.RemoveAtom(c30)
    c20, c29 = _atom(rw, "C20"), _atom(rw, "C29")
    rw.GetAtomWithIdx(c20).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rw.GetBondBetweenAtoms(c20, c29).SetBondType(Chem.BondType.DOUBLE)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def nor30_hydroxy20(mol: Chem.Mol) -> Chem.Mol:
    """30-nor + 20-OH: drop C-30, hydroxylate C-20 (configuration left unspecified)."""
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(_atom(rw, "C30"))
    c20 = _atom(rw, "C20")
    rw.GetAtomWithIdx(c20).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    o = rw.AddAtom(Chem.Atom(8))
    rw.GetAtomWithIdx(o).SetProp(SITE, "O20H")
    rw.AddBond(c20, o, Chem.BondType.SINGLE)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def _attach_sugar(mol: Chem.Mol, oxygen_site: str, sugar: str) -> Chem.Mol:
    """Bond the anomeric carbon of `sugar` to the oxygen labelled `oxygen_site`.

    The donor fragment carries a dummy atom (`*`) on its anomeric carbon so the
    carbon never passes through a radical state: the real bond is made first, the
    dummy removed afterwards, and only then is the molecule sanitised.
    """
    rw = Chem.RWMol(mol)
    oxy = _atom(rw, oxygen_site)
    frag = Chem.MolFromSmiles(SUGARS[sugar])
    offset = rw.GetNumAtoms()
    rw.InsertMol(frag)
    dummy = next(a.GetIdx() for a in rw.GetAtoms()
                 if a.GetIdx() >= offset and a.GetAtomicNum() == 0)
    anomeric = rw.GetAtomWithIdx(dummy).GetNeighbors()[0].GetIdx()
    rw.AddBond(oxy, anomeric, Chem.BondType.SINGLE)
    rw.RemoveAtom(dummy)
    rw.GetAtomWithIdx(_atom(rw, oxygen_site)).ClearProp(SITE)  # no longer a free OH
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    Chem.AssignStereochemistry(out, cleanIt=True, force=True)
    return out


def glycosylate_c3(mol: Chem.Mol, sugar: str) -> Chem.Mol:
    """3-O-glycoside (ether at the C-3 hydroxyl)."""
    return _attach_sugar(mol, "O3", sugar)


def esterify(mol: Chem.Mol, sugar: str, acid: str) -> Chem.Mol:
    """Glycosyl ester of the C-28 acid (acid='C28') or the C-23 acid (acid='C23')."""
    site = {"C28": "O28H", "C23": "O23H"}[acid]
    if not _has(mol, site):
        raise RuntimeError(f"no free {acid} acid to esterify")
    return _attach_sugar(mol, site, sugar)


def smi(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


# ---------------------------------------------------------------- the inventory
def build() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []

    def add(cid, name, mol, skeleton, kind, source_key, evidence, steps, note=""):
        rows.append({
            "id": cid, "literature_name": name, "skeleton": skeleton, "kind": kind,
            "source_key": source_key, "evidence": evidence, "smiles": smi(mol),
            "construction": steps, "stereo_note": note,
        })

    ole = label_parent(OLEANOLIC, "oleanane")
    bet = label_parent(BETULINIC, "lupane")

    ole_23al = oxidise_c23(ole, "al")
    ole_23oic = oxidise_c23(ole, "oic")
    agl2 = nor30_ene20(ole_23al)                 # 23-al, 30-nor, 12,20(29)-diene
    agl3 = nor30_ene20(ole)                      # 30-nor, 12,20(29)-diene
    agl4 = nor30_hydroxy20(ole_23oic)            # 23-oic, 30-nor, 20-OH
    agl5 = oxidise_c23(bet, "al")                # lupane 23-al
    agl6 = oxidise_c23(bet, "oic")               # lupane 23,28-dioic

    C23_NOTE = "C-23 vs C-24 methyl choice UNVERIFIED"

    add("AGL-1", "oleanolic acid", ole, "oleanane", "aglycone", "metwally2012",
        "isolated", "PubChem CID 10494, unmodified")
    add("AGL-2", "3beta-hydroxy-23-aldehyde-30-norolean-12,20(29)-dien-28-oic acid",
        agl2, "30-nor-oleanane", "aglycone", "salaheldine2019;gamal2022", "isolated",
        "CID 10494 -> C-23 to CHO -> 30-nor + 20(29)-ene", C23_NOTE)
    add("AGL-3", "30-norolean-12,20(29)-dien-28-oic acid (aglycone of compound 3)",
        agl3, "30-nor-oleanane", "aglycone", "salaheldine2019", "isolated",
        "CID 10494 -> 30-nor + 20(29)-ene")
    add("AGL-4", "3beta,20alpha-dihydroxy-30-nor-olean-12-ene-23,28-dioic acid",
        agl4, "30-nor-oleanane", "aglycone", "gamal2022", "isolated",
        "CID 10494 -> C-23 to COOH -> 30-nor + 20-OH",
        "C-20 configuration unspecified (reported 20alpha); " + C23_NOTE)
    add("AGL-5", "3beta-hydroxy-23-aldehyde-lup-20(29)-en-28-oic acid", agl5, "lupane",
        "aglycone", "gamal2022", "isolated",
        "betulinic acid CID 64971 -> C-23 to CHO", C23_NOTE)
    add("AGL-6", "3beta-hydroxy-lup-20(29)-ene-23,28-dioic acid", agl6, "lupane",
        "aglycone", "gamal2022", "isolated",
        "betulinic acid CID 64971 -> C-23 to COOH", C23_NOTE)
    add("AGL-7", "betulinic acid (lupane reference parent)", bet, "lupane", "reference",
        "-", "not_reported_in_species", "PubChem CID 64971, unmodified")

    add("SAP-1", "oleanolic acid 3-O-beta-D-glucopyranoside",
        glycosylate_c3(ole, "Glc"), "oleanane", "monodesmoside", "metwally2012",
        "isolated", "AGL-1 -> 3-O-Glc")
    add("SAP-2", "3-O-beta-D-glucopyranosyl-28-O-beta-D-xylopyranosyl oleanolic acid",
        esterify(glycosylate_c3(ole, "Glc"), "Xyl", "C28"), "oleanane", "bidesmoside",
        "metwally2012", "isolated", "AGL-1 -> 3-O-Glc -> 28-O-Xyl ester")
    add("SAP-3", "3beta-hydroxy-23-aldehyde-30-norolean-12,20(29)-dien-28-oic acid "
        "28-O-beta-D-glucopyranosyl ester (compound 1)",
        esterify(agl2, "Glc", "C28"), "30-nor-oleanane", "monodesmoside",
        "salaheldine2019", "isolated", "AGL-2 -> 28-O-Glc ester",
        "PPARalpha activation 2.25-fold reported for this compound")
    add("SAP-4", "3beta-O-D-galactopyranosyl-23-aldehyde-30-norolean-12,20(29)-dien-"
        "28-oic acid 28-O-beta-D-glucopyranosyl ester (compound 2)",
        esterify(glycosylate_c3(agl2, "Gal"), "Glc", "C28"), "30-nor-oleanane",
        "bidesmoside", "salaheldine2019", "isolated",
        "AGL-2 -> 3-O-Gal -> 28-O-Glc ester", C23_NOTE)
    add("SAP-5", "3beta-O-D-xylopyranosyl-30-norolean-12,20(29)-dien-28-oic acid "
        "28-O-beta-D-glucopyranosyl ester (compound 3)",
        esterify(glycosylate_c3(agl3, "Xyl"), "Glc", "C28"), "30-nor-oleanane",
        "bidesmoside", "salaheldine2019", "isolated",
        "AGL-3 -> 3-O-Xyl -> 28-O-Glc ester",
        "PPARalpha activation 1.86-fold reported for this compound")
    add("SAP-6", "3beta-hydroxy-lup-20(29)-ene-23,28-dioic acid 23-O-beta-D-"
        "glucopyranosyl ester (compound 5)", esterify(agl6, "Glc", "C23"), "lupane",
        "monodesmoside", "gamal2022", "isolated", "AGL-6 -> 23-O-Glc ester", C23_NOTE)
    add("SAP-7", "3-O-beta-D-glucuronopyranosyl-lup-20(29)-ene-23,28-dioic acid "
        "28-O-beta-D-glucopyranosyl ester (compound 6)",
        esterify(glycosylate_c3(agl6, "GlcA"), "Glc", "C28"), "lupane", "bidesmoside",
        "gamal2022", "isolated", "AGL-6 -> 3-O-GlcA -> 28-O-Glc ester", C23_NOTE)
    add("SAP-8", "3-O-beta-D-glucuronopyranosyl-lup-20(29)-ene-23-aldehyde-28-oic acid "
        "28-O-beta-D-glucopyranosyl ester (compound 7)",
        esterify(glycosylate_c3(agl5, "GlcA"), "Glc", "C28"), "lupane", "bidesmoside",
        "gamal2022", "isolated", "AGL-5 -> 3-O-GlcA -> 28-O-Glc ester", C23_NOTE)

    return rows


EXPECTED = {  # independent hand arithmetic, checked against RDKit's formula
    "AGL-1": "C30H48O3", "AGL-2": "C29H42O4", "AGL-3": "C29H44O3", "AGL-4": "C29H44O6",
    "AGL-5": "C30H46O4", "AGL-6": "C30H46O5", "AGL-7": "C30H48O3",
    "SAP-1": "C36H58O8", "SAP-2": "C41H66O12", "SAP-3": "C35H52O9", "SAP-4": "C41H62O14",
    "SAP-5": "C40H62O12", "SAP-6": "C36H56O10", "SAP-7": "C42H64O16", "SAP-8": "C42H64O15",
}


def main() -> int:
    rows = build()
    bad = 0
    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None:
            print(f"[01b] FAIL unparsable SMILES for {row['id']}")
            bad += 1
            continue
        formula = rdMolDescriptors.CalcMolFormula(mol)
        row["formula"] = formula
        row["exact_mass"] = round(Descriptors.ExactMolWt(mol), 4)
        row["mw"] = round(Descriptors.MolWt(mol), 2)
        row["rings"] = rdMolDescriptors.CalcNumRings(mol)
        row["n_stereocentres"] = len(Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False))
        want = EXPECTED.get(str(row["id"]))
        row["formula_check"] = "ok" if formula == want else f"MISMATCH want {want}"
        if formula != want:
            print(f"[01b] FAIL {row['id']}: formula {formula} != expected {want}")
            bad += 1

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "literature_name", "skeleton", "kind", "source_key", "evidence",
            "formula", "exact_mass", "mw", "rings", "n_stereocentres", "smiles",
            "construction", "stereo_note", "formula_check"]
    with OUT_CSV.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)

    QC_PNG.parent.mkdir(parents=True, exist_ok=True)
    mols, legends = [], []
    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        AllChem.Compute2DCoords(mol)
        mols.append(mol)
        legends.append(f"{row['id']} {row['formula']}")
    Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(430, 390), legends=legends,
                         returnPNG=False).save(QC_PNG)

    print(f"[01b] built {len(rows)} structures, {bad} formula failures")
    for row in rows:
        print(f"       {row['id']:7s} {row['formula']:12s} {row['exact_mass']:>9.4f}  "
              f"{row['formula_check']:10s} {str(row['literature_name'])[:58]}")
    print(f"[01b] wrote {OUT_CSV.relative_to(ROOT)} and {QC_PNG.relative_to(ROOT)}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
