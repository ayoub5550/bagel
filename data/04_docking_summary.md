# 04 — docking with positive controls

AutoDock Vina 1.2.5, exhaustiveness 24, 9 modes, seed 20261007, 16 CPUs. Receptors prepared with Open Babel at pH 7.4, ligands with Meeko from RDKit ETKDGv3 + MMFF94 geometries.

## Control check (this decides whether the numbers below mean anything)

| target | protein | control ligand | control score | RMSD to crystal pose (in place) | aligned RMSD (diagnostic only) | verdict |
|---|---|---|---|---|---|---|
| 4GQR | human pancreatic alpha-amylase | MYC | -7.598 kcal/mol | 2.29 A | 0.81 A | **valid** |
| 2P54 | human PPARalpha ligand-binding domain | 735 | -10.67 kcal/mol | 1.39 A | 1.17 A | **valid** |

Acceptance limit: the re-docked crystal ligand must land within 2.5 A of its crystallographic pose, measured **in place** (rdMolAlign.CalcRMS, no superposition). Until 2026-10-09 this table used the aligned RMSD (GetBestRMS), which superimposes the pose on the crystal ligand first and so cannot see a pose docked in the wrong place; the aligned value is shown for comparison only (CHANGELOG 2026-10-09).

## 4GQR — human pancreatic alpha-amylase

Why this target: crude extract inhibits alpha-amylase (aljoufi2022); myricetin co-crystal.

| ligand | class | Vina (kcal/mol) | heavy atoms | ligand efficiency (-kcal/mol per heavy atom) | rot. bonds | note |
|---|---|---|---|---|---|---|
| AGL-3 30-nor-oleanadienoic acid | plant triterpenoid aglycone | -9.935 | 32 | 0.31 | 1 |  |
| SAP-3 nor-saponin (PPARa 2.25-fold) | plant saponin | -9.666 | 44 | 0.22 | 4 |  |
| AGL-2 23-al-30-nor-oleanadienoic acid | plant triterpenoid aglycone | -9.311 | 33 | 0.282 | 2 |  |
| AGL-1 oleanolic acid | plant triterpenoid aglycone | -9.307 | 33 | 0.282 | 1 |  |
| AGL-4 20-OH-30-nor-oleanene-23,28-dioic acid | plant triterpenoid aglycone | -9.059 | 35 | 0.259 | 2 |  |
| SAP-6 lupane 23-O-Glc ester | plant saponin | -8.939 | 46 | 0.194 | 5 |  |
| AGL-6 lupene-23,28-dioic acid | plant triterpenoid aglycone | -8.86 | 35 | 0.253 | 3 |  |
| AGL-5 23-al-lupenoic acid | plant triterpenoid aglycone | -8.452 | 34 | 0.249 | 3 |  |
| catechin | plant phenolic | -8.435 | 21 | 0.402 | 1 |  |
| quercetin | plant phenolic | -8.196 | 22 | 0.373 | 1 |  |
| kaempferol | plant phenolic | -8.09 | 21 | 0.385 | 1 |  |
| myricetin | plant phenolic | -7.744 | 23 | 0.337 | 1 |  |
| chlorogenic acid | plant phenolic | -7.689 | 25 | 0.308 | 4 |  |
| CONTROL MYC (crystal ligand, re-docked) | positive control | -7.598 | 23 | 0.33 | 1 |  |
| REF fenofibric acid (PPARa drug) | reference drug | -7.425 | 22 | 0.337 | 5 |  |
| REF acarbose (amylase drug) | reference drug | -7.36 | 55 | 0.134 | 12 | 12 rotatable bonds: Vina score not interpretable |
| REF GW7647 (PPARa agonist) | reference drug | -6.894 | 29 | 0.238 | 13 | 13 rotatable bonds: Vina score not interpretable |
| gallic acid | plant phenolic | -5.365 | 12 | 0.447 | 1 |  |
| NEG betaine (osmolyte, expect no fit) | negative control | -3.186 | 8 | 0.398 | 2 |  |

## 2P54 — human PPARalpha ligand-binding domain

Why this target: pure saponins of this plant activate PPARalpha (salaheldine2019).

| ligand | class | Vina (kcal/mol) | heavy atoms | ligand efficiency (-kcal/mol per heavy atom) | rot. bonds | note |
|---|---|---|---|---|---|---|
| CONTROL 735 (crystal ligand, re-docked) | positive control | -10.67 | 33 | 0.323 | 7 |  |
| quercetin | plant phenolic | -7.812 | 22 | 0.355 | 1 |  |
| kaempferol | plant phenolic | -7.808 | 21 | 0.372 | 1 |  |
| myricetin | plant phenolic | -7.655 | 23 | 0.333 | 1 |  |
| chlorogenic acid | plant phenolic | -7.61 | 25 | 0.304 | 4 |  |
| REF GW7647 (PPARa agonist) | reference drug | -7.569 | 29 | 0.261 | 13 | 13 rotatable bonds: Vina score not interpretable |
| REF fenofibric acid (PPARa drug) | reference drug | -7.349 | 22 | 0.334 | 5 |  |
| catechin | plant phenolic | -7.172 | 21 | 0.342 | 1 |  |
| gallic acid | plant phenolic | -5.635 | 12 | 0.47 | 1 |  |
| REF acarbose (amylase drug) | reference drug | -4.541 | 55 | 0.083 | 12 | 12 rotatable bonds: Vina score not interpretable |
| NEG betaine (osmolyte, expect no fit) | negative control | -3.823 | 8 | 0.478 | 2 |  |
| SAP-6 lupane 23-O-Glc ester | plant saponin | -2.111 | 46 | 0.046 | 5 |  |
| AGL-5 23-al-lupenoic acid | plant triterpenoid aglycone | -1.527 | 34 | 0.045 | 3 |  |
| AGL-6 lupene-23,28-dioic acid | plant triterpenoid aglycone | -0.4386 | 35 | 0.013 | 3 |  |
| SAP-3 nor-saponin (PPARa 2.25-fold) | plant saponin | 0.6501 | 44 | -0.015 | 4 |  |
| AGL-4 20-OH-30-nor-oleanene-23,28-dioic acid | plant triterpenoid aglycone | 0.8901 | 35 | -0.025 | 2 |  |
| AGL-3 30-nor-oleanadienoic acid | plant triterpenoid aglycone | 1.112 | 32 | -0.035 | 1 |  |
| AGL-2 23-al-30-nor-oleanadienoic acid | plant triterpenoid aglycone | 2.043 | 33 | -0.062 | 2 |  |
| AGL-1 oleanolic acid | plant triterpenoid aglycone | 2.125 | 33 | -0.064 | 1 |  |

## How to read this

- A Vina score is a **ranking**, with roughly +-2 kcal/mol of error. It is not an affinity and it is not evidence of activity.
- Ligands with more than 10 rotatable bonds (the saponins, acarbose) are listed for completeness only; their scores are dominated by conformational search noise.
- The negative control (betaine) should score clearly worse than the real ligands. If it does not, the pocket definition is wrong, not the chemistry.
- A good score for a plant compound means exactly one thing: it is worth measuring. The measurement is an enzyme assay, and it costs little.
- **Compare ligand efficiency, not raw score.** Vina's function rewards size: a bulky triterpene beats a small flavonoid on raw score almost by construction. Efficiency (score divided by heavy-atom count) is the size-corrected number.
- Re-running this script reuses the saved poses in `work/dock/` instead of re-docking. Delete that folder to force a clean run.
