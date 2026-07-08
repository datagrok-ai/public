---
feature: chem
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md
migration_date: 2026-05-11
related_bugs: []
manual_only_reason: |
  Chem | Calculate menu walk (Descriptors / Properties / Toxicity Risks /
  MPO / BitBIRCH / Cluster MCS / Map Identifiers / InChI / InChI Keys /
  Biochemical Properties) requires the Chem-side Docker container
  (`descriptors`, etc.) to be provisioned and running for several of the
  calculators. On Playwright dev runs the container is NOT available, so
  the Chem | Calculate sub-menu items either dispatch errors or hang on
  backend calls. The UI walk remains valuable as manual regression
  coverage. Reclassified manual-only 2026-05-12 per user directive after
  the B2 batch modernization attempt surfaced the Docker dependency.
  Renamed calculate.md → calculate-ui.md per filename-convention rule
  (-ui.md suffix encodes manual-only).
---

# Chem | Calculate menu walk

Multi-format / multi-Calculate-menu-item happy-path walk. For each dataset in the linked-datasets
matrix (smiles, molV2000, molV3000, 2-column SMILES variants), iterate through the full
**Chem > Calculate** top-menu and verify each calculator appends its expected column(s) to the
active table view. The Descriptors calculator is the primary case, since it is backed by the Chem
Docker container; `smiles.csv` is the representative dataset.

## Setup

1. **Provision linked datasets.** All four files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
   - `System:AppData/Chem/tests/smiles_2_columns.csv` — 2-column SMILES dataset (Variant D —
     two molecule columns; exercises column-selection branches in calculator dialogs).
2. **Confirm Chem package is loaded** so that the full **Chem > Calculate** top-menu tree is
   registered (package source `Chem/src/package.ts`). For Descriptors specifically, the `chem-chem`
   Docker container must be available — failure mode bug-library `GROK-17621` (Docker timeout) is
   outside scope; the happy path assumes container is up.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each dataset is opened fresh inside
   the scenario.

## Scenarios

### Calculate menu walk per dataset

Data-driven walk: for each dataset (D1–D4) and each Calculate menu item (M1–M10), open the
dataset, run the menu item dialog, accept arbitrary defaults, click OK, verify the calculator
appends its expected column(s) without console errors.

Dataset matrix:

| id | dataset path | format |
|----|--------------|--------|
| D1 | `System:DemoFiles/chem/smiles.csv` | smiles |
| D2 | `System:AppData/Chem/mol1K.sdf` | molV2000 |
| D3 | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` | molV3000 |
| D4 | `System:AppData/Chem/tests/smiles_2_columns.csv` | smiles, 2-column |

Calculate menu items (run in this order on each dataset):

| id | top-menu path | added column / view (expected) |
|----|---------------|--------------------------------|
| M1 | `Chem > Calculate > Descriptors...` | one or more descriptor columns appended to the grid (selection-dependent; e.g. `MolWt`, `LogP`, etc.) — Docker-backed |
| M2 | `Chem > Calculate > Chemical Properties...` | one or more property columns appended (e.g. `MW`, `HBA`, `HBD`, `LogP`, `LogS`, `PSA`, `Rotatable Bonds`, `Stereo Centers`, `Molecule Charge` — selection-dependent) |
| M3 | `Chem > Calculate > Toxicity Risks...` | toxicity risk column(s) appended (e.g. `Mutagenicity`, `Tumorigenicity`, `Irritating Effects`, `Reproductive Effects` — selection-dependent) |
| M4 | `Chem > Calculate > MPO Score...` | MPO score column appended (profile-dependent; optionally desirability columns when enabled) |
| M5 | `Chem > Calculate > BitBIRCH Clustering...` | cluster-id column appended (BitBIRCH WASM, threshold/fingerprint configurable) |
| M6 | `Chem > Calculate > Cluster MCS...` | Cluster-MCS column appended (most-common-substructure per cluster) |
| M7 | `Chem > Calculate > Map Identifiers...` | identifier-mapping column(s) appended (PubChem/ChEMBL — source/target configurable) |
| M8 | `Chem > Calculate > To InchI...` | InChI column appended |
| M9 | `Chem > Calculate > To InchI Keys...` | InChI Key column appended |
| M10 | `Chem > Calculate > Biochemical Properties` | one or more biochemical-property columns appended (dialog discovers functions tagged `function_family: biochem-calculator`) |

For each `(D, M)` cell:

1. Open the dataset `D` (close any previously open views first to start from a clean state). Wait
   for the table view to render with the molecule column populated and the RDKit cell renderer
   applied.
2. From the top menu, run the Calculate menu item `M` (path per the table above). The
   corresponding dialog opens.
3. In the opened dialog, accept arbitrary values: select the auto-detected molecule column where
   prompted (Variant D — two molecule columns — pick the first molecule column; defaults are
   acceptable for the remaining inputs). Where the dialog presents per-property / per-descriptor
   checkboxes (M1, M2, M3, M10), select at least one item to ensure a non-empty calculation
   request; full defaults are acceptable when the dialog pre-selects items.
4. Click **OK** to run the calculator.
5. Verify: the expected column (per the table above) is appended to the active table view's grid.
   The grid row count is unchanged. No console errors fire during dialog open, OK click, or
   column-append.
6. Close the active view before moving to the next `(D, M)` cell — calculators append columns to
   the active table, so a fresh open keeps each cell isolated.

Implicit cross-cell invariants:

- Each Calculate menu item is invoked independently per dataset — no inter-cell state coupling.
- Across the 4 datasets × 10 menu items, dialog open + OK click never crashes the console; failed
  individual cells (e.g. M1 against a notation that the calculator does not support) surface as
  dialog-level error balloons, not as silent failures (assert per-cell that EITHER a column is
  appended OR a user-visible error balloon is shown — never silent).

## Notes

- **M10 scope.** Biochemical Properties (M10) is included in the walk even though it's a later
  addition to the Calculate menu, because it's a registered menu item and the original directive
  was to test every Calculate submenu entry.
- **Cross-cutting bug awareness — GROK-17621 (Descriptors Docker timeout).** This scenario walks
  the happy path only and assumes the Docker container is available; the container-unavailable /
  timeout failure mode has no dedicated spec yet — it is an open coverage gap, not exercised here.
- **Cross-cutting bug awareness — github-2942 (CSV export with filtering).** A row-count-mismatch
  bug when exporting a filtered table as CSV / Molecules-As-SMILES intersects this scenario via
  the notation-conversion dependency used by the Calculate menu walks, but has no dedicated
  bug-focused spec yet — this scenario neither filters nor exports.
- **MPO Score profile dependency (M4).** The MPO Score dialog requires at least one MPO profile to
  be available. The default profiles shipped with the Chem package should suffice on a fresh test
  instance; if none is available, the dialog shows a clear user-visible message rather than
  failing silently.
- **Map Identifiers dependency (M7).** Map Identifiers requires the Chembl package's
  identifier-mapping query surface; on environments without Chembl connectivity, the dialog shows
  a clear user-visible error.
- **Biochemical Properties dependency (M10).** The Biochemical Properties dialog discovers
  functions tagged `function_family: biochem-calculator`, which other packages may register. On
  environments with none registered, the dialog may show an empty / informational state — assert
  per-cell that either columns are appended or a clear message is shown, never silent.
- **No sibling spec.** There is no automated spec covering this Calculate menu walk yet.
