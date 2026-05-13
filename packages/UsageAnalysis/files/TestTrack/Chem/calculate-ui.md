---
feature: chem
sub_features_covered: [chem.calculate.descriptors, chem.calculate.descriptors.top-menu, chem.calculate.descriptors.transform, chem.calculate.descriptors.tree, chem.calculate.properties, chem.calculate.properties.top-menu, chem.calculate.toxicity-risks, chem.calculate.toxicity-risks.top-menu, chem.mpo, chem.mpo.top-menu, chem.calculate.bitbirch, chem.calculate.bitbirch.top-menu, chem.calculate.cluster-mcs, chem.calculate.cluster-mcs.top-menu, chem.calculate.map-identifiers, chem.calculate.map-identifiers.top-menu, chem.notation.inchi, chem.notation.inchi-keys, chem.calculate.biochemical-properties]
target_layer: manual-only
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/calculate.md
migration_date: 2026-05-11
migration_report: calculate-migration-report.md
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
active table view. Realizes atlas critical path `chem.cp.calculate-descriptors-docker` (p1) for
the Descriptors anchor.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: medium`, `pyramid_layer: integration`, `target_layer: playwright`, strategy
`simple`. Multi-format × multi-menu-item — Variant C representative source: smiles.csv;
Descriptors (Docker-backed) is the anchored test case, the remainder of the Calculate menu is the
implicit walk now enumerated below. UI coverage owned (`ui_coverage_delegated_to: null`) over
Descriptors / Properties / Toxicity Risks / MPO Score / BitBIRCH Clustering / Cluster MCS / Map
Identifiers / To InchI / To InchI Keys / Biochemical Properties.

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
   registered (per atlas `chem.calculate.*` ids; package source `Chem/src/package.ts`). For
   Descriptors specifically, the `chem-chem` Docker container must be available — failure mode
   bug-library `GROK-17621` (Docker timeout) is outside scope; the happy path assumes container
   is up.
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

- **`coverage_type: regression`** — multi-format × multi-menu-item happy-path walk. Not `smoke`
  (the section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`); not `edge` / `perf` (no specific failure-mode invariant or
  threshold being asserted). The Docker-unavailable / timeout failure-mode invariant for
  `chem.calculate.descriptors` is owned by the cross-cutting bug-focused spec
  `chem-grok-17621-spec.ts` (chain `bug_focused_candidates[]`); this scenario assumes container
  is up.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-calculate-descriptors`, `chem-calculate-properties`, `chem-calculate-toxicity-risks`,
  `chem-calculate-mpo-score`, `chem-calculate-bitbirch-clustering`,
  `chem-calculate-cluster-mcs`, `chem-calculate-map-identifiers`, `chem-calculate-to-inchi`,
  `chem-calculate-to-inchi-keys`) is exercised via UI driving — top-menu walk + dialog interaction
  + OK click. The Biochemical Properties row (M10) extends the menu walk beyond the chain
  `ui_coverage_responsibility` list because atlas `chem.calculate.biochemical-properties` is a
  registered Calculate menu item (`Chem | Calculate | Biochemical Properties`) and the original
  scenario directive ("Do the same for each section in the Calculate menu") enumerates every
  Calculate-menu entry — see migration report Decisions.
- **Cross-cutting bug awareness — GROK-17621 (Descriptors Docker timeout).** Per chain
  `bug_focused_candidates[]`, the Docker container unavailability / timeout invariant for
  `chem.calculate.descriptors` is owned by the dedicated `chem-grok-17621-spec.ts` candidate spec
  (spans `calculate.md:Step 1` per chain). This scenario walks the happy path only — container
  assumed available; the failure-mode invariant is parallel-coverage, not in scope here.
  `related_bugs: []` because the bug-focused candidate owns the invariant.
- **Cross-cutting bug awareness — github-2942 (CSV export with filtering).** Per chain
  `bug_focused_candidates[]`, the filter-active + `Export As CSV` + `Molecules As Smiles`
  row-count-mismatch invariant intersects `calculate.md` via the `convert-column` /
  `convertMoleculeNotation` dependency invoked by Calculate menu walks (spans
  `calculate.md:Step 1` + `filter-panel.md:Step 1` per chain). The dedicated
  `chem-github-2942-spec.ts` candidate spec owns the filter+export combination; this scenario
  does not filter or export. Awareness only.
- **MPO Score profile dependency (M4).** The MPO Score dialog requires at least one MPO profile
  to be available (atlas `chem.mpo.profile-app` / `chem.mpo.profile-tree-browser`). On a fresh
  test instance the default profiles shipped with the Chem package should suffice; if no profile
  is available the dialog will surface a clear user-visible message rather than silent failure
  (per the cross-cell invariant above).
- **Map Identifiers dependency (M7).** Map Identifiers (`getMapIdentifiers`) requires the Chembl
  package's `textToSmiles` / identifier-mapping query surface; on environments without Chembl
  connectivity the dialog will surface a clear user-visible error.
- **Biochemical Properties dependency (M10).** The Biochemical Properties dialog
  (`biochemPropsWidget`, atlas `chem.calculate.biochemical-properties`) discovers functions
  tagged `function_family: biochem-calculator`. On environments without registered biochem
  calculators (other packages may register them), the dialog may surface an empty / informational
  state — assert per-cell that either columns are appended OR a clear user-visible message is
  shown, never silent.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/calculate-spec.ts` (per
  `existing-test-index.yaml` line 32409: `test_name: "Chem: Calculate Descriptors"`,
  `category: Chem`, `layer: playwright`, helpers `spec-login`, `features_covered:
  [chem.calculate]`). The existing spec covers Descriptors only; the Automator will extend it to
  cover the full enumerated Calculate menu walk per this migrated scenario.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts the
  Calculate menu walk pattern; candidate helpers surfaced in migration report Decisions.
- **Order in chain.** `order: 3` per source JSON; tie-broken lexicographically with
  `r-group-analysis.md` (also order 3) — `calculate.md` runs after.
- **Source-text fixes silently applied during migration.** The original body contained two typos
  / defects ("asme" → "same"; "Calculate menu ()" → "Calculate menu" — drop empty parens) and a
  TODO note for `smiles_2_columns.csv` that has been resolved (the file exists at
  `System:AppData/Chem/tests/smiles_2_columns.csv` and is provisioned in Setup as Variant D).
  Step numbering renumbered cleanly (original collided on "1." between the dataset-open block and
  the Descriptors-walk block). See migration report Decisions § Source-text fixes.
