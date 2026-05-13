---
feature: chem
sub_features_covered: [chem.analyze.elemental, chem.analyze.elemental.top-menu, chem.analyze.elemental.run]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis.md
migration_date: 2026-05-11
migration_report: elemental-analysis-migration-report.md
related_bugs: []
---

# Chem | Analyze | Elemental Analysis

Multi-format happy-path walk of the **Chem | Analyze | Elemental Analysis** calculator. For each
linked dataset variant (smiles / molV2000 / molV3000), open the dataset, run the Elemental
Analysis top-menu item, toggle every per-element checkbox on in the dialog, click **OK**, and
verify that the per-element atom-count columns are appended to the active table view's grid.

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: simple`, `pyramid_layer: integration`, `target_layer: playwright`, strategy
`simple`. Multi-format implied — the original body promises "smiles, molV2000 and molV3000
formats" but the JSON footer declared only `smiles.csv`; chain rev 2 directive (footer note (i),
Olena 2026-05-11) instructs Migrator to extend the dataset list during migration. UI coverage
owned (`ui_coverage_delegated_to: null`) over `chem-elemental-analysis-top-menu` and
`chem-elemental-analysis-checkboxes` per chain `ui_coverage_responsibility`.

## Setup

1. **Provision linked datasets.** All three files are bundled in the platform's System file shares
   (no external provisioning required):
   - `System:DemoFiles/chem/smiles.csv` — SMILES-notation molecule column (Variant A — canonical
     SMILES).
   - `System:AppData/Chem/mol1K.sdf` — molV2000 SDF (Variant B — V2000 molblock notation).
   - `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` — molV3000 SDF (Variant C — V3000 molblock
     notation).
2. **Confirm Chem package is loaded** so that the **Chem | Analyze | Elemental Analysis...**
   top-menu item (`elementalAnalysis`, source `public/packages/Chem/src/package.ts#L961`) is
   registered. The Charts and PowerGrid packages should also be loaded so that the optional Radar
   viewer + per-row Radar grid-cell follow-ons are available; their absence surfaces a clear
   user-visible warning (per atlas `chem.analyze.elemental.top-menu` description) but does not
   block the per-element atom-count column append covered by this scenario.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. Each dataset is opened fresh inside
   the scenario.

## Scenarios

### Elemental Analysis per format

Data-driven walk: for each dataset variant (D1–D3), open the dataset, run Elemental Analysis with
all per-element checkboxes toggled on, click **OK**, and verify per-element atom-count columns
appear in the grid.

Dataset matrix:

| id | dataset path | format |
|----|--------------|--------|
| D1 | `System:DemoFiles/chem/smiles.csv` | smiles |
| D2 | `System:AppData/Chem/mol1K.sdf` | molV2000 |
| D3 | `System:DemoFiles/chem/sdf/ApprovedDrugs2015.sdf` | molV3000 |

For each dataset `D`:

1. Open the dataset `D` (close any previously open views first to start from a clean state). Wait
   for the table view to render with the molecule column populated and the RDKit cell renderer
   applied.
2. From the top menu, run **Chem > Analyze > Elemental Analysis...**. The Elemental Analysis
   dialog opens with the molecule column auto-detected and a list of per-element checkboxes.
3. Toggle every per-element checkbox in the dialog **on** (select all elements).
4. Click **OK** to run the calculator with these parameters (default behavior for any remaining
   inputs).
5. Verify: per-element atom-count columns are appended to the active table view's grid (one
   column per element that was checked, with integer atom counts per row). The grid row count is
   unchanged. No console errors fire during dialog open, OK click, or column-append.
6. Close the active view before moving to the next dataset — Elemental Analysis appends columns
   to the active table, so a fresh open keeps each variant isolated.

Implicit cross-variant invariants:

- Elemental Analysis is invoked independently per dataset — no inter-variant state coupling.
- Across the 3 variants, dialog open + OK click never crashes the console; if the calculator
  cannot parse a notation variant the dialog surfaces a user-visible error balloon, not a silent
  failure (assert per-variant that EITHER per-element columns are appended OR a user-visible
  error balloon is shown — never silent).

## Notes

- **`coverage_type: regression`** — multi-format happy-path walk over `chem.analyze.elemental`.
  Not `smoke` (the section's smoke is `Advanced/scaffold-tree-functions.md` per chain
  `ui_coverage_plan.smoke_scenario`; chain rev 1 notes a deviation from strict Rule 1 — strict
  step-count tie-break would pick `elemental-analysis.md` (3 steps) but the analyzer prefers the
  Scaffold Tree viewer smoke because Elemental Analysis is a Chem|Analyze calculator-style flow,
  not a viewer create+configure smoke). Not `edge` / `perf` (no specific failure-mode invariant
  or threshold being asserted — the happy path appends per-element columns; missing-Charts /
  missing-PowerGrid warning behavior is implicit and not exercised here). `regression` matches:
  protects against regression-of-the-set across the three supported notation formats.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility`
  (`chem-elemental-analysis-top-menu`, `chem-elemental-analysis-checkboxes`) is exercised via UI
  driving — top-menu walk + per-element checkbox interaction + OK click. The atlas
  `chem.analyze.elemental.run` (`runElementalAnalysis(table, molecules)`, role `transform`) is
  the JS API surface for the calculation, but invoking it directly would defeat the dialog-driven
  workflow (per-element checkbox toggling) that the scenario asserts. Per `pyramid_layer:
  integration`.
- **Optional Radar viewer / per-row Radar grid cells.** Atlas `chem.analyze.elemental` notes
  Charts (Radar viewer) and PowerGrid (per-row Radar grid cells) are optional add-on surfaces
  enabled when those packages are present. The Migrator preserves this awareness but the
  per-element atom-count column append (step 5) is the canonical verification — the optional
  Radar surfaces are NOT exercised in this scenario (out of scope of the original body, which
  asserts only that the calculator runs and adds columns).
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/elemental-analysis-spec.ts` (per
  `existing-test-index.yaml` line 32462: `test_name: "Chem: Elemental Analysis"`,
  `category: Chem`, `layer: playwright`, helpers `spec-login`, `features_covered:
  [chem.elemental-analysis]`, patterns `uses-fixture-runners` / `uses-grok.dapi` /
  `uses-grok.shell` / `uses-page.evaluate` / `uses-page.locator` / `uses-playwright-test`). The
  existing spec likely covers a single format; the Automator extends it to the full 3-format
  enumerated walk per this migrated scenario.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts the
  Elemental Analysis top-menu walk + per-element checkbox toggle pattern; candidate helpers
  surfaced in migration report Decisions.
- **Order in chain.** `order: 4` per source JSON; tie-broken lexicographically with
  `Advanced/structure-filter.md` (also order 4) — `elemental-analysis.md` runs after.
- **Source-text fixes silently applied during migration.** The original body had a step-numbering
  glitch (`1.` / `2.` / `2.` / `3.` — second `2.` should have been `3.`, then original `3.`
  should have been `4.`). Renumbered cleanly. The JSON dataset list was extended from a single
  `smiles.csv` entry to the three-format set (smiles + molV2000 + molV3000) per chain rev 2
  directive (footer note (i), Olena 2026-05-11) to align with the body's "smiles, molV2000 and
  molV3000 formats" coverage promise. An explicit per-element atom-count column verification
  step was added (original step 4 said "Make sure the result is correct" — this migration
  surfaces "verify per-element atom-count columns appended to the active table view's grid" as
  the concrete verification). See migration report Decisions § Source-text fixes.
