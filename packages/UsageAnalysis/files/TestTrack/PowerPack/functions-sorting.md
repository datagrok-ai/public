---
feature: powerpack
sub_features_covered:
  - powerpack.dialogs.add-new-column
  - powerpack.dialogs.add-new-column-func
target_layer: playwright
coverage_type: smoke
pyramid_layer: ui-smoke
ui_coverage_responsibility:
  - add-new-column-functions-panel
  - add-new-column-functions-sort-by-name
ui_coverage_delegated_to: add-new-column.md
ui_companion: functions-sorting-ui.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/PowerPack/AddNewColumn/functions-sorting.md
migration_date: 2026-05-20
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities: []
scope_reductions: []
related_bugs: []
realized_as:
  - functions-sorting-spec.ts
gate_verdicts:
  d:
    verdict: PASS
    cycle_id: cycle-2026-05-20-functions-sorting
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
  e:
    verdict: SCOPE_REDUCTION
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T11:05:00Z
    failure_keys: []
  a:
    verdict: PASS
    cycle_id: cycle-2026-05-20-functions-sorting
    timestamp: 2026-05-20T00:00:00Z
    failure_keys: []
    review_round: 1
  b:
    verdict: PASS
    cycle_id: 2026-05-28-powerpack-automate-02
    timestamp: 2026-05-28T10:48:30Z
    spec_runs:
      - spec: functions-sorting-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 30
        failure_keys: []
---

# Add New Column dialog — functions-panel sorting

Focused UI scenario verifying the functions-panel sorting behavior in
the Add New Column dialog: (a) selecting a column re-sorts the functions
list so functions whose input parameter type matches the selected
column appear at the top (Molecule-input functions for chem columns,
numeric-input functions for numeric columns); (b) toggling the sort
icon to "By name" re-sorts the list alphabetically; (c) once "By name"
is active, clicking different columns does NOT change the function order.

`pyramid_layer: ui-smoke` per `scenario-chains/powerpack.yaml` rev 1
(Rule 1 fallback — single short UI flow over one source `SPGI`,
`simple` classification, no related bugs). NOT the chain's elected
smoke witness; `coverage_type: regression` is used here (the chain
reserves `coverage_type: smoke` for the elected smoke scenario at
top-level `add-new-column.md`, which exercises the dialog UI
holistically).

`ui_coverage_delegated_to: add-new-column.md` — this scenario owns
the functions-panel sort-by-type / sort-by-name specialty flows on
its responsibility list and delegates the basic dialog-open-and-add
flow to the smoke witness.

`related_bugs: []` — no curated bug in
`bug-library/powerpack.yaml` intersects this scenario's
`sub_features_covered` for the functions-panel-sorting surface.
(GROK-17109 and GROK-17004 affect the same `powerpack.dialogs.add-new-column`
sub-feature but their reproduction surfaces — formula-recalc-on-rename
across save-with-datasync, and complex-paste-handler crash — are not
exercised by the functions-list sort behaviour. See chain rev 1
`bug_focused_candidates[]` for the canonical cross-cutting spec
proposals.)

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining (`depends_on: []` per chain rev 1).

## Scenarios

### Functions panel re-sorts by column-input type, then by name, then sticks

1. Open the SPGI dataset (`System:DemoFiles/chem/SPGI.csv`) — for
   example via **Browse** > **Files**, **File** > **Open**, or the
   equivalent JS-API loader. **Verify:** the SPGI grid renders with
   its chem `Structure` column (semType `Molecule`) plus numeric
   and string columns (e.g. `Chemical Space X`, `Chemist`).
2. Open the **Add New Column** dialog — click the `Add New Column`
   icon on the table view toolbar (or **Edit** > **Add New Column**
   from the top menu). **Verify:** the Add New Column dialog opens
   with its formula editor (CodeMirror), columns list, functions
   panel, and preview grid visible. The functions panel starts in
   its default sort mode ("By relevance").
3. Click the `Structure` column in the dialog's columns list.
   **Verify:** the functions list on the right re-sorts so that
   functions with `Molecule`-type input parameters appear at the
   top (for example, `canonicalize(molecule)`,
   `convertMolNotation(molecule, ...)`, `convertMoleculeNotation(molecule, ...)`,
   `getCLogP(smiles)`, `getDescriptors(molecules, ...)` precede
   functions whose first input is numeric or string).
4. Click a numeric column (e.g. `Chemical Space X`, type `double`)
   in the columns list. **Verify:** the functions list re-sorts so
   functions whose first input parameter is numeric appear at the
   top (for example, `Abs(x)`, `Acos(x)`, `Asin(x)`, `Atan(x)`,
   `Atan2(a, b)` precede chem-input or string-input families).
   Repeat for any other column type as available and **verify:**
   the matching-parameter family is on top in each case.
5. Click the sort icon on the top right corner of the functions
   panel (`[name="icon-sort-alt"]` — the two blue arrows). **Verify:**
   a popup menu appears with options "By name" and "By relevance".
   Select "By name". **Verify:** the functions list re-sorts
   alphabetically (for example, `Abs`, `Acos`, `Add`, `And`, `Asin`,
   `Atan`, `Atan2`, `Avg`, `BinByDateTime`, `BinBySpecificLimits`,
   `Boolean`, `Call`, ... appear in alphabetical order from the top).
6. With "By name" sort still active, click several different columns
   in the columns list in succession (for example `Structure`, then
   `Chemical Space X`, then `Chemist`). **Verify:** the alphabetical
   function order does NOT change — clicking columns no longer
   re-orders the functions list once "By name" sort is selected
   (sticky-sort contract).

## Notes

- **target_layer: playwright** — chosen because a sibling
  `functions-sorting-spec.ts` already exists at the playwright layer
  (per `existing-test-index.yaml` line 33452, layer `playwright`).
  All other Add-New-Column siblings (`add-new-column-spec.ts`,
  `autocomplete-spec.ts`, `hints-spec.ts`, `highlight-spec.ts`,
  `input-functions-spec.ts`, `formula-refreshing-spec.ts`) are also
  Playwright specs, so the house style is consistent.
- **coverage_type: regression** — per the chain's smoke-scenario
  election (`ui_coverage_plan.smoke_scenario: add-new-column.md`),
  the top-level `add-new-column.md` owns the `smoke` slot. This
  narrower per-flow scenario lands in `regression`. (Per the
  migration prompt's `coverage_type` heuristic: `smoke` for one
  fast happy-path test per feature; `regression` is the default
  for the rest.)
- **Step decomposition note.** The original scenario's Step 3 packs
  two sub-cases into one numbered step: (a) Structure column →
  Molecule-input functions on top; (b) change to a non-Structure
  column and verify the matching-input family is on top. The
  migrated body splits these into Step 3 (Structure / chem) and
  Step 4 (numeric column). Both sub-cases of the original Step 3
  are preserved verbatim — no silent drop. The chain's
  `ui_coverage_responsibility:` for this scenario already
  enumerates both `add-new-column-functions-sort-by-type` and
  `add-new-column-functions-sort-by-name`, so the split aligns
  with the chain's coverage intent.
- **Helpers (already in registry, available for downstream
  Automator):** `softStep`, `loginToDatagrok`, `specTestOptions`,
  `stepErrors` from
  `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts` —
  used by the existing `functions-sorting-spec.ts` per the sibling
  test index entry's `helpers_called: [spec-login]`.
- **Bug-library status:** consulted —
  `bug-library/powerpack.yaml` exists and was scanned; no curated
  bug intersects this scenario's `sub_features_covered` for the
  functions-panel sort surface. Cross-cutting candidates GROK-17109
  and GROK-17004 are emitted at the chain level
  (`bug_focused_candidates[]`) but their `spans` reference
  `AddNewColumn/add-new-column.md`, `AddNewColumn/formula-refreshing.md`,
  and `AddNewColumn/highlight.md` — NOT this scenario.
- **Decision log status:** queried — no
  `failed_attempts WHERE feature == powerpack` entries that touch
  the `AddNewColumn/functions-sorting.md` migration. No retry-skip
  approaches apply.
- **Atlas linkage (`derived_from:` provenance):**
  - `powerpack.dialogs.add-new-column-func` interactions "click
    Add New Column toolbar icon → opens AddNewColumnDialog" and
    "Edit | Add New Column top-menu → opens AddNewColumnDialog"
    are code-derived from
    `public/packages/PowerPack/src/package.ts#L405`.
  - `powerpack.dialogs.add-new-column` is the dialog class with
    the columns list + functions panel + preview grid; source
    anchor `public/packages/PowerPack/src/dialogs/add-new-column.ts#L98`.
    The functions-panel sort-by-type behaviour (matching first
    input parameter type to the selected column) and the sort-icon
    menu ("By name" / "By relevance") are surfaces of this class
    not currently captured as atlas `interactions[]` entries —
    surfaced as a curator candidate for the next atlas regen, but
    not blocking for this migration.
- **Edge cases.** The original scenario implies an edge case in
  Step 5 ("Try to click on columns — function order should not
  change"): the alphabetical-sort sticky-contract. Preserved as
  Step 6 of the migrated body. Step 3's "Change column to some
  other and check that function are sorted by matching parameters"
  generalises across column types; the migrated Step 4 picks a
  representative numeric column (`Chemical Space X`) and notes the
  pattern applies to any column type — matches the
  `functions-sorting-run.md` evidence, which records `string`
  (`Chemist`) as covered by the same assertion.

---
{
  "order": 5
}
