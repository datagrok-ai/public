---
feature: powerpack
target_layer: playwright
coverage_type: smoke
priority: p1
realizes: []
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
scope_reductions:
  - id: SR-01
    check: E-SCENARIO-RUNTIME-ALIGNMENT
    rationale: |
      The function panel re-sorts on column change, but the scenario-cited
      example function families do not reorder deterministically on the chosen
      columns. The assertion is made at the order-changed level (the top-5
      differs); the "matching family on top" expectation is log-only, not a
      hard assertion.
    verdict_status: SCOPE_REDUCTION
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

The dialog's basic open-and-add flow is covered by `add-new-column.md`.

## Setup

A clean Datagrok session is the only shared setup. The scenario opens
its own dataset; no fixture chaining is required.

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

- **Edge cases.** The sort-by-type behavior generalizes to any
  column type; string columns (e.g. `Chemist`) follow the same
  pattern as the Molecule and numeric cases exercised in Steps 3-4.

---
{
  "order": 5
}
