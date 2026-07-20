---
feature: chem
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
priority: p1
realizes_atlas: [github-3004]
realizes: []
realized_as:
  - chem-github-3004-spec.ts
produced_from: atlas-driven
produced_for: chem-github-3004-spec.ts
authored_date: 2026-05-11
authored_by: claude-code-test-designer
related_bugs: [github-3004]
---

# Chem | Scaffold Tree multi-table binding (github-3004 regression-lock)

Bug-focused scenario locking the invariant: **with ≥2 table views open, Top-menu
`Chem | Analyze | Scaffold Tree` MUST bind the new Scaffold Tree viewer to the
ACTIVE table view (`grok.shell.tv.dataFrame`), NOT to the first-opened table**.
This is Datagrok's "active table view" convention — Scaffold Tree must follow
the same default-binding behavior as other Chem viewers. Binding the new
viewer to a non-active table (the literal github-3004 surface) is the
regression this scenario catches.

This is a parallel, narrower scenario alongside the single-table Scaffold Tree walk
(`Advanced/scaffold-tree.md`) and the create+configure smoke walk
(`Advanced/scaffold-tree-functions.md`) — this one is scoped specifically to the
multi-table binding behavior.

Bug reference: github-3004 — *Two tables with molecules: scaffold tree is based on
wrong table by default in some cases* (priority p2, fixed in 1.21.0).

## Setup

1. **No external provisioning.** The scenario constructs the multi-table state
   inline — opens `System:AppData/Chem/tests/smiles-50.csv` twice with
   distinct table names (50 rows each, `canonical_smiles` auto-detects
   `semType: Molecule`). Two table views, each with one molecule column —
   sufficient to surface the binding invariant.
2. **Login** — standard test user. No special privileges required.
3. **No fixture consumed.** Independent scenario (`depends_on: []`).

## Scenarios

### Block A — Scaffold Tree from active table binds to active table, not first-opened

Locks the invariant: when the user has TableView A and TableView B open with
TableView B as the active view (`grok.shell.tv === TableView B`), invoking
`Chem | Analyze | Scaffold Tree` from the top menu MUST add the Scaffold Tree
viewer to TableView B's dataframe (B's dataframe is `viewer.dataFrame`), NOT
TableView A's.

1. **Open TableView A** — read `System:AppData/Chem/tests/smiles-50.csv` into
   `dfA`, rename `dfA.name = 'tableA'`, add the table view. Wait for
   semantic-type detection (`canonical_smiles` becomes `semType: Molecule`).
2. **Open TableView B** — read `System:AppData/Chem/tests/smiles-50.csv` into
   a NEW DataFrame `dfB`, rename `dfB.name = 'tableB'`, add the table view.
   Wait for semantic detection. The active table view at this point becomes
   TableView B (per Datagrok's default — newly-added view becomes active).
3. **Assert active table is TableView B.** `grok.shell.t.name === 'tableB'`
   AND `grok.shell.tv.dataFrame.name === 'tableB'`. This validates the
   pre-condition before invoking Scaffold Tree — if active table is NOT
   tableB at this step, the test setup itself failed (not the bug).
4. **Invoke Scaffold Tree from active TableView B via top menu.** Dispatch a
   click on `[name="div-Chem"]` (Chem top-menu) → wait for sub-menu render →
   find `.d4-menu-item-label` with text `"Scaffold Tree"` → click its parent
   `.d4-menu-item` element. Per chem.md note: Scaffold Tree menu item has NO
   `...` suffix and opens directly (no parameter dialog).
5. **Wait for Scaffold Tree viewer to materialize.** Up to 10 seconds. Look
   for `[name="viewer-Scaffold-Tree"]` in the DOM AND for a
   `"Scaffold Tree"`-type viewer in `Array.from(grok.shell.tv.viewers)` on
   TableView B.
6. **Assert regression-lock invariant — binding to active table.** Locate the
   Scaffold Tree viewer instance via
   `Array.from(grok.shell.tv.viewers).find(v => v.type === 'Scaffold Tree')`.
   Verify `viewer.dataFrame.name === 'tableB'`. If `viewer.dataFrame.name`
   equals `'tableA'` instead — the github-3004 regression has landed (viewer
   bound to first-opened table, not active).
7. **Verify the wrong-table-view is unaffected.** Verify that the OTHER
   table view (TableView A) does NOT have a Scaffold Tree viewer attached
   (no spurious binding to the inactive table either). Scope by switching
   to TableView A: `grok.shell.v = <TableView A>` → check
   `Array.from(grok.shell.tv.viewers).filter(v => v.type === 'Scaffold Tree').length === 0`.

## Notes

- **Why the top-menu vector.** The bug lives in the top-menu-driven default-table-binding
  logic, so the scenario drives Scaffold Tree via the top menu rather than the equivalent
  JS API call (`Chem:addScaffoldTree`) — the JS API path shares the same underlying
  function but isn't separately asserted here.
