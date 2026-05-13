---
feature: chem
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.viewer, chem.analyze.scaffold-tree.add]
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2
`bug_focused_candidates[chem-github-3004-spec.ts]`): independent scenario
(`depends_on: []`), `pyramid_layer: bug-focused`, `target_layer: playwright`,
strategy `simple`. Parallel-coverage on `chem.analyze.scaffold-tree.*` with
`Advanced/scaffold-tree.md` (single-table Scaffold Tree walk,
`coverage_type: regression`) and `Advanced/scaffold-tree-functions.md`
(JS API function walk). This spec owns the multi-table state-tracking
invariant.

Bug-library reference: `references/bug-library/chem.yaml :: github-3004` —
title *Two tables with molecules: scaffold tree is based on wrong table by
default in some cases*, priority p2, status `fixed`, `fixed_in: '1.21.0'`,
`test_coverage: needed`. NIBR customer signal (NX label).

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

- **`coverage_type: edge`** — multi-table state-tracking edge case; this
  scenario asserts the active-table binding invariant for
  `chem.analyze.scaffold-tree.add`. Happy-path single-table Scaffold Tree
  flows are owned by parallel `Advanced/scaffold-tree.md`
  (`coverage_type: regression`, `pyramid_layer: integration`).
- **Top-menu vector choice.** Top menu `Chem | Analyze | Scaffold Tree` is
  the canonical bug-report path. JS API alternative
  (`grok.functions.call('Chem:addScaffoldTree', {})`) also adds the viewer
  to the active table view — but the bug surface is specifically the
  top-menu-driven default-binding decision in the package's `addScaffoldTree`
  implementation. Test exercises the UI vector to confirm the user-facing
  fix; the JS API path is implicitly covered (both call the same package
  function under the hood) but not separately asserted here.
- **Dataset construction — two reads of same file.** Both tables use
  `smiles-50.csv` (50 rows, `canonical_smiles` Molecule column). Distinct
  table names (`tableA` / `tableB`) ensure the bug's "which DataFrame is
  current" invariant is exercised — the multi-table state error is about
  DataFrame identity, not data content. Same source file simplifies setup
  while still surfacing the binding error.
- **Active table view assumption.** Test relies on Datagrok's default
  behavior where the newly-added table view becomes the active view (per
  step 2 above). If the platform later changes this default, Step 3's
  pre-condition assert catches the setup regression cleanly before the
  Scaffold Tree invocation.
- **MCP-validated current behavior (2026-05-11 dev.datagrok.ai recon, qa-pw
  user, Chem package v1.17.6).** Two table views open, smiles-50-copy as
  active. Top-menu Chem | Analyze | Scaffold Tree → viewer added → bound
  to smiles-50-copy (active table). github-3004 fix verified on current
  build. Spec exercises the post-fix invariant as a regression-lock.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` from
  `helpers-registry.yaml` + `warmChemPackage` (env-bootstrap helper,
  precedent batch-1 + batch-2). Inline patterns for dataset open +
  top-menu dispatch — no new helper authored under reuse threshold
  (single-use pattern).
- **No JS API substitution.** Step 4 exercises the top-menu UI vector
  that the bug lives in. JS API invocation
  (`Chem:addScaffoldTree`) would bypass the top-menu-driven default-table
  resolution that the bug surfaces — UI driving is required to exercise
  the bug's exact code path.
- **Cold-Playwright readiness (resolved 2026-05-12 cycle).** Original
  authoring proactively applied `test.fixme()` under the env-pending
  Chem-package-cold-start hypothesis. Post-breakthrough cycle
  (chem-new-2026-05-12) confirmed that the breakthrough pattern
  (split-evaluate setup + 30s outside-evaluate post-addTableView settle)
  resolves the cold-start blocker. test.fixme REMOVED; Validator Gate B
  PASS verified across multiple 3-run sequences. One Round-1 test-bug
  surfaced (`grok.shell.tableViews` Iterable vs Array); fixed by
  wrapping in `Array.from(...)`.
- **Order in chain.** Bug-focused candidate; placed parallel to
  `Advanced/scaffold-tree.md` and `Advanced/scaffold-tree-functions.md`.
  No `depends_on:` ordering — independent scenario.
