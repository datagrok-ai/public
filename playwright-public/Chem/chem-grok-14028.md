---
feature: chem
sub_features_covered: [chem.search.substructure.filter, chem.search.substructure, chem.sketcher, chem.actions.copy-smiles]
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
produced_from: atlas-driven
produced_for: chem-grok-14028-spec.ts
authored_date: 2026-05-11
authored_by: claude-code-test-designer
related_bugs: [GROK-14028]
---

# Chem | Filter Panel Reset 3-layer cleanup invariant (GROK-14028 regression-lock)

Bug-focused scenario locking the **3-layer state cleanup invariant** for the
substructure-filter Reset action: any state transition that "clears" the filter
MUST clear ALL THREE layers — (1) sketcher UI input value, (2) underlying
BitSet (filter.trueCount === rowCount), (3) any RDKit-renderer visual
highlighting / column tags.

Per chain YAML (`scenario-chains/chem.yaml` rev 2 `bug_focused_candidates[chem-grok-14028-spec.ts]`):
independent scenario (`depends_on: []`), `pyramid_layer: bug-focused`,
`target_layer: playwright`, strategy `simple`. Parallel-coverage on
`chem.search.substructure.*` with the broader walk in `filter-panel.md`
(`coverage_type: regression`, multi-entry-point walk including sketcher draw +
match-mode tabs + Use-as-filter + hamburger). This spec owns the cross-cutting
defense-in-depth Reset cleanup invariant.

Bug-library reference: `references/bug-library/chem.yaml :: GROK-14028` —
title *Chem | Filter Panel: Reset filter doesn't clear the input line in sketcher*,
priority p2, status `regression-risk`, `test_coverage: needed`. The
chain analyzer flagged this as cross-cutting (3 affecting scenarios: filter-panel.md
Step 2 + structure-filter.md Step 1 + sketcher.md Step 6) — none of the 3
exercise the 3-layer cleanup invariant.

## Setup

1. **Provision dataset.** `System:AppData/Chem/tests/spgi-100.csv` —
   molecule column on which substructure filter can be applied with
   deterministic non-trivial filtered count.
2. **Login** — standard test user.
3. **No fixture consumed** (independent scenario).

## Scenarios

### Block A — Apply substructure filter, then Reset, assert 3-layer cleanup

Locks the invariant: after Reset, sketcher UI is empty AND BitSet is all-true
AND no substructure-highlight visual cache remains.

1. **Open spgi-100.csv** and wait for the table view + RDKit cell renderer to
   render the molecule column.
2. **Open the Filter Panel** via `grok.shell.tv.getFiltersGroup()`. Verify the
   substructure filter widget is present (filter with Structure / molecule
   column header).
3. **Apply a substructure filter** by sketching / pasting a deterministic SMILES
   (`c1ccccc1` — benzene) into the substructure-filter sketcher dialog and
   clicking OK. Verify the filter is applied — `grok.shell.t.filter.trueCount`
   strictly less than `grok.shell.t.rowCount` and strictly greater than 0.
4. **Trigger Reset.** Invoke the Reset filter action — UI vector via the
   substructure-filter hamburger / context menu Reset entry (Datagrok standard
   `Reset filter` action). The Reset action MUST be the same surface the user
   would click — JS API substitution for the trigger itself is forbidden by
   the bug-focused slice (the bug is in the Reset handler, not BitSet).
5. **Assert the 3-layer cleanup invariant.** After Reset, ALL THREE MUST hold:
   - **(L1) BitSet cleared.** `grok.shell.t.filter.trueCount === grok.shell.t.rowCount`
     (all rows visible — no filter applied).
   - **(L2) Sketcher UI cleared.** The substructure-filter widget's sketcher
     input / canvas reads as empty — either the inline sketcher value query
     returns empty / null SMILES, OR the .d4-filter for the molecule column
     no longer carries a substructure-highlight tag.
   - **(L3) No leaked column tags / visual cache.** No
     `~<col>.substructure-highlight` or `chem-substructure-filter` tag-column
     remains attached to the dataframe under a stale state.

   If any of (L1)/(L2)/(L3) FAILS — Reset cleanup is incomplete and the bug is
   present.

6. **Verify no console errors.** No uncaught Error during the Reset window
   (defensive guard).

## Notes

- **`coverage_type: edge`** — 3-layer state cleanup edge case for
  `chem.search.substructure.filter`. Happy-path filter usage is owned by parallel
  `filter-panel.md` (`coverage_type: regression`, `pyramid_layer: integration`).
- **3-layer invariant.** The bug-library entry is explicit: "any state
  transition that 'clears' filter must verify all three layers — sketcher UI
  clear, BitSet null, no row highlighting". Block A asserts all three after
  the same Reset action that the bug-library repro names.
- **No JS API substitution for the Reset trigger.** The bug LIVES in the Reset
  action handler. Triggering Reset via direct BitSet manipulation
  (`df.filter.setAll(true)`) would bypass the bug and produce a false-PASS.
  JS API IS permitted for setup (open table, open filter panel, programmatic
  SMILES input) and for assertions (read `.trueCount`, `.tags`, etc.) — but
  the Reset itself MUST be the UI-level Reset filter action.
- **L2 selector strategy.** Substructure-filter widget exposes its current
  SMILES via `.chem-canvas-wrapper` text content / data attribute, OR via
  `grok.shell.tv.dataFrame.col(<molcol>).getTag('chem-substructure-filter')`
  (atlas-documented tag). Spec checks both.
- **L3 column tag list.** Default spgi-100.csv columns. Post-Reset spec walks
  `grok.shell.t.columns.names()` to confirm no `~*substructure*` or
  `~*highlight*` virtual columns linger.
- **MCP recon NOT performed pre-author** for this scenario — author relies on
  filter-panel-spec.ts patterns + chain rationale. Section 7 hypothesis-
  investigation handles any selector / Reset-action surface that doesn't match
  assumption at Validator Gate B.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` from
  `helpers-registry.yaml`. No new helper authored.
- **Cross-cutting scope.** Per chain spans, this spec lives at the intersection
  of `filter-panel.md`, `structure-filter.md`, and `sketcher.md`. Block A
  exercises the Reset path that defines the bug invariant; sketcher copy-SMILES
  + Filter Panel insert from the bug repro Step 1-2 are setup, NOT the test
  slice (the cross-cutting invariant is the Reset handler).
- **Order in chain.** Bug-focused candidate; independent scenario
  (no `depends_on:`).
