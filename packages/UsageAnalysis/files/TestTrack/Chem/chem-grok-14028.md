---
feature: chem
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
priority: p1
realizes_atlas: [GROK-14028]
realizes: []
realized_as:
  - chem-grok-14028-spec.ts
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

This spec is narrower than, and runs parallel to, the broader `filter-panel.md`
walk (sketcher draw, match-mode tabs, Use-as-filter, hamburger menu) — it owns
specifically the Reset-cleanup invariant.

Bug reference: GROK-14028 — *Chem | Filter Panel: Reset filter doesn't clear the
input line in sketcher* (priority p2). Three other scenarios touch this filter
surface (`filter-panel.md`, `structure-filter.md`, `sketcher.md`) but none of them
exercise the full 3-layer cleanup invariant — this spec is the one that does.

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

- **Why the Reset trigger must be a UI click.** The bug lives in the Reset action's
  handler. Triggering Reset via direct BitSet manipulation (`df.filter.setAll(true)`)
  would bypass the bug and produce a false pass. JS API calls are fine for setup
  (open table, open filter panel, programmatic SMILES input) and for assertions
  (reading `.trueCount`, tags, etc.) — but the Reset action itself must be the
  UI-level control.
- **Cross-cutting scope.** This spec sits at the intersection of `filter-panel.md`,
  `structure-filter.md`, and `sketcher.md`. Block A exercises only the Reset path
  that defines the bug invariant — the sketcher copy-SMILES and Filter Panel
  insertion steps are setup, not the tested slice.
