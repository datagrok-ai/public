---
feature: chem
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
priority: p1
realizes_atlas: [GROK-12758]
realizes: []
produced_from: atlas-driven
produced_for: chem-grok-12758-spec.ts
authored_date: 2026-05-12
authored_by: claude-code-test-designer
related_bugs: [GROK-12758]
---

# Chem | Scaffold Tree node opened in Sketcher + checkbox filter (GROK-12758 regression-lock)

Bug-focused scenario locking the **Scaffold-Tree × Sketcher × Substructure
state-isolation invariant**: opening a scaffold tree node in the Sketcher
(via the per-node `Edit scaffold` pencil icon) MUST NOT corrupt internal
substructure-search state. Subsequent filter checkbox click on the SAME node
MUST apply the substructure filter cleanly — no console errors from
`searchSubstructure` ("Search pattern cannot be set"), no
"crossed-out" molecule rendering.

This scenario cuts across three Chem subsystems — the Scaffold Tree viewer and
filter, substructure search, and the sketcher cell editor. It runs parallel to
`Advanced/scaffold-tree.md` (single-vector tree-population walk) and
`Advanced/scaffold-tree-functions.md` (create+configure smoke walk).

Bug reference: GROK-12758 — *Chem | Scaffold Tree: incorrect molecule rendering
after opening some scaffold tree elements* (priority p1).

## Setup

1. **Provision dataset.** `System:AppData/Chem/tests/spgi-100.csv` — has
   Structure column (Molecule semType, 100 rows) — sufficient for non-trivial
   scaffold tree with several visible nodes post-Generate.
2. **Login** — standard test user.
3. **No fixture consumed** (independent scenario).

## Scenarios

### Block A — Open scaffold node in sketcher, then apply as filter, assert clean state

Locks the invariant: Edit-Scaffold dialog dismiss + per-node filter checkbox
click on the SAME node = clean substructure-filter application with no
console errors and no crossed-out rendering.

1. **Open spgi-100.csv** and wait 30s for Chem autostart cascade (Structure
   detected as Molecule, RDKit renderer + Scaffold Tree function registered).
2. **Open Scaffold Tree viewer** via top-menu `Chem | Analyze | Scaffold Tree`.
   Verify `[name="viewer-Scaffold-Tree"]` is present.
3. **Click magic wand** (`[aria-label="Generate"]`) to auto-generate scaffold
   tree from Structure column. Wait ≥30s for tree population (empirically
   a few seconds on dev for spgi-100.csv; toolbar drops `empty-tree` class
   when ready).
4. **Verify tree populated** — `.d4-tree-view-node` count > 4 (need at least
   4 visible scaffold nodes to exercise "4th node" repro step).
5. **Hook console.error** to capture `Search pattern cannot be set` /
   `searchSubstructure` error patterns before any node interaction.
6. **Click `Edit scaffold` pencil icon on the 4th visible node** —
   `i[aria-label="Edit scaffold"]` inside the 4th `.d4-tree-view-node`
   carrying a `canvas.chem-canvas`. The `Edit Scaffold...` dialog opens
   (`.d4-dialog` containing `.grok-sketcher`).
7. **Click CANCEL** to dismiss without edits (the bug fires regardless of
   whether the sketcher contents were modified — the open-then-close action
   alone corrupts internal state per bug-library).
8. **Wait 2s for state settle**, then **click the filter checkbox** on the
   SAME 4th `.d4-tree-view-node` (`input[type="checkbox"]` direct child of
   the node). This is the "Select the checkbox next to the 4th structure to
   apply it as a filter" step in the bug repro.
9. **Wait 4s for filter application**.
10. **Assert clean filter state.** ALL THREE MUST hold:
    - **(A1) No `searchSubstructure` error class** — console.error hook
      captures zero entries matching
      `/Search pattern cannot be set|searchSubstructure/i`.
    - **(A2) Filter applied cleanly** — `df.filter.trueCount` strictly
      between 0 and `df.rowCount` (some rows masked out by the scaffold
      filter, but not all-zero / all-one).
    - **(A3) No selection corruption** — `df.selection.trueCount === 0`
      (bug surface: "molecules crossed out instead of clean filter" implies
      selection / strike-through visual state pollution; check selection
      BitSet is untouched).

## Notes

- **Mixed JS API + UI vector.** Setup (open table, add table view, generate tree, count
  nodes) uses the JS API for reliability. The bug-triggering actions — open Edit Scaffold
  dialog, CANCEL, click checkbox — must be real UI clicks, because the bug is in the UI
  handler chain (sketcher state save/restore and the substructure-filter handler share
  internal state); a direct JS API call would bypass the bug.
- **Why the 4th node.** The bug-library repro specifically names the 4th visible scaffold
  node as the trigger. If a scaffold tree has fewer than 4 visible nodes, the spec falls
  back to the last available node instead.
- **Bug status caveat.** GROK-12758 is tracked as `regression-risk` rather than
  definitively fixed — if this scenario starts failing, treat it as a live bug
  reproduction rather than a flaky test.
