---
feature: chem
sub_features_covered: [chem.analyze.scaffold-tree, chem.analyze.scaffold-tree.viewer, chem.analyze.scaffold-tree.filter, chem.search.substructure, chem.search.substructure.api, chem.sketcher.cell-editor]
target_layer: playwright
pyramid_layer: bug-focused
coverage_type: edge
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2
`bug_focused_candidates[chem-grok-12758-spec.ts]`): independent scenario
(`depends_on: []`), `pyramid_layer: bug-focused`, `target_layer: playwright`,
strategy `simple`. Cross-cutting on three Chem subsystems —
`chem.analyze.scaffold-tree` (viewer + filter), `chem.search.substructure`
(API + UI), `chem.sketcher.cell-editor` (scaffold edit dialog). Parallel-
coverage with `Advanced/scaffold-tree.md` (single-vector tree-population
smoke) and `Advanced/scaffold-tree-functions.md` (JS API walk).

Bug-library reference: `references/bug-library/chem.yaml :: GROK-12758` —
title *Chem | Scaffold Tree: incorrect molecule rendering after opening some
scaffold tree elements*, priority p1, status `regression-risk`,
`test_coverage: needed`.

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

- **`coverage_type: edge`** — cross-cutting state-isolation edge case across
  Scaffold Tree + Sketcher + Substructure Search. Happy-path scaffold tree
  smoke is owned by `Advanced/scaffold-tree.md`.
- **Mixed JS API + UI vector.** Setup (open table, addTableView, generate
  tree, count nodes) uses JS API for reliability. The bug-triggering actions
  (open Edit Scaffold dialog, CANCEL, click checkbox) MUST be UI clicks —
  the bug is in the UI handler chain (sketcher state save/restore +
  substructure-filter handler share internal state); direct JS API
  invocations would bypass the bug.
- **30s post-Generate budget.** Empirical on dev.datagrok.ai 2026-05-12:
  scaffold tree on spgi-100.csv populates in a few seconds, but reaching the
  stable settled state (toolbar `empty-tree` dropped, all visible nodes
  rendered, canvases painted) is generously budgeted at ≥30s.
- **Selector strategy — visible-node filter.** `.d4-tree-view-node` count
  includes a hidden empty root wrapper. Use the canvas presence as a
  "visible node" filter:
  `Array.from(nodes).filter(n => n.querySelector('canvas.chem-canvas'))`.
- **4th visible node choice** — bug-library repro specifies the 4th node
  as the trigger. Some scaffold trees may have <4 visible nodes initially
  (e.g. small heterogeneous datasets); if so, the spec degrades gracefully
  by using `Math.min(3, visibleNodes.length - 1)`. spgi-100.csv typically
  yields enough visible nodes to satisfy the 4-node threshold; if not, the
  fallback kicks in.
- **Browser hang during recon (2026-05-12).** MCP recon of this surface
  triggered a `Runtime.evaluate` timeout / Accessibility tree hang after
  the Edit-Scaffold dialog dismiss + checkbox click sequence. This may
  itself be a GROK-12758 symptom (sketcher resource leak corrupts page
  state). Spec accommodates via softStep boundary: if checkbox click hangs,
  Validator surfaces a per-step timeout error rather than a clean
  invariant failure. Cross-reference: chem.md "GROK-12758 known fragility".
- **regression-risk caveat.** Bug-library lists GROK-12758 as
  `regression-risk` with `fixed_in: ''`. If spec FAILS on first Validator
  run, mark `test.fixme()` with bug-state-dependent rationale (mirror
  github-2942 pattern). If spec PASSES, bug may have been silently fixed —
  update bug-library accordingly.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` from
  `helpers-registry.yaml`. No new helper authored.
- **Order in chain.** Bug-focused candidate; independent scenario
  (no `depends_on:`).
