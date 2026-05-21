# Features Table Layout — Mockup Comparison

Feature: Post-Stage-4 ("Features") view of the Consensus Pharmacophore app.
The question motivating this: "why does the table show 25 rows when I only gave it 5 PDBs?"

## Option A — Per-PDB summary in grid; interaction detail in context panel

The grid is replaced with a 5-row aggregation: one row per input PDB. Count columns
(`n_hbd`, `n_hba`, `n_hyd`, `n_pi`, `n_xb`, `n_total`) are colored by family code. Clicking
a row opens the full 9-row (or n-row) per-interaction detail in the Datagrok context panel on
the right. The 3D viewer always shows all 25 chain-F spheres regardless of which row is selected.

Trade-offs: immediately answers the "5 inputs = 5 rows" question; colored count columns let the
user compare family coverage across ligands at a glance. The cost is that the raw `ligand_features`
DataFrame is no longer the live grid object — it must be cached and rendered only in the context
panel. Power users who want to sort/filter/export individual interaction rows lose direct grid
access. Implementation requires a JS aggregation pass and a custom context-panel widget. The
context panel is the standard Datagrok idiom for entity detail (see `add-info-panel` skill), so
the pattern is familiar.

## Option B — Keep 25-row detail grid; add per-PDB summary above buttons in left rail

The grid is unchanged (`ligand_features`, 25 rows). A new "SUMMARY" accordion block is inserted
in the left input rail above the preview buttons, showing one compact entry per PDB: `pdb_id ·
ligand_comp_id · N interactions` plus a row of color-coded family chips (`D ×1`, `A ×1`, `H ×3`,
etc.). A warning icon flags PDBs with `skip_reason` rows. Clicking a PDB entry in the summary
applies a Datagrok filter on `pdb_id` and shows a filter pill in the grid toolbar; clicking the
pill's × clears it.

Trade-offs: zero change to the data model — the 25-row DataFrame is preserved exactly, so Stage
5a still receives it unchanged and sorting/filtering/exporting remains full-featured. The summary
answers the "why 25?" question by surfacing the per-PDB breakdown in the UI without replacing the
grid. Implementation cost is low: a JS loop over the existing DataFrame produces the HTML. The
rail becomes taller and may require scrolling on small screens to reach BUILD. The filter-on-click
interaction is additive convenience, not required for the feature to work.

## Option C — Current state (baseline, no changes)

Grid = 25-row `ligand_features` as ProLIF returns it. Left rail = PDB textarea + options + 5
preview buttons + BUILD. Mol* right pane = chain-F spacefill overlay.

Trade-offs: no implementation work required. Every interaction is directly visible and sortable.
The 25-row count is correct and expected by Stage 5a downstream. However, without any summary,
a new user confronted with the grid naturally asks "why 25 rows?" — the connection between 5
PDBs and 25 interactions is not surfaced anywhere in the UI. Repeated `pdb_id` and `ligand_comp_id`
values across consecutive rows make the table visually noisy.

---

## Recommendation: Option B

Option B is recommended. It preserves the `ligand_features` DataFrame intact (required for Stage
5a k-means and for any future export), adds actionable context directly in the left rail where
the user is already looking at their inputs, and answers the "why 25 rows?" question the moment
the stage completes. The filter-on-click interaction is a natural Datagrok pattern (grid + filter
pill) and can be implemented as a pure JS loop + `df.filter` call with minimal surface area.
Option A is more visually elegant but requires a non-trivial context-panel widget, a caching
strategy for the raw DataFrame, and means the grid object no longer holds the stage output
directly — which breaks the current `view.dataFrame = newDf` swap pattern the orchestrator
relies on. Option C preserves the status quo but leaves the core UX confusion unaddressed.

### Open questions for reviewer / user

1. Should the summary block in Option B be collapsible (accordion)? The block is ~130 px tall for
   5 PDBs; on a 800-px-tall window the BUILD button may fall off-screen on first load.
2. In Option B, should clicking a PDB chip in the summary highlight the corresponding spheres in
   Mol* (e.g. dim all non-matching PDB's spheres)? This would require Mol* state management.
3. In Option A, if the user re-runs the pipeline with a new PDB list, should the context panel
   stay open with stale data, or be auto-closed? The Datagrok context panel pattern auto-updates
   on selection, but the aggregated DataFrame is synthetic.

### Precedents consulted

- `.claude/mockups/consensus-pharmacophore-3d/variant-a-state-result.html` — left-rail layout,
  stage badge bar, grid tab bar, Mol* pane proportions, property-grid spacing.
- `.claude/mockups/consensus-pharmacophore-3d/variant-a-app.html` — input panel accordion
  structure, demo badge, preview button style.
- `.claude/skills/ui/SKILL.md` — Datagrok UI idioms: context panel pattern, `ui.accordion()`,
  `DG.Viewer.grid()`, `ui.splitH`.
- `.claude/skills/add-info-panel/SKILL.md` (name only) — confirms context panel is the platform
  idiom for entity detail.
- `src/orchestrator.ts` (lines 67–218) — actual input panel construction, stage badge HTML,
  `view.dataFrame = newDf` swap strategy.
