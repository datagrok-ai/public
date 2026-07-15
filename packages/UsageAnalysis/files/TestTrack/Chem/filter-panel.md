---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [chem.cp.substructure-search-with-filter]
realizes: [chem.substructure-filter]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel.md
migration_date: 2026-05-11
source_text_fixes:
  - step-numbering-chaos-resolved-with-block-labels
  - typo-sketcherr-sketcher
  - drag-n-drop-drag-and-drop-canonical-phrasing
candidate_helpers:
  - helpers.playwright.chem.toggleFilterPanel
  - helpers.playwright.chem.drawSubstructureInFilterPanelSketcher
  - helpers.playwright.chem.cycleSubstructureFilterSettingsTabs
  - helpers.playwright.chem.useAsFilter
  - helpers.playwright.chem.openColumnHamburgerFilter
  - helpers.playwright.chem.dragColumnHeaderToFilterPanel
unresolved_ambiguities:
  - sketch-a-substructure-which-specific-substructure
  - check-all-filter-settings-exact-verification-granularity-per-tab
  - use-as-filter-molecule-cell-selection
  - hamburger-menu-filter-exact-menu-label
  - add-filter-button-label-in-sketcher-dialog
  - check-the-filtering-of-chemical-columns-with-any-other-filter-which-other-filter
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: [GROK-14028]
---

# Chem | Filter Panel deep-dive — substructure filter draw / settings / Use as filter / hamburger / drag-drop / cross-filter / sketcher↔filter sync

Filter Panel walk over the substructure-filter widget on `System:AppData/Chem/tests/spgi-100.csv`. Exercises every
entry point into the substructure filter — sketcher draw inside the Filter Panel, the four
substructure-match settings tabs (Contains / Included in / Exact / Similar), the cell action
**Current Value > Use as filter**, the column hamburger-menu **Filter** entry, **Add filter** from the
sketcher dialog, drag-and-drop of the column header onto the Filter Panel, and the cross-filter
interaction with another (non-molecule) filter type — plus the filter↔sketcher synchronization
invariant when the sketcher is modified through the hamburger-menu route.

Bug GROK-14028 (Reset filter doesn't clear the sketcher input) intersects this scenario at the
Use-as-filter and hamburger-menu-modify surfaces; defense-in-depth coverage of the full three-layer
cleanup invariant (UI clear + BitSet null + visual cleanup) is owned by a dedicated bug-focused spec,
not by this walk.

## Setup

1. **Provision dataset.** The linked dataset is bundled in the platform's System file shares (no
   external provisioning required):
   - `System:AppData/Chem/tests/spgi-100.csv` — SMILES-notation molecule column + numeric activity / properties.
   Order `11` per source JSON footer.
2. **Confirm Chem package is loaded** so the substructure filter widget is registered as the primary
   filter for `semType: Molecule` columns (package source `Chem/src/package.ts`). The OCL sketcher
   backend must be available as the platform fallback so the in-Filter-Panel sketcher renders without
   an external sketcher package.
3. **No fixture consumed.** Per chain YAML `depends_on: []`. spgi-100.csv is opened fresh inside the
   scenario.

## Scenarios

### Filter Panel deep-dive on spgi-100.csv

Single-pass walk on a fresh spgi-100.csv table view, organized into five sequential blocks. Each block
opens / inspects / interacts with the substructure filter through a different entry point. Blocks
share the same active table view (no per-block reload); inter-block invariants are stated where they
apply.

Open the dataset:

1. Open `System:AppData/Chem/tests/spgi-100.csv` (close any previously open views first to start from a clean
   state). Wait for the table view to render with the molecule column populated and the RDKit cell
   renderer applied.

#### Block A — Filter Panel basics: draw a substructure query

2. Open the Filter Panel (toggle the filter icon on the table view toolbar). Verify the Filter Panel
   appears as a side pane with the molecule-column substructure filter widget rendered as a sketcher
   tile.
3. In the Filter Panel's substructure filter, sketch a substructure (any deterministic small fragment
   — e.g. a benzene ring or pyridine — drawn via the OCL sketcher inline canvas). Verify the table
   filters down to rows whose molecule column contains the sketched substructure; matching atoms are
   visually highlighted in the molecule cells (substructure highlight overlay rendered by the RDKit
   cell renderer).
4. Open the filter settings (the gear / context menu on the substructure filter widget). Cycle
   through the four match-mode tabs in turn — **Contains**, **Included in**, **Exact**, **Similar** —
   and for each tab verify that the filter re-evaluates and the visible-row set updates accordingly:
   - **Contains** — rows whose molecule contains the sketched substructure (default).
   - **Included in** — rows whose molecule is fully included as a substructure of the sketched
     query (inverse direction).
   - **Exact** — rows whose molecule is an exact structural match to the sketched query.
   - **Similar** — rows whose molecule has Tanimoto similarity above the configured threshold to
     the sketched query.
   Each tab switch must complete without console errors and produce a deterministically-different
   filtered row set (or, where the set legitimately overlaps, the same set with the mode selector
   visibly toggled).

#### Block B — Use as filter from the molecule cell context menu

5. Close the Filter Panel.
6. On the grid, right-click a molecule cell and select **Current Value > Use as filter**. Verify the
   Filter Panel reopens with the substructure filter populated with the clicked molecule's structure
   (orientation preserved when the source is SMILES), and the table filters to rows containing that
   substructure.

#### Block C — Filter from the column header hamburger menu

7. Close the Filter Panel.
8. Open the molecule column's header hamburger menu and select **Filter**. Verify a sketcher dialog
   opens (the editor function prompts for a query when invoked from the column hamburger menu route).
9. Draw another structure in the dialog and click **Add filter**. Verify the substructure filter is
   added to the Filter Panel with the newly-drawn structure, and the table filters accordingly.
10. Remove the structure filter from the Filter Panel (the filter widget's remove control / context
    menu **Remove filter** entry). Verify the substructure filter widget is removed from the Filter
    Panel and the table reverts to no-substructure-filter state.

#### Block D — Drag-and-drop the column header onto the Filter Panel

11. Drag-and-drop the molecule column's header onto the Filter Panel. Verify the substructure filter
    is added back to the **top** of the Filter Panel (insertion position: top of the panel, not
    bottom) and the filter is empty (no pre-populated query) until the user draws one.

#### Block E — Cross-filter with another column type + sketcher↔filter sync via hamburger menu

12. Apply any additional non-molecule column filter (e.g. add a numeric range filter on a numeric
    spgi-100 column, or a categorical filter on a string column). Verify the substructure filter and the
    additional filter compose (logical AND) — the visible-row set is the intersection of both
    filters; toggling either filter independently updates the row set accordingly without console
    errors. Confirms substructure filter participates correctly with other filter types.
13. With the substructure filter active in the Filter Panel, open the molecule column's header
    hamburger menu again and select **Filter** (re-enters the sketcher dialog with the current filter
    structure pre-loaded). Modify the structure in the sketcher (e.g. add or remove an atom / bond)
    and apply the change (click **Add filter** / **OK** / equivalent commit control). Verify the
    Filter Panel's substructure filter updates to match the modified structure AND the table re-
    filters accordingly. **In the filter widget and in the sketcher dialog we must always observe the
    same structure** — the filter↔sketcher synchronization invariant: modifying the sketcher through
    the hamburger-menu route must propagate to the filter widget (visible same structure) and the
    filter re-applies to the data (table re-filters).

## Notes

- **Sibling spec.** A Playwright spec already exists at `filter-panel-spec.ts` ("Chem: Filter
  Panel"); it is extended to cover the full Block-A-through-E walk described here.
