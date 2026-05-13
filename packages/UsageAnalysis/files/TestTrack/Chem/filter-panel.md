---
feature: chem
sub_features_covered: [chem.search.substructure, chem.search.substructure.filter, chem.search.substructure.editor, chem.search.use-as-filter, chem.sketcher, chem.sketcher.ocl]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel.md
migration_date: 2026-05-11
migration_report: filter-panel-migration-report.md
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario (`depends_on: []`),
`classification: medium`, `pyramid_layer: integration`, `target_layer: playwright`, strategy `simple`.
UI coverage owned (`ui_coverage_delegated_to: null`) over the 11 entries in
`ui_coverage_responsibility` — `chem-substructure-filter-draw`,
`chem-substructure-filter-settings-contains` / `-included-in` / `-exact` / `-similar`,
`chem-substructure-use-as-filter`, `chem-substructure-hamburger-menu-add-filter`,
`chem-substructure-filter-remove`, `chem-substructure-filter-drag-drop-to-panel`,
`chem-substructure-filter-cross-filter-other-types`,
`chem-substructure-filter-sketcher-modify-from-hamburger`. Bug `GROK-14028` (Reset filter doesn't
clear sketcher input) intersects this scenario at the Block B / Block E sketcher-state-cleanup
surfaces and is cited in `related_bugs`; defense-in-depth coverage of the three-layer cleanup
invariant (UI clear + BitSet null + visual cleanup) is owned by the dedicated bug-focused spec
`chem-grok-14028-spec.ts` per chain `bug_focused_candidates[]`.

## Setup

1. **Provision dataset.** The linked dataset is bundled in the platform's System file shares (no
   external provisioning required):
   - `System:AppData/Chem/tests/spgi-100.csv` — SMILES-notation molecule column + numeric activity / properties.
   Order `11` per source JSON footer.
2. **Confirm Chem package is loaded** so the substructure filter widget is registered as the primary
   filter for `semType: Molecule` columns (atlas `chem.search.substructure.filter`, package source
   `Chem/src/package.ts`). The OCL sketcher backend (atlas `chem.sketcher.ocl`) must be available as
   the platform fallback so the in-Filter-Panel sketcher renders without an external sketcher package.
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
   tile (atlas `chem.search.substructure.filter`).
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
6. On the grid, right-click a molecule cell and select **Current Value > Use as filter** (atlas
   `chem.search.use-as-filter`). Verify the Filter Panel reopens with the substructure filter
   populated with the clicked molecule's structure (orientation preserved when the source is SMILES,
   per atlas), and the table filters to rows containing that substructure.

#### Block C — Filter from the column header hamburger menu

7. Close the Filter Panel.
8. Open the molecule column's header hamburger menu and select **Filter**. Verify a sketcher dialog
   opens (atlas `chem.search.substructure.editor` — the editor function prompts for a query when
   invoked from the column hamburger menu route).
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

- **`coverage_type: regression`** — multi-entry-point × settings-tabs × cross-filter × filter↔
  sketcher-sync walk across the substructure filter surface. Not `smoke` (the section's smoke is
  `Advanced/scaffold-tree-functions.md` per chain `ui_coverage_plan.smoke_scenario`); not `edge` /
  `perf` (no specific failure-mode invariant or threshold being asserted — Reset-filter-doesn't-clear-
  sketcher-input invariant per `GROK-14028` is owned by the dedicated bug-focused spec
  `chem-grok-14028-spec.ts` per chain `bug_focused_candidates[]`, NOT this happy-path walk; the
  `MAX_SUBSTRUCTURE_SEARCH_ROW_COUNT` row-count bound per atlas `chem.search.substructure` is not
  exercised — spgi-100.csv sits well below threshold). `regression` is the natural fit for a multi-entry-
  point × multi-interaction walk that guards against regression-of-the-set across the substructure
  filter surface.
- **No JS API substitution.** Every entry in chain `ui_coverage_responsibility` is exercised via UI
  driving — Filter Panel toggle, in-panel sketcher draw, settings-tabs cycle, right-click context
  menu **Use as filter**, column hamburger-menu **Filter**, sketcher dialog draw + **Add filter**,
  filter-widget remove control, drag-and-drop of the column header, additional non-molecule filter
  composition, and the hamburger-menu re-edit propagation to the filter widget. Substituting any of
  these via `grok.dapi.*` or programmatic BitSet manipulation would defeat the `pyramid_layer:
  integration` purpose and lose the filter↔sketcher synchronization invariant verification.
- **Existing sibling spec.** A test file already exists at
  `public/packages/UsageAnalysis/files/TestTrack/Chem/filter-panel-spec.ts` (per
  `existing-test-index.yaml` line 32480: `test_name: "Chem: Filter Panel"`, `category: Chem`,
  `layer: playwright`, `test.skip` pattern present, helpers `spec-login`, `features_covered:
  [chem.filter-panel]`). The Automator will extend it to cover the full block-A-through-E walk per
  this migrated scenario.
- **Bug `GROK-14028` awareness.** The Reset-filter cleanup gap (Reset clears BitSet but sketcher UI
  input retains the structure) intersects Blocks B (Use as filter populates the sketcher; subsequent
  Reset must clear it) and E (hamburger-menu Re-edit modifies the filter; Reset must clear all three
  layers). This scenario does NOT exercise the Reset action explicitly — defense-in-depth coverage
  of the three-layer cleanup invariant (UI clear + BitSet null + visual cleanup) is owned by the
  dedicated bug-focused spec `chem-grok-14028-spec.ts` per chain `bug_focused_candidates[]`.
- **Helpers usage.** Standard `loginToDatagrok` + `softStep` + `closeAllViews` from
  `helpers-registry.yaml` cover the section harness. No registered helper currently abstracts the
  Filter Panel toggle / in-panel-sketcher draw / settings-tabs cycle / drag-and-drop column header
  pattern; candidate helpers surfaced in migration report Decisions.
- **Order in chain.** `order: 11` per source JSON footer.
- **Source-text fixes silently applied during migration.** The original body had severe step-
  numbering chaos (numerals `1.1.1.2.1.2.1.1.2.2.2.3.4.` interleaved without semantic grouping — the
  markdown renderer auto-renumbers display but the source is unintelligible). Renumbered with
  explicit block labels (Block A — Filter Panel basics, Block B — Use as filter, Block C — Hamburger
  menu Filter, Block D — Drag-and-drop column header, Block E — Cross-filter + sketcher↔filter sync)
  per chain rev 2 footer note (c) (Olena 2026-05-11). Typo `sketcherr` → `sketcher`. `Drag'n'drop` →
  `Drag-and-drop` (canonical phrasing). See migration report Decisions § Source-text fixes.
