---
feature: chem
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [chem.cp.substructure-search-with-filter]
realizes: [chem.substructure-filter]
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/structure-filter.md
migration_date: 2026-05-11
source_text_fixes: []
candidate_helpers:
  - helpers.playwright.chem.openFilterPanel
  - helpers.playwright.chem.drawStructureFilter
  - helpers.playwright.chem.useAsFilterFromCell
  - helpers.playwright.chem.hamburgerMenuFilter
  - helpers.playwright.views.cloneCurrentView
unresolved_ambiguities:
  - set-a-filter-blocks-2-3-4-preamble-entry-point-freedom
  - filter-ui-is-correct-block-2-step-3-verification-depth
  - some-structure-cell-selection-block-5-step-4-which-row
  - block-4-4-4-numbering-interpretation
scope_reductions:
  - id: SR-01
    check: A-STRUCT-02
    rationale: "A-STRUCT-02 carryforward (chain-level edge/perf coverage)"
    verdict_status: SCOPE_REDUCTION
related_bugs: []
---

# Structure Filter

End-to-end exercise of the Chem Substructure Filter on `spgi-100.csv` across
five invariants. Each `### <name>` sub-scenario is an independent
slice of the substructure-filter state machine:

- **Disable / re-enable through panel close+reopen** — drawing a
  structure filter survives toggling the filter off, closing the panel,
  re-opening, and toggling back on (filter restored intact).
- **Cell context-menu Use as filter after panel close** — running
  *Current value > Use as filter* with the filter panel closed builds
  a correct filter UI.
- **Hamburger-menu Filter + draw + Add filter** — combining an existing
  filter with a new sketcher-drawn structure via the column hamburger
  menu's Filter action and the *Add filter* button.
- **Remove + Use as filter ordering** — after removing a substructure
  filter, *Current value > Use as filter* sketches the molecule
  correctly AND positions the new Structure filter first on the
  Filter Panel.
- **Clone view + cross-view filter sync** — a substructure filter
  applied through *Current value > Use as filter* on one view is
  synchronized across both views of the same table.

Exercises the substructure search + filter workflow end-to-end.

## Setup

1. Open `spgi-100.csv` (`System:AppData/Chem/tests/spgi-100.csv`) — the
   linked dataset for this scenario per the original JSON footer
   (`{ "order": 4, "datasets": ["System:AppData/Chem/tests/spgi-100.csv"] }`).
   Confirm the molecule column auto-detects as semantic type
   `Molecule`.

## Scenarios

### Disable filter + close+reopen panel + re-enable

Filter draw → disable → panel close → panel reopen → re-enable;
verify each transition preserves the filter shape end-to-end.

1. Open the **Filters** panel and add a Substructure filter for the
   molecule column — either draw a structure in the filter sketcher
   OR right-click a molecule cell and choose **Current value > Use
   as filter**. The substructure filter is added to the Filter Panel
   and the table is filtered to matching rows.
2. Disable the Structure filter (toggle the filter off in the Filter
   Panel). The filter is inactive; the table un-filters; the sketcher
   UI on the filter widget retains the drawn structure.
3. Close the Filter Panel. The panel collapses; filter state is
   preserved.
4. Open the Filter Panel. The Structure filter is still present, still
   disabled, sketcher UI still shows the previously-drawn structure;
   no console errors.
5. Turn the Structure filter back on. The filter re-applies; the table
   re-filters to matching rows identical to the pre-disable state; no
   console errors.

### Set filter then Use as filter after panel close

Verify the cell context-menu *Current value > Use as filter* action
constructs a correct filter UI even when invoked with the filter
panel closed.

1. Set a substructure filter (draw a structure in the filter sketcher
   or via **Current value > Use as filter**). The table is filtered.
2. Close the Filter Panel.
3. Right-click a molecule cell and choose **Current value > Use as
   filter**. The Filter Panel re-opens (or the filter applies behind
   the closed panel); the Substructure filter UI shows the molecule
   from the clicked cell sketched correctly; the table is re-filtered
   accordingly; no console errors.

### Hamburger menu Filter + Draw + Add filter

Verify combining an existing filter with a new sketcher-drawn
structure via the column hamburger menu and the *Add filter*
button.

1. Set a substructure filter (draw a structure or via
   **Current value > Use as filter**). The table is filtered.
2. Close the Filter Panel.
3. On the molecule column header, open the hamburger menu and choose
   **Filter**. A sketcher dialog (or in-place filter sketcher widget)
   opens for the molecule column.
4. Draw another structure in the sketcher.
5. Click **Add filter**. The new structure is combined with the
   existing filter (added as an additional Substructure filter on
   the Filter Panel); the table re-filters to rows matching the
   combined predicate; no console errors.

### Remove + Use as filter ordering

Verify that *Current value > Use as filter* after the Structure filter
has been removed re-creates the filter correctly AND positions it
first on the Filter Panel.

1. Remove the Structure filter from the Filter Panel (delete / close
   the filter entry).
2. Right-click a molecule cell and choose **Current value > Use as
   filter**. Verify:
   - the molecule is sketched correctly in the new Structure filter
     sketcher (orientation preserved when the source column is
     SMILES),
   - the Structure filter is the first entry on the Filter Panel
     (positioned ahead of any other filters present).
3. Confirm no console errors fire during the remove → re-add
   transition.

### Clone view + cross-view filter sync

Verify a substructure filter applied through *Current value > Use as
filter* on one view is synchronized across both views of the same
table.

1. Open the linked dataset (`spgi-100.csv`) and open the Filter Panel.
   Confirm at least one structure column is added to the filter
   panel.
2. Add another view onto the same table (e.g. clone the current view
   via the view tab context menu / *Clone view* action). The cloned
   view opens as a second tab/view on the same DataFrame.
3. In the cloned view, open the Filter Panel and confirm the same
   structure column is added to the filter panel in the 2nd view.
4. In one of the two views, right-click a molecule cell and choose
   **Current value > Use as filter**. The substructure filter applies
   and the filtered row set is synchronized across both views — both
   views show the same filtered subset.

## Notes

- **Cross-cutting bug awareness — GROK-14028.** The Filter Panel Reset cleanup gap (Reset clears the
  filter's BitSet but the sketcher UI input retains the drawn structure) is owned by a dedicated
  bug-focused spec, not this scenario. Block 1 (disable + close panel + re-open) is a defense-in-depth
  sister check on the same UI-state-survival question — does panel re-rendering preserve filter state? —
  but is distinct from the Reset-clears-BitSet-not-UI invariant that GROK-14028 targets.
- **Verification scope.** "Sketched correctly" (Block 4) means the sketcher renders the molecule from
  the clicked cell without falling back to an empty/placeholder state — not necessarily an exact SMILES
  round-trip comparison.
- **Ordering invariant.** Block 4's assertion that the newly re-added Structure filter appears first on
  the Filter Panel is a presentation-order rule specific to this scenario; it is not documented as a
  general Filter Panel layout rule elsewhere.
- **Why cross-view sync should hold.** Filter state is a property of the underlying DataFrame rather
  than of an individual view, so both views should always show the same filtered set; Block 5 makes this
  behavior observable explicitly.
