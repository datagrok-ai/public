---
feature: chem
sub_features_covered: [chem.search.substructure, chem.search.substructure.filter, chem.search.substructure.top-menu, chem.search.use-as-filter, chem.sketcher]
target_layer: playwright
coverage_type: regression
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Chem/Advanced/structure-filter.md
migration_date: 2026-05-11
migration_report: structure-filter-migration-report.md
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

Per chain YAML (`scenario-chains/chem.yaml` rev 2): independent scenario
(`depends_on: []`), `classification: medium`, `pyramid_layer: integration`,
`target_layer: playwright`, strategy `simple`. Realizes atlas critical
path `chem.cp.substructure-search-with-filter` (p0).

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
     SMILES per atlas `chem.search.use-as-filter`),
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

- **Dataset.** spgi-100.csv is the canonical chem test dataset for the
  substructure-filter cycle. Provisioning path: `System:AppData/Chem/tests/spgi-100.csv`.
- **Filter-panel toggle.** Datagrok exposes the Filter Panel via the
  `.grok-icon-filter` toggle on the Table View toolbar (per MEMORY.md
  Datagrok UI Selector Notes — `table_view.dart:647`). Container is
  `.d4-filter-panel`. Spec selector pinning is the Automator's
  responsibility.
- **Use as filter cell action.** The `Current value > Use as filter`
  cell action is provided by atlas `chem.search.use-as-filter`
  (`Chem/src/package.ts#L1957`) — preserves orientation when the
  source column is SMILES.
- **Cross-cutting bug awareness — GROK-14028 (Filter Panel Reset
  cleanup gap).** GROK-14028 affects `chem.search.substructure.filter`
  (this scenario's coverage). The bug invariant — Reset filter clears
  the BitSet but the sketcher UI input retains the SMILES — is owned
  by the **bug-focused candidate spec** `chem-grok-14028-spec.ts`
  (chain `bug_focused_candidates[]` entry; spans
  `filter-panel.md:Step 2` + `Advanced/structure-filter.md:Step 1`
  + `sketcher.md:Step 6`). Per `pattern 2` (bug-focused slice — no
  separate `.md` authored at this stage), `related_bugs` is intentionally
  empty here and the bug invariant lives in the dedicated bug-focused
  spec. Block 1 (disable + close panel + re-open) IS a defense-in-depth
  sister check against the same UI-state-survival invariant — does
  panel re-rendering preserve filter state? — but distinct from the
  Reset-clears-BitSet-not-UI invariant of GROK-14028 itself.
- **No JS API substitution for UI flows.** All 8 entries of
  `ui_coverage_responsibility`
  (`chem-structure-filter-draw`, `chem-substructure-use-as-filter`,
  `chem-substructure-disable-filter`, `chem-substructure-remove-filter`,
  `chem-substructure-hamburger-menu-add-filter`,
  `filter-panel-close-open`, `clone-view`, `filter-cross-view-sync`)
  MUST be exercised through Playwright DOM driving — Filter Panel
  toggle clicks, sketcher draw, cell context-menu navigation,
  hamburger-menu navigation, view-clone action, and cross-view
  inspection are UI behaviors the scenario contracts on. The atlas
  also exposes `grok.chem.searchSubstructure(col, query)` for the
  underlying compute, but `pyramid_layer: integration` plus
  `ui_coverage_delegated_to: null` mandate UI driving for this
  scenario.
- **"Use as filter" sketcher rendering correctness.** Block 4
  asserts the molecule is sketched correctly when re-added via
  *Current value > Use as filter* after a prior Remove. The
  Automator's verification surface is bounded to "the sketcher
  widget renders the molecule from the clicked cell without
  fallback / placeholder / empty state" — exact SMILES round-trip
  comparison is spec-time discretion.
- **First-on-panel ordering invariant.** Block 4 asserts the
  Structure filter is the first entry on the Filter Panel after the
  *Use as filter* re-add. This is a presentation-order invariant
  on the Filter Panel layout (newly-added Structure filter on top,
  ahead of pre-existing non-Structure filters). The atlas does not
  enumerate the exact ordering rule; the original scenario's
  assertion governs.
- **Cross-view sync invariant.** Block 5 verifies that the filtered
  row set on the underlying DataFrame is synchronized across two
  views — both views witness the same `filter.selectedIndexes`.
  Datagrok's filter state is a property of the underlying DataFrame
  rather than per-view, so the sync should be implicit; the
  scenario explicitly contracts on observable cross-view
  consistency.
