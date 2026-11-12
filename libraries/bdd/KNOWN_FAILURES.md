# BDD known-failure audit

Rechecked on 2026-09-13 on macOS, against the local stand at `http://localhost:8889`, starting
from public master `9bd016909a` and core master `c6f235f620`. Twelve scenarios were tagged;
one tag was stale and eleven still describe product defects. Each retained tag was checked
against its actual failing step and the implementation; an expected failure elsewhere in the
scenario is not evidence for the named bug.

## Fixed test and shared bindings

- [Line-chart selection](../../packages/UsageAnalysis/bdd/features/viewers/line-chart/line-chart-selection.feature):
  choosing a checkbox menu item deliberately leaves the menu open. Area selection ignores a
  drag while that menu is open. Closing it makes the same drag select the expected 82 rows.
  Removed the tag and also closed the menu earlier in the journey, where it made an empty-drag
  assertion pass without exercising selection.
- Tooltip text assertions now count visible matching tooltips. A hidden tooltip can retain its
  old text, and a missing tooltip must satisfy a negative text check. A Chromium unit test covers
  hidden, visible and removed tooltip elements, including assertions that must fail.
- An area hover now settles the viewer before returning. Grid cell tooltip requests use
  `grid.debounced` in `core/client/d4/lib/src/viewers/grid/features/grid_tooltip.dart`, so the
  existing `isRenderPending` signal covers their delay. Without that signal, the corrected
  absence assertion could pass before the broken correlation tooltip appeared.
- Range-check failures now report the sampled value through the assertion. The previous error
  interpolated the initial `NaN` before polling, hiding a valid but out-of-range row height.
- The full regression exposed a cleanup defect in Spaces. Its stage lookup used an ID filter
  which returned no matches, falsely declaring existing fixtures deleted. Subsequent runs failed
  to create spaces because those names already existed. Spaces cleanup now verifies captured IDs
  against the complete listing, as its name lookup already did; unrelated spaces are retained.
  A unit test checks both pre-run and post-run cleanup with a server that returns empty filtered
  listings. The underlying Spaces endpoint does not decode the smart filter as Projects does.
  Setup also uses the open Browse panel's Refresh action: API deletion does not invalidate its
  cached nodes, which otherwise coexist with a newly created fixture of the same name.
- [Space file operations](../../packages/UsageAnalysis/bdd/features/spaces/spaces-entity-ops.feature)
  now use Rule backgrounds to require the fixture space to be current. Previously, a failed copy
  setup left the Demo gallery open and subsequent scenarios could rename the original file.
  The affected TSLA fixture was verified against Git and restored unchanged. A view assertion
  also guards the return from the copy space before deleting the original test copy.
- [Box-plot statistics](../../packages/UsageAnalysis/bdd/features/viewers/box-plot/box-plot-statistics.feature)
  uses the existing `user presses t in box plot viewer` binding for both toggles and explicitly
  closes the preceding checkbox menu. In one full run, the table's delayed initial grid focus
  landed between two unscoped key presses, sending the second to the grid. The trace placed the
  presses 25 ms apart at the one-second autofocus boundary (`xamgle/lib/src/views/table_view.dart`).
  Targeting the viewer tests its shortcut independently of that global focus change; the core
  autofocus behavior itself is unchanged.
- [Pivot table links](../../packages/UsageAnalysis/bdd/features/viewers/pivot-table/pivot-table-links.feature)
  (2026-09-21): a filter `DIS_POP in [AS]` appeared after Control-clicks on the aggregated grid's
  row headers, by a group no click named. Cause in the core: the grid's key handler
  (`d4/lib/src/viewers/grid/features/grid_keyboard_navigation.dart`) read `currentPos` on every
  key, and that getter makes row 0 current when no cell is — so the Control key itself, pressed
  in a focused grid, made row 0 current, and a pivot with Filtering Enabled filtered its source by
  the first group. The handler now reads `currentGridCell?.pos`. The grid status
  (`grid_status.dart`) had the same getter behind its `current row` reading, which reset the
  current cell on every read and hid the defect; it now reports the table's current row and
  column, which is what the grid highlights. The scenario is untagged.

## Product defects fixed on 2026-09-21

Every `@known-failure` tag of that date was retired the same day: each defect below is fixed at
its cause and the scenario that carried the tag passes untagged. Core paths are relative to
`core/`; public source links are relative to this file. Expected values were not changed.

| Scenario | Was | Fix |
| --- | --- | --- |
| [Correlation tooltip disabled](../../packages/UsageAnalysis/bdd/features/viewers/correlation-plot/correlation-plot-cells.feature) | A tooltip with `Pearson R` with `showTooltip=false`. | `client/d4/lib/src/viewers/correlation_plot/correlation_plot_core.dart`: the cell tooltip handler reads the flag and suppresses the default tooltip too. |
| [Correlation color scale](../../packages/UsageAnalysis/bdd/features/viewers/correlation-plot/correlation-plot.feature) | 0.0648 and 0.4124 painted the same full red. | `client/d4/lib/src/common/color_coding.dart`: the numerical branch honours a grid column's own `minScale`/`maxScale` (the plot's -1..1). |
| [Forms after header sort](../../packages/UsageAnalysis/bdd/features/viewers/forms/forms-core.feature), [Forms with Show Current Row off](../../packages/UsageAnalysis/bdd/features/viewers/forms/forms-interactions.feature), [PowerGrid twin](../../packages/PowerGrid/bdd/features/viewers/forms/forms-interactions.feature) | 0 cards: the leading card for row -1 was zero pixels tall and the virtual view laid out nothing. | [Forms](../utils/src/viewers/forms-viewer.ts): the card for no row is a placeholder with a real card's size and no handlers. |
| [Heatmap Colors off — GROK-20619](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map-colors.feature) | The AGE band changed by 0 pixels. | `client/d4/lib/src/viewers/grid/grid_core.dart`: the dense heat-map path shares the cell renderer's `heatmapColors` gate. |
| [Returning to heat-map mode](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map-navigation.feature) | Row height stayed at 28 px. | `client/d4/lib/src/viewers/grid/grid_look.dart`: the mode switch rebuilds the vertical range. |
| [Heat-map Row Height disabled](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map.feature) | The Row Height row was not disabled. | `client/libs/property_grid/lib/property_grid_lib.dart`: the initial dependency pass runs when the controlling property exists on the source, not only when it has an editor. |
| [Heat-map column cap before settings](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map.feature) | 11 columns visible after setting the cap to 3. | `grid_core.dart`: the grid binds `look.viewer` on every look it is given, so a setter's `refreshGrid()` works before the property panel opens. |
| [Second scatter plot's filter after a layout — GROK-20896](../../packages/UsageAnalysis/bdd/features/viewers/scatter-plot/scatter-plot-legend.feature) | The second viewer showed every row. | `client/d4/lib/src/viewer_base/data_frame_viewer.dart`: a viewer that finds another one computing the same formula column looks again shortly instead of waiting for nothing. |
| [Filtered group comparison — GROK-20795](../../packages/EDA/bdd/features/analyze/filtered-group-comparison.feature) | Full-table counts 157/5266/354 instead of 104/2823/279. | [Control comparisons](../../packages/EDA/src/control-comparisons/control-comparisons-ui.ts) clone the columns through the table's filter before factorization. |
| [Empty Pareto objective](../../packages/EDA/bdd/features/pareto-front-objectives.feature) | 17 columns offered, the all-null one included. | `shared/ddt/lib/src/data_frame/column_filter.dart` combines conditions (`numerical; not empty`, `Column.matches` agrees); the [Pareto viewer](../../packages/EDA/src/pareto-optimization/pareto-front-viewer.ts) asks for it. |
| [Similarity leaves no cell blank](../../packages/Bio/bdd/features/calculate/scoring.feature) | Blank in every row whose length differed from the reference's. | `libraries/bio/src/monomer-works/monomer-utils.ts`: similarity scores the reference's positions, as identity does. |
| [Centroid linkage — GROK-19595](../../packages/Dendrogram/bdd/features/clustering/chem-dialog.feature) | No tree: the cluster matrix mixed two labelings. | `libraries/math/.../fastcluster.cpp`: centroid takes the same dendrogram path as median (wasm rebuilt). |
| [Global permissions of a role — GROK-20902](../../packages/UsageAnalysis/bdd/features/users-groups-roles/roles-assignment.feature) | The pane listed every group's grants. | `server/datlas/lib/src/services/privileges_service.dart`: `getPermissions(groupId, global: true)` returns the group's own grants. |
| [Deleting a role keeps its grants — GROK-20904](../../packages/UsageAnalysis/bdd/features/users-groups-roles/roles-assignment.feature) | The delete violated `permissions_user_group_id_fkey`. | `server/db/db_up/20260922_0_permissions_group_cascade.sql` (and `init_db.sql`): the key cascades, so a group's grants are deleted with it; a group delete already drops the permission caches. |

## Tags retired on 2026-09-22

The `bdd/grid-gaps` branch (2026-09-14) carried two tags for defects that master fixed on
2026-09-15, before the branch was merged; both scenarios passed untagged in the review run.

| Scenario | Was | Fix |
| --- | --- | --- |
| [Column tooltip set to Columns — GROK-20890](../../packages/UsageAnalysis/bdd/features/viewers/grid/grid-viewer.feature) | `Invalid argument (index): null` from a pointer on the header of the column. | `client/d4/lib/src/viewers/grid/features/grid_tooltip.dart`: the Columns branch is gated on a data cell. |
| [Tags column — GROK-20888](../../packages/PowerGrid/bdd/features/grid/summary-columns.feature) | The grid body painted over in the chip colour. | [Tags renderer](../../packages/PowerGrid/src/cell-types/tags-cell-renderer.ts): `beginPath` before `roundRect`. The scenario, kept apart for the tag, is now the last one of the summary-columns feature. |

No `@known-failure` tag remained in any project after that date.

## Tags retired on 2026-09-23

Eleven scenarios, checked one by one on 2026-09-23 against master with the core and packages of
that day. Nine were product defects, now fixed. Two tagged claims made no sense: GROK-16111 wanted
a balloon for an empty current row, and GROK-20964 was an API check inside a UI feature.

| Scenario | Was | Resolution |
| --- | --- | --- |
| [A band title that does not fit — GROK-20388](../../packages/PowerPack/bdd/features/annotation-regions/annotation-region-titles.feature), [Auto Layout off](../../packages/PowerPack/bdd/features/annotation-regions/annotation-region-titles.feature), [A resized line chart keeps its band title](../../packages/PowerPack/bdd/features/annotation-regions/annotation-regions-2d-viewers.feature) | A const band's title vanished when its strip was not reserved. | `client/d4/lib/src/common/polygon.dart`: `containsPoint` counted a vertex level with the point on both of its edges; a band's sampled edge has one at every pixel row, so the title's cached anchor read as outside and the in-data fallback drew nothing. `annotation_regions_mixin.dart`: a layout that reserved a strip before the band's polygon existed lays out once more; the line chart, which maps per chart, keeps reserving unmeasured (`mapsDuringLayout`). |
| [Similarity, Diversity and Activity Cliffs on an empty sequence — GROK-16111](../../packages/Bio/bdd/features/analyze/empty-current-row.feature) | A balloon was expected for an empty current row. | Won't fix: a balloon on every click on an empty row would be noise. The feature claims what the viewers report: Similarity Search targets the empty row and lists 11 rows, Diversity Search picks 10, Activity Cliffs finds 2 cliffs; no balloon. The "progress entry never ends" note was the dead shell-reset wait printing without waiting. |
| [The identity function scores a HELM sequence — GROK-20964](../../packages/Bio/bdd/features/calculate/scoring.feature) | `Bio:seqIdentity` throws on HELM. | Removed: a function call is an API test (seqIdentity was never meant for HELM). |
| [The Dendrogram box and the tree — GROK-20640](../../packages/Peptides/bdd/features/sar/from-top-menu.feature) | The box reopened unchecked, and unchecking it left the tree. | [Peptides model](../../packages/Peptides/src/model.ts) finds the tree the Dendrogram package attaches beside the grid by its close button (`dendrogramCloseButton`), for the box and for removing it. |
| [Removing a scaffold under a colored one — GROK-18286](../../packages/Chem/bdd/features/scaffold-tree/colors-and-limits.feature) | The scatter plot colored by "Structure colors" lost its coloring. | [Scaffold Tree](../../packages/Chem/src/widgets/scaffold-tree.ts) `removeColorCoding` removed the colored ancestor's entry for a child that only inherits its color, so the column was dropped and rebuilt, unbinding the plot. The old claim (the last colored scaffold removed, the plot still colored) contradicted the design and became two scenarios. |
| [Names beside an existing canonical_smiles — GROK-20955](../../packages/Chem/bdd/features/transform/names-to-smiles.feature) | "Column named 'canonical_smiles' already exists", nothing written. | [Chem](../../packages/Chem/src/package.ts) `namesToSmiles` takes an unused name: canonical_smiles (2). |

## Tags retired on 2026-10-12

Every tag was pinned to the step it describes (`@known-failure:<line>`, compile.ts / harness.ts): a failure
anywhere else now fails the test. Pinning showed one tag green for the wrong reason.

| Scenario | Was | Resolution |
| --- | --- | --- |
| [Open table residues gives compId — GROK-21118](../../packages/UsageAnalysis/bdd/features/biostructure-viewer/file-open-and-preview.feature) | The residue names came out empty. | [pdb-helper](../../packages/BiostructureViewer/src/utils/pdb-helper.ts): `compId` is a string column (it was made an int column). |
| [A Biostructure viewer takes the structure column by itself — GROK-21119](../../packages/UsageAnalysis/bdd/features/biostructure-viewer/file-open-and-preview.feature) | The viewer never picked a column. | [Mol* viewer](../../packages/BiostructureViewer/src/viewers/molstar-viewer/molstar-viewer.ts): detection returned early on `dataJson`, whose empty value is the string `'null'`. |
| [The context panel title follows the new name — GROK-21137](../../packages/UsageAnalysis/bdd/features/layouts/layout-names-tags.feature) | The title kept the old layout name. | `client/xamgle/lib/src/meta/view_layout_meta.dart`: the title follows a save of the same layout. |
| [Edit tags... keeps the name typed into the context panel — GROK-21138](../../packages/UsageAnalysis/bdd/features/layouts/layout-names-tags.feature) | The card saved its stale copy and put the old name back. | `client/xamgle/lib/src/features/viewer_suggestions.dart`: a save rebuilds the card on the saved entity. |
| [The rows a molecule used as a filter lets through contain that molecule](../../packages/UsageAnalysis/bdd/features/viewers/nx/nx-chain.feature) | Rows without the molecule passed. | [Chem substructure filter](../../packages/Chem/src/widgets/chem-substructure-filter.ts): Remove All deactivated the card before removing it, and the broadcast stopped the peers in the other views; a removed card now releases them. The step's own check matched a Kekulé query (`get_qmol` → `get_mol`, `bindings/nx.ts`). |
| [After a layout saved before a rename, the formula names only existing columns](../../packages/UsageAnalysis/bdd/features/viewers/nx/nx-chain.feature) | A layout wrote its saved `formula` tag over the column's own. | `client/xamgle/lib/src/views/table_view.dart` `applyUserDataTags`; Dart test in "Core: ddt: Calculated columns: Save/Apply". |
| [Converting the same column twice keeps both columns — GROK-20806](../../packages/UsageAnalysis/bdd/features/sequence-translator/oligo-toolkit-translator-structure.feature) | Green for the wrong reason: the scenario failed at its first claim, a column name Convert no longer writes (it writes HELM). | The bug is fixed: with the name corrected, both columns stay, the second named "(2)". |

## Tags retired on 2026-10-18

| Scenario | Was | Resolution |
| --- | --- | --- |
| [Shift-click on a network diagram node adds its rows](../../packages/UsageAnalysis/bdd/features/viewers/network-diagram/network-diagram-selection.feature) | Shift-click replaced the selection (6 rows, not 554). | `client/d4/lib/src/viewers/network_diagram/network_diagram_core.dart`: the click took its modifiers from the last key event the viewer's root received, never sent while the focus is outside it; it reads them off the click's own pointer event (vis's `event.srcEvent`). The scenario's Control-click claim said "only the M rows" (447) for the one row left; it now says no F row is selected. |
| [Save SDF with an antisense strand that cannot be converted is refused by name — GROK-20959](../../packages/UsageAnalysis/bdd/features/sequence-translator/oligo-toolkit-translator-structure.feature) | The refusal balloon does not name the strand. | The product names it ("Unable to detect the format of the sequence 'NOTASEQUENCE'"); the scenario claimed the one balloon twice, and the first claim takes it. One claim now. |

## Fixed instead of tagged on 2026-11-08

PRs 4479, 4488, 4516 and 4541 arrived with 35 known-failure scenarios, 13 of them in five temporary
`<folder>-known-failures.feature` files. None was merged as a known failure. The defects below were fixed at the
cause and are claimed by ordinary steps or tests; the rest are one `## Still open` line each in their folder's
MISSING.md, with the claim to add once fixed. Core paths are relative to `core/`.

| Defect | Fix | Claimed by |
| --- | --- | --- |
| core#1159, public#2548 (master's tag, retired): Up/Down no longer moved a legend's sole selection | `client/d4/lib/src/widgets/legend/legend_interaction.dart` `wireKeys`, wired by `legend_mixin.dart` and `legend.dart` (the handler "Legend: Refactor" dropped) | `viewers/legend/legend-interaction.feature` |
| core#1174: a removed column left its rows filtered (Combined Boolean, Expression and free-text rules — which also threw — the bar chart's click filter, the sunburst's and tree's) | `bool_columns_filter.dart`, `free_text_filter.dart` (a renamed column is followed), `bar_chart_core.dart` `onLookChanged`; Charts `sunburst-viewer.ts`, `tree-viewer.ts` | xamgle `Core: d4: Filters: Removed column`, `Fixed issues`; `sunburst`, `tree-on-click-filter` |
| core#1223: a hierarchical card kept a removed level, also from a layout | `hierarchical_filter.dart` `refresh`, `applyState` | xamgle `Core: d4: Filters: Removed column` |
| core#1194, core#1195: Shift+Enter selected rows with Allow Row Selection off (and opened the cell editor); a blank group name made a group | `grid_row_selection.dart`, `grid_editors.dart` `isEditButton`; xamgle `column_commands.dart` | `viewers/grid/grid-rows`, `grid-column-groups` |
| core#1137, core#1191: Unpin Column threw with Frozen Columns above the column count; pinned values past the 10-row cap came back from a layout | `grid_core.dart` (`_clampFrozenColumns`, `pinRowsByValues`), `grid_popup_menu.dart` | xamgle `Core: d4: Viewers: Grid` |
| core#1213: a viewer's Filter formula kept a renamed column's old name | `client/d4/lib/src/viewer_base/data_frame_viewer.dart` `_renameColumnInLook` | ApiTests `AI: Viewers: Shared rules` |
| core#1252, core#980, core#1180: a layout renamed to nothing in its Details pane; a Share dialog's OK with nobody picked; a profile name empty or over 64 characters | xamgle `view_layout_meta.dart`, `share_dataset.dart` + `permissions_browser.dart`, `user_profile.dart` | `layouts/layout-names-tags`, `sharing/sharing-dialog`, `general/profile` |
| public#4473, public#4475: MCL on number columns never finished; Compare Columns over two tables was false in every row | libraries/ml `mcl-editor.ts`; PowerPack `compare-columns.ts` | EDA `analyze/clustering`, PowerPack `compare-columns` |
| core#641, core#1246: Run query... on an unsaved query; Data to Curves and Compare Columns missing from the Data menu | fixed in core before the round (8fad96d965, 2f1fee3171) | `queries/parameterized-queries`, Curves and PowerPack features |
| The Recent view's calendar offered a placeholder day for 500 ms (browse-my-stuff red on master after the Dart 3 client); its feed failed with a server RangeError when a package's log message was among the user's events | xamgle `recent_activity_view.dart`: no placeholder, the dates load at once; `server/datlas/lib/src/services/action_logger_service.dart` `getActivity` skips an event that names fewer than two entities | `browse/browse-my-stuff` |
| The Move entity dialog of a drop on a project or space opened before its rows, which wait for each entity to be fetched again: under two workers "selects Move" found no choice (projects-move, projects-complex, now and then) | xamgle `meta/project_meta.dart` `entityProjectMove` shows the dialog once the move widget is built | xamgle `Core: xamgle: Projects/Layouts: Projects` 'The Move entity dialog opens with its choices'; the nine features that use the dialog |
| A connection a local docker spawner serves (the ChEMBL container, PostgresDocker) failed: grok_connect, itself in docker, was handed `localhost:<port>` (the 2026-10-29 fix, lost in the Dart 3 merge) | `server/datlas/lib/src/services/connectors/connectors_service.dart` `startConnectionDocker` | `viewers/trellis-plot/trellis-plot-keyboard`; `connections/connections-schema`, which skipped at its gate before |

## Open on 2026-10-12

| Scenario | Pinned at | Was |
| --- | --- | --- |
| [Converted molecules keep their stereochemistry — GROK-20956](../../packages/Chem/bdd/features/transform/notation.feature) | line 58, `every molecule of "canonical_smiles_molblock" column should be the same as in "canonical_smiles" column` | Rows 478, 487 and 488: RDKit's own SMILES → molblock → SMILES round trip flips the vinyl stereocentre of the quinuclidines in Chem's RDKit_minimal 1.2.23 (2024.09) and keeps it in 2026.03; the tag goes with the upgrade. CoordGen also flips the C=N geometry of row 833 — a separate defect the next line would report once the first is fixed. |

## Validation

Final complete suites on four workers, 2026-09-13:

| Suite | Passed | Time |
| --- | ---: | ---: |
| UsageAnalysis, including all viewers and Spaces | 135 | 129.4 s |
| EDA | 13 | 18.1 s |
| Bio | 18 | 26.2 s |
| DiffStudio | 9 | 31.3 s |
| Peptides | 1 | 4.7 s |
| PowerGrid | 5 | 10.4 s |
| U2Demo | 64 | 14.6 s |
| BDD library platform features | 3 | 5.9 s |
| **Total** | **248** | |

The JSON reports contain no unexpected failures, skipped tests, flakes or retries. Eleven
scenario assertions failed as expected, each at the intended step listed above; those expected
failures are included in the passing journey totals.

The UsageAnalysis known-failure journeys plus the corrected line-chart journey also passed two
serial repetitions (18 tests), three repetitions on four workers (27), and a headed run together
with all Spaces features (19). Space file operations and drag/drop passed two serial repetitions
(4). Box-plot statistics passed two serial repetitions, eight on four workers, and one headed
run (11). These additional runs had no unexpected failures, skips, flakes or retries.

Library build, generated-spec checks, type checking and all 64 unit tests passed, including the
Chromium regression tests. The legacy d4 analyzer reported 16 existing hints and no errors or
warnings. After the runs, no named Spaces fixtures remained on the server, and the TSLA and
acidiq demo files matched Git.

## Rechecking

From a package directory, after building the library and serving the matching core:

```bash
export DATAGROK_URL=http://localhost:8889
PLAYWRIGHT_JSON_OUTPUT_FILE=/tmp/bdd-known-failures.json npx grok-bdd run --workers=1 --grep @known-failure --reporter=list,json
```
(Or whatever the port is on that stand (8889 or 8888))
Inspect the tagged scenario's step error in the JSON report. A journey may contain several tagged
scenarios, so test totals and known-failure totals differ. The corrected line-chart journey no
longer matches this grep; run `generated/viewers/line-chart/line-chart-selection` separately.
