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

No `@known-failure` tag remains in any project.

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
