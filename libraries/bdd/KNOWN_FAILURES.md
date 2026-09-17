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

## Retained product defects

Core paths below are relative to `core/client/`; public source links are relative to this file.
Expected values remain unchanged.

| Scenario | Observed failure | Cause |
| --- | --- | --- |
| [Correlation tooltip disabled](../../packages/UsageAnalysis/bdd/features/viewers/correlation-plot/correlation-plot-cells.feature) | A visible tooltip still contains `Pearson R` with `showTooltip=false`. | `d4/lib/src/viewers/correlation_plot/correlation_plot_core.dart` supplies a custom cell tooltip without checking the flag. The nested grid must also suppress its default fallback when the tooltip is disabled. |
| [Correlation color scale](../../packages/UsageAnalysis/bdd/features/viewers/correlation-plot/correlation-plot.feature) | Correlations approximately 0.0648 and 0.4124 have the same full-red color. | The plot sets column bounds to -1 and 1, but the numeric path in `d4/lib/src/common/color_coding.dart` passes null bounds to linear coloring. |
| [Forms after header sort](../../packages/UsageAnalysis/bdd/features/viewers/forms/forms-core.feature) | 0 visible cards instead of 7. | Sorting clears the current row. [Forms](../utils/src/viewers/forms-viewer.ts) creates an empty leading card; `d4/lib/src/widgets/virtual_item_view.dart` measures its zero height and computes zero layout columns. |
| [Forms with Show Current Row off](../../packages/UsageAnalysis/bdd/features/viewers/forms/forms-interactions.feature) | 0 visible cards instead of 6. | The same layout defect: a blank mouseover card becomes the first measured item. Both tags remain because they cover different user actions. |
| [Heatmap Colors off — GROK-20619](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map-colors.feature) | The AGE band changes by 0 pixels. | The dense path at row height ≤5 in `d4/lib/src/viewers/grid/grid_core.dart` bypasses the normal cell renderer's `heatmapColors` gate and colors cells directly. The normal renderer does honor the flag. |
| [Returning to heat-map mode](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map-navigation.feature) | Row height stays near the grid's 28 px instead of fitting the table within 0–8 px. | `GridLook.refreshGrid` and `GridCore.onLookChanged` preserve the grid's vertical viewport instead of restoring the full heat-map range. |
| [Heat-map Row Height disabled](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map.feature) | The Row Height property exists but is not disabled. | `libs/property_grid/lib/property_grid_lib.dart` skips initial dependency evaluation when the controlling property has no editor. `isGrid` is non-editable, although its value is available on the look. |
| [Heat-map column cap before settings](../../packages/UsageAnalysis/bdd/features/viewers/heat-map/heat-map.feature) | Setting the cap to 3 leaves 11 columns visible. | The API does invoke the Dart setter. `GridLook.refreshGrid` returns because the fresh look has no `viewer` reference; opening settings binds that reference, after which the same API write works. The feature's earlier API-versus-setter explanation was corrected. |
| [PC plot selection after reset — GROK-17306](../../packages/UsageAnalysis/bdd/features/viewers/pc-plot/pc-plot-transformation.feature) | 0 selected rows instead of the prior 10. | Reset rebuilds the transformation and detaches the old frame link. Link cleanup clears the source selection, and the replacement link does not restore it. The assertion reads the original table. |
| [Filtered group comparison — GROK-20795](../../packages/EDA/bdd/features/analyze/filtered-group-comparison.feature) | The first result count is 157 instead of 104; full-table counts are 157/5266/354 instead of filtered 104/2823/279. | [Control comparisons](../../packages/EDA/src/control-comparisons/control-comparisons-ui.ts) passes the original columns to factorization without applying the table's filter. Fixture counts were checked against demog, including missing AGE values and the excluded control group. |
| [Empty Pareto objective](../../packages/EDA/bdd/features/pareto-front-objectives.feature) | The picker offers 17 columns instead of 16, including the all-null integer column. | [Pareto viewer](../../packages/EDA/src/pareto-optimization/pareto-front-viewer.ts) filters the picker by numerical type only. Its separate initialization check for nonempty columns is not applied to the picker. |

The retained defects have implementation causes independent of keyboard differences between
macOS and Windows. Their tags remain assertions of the desired behavior, not skips.

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
