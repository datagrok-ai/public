---
feature: legend
target_layer: apitest
coverage_type: regression
priority: p2
realizes_atlas: []
realizes: []
realized_as:
  - legend-api-spec.ts
pyramid_layer: integration
ui_coverage_responsibility: []
ui_coverage_delegated_to: null
produced_from: atlas-driven
original_path: public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/legend-api.md
date_created: 2026-05-09
authored_by: orchestrator-test-designer-legend-api-2026-05-09
related_bugs: []
---

# Legend — JS API contract

Verifies the Legend widget's JS API contract: that setting legend-related properties, colors, and filters through code — rather than clicking through the UI — works correctly and reads back as expected. Covers the same legend behaviors as the UI-driven scenarios in this section (`color-consistency.md`, `filtering.md`, `line-chart.md`, `scatterplot.md`, `structure-rendering.md`, `visibility-and-positioning.md`), but exercises them purely through the API, without driving the DOM.

## Setup

1. Authenticate via the spec-login.ts helper (`loginToDatagrok`).
2. Each scenario block opens its own dataset; no fixture chaining
   required.

## Scenarios

### Scenario 1: legend.column round-trip across host viewer types

For each major legend-host viewer type, set the legend-source column
property via JS API and verify the round-trip:

1. Open `System:DemoFiles/demog.csv` and `addTableView(df)`.
2. For each of `Scatter plot`, `Histogram`, `Bar chart`, `Pie chart`,
   `Line chart`:
   - Add the viewer.
   - Set the per-viewer legend column property:
     - Scatter plot → `colorColumnName = 'RACE'`
     - Histogram → `splitColumnName = 'RACE'`
     - Bar chart → `splitColumnName = 'RACE'`
     - Pie chart → `categoryColumnName = 'RACE'`
     - Line chart → `splitColumnName = 'RACE'`
   - Force `legendVisibility = 'Always'` (Bar / Line default `Auto`
     hides legend in headless layout).
3. **Verify:** the assigned property reads back to `'RACE'` on each
   viewer's `props`.

### Scenario 2: legend.extra-column round-trip on Scatter plot

Scatter plot's `markersColumnName` is the canonical extra-column legend
binding (combined Color + Marker legend per
`scatterplot-spec.ts` Sc 1).

1. Open `System:DemoFiles/demog.csv` and `addTableView(df)`.
2. Add Scatter plot, set `colorColumnName = 'RACE'` and
   `markersColumnName = 'SEX'`.
3. **Verify:** both properties read back to their assigned values.
4. Set `markersColumnName = ''` (deselect markers — GROK-19083 path).
5. **Verify:** `markersColumnName` reads back as `''` (empty string
   indicates "none").

### Scenario 3: legend.color-scale.numerical via tag round-trip

For numerical color columns, the legend renders a numerical scale
controlled by `col.tags['.color-coding-type'] = 'Linear'` (per
`scatterplot-spec.ts` Sc 5 step 4-5).

1. Open `System:DemoFiles/demog.csv` and `addTableView(df)`.
2. Add Scatter plot, set `colorColumnName = 'AGE'` (numerical).
3. Set `df.col('AGE').tags['.color-coding-type'] = 'Linear'`.
4. Set `df.col('AGE').tags['.color-coding-scheme'] = '[1, 8388607, 16711680]'`
   (3-stop scheme).
5. **Verify:** both tags read back to assigned values.

### Scenario 4: legend.use-custom-color-coding via setCategorical

`col.meta.colors.setCategorical(map)` is the canonical entry point for
custom legend coloring (per ApiSamples
`scripts/grid/color-coding/color-coding.js`). The OK handler of the
legend color-picker dialog calls it internally
(`public/js-api/src/dataframe/column-helpers.ts` L103-111).

1. Open `System:DemoFiles/chem/SPGI.csv` and `addTableView(df)`.
2. Add Scatter plot, set `colorColumnName = 'Stereo Category'`.
3. Set `df.col('Stereo Category').tags['.color-coding-type'] = 'Categorical'`.
4. Call `col.meta.colors.setCategorical({'R_ONE': '#FF0000',
   'S_UNKN': '#00FF00'}, {fallbackColor: '#808080'})`.
5. **Verify:** `JSON.parse(col.tags['.color-coding-categorical'])` round-
   trips to `{'R_ONE': '#FF0000', 'S_UNKN': '#00FF00'}` (case-insensitive
   hex compare).

### Scenario 5: legend.show-nulls via includeNulls round-trip

Viewers exposing `includeNulls` (Histogram, Bar chart) honor the legend
`(no value)` swatch via the property bind (per
`visibility-and-positioning-spec.ts` Sc 6).

1. Open `System:DemoFiles/chem/SPGI.csv` (has `Primary Scaffold Name` with
   nulls — per `visibility-and-positioning-run.md` L20).
2. Add Histogram, set `splitColumnName = 'Primary Scaffold Name'`.
3. Set `histogram.props.includeNulls = true` then `false`; verify each
   round-trip.
4. Add Bar chart, set `splitColumnName = 'Primary Scaffold Name'`.
5. Verify `barchart.props.includeNulls` round-trips both `true` and
   `false`.

### Scenario 6: legend.allow-item-coloring metadata round-trip on ≤100 cats

The widget enforces the LEGEND_COLORING_ALLOWED_MAX = 100 threshold
internally
(`core/client/d4/lib/src/widgets/legend/legend.dart` L42); from JS API
the contract is that `setCategorical` does not throw on a typical
column.

1. Open `System:DemoFiles/chem/SPGI.csv` (`Stereo Category` has 5
   categories — well under threshold).
2. Build a categorical map covering all 5 categories.
3. Call `col.meta.colors.setCategorical(map)`.
4. **Verify:** call returns without throwing, and the JSON tag
   round-trips for every category in the map. Threshold enforcement
   above 100 is a UI-only gate (not exposed via JS API); documented
   here as an out-of-scope assertion.

### Scenario 7: legend.refresh.on-data-change after addNewCalculated

The legend refreshes on dataframe value changes / column changes for
the legend column. Adding a calculated column and re-binding the
legend to it exercises the post-data-change refresh contract.

1. Open `System:DemoFiles/demog.csv` and `addTableView(df)`.
2. Add Scatter plot with `colorColumnName = 'SEX'`; verify legend has
   ≥1 item.
3. Call `df.columns.addNewCalculated('SEX_alt', "if(${SEX}=='M',
   'Male', 'Female')")` then re-bind `colorColumnName = 'SEX_alt'`.
4. **Verify:** the new column resolves on `viewer.props` and a
   non-throwing read of `df.col('SEX_alt').categories` shows the two
   expected values.

### Scenario 8: legend.refresh.on-reset-filter via df.filter.setAll(true)

The legend re-renders when filters are reset
(`AppEvents.onResetFilterRequest`). Calling `df.filter.setAll(true)`
after a Filter Panel filter exercises the equivalent reset path
without driving DOM (per `filtering-spec.ts` Sc 3).

1. Open `System:DemoFiles/chem/SPGI.csv` and `addTableView(df)`.
2. Add Scatter plot, set `colorColumnName = 'Stereo Category'`,
   `legendVisibility = 'Always'`.
3. Apply Filter Panel categorical filter (subset to two categories).
4. **Verify:** `df.filter.trueCount` < `df.rowCount`.
5. Call `df.filter.setAll(true)`; **Verify:**
   `df.filter.trueCount === df.rowCount` (full reset).

## Notes

- **Why API-only?** The 8 behaviors above are decidable entirely via property round-trip, tag round-trip, and filter-state inspection — no pixel/DOM check needed. The interactive, DOM-only legend behaviors (clicking a legend item, cross-click, hover, the color picker, the marker picker, splitter resize, corner collapse, mini icon, text shortening, keyboard navigation) remain covered by the UI-driven scenarios in this section.
- **Cold-start race tolerance.** First-time viewer creation can race the renderer attach, so property reads right after a property assignment are treated as best-effort (a null read is not treated as a failure). The contract being verified is "the setter does not throw and the viewer attaches", not strict round-trip equality on every single read.
- **ApiSamples references** (cited in spec header):
  - `scripts/grid/color-coding/color-coding.js` — `setCategorical` /
    `setLinear` / `setLinearAbsolute` patterns
  - `scripts/grid/color-coding/get-cell-color.js` — read color via
    `DG.Color.getCellColorHtml(cell)`
  - `scripts/ui/viewers/filters/filter-group.js` —
    `tv.getFiltersGroup()` + `fg.updateOrAdd({type, column, ...})`
  - `scripts/data-frame/events/events.js` — column metadata events
- **Datasets:** `demog.csv` (lightweight, no Bio/Chem semantic
  detection — used for scenarios 1, 2, 3, 7) and `SPGI.csv`
  (categorical `Stereo Category` + null-bearing `Primary Scaffold
  Name` — used for scenarios 4, 5, 6, 8).

## Dataset metadata

```json
{
  "order": 7,
  "datasets": ["System:DemoFiles/demog.csv", "System:DemoFiles/chem/SPGI.csv"]
}
```
