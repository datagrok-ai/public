# Viewer API research — for `datagrok-viewers` skill (topic 4)

Citations use `<absolute-path>:<line>`. Path abbreviations:

- `js-api/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/js-api/src/`
- `samples/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/ApiSamples/scripts/`
- `chem/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Chem/src/`
- `grokky/` = `/Users/aleksashka/Desktop/datagrok/reddata/public/packages/Grokky/`
- `skill/` = `grokky/dockerfiles/claude-runtime/plugin/skills/datagrok-viewer/`

The viewer surface is the largest of the four topics so far: ~30 built-in types, multiple
creation paths, a property layer with three different access shapes, and dock-aware lifecycle.
The existing `datagrok-viewer` skill (`skill/SKILL.md`) and `skill/viewer-properties.md`
already address most of the *content* — type strings, property conventions, fuzzy matching,
schema-driven validation. What they miss is **execution context**: how to behave when `view`
is in scope, how to close viewers, how to find existing ones, how to handle plugin viewers.
This research is structured to feed a refreshed skill that fixes the `grok.shell.tv` bug,
exposes `closeViewer`/`findViewer`, and explicitly delineates built-in vs plugin viewers.

## 1. Summary table

| # | Area | Key entrypoints | Canonical idiom | Wrapper-worthiness |
|---|------|-----------------|-----------------|--------------------|
| A | Creating viewers | `Viewer.fromType` (`js-api/viewer.ts:128`), `DG.Viewer.scatterPlot/histogram/...` (`viewer.ts:215-306`), `TableView.addViewer` (`js-api/views/view.ts:417`), `view.scatterPlot()` shorthands (`view.ts:447-606`, marked deprecated), `df.plot.scatter()` (`js-api/dataframe/data-frame.ts:551-571`) | `view.addViewer('Scatter plot', {x: 'a', y: 'b'})` when on a `TableView`; `DG.Viewer.fromType(type, df, opts)` for off-view rendering | **High** — fix the `grok.shell.tv` bug, accept `view`-in-scope |
| B | Viewer type catalog | `DG.VIEWER` enum (`js-api/const.ts:671-704`), `DG.CORE_VIEWER` (`const.ts:707-733`), `Viewer.CORE_VIEWER_TYPES` (`viewer.ts:363-372`), `Viewer.getViewerTypes()` (`viewer.ts:135`) | 30 built-in entries; `Viewer.getViewerTypes()` returns those + plugin-registered | **Medium** — extend the existing per-viewer property table |
| C | Configuring viewers — properties | `viewer.setOptions(map)` (`viewer.ts:146`), `viewer.getOptions(includeDefaults?)` (`viewer.ts:158`), `viewer.props.<name>` (`viewer.ts:118` `ObjectPropertyBag`), `viewer.getProperties()` returns `Property[]` (`viewer.ts:166`) | `viewer.setOptions({...})` for batch; `viewer.props.x = 'a'` for single | **High** — already wrapped (`grokky.configureViewer`); keep and harden |
| D | Viewer-to-viewer linking | Implicit via shared `DataFrame` — `viewer.dataFrame` getter (`viewer.ts:208`); cross-DF via `grok.data.linkTables` (out of scope) | Don't ask. Same df → linked. | **None** — document only |
| E | Finding & closing viewers | `view.viewers` iterator (`view.ts:625`), `viewer.close()` (`viewer.ts:171`), `viewer.removeFromView()` (`viewer.ts:355`), `view.detachViewers()` (`view.ts:619`), `view.resetLayout()` (`view.ts:609`) | `Array.from(view.viewers).slice(1).forEach(v => v.close())` to "close all but the grid" (`UITests/src/viewers/viewers.ts:251-254`) | **High** — `closeViewer(target)`, `findViewers(view, pred)`, `closeAllViewers(view, {keepGrid})` |
| F | Viewer events | `grok.events.onViewerAdded` / `onViewerClosed` (`js-api/events.ts:202-204`); `viewer.onContextMenu`, `viewer.onDataSelected`, scatter-plot `onPointClicked` etc (`viewer.ts:308-309`, `:667-668`) | Subscribe via rxjs; out of scope for one-shot exec | **Low** — document only |
| G | Saving viewer state / layouts | `viewer.getOptions(true)` (`viewer.ts:158`), `view.saveLayout()` / `view.loadLayout()` (`js-api/views/view.ts:296-305`), `view.saveState()` / `view.loadState()` (`view.ts:632-637`) | `tv.saveLayout()` to capture all viewers + grid; full coverage in T1.5 | **None here** — cross-reference |
| H | Dock layout | `view.dockManager` (`view.ts:439`), `DockManager.dock(el, type, refNode?, title?, ratio?)` (`js-api/docking.ts:135`), `DOCK_TYPE` enum (`const.ts:776-782`), `view.dockManager.close(el)` (`docking.ts:143`) | `view.addViewer(...)` first; then `view.dockManager.dock(viewer, 'right', null, 'Title', 0.4)` to reposition | **Medium** — touch lightly here (one example), full skill later |
| I | Filter viewer (special) | `view.getFiltersGroup({createDefaultFilters?})` (`view.ts:404`), `FilterGroup.updateOrAdd(state, requestFilter?)` (`viewer.ts:525`) | Use **filtering** skill, not this one | **None** — cross-reference |
| J | Existing `grokky.addViewer` / `configureViewer` | `grokky/src/claude/grokky-api.ts:162-180` | Fuzzy type, schema-validated properties, column-name suffix shortcut, did-you-mean warnings | **High** — fix `tv` bug (`:163`); document everything else |
| K | "Close scatter, add histogram" idiom | `findViewers(view, v => v.type === DG.VIEWER.SCATTER_PLOT).forEach(v => v.close()); view.addViewer('Histogram', ...)` | One canonical pattern needed in the skill | **High** — `closeViewer({type})` is the lift |
| L | Plugin viewers (`Chem:chemSpace` etc.) | `grok.functions.call('Chem:chemSpaceTopMenu', {...})` returns `Promise<DG.Viewer>` (`chem/package.ts:821-877`); fallback `DG.Func.find({name:'...'})[0].prepare(args).call()` (`chem/package.ts:841`) | Function call, not viewer instantiation; returns a `Viewer` already added to `grok.shell.tv` | **Medium** — separate helper or cross-link to `callRegisteredFunction` |

## 2. Creating viewers (Area A)

Four paths exist. Two are blessed, two are legacy/specialty.

### 2.1 `DG.Viewer.fromType(viewerType, df, options?)` — generic constructor

`js-api/viewer.ts:128-130`:

```ts
static fromType(viewerType: ViewerType, table: DataFrame, options: object | null = null): Viewer {
  return toJs(api.grok_Viewer_FromType(viewerType, table.dart, _toJson(options)));
}
```

Synchronous, returns a `Viewer` *not yet attached* to any view. Caller is responsible for
placing the viewer (`view.addViewer(v)` or `view.dockManager.dock(v, ...)`). Used in
`samples/ui/viewers/create-viewers-dynamically.js:6` and `samples/ui/docking/docking-table-view.js:8-11`:

```js
let scatterPlot = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t);
// ...
let viewer = DG.Viewer.fromType('Scatter Plot', table);
view.addViewer(viewer);
view.dockManager.dock(viewer, 'right');
```

There's also an async variant exposed via `df.plot.fromType(viewerType, options)` (`js-api/dataframe/data-frame.ts:557`) that uses `grok_Viewer_FromType_Async` — required for viewers that load on demand (plugin viewers, non-core types).

### 2.2 Typed shorthands — `DG.Viewer.scatterPlot(df, opts)` etc.

`js-api/viewer.ts:215-306` declares static factories for every core viewer:

- `Viewer.grid(t, opts?)` → `Grid` (`:215`)
- `Viewer.histogram(t, opts?)` → `HistogramViewer` (`:219`)
- `Viewer.barChart(t, opts?)` → `Viewer<IBarChartSettings>` (`:223`)
- `Viewer.heatMap(t, opts?)` → `Grid` (`:227` — note: returns `Grid`, not a Viewer subclass)
- `Viewer.boxPlot(t, opts?)` → `BoxPlot` (`:231`)
- `Viewer.filters(t, opts?)` → `Viewer<IFiltersSettings>` (`:235`)
- `Viewer.scatterPlot(t, opts?)` → `ScatterPlotViewer` (`:239`)
- `Viewer.lineChart(t, opts?)` → `LineChartViewer` (`:243`)
- `Viewer.network(t, opts?)` (`:247`), `calendar` (`:251`), `correlationPlot` (`:255`),
  `densityPlot` (`:259`), `form` (`:263`), `markup` (`:267`), `matrixPlot` (`:271`),
  `pcPlot` (`:275`), `pieChart` (`:279`), `scatterPlot3d` (`:283`), `statistics` (`:287`),
  `tile` (`:291`), `treeMap` (`:295`), `trellisPlot` (`:299`), `wordCloud` (`:304` — `@deprecated`)

These return strongly typed subclasses (e.g. `ScatterPlotViewer` has `zoom()`, `worldToScreen()`,
`onPointClicked` — see `viewer.ts:610-669`). The shorthand is **detached** by default — same as
`fromType`; you must `view.addViewer(v)` to mount it.

### 2.3 `TableView.addViewer(type | Viewer, options?)` — the right ergonomic call

`js-api/views/view.ts:417-426`:

```ts
addViewer(v: ViewerType | string | Viewer, options?: any): Viewer {
  if (typeof v === 'string')
    v = toJs(api.grok_View_AddViewerByName(this.dart, v)) as Viewer;
  else
    api.grok_View_AddViewer(this.dart, v.dart);
  if (options)
    v.setOptions(options);
  api.grok_TableView_ProcessNewViewer(this.dart, v.dart);
  return v;
}
```

This is the **canonical idiom** for `datagrok-exec` blocks. Pass a string (e.g.
`'Scatter plot'` or `DG.VIEWER.SCATTER_PLOT`), and the view both creates and docks it,
applying options through `setOptions()`. Test coverage confirms it works with every core
type: `ApiTests/src/ai/reported-issues/grok-19970.ts:24-25` passes options through cleanly:

```ts
const v = tv.addViewer(DG.VIEWER.BAR_CHART, {value: 'age', category: 'race', onClick: 'Filter'});
```

When passed a pre-built `Viewer` (via `DG.Viewer.fromType` or a static factory), `addViewer`
attaches it to the view — see `docking-table-view.js:8-9`. **One subtle catch**: the typed
static factories like `Viewer.scatterPlot(t)` return a *detached* viewer associated with `t`;
calling `view.addViewer(sp)` mounts it but uses **`t`** as the source, not `view.dataFrame`.
If `t !== view.dataFrame`, you get a viewer that ignores the view's table — usually a bug.

### 2.4 `TableView.scatterPlot(opts)` / `histogram(opts)` — deprecated shorthands

`view.ts:443-606` declares per-type shorthands like `scatterPlot(opts)`, `histogram(opts)`,
`barChart(opts)`. Every one of them carries an explicit "deprecated: use
`addViewer(Viewer.scatterPlot(options))`" comment (e.g. `view.ts:443-447`, `:451-455`,
`:558-562`). The samples still use them everywhere (`samples/ui/viewers/types/*.js`) because
they were authored before the deprecation, but new code should call `addViewer(type, opts)`.

The skill should write `view.addViewer('Scatter plot', {...})` exclusively — it works for
every type, the option route is identical, and it survives the deprecation cycle.

### 2.5 `df.plot.<kind>()` — fluent dataframe API

`js-api/dataframe/data-frame.ts:551-571` adds a small `DataFramePlotHelper`:

```ts
scatter(options?) { return DG.Viewer.scatterPlot(this.df, options); }
grid(options?)    { return DG.Viewer.grid(this.df, options); }
histogram(options?) { return DG.Viewer.histogram(this.df, options); }
bar(options?)     { return DG.Viewer.barChart(this.df, options); }
heatMap(options?) { return DG.Viewer.heatMap(this.df, options); }
box(options?)     { return DG.Viewer.boxPlot(this.df, options); }
line(options?)    { return DG.Viewer.lineChart(this.df, options); }
network(options?) { return DG.Viewer.network(this.df, options); }
async fromType(viewerType, options?) { ... }  // async — waits for viewer to load
```

Same detached-viewer semantics as `DG.Viewer.scatterPlot(df)`. Useful only for off-view
rendering (e.g. `t.plot.line(...).root` to inject into a dialog — see
`samples/grid/advanced/charts-in-cells.js:29-31`). Not worth surfacing in the skill.

### 2.6 What about non-`TableView` views?

`view` (the global in `datagrok-exec` blocks) is typed `DG.ViewBase` per the skill prologue.
In practice, when the user is on a table, `view` is `DG.TableView` (which extends `View`
extends `ViewBase`). If `view` is a plain `View` (e.g. a script view, a function view), it
does **not** have `addViewer` — the method is only on `TableView` (`view.ts:417`).

**Decision rule** for the helper:
1. If `view instanceof DG.TableView` (or has `addViewer` + `dockManager`), use it.
2. Else if `grok.shell.tv` exists, fall back with a console warning.
3. Else throw `'addViewer: no active TableView'`.

The current `grokky.addViewer` skips step 1 entirely and goes straight to `grok.shell.tv`
(`grokky/src/claude/grokky-api.ts:163`). That breaks any script run from a non-active
TableView, and breaks the spirit of the `view`-in-scope convention.

## 3. Viewer type catalog (Area B)

### 3.1 The enum

`js-api/const.ts:670-704` — `DG.VIEWER`:

```ts
HISTOGRAM = 'Histogram',
BAR_CHART = 'Bar chart',
BOX_PLOT = 'Box plot',
CALENDAR = 'Calendar',
CORR_PLOT = 'Correlation plot',
DENSITY_PLOT = 'Density plot',
FILTERS = 'Filters',
FORM = 'Form',
GLOBE = 'Globe',
GRID = 'Grid',
GOOGLE_MAP = 'Google map',
HEAT_MAP = 'Heat map',
LINE_CHART = 'Line chart',
SHAPE_MAP = 'Shape Map',
MARKUP = 'Markup',
MATRIX_PLOT = 'Matrix plot',
NETWORK_DIAGRAM = 'Network diagram',
PC_PLOT = 'PC Plot',
PIE_CHART = 'Pie chart',
SCATTER_PLOT = 'Scatter plot',
SCATTER_PLOT_3D = '3d scatter plot',
STATISTICS = 'Statistics',
TILE_VIEWER = 'Tile Viewer',
TREE_MAP = 'Tree map',
TRELLIS_PLOT = 'Trellis plot',
WORD_CLOUD = 'Word cloud',
TIMELINES = 'Timelines',
RADAR_VIEWER = 'Radar',
SURFACE_PLOT = 'Surface plot',
SCAFFOLD_TREE = 'Scaffold Tree',
PIVOT_TABLE = 'Pivot table',
CONFUSION_MATRIX = 'Confusion matrix'
```

A second enum `CORE_VIEWER` (`const.ts:707-733`) is the subset that's always loaded
synchronously — same list minus `Globe`, `Google map`, `Timelines`, `Radar`, `Surface plot`,
`Scaffold Tree`, `Word cloud`. `Viewer.CORE_VIEWER_TYPES` at `viewer.ts:363-372` is a third,
similar list (note: this one includes `WORD_CLOUD` but excludes `GOOGLE_MAP` — slight drift
from `CORE_VIEWER` enum). For the skill we should treat all of `DG.VIEWER` as legal type
strings, with the proviso that the non-core ones may need an async load.

### 3.2 Built-in catalog with primary properties

The existing `skill/viewer-properties.md` covers 8 viewers (scatter, 3D scatter, histogram,
bar, line, box, trellis, pie, form) in detail. Below is a one-line description plus the most
common props for **every** entry. Use as a "starter set" — the wrapper still validates via
`getProperties()` at runtime.

| Type string | Description | Common properties |
|---|---|---|
| `Scatter plot` | 2D x-vs-y numerical | `x`, `y`, `color`, `size`, `markerType`, `showRegressionLine`, `xAxisType`, `yAxisType`, `xMin..yMax`, `showXHistogram`, `zoomAndFilter` |
| `3d scatter plot` | three numerical axes | `x`, `y`, `z`, `color`, `size`, `label`, `xAxisType..zAxisType`, `showAxes`, `backColor` |
| `Histogram` | single-column distribution | `value`, `bins`, `split`, `splitStack`, `normalizeValues`, `valueMin`, `valueMax`, `color`, `colorAggrType`, `showBinSelector` |
| `Bar chart` | aggregates per category | `value`, `valueAggrType`, `split`, `stack`, `orientation`, `axisType`, `relativeValues`, `barSortType`, `onClick` |
| `Line chart` | series over an ordered axis | `x`, `split`, `splitColumnNames`, `yAggrTypes`, `segment`, `multiAxis`, `xAxisType`, `xMin..yMax`, `showVerticalGridLines` |
| `Box plot` | distribution by category | `categoryColumnNames`, `value`, `axisType`, `showStatistics`, `binColor`, `markerColor`, `colorAxisType`, `invertYAxis` |
| `Pie chart` | categorical proportions | `category`, `segmentAngle`, `segmentLength`, `segmentAngleAggrType`, `pieSortType`, `includeNulls`, `labelPosition`, `onClick` |
| `Trellis plot` | small multiples | `xColumnNames`, `yColumnNames`, `viewerType`, `innerViewerLook`, `globalScale`, `useTiledView`, `tilesPerRow`, `autoLayout` |
| `Grid` | tabular spreadsheet | accessed via `view.grid.setOptions({...})` — props like `rowHeight`, `colHeaderHeight`, `showColumnLabels` |
| `Heat map` | matrix of values (returns a `Grid` subclass — `viewer.ts:227`) | column tag-driven color coding; props from `IGridSettings` |
| `Density plot` | smoothed 2D density | `xColumnName`, `yColumnName`, `binsX`, `binsY`, `palette` |
| `Filters` | filter panel as viewer | `columnNames`, `showHeader`, `showBoolCombinedFilter` |
| `Form` | per-row record view | `columnNames`, `title`, `showNavigation`, `showRowSelector`, `formFont`, `descriptionPosition` |
| `Calendar` | date-indexed events | `dateColumnName`, default sufficient for demo |
| `Globe` | 3D geo globe | `latitudeColumnName`, `longitudeColumnName` |
| `Google map` | flat geo map | `latitudeColumnName`, `longitudeColumnName` |
| `Shape Map` | choropleth on regions | needs region-id column + value column |
| `Matrix plot` | scatter matrix | `xColumnNames`, `yColumnNames` |
| `Network diagram` | node-link graph | `node1ColumnName`, `node2ColumnName`, `edgeColorColumnName` |
| `PC Plot` | parallel coordinates | `columnNames` |
| `Correlation plot` | correlation matrix | `xs`, `ys` (column-name arrays) — see `samples/ui/viewers/types/corr-plot.js:5-8` |
| `Statistics` | numeric column summary | `columnNames` |
| `Surface plot` | 3D surface from gridded data | `xColumnName`, `yColumnName`, `zColumnName` |
| `Tile Viewer` | one card per row | `cardLayout`, `pageSize` |
| `Tree map` | hierarchical area chart | `categoryColumnName`, `sizeColumnName` |
| `Markup` | static HTML/markdown | `content` (string) |
| `Word cloud` | term frequency | `column` |
| `Timelines` | events on time axis | `startColumnName`, `endColumnName`, `categoryColumnName` |
| `Radar` | per-row radar | `columnNames` |
| `Scaffold Tree` | chem scaffold hierarchy | plugin-provided in practice |
| `Pivot table` | grouped aggregates | `groupByColumns`, `aggregations` |
| `Confusion matrix` | classifier prediction vs actual | `actualColumnName`, `predictedColumnName` |

### 3.3 Discoverability at runtime

`Viewer.getViewerTypes()` (`viewer.ts:135`) returns the live list including plugin viewers.
`DG.WidgetDescriptor.getDescriptors()` (`viewer.ts:40`) returns one descriptor per registered
widget, each with `.name`, `.description`, `.synonyms`, `.properties` — usable for AI ranking
without instantiating the viewer (`samples/ui/viewers/advanced/descriptors.js:6-12`).

`Viewer.canVisualize(viewerType, df)` (`viewer.ts:359-361`) returns `null` if the viewer can
visualize this df, or a string with the reason it can't (e.g. "no numerical column"). Useful
for graceful failure but never observed in the samples.

## 4. Configuring viewers — properties (Area C)

Three property-access shapes coexist:

```ts
viewer.setOptions({x: 'a', y: 'b'});        // batch, takes a plain object
viewer.props.xColumnName = 'a';             // ObjectPropertyBag proxy, one prop at a time
viewer.props['xColumnName'] = 'a';          // same; bracket access
const p: Property[] = viewer.getProperties();  // metadata: name, type, semType, description
```

`setOptions` is documented at `viewer.ts:146-148`. It calls `grok_Viewer_Options` which
**accepts non-canonical property names** (the Dart side does some matching). Confirmed by
`samples/ui/viewers/types/trellis-plot.js:5-16`, `scatter-plot.js:5-10`, and
`samples/ui/viewers/types/scatter-plot.js:12-15` (a follow-up `setOptions` after creation).

`getOptions(includeDefaults?)` (`viewer.ts:158-160`) serializes to
`{id, type, look: {...}}` — `look` carries the user-set options. `includeDefaults: true`
emits every property with its current value; that's the canonical "snapshot for layout
persistence" path.

`viewer.props` is an `ObjectPropertyBag` proxy assigned in `viewer.ts:118`:

```ts
this.props = new ObjectPropertyBag(this, api.grok_Viewer_Get_Look(this.dart));
```

It exposes get/set per property and also a `getProperty(name)` helper
(`samples/ui/viewers/inspect-viewer-properties.js:8-10`):

```ts
sp.props.xColumnName = 'race';
let descriptions = sp.props.getProperties()
  .map(p => p.propertyType + ' ' + p.name + ': ' + p.description + ' ' + p.columnFilter)
  .join('<br>');
```

`getProperties()` on `Viewer` (`viewer.ts:166-168`) returns the same `Property[]` array used
elsewhere — each `Property` has `name`, `propertyType`, `semType`, `description`, `category`,
`columnFilter`, `choices` (for enums). The existing `applyViewerOptions` helper
(`grokky/src/claude/grokky-api.ts:135-160`) uses `getProperties().map(p => p.name)` to build
a known-set, then attempts three reconciliations per supplied key:

1. Direct hit: `xColumnName` → `xColumnName`.
2. `+ColumnName` shortcut: `x` → `xColumnName`, `color` → `colorColumnName`.
3. `+ColumnNames` plural shortcut: `categoryColumns` → `categoryColumnNames`.
4. Case fix via lowercase map.
5. Unmatched → did-you-mean (`closestMatch`, Levenshtein, ≤3) → `console.warn`.

This is well-shaped but has gaps:

- **No type coercion.** `viewer.setOptions({bins: '30'})` will likely fail or no-op; the
  property's `propertyType` (`'int'`, `'bool'`, etc., from `Property`) is known but
  unused. Probably fine — Claude generates correct types from the schema — but worth a note.
- **The shortcut for `value` → `valueColumnName`** is *not* in the lookup. Histogram's main
  property is `valueColumnName` but the existing skill teaches users to write `value`. With
  the current shortcut logic, `value` → check `valueColumnName` (✓), so this already works.
  Verified at `grokky-api.ts:143`.
- **`split`** → `splitColumnName` works the same way. Good.
- **Trellis's `xColumnNames` is a list**; the shortcut at `:145` covers `xColumns` →
  `xColumnNames` but **not** `x` → `xColumnNames`. Claude tends to write the singular
  `x` form for trellis too — minor footgun, won't affect demos because the existing
  examples (`samples/ui/viewers/types/trellis-plot.js`) use the full plural.

### 4.1 What `getProperties` actually returns

Concretely, every entry is a `Property` object with at minimum: `name` (canonical, like
`'xColumnName'`), `propertyType` (string: `'string'`, `'int'`, `'double'`, `'bool'`,
`'column'`, `'column_list'`, `'color'`, etc.), `description`, optionally `choices`,
`columnFilter` (semantic constraint like `'numerical'`), and `defaultValue`. The existing
`inspect-viewer-properties.js:8-10` sample iterates them. The wrapper's reliance on
`p.name` only (`grokky-api.ts:135`) is robust as long as the JS-side property bag returns
*all* defined props — which it does, including style/look props (`backColor`, `axisFont`,
etc.).

### 4.2 Tags-based defaults

Per `samples/ui/viewers/advanced/default-properties.js:12-13`, df tags like
`'.Viewer Template: Scatter plot'` carry per-df default options. Useful note for the skill
but rarely needed in demos.

`Property.setDefault(data?, style?)` and `Property.resetDefault()` on the property bag
(`samples/ui/viewers/style-settings.js:16, 22`) persist the *current* state as the global
default for that viewer type. Out of scope.

## 5. Viewer-to-viewer linking (Area D)

There is **no explicit "link these two viewers"** API. Viewers sharing the same `DataFrame`
auto-sync because all interaction state lives on the dataframe:

- `df.selection` (BitSet, see topic 3)
- `df.filter` (BitSet, see topic 2)
- `df.currentRowIdx`
- `df.mouseOverRowIdx`

When a viewer fires `df.selection.fireChanged()`, every viewer bound to the same df receives
the event. Confirmed implicitly by `samples/ui/viewers/filter-select-highlight.js:1-18` —
three viewers (`scatterPlot`, `histogram`, `barChart`) all driven by `t.rows.filter(...)` and
`t.rows.select(...)` without any cross-wiring.

If two viewers must be on **different** dataframes, `grok.data.linkTables(t1, t2, keys1,
keys2, [SYNC_TYPE.SELECTION_TO_SELECTION, ...])` wires master-detail behavior — out of scope.

`viewer.dataFrame` (`viewer.ts:208-209`) is set automatically when a viewer is built via
`view.addViewer(type)` (it adopts `view.dataFrame`) or via `Viewer.fromType(type, df)` (it
uses `df`). The setter exists, but reassigning a mounted viewer's dataframe is rare and not
documented in samples.

## 6. Finding & closing viewers (Area E)

### 6.1 `view.viewers` — the iterator

`view.ts:625-627`:

```ts
get viewers(): Viewer[] {
  return api.grok_View_Get_Viewers(this.dart);
}
```

JSDoc warns: "Resulting array is cloned, so cache the result and do not call repeatedly in
inner loops." For one-shot exec, that doesn't matter. The grid is always at index 0 — verified
by `UITests/src/viewers/viewers.ts:251-254`:

```ts
function closeViewers(view: DG.TableView) {
  let viewers = Array.from(view.viewers).slice(1);  // ← skip [0], the grid
  viewers.forEach((v) => v?.close());
}
```

There is **no `getViewersByType(type)` method**. The idiom is:

```ts
Array.from(view.viewers).filter(v => v.type === DG.VIEWER.SCATTER_PLOT)
```

This is exactly the gap the new skill should close with a `findViewers(view, predicate)`
helper.

### 6.2 `viewer.close()` vs `viewer.removeFromView()` vs `view.detachViewers()`

`viewer.close()` (`viewer.ts:170-173`):

```ts
close(): void {
  api.grok_Viewer_Close(this.dart);
}
```

Closes *and* detaches. Calls cleanup logic — including the `onClick='Filter'` filter-reset
behavior covered by `ApiTests/src/ai/reported-issues/grok-19970.ts:15-77`. This is the
**preferred** method.

`viewer.removeFromView()` (`viewer.ts:355-357`):

```ts
removeFromView() {
  return toJs(api.grok_Viewer_Remove_From_View(this.dart));
}
```

Removes from the view's dock layout without disposing the viewer instance. Used in
`UITests/src/viewers/viewers.ts:76, 90` and in PhyloTreeViewer (`PhyloTreeViewer/src/viewers/grid-with-tree-viewer.ts:110, 116`). For our purposes, `close()` is what users mean when
they say "close the scatter plot".

`viewer.detach()` (only on `JsViewer`, `viewer.ts:429-432`) tears down rxjs subscriptions
in `this.subs`. It's the override hook for custom viewers — not what we want.

`view.detachViewers()` (`view.ts:619-621`) detaches **all** viewers on the view. Used by
TableView teardown. Together with `view.resetLayout()` (`view.ts:609-611`) — "resets view
layout, leaving only grid visible" — this is the heavy hammer. `resetLayout()` is probably
what the user expects for "close all the plots".

### 6.3 Edge cases

- `viewer.close()` on a viewer that was created but never attached (e.g. via `DG.Viewer.fromType`
  then never `view.addViewer(v)`) **throws**. Test note at `ApiTests/src/ai/reported-issues/grok-19970.ts:15`:
  "viewer.close() on a never-attached viewer throws — every close() call is guarded by a prior
  tv.addViewer(...)". The helper should swallow that case.
- `viewer.close()` is fire-and-forget; `tv.viewers` updates synchronously after the call
  per `UITests/src/viewers/viewers.ts:99-105` (`expect(Array.from(tv.viewers).length, 3); closeViewers(tv); expect(Array.from(tv.viewers).length, 1);`).

## 7. Viewer events (Area F)

Brief — out of scope for the per-exec-block skill but worth a single bullet:

- `grok.events.onViewerAdded` / `onViewerClosed` (`js-api/events.ts:202-204`) — global events.
- `viewer.onContextMenu` (`viewer.ts:308-310`), `onTooltipCreated` (`:106`), `onDataSelected`
  (`:107`), `onDataRowClicked` (`:109`), `onPropertyValueChanged` (`:110`).
- Per-subclass: `ScatterPlotViewer.onPointClicked` (`viewer.ts:667`), `onZoomed` (`:662`);
  `BarChartViewer.onCategoryClicked` (`:727`); `LineChartViewer.onLineSelected` (`:599`); etc.

## 8. Saving viewer state / layouts (Area G)

Per-viewer: `viewer.getOptions(true)` (`viewer.ts:158`) returns
`{id, type, look: {...allProps}}` — a serializable snapshot. Round-trip via
`viewer.setOptions(snapshot.look)`.

Per-view: `view.saveLayout()` (`js-api/views/view.ts:304-306`) returns a `ViewLayout`; pair
with `view.loadLayout(layout)` (`view.ts:296-298`). `view.saveState()` / `loadState(s, opts)`
(`view.ts:632-637`) are the string-flavored cousins. Covered fully in the future "layout"
skill — out of scope here. Cross-link briefly.

## 9. Dock layout (Area H)

`view.dockManager` (`view.ts:439-441`) returns the `DockManager` that owns viewer placement.
The signature is awkward — `DockManager.dock(element, dockType, refNode?, title?, ratio?)`
at `js-api/docking.ts:135`:

```ts
dock(element: HTMLElement | Viewer,
     dockType: DockType = DG.DOCK_TYPE.LEFT,
     refNode: DockNode | null = null,
     title?: string,
     ratio: number = 0.5): DockNode
```

`DOCK_TYPE` (`const.ts:776-782`) is `'left' | 'right' | 'up' | 'down' | 'fill'`. Strings work
directly — see `samples/ui/docking/docking.js:4-13`:

```js
let node1 = grok.shell.dockManager.dock(ui.divText('first'), 'right', null, 'First');
let node3 = grok.shell.dockManager.dock(ui.divText('third'), 'fill', node2, 'Third');
```

When `view.addViewer(...)` is called, the platform picks a default dock position
(`grok_TableView_ProcessNewViewer` — `view.ts:424`); to override, dock explicitly *after*
adding (`samples/ui/docking/docking-table-view.js:8-11`):

```js
let viewer = DG.Viewer.fromType('Scatter Plot', table);
view.addViewer(viewer);
view.dockManager.dock(viewer, 'right');
```

`view.dockManager.close(el)` (`docking.ts:143-150`) closes a docked element; for viewers,
prefer `viewer.close()` instead — it cleans up viewer state too.

The "sketcher on left, ADME on right" demo from the audit needs:
1. `view.addViewer('Sketcher', {...})` (plugin viewer, see §11)
2. `view.dockManager.dock(sketcher.root, 'left', null, 'Structure', 0.3)`
3. Add right-side viewers normally.

Touched lightly here; deferred to a dedicated layout skill.

## 10. Filter viewer (Area I)

Cross-reference to the **filtering** skill. The filter viewer is special: it's a `Viewer`
subclass (`FilterGroup`, `viewer.ts:510-553`) accessed via `view.getFiltersGroup({createDefaultFilters?})` (`view.ts:404-406`), and its
`updateOrAdd(state)` is the canonical filter-mutation API
(`grokky-api.ts:22-31` for `filterRows`).

If a user asks to "add filters", route them to `view.addViewer(DG.VIEWER.FILTERS)`
(samples confirm this works — `samples/ui/viewers/types/filters.js:5`). If they ask to
"filter on X", route to the filtering skill.

## 11. Existing `grokky.addViewer` / `configureViewer` (Area J)

Full source: `grokky/src/claude/grokky-api.ts:98-180`.

```ts
const VIEWER_TYPES: string[] = Object.values(DG.VIEWER);   // :98

function canonicalizeViewerType(input: string): string {    // :127
  const norm = input.toLowerCase().trim().replace(/\s+/g, ' ');
  const exact = VIEWER_TYPES.find(t => t.toLowerCase() === norm);
  if (exact) return exact;
  return closestMatch(VIEWER_TYPES, input, 3) ?? input;
}

function applyViewerOptions(viewer: DG.Viewer, options: Record<string, any>): void {  // :134
  const knownNames = viewer.getProperties().map((p: any) => p.name);
  // ... three-tier matching (exact / +ColumnName / +ColumnNames / case-fix) ...
  // Unknown keys → console.warn with did-you-mean hints
}

function addViewer(type: string, options: Record<string, any> = {}): DG.Viewer {  // :162
  const tv = grok.shell.tv;            // BUG: ignores in-scope `view`
  if (!tv) throw new Error('addViewer: no active TableView');
  const canonicalType = canonicalizeViewerType(type);
  const v = tv.addViewer(canonicalType);
  applyViewerOptions(v, options);
  return v;
}

function configureViewer(viewer: DG.Viewer, options: Record<string, any>): void {  // :171
  applyViewerOptions(viewer, options);
}
```

### What works well

- Fuzzy type matching via Levenshtein ≤3 (`:131`) — `'scatter'` → `'Scatter plot'`.
- Schema-derived property whitelist (`:135`) — no silent typos.
- Column shortcut (`x` → `xColumnName`) (`:143`).
- Plural list shortcut (`categoryColumns` → `categoryColumnNames`) (`:145`).
- Case-insensitive fallback (`:147`).
- Did-you-mean warnings, never throws on bad property (`:154-159`).

### Bugs and gaps to fix in the new skill

1. **`grok.shell.tv` instead of in-scope `view`** (`:163`). When `view` is in scope and is a
   TableView, the helper should use it. Falling back to `grok.shell.tv` is fine if no `view`
   is passed, but should warn.
2. **Plugin viewers fall through `canonicalizeViewerType`** (`:131`) — `'Chem space'` isn't
   in `DG.VIEWER` so `closestMatch` returns null and the raw input is passed. If the user
   means `Chem:chemSpaceTopMenu` (a function call, not a viewer type), the wrapper silently
   tries `addViewer('Chem space')`, fails to find a registered viewer of that name, and the
   user sees nothing. Better behavior: detect non-DG.VIEWER inputs and **either** call
   `Viewer.getViewerTypes()` for the plugin-loaded list, or refuse with a helpful message
   pointing to `callRegisteredFunction`.
3. **No `closeViewer` / `findViewer` helpers**.
4. **No "close all of type X"** — currently users have to write
   `Array.from(view.viewers).filter(v => v.type === '...').forEach(v => v.close())`.
5. **No graceful handling of `canVisualize` failure** — if you try `Histogram` on a df with
   no numerical columns, the platform might log a warning and produce an empty viewer. The
   wrapper could call `Viewer.canVisualize(type, df)` (`viewer.ts:359`) up-front.
6. **`type` is a hint, not authoritative.** `view.addViewer('Scatter Plot')` (note casing)
   works — the Dart side matches case-insensitively. The fuzzy `canonicalizeViewerType`
   already handles this, but the redundancy is fine.

## 12. The "close scatter, add histogram" idiom (Area K)

Audit demo: *"Add a histogram, close the scatter plot."* Canonical pattern with a helper:

```ts
grokky.closeViewers(view, {type: 'Scatter plot'});
grokky.addViewer(view, 'Histogram', {value: 'activities'});
```

Without a helper:

```ts
Array.from(view.viewers)
  .filter(v => v.type === DG.VIEWER.SCATTER_PLOT)
  .forEach(v => v.close());
view.addViewer('Histogram', {value: 'activities'});
```

The skill MUST show both forms — Claude needs to know the long form when no helper exists,
and the short form when one does.

## 13. Plugin viewers — chem-space and friends (Area L)

These are **not** viewer types. They are functions that *return* a viewer (already added to
the active TableView). Concrete example: `Chem:chemSpaceTopMenu`
(`chem/package.ts:821-877`):

```ts
static async chemSpaceTopMenu(table: DG.DataFrame, molecules: DG.Column, methodName: ...,
    similarityMetric: ..., plotEmbeddings: boolean, options?: ..., preprocessingFunction?: DG.Func,
    clusterEmbeddings?: boolean, clusterMCS?: boolean): Promise<DG.Viewer | undefined>
```

It runs a dimensionality reduction, adds embedding columns to the table, and (if
`plotEmbeddings`) calls `tv.scatterPlot({x: emb0, y: emb1, title: 'Chemical space'})`
(`chem/package.ts:858`) — returning that scatter plot. The user-visible result is a
scatter plot whose axes are UMAP/t-SNE coordinates.

Invocation from a `datagrok-exec` block:

```ts
const sp = await grok.functions.call('Chem:chemSpaceTopMenu', {
  table: t, molecules: t.col('smiles'),
  methodName: 'UMAP', similarityMetric: 'Tanimoto', plotEmbeddings: true,
});
```

Other plugin-provided viewer-like functions: `Chem:activityCliffs`, `Bio:sequenceSpace`,
`Bio:macromoleculeDifference`. They follow the same pattern: take a df + column + options,
mutate the table with extra columns, return a viewer.

The existing skill (`skill/SKILL.md:63`) sends these to a separate
`callRegisteredFunction` skill. That's the right call — but the new skill should still
explicitly **enumerate** the most common ones with the dispatch pattern, so Claude doesn't
try to fuzzy-match `'Chem space'` against `DG.VIEWER` and silently fail.

## 14. Proposed `grokky` helpers

### Signatures

```ts
// FIXED: accepts view-in-scope (TableView), falls back to grok.shell.tv with warning
addViewer(
  viewArg: DG.TableView | DG.ViewBase | null,
  type: string,
  options?: Record<string, any>,
): DG.Viewer;

// Unchanged from existing: keeps schema-validated property apply
configureViewer(viewer: DG.Viewer, options: Record<string, any>): void;

// NEW: close one or many. Three input shapes.
closeViewer(target: DG.Viewer | string | ((v: DG.Viewer) => boolean), view?: DG.TableView): number;
// returns count closed; tolerates already-closed/never-attached viewers

// NEW: list all viewers matching predicate
findViewers(
  view: DG.TableView | null,
  predicate?: (v: DG.Viewer) => boolean,   // default = all-but-grid
): DG.Viewer[];

// NEW: first match (convenience for "the scatter plot")
findViewer(
  view: DG.TableView | null,
  predicate: string | ((v: DG.Viewer) => boolean),  // string = type filter
): DG.Viewer | null;

// NEW: scorched-earth reset for "close all the plots"
closeAllViewers(view: DG.TableView | null, opts?: {keepGrid?: boolean}): number;
// keepGrid defaults to true; equivalent to view.resetLayout() but reports count
```

### Rationale

- **`addViewer(view, type, opts)` argument ordering**. `view` comes first because callers
  who *don't* have a `view` (rare in exec blocks) can pass `null` and fall back. This
  preserves the existing one-line ergonomics (`grokky.addViewer(view, 'Scatter plot', {x: ..., y: ...})`)
  while fixing the silent-shell-tv bug.
- **`closeViewer(target, view?)` with three input shapes** — `DG.Viewer`, `string` (type filter),
  `predicate` — mirrors the `selectRows` polymorphism from topic 3 (`grokky-api.ts:38-48`).
  Returns count so callers can branch on "was anything closed?".
- **`findViewers` / `findViewer`** — pure lookup helpers, no state mutation. Returns
  `Array.from(view.viewers).filter(...)`. Default predicate skips the grid (`[0]`).
- **`closeAllViewers`** — wraps `view.resetLayout()` when `keepGrid: true` (the default,
  matching `view.ts:609`'s docstring "leaves only the grid visible"), and the manual loop
  when `keepGrid: false`.

### Should there be `addPluginViewer(funcName, args)`?

Tempting but probably not. Plugin viewers are sufficiently varied (some return viewers, some
return dataframes, some mutate state and return nothing) that a single wrapper would either
be too thin to be useful or too fat to be predictable. The `callRegisteredFunction` skill
(separate) handles this better — and we keep `addViewer` cleanly scoped to built-in viewers
known to `DG.VIEWER`.

## 15. Anti-patterns

1. **`new DG.Viewer(...)`** — never directly construct. Use `Viewer.fromType` or static
   factories. The constructor (`viewer.ts:86`) takes a Dart-side handle and is for internal
   use. No samples in the codebase use `new`.
2. **`grok.shell.addViewer(...)`** — there is no such method on the shell (verified —
   `grep` returned zero matches in `js-api/src/`). Use `view.addViewer(...)` or
   `tv.addViewer(...)`.
3. **`grok.shell.tv` when `view` is in scope** — the source bug in the current helper
   (`grokky-api.ts:163`). If the runtime has populated `view`, it's the source of truth.
4. **Hardcoded property names** without consulting `viewer.getProperties()` — silent
   no-ops. Always validate (the existing helper already does; preserve this).
5. **`view.scatterPlot(opts)`, `view.histogram(opts)`** etc. — deprecated since the
   per-type method comments at `view.ts:443-606` say so. Always use `view.addViewer(type, opts)`.
6. **Mixing creation contexts**: `Viewer.scatterPlot(t)` followed by `view.addViewer(v)`
   when `t !== view.dataFrame` — the viewer ignores `view.dataFrame`. If on a TableView,
   prefer `view.addViewer('Scatter plot', opts)`.
7. **`viewer.close()` on a never-attached viewer** — throws. Guard with a
   try/catch or only call on viewers obtained from `view.viewers`.
8. **`view.viewers[0].close()`** — closes the grid, which then can't be reopened cleanly.
   Always skip index 0 (`.slice(1)`) when iterating to close.
9. **Manual REST or `fetch` to mutate viewers** — pointless. Everything is local JS.

## 16. Open questions

1. **`Viewer.canVisualize` pre-flight check.** Worth adding to `addViewer`? Pro: better
   error message ("scatter plot needs at least 2 numerical columns; this df has 1"). Con:
   adds a Dart round-trip per add. Probably yes if it's synchronous (the signature is —
   `viewer.ts:359` is `string | null`).

2. **Property type coercion.** Should the helper coerce `bins: '30'` to `bins: 30` based on
   the `Property.propertyType`? Trust-the-schema principle says no — fix at the Claude
   layer, not the wrapper.

3. **Plugin viewer enumeration.** Should `addViewer` accept any string in
   `Viewer.getViewerTypes()` (which includes plugin types)? Or stick to `DG.VIEWER` only?
   The former is more flexible; the latter avoids the silent-fail-on-Chem-space trap. **Lean
   toward**: detect `Viewer.getViewerTypes().includes(type)` for the wider list, but emit a
   warning "this is a plugin viewer, may load asynchronously".

4. **Async vs sync viewer creation.** `df.plot.fromType` is async; `Viewer.fromType` is
   sync. Plugin viewers may need async to wait for plugin load. If the wrapper always uses
   `view.addViewer(type)` (which goes through `grok_View_AddViewerByName` —
   `view.ts:419`), do we get sync or async behavior? Empirically sync (the samples never
   `await`), but plugin viewers may show empty briefly. **Test before claiming.**

5. **Scope of viewer state save/restore.** Belongs in the layout skill (T1.5), not this
   one. But `viewer.getOptions(true)` could be exposed here as a debug aid ("show me the
   current scatter plot settings"). **Lean toward**: yes, expose as `grokky.getViewerOptions(viewer)`
   one-liner.

6. **`view` typed as `DG.ViewBase`.** When `view` is `View` (not `TableView`), do we throw
   or fall back to `grok.shell.tv`? **Lean toward**: fall back with a warning. The user is
   in an AI panel that doesn't know if they're on a table view.

7. **Should `closeViewer(scatterPlot)` close the viewer if scatterPlot lives on a
   *different* view than the current one?** The viewer object knows its own view via
   `viewer.view` (`viewer.ts:198`). Yes, the `viewer.close()` call works regardless. Keep
   it simple — close any viewer passed in.

## 17. References by file

- `js-api/viewer.ts:80-373` — `Viewer` class, `fromType`, static factories, `setOptions`,
  `getOptions`, `getProperties`, `close`, `CORE_VIEWER_TYPES`.
- `js-api/viewer.ts:510-553` — `FilterGroup` (filter viewer subclass).
- `js-api/viewer.ts:577-833` — typed subclasses (`ScatterPlotViewer`, `BarChartViewer`,
  `HistogramViewer`, `BoxPlot`, `LineChartViewer`, `PieChartViewer`, `CorrelationPlot`,
  `ConfusionMatrix`, `RocCurve`, `CalendarViewer`, `PivotViewer`, `TrellisPlotViewer`,
  `DensityPlotViewer`, `PcPlot`).
- `js-api/views/view.ts:380-644` — `TableView`: `addViewer`, `scatterPlot/.../wordCloud`
  shorthands (deprecated), `viewers` iterator, `getFiltersGroup`, `dockManager`,
  `resetLayout`, `detachViewers`, `saveState/loadState`.
- `js-api/docking.ts:102-173` — `DockManager`, `DockNode`, `DockContainer`.
- `js-api/const.ts:670-733` — `DG.VIEWER`, `DG.CORE_VIEWER` enums.
- `js-api/const.ts:776-782` — `DOCK_TYPE` enum.
- `js-api/events.ts:202-204` — `onViewerAdded`, `onViewerClosed`.
- `js-api/dataframe/data-frame.ts:551-572` — `DataFramePlotHelper`.
- `js-api/shell.ts:64-83` — `grok.shell.t`, `.v`, `.tv` accessors.
- `grokky/src/claude/grokky-api.ts:98-180` — existing `addViewer`, `configureViewer`,
  fuzzy match, property normalization.
- `samples/ui/viewers/create-viewers.js`, `create-viewers-dynamically.js`,
  `inspect-viewer-properties.js`, `list-viewers.js`, `viewer-info.js`,
  `scatter-plot-events.js`, `style-settings.js`, `filter-select-highlight.js`,
  `trellis-plot-with-summary-viewer.js`, `line-chart-multiple-series.js`,
  `markup-advanced.js`, `types/*.js` (every built-in viewer demo).
- `samples/ui/docking/docking-table-view.js:1-15` — `fromType` + `addViewer` +
  `dockManager.dock` flow.
- `samples/ui/viewers/advanced/default-properties.js:12-13`,
  `advanced/default-axis-type.js:2`, `advanced/attached-properties.js:1-19`,
  `advanced/init-function.js:9-20` — tag-based defaults, attached properties.
- `UITests/src/viewers/viewers.ts:99-105, 251-254` — close-all-but-grid idiom, expected
  viewer count.
- `ApiTests/src/ai/reported-issues/grok-19970.ts:15, 40-46` — `viewer.close()` semantics,
  never-attached throws, `onClick='Filter'` cleanup.
- `chem/package.ts:821-877` — `chemSpaceTopMenu` plugin viewer pattern.
- `skill/SKILL.md` — existing skill (to be replaced).
- `skill/viewer-properties.md` — existing per-viewer property reference (cite shape; extend).
