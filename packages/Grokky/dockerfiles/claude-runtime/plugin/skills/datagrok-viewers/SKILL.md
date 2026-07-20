---
name: datagrok-viewers
description: Add a viewer, configure a viewer, change viewer options, find viewer, close viewer, view a scatter plot, bar chart, histogram, line chart, box plot, pie chart, heat map, correlation plot, 3D scatter, trellis, density plot, statistics, on a Datagrok TableView inside a datagrok-exec block. Use whenever the user asks to plot, chart, visualize, show a graph, draw a distribution, color by a column, swap a viewer's axis, toggle a legend / regression line / log scale, replace one viewer with another, close every chart, reset the view to just the grid, or find an existing viewer by type. Plugin viewers like "Chem space", "sequence space", "activity cliffs" are NOT viewer types — they're registered functions — route those to `grok.functions.call`. Does NOT cover filtering (separate skill `datagrok-filtering`), selection (`datagrok-selection`), grid cell rendering (`datagrok-grid-customization`), layout save/restore, or custom-viewer authoring.
---

# datagrok-viewers

Add, configure, find, and close viewers on a `DG.TableView` from inside a
`datagrok-exec` block. Globals injected by the runtime: `grok`, `ui`, `DG`,
`view`, `t` (the current `DG.DataFrame`, when the view is a `TableView`).

## Quick reference

| Need                                 | Call                                                                  |
|--------------------------------------|-----------------------------------------------------------------------|
| Add a viewer                         | `view.addViewer(DG.VIEWER.SCATTER_PLOT, {xColumnName: 'a', ...})`     |
| Configure an existing viewer         | `viewer.setOptions({colorColumnName: 'logP'})`                         |
| All non-grid viewers                 | `Array.from(view.viewers).slice(1)`                                    |
| First viewer of a type               | `Array.from(view.viewers).slice(1).find((v) => v.type === DG.VIEWER.SCATTER_PLOT)` |
| Close one viewer                     | `viewer.close()` (wrap in try/catch — throws if never attached)        |
| Close every non-grid viewer          | `Array.from(view.viewers).slice(1).forEach((v) => { try { v.close(); } catch {} })` |
| Inspect a viewer's option schema     | `viewer.getProperties()` → `DG.Property[]` (each has `.name`)          |
| Standalone (detached) viewer         | `DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t, {...})`                 |

## THE viewer footgun: column-binding properties always end in `ColumnName(s)`

Users say "x is mw, y is logp". You must emit `xColumnName` / `yColumnName` —
never bare `x` / `y` / `color`. Bare aliases are inconsistent across viewer
types and may silently do nothing for the canonical `addViewer`/`setOptions`
path. Always use the full `*ColumnName(s)` form.

| User says           | Emit (single column)                | Emit (multiple columns)         |
|---------------------|-------------------------------------|---------------------------------|
| `x is mw`           | `xColumnName: 'mw'`                 | `xColumnNames: ['mw', ...]`     |
| `y is logp`         | `yColumnName: 'logp'`               | `yColumnNames: ['logp', ...]`   |
| `color by activity` | `colorColumnName: 'activity'`       | —                               |
| `size by mw`        | `sizeColumnName: 'mw'`              | —                               |
| `label by smiles`   | `labelColumnName: 'smiles'`         | —                               |
| `split by class`    | `splitColumnName: 'class'`          | `splitColumnNames: ['class']`   |
| `value is ic50`     | `valueColumnName: 'ic50'`           | `valueColumnNames: ['ic50']`    |
| `category is class` | `categoryColumnName: 'class'`       | `categoryColumnNames: ['class']`|
| `columns are a, b`  | —                                   | `columnNames: ['a', 'b']`       |

- **Single column** properties: singular `<role>ColumnName` (string).
- **Multi-column** properties: plural `<role>ColumnNames` (array of strings).
- Some viewers take both forms for different roles (Line chart has
  `xColumnName` for the x axis and `yColumnNames` for the series).
- When in doubt, inspect the live schema: `viewer.getProperties()` returns
  `DG.Property[]`. Each has `.name` (canonical), `.propertyType`, `.choices`.

Other naming conventions across viewers:

- **Boolean toggles start with `show`**: `showRegressionLine`, `showXAxis`,
  `showLegend`, `showStatistics`, `showColorSelector`. Not `regressionLine`
  / `legend` / `statistics`.
- **Axis scale**: `xAxisType: 'linear' | 'logarithmic'` (and `y`, `z`,
  `colorAxisType`). Not `logX`, not `logScale`.
- **Axis bounds**: `xMin`, `xMax`, `yMin`, `yMax`. Not `xRange: [a, b]`.

## The viewers array

`view.viewers` is a synchronous array snapshot:

- `view.viewers[0]` is **always the grid**. Closing it breaks the table widget.
  Every iteration must `.slice(1)` (or filter out `v.type === DG.VIEWER.GRID`).
- Every other entry is a non-grid viewer. Each has `.type` (a string equal to
  one of the `DG.VIEWER.*` constants) and `.dataFrame`.

## Viewer types catalog

Always use the `DG.VIEWER.*` constants, never literal strings:

| Constant                       | String value         | Primary properties                                                                |
|--------------------------------|----------------------|-----------------------------------------------------------------------------------|
| `DG.VIEWER.SCATTER_PLOT`       | `'Scatter plot'`     | `xColumnName`, `yColumnName`, `colorColumnName`, `sizeColumnName`, `showRegressionLine` |
| `DG.VIEWER.HISTOGRAM`          | `'Histogram'`        | `valueColumnName`, `splitColumnName`, `bins`                                      |
| `DG.VIEWER.LINE_CHART`         | `'Line chart'`       | `xColumnName`, `yColumnNames` (array)                                             |
| `DG.VIEWER.BAR_CHART`          | `'Bar chart'`        | `splitColumnName`, `valueColumnName`, `valueAggrType`                             |
| `DG.VIEWER.PIE_CHART`          | `'Pie chart'`        | `categoryColumnName`                                                              |
| `DG.VIEWER.BOX_PLOT`           | `'Box plot'`         | `valueColumnName`, `categoryColumnName`                                           |
| `DG.VIEWER.HEAT_MAP`           | `'Heat map'`         | (returns a Grid subclass — column tag-driven color coding)                        |
| `DG.VIEWER.STATISTICS`         | `'Statistics'`       | `columnNames`                                                                     |
| `DG.VIEWER.CORR_PLOT`          | `'Correlation plot'` | `columnNames` (or `xs` / `ys`)                                                    |
| `DG.VIEWER.DENSITY_PLOT`       | `'Density plot'`     | `xColumnName`, `yColumnName`                                                      |
| `DG.VIEWER.SCATTER_PLOT_3D`    | `'3d scatter plot'`  | `xColumnName`, `yColumnName`, `zColumnName`, `colorColumnName`                    |
| `DG.VIEWER.TRELLIS_PLOT`       | `'Trellis plot'`     | `xColumnNames`, `yColumnNames`, `viewerType`                                      |
| `DG.VIEWER.PC_PLOT`            | `'PC Plot'`          | `columnNames`                                                                     |
| `DG.VIEWER.TREE_MAP`           | `'Tree map'`         | `splitColumnNames`, `sizeColumnName`                                              |
| `DG.VIEWER.MATRIX_PLOT`        | `'Matrix plot'`      | `xColumnNames`, `yColumnNames`                                                    |
| `DG.VIEWER.FILTERS`            | `'Filters'`          | `columnNames` (see `datagrok-filtering` for the real filtering API)               |
| `DG.VIEWER.FORM`               | `'Form'`             | `columnNames`                                                                     |
| `DG.VIEWER.MARKUP`             | `'Markup'`           | `content` (string)                                                                |
| `DG.VIEWER.NETWORK_DIAGRAM`    | `'Network diagram'`  | `node1ColumnName`, `node2ColumnName`                                              |
| `DG.VIEWER.WORD_CLOUD`         | `'Word cloud'`       | column                                                                            |
| `DG.VIEWER.TILE_VIEWER`        | `'Tile Viewer'`      | per-row cards                                                                     |
| `DG.VIEWER.PIVOT_TABLE`        | `'Pivot table'`      | grouped aggregates                                                                |
| `DG.VIEWER.SHAPE_MAP`          | `'Shape Map'`        | choropleth                                                                        |
| `DG.VIEWER.CALENDAR`           | `'Calendar'`         | `dateColumnName`                                                                  |
| `DG.VIEWER.GLOBE`              | `'Globe'`            | `latitudeColumnName`, `longitudeColumnName`                                       |
| `DG.VIEWER.GOOGLE_MAP`         | `'Google map'`       | `latitudeColumnName`, `longitudeColumnName`                                       |
| `DG.VIEWER.TIMELINES`          | `'Timelines'`        | `startColumnName`, `endColumnName`                                                |
| `DG.VIEWER.RADAR_VIEWER`       | `'Radar'`            | `columnNames`                                                                     |
| `DG.VIEWER.GRID`               | `'Grid'`             | (the main grid — `view.grid` is the typed accessor)                               |

Full canonical list at runtime: `Object.values(DG.VIEWER)`.

## Adding viewers

Use the in-scope `view` global, not `grok.shell.tv` (the *globally* active
view, which may be a different tab when the user has several open).

```datagrok-exec
// Scatter plot of MW vs LogP, colored by activity, with regression line.
view.addViewer(DG.VIEWER.SCATTER_PLOT, {
  xColumnName: 'MW',
  yColumnName: 'LogP',
  colorColumnName: 'activity',
  showRegressionLine: true,
});
```

```datagrok-exec
// Line chart with multiple y series — yColumnNames is the plural form (array).
view.addViewer(DG.VIEWER.LINE_CHART, {
  xColumnName: 'date',
  yColumnNames: ['revenue', 'expenses', 'profit'],
});
```

```datagrok-exec
// Bar chart of count by category.
view.addViewer(DG.VIEWER.BAR_CHART, {
  splitColumnName: 'category',
  valueColumnName: 'count',
  valueAggrType: 'sum',
});
```

## Configuring an existing viewer

`viewer.setOptions({...})` applies a batch of property changes and re-renders
once. Same property-name rules as `addViewer`.

```datagrok-exec
// Recolor the existing scatter plot by a different column.
const sp = Array.from(view.viewers).slice(1)
  .find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
if (sp)
  sp.setOptions({colorColumnName: 'logP'});
```

`view.addViewer(...)` returns the freshly attached viewer, so you can chain
create + configure:

```datagrok-exec
const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT, {
  xColumnName: 'height',
  yColumnName: 'weight',
  showRegressionLine: true,
});
sp.setOptions({xMin: 150, xMax: 200, colorColumnName: 'age'});
```

To discover property names, ask the viewer:

```datagrok-exec
const v = view.addViewer(DG.VIEWER.SCATTER_PLOT);
const names = v.getProperties().map((p) => p.name);
return ui.divText(names.join(', '));
```

## Finding viewers

`view.viewers` is an array; **index 0 is always the grid**, so iteration always
starts with `.slice(1)`.

```datagrok-exec
// First scatter plot on the view (or undefined).
const sp = Array.from(view.viewers).slice(1)
  .find((v) => v.type === DG.VIEWER.SCATTER_PLOT);
```

```datagrok-exec
// All histograms.
const hs = Array.from(view.viewers).slice(1)
  .filter((v) => v.type === DG.VIEWER.HISTOGRAM);
```

## Closing viewers

`viewer.close()` closes and detaches the viewer. It **throws** if the viewer
was never attached to a view, so wrap in try/catch.

```datagrok-exec
// Close every scatter plot on the view.
Array.from(view.viewers).slice(1)
  .filter((v) => v.type === DG.VIEWER.SCATTER_PLOT)
  .forEach((v) => { try { v.close(); } catch {} });
```

```datagrok-exec
// Close every non-grid viewer (reset the view to just the grid).
Array.from(view.viewers).slice(1)
  .forEach((v) => { try { v.close(); } catch {} });
```

Close-and-replace pattern for "show me X instead of Y":

```datagrok-exec
// Close all scatter plots, then add a histogram.
Array.from(view.viewers).slice(1)
  .filter((v) => v.type === DG.VIEWER.SCATTER_PLOT)
  .forEach((v) => { try { v.close(); } catch {} });
view.addViewer(DG.VIEWER.HISTOGRAM, {valueColumnName: 'activity'});
```

## Plugin viewers — they are functions, not types

Some "viewers" users name are actually package-registered functions that
*create* a viewer as a side effect. They aren't in `DG.VIEWER`. Route via
`grok.functions.call`:

| User phrase                  | Function call                                                                                                                  |
|------------------------------|--------------------------------------------------------------------------------------------------------------------------------|
| "chem space" / "chemical space" | `grok.functions.call('Chem:chemSpaceTopMenu', {table: t, molecules: t.col('smiles'), methodName: 'UMAP', similarityMetric: 'Tanimoto', plotEmbeddings: true})` |
| "sequence space"             | `grok.functions.call('Bio:sequenceSpaceTopMenu', {table: t, molecules: t.col('sequence'), methodName: 'UMAP', similarityMetric: 'Levenshtein', plotEmbeddings: true})` |
| "activity cliffs"            | `grok.functions.call('Chem:activityCliffs', {table: t, molecules: t.col('smiles'), activities: t.col('activity'), similarity: 80, methodName: 'UMAP'})` |
| "elemental analysis"         | `grok.functions.call('Chem:elementalAnalysis', {table: t, molecules: t.col('smiles')})`                                       |

Each mutates the DataFrame and attaches a viewer to the current TableView.
Only the *Space functions return the freshly-attached viewer (`Chem:chemSpaceTopMenu`,
`Bio:sequenceSpaceTopMenu`). `Chem:activityCliffs` and `Chem:elementalAnalysis`
return `void`, so to keep configuring you must find the new viewer via
`view.viewers`:

```datagrok-exec
const sp = await grok.functions.call('Chem:chemSpaceTopMenu', {
  table: t, molecules: t.col('smiles'),
  methodName: 'UMAP', similarityMetric: 'Tanimoto', plotEmbeddings: true,
});
if (sp)
  sp.setOptions({colorColumnName: 'activity'});
```

## Standalone (detached) viewers

Use `DG.Viewer.fromType` **only** when embedding a viewer in a custom UI
element you return from the block. The caller is responsible for placing the
result. For everything else, prefer `view.addViewer`.

```datagrok-exec
// Build a scatter plot of t and embed its DOM root in the chat.
const sp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, t, {
  xColumnName: 'mw', yColumnName: 'logp',
});
return sp.root;
```

## Anti-patterns

1. **Bare `x` / `y` / `color` in options** — `{x: 'mw', y: 'logp'}` silently
   does nothing. The property is `xColumnName`.
2. **Literal strings instead of `DG.VIEWER.*`** — use `DG.VIEWER.SCATTER_PLOT`,
   not `'scatter plot'`.
3. **`view.scatterPlot(opts)` / `view.histogram(opts)` / `view.barChart(opts)`** —
   per-type shorthands are deprecated. Use `view.addViewer(DG.VIEWER.*, opts)`.
4. **`new DG.Viewer(...)`** — internal constructor over a Dart handle. Use
   `view.addViewer(type, opts)` or `DG.Viewer.fromType(type, df, opts?)`.
5. **`grok.shell.addViewer(...)`** — no such method. Use `view.addViewer(...)`.
6. **`grok.shell.tv` when `view` is in scope** — `view` is the in-scope
   TableView. `grok.shell.tv` is the globally active view and may be a
   different tab.
7. **Treating "chem space" / "sequence space" / "activity cliffs" as viewer
   types** — they're functions called via `grok.functions.call`.
8. **Iterating `view.viewers` without `.slice(1)`** — index 0 is the grid.
   `view.viewers.forEach(v => v.close())` closes the grid.
9. **`viewer.close()` on a never-attached viewer** — throws. Guard with
   try/catch.
10. **`await view.addViewer(...)`** — synchronous for built-in types. Only
    plugin viewers via `grok.functions.call(...)` are async.

## Out of scope

- **Row filtering** (`df.filter`) — `datagrok-filtering`.
- **Grid customization** (visibility, widths, pinning) —
  `datagrok-grid-customization`; column color coding —
  `datagrok-df-and-columns`.
- **The Filters viewer (`DG.VIEWER.FILTERS`)** — `view.getFiltersGroup()`
  is the right API and is covered by `datagrok-filtering`. Only add a
  Filters viewer here if the user explicitly says "add a filter panel as a
  viewer".
