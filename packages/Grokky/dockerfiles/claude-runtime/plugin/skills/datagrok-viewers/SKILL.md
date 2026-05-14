---
name: datagrok-viewers
description: Add, configure, find, and close viewers (scatter plot, histogram, line/bar chart, box plot, pie chart, trellis, heat map, correlation plot, 3D scatter, statistics, density, etc.) on a Datagrok TableView inside a datagrok-exec block. Use whenever the user asks to plot, chart, visualize, show a graph, draw a distribution, color by a column, swap a viewer's axis, toggle a legend / regression line / log scale, replace one viewer with another, close every chart, reset the view to just the grid, or find an existing viewer by type. Plugin viewers like "Chem space", "sequence space", "activity cliffs" are NOT viewer types — they're registered functions — and this skill routes those to `grok.functions.call`. Does NOT cover filtering (separate skill `datagrok-filtering`), selection (`datagrok-selection`), grid cell rendering, layout save/restore, or custom-viewer authoring.
---

# datagrok-viewers

Use the `grokky.*` viewer helpers inside a `datagrok-exec` block. They wrap
`TableView.addViewer`, the per-viewer `Property[]` schema (`getProperties()`),
and the `view.viewers` iterator so Claude doesn't have to remember that
`new DG.Viewer(...)` is internal-only, that the typed `view.scatterPlot(opts)`
shorthands are deprecated, that `view.viewers[0]` is always the grid (and
closing it breaks the view), that `viewer.close()` throws on a never-attached
viewer, or that "Chem space" is a function call — not a viewer type.

## What this skill covers

Creating and configuring built-in viewers on a `DG.TableView`: pick the type
(fuzzy-matched against `DG.VIEWER`), set column-binding properties
(`xColumnName`, `yColumnName`, `colorColumnName`, ...), tweak look (axis
log scale, regression line, legend visibility), find existing viewers
(predicate-based, type-based), close one viewer / all viewers of a type /
every viewer (keep grid), and the canonical "close X then add Y" replace
pattern.

**Out of scope.** Row filtering (`df.filter`) lives in `datagrok-filtering`.
Row selection (`df.selection`) lives in `datagrok-selection`. Custom grid
cell rendering, pinned columns, and grid-specific UI are in the (future)
`datagrok-grid` skill. Saving / restoring viewer state across sessions
(`view.saveLayout()`, `viewer.getOptions(true)`) belongs to the layout skill.
Authoring a custom `JsViewer` subclass is the package-development path, not
exec-block territory. The Filters viewer is handled by `datagrok-filtering`
(`view.getFiltersGroup`) — only add `DG.VIEWER.FILTERS` from here if the
user explicitly wants a filter panel as a viewer.

## Quick reference

| Helper                                          | One-liner                                                                          |
|-------------------------------------------------|------------------------------------------------------------------------------------|
| `grokky.addViewer(view, type, opts?)`           | Attach a built-in viewer. Fuzzy type, schema-validated options, did-you-mean hints.|
| `grokky.configureViewer(viewer, opts)`          | Set/update properties on an existing viewer. Same schema validation.               |
| `grokky.findViewer(view, pred)`                 | First non-grid viewer matching predicate, or `null`. String pred = type filter.    |
| `grokky.findViewers(view, pred?)`               | All non-grid viewers matching predicate. No pred = every non-grid viewer.          |
| `grokky.closeViewer(target, view?)`             | Close one viewer, all viewers of a type, or all matching a predicate. Tolerant.    |
| `grokky.closeAllViewers(view, opts?)`           | Scorched earth. `{keepGrid: true}` by default — view.viewers[0] (grid) is spared.  |

Globals available inside every `datagrok-exec` block: `grok`, `ui`, `DG`,
`view`, `t` (the current `DG.DataFrame`, when the view is a TableView),
`grokky`.

## The mental model

`view.viewers` is a synchronous array snapshot:

- `view.viewers[0]` is **always the grid**. Always. Closing it breaks the
  view's main table widget. Every iteration must `.slice(1)` (or use the
  helpers, which do that for you).
- Every other entry is a viewer the user (or previous code) added — scatter
  plot, histogram, filter panel, etc. Each has `.type` (one of the
  `DG.VIEWER.*` strings) and `.dataFrame` (typically `view.dataFrame`).

Three viewer creation paths exist in the JS API. **Only one is right for
exec blocks**:

| Path                                     | Verdict for this skill                                                |
|------------------------------------------|-----------------------------------------------------------------------|
| `view.addViewer('Scatter plot', opts)`   | **Canonical**. Creates + docks in one call. Use via `grokky.addViewer`.|
| `DG.Viewer.fromType(type, df, opts?)`    | Detached viewer; caller must dock. Off-view rendering only.            |
| `view.scatterPlot(opts)` / `histogram()` | **Deprecated**. Source comments say so. Never use.                     |
| `new DG.Viewer(...)`                     | Internal constructor over a Dart handle. Never use.                    |

Viewers cross-talk automatically when they share a `DataFrame`. Selection,
filter, current row, mouse-over — all live on the DF, and every viewer bound
to that DF reacts. There is no `linkViewers()` to call.

## Viewer types catalog

These are the top demo viewers and their primary properties. Full list lives
in `DG.VIEWER.*` (~30 entries). The wrapper fuzzy-matches type names with
Levenshtein distance ≤ 3, case-insensitive — `'scatter'`, `'Scatter Plot'`,
`'scatterplot'` all resolve to `'Scatter plot'`.

| Type string          | Use for                          | Primary axes / properties                                                              |
|----------------------|----------------------------------|----------------------------------------------------------------------------------------|
| `'Scatter plot'`     | XY relationship                  | `xColumnName`, `yColumnName`, `colorColumnName`, `sizeColumnName`, `showRegressionLine`|
| `'Histogram'`        | Distribution of a single column  | `valueColumnName`, `splitColumnName`, `bins`                                           |
| `'Line chart'`       | Trend / time-series              | `xColumnName`, `yColumnNames` (array)                                                  |
| `'Bar chart'`        | Categorical counts / aggregates  | `splitColumnName`, `valueColumnName`, `valueAggrType`                                  |
| `'Pie chart'`        | Categorical proportions          | `categoryColumnName`                                                                   |
| `'Box plot'`         | Numeric by category              | `valueColumnName`, `categoryColumnName`                                                |
| `'Heat map'`         | Cell-color matrix                | (returns a Grid subclass — unusual)                                                    |
| `'Statistics'`       | Per-column summary               | `columnNames`                                                                          |
| `'Correlation plot'` | Pairwise correlations            | `columnNames` (or `xs` / `ys`)                                                         |
| `'Density plot'`     | Smoothed 2D density              | `xColumnName`, `yColumnName`                                                           |
| `'3d scatter plot'`  | Three-axis numerical             | `xColumnName`, `yColumnName`, `zColumnName`, `colorColumnName`                         |
| `'Trellis plot'`     | Small multiples                  | `xColumnNames`, `yColumnNames`, `viewerType`                                           |
| `'PC Plot'`          | Parallel coordinates             | `columnNames`                                                                          |
| `'Tree map'`         | Hierarchical area                | `splitColumnNames`, `sizeColumnName`                                                   |
| `'Matrix plot'`      | Scatter matrix                   | `xColumnNames`, `yColumnNames`                                                         |
| `'Filters'`          | Filter panel as a viewer         | `columnNames` (see `datagrok-filtering` instead)                                       |
| `'Form'`             | Per-row record view              | `columnNames`                                                                          |
| `'Markup'`           | Static HTML/markdown             | `content` (string)                                                                     |

For the full canonical list at runtime: `Object.values(DG.VIEWER)` (every
key surfaced by `DG.VIEWER.*`). Common variants `'scatter plot'`,
`'scatter'`, `'Scatterplot'`, `'Scatter Plot'` all resolve through the
helper. The property table above is a starter set — `getProperties()` is
the authoritative schema and the wrapper validates every key against it,
emitting a did-you-mean warning on unknown keys (`'xColumn'` →
`'xColumnName'`).

### Boolean and axis property conventions

These conventions hold across nearly all viewers — use them as a writing
prior, then trust the schema-validation warnings to correct edge cases:

- **Column references** end in `ColumnName` or `ColumnNames`:
  `xColumnName`, `colorColumnName`, `categoryColumnNames`. The wrapper
  accepts the short forms (`x`, `color`, `categoryColumns`) and rewrites
  them. Never `xColumn` or `colorCol`.
- **Boolean toggles start with `show`**: `showRegressionLine`, `showXAxis`,
  `showLegend`, `showStatistics`, `showXSelector`, `showColorSelector`.
  Never `regressionLine` / `legend` / `statistics`.
- **Axis scale**: `xAxisType: 'linear' | 'logarithmic'` (and `y`, `z`,
  `colorAxisType`). Not `logX`, not `logScale`.
- **Axis bounds**: `xMin`, `xMax`, `yMin`, `yMax`. Not `xRange: [a, b]`.

## Plugin viewers — they are functions, not types

A handful of "viewers" users will name are actually package-registered
functions that *create* a viewer as a side effect. Don't try to fuzzy-match
them against `DG.VIEWER` — `grokky.addViewer` will fail to find them in the
canonical list. Route via `grok.functions.call(...)` instead:

| User phrase                              | Function call                                                                                                                       |
|------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------|
| "chemical space" / "chem space"          | `grok.functions.call('Chem:chemSpaceTopMenu', {table: t, molecules: t.col('smiles'), methodName: 'UMAP', similarityMetric: 'Tanimoto', plotEmbeddings: true})` |
| "sequence space"                         | `grok.functions.call('Bio:sequenceSpaceTopMenu', {table: t, molecules: t.col('sequence'), methodName: 'UMAP', similarityMetric: 'Levenshtein', plotEmbeddings: true})` |
| "activity cliffs"                        | `grok.functions.call('Chem:activityCliffs', {table: t, molecules: t.col('smiles'), activities: t.col('activity'), similarity: 80, methodName: 'UMAP', ...})` |
| "elemental analysis"                     | `grok.functions.call('Chem:elementalAnalysis', {table: t, molecules: t.col('smiles')})`                                            |

Each of these mutates the DataFrame (adds embedding columns, transformations)
and returns the freshly-attached viewer. Common pattern:

```datagrok-exec
const sp = await grok.functions.call('Chem:chemSpaceTopMenu', {
  table: t, molecules: t.col('smiles'),
  methodName: 'UMAP', similarityMetric: 'Tanimoto', plotEmbeddings: true,
});
```

You **can** further `grokky.configureViewer(sp, {...})` afterwards — the
returned viewer obeys the regular property schema. But don't try
`grokky.addViewer(view, 'Chem space', ...)` — there's no such viewer type.

## Adding viewers

```ts
grokky.addViewer(
  view: DG.TableView | DG.ViewBase | null,
  type: string,
  options?: Record<string, any>,
): DG.Viewer;
```

Pass `view` first. If `view` is a `TableView`, the viewer is attached there.
If `view` is `null` or a non-TableView (script view, function view, ...), the
helper falls back to `grok.shell.tv` with a `console.warn`. Throws only if
no active TableView exists anywhere.

The wrapper:

- Fuzzy-matches `type` against `DG.VIEWER` strings (Levenshtein ≤ 3,
  case-insensitive). `'scatter'` → `'Scatter plot'`.
- Reads the viewer's live property schema via `viewer.getProperties()`.
- Applies your `options` keys with three rewrites: direct hit (`xColumnName`),
  column-name suffix (`x` → `xColumnName`), plural list (`categoryColumns`
  → `categoryColumnNames`), then a case-insensitive fallback.
- Unknown keys → `console.warn` with a did-you-mean (`'regressionLine'` →
  `'showRegressionLine'?`). Never throws on typos. The valid subset is
  still applied; the unknown keys are dropped.

```datagrok-exec
// Scatter plot of height vs weight, colored by age, with regression line.
grokky.addViewer(view, 'Scatter plot', {
  xColumnName: 'height',
  yColumnName: 'weight',
  colorColumnName: 'age',
  showRegressionLine: true,
});
```

```datagrok-exec
// Histogram of activities per compound class — `split` is the column we
// overlay one bar stack per category on.
grokky.addViewer(view, 'Histogram', {
  valueColumnName: 'activity',
  splitColumnName: 'compound_class',
});
```

```datagrok-exec
// Line chart with multiple y series. `yColumnNames` is plural — an array.
grokky.addViewer(view, 'Line chart', {
  xColumnName: 'date',
  yColumnNames: ['revenue', 'expenses', 'profit'],
});
```

Short-form column properties also work — the wrapper rewrites to canonical:

```datagrok-exec
// Equivalent to the first example. `x` → `xColumnName`, `y` → `yColumnName`.
grokky.addViewer(view, 'Scatter plot', {x: 'height', y: 'weight', color: 'age', showRegressionLine: true});
```

### View scope

Old behavior used `grok.shell.tv` (the global active TableView), which broke
if the user was on a different view at the moment of execution. The fixed
helper prefers in-scope `view` and only falls back if it has to. Always
pass `view` when you have it.

```datagrok-exec
// Wrong: ignores `view`, hits the global. If user just switched tabs,
// the viewer lands on the wrong table.
const tv = grok.shell.tv;
const sp = tv.addViewer('Scatter plot');

// Right: uses the view you were given.
grokky.addViewer(view, 'Scatter plot', {x: 'height', y: 'weight'});
```

## Configuring an existing viewer

```ts
grokky.configureViewer(viewer: DG.Viewer, options: Record<string, any>): void;
```

Same schema validation and did-you-mean logic as `addViewer`. Use to tweak
already-attached viewers:

```datagrok-exec
// Find the scatter plot and recolor it by a different column.
const sp = grokky.findViewer(view, (v) => v.type === DG.VIEWER.SCATTER_PLOT);
if (sp)
  grokky.configureViewer(sp, {colorColumnName: 'logP'});
```

```datagrok-exec
// Flip a histogram's y-axis to logarithmic.
const h = grokky.findViewer(view, (v) => v.type === DG.VIEWER.HISTOGRAM);
if (h)
  grokky.configureViewer(h, {yAxisType: 'logarithmic'});
```

`configureViewer` calls `viewer.setOptions(...)` under the hood, which
notifies the viewer to re-render in one batch. Equivalent to setting each
property via `viewer.props.xxx = yyy` — preferred because it's a single
notification.

## Finding viewers

| Need                          | Helper                                                                |
|-------------------------------|-----------------------------------------------------------------------|
| The first scatter plot        | `findViewer(view, v => v.type === DG.VIEWER.SCATTER_PLOT)`            |
| All histograms                | `findViewers(view, v => v.type === DG.VIEWER.HISTOGRAM)`              |
| All non-grid viewers          | `findViewers(view)`                                                   |
| First match by type string    | `findViewer(view, 'Scatter plot')`  *(string pred → type filter)*    |
| Count of viewers in the view  | `findViewers(view).length`                                            |

Both helpers **skip `view.viewers[0]` (the grid)** — closing the grid is
almost certainly a bug, and "the viewers" in user speech almost always means
the non-grid ones. Need the grid? `view.grid` is the typed accessor.

```datagrok-exec
// Count non-grid viewers currently attached.
return {count: grokky.findViewers(view).length};
```

```datagrok-exec
// First scatter plot (or null) — same predicate, two phrasings.
const sp1 = grokky.findViewer(view, (v) => v.type === DG.VIEWER.SCATTER_PLOT);
const sp2 = grokky.findViewer(view, 'Scatter plot');
return {handle: sp1 === sp2};
```

## Closing viewers

```ts
grokky.closeViewer(
  target: DG.Viewer | string | ((v: DG.Viewer) => boolean),
  view?: DG.TableView,
): number;
```

Returns the count closed. Tolerates the "never-attached viewer throws on
`close()`" edge case (try/catch internally — failures are logged, not
propagated). Polymorphic by target type:

```datagrok-exec
// Close a specific viewer handle.
const sp = grokky.findViewer(view, 'Scatter plot');
if (sp)
  grokky.closeViewer(sp);  // returns 1
```

```datagrok-exec
// Close every scatter plot on the view. View is required for string/pred
// inputs because the helper has no other way to find candidates.
grokky.closeViewer('Scatter plot', view);
```

```datagrok-exec
// Close every histogram whose value column has been dropped from the DF.
grokky.closeViewer((v) => {
  if (v.type !== DG.VIEWER.HISTOGRAM) return false;
  const colName = v.props.valueColumnName;
  return colName && !t.col(colName);
}, view);
```

### Close-all (keep the grid)

```ts
grokky.closeAllViewers(view: DG.TableView | null, opts?: {keepGrid?: boolean}): number;
```

`keepGrid` defaults to `true` — `view.viewers[0]` is preserved. The helper
iterates `view.viewers.slice(1)` and closes each, returning the count.

```datagrok-exec
// Reset the view to just the grid.
const n = grokky.closeAllViewers(view);
return {closed: n};
```

## Close-and-replace pattern

Two-step idiom for "show me X instead of Y". Close the old viewer, add the
new one.

```datagrok-exec
// Close the scatter plot, then add a histogram of activities.
grokky.closeViewer('Scatter plot', view);
grokky.addViewer(view, 'Histogram', {valueColumnName: 'activity'});
```

Equivalent without the helpers (so Claude knows the long form when no
wrapper exists):

```datagrok-exec
// Long form: iterate view.viewers, filter to type, close each (skipping the grid).
Array.from(view.viewers)
  .slice(1)
  .filter((v) => v.type === DG.VIEWER.SCATTER_PLOT)
  .forEach((v) => v.close());
view.addViewer(DG.VIEWER.HISTOGRAM, {valueColumnName: 'activity'});
```

## Viewer cross-talk

Viewers sharing a DataFrame auto-sync. Click a bar in a bar chart → the
selection BitSet updates on `df.selection` → every other viewer bound to
`df` reacts. Same for filter (`df.filter`), current row, mouse-over. No
explicit linking call.

If two viewers must coordinate across **different** dataframes, see
`grok.data.linkTables(t1, t2, keys1, keys2, [SYNC_TYPE.SELECTION_TO_SELECTION, ...])`
— rare and not handled here.

## Worked demo recipes

### Recipe 1 — Iterative scatter plot tweaking

> "Add a scatter plot, show height vs weight, show the regression line,
> zoom in, color by age."

```datagrok-exec
const sp = grokky.addViewer(view, 'Scatter plot', {
  xColumnName: 'height',
  yColumnName: 'weight',
  showRegressionLine: true,
});
// Zoom and recolor on the existing viewer instead of recreating it.
grokky.configureViewer(sp, {
  xMin: 150, xMax: 200,
  yMin: 50,  yMax: 110,
  colorColumnName: 'age',
});
```

### Recipe 2 — Distribution by category

> "Show me the distribution of activities per compound class."

```datagrok-exec
// One histogram, overlaid by category.
grokky.addViewer(view, 'Histogram', {
  valueColumnName: 'activity',
  splitColumnName: 'compound_class',
});
```

If the user later wants the categories side-by-side instead of overlaid:

```datagrok-exec
const h = grokky.findViewer(view, 'Histogram');
if (h)
  grokky.configureViewer(h, {splitStack: true});
```

### Recipe 3 — Chemical space (plugin viewer)

> "Show me the chemical space of the current molecule column."

```datagrok-exec
const sp = await grok.functions.call('Chem:chemSpaceTopMenu', {
  table: t,
  molecules: t.col('smiles'),
  methodName: 'UMAP',
  similarityMetric: 'Tanimoto',
  plotEmbeddings: true,
});
// The returned `sp` is a scatter plot whose axes are UMAP coords; you can
// still configure it normally.
if (sp)
  grokky.configureViewer(sp, {colorColumnName: 'activity'});
```

## Anti-patterns

1. **`new DG.Viewer(...)`** — internal constructor over a Dart handle. Use
   `view.addViewer(type, opts)` (via `grokky.addViewer`) or
   `DG.Viewer.fromType(type, df, opts?)` for detached viewers.
2. **`view.scatterPlot(opts)` / `view.histogram(opts)` / `view.barChart(opts)`** —
   the per-type shorthands are explicitly **deprecated** in the source
   (`view.ts`). Use `view.addViewer('Scatter plot', opts)` (or
   `grokky.addViewer(view, 'Scatter plot', opts)`).
3. **`grok.shell.tv` when `view` is in scope** — the fixed helper takes
   `view` as the first argument precisely because `grok.shell.tv` is the
   *globally active* view, which may not be where the user clicked.
4. **Hardcoded property typos** — `xColumn`, `colorCol`, `regressionLine`.
   The wrapper fuzzy-matches but emits a console warning. Cite canonical
   names from `viewer.getProperties()` to avoid the noise.
5. **`grok.shell.addViewer(...)`** — no such method on the shell. Use
   `view.addViewer(...)` or `tv.addViewer(...)`.
6. **Treating "chem space" / "sequence space" / "activity cliffs" as viewer
   types** — they're functions. `grokky.addViewer(view, 'Chem space', ...)`
   silently fails with a did-you-mean miss. Use
   `grok.functions.call('Chem:chemSpaceTopMenu', {...})`.
7. **Iterating `view.viewers` without `.slice(1)`** — index 0 is the grid.
   `view.viewers.forEach(v => v.close())` closes the grid. Always skip it,
   or use `grokky.findViewers(view)` / `closeAllViewers(view)` which do.
8. **Mixing creation contexts** — `Viewer.scatterPlot(otherDf)` then
   `view.addViewer(sp)` where `otherDf !== view.dataFrame`. The viewer
   ignores `view.dataFrame`. Use `view.addViewer('Scatter plot', opts)`
   when on a TableView.
9. **`viewer.close()` on a never-attached viewer** — throws. Always either
   guard with try/catch or only close viewers you got from `view.viewers` /
   `grokky.findViewer*`. The helpers do this for you.
10. **`view.detachViewers()` / `view.resetLayout()` to "close everything"** —
    they work, but report no count. Use `grokky.closeAllViewers(view)` if
    you want to know how many were closed.
11. **Wrapping every `addViewer` call in `await`** — `view.addViewer(type)`
    is synchronous for built-in types. Plugin viewers via
    `grok.functions.call(...)` are async — that's the only `await` needed.

## Out of scope

- **Dock layout deep dive.** `view.dockManager.dock(viewer, 'right', null,
  'Title', 0.4)` repositions an attached viewer. Brief usage here is fine;
  full coverage (split sketcher / panel layouts, "sketcher on left, ADME on
  right" demos) belongs to the future `datagrok-layout` skill.
- **Viewer state save/restore.** `viewer.getOptions(true)` snapshots
  current props as `{id, type, look}`. `view.saveLayout()` and
  `view.loadLayout(layout)` round-trip the whole view. Lives in the layout
  skill.
- **Custom viewer authoring.** Subclassing `DG.JsViewer`, implementing
  `onTableAttached`, `onPropertyChanged`, `detach` — that's package
  development, not exec blocks.
- **Grid-specific cell rendering and pinned columns.** Future
  `datagrok-grid` skill.
- **The Filters viewer (`DG.VIEWER.FILTERS`)** — `view.getFiltersGroup()`
  is the right API and is covered by `datagrok-filtering`. Only add a
  Filters viewer here if the user explicitly says "add a filter panel as a
  viewer".
- **Cross-table linking** — `grok.data.linkTables` and `SYNC_TYPE` lives in
  a future cross-table skill. Single-table viewer cross-talk is automatic
  via shared DataFrame.
