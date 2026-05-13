---
name: manipulate-viewers
version: 0.1.0
description: |
 Drive Datagrok viewers (scatter plot, histogram, bar chart, custom)
 from package or app code — attach to a `TableView`, set options,
 read property metadata, dock to a panel side, and persist custom
 on-attach behavior across layout restores. For plugin authors who
 need to script a chart instead of letting users add it from the UI.
 Use when asked to "open a chart on the current table from code",
 "lay out plots in a panel programmatically", or "make a viewer
 re-draw the same overlay after the layout reloads".
triggers:
 - open a plot from code
 - configure a chart programmatically
 - dock a chart in a panel
 - read viewer properties from a package
 - persist a custom overlay across layouts
 - drive native and custom charts from a function
allowed-tools:
 - Read
 - Edit
 - Bash
harness-authored: true
---

# manipulate-viewers

## When to use

Scripting a viewer from a package, app, or function — attach a chart
to a `TableView`, dock it next to another, read its property state,
or wire an overlay that re-runs on every layout restore.

## Prerequisites

- A package scaffold (`grok create <Name>`); code lives under `src/`.
- `datagrok-api` imports — the article omits these; snippets won't
 compile without them:
 ```typescript
 import * as grok from 'datagrok-api/grok';
 import * as DG from 'datagrok-api/dg';
 ```
- A `DG.DataFrame` (`grok.data.demo.demog` for demos) and a
 `TableView` (`grok.shell.addTableView(df)` or `grok.shell.tv`).

## Steps

1. **Use the universal attach path: `view.addViewer(...)`.**
 First arg is a string, `DG.VIEWER.<TYPE>` constant, or a `Viewer`
 instance — works for ALL viewers, native and custom (`DG-FACT-198`).
 The per-type wrappers `view.histogram(opts)`, `view.scatterPlot(opts)`,
 `view.barChart(opts)` featured in the article are JSDoc-marked
 `@deprecated use addViewer(...)` and slated for removal:
 ```typescript
 const data = grok.data.demo.demog;
 const view = grok.shell.addTableView(data);
 view.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'});
 view.addViewer('Leaflet', { // custom viewers too
 latitudeColumnName: 'lat', longitudeColumnName: 'lng',
 renderType: 'heat map',
 });
 ```
 Expected: viewer mounts in the table view. Pattern matches
 `packages/Chem/src/package.ts:1207-1221`.

2. **Choose the right factory when you need the instance first.**
 `DG.Viewer.fromType(type, df, opts?)` is SYNC and only works for
 core viewers in `DG.CORE_VIEWER` (`DG-FACT-199`, `DG-FACT-206`).
 `df.plot.fromType(type, opts?)` is ASYNC — required for
 package-shipped types like `GLOBE`, `GOOGLE_MAP`, `WORD_CLOUD`,
 `TIMELINES`, `RADAR_VIEWER`, `SURFACE_PLOT`, `SCAFFOLD_TREE`:
 ```typescript
 const sp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, data); // sync
 view.addViewer(sp);
 const wc = await data.plot.fromType(DG.VIEWER.WORD_CLOUD); // async
 view.addViewer(wc);
 ```
 `df.plot` also has typed sync shortcuts (`scatter`, `histogram`,
 `bar`, `box`, `line`, `heatMap`, `network`, `grid`, `tile`, `form`)
 returning the typed subclass (`DG-FACT-200`). Names differ from
 `DG.Viewer` statics — `df.plot.bar` not `barChart`.

3. **Pass / read options as a flat map; values live under `.look`.**
 `setOptions(map)` takes a flat bag (no `{look:...}` wrapper).
 `getOptions(includeDefaults = false)` returns `{id, type, look:{...}}` —
 pass `true` for an exhaustive snapshot, `false` (default) keeps the
 payload small for serialization (`DG-FACT-201`):
 ```typescript
 const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT) as DG.ScatterPlotViewer;
 sp.setOptions({xColumnName: 'weight', yColumnName: 'height'});
 const opts = sp.getOptions(true); // {id, type, look:{xColumnName:...,...}}
 const xCol = opts.look['xColumnName'];
 ```
 Pattern: `packages/Chem/src/package.ts:1253` reads
 `sp.getOptions.look['xColumnName']`.

4. **Inspect available properties via `getProperties` — use `columnTypeFilter`.**
 Each `Property` exposes `name`, `propertyType`, `semType`,
 `description`, `defaultValue`, `choices`, and `columnTypeFilter`
 (`DG-FACT-202`). The article's `p.columnFilter` is stale; no
 backward-compat alias:
 ```typescript
 const bc = view.addViewer(DG.VIEWER.BAR_CHART);
 for (const p of bc.getProperties)
 console.log(`${p.propertyType} ${p.name}: ${p.description} [${p.columnTypeFilter}]`);
 ```
 `columnTypeFilter` is `'numerical' | 'categorical' | ColumnType | null`;
 compare with `DG.COLUMN_TYPE.*` constants, not raw strings.

5. **Dock through the right manager — never the literal `'top'`.**
 `DockManager.dock(viewer, dockType?, refNode?, title?, ratio?)`
 defaults to `LEFT` (`DG-FACT-204`). Pick by intent:
 `view.dockManager` keeps the viewer INSIDE the view; the
 `grok.shell.dockManager` floats it independently. `DG.DOCK_TYPE.TOP`
 resolves to `'up'` (NOT `'top'`) — the article's "top" misleads
 (`DG-FACT-203`, ). Always go through the enum:
 ```typescript
 const wl = view.addViewer('WebLogo', {sequenceColumnName: col.name});
 view.dockManager.dock(wl, DG.DOCK_TYPE.DOWN, null, 'Composition', 0.25);
 ```
 Pattern: `packages/Bio/src/package.ts:1050-1051`,
 `packages/PowerPack/src/search/power-search.ts:147,397`.

6. **Persist custom init via `initializationFunction` (re-runs on restore).**
 `initializationFunction` is a built-in viewer property whose value
 is a registered Datagrok function name. The function takes a single
 `viewer` input typed as the right subclass and runs every time the
 viewer is constructed — including on layout/project restore
 (`DG-FACT-205`):
 ```typescript
 export class PackageFunctions {
 @grok.decorators.func({name: 'highlightInitFunction'})
 static async highlightInitFunction(
 @grok.decorators.param({type: 'viewer'}) sp: DG.ScatterPlotViewer,
 ): Promise<void> {
 sp.onAfterDrawScene.subscribe( => {
 const ctx = sp.canvas.getContext('2d')!;
 ctx.fillStyle = 'red'; ctx.fillRect(100, 100, 50, 50);
 });
 }
 }
 view.addViewer(DG.VIEWER.SCATTER_PLOT, {
 xColumnName: 'weight', yColumnName: 'height',
 initializationFunction: 'highlightInitFunction',
 });
 ```
 Pattern: `packages/Chem/src/package.ts:1218,1238-1242`.

## Common failure modes

- **Package viewer attaches but `setOptions` silently no-ops.** You
 used `DG.Viewer.fromType` (sync, core-only) for a packaged viewer,
 or `view.addViewer('CustomType', opts)` returned before the package
 loaded. Switch to `await df.plot.fromType('Type')` and configure on
 the resolved instance (`DG-FACT-199`, `DG-FACT-206`).
- **`view.dockManager.dock(viewer, 'top')` throws or no-ops.** The
 position string is `'up'`. Use `DG.DOCK_TYPE.TOP` (`DG-FACT-203`,
 ).
- **`p.columnFilter` is `undefined` in `getProperties`.** Accessor
 was renamed to `columnTypeFilter`; no alias.
- **`view.histogram(...)` / `view.scatterPlot(...)` flagged deprecated
 by `tsc`.** All `View.<viewer>` wrappers carry `@deprecated`;
 rewrite as `view.addViewer(DG.VIEWER.<TYPE>, {...})`.
- **Custom overlay runs once but vanishes on layout reload.** You
 subscribed after `addViewer` instead of passing
 `initializationFunction`. Move the body into a registered package
 function and pass its name in viewer options (`DG-FACT-205`).

## Verification

- `npm run build` (or `grok check`) exits `0` with no `deprecated`
 warnings on `view.<type>(...)` — every attach goes through `addViewer`.
- In Datagrok, the viewer mounts; `viewer.getOptions(true).look`
 returns your options plus defaults.
- `viewer.getProperties.every(p => 'columnTypeFilter' in p)` is `true`.
- Save the layout, reopen — the viewer reappears AND the
 `initializationFunction` behavior re-runs.
- `view.dockManager.dock(viewer, DG.DOCK_TYPE.TOP)` actually docks at
 the top (proves you didn't paste the literal `'top'`).

## See also

- Source: `help/develop/how-to/viewers/manipulate-viewers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
 `DG-FACT-198`…`DG-FACT-206`/`079`/`080`.
- Reference packages:
 - `packages/Chem/src/package.ts:1207-1221,1238-1242` — scatter plot
 attach + `initializationFunction: 'activityCliffsInitFunction'` and
 the matching decorator-form init function.
 - `packages/Bio/src/package.ts:1050-1051` — WebLogo docked via
 `dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'Composition', 0.25)`.
 - `packages/PowerPack/src/viewers-gallery.ts:387` —
 `view.addViewer(DG.Viewer.fromType(viewer, table))` instance form.
 - `packages/PowerPack/src/search/power-search.ts:147,397` —
 `grok.shell.dockManager.dock(...)` for view-independent docks.
- Related: `develop-custom-viewer` (defines the `JsViewer` this skill
 attaches and configures).
