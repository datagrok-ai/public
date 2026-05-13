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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# manipulate-viewers

## When to use

Scripting a viewer from a package, app, or function — attach a chart
to a `TableView`, dock it next to another, read its property state,
or wire an overlay that re-runs on every layout restore.

## Prerequisites

- A `DG.DataFrame` (`grok.data.demo.demog` for demos) and a
 `TableView` (`grok.shell.addTableView(df)` or `grok.shell.tv`).

## Steps

1. **Use the universal attach path: `view.addViewer(...)`** — accepts string, `DG.VIEWER.<TYPE>` constant, or `Viewer` instance; works for native and custom (`DG-FACT-198`). Per-type wrappers like `view.histogram` are deprecated.
 ```typescript
 const data = grok.data.demo.demog;
 const view = grok.shell.addTableView(data);
 view.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'});
 view.addViewer('Leaflet', { // custom viewers too
 latitudeColumnName: 'lat', longitudeColumnName: 'lng',
 renderType: 'heat map',
 });
 ```

2. **Choose the right factory when you need the instance first.** `DG.Viewer.fromType` is sync and core-only; `df.plot.fromType` is async and required for package viewers (`DG-FACT-199`, `DG-FACT-206`). Typed `df.plot.{scatter,bar,...}` shortcuts in `DG-FACT-200`.
 ```typescript
 const sp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, data); // sync
 view.addViewer(sp);
 const wc = await data.plot.fromType(DG.VIEWER.WORD_CLOUD); // async
 view.addViewer(wc);
 ```

3. **Pass / read options as a flat map; values live under `.look`** (`DG-FACT-201`). `setOptions` is flat; `getOptions(includeDefaults?)` returns `{id, type, look:{...}}`.
 ```typescript
 const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT) as DG.ScatterPlotViewer;
 sp.setOptions({xColumnName: 'weight', yColumnName: 'height'});
 const opts = sp.getOptions(true); // {id, type, look:{xColumnName:...,...}}
 const xCol = opts.look['xColumnName'];
 ```

4. **Inspect available properties via `getProperties` — use `columnTypeFilter`** (article's `p.columnFilter` is stale, no alias; `DG-FACT-202`).
 ```typescript
 const bc = view.addViewer(DG.VIEWER.BAR_CHART);
 for (const p of bc.getProperties)
 console.log(`${p.propertyType} ${p.name}: ${p.description} [${p.columnTypeFilter}]`);
 ```

5. **Dock through the right manager — never the literal `'top'`.** `DG.DOCK_TYPE.TOP` resolves to `'up'`; always go through the enum (`DG-FACT-203`, `DG-FACT-204`). Use `view.dockManager` for in-view, `grok.shell.dockManager` to float.
 ```typescript
 const wl = view.addViewer('WebLogo', {sequenceColumnName: col.name});
 view.dockManager.dock(wl, DG.DOCK_TYPE.DOWN, null, 'Composition', 0.25);
 ```

6. **Persist custom init via `initializationFunction`** — built-in viewer property whose value is a registered function name; re-runs on layout/project restore (`DG-FACT-205`).
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

## Common failure modes

- **Package viewer `setOptions` silently no-ops** — used sync `DG.Viewer.fromType` for a package type. Switch to `await df.plot.fromType` (`DG-FACT-199`, `DG-FACT-206`).
- **`dock(viewer, 'top')` no-ops** — position is `'up'`; use `DG.DOCK_TYPE.TOP` (`DG-FACT-203`).
- **`p.columnFilter` undefined** — renamed to `columnTypeFilter` (`DG-FACT-202`).
- **`view.histogram(...)` deprecated** — rewrite as `view.addViewer(DG.VIEWER.HISTOGRAM, {...})`.
- **Custom overlay vanishes on layout reload** — use `initializationFunction` (`DG-FACT-205`).

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
