---
name: manipulate-viewers
description: Attach, configure, dock, and inspect Datagrok viewers from a package using the JS API
---

# manipulate-viewers

## When to use

You are scripting a viewer from a package or app — adding a chart to a
`TableView`, docking it next to another, reading its current property
state, or wiring a custom initializer that re-runs on layout restore.
Triggers: "open a scatter plot with these columns", "dock a histogram
on the right", "attach a custom viewer from my package".

## Prerequisites

- A package scaffold (`grok create <Name>`); commands run from the
  package root, code lives under `src/`.
- `datagrok-api` imports available — the article omits these and
  snippets won't compile as written:
  ```typescript
  import * as grok from 'datagrok-api/grok';
  import * as DG   from 'datagrok-api/dg';
  ```
- A `DG.DataFrame` to plot (`grok.data.demo.demog()` for demos) and a
  `TableView` (`grok.shell.addTableView(df)` or `grok.shell.tv`).

## Steps

1. **Use the universal attach path: `view.addViewer(...)`.**
   First arg is a string, `DG.VIEWER.<TYPE>` constant, or a `Viewer`
   instance — works for ALL viewers, native and custom
   (`DG-FACT-198`). The per-type wrappers (`view.histogram(opts)`,
   `view.scatterPlot(opts)`, `view.barChart(opts)`) the article
   features are JSDoc-marked `deprecated: use addViewer(...)` and
   slated for removal in 1.21 (`DG-FACT-DRIFT-080`):
   ```typescript
   const data = grok.data.demo.demog();
   const view = grok.shell.addTableView(data);
   view.addViewer(DG.VIEWER.HISTOGRAM, {value: 'age'});
   view.addViewer('Leaflet', {                          // also custom viewers
     latitudeColumnName: 'lat', longitudeColumnName: 'lng',
     renderType: 'heat map',
   });
   ```
   Expected: viewer mounts in the table view; pattern matches
   `packages/Chem/src/package.ts:1207-1221`.

2. **Choose the right factory when you need the instance first.**
   `DG.Viewer.fromType(type, df, options?)` is SYNC and only works for
   the 25 viewers in `DG.CORE_VIEWER` (`DG-FACT-199`, `DG-FACT-206`).
   `df.plot.fromType(type, options?)` is ASYNC and is required for
   the 7 package-shipped types: `GLOBE`, `GOOGLE_MAP`, `WORD_CLOUD`,
   `TIMELINES`, `RADAR_VIEWER`, `SURFACE_PLOT`, `SCAFFOLD_TREE`.
   ```typescript
   // CORE viewer — sync
   const sp = DG.Viewer.fromType(DG.VIEWER.SCATTER_PLOT, data);
   view.addViewer(sp);

   // package viewer — async; await before configuring
   const wc = await data.plot.fromType(DG.VIEWER.WORD_CLOUD);
   view.addViewer(wc);
   ```
   `df.plot` ALSO has 10 sync shortcuts (`scatter`, `histogram`, `bar`,
   `box`, `line`, `heatMap`, `network`, `grid`, `tile`, `form`) that
   return the typed `Viewer` subclass (`DG-FACT-200`). Names differ
   from `DG.Viewer` statics — `df.plot.bar` not `barChart`.

3. **Pass / read options as a flat map; values live under `.look`.**
   `setOptions(map)` accepts a flat property bag (no `{look: ...}`
   wrapper). `getOptions(includeDefaults = false)` returns
   `{id, type, look: {...}}` — pass `true` for an exhaustive snapshot,
   default `false` keeps the payload small for serialization
   (`DG-FACT-201`):
   ```typescript
   const sp = view.addViewer(DG.VIEWER.SCATTER_PLOT) as DG.ScatterPlotViewer;
   sp.setOptions({xColumnName: 'weight', yColumnName: 'height'});

   const opts = sp.getOptions(true);     // {id, type, look:{xColumnName:..., ...}}
   const xCol = opts.look['xColumnName'];
   ```
   Pattern: `packages/Chem/src/package.ts:1253` reads
   `sp.getOptions().look['xColumnName']`.

4. **Inspect available properties via `getProperties()` — use `columnTypeFilter`.**
   Each `Property` exposes `name`, `propertyType`, `semType`,
   `description`, `defaultValue`, `choices`, and `columnTypeFilter`
   (`DG-FACT-202`). The article's `p.columnFilter` is STALE — the
   accessor was renamed and there is no backward-compat alias
   (`DG-FACT-DRIFT-078`):
   ```typescript
   const bc = view.addViewer(DG.VIEWER.BAR_CHART);
   for (const p of bc.getProperties())
     console.log(`${p.propertyType} ${p.name}: ${p.description} [${p.columnTypeFilter}]`);
   ```
   `columnTypeFilter` is `'numerical' | 'categorical' | ColumnType | null`;
   compare with `DG.COLUMN_TYPE.*` constants, not raw strings.

5. **Dock through the right manager — and never the literal `'top'`.**
   `DockManager.dock(viewer, dockType?, refNode?, title?, ratio?)`
   defaults to `LEFT` (`DG-FACT-204`). Two managers — pick by intent:
   `view.dockManager` keeps the viewer INSIDE the view (moves with
   it); `grok.shell.dockManager` floats it independently.
   `DG.DOCK_TYPE.TOP` resolves to the string `'up'` (NOT `'top'`) —
   the article's "left | right | top | down | fill" misleads
   (`DG-FACT-203`, `DG-FACT-DRIFT-079`). Always go through the enum:
   ```typescript
   const wl = view.addViewer('WebLogo', {sequenceColumnName: col.name});
   view.dockManager.dock(wl, DG.DOCK_TYPE.DOWN, null, 'Composition', 0.25);
   ```
   Pattern: `packages/Bio/src/package.ts:1050-1051`,
   `packages/PowerPack/src/search/power-search.ts:147,397`.

6. **Persist custom init via `initializationFunction` (re-runs on restore).**
   `initializationFunction` is a built-in viewer property whose value
   is a registered Datagrok function name (`<Package>:<funcName>`).
   The function takes a single `viewer` input typed as the right
   subclass and runs every time the viewer is constructed — including
   on layout/project restore (`DG-FACT-205`):
   ```typescript
   // src/package.ts — canonical decorator form
   export class PackageFunctions {
     @grok.decorators.func({name: 'highlightInitFunction'})
     static async highlightInitFunction(
       @grok.decorators.param({type: 'viewer'}) sp: DG.ScatterPlotViewer,
     ): Promise<void> {
       sp.onAfterDrawScene.subscribe(() => {
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

- **Custom/package viewer attaches but `setOptions` silently no-ops.**
  You used `DG.Viewer.fromType` (sync, core-only) for a packaged
  viewer, or `view.addViewer('CustomType', opts)` returned before the
  package loaded. Switch to `await df.plot.fromType('Type')` and
  configure on the resolved instance (`DG-FACT-199`, `DG-FACT-206`).
- **`view.dockManager.dock(viewer, 'top')` throws / no-ops.**
  The position string is `'up'`. Use `DG.DOCK_TYPE.TOP`
  (`DG-FACT-203`, `DG-FACT-DRIFT-079`).
- **`p.columnFilter` is `undefined` while iterating
  `getProperties()`.** Accessor was renamed to `columnTypeFilter`;
  no alias remains (`DG-FACT-DRIFT-078`).
- **`view.histogram(...)` / `view.scatterPlot(...)` flagged deprecated
  by the type checker.** All `View.<viewer>()` wrappers carry
  `@deprecated`; rewrite as `view.addViewer(DG.VIEWER.<TYPE>, {...})`
  (`DG-FACT-DRIFT-080`).
- **Custom init runs once but is lost when the layout is reloaded.**
  You subscribed after `addViewer` instead of passing
  `initializationFunction`. Move the body into a registered package
  function and pass its name in viewer options (`DG-FACT-205`).

## Verification

- TypeScript build (`npm run build` or `grok check`) exits `0` and
  emits no `deprecated` warnings on `view.<type>(...)` calls — every
  attach goes through `addViewer`.
- In Datagrok, opening the view shows the viewer immediately;
  `viewer.getOptions(true).look` returns the options you passed plus
  platform defaults.
- `viewer.getProperties().every(p => 'columnTypeFilter' in p)` is
  `true` (sanity-checks the API surface you're coding against).
- Save the layout, reopen the project — the viewer reappears AND the
  `initializationFunction`-registered behavior re-runs.
- `view.dockManager.dock(viewer, DG.DOCK_TYPE.TOP)` actually docks at
  the top (proves you didn't paste the literal `'top'`).

## See also

- Source articles:
  - `help/develop/how-to/viewers/manipulate-viewers.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-198` … `DG-FACT-206` and drifts `DG-FACT-DRIFT-078`,
  `DG-FACT-DRIFT-079`, `DG-FACT-DRIFT-080`.
- Reference packages:
  - `packages/Chem/src/package.ts:1207-1221,1238-1242` —
    `view.addViewer(DG.VIEWER.SCATTER_PLOT, {...
    initializationFunction: 'activityCliffsInitFunction' ...})` and
    the matching typed decorator-form init function.
  - `packages/Bio/src/package.ts:1050-1051` — `WebLogo` docked via
    `dockManager.dock(viewer, DG.DOCK_TYPE.DOWN, null, 'Composition', 0.25)`.
  - `packages/PowerPack/src/viewers-gallery.ts:387` —
    `view.addViewer(DG.Viewer.fromType(viewer, table))` (pre-built
    instance form).
  - `packages/PowerPack/src/search/power-search.ts:147,397` —
    `grok.shell.dockManager.dock(...)` for view-independent docks.
- Related skills:
  - `develop-custom-viewer` (sibling — defines the `JsViewer` this
    skill attaches and configures).
