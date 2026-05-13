---
name: develop-custom-viewer
version: 0.1.0
description: |
  Subclass `DG.JsViewer` to ship a DataFrame-bound chart the platform
  doesn't have built-in (D3 bar chart, custom domain plot) — or declare
  a Python/R/Julia scripting viewer for a server-rendered visualization.
  Produces a class wired to dataframe filter/selection/size events plus
  the package-function registration that surfaces it in the *Add Viewer*
  menu and lets users persist it in layouts.
  Use when asked to "add a new chart kind to the Add menu",
  "ship a custom D3 visualization for a dataframe", or "register a
  domain-specific plot type the platform doesn't have built-in".
triggers:
  - new chart kind for dataframe
  - add custom d3 visualization
  - register chart in add menu
  - domain-specific plot type
  - ship a python scripting plot
  - dataframe-bound visualization
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# develop-custom-viewer

## When to use

A built-in viewer (scatter, bar, histogram, line, network, …) can't
express the picture your data needs — a D3 bar chart with custom
interactions, a domain plot (radar, sankey), or a server-rendered
matplotlib scene. Output must live in *Add Viewer*, persist in
layouts, and react to filter/selection/current row.

## Steps

1. **Subclass `DG.JsViewer`; declare properties in the constructor.**
   Base class at `js-api/src/viewer.ts:379` (`DG-FACT-186`). Each typed
   helper — `this.int`, `this.float`, `this.string`, `this.stringList`,
   `this.bool`, `this.dateTime`, `this.columnList` (`viewer.ts:457-499`,
   `DG-FACT-187`) — both registers the property in the context panel
   AND returns its initial value. Column-bound properties MUST end with
   `ColumnName`; that suffix flags them as Data (`viewer.ts:455-459`,
   `DG-FACT-188`). Properties auto-group into tabs by name pattern
   (Data / Colors / Axes / Legend / Margins / Markers / Description /
   Misc — `DG-FACT-189`); override via `{category: 'X'}`.

   ```typescript
   import * as DG from 'datagrok-api/dg';
   import * as grok from 'datagrok-api/grok';

   @grok.decorators.viewer({name: 'Awesome', icon: 'icons/awesome.svg'})
   export class AwesomeViewer extends DG.JsViewer {
     splitColumnName = this.string('splitColumnName', 'site');
     valueColumnName = this.string('valueColumnName', 'age');
     valueAggrType = this.string('valueAggrType', 'avg', {choices: ['avg', 'count', 'sum']});
     color = this.string('color', 'steelblue', {choices: ['darkcyan', 'seagreen', 'steelblue']});
     initialized = false;
   }
   ```
   Expected: properties appear in the context panel under their
   auto-derived tabs when the viewer is attached.

2. **Wire lifecycle hooks; push every subscription onto `this.subs`.**
   `onTableAttached()` runs when a DataFrame binds — do data-dependent
   init here (`DG-FACT-190`, `viewer.ts:412-444`). Push every
   `selection.onChanged` / `filter.onChanged` / `ui.onSizeChanged(...)`
   subscription onto `this.subs`; default `detach()` unsubscribes them
   (`DG-FACT-191`, `viewer.ts:397-398,429-432`). Wrap noisy streams
   with `DG.debounce(obs, 50)`. Don't override `detach()` without
   `super.detach()` — leaks.

   ```typescript
   onTableAttached() {
     this.init();
     this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50).subscribe(_ => this.render()));
     this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50).subscribe(_ => this.render()));
     this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50).subscribe(_ => this.render(false)));
     this.render();
   }
   onPropertyChanged(property: DG.Property) {
     super.onPropertyChanged(property);
     if (this.initialized) this.render();
   }
   ```

3. **Render against the live filter; gate compute on a flag.**
   In `render(computeData = true)` recompute the aggregated frame ONLY
   when `computeData` (resize events skip the work). Chain
   `groupBy([col]).whereRowMask(this.dataFrame.filter).add(...).aggregate()`
   so the viewer respects user filters — skipping `whereRowMask`
   aggregates over the full frame regardless (`DG-FACT-461`,
   `stats.ts:350-354`). Clear `this.root` and re-draw to its size.

   ```typescript
   render(computeData = true) {
     if (computeData) {
       this.aggregatedTable = this.dataFrame
         .groupBy([this.splitColumnName])
         .whereRowMask(this.dataFrame.filter)
         .add(this.valueAggrType, this.valueColumnName, 'result')
         .aggregate();
       /* populate this.data from aggregatedTable.getCol('result') */
     }
     $(this.root).empty();
     /* ... d3 draw using this.data, this.color, this.root size ... */
   }
   ```

4. **Bridge mouse events into platform state.**
   `ui.tooltip.showRowGroup(dataFrame, predicate, x, y)` pops the
   standard row-group tooltip; `dataFrame.selection.handleClick(predicate,
   event)` extends/replaces selection per Ctrl/Shift/Meta — pass the raw
   `MouseEvent` (`DG-FACT-460`, `ui.ts:1606-1608`, `bit-set.ts:220-222`).

   ```typescript
   bars.on('mouseover', (event, d) => ui.tooltip.showRowGroup(this.dataFrame,
       i => d.category === this.dataFrame.getCol(this.splitColumnName).get(i),
       event.x, event.y))
     .on('mouseout', () => ui.tooltip.hide())
     .on('mousedown', (event, d) => this.dataFrame.selection.handleClick(
       i => d.category === this.dataFrame.getCol(this.splitColumnName).get(i),
       event));
   ```

5. **Register so it appears in *Add Viewer*.** Two paths (`DG-FACT-193`):
   the `@grok.decorators.viewer({name, description?, icon?, toolbox?,
   trellisable?, viewerPath?})` class decorator (canonical:
   `packages/Charts/src/viewers/sankey/sankey.ts:53-57` →
   `packages/Charts/src/package.g.ts:129-136`); or an annotated factory
   in `package.ts` returning `new AwesomeViewer()`. Either way the build
   emits `//name:` + `//meta.role: viewer` + `//output: viewer result`
   into `src/package.g.ts` — commit it (NOT gitignored). Metadata
   (`DG-FACT-194`, `decorators/functions.ts:36-45`): `meta.icon`,
   `meta.toolbox: true`, `meta.trellisable: true` (inner viewer of
   `DG.VIEWER.TRELLIS_PLOT`), `meta.viewerPath: 'Cat | Friendly Name'`
   (overrides default `Add > JavaScript Viewers > <Pkg> > <Name>`).
   `meta.viewerPosition` (`top|bottom|left|right|fill|auto`,
   `DG-FACT-195`) is HEADER-FORM ONLY — patch into `package.g.ts`
   post-build (see `packages/PowerGrid/src/package.g.ts:225`).

6. **Alternative: scripting viewer (server-rendered).** For Python/R/
   Julia, drop a `.py`/`.r`/`.jl` under `<pkg>/scripts/` with
   `# language: python`, `# tags: viewers` (registers in Script Browser
   at `/scripts?q=%23viewers`), typed inputs (e.g. `# input: column
   splitColumnName {type: categorical}`), and `# output: graphics`
   — REQUIRED (`DG-FACT-196`). Platform auto-adds `Refresh on Filter`,
   `Title`, `Description`, `Description Position`, `Description
   Visibility Mode` (`DG-FACT-197`).

## Common failure modes

- **Viewer doesn't appear in *Add Viewer* menu.** `src/package.g.ts`
  missing `//meta.role: viewer` or `//output: viewer result` — rebuild
  after fixing the decorator; confirm `package.g.ts` is committed
  (`DG-FACT-193`).
- **Memory leak: closed viewer keeps re-rendering.** A subscription
  wasn't pushed onto `this.subs`, so default `detach()` can't
  unsubscribe it (`DG-FACT-191`).
- **Chart shows stale rows after the user filters.** Aggregation
  skipped `whereRowMask(this.dataFrame.filter)` and ran over the full
  frame (`DG-FACT-461`).
- **Property doesn't show in the Data tab.** Name lacks the
  `ColumnName` suffix — rename `splitCol` → `splitColumnName`
  (`DG-FACT-188`).
- **Scripting viewer fails to render.** Header omits `# output: graphics`
  or `# tags: viewers` (`DG-FACT-196`).

## See also

- Source: `help/develop/how-to/viewers/develop-custom-viewer.md`;
  related `help/develop/how-to/viewers/manipulate-viewers.md`,
  `help/visualize/viewers/viewers.md`; API `js-api/src/viewer.ts:376-499`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-186` … `DG-FACT-197`, `DG-FACT-460`, `DG-FACT-461`.
- Reference: `packages/Charts/src/viewers/sankey/sankey.ts:53-100`
  (decorator + JsViewer); `.../radar/radar-viewer.ts:20-240` (full
  lifecycle); `packages/Charts/src/package.g.ts:129-136` (auto-emitted
  header); `packages/PowerGrid/src/package.g.ts:220-230` (viewerPosition).
- Related skills: `manipulate-viewers`, `custom-filters`, `custom-cell-renderers`.
