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
   Typed `this.int/float/string/bool/dateTime/columnList/…` helpers
   register-and-return the initial value (see `DG-FACT-186`, `187`).
   Column-bound names must end `ColumnName` (`DG-FACT-188`); tab
   grouping is by name pattern (`DG-FACT-189`).
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

2. **Wire lifecycle hooks; push every subscription onto `this.subs`.**
   Default `detach()` unsubscribes the bag for you (see `DG-FACT-190`,
   `DG-FACT-191`).
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

3. **Render against the live filter; gate compute on a flag.** Always
   chain `.whereRowMask(this.dataFrame.filter)` (see `DG-FACT-461`).
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

4. **Bridge mouse events to platform state** (see `DG-FACT-460`).
   ```typescript
   bars.on('mouseover', (event, d) => ui.tooltip.showRowGroup(this.dataFrame,
       i => d.category === this.dataFrame.getCol(this.splitColumnName).get(i),
       event.x, event.y))
     .on('mouseout', () => ui.tooltip.hide())
     .on('mousedown', (event, d) => this.dataFrame.selection.handleClick(
       i => d.category === this.dataFrame.getCol(this.splitColumnName).get(i),
       event));
   ```

5. **Register so it appears in *Add Viewer*.** Use
   `@grok.decorators.viewer({...})` on the class, or an annotated
   factory (see `DG-FACT-193`, `DG-FACT-194`). Commit `package.g.ts`
   (not gitignored). `meta.viewerPosition` is header-form only — patch
   into `package.g.ts` post-build (see `DG-FACT-195`).

6. **Alternative: scripting viewer (server-rendered).** Python/R/Julia
   under `<pkg>/scripts/` with `# language:`, `# tags: viewers`,
   `# output: graphics` required (see `DG-FACT-196`, `DG-FACT-197`).

## Common failure modes

- Not in *Add Viewer* — `package.g.ts` missing `//meta.role: viewer` /
  `//output: viewer result`; rebuild + commit (`DG-FACT-193`).
- Memory leak on close — subscription wasn't pushed onto `this.subs`
  (`DG-FACT-191`).
- Stale rows after filter — missing `whereRowMask(...filter)`
  (`DG-FACT-461`).
- Property absent from Data tab — name lacks `ColumnName` suffix
  (`DG-FACT-188`).
- Scripting viewer blank — header missing `# output: graphics` or
  `# tags: viewers` (`DG-FACT-196`).

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
