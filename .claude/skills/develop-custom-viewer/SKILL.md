---
name: develop-custom-viewer
description: Author a custom JavaScript viewer in a Datagrok package — extend DG.JsViewer, declare properties, wire lifecycle, register via decorator
---

# develop-custom-viewer

## When to use

Your package needs a chart Datagrok doesn't ship (Sankey, word cloud,
WebLogo, plate map, …) bound to the current dataframe — properties on
the context panel, react to filter / selection / resize, persist with
layout. Python/R chart → *scripting viewer* (last step). EXISTING
viewer → `manipulate-viewers`. File-share preview →
`create-custom-file-viewers`.

## Prerequisites

- Package scaffold (`grok create <Name> --ts`); run from package root.
- `datagrok-tools` ≥ `4.12.x` for `@grok.decorators.viewer` (`DG-FACT-193`).
- `DG.JsViewer` (base at `js-api/src/viewer.ts:379`): no-arg ctor,
  `this.root` is `ui.box()`, `this.subs: Subscription[]` pre-initialized
  (`DG-FACT-186`, `DG-FACT-191`).

## Steps

1. **Scaffold the viewer file.**
   ```bash
   grok add viewer AwesomeViewer
   ```
   Expected: `src/awesome-viewer.ts` with a `DG.JsViewer` subclass; class
   names end in `Viewer` by convention (`DG-FACT-186`).

2. **Declare properties with the typed helpers.**
   Helpers on `DG.JsViewer`: `int`, `float`, `string`, `stringList`,
   `bool`, `dateTime`, plus `column(...)` (`DG-FACT-187`). **Column-bound
   properties are always `string`** — they store the column NAME, not
   the column dtype (`DG-FACT-188`). The article's worked example
   contradicts this with `this.int('valueColumnName', 'age')` — drift
   `DG-FACT-DRIFT-076`; do not copy.
   ```typescript
   import * as DG from 'datagrok-api/dg';   // + ui, grok as needed
   @grok.decorators.viewer({
     name: 'Awesome', description: 'Aggregated bar chart',
     icon: 'icons/awesome.svg', toolbox: true,
   })
   export class AwesomeViewer extends DG.JsViewer {
     splitColumnName!: string;
     valueColumnName!: string;     // string, not int (DG-FACT-DRIFT-076)
     valueAggrType!: string;
     color!: string;
     initialized = false;

     constructor() {
       super();
       this.splitColumnName = this.string('splitColumnName', null,
         {nullable: false, columnTypeFilter: DG.TYPE.STRING});
       this.valueColumnName = this.string('valueColumnName', null,
         {nullable: false, columnTypeFilter: DG.TYPE.NUMERICAL});
       this.valueAggrType = this.string('valueAggrType', 'avg',
         {choices: ['avg', 'count', 'sum']});
       this.color = this.string('color', 'steelblue',
         {choices: ['darkcyan', 'seagreen', 'steelblue']});
       this.addRowSourceAndFormula();   // standard rowSource (DG-FACT-192)
     }
   }
   ```
   Expected: properties auto-group into context-panel tabs by name
   pattern — `*ColumnName` → **Data**, `*color` → **Colors**, `*axis*`
   → **Axes**, `legend*` → **Legend**, `*margin*` → **Margins**,
   `*marker*` → **Markers**, `title|description` → **Description**,
   else **Misc** (`DG-FACT-189`). Override with `{category: '<Tab>'}`.

3. **Wire lifecycle and subscriptions in `onTableAttached`.**
   Push every subscription into `this.subs`; default `detach()`
   unsubscribes them — only override `detach` if you super-call
   (`DG-FACT-190`, `DG-FACT-191`). Wrap noisy streams in `DG.debounce`.
   ```typescript
   onTableAttached() {
     this.init();
     this.subs.push(DG.debounce(this.dataFrame.selection.onChanged, 50)
       .subscribe(() => this.render()));
     this.subs.push(DG.debounce(this.dataFrame.filter.onChanged, 50)
       .subscribe(() => this.render()));
     this.subs.push(DG.debounce(ui.onSizeChanged(this.root), 50)
       .subscribe(() => this.render(false)));
     this.render();
   }
   onPropertyChanged(p: DG.Property) {
     super.onPropertyChanged(p);          // always super-call first
     if (this.initialized) this.render();
   }
   ```
   Expected: closing the viewer triggers default `detach` → all three
   subscriptions unsubscribe; no leaked listeners on the dataframe.

4. **Implement `render(computeData)` and respect the filter.**
   Split data prep from drawing so resize doesn't re-aggregate; use
   `whereRowMask(this.dataFrame.filter)` to keep aggregation in sync.
   Add `ui.tooltip.showRowGroup(...)` on hover and
   `dataFrame.selection.handleClick(...)` on mousedown so the viewer
   highlights and selects across sibling viewers.
   ```typescript
   render(computeData = true) {
     if (computeData) {
       this.aggregated = this.dataFrame.groupBy([this.splitColumnName])
         .whereRowMask(this.dataFrame.filter)
         .add(this.valueAggrType, this.valueColumnName, 'result')
         .aggregate();
     }
     $(this.root).empty();          // draw into this.root w/ d3/echarts/svg
   }
   ```
   Expected: filter toggle re-aggregates; resize re-draws without
   recomputing; hover highlights matching rows in sibling viewers.

5. **Confirm registration form and (if needed) docking position.**
   The decorator emits a `package.g.ts` wrapper with `//name:`,
   `//tags: viewer`, `//output: viewer result`, plus
   `//meta.icon|toolbox|trellisable|viewerPath` (`DG-FACT-194`). Default
   top-menu path: `Add > JavaScript Viewers > <Pkg> > <Name>`; override
   with `viewerPath: 'Subcategory | Friendly Name'`. **`viewerPosition`
   is NOT a typed decorator option** (`DG-FACT-DRIFT-077`) — annotate
   the generated wrapper header-form:
   `//meta.viewerPosition: bottom` (`top|bottom|left|right|fill|auto`,
   `DG-FACT-195`). Use `grok.shell.registerViewer(...)` (`DG-FACT-193`)
   ONLY from a `meta.role: autostart` function when you need a
   synchronous handle inside package code.

6. **Build and smoke-test.**
   ```bash
   npm install && npm run build
   grok publish dev          # add --release once stable
   ```
   In the Datagrok console:
   ```javascript
   grok.shell.addTableView(grok.data.demo.demog()).addViewer('Awesome');
   ```
   Expected: `package.g.ts` contains
   `export function _AwesomeViewer() { return new AwesomeViewer(); }`
   with `//meta.role: viewer` + `//tags: viewer`; the viewer renders;
   context panel shows properties under Data/Colors/Misc.

7. **(Alternative) Scripting viewer in Python / R / Julia.**
   `scripts/<name>.py` with `# tags: viewers` and `# output: graphics`
   registers the script as a viewer (`DG-FACT-196`). Platform auto-adds
   `Refresh on Filter`, `Title`, `Description`, `Description Position`,
   `Description Visibility Mode` (`DG-FACT-197`).

## Common failure modes

- **`Argument of type 'string' is not assignable to parameter of type
  'number'` on `this.int('valueColumnName', 'age')`.** Drift
  `DG-FACT-DRIFT-076`. Column-name properties are STRINGS regardless of
  dtype — use `this.string('valueColumnName', 'age', {columnTypeFilter:
  DG.TYPE.NUMERICAL})`.
- **Memory leak / "viewer still firing after close".** Subscription not
  pushed into `this.subs`, or `detach()` overridden without
  `super.detach()`. Default `detach` only unsubscribes what's in
  `this.subs` (`DG-FACT-190`, `DG-FACT-191`).
- **`viewerPosition` set on the decorator is ignored.** Typed signature
  exposes only `name|description|icon|toolbox|trellisable|viewerPath`
  (`DG-FACT-DRIFT-077`). Annotate generated `package.g.ts` with
  `//meta.viewerPosition: <pos>` instead.
- **Property lands in the wrong context-panel tab.** Auto-grouping is
  by name pattern (`DG-FACT-189`), not type. Rename or pass
  `{category: 'Colors'}` explicitly.
- **`view.addViewer('Awesome')` returns blank.** Non-core viewer
  resolution is async until the package loads (`DG-FACT-193`).
  Pre-register via `grok.shell.registerViewer(...)` from a
  `meta.role: autostart` function for synchronous attach.
- **`super.onPropertyChanged(p)` skipped.** Layout serialization and
  undo break silently. Always super-call first (`DG-FACT-190`).

## Verification

- `npm run build` exits `0`; regenerated `src/package.g.ts` contains
  `//meta.role: viewer` + `//tags: viewer` for the wrapper.
- `Add → JavaScript Viewers → <Pkg> → Awesome` opens the viewer;
  properties grouped under Data / Colors / Misc.
- Toggle a filter — chart re-aggregates; resize — re-draws without
  recomputing (log in `render(computeData)`). Close viewer; dataframe's
  `selection.onChanged` listener count drops back.

## See also

- Source articles:
  - `help/develop/how-to/viewers/develop-custom-viewer.md`
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-186`…`DG-FACT-197` and drifts `DG-FACT-DRIFT-076`,`-077`.
- Reference packages:
  - `packages/Charts/src/viewers/sankey/sankey.ts:53-130` — decorator
    form, `addRowSourceAndFormula`, `columnTypeFilter`.
  - `packages/Charts/src/viewers/word-cloud/word-cloud-viewer.ts:15-90`
    — typed property bag with `{choices}` and `{min}`.
  - `packages/PowerGrid/src/package.g.ts:225` — header-form
    `//meta.viewerPosition: bottom` workaround for `DG-FACT-DRIFT-077`.
- Related skills:
  - `manipulate-viewers` (sibling — attach/configure EXISTING viewers).
  - `create-custom-file-viewers` (sibling — file-share preview viewers).
