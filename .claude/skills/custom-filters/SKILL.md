---
name: custom-filters
description: Implement and register a custom Datagrok column filter that participates in the platform's collaborative filter chain
---

# custom-filters

## When to use

A built-in filter (categorical, histogram, search, sketcher) does not
fit the column's domain — radio-button single-choice picker,
multi-value checkbox list against `' | '`-joined cells, substructure
search, date-range with custom semantics — and the filter must compose
with every other filter in the *Filters* viewer.

## Prerequisites

- A package scaffold (`grok create <Name>`); run from the package root.
- `datagrok-api` imports (`* as DG`, `* as ui`, `* as grok`).
- `datagrok-tools ^4.12.x` for the `@grok.decorators.filter` shortcut
  (`DG-FACT-267`); older toolchains use the function-role surface.

## Steps

1. **Subclass `DG.Filter`; implement the two abstract members.**
   `DG.Filter` extends `Widget`; `filterSummary` (caption text) and
   `applyFilter()` (writes the mask) are abstract and MUST be
   overridden (`DG-FACT-265`). Constructor calls `super()` (no args
   — base assigns `this.root = ui.div()`); subclass reassigns
   `this.root` to its own container and resets `this.subs = []`.

   ```typescript
   // src/filters/radio-button-filter.ts
   import * as ui from 'datagrok-api/ui';
   import * as DG from 'datagrok-api/dg';

   export class RadioButtonFilter extends DG.Filter {
     constructor() {
       super();
       this.root = ui.divV([], 'd4-radio-button-filter');
       this.subs = [];
     }

     get isFiltering()   { return true; }                                  // step 4
     get filterSummary() { return this.column!.getCategory(this.checkedId); }

     attach(dataFrame: DG.DataFrame) {
       super.attach(dataFrame);                                            // DG-FACT-272: super FIRST
       this.column     = DG.Utils.firstOrNull(dataFrame.columns.categorical);
       this.columnName = this.column!.name;
       this.render();
     }

     applyState(state: any) { super.applyState(state); this.render(); }    // DG-FACT-273
     // applyFilter(), render() — steps 3 and 4
   }
   ```
   Expected: TS compiles; `applyFilter is abstract` errors gone.

2. **Register the filter as a package function with `role: filter`.**
   Datagrok discovers filters by *function role*, not by class name —
   the function MUST declare `meta.role: filter` AND output type
   `filter` (`DG-FACT-266`). Function-decorator form (canonical for
   Widgets):

   ```typescript
   // src/package.ts
   export class PackageFunctions {
     @grok.decorators.func({
       name: 'Single Choice',
       outputs: [{type: 'filter', name: 'result'}],
       meta: {role: 'filter'},
     })
     static radioButtonFilter() { return new RadioButtonFilter(); }
   }
   ```

   The class-decorator shortcut `@grok.decorators.filter({name?,
   description?, semType?})` implicitly applies `meta.role: filter`
   plus the `filter` output (`DG-FACT-267`); no `package.ts` factory
   needed. Do NOT hand-author the bare-function `//name:` /
   `//meta.role:` header block shown in the article — that is the
   auto-emitted `package.g.ts` shape (`DG-FACT-DRIFT-CF-001`).
   `package.g.ts` is generated on `npm run build` and IS committed.

3. **Filter additively in `applyFilter`; never set bits TRUE.**
   Multiple filters compose by mask AND. Each filter MUST only flip
   bits to FALSE. Iterate `dataFrame.rowCount`, evaluate the
   predicate, and call `filter.set(i, <pass>, false)` (third arg =
   `notify=false`); when the loop is done, call `filter.fireChanged()`
   once (`DG-FACT-270`). Base-class JSDoc: "should disregard false
   values (these are filtered out already by other filters), and
   should filter out corresponding indexes."

   ```typescript
   applyFilter() {
     const indexes = this.column!.getRawData();
     const filter  = this.dataFrame!.filter;
     for (let i = 0; i < this.dataFrame!.rowCount; i++)
       filter.set(i, indexes[i] === this.checkedId, false);
     this.dataFrame!.filter.fireChanged();
   }
   ```

4. **Trigger re-filtering via `requestFilter`, NOT a direct call.**
   When the UI changes, request a chain re-run with
   `this.dataFrame!.rows.requestFilter()` — this fires
   `onRowsFiltering`, and the platform invokes `applyFilter()` on
   every filter in the group in order. Calling `applyFilter()`
   directly bypasses the chain and produces wrong results when other
   filters are active (`DG-FACT-271`):

   ```typescript
   radio.on('change', () => this.dataFrame!.rows.requestFilter());
   ```

   Optionally override `get isFiltering()` — return `super.isFiltering
   && <your predicate>` to short-circuit when nothing is selected. The
   base getter checks the `d4-filter-disabled` / `d4-filters-disabled`
   CSS classes set by *FilterGroup* when the user toggles the panel
   checkbox; chaining through `super` preserves that behaviour
   (`DG-FACT-274`).

5. **Build, publish, attach to a TableView.**
   ```bash
   npm install && npm run build && grok publish <host>
   ```
   Set `type: '<Package>:<funcName>'` on the *Filters* viewer
   (`DG-FACT-268`). When the filter lives in a different package than
   the caller, `await grok.functions.call('<Package>:<funcName>')`
   FIRST so the function is loaded before the synchronous
   `view.filters(...)` call (`DG-FACT-269`). Prefer
   `addViewer(Viewer.filters(...))` over the `@deprecated`
   `view.filters(...)` (`DG-FACT-DRIFT-CF-002`):

   ```typescript
   await grok.functions.call('Widgets:radioButtonFilter');
   const tv = grok.shell.addTableView(grok.data.demo.demog());
   tv.addViewer(DG.Viewer.filters({
     filters: [{type: 'Widgets:radioButtonFilter', columnName: 'race'}],
   }));
   ```

## Common failure modes

- **Filter never appears in *Add filter…* picker.** Function role is
  wrong or missing — inspect `src/package.g.ts`: the entry MUST
  contain `//meta.role: filter` AND `//output: filter result`. Both
  tokens are case-sensitive (`DG-FACT-266`).
- **Other filters' selections clear when this filter runs.**
  `applyFilter` is setting bits to TRUE (e.g. `filter.set(i, true,
  false)`) instead of only flipping its own failing rows. Never set
  TRUE — only `filter.set(i, <predicate>, false)` so already-FALSE
  bits remain FALSE (`DG-FACT-270`).
- **UI changes don't refilter.** Change handler calls
  `this.applyFilter()` directly. Switch to
  `this.dataFrame!.rows.requestFilter()` so the whole chain re-runs
  (`DG-FACT-271`).
- **`Cannot read properties of undefined` in `attach`.** Subclass
  read `this.dataFrame` / `this.column` before calling
  `super.attach(dataFrame)` — call `super.attach` FIRST, then assign
  `this.column` / `this.columnName`, then build UI (`DG-FACT-272`).
- **Cross-package filter throws "function not found".** Caller
  invoked `view.filters({...})` synchronously without first loading
  the filter's package. Add `await grok.functions.call('<Package>:<func>')`
  before the call (`DG-FACT-269`).
- **State doesn't restore on layout reload.** `applyState` not
  overridden, or overridden without `super.applyState(state)`.
  Override to call `super.applyState(state)` then re-render
  (`DG-FACT-273`); platform calls `applyState` AFTER `attach`.

## Verification

- `npm run build` and `grok publish <host>` exit `0`;
  `src/package.g.ts` contains a wrapper with `//meta.role: filter`
  and `//output: filter result`.
- In Datagrok: open the *Filters* viewer, click *Add filter* — your
  filter appears in the menu and paints its UI. Toggling it reduces
  visible rows; toggling another filter further reduces them
  (composition works) instead of resetting them.

## See also

- Source: `help/develop/how-to/viewers/custom-filters.md` (mirror at
  `docs/_internal/articles-mirror/...`); related
  `help/visualize/viewers/filters.md`; API `js-api/src/widgets/filter.ts`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` — facts
  `DG-FACT-265` … `DG-FACT-275` and drifts `DG-FACT-DRIFT-CF-001` …
  `DG-FACT-DRIFT-CF-003`.
- Reference packages:
  - `packages/Widgets/src/filters/radio-button-filter.ts:12-78` —
    categorical, single-choice, constant `isFiltering`.
  - `packages/Widgets/src/filters/multi-value-filter.ts:6-93` —
    multi-value parsing, computed `isFiltering`, checkbox UI.
  - `packages/Widgets/src/package.ts:14-33` + `package.g.ts:4-18` —
    function-decorator registration and auto-emitted header form.
  - `packages/ApiSamples/scripts/ui/viewers/filters/custom-filters.js` —
    cross-package invocation pattern.
- Related skills: `custom-cell-renderers` (same role-based
  registration pattern, different role token).
