---
name: custom-filters
version: 0.1.0
description: |
  Implement and register a Datagrok column-slicing widget that
  participates in the platform's collaborative row-selection chain —
  for domains where the built-in slicers (categorical, histogram,
  search, sketcher) cannot express the predicate (fingerprint
  similarity, multi-value cells, custom ranges, domain-specific
  pickers). Produces a `DG.Filter` subclass plus the package-function
  registration that makes it appear in the *Add filter…* menu and
  compose with every other row predicate on the view.
  Use when asked to "narrow rows by a predicate the default slicers
  can't express", "add a column slicing widget that composes with
  others", or "build a domain-specific row picker for fingerprint or
  multi-value columns".
triggers:
  - narrow rows by custom predicate
  - column slicing widget
  - domain-specific row picker
  - chained row selection
  - substructure or similarity slicer
  - multi-value column slicer
allowed-tools:
  - Read
  - Write
  - Edit
  - Bash
harness-authored: true
---

# custom-filters

## When to use

A built-in slicer (categorical, histogram, search, sketcher) cannot
express your column's predicate — single-choice radio, multi-value
checkbox list, substructure / similarity search, custom date range —
and the widget must compose with every other selector in *Filters*.

## Prerequisites

- A package scaffold (`grok create <Name>`); run from the package root.
- `datagrok-api` imports (`* as DG`, `* as ui`, `* as grok`).
- `datagrok-tools ^4.12.x` for the `@grok.decorators.filter` shortcut
  (`DG-FACT-267`); older toolchains use the function-role surface.

## Steps

1. **Subclass `DG.Filter`; implement the two abstract members.**
   `DG.Filter` extends `Widget`; `filterSummary` (caption text) and
   `applyFilter()` (writes the mask) are abstract and MUST be
   overridden (`DG-FACT-265`). Constructor calls `super()` no-args
   (base assigns `this.root = ui.div()`); the subclass reassigns
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
     get isFiltering()   { return true; }                              // step 4
     get filterSummary() { return this.column!.getCategory(this.checkedId); }
     attach(dataFrame: DG.DataFrame) {
       super.attach(dataFrame);                                        // DG-FACT-272: super FIRST
       this.column     = DG.Utils.firstOrNull(dataFrame.columns.categorical);
       this.columnName = this.column!.name;
       this.render();
     }
     applyState(state: any) { super.applyState(state); this.render(); } // DG-FACT-273
   }
   ```
   Expected: TS compiles; `applyFilter is abstract` errors gone.

2. **Register the class as a package function with `role: filter`.**
   The platform discovers filters by *function role*, not class name —
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
   plus the `filter` output (`DG-FACT-267`) — no `package.ts` factory
   needed. Do NOT hand-author the bare-function `//name:` / `//meta.role:`
   header from the article: that is the auto-emitted `package.g.ts`
   shape and IS committed (not gitignored).

3. **Filter additively in `applyFilter`; never set bits TRUE.**
   Filters compose by mask AND; each filter MUST only flip bits to
   FALSE. Iterate `dataFrame.rowCount`, call `filter.set(i, <pass>,
   false)` (third arg = `notify=false`), and call `filter.fireChanged()`
   once at the end (`DG-FACT-270`). Base JSDoc: "should disregard
   false values (these are filtered out already by other filters), and
   should filter out corresponding indexes".

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
   On every UI change call `this.dataFrame!.rows.requestFilter()`;
   it fires `onRowsFiltering`, prompting the platform to invoke
   `applyFilter()` on every filter in the group. Direct calls bypass
   the chain (`DG-FACT-271`). Optionally override `get isFiltering()`
   as `super.isFiltering && <predicate>` to short-circuit when
   nothing is selected; chaining through `super` preserves the
   `d4-filter-disabled` checkbox behaviour (`DG-FACT-274`).

   ```typescript
   radio.on('change', () => this.dataFrame!.rows.requestFilter());
   ```

5. **Build, publish, attach to a TableView.**
   ```bash
   npm install && npm run build && grok publish <host>
   ```
   Pass `type: '<Package>:<funcName>'` to the *Filters* viewer
   (`DG-FACT-268`). For cross-package use, call
   `await grok.functions.call('<Package>:<funcName>')` FIRST so the
   function loads before the synchronous viewer call (`DG-FACT-269`).
   Prefer `addViewer(Viewer.filters(...))` over deprecated
   `view.filters(...)`:

   ```typescript
   await grok.functions.call('Widgets:radioButtonFilter');
   const tv = grok.shell.addTableView(grok.data.demo.demog());
   tv.addViewer(DG.Viewer.filters({
     filters: [{type: 'Widgets:radioButtonFilter', columnName: 'race'}],
   }));
   ```

## Common failure modes

- **Filter never appears in *Add filter…* picker.** Function role is
  wrong or missing — `src/package.g.ts` MUST contain
  `//meta.role: filter` AND `//output: filter result`, case-sensitive
  (`DG-FACT-266`).
- **Other filters' selections clear when this filter runs.**
  `applyFilter` is setting bits to TRUE; only flip failing rows to
  FALSE with `filter.set(i, <predicate>, false)` (`DG-FACT-270`).
- **UI changes don't refilter.** Change handler calls
  `this.applyFilter()` directly. Switch to
  `this.dataFrame!.rows.requestFilter()` (`DG-FACT-271`).
- **`Cannot read properties of undefined` in `attach`.** Subclass
  read `this.dataFrame` / `this.column` before `super.attach(...)` —
  call `super.attach` FIRST, then assign `this.column` /
  `this.columnName`, then build UI (`DG-FACT-272`).
- **Cross-package filter throws "function not found".** Add
  `await grok.functions.call('<Package>:<func>')` before
  `view.filters({...})` (`DG-FACT-269`).
- **State doesn't restore on layout reload.** Override `applyState`
  to call `super.applyState(state)` then re-render (`DG-FACT-273`);
  platform calls `applyState` AFTER `attach`.

## Verification

- `npm run build` + `grok publish <host>` exit `0`; `src/package.g.ts`
  contains a wrapper with `//meta.role: filter` and `//output: filter result`.
- In Datagrok: open *Filters* → *Add filter* — your widget appears and
  paints its UI; toggling it reduces visible rows and composes (rather
  than resets) other active filters.

## See also

- Source: `help/develop/how-to/viewers/custom-filters.md`; related
  `help/visualize/viewers/filters.md`; API `js-api/src/widgets/filter.ts`.
- Knowledge: `docs/_internal/knowledge/knowledge-graph.md` —
  `DG-FACT-265` … `DG-FACT-275`.
- Reference packages:
  `packages/Widgets/src/filters/radio-button-filter.ts:12-78` (single-choice);
  `packages/Widgets/src/filters/multi-value-filter.ts:6-93` (multi-value);
  `packages/Widgets/src/package.ts:14-33` + `package.g.ts:4-18`
  (decorator + auto-emitted header);
  `packages/ApiSamples/scripts/ui/viewers/filters/custom-filters.js`
  (cross-package invocation).
- Related skills: `custom-cell-renderers` (role-based registration).
