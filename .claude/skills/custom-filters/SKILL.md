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
---

## Cited facts

See [`facts.yaml`](./facts.yaml) — concrete API references for the `DG-FACT-NNN` citations used below.

# custom-filters

## When to use

A built-in slicer (categorical, histogram, search, sketcher) cannot
express your column's predicate — single-choice radio, multi-value
checkbox list, substructure / similarity search, custom date range —
and the widget must compose with every other selector in *Filters*.

## Prerequisites

- `datagrok-tools ^4.12.x` for the `@grok.decorators.filter` shortcut
  (`DG-FACT-267`); older toolchains use the function-role surface.

## Steps

1. **Subclass `DG.Filter`; implement the two abstract members.**
   `filterSummary` and `applyFilter()` MUST be overridden — see
   `DG-FACT-265` for constructor/super contract.

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

2. **Register the class as a package function with `role: filter`.**
   Platform discovers filters by function role (`DG-FACT-266` — needs
   `meta.role: filter` AND `filter`-typed output). Class-decorator
   shortcut `@grok.decorators.filter({name?, description?, semType?})`
   applies both implicitly — no `package.ts` factory needed
   (`DG-FACT-267`).

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

3. **Filter additively in `applyFilter`; never set bits TRUE.**
   See `DG-FACT-270` — only flip rows to FALSE via `filter.set(i,
   <pass>, false)`, finish with `filter.fireChanged()`.

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
   See `DG-FACT-271`. Optionally chain `get isFiltering()` through
   `super.isFiltering` to honor user-disable (`DG-FACT-274`).

   ```typescript
   radio.on('change', () => this.dataFrame!.rows.requestFilter());
   ```

5. **Build, publish, attach to a TableView.**
   ```bash
   npm install && npm run build && grok publish <host>
   ```
   Pass `type: '<Package>:<funcName>'` to the *Filters* viewer
   (`DG-FACT-268`). For cross-package use, `await grok.functions.call(
   '<Package>:<funcName>')` first so the function loads before the
   synchronous viewer call (`DG-FACT-269`). Prefer
   `addViewer(Viewer.filters(...))` over deprecated `view.filters(...)`.

   ```typescript
   await grok.functions.call('Widgets:radioButtonFilter');
   const tv = grok.shell.addTableView(grok.data.demo.demog());
   tv.addViewer(DG.Viewer.filters({
     filters: [{type: 'Widgets:radioButtonFilter', columnName: 'race'}],
   }));
   ```

## Common failure modes

- Filter missing from *Add filter…* — check `package.g.ts` has `//meta.role: filter` AND `//output: filter result` (`DG-FACT-266`).
- Other filters' selections wiped — `applyFilter` is setting TRUE; only flip to FALSE (`DG-FACT-270`).
- UI changes don't refilter — replace `this.applyFilter()` call with `rows.requestFilter()` (`DG-FACT-271`).
- `Cannot read properties of undefined` in `attach` — call `super.attach` FIRST (`DG-FACT-272`).
- Cross-package "function not found" — `await grok.functions.call(...)` before `view.filters(...)` (`DG-FACT-269`).
- State doesn't restore on layout reload — override `applyState` (`DG-FACT-273`).

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
