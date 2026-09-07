# JS API conventions

How members of `datagrok-api` are shaped, named, documented, deprecated and tested. Distilled from the parts of the API that already work well (the audit behind it lives in the core repo under `core/docs/reviews/js-api-audit/`).

## 1. Shape of a member

- **One way to do a thing.** Before adding a method, search for the concept; extend the existing one with an
  options field rather than adding a sibling (`form`, not `narrowForm`/`wideForm`/`forms.condensed`).
- **Options object after two parameters.** `f(required, required, options?: IFooOptions)`. Never a third
  positional boolean. Options interfaces are exported, named `I<Method>Options` or `<Noun>Options` (one convention;
  the codebase currently has both — pick `I…Options` for new code, it is the majority).
- **Booleans are named.** `{recursive: true}` not `list(path, true)`.
- **Argument order is (row, column) everywhere** a cell is addressed (`cell(row, col)`); when a method must take
  (column, row) for a historical reason, its doc's first line says so.
- **Fluent builders are immutable** and return a new instance (`where().orderBy().page()`), or the class doc says
  "mutates and returns this" in its first line. Never mix within one class hierarchy.
- **Async methods are named like sync ones.** No `Ex`/`Sync`/`Async` suffix when only one exists. If both exist,
  the sync one is the plain name and the async one ends in `Async` (`toCsv` / `toCsvAsync`), never `Ex`.
- **Return the most specific wrapper type** (`DateTimeColumn`, `ScatterPlotViewer`, `Column<string>`), and `null`
  in the type when the runtime can return null. A `find*` method returns `T | null`; a `get*` method returns `T`
  or throws.
- **No `any` on the public surface.** Use generics with defaults (`Column<T = any>`), `unknown` for opaque payloads,
  `{[key: string]: string}` for tag bags, a callback type for callbacks. `Function`, `String`, `Object`, `Boolean`,
  `Number` are banned (lint).
- **Enums for closed sets, `(string & {})` for open ones.** `type ViewerType = \`${VIEWER}\` | (string & {})` keeps
  autocomplete and accepts plugin viewers.
- **Generics carry through.** If a class is generic, every factory that knows the type parameter states it
  (`Column.fromStrings(): Column<string>`), and every method that consumes the settings type uses it
  (`Viewer<TSettings>.setOptions(options: Partial<TSettings>)`).

## 2. Naming

- Classes: `PascalCase` nouns; viewer classes end in `Viewer` (`BoxPlotViewer`), data sources in `DataSource`,
  clients in `Client`, options in `Options`, event args in `Args`.
- Members: `camelCase`; getters are nouns (`rowCount`), predicates start with `is/has/can`, event streams start with
  `on` and are `Observable<TypedArgs>`, never `Observable<any>`.
- No abbreviations on the main path (`grok.shell.currentTable`, with `t` kept as a documented alias); no one-letter
  members in new code.
- Constants: one enum per concept. Subset enums (`CORE_VIEWER`, `STATS`) are replaced by a predicate or a
  `readonly` array on the owning class.
- Public members never start with `_`. Anything that must be reachable but is not API is `@internal` and lives in a
  file the barrels do not re-export.

## 3. Documentation template

Every exported class, function and public member gets a JSDoc block. Minimum: one sentence that says what it
*does or returns*, in the imperative or as a noun phrase. Then, in this order, only what applies:

```ts
/** Sorts rows by the specified columns.
 * Passing `orders` shorter than `columns` sorts the remaining columns ascending.
 * @param columns  Column names or objects; the first column is the primary key.
 * @param orders   `true` = ascending, `false` = descending; defaults to ascending.
 * @returns this grid (fluent).
 * @example
 * grid.sort(['age', 'sex'], [false, true]);
 * @see {@link https://public.datagrok.ai/js/samples/grid/order/order-rows} */
```

Rules:
- `/** */` only. `///` is a Dart habit and is invisible to TypeDoc, IDE hover and AI tooling. Trailing `// comment`
  on an interface field is also invisible — put it in a `/** */` above the field.
- No `@param {type}` braces — the type is in the signature. Name the parameter and say what a non-obvious value does.
- No `@constructs`, `@type`, `@returns {string}` boilerplate. No `@param dart - The underlying Dart object`.
- Contracts go in the first paragraph: nullability, mutation, sync/async, what throws, side effects
  (`Adds the view to the workspace and makes it current.`).
- One `@example` (2–5 lines, runnable in the platform console) on every class and every method with more than one
  parameter or a non-obvious result. Reuse the sample link when a sample exists instead of duplicating it.
- Sample links use one form, `@see {@link https://public.datagrok.ai/js/samples/<path>}`, where `<path>` is the
  file path under `ApiSamples/scripts` without extension. CI checks the link resolves (`data/check-sample-links.cjs`).
- A deprecated member has `@deprecated Use {@link Replacement} — <one reason>.` and nothing else changes. Prose
  markers ("Obsolete", "softly deprecated", "@Obsolete") are not allowed.
- Comments explain *why* only when the code cannot (`Row`'s Proxy note, the `.js` import suffix note in `base.ts`).
  They never narrate what the next line does, and never editorialise.
- Cross-cutting behaviour that a caller must know (`grok.functions.call` result shape, `HttpDataSource` statefulness,
  `FilterGroup.filters` mixed content) is documented **on the member**, not only in `CLAUDE.md`. `CLAUDE.md` links
  to the member.

## 4. Module layout

- One concept per file; a file over ~800 lines is split. Barrels (`dataframe.ts`, `widgets.ts`, `entities.ts`,
  `ui.ts`) re-export and hold no logic.
- `ui.ts` contains only DOM builders. `ObjectHandler` lives in `src/object-handler.ts`; inputs in `src/ui/inputs.ts`;
  the form auto-layout engine in `src/ui/form-layout.ts`; `css` enums in `src/ui/css.ts`.
- Interop plumbing (`wrappers*`, `proxies`, `utils_convert`) is imported by name where needed and marked
  `@internal`; `dg.ts` re-exports `toJs`/`toDart` explicitly and nothing else from those files.
- Test, benchmark and analytics helpers (`Utils.executeTests`, `Shell.reportTest`, `ClickUtils`, `Test`) move to
  `@datagrok-libraries/test` or an `@internal` `src/testing.ts`.
- No `declare let grok: any` / `DG: any` / `ui: any` in new code; import the module. Where a cycle forces a global,
  a one-line comment names the cycle.
- No import-time side effects other than registering the module. Prototype patches on DOM globals are removed or
  made explicit opt-in (`DG.installCanvasExtensions()`).

## 5. Deprecation protocol

1. Add the replacement in the same release, with a sample and a test.
2. Mark the old member `@deprecated Use {@link New}.` and make it delegate to the new one.
3. Add a `CHANGELOG.md` line under `v.next` that names both.
4. Run `rg` across `public/packages` and file one ticket per package that still uses the old member.
5. Remove in the next major (`versioning-policy.md`, rule 6). Keep a `MIGRATION.md` at the package root listing
   every removal since the previous major.

## 6. Tests and samples

- A new public member ships with an `ApiTests` test (category named after the class) and, if it has more than one
  parameter or returns a UI element, an `ApiSamples` sample linked from its JSDoc
  (`/grok-add-api-coverage` scaffolds both).
- Samples live under a stable path; moving one requires updating every `@see` that points at it (the link checker
  fails the build otherwise).
- `data/inventory.cjs` runs in CI and fails when the documented share of hand-written public members drops, or when a
  new export has no JSDoc.

## 7. AI-specific

- Prefer explicit option names that read as English (`createDefaultFilters: false`) over positional flags; AI code
  copies whatever the signature suggests.
- Put the one-sentence rule that prevents the common mistake in the first line of the doc: "Result is the single
  output value, or an object keyed by output name when the function declares several outputs."
- Keep `public/js-api/CLAUDE.md` as the *map* (where things live, the five traps) and the JSDoc as the *truth*; a
  fact that lives only in `CLAUDE.md` is a bug.
- Every namespace root (`grok`, `ui`, `DG`) exposes a discoverable index: `grok.ts` doc comments name each sub-object
  and its three most common calls.
