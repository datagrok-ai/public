# ems/dg — the platform layer

## Domain controls: `domains.*`

The platform adapter of the domain (EMS) seam (`src/dg/`; recipes: `docs/recipes/crud-app.md`
(zero code), `spec-app.md` (configuration), `custom-app.md` (code); what both backends must agree
on: `docs/domain-backend-contract.md`). `emsBackends.domain` is filled
at dg import with `DgDomainBackend` over `grok.dapi.domains`, so a `DomainSource` (core) built
anywhere in the platform loads real rows; the same source over `MemoryDomainBackend` — a frame
host too, over `MemoryFrame` — is what the headless tests and the gallery run.

| Export | What it does |
|---|---|
| `src/dg/styles.js` | a side-effect import: every stylesheet the domain stack paints with (tokens, the inputs and their editors, the list and the grid, the breadcrumb, the tabs, the section, the filter box, the dialogs and balloons). An app built on `domains.*` imports this one module instead of a hand-kept list of sheets; a package rendering other u2 controls adds those sheets itself |
| `domains.table<TRow>(address)` | one await → `DomainTable<TRow>`: `properties`, `info`, `access` (an `Access`), `handler` (the platform's `DomainObjectHandler` for the rows — never registered by u2), and the registries an app fills once: `actions` (`DomainAction {name, icon?, requires?, when?(row), run(row)}`), `validators.add(column, (value, row) => message \| null)`, `derived.add({from, to, derive(row)})` (columns the browser computes from other columns, waited for by Save and read-only in a form; a failure — the thrown message, or "did not finish in 60 s" — is shown on the field of its first input as soon as the run ends and refuses Save until the next good run; `deriveRows(rows)` runs them over plain rows), `operators.add({id, label, semType, parameter?, rows(condition, signal)})` (a filter operator the app answers — offered on columns of `semType`, sent to the server as the ids `rows` returns), `renderer` (an `ObjectRenderer<RowView<TRow>>`; the handler's rendering, the name column as caption — in its display `format`, so an autonumbered key reads `DG-000123`). `TRow` is the app's own row interface (`IssueRow` from the generated `db.ts`; any `{id: string}` qualifies): rows, `defaults`, `newRow` values and the callbacks of both registries are typed by it, `validators.add` offers its keys, and a `DomainTable<IssueRow>` goes wherever a `DomainTable` is expected — the controls take any. Default: `DomainRowLike`, every column `unknown` |
| `table.source(options)` / `table.draft(values)` | a started `DomainSource<TRow>` (`options.defaults`, `src.newRow(values)`: `Partial` of the typed row; `src.readScope()` what every read sends — the query with host conditions answered, the search, the trash mode); a `draft: true` source of one pristine draft. Both are tracked, so the controls find the handle (`DomainTable.of(src)`); `src.session` is the unit of work Save drives, `src.activate` the signal a list bumps on Enter |
| `table.app(options?)` | the table as a platform view (`grok.shell.addView(issues.app())`): a `DomainApp` under its own session in an `appView` — list ⇄ entity page, ribbon, status, URL, the unsaved gates. `options.app` names a `DomainApp` subclass to build instead (see below); `name`, `path`, and the app's own options pass through |
| `domains.app(options)` | the `DomainApp` control alone (`{table, base, query?, sort?, sorts?, pageSize?, mode?, filters?, include?, immutable?, sections?, children?, history?, shortcuts?, context?}` — `filters: {mode: 'builder', placement: 'left'}` puts the condition builder in a pane on the list's left) for a hand-built view with its own chrome — see "Domain apps" below for what it does and how it is subclassed |
| `domains.form(target, options?)` | `DomainForm` over a `DomainSource`, its `currentRow` signal, or a bare row (+ `source`): a fresh `propertyForm` per row under the row's access — the system columns left out (`system: 'footer'` shows them as a muted block: Id/Version/Created/Updated/Author, local times, the author resolved, a draft's id "assigned on save"), a draft's business key the caller may not type (`fields[col]` readonly — an autonumber) reading "assigned on save" until the insert lands, a readonly reference showing its target's name as a link that opens the row (`DomainPick.open`), `immutable: [columns]` — columns set on a new row only, an input on a draft and read-only once saved (a batch's compound), a pristine draft's verdicts held back until a field is edited or left (`validate()` shows them all), an edited field marked with an amber edge, `focus()` on the first field without scrolling the page; validity = schema rules + `table.validators` + a failed derivation from the field, shown as the derivation ends (`table.derivationProblem(src, id, column)`) + the writer's per-cell errors (over the platform the editor maps the server's column errors onto the cells and resolves a 409 itself). Paired through the source: it guards `src.save()` (a refusal is "Cannot save: <Caption> is required" whichever path ran — the button, Ctrl+S, `cmd:save`), takes the focus back on the session's `onSaved`/`onDiscarded` and on `src.activate`; a refused save is one balloon, a 403 also drops the access caches |
| `domains.list(src, options?)` | `DomainList`: a `VirtualList` over `src.rows`, `brief` or `cards` (the recipe card: title, one muted description line, the creation time; a handler's own card where it defines one), Open (saved rows) + Delete (`requires: 'delete'`) + the table's and the list's own actions per row; a row marked deleted stays struck through with Restore until the save; selection = `src.currentRow` (a loaded list starts on its first row), Enter bumps `src.activate` (the paired form takes the focus), Delete deletes, the selected row's actions rove into the tab order; the next page near the bottom, loading/empty/error under the rows |
| `domains.pick(table, options?)` | `DomainPick`, an `Input<string \| null>` face over a `TypeAhead` that queries the target table by its name column (else the first business-key column, else the id), with an open control beside a saved row; `PickInput` is the same face over any type-ahead — the `User`/`Group` editors ride it. `DomainPick.open(table, id)` (the open control, a readonly reference's link) brings forward the app open over that table, else opens the row in the view `domains.views` made for the table (`DomainViews.openRow`, which `table.open(row)` — a list's and a tree's Open — runs too), else the platform's entity view; it runs after the click, since the view a click lands in makes itself current as the click bubbles |
| `domains.grid(src, options?)` | `DomainGrid`: the platform grid over the source's frame with the source's writer attached (`Grid.attachEditor`), the table's handler decorating the columns; `options.look` binds `IGridSettings` |
| `domains.search(src, options?)` | `DomainSearch`, a ribbon search box two-way with `src.search` — the schema's searchable columns, AND-ed with the query; written after a pause (`debounceMs`, 300), at once on Enter, cleared by Escape; every search the user makes goes through `confirmDiscard` (H6) — refused while a batch is being written back, asked about while changes are pending |
| `domains.filters(src, options?)` | `DomainFilters`, the query box (`mode: 'query'`, completion over the table's columns and, over the platform, their values) or the condition builder (`'builder'`), two-way with `src.query`; every change the user makes goes through `confirmDiscard` first; text stays text, a tree stays a tree. The table's own columns are offered first and the system ones last, so a new condition starts on the first declared column; the values offered to pick hold no empty cell. The builder applies its complete rows; a row whose value is wrong (an emptied, mistyped or out-of-range parameter — not one still being written) holds the filter applied before it and is `problem`. A host operator is offered on the columns of a table that answers it — alone there, beside the null tests — and every other table keeps the core operators of the column |
| `domains.history(src, row?)` | `DomainHistory`, the entity-page pane over `table.audit(id)`: the current row's (or `row`'s) trail newest first, `caption: before → after` per updated column, system columns skipped, refreshed when the session saves; a draft says "Not saved yet" |
| `domains.children(src, options?)` | `DomainChildren`, the entity-page pane: one tab per table referring to the source's rows (`tables` narrows), a child `DomainSource` per tab in the parent's session with the FK defaulted to the parent's id (a draft parent → `~new:`, one transaction), `mode: 'grid'` (default) or `'list'` beside a form; rebuilt when the parent's id changes, never mid-save. New adds a row in the tab; `onNew(entry, values)` routes it elsewhere on a saved parent — the child app's create page, `DomainApp.createAt(base, values)` — while a draft parent's children stay in the tab |
| `domains.saveButton(src \| session, {blocked?})` / `domains.discardButton(src \| session)` / `domains.newButton(src, values?)` | buttons over a `DomainSession` (a source stands for its own): disabled while clean or saving, Save also while `blocked` (a signal) says why — its tooltip; an app's Save is blocked while its form shows a verdict other than "required", and its New while it is on a create page; the session balloons "<Singular> saved" once and emits `onSaved`, so no form is passed around; New adds a pristine draft (`src.newRow`), hidden without `insert` |
| `installEms(reg?)` | installs ems into a registry (the global one by default) as one u2 plugin (`emsPlugin`): `u2-domain-source`, the `u2-domain-*` controls — `form`, `list`, `pick`, `grid`, `data-table`, `tree`, `search`, `filters`, `history`, `children` — and `u2-manifest-editor`, plus the spec ambient that joins every source of one spec into one session. Once per registry, whoever calls it again; answers the uninstall. u2's `registerAll`/`registerPlatformComponents` register no domain tag, so a host rendering specs with domain tags calls it before the first render — the placeholder of a missing tag names it. A control binds to its source whole: `bind: {source: "$.issues.source"}`. `registerEmsComponents(reg)` is its deprecated alias; `registerDomainComponents(reg)` adds the controls alone |

The classes (`DomainForm`, `DomainList`, `DomainApp`, …) and the option types are exported beside
`domains` for subclassing and typing; the factories live on `domains` alone.

### Domain apps

Three tiers, each a complete app, each building on the one below (the recipes walk one reference
app per tier):

| Tier | What the author writes | Reference |
|---|---|---|
| Zero code | `schema.json` and `(await domains.table(address)).app()` — the app is the schema | Stockroom (`packages/Stockroom`), `docs/recipes/crud-app.md` |
| Configuration | a `dg-ui/1` spec: `u2-domain-source` on the tray + the `u2-domain-*` tags, mounted with `renderSpec` after `installEms()`, edited in the designer; phase 4 stores it on the table | Stockroom's `app.spec.json`, `docs/recipes/spec-app.md` |
| Code | the table's registries (`actions`, `validators`, `renderer`) and a `DomainApp` subclass passed as `app({app})` | Grit (`packages/Grit`) over the generated `getGritDb()`, `docs/recipes/custom-app.md` |

`DomainTable.app(options?: DomainAppViewOptions): DG.ViewBase` builds a `DomainApp` under its own
`SharedSession` (`SharedSession.runWith`) inside an `appView` and returns the view:

| Option | Default | What it does |
|---|---|---|
| `name` | `info.pluralName` | the view's name |
| `path` | `/domains/<schema>/<table>` | the base of the app's paths; `acceptsPath`/`handlePath` route deep links under it |
| `app` | `DomainApp` | the class to build — a subclass with its own ribbon, presets and shortcuts |
| `query`, `pageSize`, `mode` | `''`, 50, `'brief'` | the list page's source and row shape |
| `sort`, `sorts` | `sorts[0]`, none | the list's order (`'col,!col'`, the table's own without either) and the orders a ⋯ Sort submenu offers as `[label, sort]` — no submenu without them |
| `context` | the row's `DomainRow` | what the platform's context panel shows (`grok.shell.o`) for the row the user moves to in the list and for the entity page's row — again when a reference of it changes, so a create page shows it too; null shows nothing; the row the list opens on is not a move |
| `include` | every column | the entity form's columns, in this order |
| `sections` | none | `[{title, include, expanded?, hint?}]`: those fields under a heading after the others, in the same form (one label column); `expanded: false` folds it |
| `children` | `true` | the children pane: `false`, or `DomainChildrenOptions` (`tables`, `mode`, `onNew`) |
| `history` | `true` | the history pane |
| `shortcuts` | `{}` | key → action name, see below |

What the app does (`app.ts`): `page` (`'list' \| 'entity'`) and `entity` (the row id, `'new'` for a
draft) are read-only signals — `goTo(page, id?)` is their only writer and `open(path?)` the way in
from the address bar (`?entity=<id>`, `?q=<query>`; a path given is authoritative, `open()` with
none reads `location.search`, which is how a cold deep link lands); `path` mirrors them for `appView` and
`summary` is the session's while dirty, else the page's source's. The list page is a `DomainList`
over `listSource` (Enter → the entity page); the entity page a source of one under the same
session, `DomainForm` with `system: 'footer'`, then `panes`; `DomainTable.open(row)` activates a
live app (`DomainApp.activate(base, id)`) instead of the platform's row view. `create(values?)` opens
the create page over the app's `defaults` and `values`; `DomainApp.createAt(base, values)` does it on
the live app under `base`. The entity page's ⋯ Delete marks the row: a banner says it goes with the
next Save and offers Restore (the list's strike-through and Restore), and once the deletion is saved
the page goes back to the list. With `filters: {placement: 'left'}` the builder's pane is folded to a
strip while it holds no condition and opens when a filter is applied from the address or a preset.

Override points for a subclass (`class IssuesApp extends DomainApp`):

- `ribbon(): (Control \| HTMLElement)[][]` — `[[New, Save, Discard], [search, filters]]`, built
  once and owned by the app; extend it with `[...super.ribbon(), [...]]`. New opens `?entity=new`
  and is hidden without `insert`; search and filters show on the list page only.
- `presets(...entries: [label, query][]): ButtonGroup` — a single-toggle switch writing the list's
  query; `$me` binds to the current user's id (`Filters.bind`), the preset equal to the query in
  force is shown pressed, a press while dirty goes through the gate; list page only.
- `shortcuts: Record<string, string>` — `{'Ctrl+Shift+C': 'Close'}` (modifiers Ctrl, Alt, Shift in
  that order; Cmd counts as Ctrl) onto the names of `table.actions`, run over the current row of
  the page while the app has the focus, under the row's access; a subclass field or the option.

The gates — the same `confirmDiscard(session, {action})` dialog (Save / Discard / Cancel; a save
in flight → a warning and no move, checked before dirtiness, since a batch reads clean for the
length of its write-back) in front of every way to drop unsaved changes:

| Way out | Gate |
|---|---|
| a row, Back, `open(path)`, New | `goTo` → `confirmDiscard(session, {action: 'leave this page'})`; `open(path)` on the list page gates the `?q=` it applies too |
| the filter box, a preset | `DomainFilters` / `presets()` → `{action: 'change the filter'}`; the input is put back on cancel |
| the search box | `DomainSearch` → `{action: 'change the search'}`; the text is put back on cancel |
| the view's ✕ | `grok.events.onViewRemoving` — `preventDefault`, then `view.close()` once confirmed (`{action: 'close the view'}`) |
| the browser tab | `guardUnload()` — `beforeunload` armed while `session.isDirty` or `session.isSaving` |

Editor rules registered at import: a property whose `semType` is a table address
(`<schema>.<table>`) gets `domains.pick` over that table; `User` → the user picker, `Group` → the
group picker — everywhere `propertyForm` is generated, not only in domain forms.

### Access

`Access` (core) is the one protocol every control consults; the server answers it in exactly that
shape (`DomainTableClient.access()` → `{can, fields}`) and nothing asks the server twice.

- `can(capability)`: `view`, `insert`, `edit`, `delete`, `share` (plus a schema's custom names).
  `Action.requires` names one; a denied action is not rendered.
- `field(name)`: `hidden` where the server did not list the column, `editable` where it listed
  it so AND the row may be written, `readonly` otherwise. The write capability is `edit` for an
  existing row and `insert` for a draft (`access.forDraft()`) — a listed `editable` only means
  "not column-restricted", never "this caller may write this row".
- `row(r)`: the per-row view. A `withAccess` query adds `~can_edit`, `~can_delete`, `~can_share`
  to every row (`Access.ROW_COLUMNS`, the documented server set) — computed by the same predicates
  the server enforces on write; `row()` reads any `~can_<name>` the table has a capability for. A
  boolean the row carries REPLACES the table's answer for that capability either way: a row-mode
  table's table-level flags are false negatives, so a row may be editable under a table that says
  `edit: false`, and read-only under one that says `true`. A column that is absent or not a
  boolean is not carried — off row mode the server sends `~can_share` as null and omits it from
  a frame. `insert` is never per row.
- The controls pick the view for you: `domains.form` renders a draft under `forDraft()` and an
  existing row under `row(r)`; `domains.list` filters each row's actions under `row(r)` (a draft
  under the table's access — its `~can_*` cells are the frame's defaults, not truth). The frame
  path hides the three columns from every export the way the editor hides its own service columns.

### What the frame host honours (STATE-CONTRACT)

Over the platform backend a source is a frame with the js-api `DomainFrameEditor` attached as its
single writer (`EditorEditState`); over the memory backend a `MemoryFrame` with `MemoryEditState`.
Rows are keyed the way `FrameRows` keys them — the `id` cell, `Rows.draftKey(index)` (`~row:<index>`)
for a draft that has no id yet; `Rows` (core) is the one holder of the `~` conventions (`Rows.STATE`,
`isService`, `isDraft`). What that means for an app:

- **H6** — a programmatic re-query (a bound `query` changing) is skipped while the source is dirty;
  a user-initiated one must go through the unsaved-changes gate. A refresh drops the batch.
- **H7** — never export or hand off the source's frame: it carries the editing state and the
  access columns (export-tagged, but that is the safety net). Query a fresh frame.
- **H8** — a re-query replaces the frame and detaches the editor from the old one; the next page
  (`loadMore`) is appended into the same frame, pending edits kept.
- **H10** — everything mirrored (`isDirty`, `changeCount`, `validity`, `errorOf`) follows the
  editor's `onChanged`, never frame events; a direct `df.set` is invisible to save.
