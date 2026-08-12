# @datagrok-libraries/domain-ui

Composable UI over [entity-mapped domain tables](https://datagrok.ai/help/develop/how-to/db/domain-schemas) —
the `grok.dapi.domains` data plane's UI counterpart.

Everything is reflective: components take columns, labels, choices, validation rules and
permissions from the runtime domain registry, so they work on any registered table with no
codegen and no hand-wired permission checks. Common operations are one or two lines.

## Install

```
npm install @datagrok-libraries/domain-ui
```

Requires `datagrok-api` ^1.27.8.

## The facade

The way in is the `domains` facade. `domains.table(...)` is the one await — it prefetches
the typed client, the registry metadata and the caller's capabilities — and every widget
factory on the resulting handle is synchronous (data loads inside the widget; await
`widget.ready` on the rare occasion you need the loaded frame):

```ts
import {domains} from '@datagrok-libraries/domain-ui';

const issues = await domains.table('grit.issue');
```

A whole browse/CRUD app is one expression — the list page (cards / brief / editable grid,
search, New, deep-linkable query), a page per row (form, child-table tabs, history,
permission-gated actions), URL round-trip, and the unsaved-changes gate on every
navigation path:

```ts
grok.shell.addView(issues.app());
```

Widget factories on the handle:

| Factory                    | Returns                | What you get                                                        |
|----------------------------|------------------------|----------------------------------------------------------------------|
| `form(options?)`           | `DomainForm`           | Property form of one row (new by default, `{row}` or `{id}` to edit): inputs generated from the registry, reference pickers with async typeahead, client + server validation mapped onto the named fields |
| `formDialog(options?)`     | `Promise<boolean>`     | The form in a dialog — OK saves (a failed save keeps it open), Cancel/Esc discards silently |
| `grid(options?)`           | `DomainGrid`           | Editable grid: batch editing with amber/red/conflict markers, one-transaction save with the 409 flows, per-column writability gating, `defaults` prefilling added rows |
| `list(options?)`           | `EntityListWidget`     | The row list (cards / brief / grid, search, New, per-item commands) as an embeddable widget |
| `listView(options?)` / `app(options?)` | `DomainAppView` | The list as a page: ribbon, `DomainQuery` ⇄ URL round trip, row pages, the gate |
| `entityView(rowOrId)`      | `DomainEntityAppView`  | The page of one row                                                  |

Data-plane passthroughs (`query`, `queryDf`, `get`, `count`) ride along, so a handle is
also enough to read back what a form just wrote. Capabilities are a snapshot taken at
acquisition — re-acquire the handle to re-gate after permission changes.

An app that works with several tables acquires the whole schema at once. Each table is a
camelCase plural property (`db.issues`, `db.issueLabels`); declared singular names remain
in `db.tables` and `db.table(name)`:

```ts
const db = await domains.db('grit');
grok.shell.addView(db.issues.app());
const saved = await db.projects.formDialog({values: {key: 'K', name: 'N'}});
```

With `grok api --ui`, the generated `db-ui.ts` types the same handles per table —
`const db = await gritUiDb(); db.issues.form({values: {priorty: 'x'}})` fails to compile.

## Composed pages and dialogs

`domains.view(widgets)` hosts any `DG.Widget`s as one page with one aggregated dirty
state, one status bar, one unsaved-changes gate, and the children's actions merged onto
the ribbon. `formView(w)` is the one-widget shorthand, `dialog(w)` the dialog host:

```ts
grok.shell.addView(domains.view([
  projects.form({row}),
  issues.grid({query: {filter: `project_id = "${row.id}"`}, defaults: {project_id: row.id}}),
], {name: row.displayName}));

grok.shell.preview = domains.formView(issues.form({values: {reporter: me.id}}));
```

A `DG.Viewer` is already a widget, so a chart composes directly and stays live on the
grid's DataFrame — in-grid edits reflect before saving.

One-shot openers stay string-addressable: `domains.pick(table, {anchor?})` (the lookup
picker — a dialog, or a drop-down below `anchor`), `domains.create(table)`,
`domains.edit(row)`, plus `conflictDialog`, `auditPane`, `grantsPane`, `open(path)`, and
`deepLink(row)`.

## Customizing

Override one thing and stay in the few-lines world:

```ts
// One form field, one validation:
const form = issues.form();
form.replaceInput('description', (p) => ui.input.textArea(p.caption));
form.addValidator('title', (s) => s.trim().length < 5 ? 'At least 5 characters' : null);

// One action, appended to (or filtering) the permission-gated defaults:
issues.listView({actions: (row, defaults) => [...defaults,
  {name: 'Escalate', icon: 'exclamation-triangle', run: () => escalate(row)}]});

// One renderer — the standard handler mechanism; every list, grid, picker
// and view of the table picks it up:
DG.ObjectHandler.register(new class extends DG.DomainObjectHandler {
  constructor() { super('grit.issue'); }
  renderCard(x) { return ui.card(ui.divText(`#${x.values.number} — ${x.values.title}`)); }
});
```

## The machine surface

Every widget reports itself through the platform's own protocol — no DOM scraping in
tests, and AI assistants drive the same members a keystroke does:

```ts
const form = issues.form();
form.props['title'] = 'Login page 500s';            // same path as typing
form.getWidgetStatus().inputs;                      // values, choices, validation, per input
await form.getFunctions().find((f) => f.name === 'Save')!.apply({widget: form});
```

Actions are real registered platform Funcs (`Save`, `Discard`, `AddRow`, ...), returned
by `getFunctions()` filtered to what the caller's permissions allow.

## The state model

The editing state machine ships in `datagrok-api` as `DG.DomainFrameEditor` (with
`DG.DomainGrid` as its platform grid host); this library re-exports the same classes and
builds the forms, lists and pages on top. The editor is the single writer of three
service columns on the frame:

| Column     | Content                                             |
|------------|-----------------------------------------------------|
| `~state`   | `'' \| 'new' \| 'modified' \| 'deleted'` per row    |
| `~changes` | JSON — the original value of every changed cell     |
| `~errors`  | JSON — per-cell validation / conflict messages      |

Each is tagged out of binary and CSV export, so the editing state is memory-only: a saved
project, `toCsv()`, an export or a `batch()` upload built from the frame never carry it.
Write through `setValue()` (programmatic) or `beginEdit()` + `commitEdit()` (what a grid
does) — writing a cell directly on the DataFrame bypasses tracking and is silently not
saved. Save writes the whole batch as one `/transaction`: audit rows share a `tx_id`, and
a version conflict goes through the platform's standard reload/overwrite dialog.

`refresh()` re-runs the query and rebuilds the frame from scratch — there is no merge.
Deciding whether it is safe to refresh is the caller's job (`isDirty` + the shared
`confirmDiscardChanges` prompt); the page and app classes above already do this at every
navigation path.

## Permissions

Every affordance comes from `DG.DomainTableCapabilities`: no `canEdit` gives a read-only
grid, no `canInsert` removes New/Add row, no `canDelete` removes Delete, and columns
outside `writableColumns` stay read-only and never appear in a payload.

## Learn more

* [Domain schemas](https://datagrok.ai/help/develop/how-to/db/domain-schemas) — declaring
  schemas, security modes, the data-plane API
* [Domains](https://datagrok.ai/help/govern/catalog/domains) — the user-facing guide
* [Grit](https://github.com/datagrok-ai/public/tree/master/packages/Grit) — an issue
  tracker built on this library
