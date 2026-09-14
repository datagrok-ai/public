# Recipe: a CRUD app over a domain table (zero code)

The first of the three tiers of a domain (EMS) app: the app is the schema. `schema.json` declares
the tables; `domains.table(address).app()` turns one of them into a platform view — list ⇄ entity
page, ribbon, status bar, URL, children, history and every unsaved-changes gate. Nothing fetches,
validates, tracks a change or checks a permission by hand. Reference implementation: **Stockroom**
(`packages/Stockroom`), a chemical stockroom on the GHS taxonomy; the whole app is:

```ts
import {domains} from '@datagrok-libraries/u2/src/dg/index.js';

//name: Stockroom
//tags: app
//output: view result
export async function stockroomApp(): Promise<DG.ViewBase> {
  return (await domains.table('stockroom.substance')).app();
}
```

The other two tiers build on the same handle: [spec-app](spec-app.md) lays the page out as a
designer-editable spec, [custom-app](custom-app.md) adds actions, validators, renderers and an app
subclass in code.

## What the schema declares → what the app shows

Every key below is one Stockroom declares (`databases/stockroom/schema.json`); nothing else drives
the app.

| `schema.json` | In the app |
|---|---|
| `friendlyName: "Substances"` (singular/plural derived) | the view's name, the breadcrumb, "New substance", "Substance saved", "No substances." |
| `isName: true` on `name` | the row's caption everywhere — list rows, cards, breadcrumbs, pickers, the "kept X" hint |
| `searchable: true` on `name`, `cas` | the ribbon search box: ILIKE over exactly these columns, AND-ed with the filter; a paged count that agrees |
| `filters: [{column: 'hazards.code'}, …]` | the platform's Domain View filter panel; the u2 ribbon filter box needs no declaration — it completes every column, and over the platform the values behind them |
| `required`, `choices`, `min`/`max`, `type` | the form's editors and its first-line validation: "Value can't be empty", a choice input, a bounded number, a date picker |
| `ref: "location"` (+ `filter: "site = $site"`) | a type-ahead picker over the target's name column; with `filter`, narrowed by the sibling column — "Pick a site first" while it is empty |
| `type: "user"` (`owner`, `approved_by`) | the platform's user picker; shown as the user, stored as the id |
| `editor: "textarea"`, `semType: "Molecule"` | a multi-line editor; the platform's molecule renderer where the Chem package is installed |
| `constraints: {dates: {expr: "expires >= received", message: …}}` | pre-validated in the form on every change — the message on the first named field; the same CHECK on the server |
| `autoNumber` on `label`, `number` | read-only, "assigned on save" on a draft, filled by the engine in the saved row |
| a table referring to this one (`container.substance_id`, `sds_document.substance_id`) | one tab per referring table under the form, its rows queried by the FK, New pre-filled with the parent's id, saved with the parent as one transaction |
| `audit` (on by default) | the History pane: who did what when, `caption: before → after` per updated column |
| `permissions: {approve: …}`, the five built-ins | `~can_<name>` per row; a control the caller may not use is not rendered (the code tier binds actions to a custom name) |
| `securityMode: "master"`, `delegate` | the child's access is its parent's — the tabs under a substance edit under the substance's rights |

## What `app()` gives

`table.app(options?)` returns a `DG.ViewBase` for `grok.shell.addView`. Inside it a `DomainApp`
under its own `SharedSession`:

- **Two pages, one session.** The list page (`domains.list` over a paged source; Enter or Open →
  the entity page) and the entity page (breadcrumbs, `domains.form` with the system columns as a
  footer, the children tabs, the history pane). A New in the ribbon opens the entity page over a
  draft; once saved, the page is the row's.
- **The URL follows.** `base?q=<query>` on the list page, `base?entity=<id>` (`?entity=new` for a
  draft) on the entity page; `base` is `/domains/<schema>/<table>` unless `path` says otherwise.
  Deep links land through `acceptsPath`/`handlePath`; a saved row opened from anywhere
  (`DomainTable.open(row)`) activates the running app instead of a second view.
- **The ribbon** is `[[New, Save, Discard], [search, filters]]`: New hidden without `insert`, Save
  (Ctrl+S; Tab from the form's last field lands on it) and Discard disabled while there is nothing
  to save, search and filters shown on the list page only.
- **The status bar** is the session's summary while it is dirty ("3 unsaved changes in 2 tables"),
  else the page's ("50 of 2,077", "New substance", "Substance saved").
- **Three gates in front of dropping unsaved changes**, each the same Save / Discard / Cancel
  dialog: in-app navigation (a row, Back, a filter or preset change), the view's ✕
  (`onViewRemoving`), the browser's unload (`beforeunload`, armed exactly while dirty).

| `app()` option | What it does |
|---|---|
| `name`, `path` | the view's name (the plural name by default) and its base path |
| `query`, `pageSize`, `mode: 'brief' \| 'cards'` | the list page's initial query, page size and row shape |
| `include: string[]` | the form's columns, in this order |
| `children: false \| {tables, mode: 'grid' \| 'list'}`, `history: false` | the entity-page panes (both on by default) |
| `shortcuts`, `app` | keys onto the table's actions, and a `DomainApp` subclass — the [code tier](custom-app.md) |

## Rules the controls follow (so the app does not re-implement them)

- **Permission ⇒ hidden, state ⇒ disabled.** A column the caller may not see is absent; a row
  action they may not run is absent; Save is disabled while there is nothing to save.
- **Readonly is text.** A column the caller may not write — the table says so, the row does
  (`~can_edit`), or a draft under a caller without `insert` — is caption + value, never a dead
  input, and never in the payload.
- **One writer, one Save.** Every edit goes through the source's edit state; every Save path —
  the button, Ctrl+S, the gate's SAVE — is `session.save()`: the dirty sources plus every source
  holding a draft they refer to, as one transaction. A refusal is one balloon ("Cannot save:
  Name is required"); the server's column errors land on the cells; a 409 is resolved through
  the platform's reload/overwrite dialog.
- **Ownership follows the tree.** The view owns the app, the app owns its sources and controls;
  closing the view disposes everything.

## When to reach for configuration or code

- The page needs a different layout, a second table side by side, or an admin should edit it in
  the designer → [spec-app](spec-app.md).
- A row needs an action ("Approve", "Assign to me"), a cross-column rule the grammar cannot say,
  a card of its own, presets ("Mine"), keyboard shortcuts → [custom-app](custom-app.md).

## Anti-patterns

- A `grok.dapi.domains.table(...).update(...)` behind a Save button — `app()` already saves every
  pending row as one transaction.
- Re-declaring in code what the schema says (`if (!access.can.edit) input.enabled = false`, a
  hand-written "required" check) — declare it in `schema.json` and every tier gets it.
- A second view over the same table for "create" — New is the entity page over a draft.
