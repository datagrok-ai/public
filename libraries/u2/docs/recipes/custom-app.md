# Recipe: a domain app in code (actions, validators, renderers, an app subclass)

The third tier of a domain (EMS) app: everything the [zero-code](crud-app.md) and
[spec](spec-app.md) tiers give, plus what only code can say — an action on a row, a rule across
columns, a card of your own, presets, shortcuts. App-specific code stays app-specific: one
registry each on the table handle, nothing re-implemented. Reference implementation: **Grit**
(`packages/Grit`), the issue tracker, over the typed handles `grok api --ui` generates.

## The shape

```ts
import * as grok from 'datagrok-api/grok';
import {badge, divV, span, timestamp} from '@datagrok-libraries/u2';
import type {Control} from '@datagrok-libraries/u2';
import {DomainApp} from '@datagrok-libraries/u2/src/dg/index.js';
import {gritDb} from './generated/db';                     // grok api: typed data clients
import {getGritDb} from './generated/db-ui';               // grok api --ui: typed u2 handles

const CLOSED = (await gritDb.statuses.getByKey({name: 'closed'}))?.id ?? null;   // a lookup row id, resolved once

class IssuesApp extends DomainApp {
  shortcuts = {'m': 'Assign to me', 'c': 'Close'};        // plain keys fire outside text entry; chords (Ctrl+Shift+C) everywhere
  ribbon(): (Control | HTMLElement)[][] {
    return [...super.ribbon(), [this.presets(['Mine', 'assignee = $me'], ['Open', 'status_id.name != "closed"']   // a FK path, no id needed)]];
  }
}

//name: Issues
//tags: app
//output: view result
export async function issuesApp(): Promise<DG.ViewBase> {
  const db = await getGritDb();                            // one await: schema, data clients, u2 handles
  const {issues} = db.tables;                              // DomainTable<IssueRow>
  const me = grok.shell.user.id;

  issues.actions.add({name: 'Assign to me', icon: 'user', requires: 'edit',
    when: (r) => r.assignee !== me, run: (r) => r.assignee = me});
  issues.actions.add({name: 'Close', icon: 'check', requires: 'edit',
    when: (r) => r.status_id !== CLOSED, run: (r) => r.status_id = CLOSED});
  issues.validators.add('status_id', (value, r) =>
    value === CLOSED && !r.assignee ? 'Assign the issue before closing it' : null);
  issues.renderer = {...issues.renderer, card: (r) => divV([
    span(r.title, 'u2-domain-card-title'), badge(r.status_id ?? 'new'), timestamp(r.created_on)])};

  return issues.app({app: IssuesApp, path: '/apps/Grit/Issues', children: {tables: ['comment'], mode: 'list'}});
}
```

Every callback is typed by `IssueRow`: `r.assignee` is `string | undefined`, `validators.add`
offers the row's columns, `run` writes through the row (the source's edit state — the change is
pending, Save lands it). The registries are filled once per handle; `issues.app()` and every
control over a source of `issues` read them.

## The pieces

| Piece | What it is |
|---|---|
| `getGritDb(): Promise<GritDb>` | generated (`grok api --ui`): `{schema, data, tables}` — `data` the typed dapi clients (`gritDb.issues.query(...)`), `tables` one `DomainTable<XRow>` per table, opened in parallel. Without the generator: `await domains.table<IssueRow>('grit.issue')` |
| `table.actions.add({name, icon?, requires?, enabled?, when?(row), run(row)})` | the ribbon, the row's hover actions, the context menu, the shortcuts. `requires` names a capability (`'edit'`, a schema `permissions` name such as `'approve'`) — denied, the action is not rendered; `when` narrows it per row; returns the unregister |
| `table.validators.add(column, (value, row) => message \| null)` | a rule the form checks on every change, beside the schema's and the server's; the message shows on the column's field and refuses Save |
| `table.renderer: ObjectRenderer<RowView<TRow>>` | `caption`, `icon?`, `listItem?`, `markup?`, `tooltip?`, `card?` — lists, cards, pickers, chips, breadcrumbs; spread the default to override one |
| `class X extends DomainApp` | `ribbon()` returns `[[New, Save, Discard], [search, filters]]` — extend it; `presets(...[label, query])` is a query switch (`$me` = the current user's id; the one matching the query in force is pressed; hidden off the list page; a press while dirty goes through the gate); `shortcuts` maps `'Ctrl+Shift+C'`-style keys onto action names, run over the current row while the app has the focus |
| `table.app({app: X, ...})` | the view: `name`, `path`, `query`, `pageSize`, `mode`, `include`, `children`, `history`, `shortcuts` — see [crud-app](crud-app.md) |
| `domains.app({table, base, ...})` | the `DomainApp` control alone, for a hand-built `appView` with its own chrome |

## Master–detail in one session

The entity page already is one: `domains.children(parent)` makes a child source per referring
table inside the parent's session, the FK defaulted to the parent's id. Written by hand — a
project form over a grid of its issues, one Save:

```ts
import {SharedSession, Splitter, computed} from '@datagrok-libraries/u2';
import {appView, domains} from '@datagrok-libraries/u2/src/dg/index.js';

const {projects, issues} = (await getGritDb()).tables;
const session = new SharedSession();                       // or build under SharedSession.runWith(session, …)
const project = projects.source({query: `id = "${id}"`, pageSize: 1, session});
const kids = issues.source({query: `project_id = "${id}"`, defaults: {project_id: id}, session});

grok.shell.addView(appView({
  name: 'Project',
  content: new Splitter([domains.form(project), domains.grid(kids)], {sizes: [40, 60]}),
  own: [project, kids],
  ribbon: [[domains.newButton(kids), domains.saveButton(session), domains.discardButton(session)]],
  status: computed(() => session.summary.value || kids.summary.value),
}));
```

A draft parent works the same way — its child holds the draft id, the server orders the inserts:

```ts
const draft = project.newRow({key: 'S7', name: 'Study 7'});   // id = '~new:…'
kids.newRow({project_id: draft.id, title: 'First issue'});     // '$~new:…' in the one transaction
await session.save();                                          // both tables; the child's FK re-points to the real id
```

A source made without `session` is a session of one; `app()` and `renderSpec` build under an
ambient session, so their sources join it without the option.

## Rules the code follows

- **Declare, don't wire.** An action is registered, not attached to a button: the ribbon, the
  row, the menu and the shortcut all read the registry, under the row's access.
- **Write through the row.** `run` mutates `r`; never `grok.dapi.domains.table(...).update` —
  the change joins the batch and shows as pending until Save.
- **A validator returns a message or null**, synchronously, from the row it is given; a rule the
  grammar can say (`end_date >= start_date`) belongs in `schema.json` `constraints` instead, so
  the server checks it too.
- **Presets are queries**, bound through `Filters.bind`; a preset the user cannot express in the
  filter box (`$me`) is exactly what they are for.

## Anti-patterns

- `class IssuesApp extends DomainApp { constructor() { … this.root.append(myToolbar) } }` — the
  ribbon is `ribbon()`; the shell owns the chrome.
- A validator that fetches (`await dapi…`) — validators are synchronous; a server-side rule is a
  constraint or a server-mapped column error.
- `issues.actions.add({name: 'Approve', run})` without `requires: 'approve'` — the permission is
  declared in the schema for exactly this; without it the action shows to everyone and 403s.
- A second `domains.table('grit.issue')` to reach the registries from another module — a handle
  per call, so the second one's registries are empty; hand the first handle around (the
  generated `db.tables` is that one place).
