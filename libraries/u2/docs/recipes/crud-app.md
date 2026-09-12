# Recipe: a CRUD app over a domain table

The phase-1 shape of an EMS app: one table handle, sources over it, the domain controls on top,
the shell's own chrome around them. Nothing here fetches, validates, tracks changes, checks a
permission or wires a button to a form by hand — the table's schema and the caller's access
drive all of it, and the controls pair through the source.

## The pieces

| Piece | What it is |
|---|---|
| `domains.table(address)` | one await: the table's schema, its display identity and the caller's `Access`; the per-table registries (`actions`, `validators`, `renderer`) live on the handle |
| `table.source({query, pageSize, defaults})` | a started `DomainSource`: `rows`, `currentRow`, `state`, `total`, `access`, `edit`, `isDirty`, `summary`, `session`; `save()` writes every pending row as ONE transaction through the session, `discard()` drops them |
| `table.draft(values)` | a `draft: true` source: loads nothing, holds one pristine draft — what a create form binds to; `save()` inserts it |
| `src.session` | the unit of work Save and Discard drive (`isDirty`, `isSaving`, `save`, `discard`, `onSaved`, `onDiscarded`) — a source's own session of one today; the phase-2 session over several tables has the same shape |
| `domainForm(src)` / `domainForm(src.currentRow)` | `propertyForm` over the current row: editors from the schema, text for what the caller may only read, writes through the source; it guards the source's Save (a refusal names the field), takes the focus back after a save, a discard, or a list's Enter |
| `domainList(src, {mode})` | a virtual list over `src.rows` with Open / Delete / the table's actions per row, selection = `src.currentRow`, Enter hands the row to the paired form, the next page on scroll |
| `domainPick(table)` | a type-ahead over another table's name column; the value is the row id — what a `ref` column holds (the form uses it by itself) |
| `saveButton(src)` / `discardButton(src)` / `newButton(src, values?)` | ribbon buttons over the source's session (a session may be passed instead): Save and Discard follow its state (disabled while clean or saving), Save announces "<Singular> saved"; New adds a pristine draft and is hidden without `insert` |

## The shape

```ts
import {Splitter} from '@datagrok-libraries/u2';
import {appView, domains, domainForm, domainList, newButton, saveButton, discardButton}
  from '@datagrok-libraries/u2/src/dg/index.js';

const issues = await domains.table('grit.issue');
issues.validators.add('title', (v) => String(v ?? '').trim().length < 5 ? 'At least 5 characters' : null);
issues.actions.add({name: 'Escalate', icon: 'arrow-up', requires: 'edit',
  when: (r) => r.priority !== 'p0', run: (r) => r.priority = 'p0'});

// list ⇄ form over one source; the view owns the source, the splitter owns the controls
const src = issues.source({query: 'status = "open"', pageSize: 50});
grok.shell.addView(appView({
  name: 'Issues',
  content: new Splitter([domainList(src, {mode: 'cards'}), domainForm(src)], {sizes: [40, 60]}),
  own: [src],
  ribbon: [[newButton(src, {status: 'open'}), saveButton(src), discardButton(src)]],
  status: src.summary,                                  // "50 of 2,077" · "3 unsaved changes"
}));

// a create form: a pristine draft, saved as an insert
const draft = issues.draft({reporter: grok.shell.user.id});
grok.shell.addView(appView({name: 'New issue', content: domainForm(draft), own: [draft],
  ribbon: [[saveButton(draft), discardButton(draft)]], status: draft.summary}));
```

The same page as a spec — a control binds to its source as a whole through the `source` step,
and `cmd:issues.save` goes through the source's session:

```json
{"$schema": "dg-ui/1",
 "components": [{"tag": "u2-domain-source", "name": "issues", "props": {"table": "grit.issue"}}],
 "root": {"tag": "u2-splitter", "children": [
   {"tag": "u2-domain-list", "bind": {"source": "$.issues.source"}, "props": {"mode": "cards"}},
   {"tag": "u2-domain-form", "bind": {"source": "$.issues.source"}},
   {"tag": "u2-button", "props": {"text": "Save"}, "on": {"click": "cmd:issues.save"}}]}}
```

## Rules the controls follow (so the app does not re-implement them)

- **Permission ⇒ hidden, state ⇒ disabled.** A column the caller may not see is not rendered;
  a row action the caller may not run is not rendered; Save is disabled while there is nothing
  to save. Nothing shows a 403 after a click.
- **Readonly is text.** A column the caller may not write — the table says so, or the row does
  (`~can_edit`), or a draft under a caller without `insert` — is caption + value, never a dead
  input, and never in the payload.
- **One writer.** Every edit goes through the source's edit state (the row proxy's `set`), so
  dirty state, validation and the transaction cannot disagree. A direct write to the frame is
  silently not saved.
- **One Save.** Every path — the button, Ctrl+S in the form, `cmd:issues.save` — is
  `src.session.save()`. The form registers a guard on the source, so a refusal is "Cannot save:
  Title is required" whichever path ran; the session announces a landed batch once and the form
  takes the focus back.
- **Validation is the union** of the schema (`nullable`, choices, ranges), `table.validators`
  and what the writer reports per cell. Over the platform the editor maps the server's column
  errors onto the cells and resolves a 409 through the platform's reload/overwrite dialog itself;
  a 403 re-gates the affordances.
- **Ownership follows the tree.** A `Control` handed to a container (`Splitter`, `appView`'s
  ribbon, toolbox and status) is disposed with it; a source is owned by the view through
  `appView({own})`. Nothing needs an `own(() => x.dispose())` line.

## Anti-patterns

- `grok.dapi.domains.table(...).update(...)` from a form's Save button — the source already
  batches every pending row into one transaction; call `src.save()`.
- A hand-written `if (!access.can.edit) input.enabled = false` — pass the access (the form
  already takes the source's) and let readonly render as text.
- Refreshing a source on a timer or a route change without looking at `isDirty` — a refresh
  drops the pending batch by design; the source skips a programmatic re-query while dirty, a
  user-initiated one goes through the unsaved-changes gate.
- Exporting or handing off the source's frame (`toCsv`, `batch`) — it carries editing state;
  query a fresh frame for that.
- Binding `$.issues` or `$.issues.currentRow` to a domain control — `$.issues` is the rows, and
  `currentRow` is a walkable step; the control needs the source: `$.issues.source`.
