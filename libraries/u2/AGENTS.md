# u2 for coding agents

The short version of `@datagrok-libraries/u2`. Behind it: [README.md](README.md) (each rule with its how and what
catches a slip), [docs/recipes/](docs/recipes/README.md), [docs/samples/](docs/samples/README.md), `manifest.json`
and [src/dg/README.md](src/dg/README.md).

## Entry points

| Import | What it holds |
|---|---|
| `@datagrok-libraries/u2` | the platform-free library: signals, `component`, inputs, forms, dialogs, menus, collections, `Format` |
| `@datagrok-libraries/u2/dg` | the platform layer: `ViewHost`, `host`, `columnInput`, `tableInput`, `propertyForm`, `Platform` |
| `@datagrok-libraries/u2/filter`, `@datagrok-libraries/u2/dg/filter` | the filter model and its platform side |
| `@datagrok-libraries/u2/testing` | the testing kit |

Never a path inside the package, and no CSS to import. `datagrok-api` is a peer dependency: the signal core is
`datagrok-api/u2core`, which the platform serves once per page as `DG.U2`.

## Ten ideas

Everything u2 builds belongs to an owner and goes when the owner goes:

1. `component((self) => { …; return {root, …api}; }, {name})`: a piece of UI with a life of its own; the builder runs
   with `self` as the owner, and what it answers becomes the control and its handle.
2. `owner.build(fn)`: runs `fn` with that owner building (`self`, a host), so what it makes is the owner's.
3. `owner.task(async (t) => …)`: async work for an owner; after an `await`, get back in with `t.build(…)`.
4. `owner.subscribe(source, handler)`: a subscription that ends with the owner.
5. `effect(fn)`: runs again when the signals it reads change, until its owner goes.
6. `owner.own(cleanup)`: a cleanup the owner runs when it goes.
7. `on(el, type, handler)`: a listener that goes with whoever built the element.
8. `Scope.root(label)`: a tree nothing else owns, started on purpose and disposed by whoever started it.
9. `ViewHost`: a platform view as an owner, `ViewHost.of(view)` for an open view; an app opens a new one with
   `grok.shell.addView(ViewHost.create({name, content: () => myApp()}).view)`, never a control's `.root` in a view.
10. `host(control)`: chrome that is not a view (a context-panel widget, a pane); its close disposes the control.

A value is `input.value` (followed inside an effect), `input.peek()` (not followed) or `input.signal` (to wire). An
off is `disabled`: the text the user reads, a signal of it, or `null` while on; never a boolean.

## Rules

<!-- rules: generated from src/core/rules.ts by `npm run manifest`; edit them there -->
- **Import from the entry points only.** A path inside the package is no API: its files move and change between
  releases. Write `import {TextInput} from '@datagrok-libraries/u2';`, not
  `import {TextInput} from '@datagrok-libraries/u2/dist/src/index.js';`.
- **Build inside an owner, as a function.** What an owner builds is disposed with it; what nothing owns outlives the
  screen that showed it. Write
  `component((self) => { const nameInput = new TextInput({label: 'Name'}); return {root: [nameInput], nameInput}; })`,
  not `const nameInput = new TextInput({label: 'Name'});  // at module load or after an await: nothing owns it`.
- **An input's value is `input.value`.** A slot u2 keeps live needs the signal: the plain value is read once there
  and never changes. Write `span(nameInput.signal)`, not `span(nameInput.value)`.
- **Register on the thing that lives as long.** An effect put on a long-lived control from a row is added again by
  every rebuild of the row. Write `rowInput.effect(() => paint(row));  // on a control built in that row`, not
  `self.effect(() => paint(row));  // inside an each() row`.
- **Act on the user with `onInput`.** `onChanged` also hears writes from code, so the work runs again whenever code
  sets the value. Write `new ChoiceInput({label: 'Library', items, onInput: (name) => load(name)})`, not
  `new ChoiceInput({label: 'Library', items, onChanged: (name) => load(name)})`.
- **Load a record with `reset(value)`; follow a model with `follow(read)`.** A plain write fires `onChanged` and
  keeps the touched state of the record before, so a loaded record looks edited. Write
  `nameInput.reset(record.name);`, not `nameInput.value = record.name;  // loading a record`.
- **Save goes by `validity`; what a field shows goes by `error`.** `error` stays empty until a field is edited, left
  or validated, so a Save that follows it is on for an empty form. Write
  `button('Save', save, {disabled: form.validity})`, not `button('Save', save, {disabled: nameInput.error})`.
- **A tab marks itself: give it a `validity`, or make its content a `Form`.** A marker set by hand drifts from the
  fields it stands for. Write
  `{id: 'fields', title: 'Fields', content: () => fieldsForm, validity: fieldsForm.validity}`, not
  `tabs.setTabTitle('fields', 'Fields (!)');  // marked by hand`.
- **Say why something is disabled.** A control that is off with no reason leaves the user guessing what to do. Write
  `button('Save', save, {disabled: 'Pick a file first'})`, not `button('Save', save, {disabled: true})`.
- **An action's `run` returns a promise or nothing.** u2 waits for what `run` answers: a value handed back by
  mistake, such as a `DG.Row`, is no work to wait for. Write `{text: 'Add row', run: () => { df.rows.addNew(); }}`,
  not `{text: 'Add row', run: () => df.rows.addNew()}`.
- **Forms that must line up share one `LabelColumn`.** Each form aligns its own captions, so two forms one above the
  other start their editors at different places. Write
  `const labelColumn = new LabelColumn(); const tabs = [new Form({labelColumn}), new Form({labelColumn})];`, not
  `[new Form(), new Form()]  // two tabs of one editor`.
- **Widths are options, not CSS.** A width rule over u2's classes fights u2's own sizing (hence the `!important`) and
  breaks when they change. Write `new TextInput({label: 'Name', fill: true})`, not
  `.my-panel .u2-input-editor { width: 100% !important; }`.
- **Lay a view out with `Page`; a list beside its details is a `MasterDetail`, past a few hundred items a
  `VirtualList`.** A layout written in package CSS breaks as the view narrows, and drifts from every other view's.
  Write `new Page({content: form, footer: [statusLine(form.shownProblem), saveButton], maxWidth: 820})`, not
  `.my-editor { display: flex; flex-direction: column; height: 100%; }  /* and a pinned footer by hand */`.
- **Style your own classes, never u2's.** u2's classes and attributes are its internals: they change without notice,
  and a rule over them fights u2's own. Write `.my-app { --u2-row-hover-bg: var(--blue-1); }`, not
  `.my-app .u2-table-row:hover { background: #eef; }`.
- **A value shown as an element is a `DynamicInput` row; the rail takes buttons.** The rail beside an editor is one
  row high: a drawing put there spills over the rows under it. Write
  `new DynamicInput({label: 'Structure', value, render: (s) => s === null ? null : grok.chem.drawMolecule(s, 160, 100)})`,
  not `smiles.addOptions(grok.chem.drawMolecule(value, 160, 100))`.
- **Listen with `on(el, type, handler)`, follow with `effect(fn)`, on the elements a row builds.** Registered on the
  long-lived control instead, the handlers and the element stay with it after the row is gone. Write
  `on(card, 'click', open);`, not `self.listen(card, 'click', open);  // inside an each() row`.
- **A docked editor panel is laid out by its form: columns, a status line and folds.** A panel laid out by package
  CSS breaks when the dock narrows it. Write `form.addColumns([symbol, name], {minWidth: 160})`, not
  `.my-editor .u2-form { display: grid; grid-template-columns: 1fr 1fr; }`.
- **Dock with `ViewHost.of(view).dock(control, place)`.** Docked by hand, a control outlives its panel and its view.
  Write `ViewHost.of(view).dock(panel, 'right', {title: 'Details'})`, not
  `grok.shell.dockManager.dock(panel.root, 'right')`.
- **Find a view's control from its root.** The platform hands out a new wrapper object for a view each time, so
  nothing stored on one is there the next time. Write `Control.forElement(view.root.querySelector('.my-view'))`, not
  `(view as any).editor = editor;  // read back from grok.shell.v later`.
- **A dialog comes with CANCEL; one whose work applies as it is done is `closeOnly`.** A dialog that applies its work
  as it goes has nothing to cancel: one CLOSE says so, and ✕ and Escape agree. Write
  `new Dialog('Progress', {closeOnly: true})`, not `dialog.getButton('CANCEL')?.remove();`.
- **Personal settings apply at once; shared and admin settings stage behind SAVE N CHANGES.** The user's own
  preference should take effect as it is set, while a change other people feel waits until it is saved on purpose.
  Write `propertyForm(props, settings, {onChanged: () => save(settings)})  // the user's own settings`, not
  `form.addButtons((row) => row.append(button('Save', save)))  // the user's own settings`.
- **A dialog opened from a row outlives the row.** A rebuild of the row does not close the dialog, so content the row
  built would lose its handlers under it. Write `dialog.add(() => detailsOf(row));`, not
  `dialog.add(detailsOf(row));  // built by the row`.
- **A balloon's content is the balloon's.** Content built in a row's handler keeps its handlers on the row, so its
  Undo dies with a rebuild. Write `notify.info(() => div(['Deleted. ', link('Undo', undo)]), {onClose})`, not
  `notify.info(div(['Deleted. ', link('Undo', undo)]))`.
- **An app home starts with `AppHeader`.** The same header in every app tells the user what the app does; one built
  by hand drifts from the rest. Write
  `ViewHost.create({name: 'My App', content: () => [header(), campaigns()]})  // header() builds the AppHeader`, not
  `div([h1('My App'), span('Design molecules with your team')])`.
- **Own a platform subscription with the owner's `subscribe`.** A raw subscription is never ended: its handler keeps
  running after the view or the control is gone. Write `self.subscribe(grok.events.onViewRemoved, handler)`, not
  `grok.events.onViewRemoved.subscribe(handler)`.
- **Code living in a platform view is owned by its `ViewHost`.** The view is what the user closes: what its host owns
  ends then, whichever wrapper of the view the code held. Write `ViewHost.of(tv).build(() => detailsPanel())`, not
  `grok.events.onViewRemoved.subscribe((v) => { if (v.name === tv.name) panel.dispose(); })`.
- **Ask with `confirmDialog()` before what cannot be undone.** The platform asks the same way: CANCEL first and the
  action named as a verb, busy until the work is done. Write
  `confirmDialog({title: 'Delete campaign', message: 'Delete KHL-3?', confirmText: 'Delete', onConfirm: () => remove()})`,
  not `if (window.confirm('Delete?')) remove();`.
- **Rows the user edits are a `RecordsInput`; a pick with an order is a `ChecklistInput`.** Both are inputs: they
  validate, mark their form and read back as one value. Write
  `new ChecklistInput({label: 'Run', items: functionNames, orderable: true})`, not
  `new MultiChoiceInput({label: 'Run', items: functionNames})  // the order is lost`.
- **A list of records is a `CatalogTable`.** It owns the order, the groups and the chosen columns, and keeps them
  where the app keeps them. Write
  `new CatalogTable({items, keyOf: (c) => c.code, columns, groupings, persist: 'my-app.campaigns'})`, not
  `ui.table(campaigns, (c) => [c.name, c.status])  // header listeners, rows put in by hand`.
- **A column choice is `columnInput`.** A list built from `df.columns` offers the hidden `~` columns and nameless
  ones the platform's pickers leave out. Write
  `columnInput('Column', df, {filter: (c) => c.semType === DG.SEMTYPE.MOLECULE})`, not
  `new ChoiceInput({label: 'Column', items: df.columns.names()})`.
- **A table choice is `tableInput`.** A list of the shell's tables offers a table with no name, which no picker can
  hold. Write `tableInput('Table')`, not
  `new ChoiceInput({label: 'Table', items: grok.shell.tables.map((t) => t.name)})`.
- **One text for a value: `Format`.** A value reads the same everywhere u2 shows it: a date in the user's platform
  format, in one zone. Write `span(Format.value(row.created))`, not `span(String(row.created))`.
- **Comment only what the code cannot say.** A narrative goes stale the day the code changes, and the code, the
  README and the history already say it. Write
  `// a rename fires onColumnNameChanged only, so the held name is remapped here`, not
  `/* The campaign editor: fields, stages and the save flow, reworked for the new templates. */`.
<!-- /rules -->

## Choosing a component

`manifest.json` has each component's `usage` and each kind's `whenToUse` and `notFor`; the tables of the recipes
[forms-and-details](docs/recipes/forms-and-details.md), [catalog-table](docs/recipes/catalog-table.md) and
[shell-integration](docs/recipes/shell-integration.md) pick between neighbours.

## Runnable samples

[docs/samples/](docs/samples/README.md): whole, copy-ready code for the everyday jobs, one `<id>.ts` each.
`tests/recipes.test.js` compiles every sample and every TypeScript block of the README and the recipes against u2's
sources, and runs the live ones under the testing kit.

## Testing

`@datagrok-libraries/u2/testing`, the same in `node:test` and `@datagrok-libraries/test` (an example:
[README › Testing](README.md#testing)): `u2test(body)` fails a test that leaves a scope, a timer or a listener live,
raises a usage error or lets an error reach `Scope.onUnhandledError`; `mount(build)` builds under the test;
`user.click`, `user.type`, `user.press` and `user.pick` act as a user does and refuse what a user could not do;
`expectOperable(el)` lists what is offered but cannot be operated; `byName(within, name)`, `byKind(within, kind)`
and `state(el)` find controls and read their state.

## Checks

From `public/`:

- `node libraries/u2/scripts/check-consumer.mjs packages/<Package>`: the u2 lint rules and the package's counts
  against `libraries/u2/consumers-baseline.json` (`grok check` runs the same; `--all` checks every consumer). A lint
  finding alone is a warning; a count that grows fails, and a new package starts at 0. A count only shrinks: growth
  needs `--update-baseline --allow-growth "<why>"` in the same change.
- `node libraries/u2/scripts/audit-consumers.mjs --check`: the counts alone, every tree.
- u2 reports a misuse as `u2 usage: <code>: <what to do>`; `DG.U2.Usage.mode = 'throw'` makes it throw.

## When u2 lacks something

Say so in your answer, naming the component and the job, and build the closest shape the entry points allow in your
own package. Never import from a path inside the package, copy or patch u2's sources, subclass a control to reach
its internals, style `.u2-*` classes or `[data-u2…]` attributes, or look u2's parts up by class name.
