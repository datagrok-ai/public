# Recipe: forms and detail panes

## Which component

| Situation | Use |
|---|---|
| Editing an object that has Property metadata (`DG.Property`, viewer props, domain rows) | `propertyForm(props, target)` / `objectForm(source)` (u2/dg) — editors, validation and writing valid values through are generated |
| A hand-assembled set of inputs | `Form` — shared label column, aggregated validity, Enter navigation |
| Descriptor-driven settings panel (many small typed values) | `PropertyGrid` |
| Read-only details of a selected item, in a pane of their own | `tableFromMap` — muted keys, wrapping values, a value that is an element shown as it is |
| A value shown as an element among a form's inputs (a structure drawing, a chart) | `DynamicInput` — a row of the form, in its label column, that grows to the element |
| Picking a column of a table | `columnInput(label, df, …)` (u2/dg) — what the platform's pickers list, hidden `~` columns left out; never a `ChoiceInput` over `df.columns` |
| A short list of typed rows edited in place (a template's fields, its stages) | `RecordsInput` — a column per field, add / remove / `orderable` |
| A subset where the order matters too (the functions to run, in order) | `ChecklistInput({orderable: true})` |
| Asking before what cannot be undone | `confirmDialog({title, message, confirmText: 'Delete', destructive: true})` |

The two read-only shapes do not compete: a detail pane is a `tableFromMap`; inside a form, where the rows are
inputs, a value drawn as an element is a `DynamicInput`, never an element put on an input's rail.

## Layout rules

- Labels form one left column (the form aligns them); values are **left-aligned**, starting
  right after the label column — never centered, never pushed to mid-pane.
- Long read-only values **wrap** (`tableFromMap` does this by default); inputs ellipsize.
- `required: true` marks a field to fill in: the caption's asterisk and `aria-required`, and an empty value fails
  it (blank text is empty). `nullable: false` (the model's) refuses an empty value the same way, with the same
  message, but without the asterisk — a choice that always holds a value is not a field to fill in; `required:
  false`, given, lets an empty value through over it. The message is the input's own: "Value can't be empty", for
  a file or a key "Can't be empty.", for files "At least one file is required.". The required check comes first:
  a validator never sees the empty value it refuses, so a format check need not pass on nothing.
  Never a hand-added asterisk. An input shows its error — a required one left empty, a validator's — once
  edited, left or validated, and `form.validate()` shows every error and focuses the first offender.
  Save and OK go by the verdict from the start: disabled, the reason in their tooltip. Hidden and
  disabled inputs never fail the form.
- A validator may read other inputs and signals: it re-judges its input when any of them changes, from
  the moment it is added (`to.addValidator((d) => d < from.value ? 'Ends before it starts' : null)`).
  `revalidate()` is for state that is not a signal.
- Every input takes the same contract options: `readOnly` (selectable, not editable), `required`,
  `visible`, `disabled` (the reason the user reads while it is off) — each a value or a signal, followed live —
  plus `focus({select})`, `fill` and `onEnter` / `onEscape` handlers that keep the key from the form and the
  dialog — a handler answering `false` leaves it to them.
- The rail beside an editor takes buttons and links (`addOptions(iconButton(…))`) and works on a read-only
  input: a download, or a pencil that turns editing on; `addOptions(button, {edits: true})` marks one that
  changes the value, which read-only stops. The rail is one row high.
- Entity values render as chips (`chip(entity)`), not as text.
- The form's button row is `form.addButtons((row) => row.append(…))`: under the rows, at the form's end — right-
  aligned, where a platform dialog's buttons are. `{align: 'editors'}` starts it where the editors start (as
  `ui.buttonsInput` does), `{align: 'start'}` at the left; `{host}` puts the row elsewhere (`page.footer`, the
  action row a `Page` pins under a body that scrolls) and keeps its look. A button row outside any form is
  `divH([cancel, save], {gap: 'm'})` in a host of your own.
- A form is 700px at most: `new Form({maxWidth: 820})` sets its own, and in a `Page` given a `maxWidth` its forms
  go as wide as the page. The forms on one `Page` share a label column (`labelColumn: false` keeps one's own);
  `density: 'compact'` takes the 24px rows of the compact token set. The space between a form's sections is
  `--u2-form-section-gap`, set on an element of your own (README › Styling).

## The shape (metadata-driven, the preferred path)

```ts
import {propertyForm} from '@datagrok-libraries/u2/dg';

const form = propertyForm(props, entity, {
  overrides: {smiles: {input: moleculeInput('SMILES')}},   // per-field replacement when needed
  onChanged: (name) => markDirty(name),
});
form.addButtons((row) => row.append(button('Save', async () => {
  if (!form.validate()) return;
  await grok.dapi.groups.save(entity);           // persisting stays with the caller
  form.refresh();
})));
// inspection: form.properties, form.input('name').validity, .addValidator(...), .disabled
```

## Reading and writing a value

`input.value` is the current value (followed inside an effect or a `computed`), `input.value = x` writes
it, `input.peek()` reads it without following, and `input.signal` is the signal to wire (`bind`,
`span(input.signal)`). A handle names the input for what it holds (`editor.person.value.trim()`).

```ts
const personEditor = () => component(() => {
  // built inside the builder: owned by the editor and disposed with it
  const person = new TextInput({label: 'Name', required: true});
  const greeting = computed(() => `Hello, ${person.value.trim()}`);
  return {root: divV([person, span(greeting)]), person};
}, {name: 'person-editor'});
```

## Rows, picks and confirmations

```ts
const stages = new RecordsInput<{name: string}>({
  label: 'Stages', orderable: true, addText: 'Add stage', required: true,
  columns: [{name: 'name', header: 'Name', input: (v) => new TextInput({inline: true, value: v ?? '', required: true})}],
  value: template.stages.map((name) => ({name})),
  rowName: (_, i) => `Stage ${i + 1}`,                  // "Stage 3: Value can't be empty"
});
stages.addValidator((rows) => new Set(rows.map((r) => r.name)).size < rows.length ? 'Stages must differ' : null);

const order = new ChecklistInput({label: 'Functions', items: functionNames, orderable: true,
  value: saved.order.filter((f) => !saved.hidden.includes(f))});
// persist: {order: order.value, hidden: functionNames.filter((f) => !order.value.includes(f))}

async function removeCampaign(name: string, code: string): Promise<void> {
  if (await confirmDialog({title: 'Delete campaign',
    message: `Are you sure you want to delete campaign ${name} (${code})?`,
    detail: 'This cannot be undone.', confirmText: 'Delete', destructive: true,
    onConfirm: () => deleteCampaign(code)}))    // busy until done; a failure stays in the dialog
    notify.info(`Campaign ${code} deleted`);
}
```

- Each cell is an inline input from the column's `input`; the value is the records, one new array
  per edit, the untouched records kept as they were. A cell's verdict is its red underline; the
  input's own line names the row.
- A dialog with a fixed height (`show({height})`) or `show({resizable: true})` has its body fill it.
  A dialog that belongs to a platform view is shown with u2/dg `ViewHost.of(view).show(dialog)`.

## Load, then build rows

Rows that depend on data the app loads: hold the load in an `AsyncValue` and build from its state. While it
loads, the area says so; a failure shows why, with Retry; once ready, the rows are built from the value — and
built again, owned by the new run, on every refresh.

```ts
const fieldsEditor = () => component(() => {
  const template = new AsyncValue(() => loadTemplate('Kinase hit-to-lead'));
  const area = div();
  each(area, () => [template.state.value], (s) => s.kind, (s) =>    // one row per state: rebuilt when it changes
    s.kind === 'ready' ? new Form().addAll(s.value.fields.map((f) => new TextInput({label: f, required: true}))).root :
    s.kind === 'error' ? emptyState(s.message, {error: true, actions: [button('Retry', () => template.refresh())]}) :
    skeleton());
  return area;
}, {name: 'fields-editor'});
```

- `AsyncView.owned(fetch, render)` is the same for a list (a search): it renders loading, the results, an empty
  answer and a failure with Retry itself — docs/samples/async-search.ts.
- A pane that loads once and then only edits needs no rebuild: build the form when the value is ready
  (`template.state.value.kind === 'ready'`) and load it with `input.reset(value)`.

## Dates and the text of a value

A `Date` reads in the browser's zone, the way the platform's grid shows a datetime column. `Format.value(v)` is
the text u2 shows for any value — a date in the user's platform date format (`M/d/yyyy` until they choose one),
with the time of day where it carries one; `Format.date(d, {format: 'yyyy-MM-dd'})` takes a platform pattern.
`timestamp(d)` draws the short date for a list row with the full date and time as its tooltip. A value that
stands for a day, not a moment, stored as midnight UTC (a domain `date`), reads in UTC with `{utcDates: true}`:
read in the local zone it is the day before west of Greenwich.

```ts
const measured = new Date(2026, 9, 27, 14, 30);           // a moment, in the browser's zone
const dueDay = new Date(Date.UTC(2026, 9, 27));           // a day, stored as midnight UTC
const dates = divV([
  span(Format.value(measured)),                           // 10/27/2026 14:30 — the user's format
  span(Format.date(dueDay, {format: 'MMM d, yyyy', utcDates: true})),   // Oct 27, 2026 everywhere
  timestamp(dueDay, undefined, {utcDates: true}),
]);
```

`CatalogTable` cells and group headers without a `render` go through `Format.value`; a `DataTable` or
`BasicTable` cell without one shows `String(value)`, so give a date column a `render` with `Format.date`.

## Multi-step flows

`Wizard` for steps known up front. A step's `commit` and the wizard's `onFinish` get
`{signal: AbortSignal}`: Esc or ✕ aborts the run in flight; a throw or a rejection keeps the step
and shows the message in the footer. React to `wizard.on('step' | 'finish' | 'cancel', …)`; read
`wizard.currentStep`, `busy` and `error`. In a spec the same events are `on` entries (see
[specs](specs.md#events)).

## A docked editor panel

An editor beside a grid (Bio's monomer editor, `packages/Bio/src/utils/monomer-lib/monomer-manager/`)
lives in a panel the user resizes and on screens from 720 to 1440 px tall. The shape that wastes no
space at any of those sizes:

```ts
const column = new LabelColumn();                          // one caption column for every form below
const fields = new Form({layout: 'auto', labelColumn: column})
  .addColumns([symbol, name, polymer, type], {minWidth: 176}); // two to a row, one to a row when narrow
const save = button('Save', () => onSave(), {primary: true, disabled: problem}); // the reason is the tooltip
const details = new Section({title: 'Details', hint: summary, validity: detailsProblem, chevron: 'end',
  expanded: false, persistKey: 'my-editor.details'});      // folded: what it holds, or what is wrong in it
fields.addElement(details.root).add(id, details.body)      // last in the fields: it opens right above Save
  .addGroup('Colors', [line, background, text], {captionsAfter: true, host: details.body});
const footer = div([save, statusLine(problem, computed(() => `Saves into ${library.value}`)), details.header]);
const scroller = div([fields.root], 'my-scroller');        // display: flex; flex-direction: column-reverse
panel.append(div([header, preview, scroller, footer], 'my-grid'));
```

- Lay the panel out as a grid of four rows: the header; the big control (a sketcher, a preview)
  taking whatever height the fields leave, down to a floor (a share of the panel's height, 300 px at
  least); the fields at their own height, scrolling once the big control is at its floor; and
  Save, pinned in the last row. Make the panel a size container and key the floor off it.
- With Save pinned at the bottom, whatever opens makes room by moving the rows above it up. Open a
  sub-record's fields under its chips, where users look for them, and a fold whose toggle sits in the
  pinned row (Details) above that row. Give the scroller `flex-direction: column-reverse` around one
  child, so it stays scrolled to the bottom and what just opened stays in view; the trigger moves up
  by what opened under it. (Opening everything upward keeps the trigger still, but read as
  unexpected.) Fade the scroller's top edge while rows are hidden above it.
- Every editor keeps its underline and its caption right-aligned beside it, as in every platform
  form (`quiet` hides underlines at rest, and a form people fill in then reads as half fields, half
  text). `addColumns` keeps every column's captions and editors in line, each editor filling its
  column, and `minWidth` is the narrowest a caption and its editor get — the longest choice item and a
  dozen characters of text — before the row goes one to a row.
- Save with a `statusLine` beside it says what Save will do, or why it cannot yet, without hovering.
  The fold's header (`chevron: 'end'`) sits at the end of that row with a `hint` summarizing what it
  holds (an element can carry swatches); open it where there is room (from a `ResizeObserver`, until
  `section.userSet`). A form whose panes stand on their own rows lines them up with
  `new Accordion({flush: true})`.
- A small set of sub-records (R groups) reads as chips that wrap: a `ChoiceChip` each ("R1 OH ▾"), its
  menu of choices built on every open — an item's `description` says what each stands for, the current
  one is checked — so a change is two clicks and nothing in the panel moves; `invalid` turns a chip red
  with its reason. Fields needed rarely (a cap's raw definition) open under the chips from the chip's
  `actions` ("Edit definition…"), and `active` keeps the chip pressed while they show. Build the chips
  with `each()` (docs/samples/choice-chip.ts).
- Tabs in a narrow panel hide most of the record and leave an empty band under the short tabs.
- Keep a canvas that takes the mouse wheel (a sketcher, a chart) out of any scroller.

## Anti-patterns

- Hand-building an input per property when the object carries metadata — that is
  `propertyForm`'s job, including the platform's own editors (`editors: 'auto'`).
- A `<table>`/grid of label+input pairs instead of `Form` — loses validity aggregation,
  label alignment and keyboard flow.
- Centered or right-aligned value columns in detail panes.
- The form saving to the server itself — mutation stays in the form, persistence with the
  caller.
- A column list built from `df.columns` into a `ChoiceInput`: it offers the hidden `~` columns the platform keeps
  out of every picker — `columnInput` lists what the platform lists.
