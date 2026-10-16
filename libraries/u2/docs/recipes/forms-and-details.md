# Recipe: forms and detail panes

## Which component

| Situation | Use |
|---|---|
| Editing an object that has Property metadata (`DG.Property`, viewer props, domain rows) | `propertyForm(props, target)` / `objectForm(source)` (u2/dg) — editors, validation and write-through are generated |
| A hand-assembled set of inputs | `Form` — shared label column, aggregated validity, Enter navigation |
| Descriptor-driven settings panel (many small typed values) | `PropertyGrid` |
| Read-only details of a selected item | `tableFromMap` — muted keys, wrapping values |

## Layout rules

- Labels form one left column (the form aligns them); values are **left-aligned**, starting
  right after the label column — never centered, never pushed to mid-pane.
- Long read-only values **wrap** (`tableFromMap` does this by default); inputs ellipsize.
- Required comes from the model (`nullable: false`) or the input's `required` option, not from a
  hand-added asterisk; an empty required input shows its error once edited, left or validated,
  and `form.validate()` shows every error and focuses the first offender. Hidden and disabled
  inputs never fail the form.
- Every input takes the same contract options: `readOnly` (selectable, not editable), `required`,
  `visible`, `enabled` — each a value or a signal, followed live — plus `focus()`, `fill` and
  `onEnter` / `onEscape` handlers that keep the key from the form and the dialog.
- Entity values render as chips (`chip(entity)`), not as text.

## The shape (metadata-driven, the preferred path)

```ts
import {propertyForm} from '@datagrok-libraries/u2/dg';

const form = propertyForm(props, entity, {
  overrides: {smiles: {input: moleculeInput}},   // per-field replacement when needed
  onChanged: (name) => markDirty(name),
});
form.addButtons((row) => row.append(button('Save', async () => {
  if (!form.validate()) return;
  await grok.dapi.groups.save(entity);           // persisting stays with the caller
  form.refresh();
})));
// inspection: form.properties, form.input('name').validity, .addValidator(...), .enabled
```

## Multi-step flows

`Wizard` for steps known up front. A step's `commit` and the wizard's `onFinish` get
`{signal: AbortSignal}`: Esc or ✕ aborts the run in flight; a throw or a rejection keeps the step
and shows the message in the footer. React to `wizard.on('step' | 'finish' | 'cancel', …)`; read
`wizard.currentStep`, `busy` and `error`. In a spec the same events are `on` entries (see
[specs](specs.md#events)).

## A docked editor panel

An editor beside a grid (Bio's monomer editor, `packages/Bio/src/utils/monomer-lib/monomer-manager/`)
lives in a panel the user resizes and on screens from 720 to 1440 px tall. The shape that wastes no
space at any of those sizes and moves nothing under the pointer:

```ts
const column = new LabelColumn();                          // one caption column for every form below
const fields = new Form({layout: 'auto', labelColumn: column, quiet: true, captionAlign: 'left'})
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
- With Save pinned at the bottom, whatever opens grows the fields upward. Put what a click opens
  above what was clicked — a sub-record's row above its chips, a fold's body above its toggle — and
  give the scroller `flex-direction: column-reverse` around one child, so it stays scrolled to the
  bottom: then nothing moves under the pointer, whether the big control gives way or the fields scroll.
  Fade the scroller's top edge while rows are hidden above it.
- `quiet` reads values as text: no column of underlines, a line where the pointer, the focus, an
  empty field or an error is. Left captions give the panel one left edge; `addColumns` keeps every
  column's captions and editors in line, and `minWidth` is the narrowest a caption and its editor
  get — the longest choice item and a dozen characters of text — before the row goes one to a row.
- Save with a `statusLine` beside it says what Save will do, or why it cannot yet, without hovering.
  The fold's header (`chevron: 'end'`) sits at the end of that row with a `hint` summarizing what it
  holds (an element can carry swatches); open it where there is room (from a `ResizeObserver`, until
  `section.userSet`). A form whose panes stand on their own rows lines them up with
  `new Accordion({flush: true})`.
- A small set of sub-records (R groups) reads as chips that wrap, one opened at a time into its row:
  build them with `each()`, and give their plain elements `effect()` and `on()` — both owned by the row.
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
