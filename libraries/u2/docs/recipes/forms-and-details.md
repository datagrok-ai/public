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

## Anti-patterns

- Hand-building an input per property when the object carries metadata — that is
  `propertyForm`'s job, including the platform's own editors (`editors: 'auto'`).
- A `<table>`/grid of label+input pairs instead of `Form` — loses validity aggregation,
  label alignment and keyboard flow.
- Centered or right-aligned value columns in detail panes.
- The form saving to the server itself — mutation stays in the form, persistence with the
  caller.
