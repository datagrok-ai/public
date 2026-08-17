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
- Required comes from the model (`nullable: false`), not from a hand-added asterisk;
  validation messages render under the input, and `form.validate()` focuses the first
  offender.
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

## Anti-patterns

- Hand-building an input per property when the object carries metadata — that is
  `propertyForm`'s job, including the platform's own editors (`editors: 'auto'`).
- A `<table>`/grid of label+input pairs instead of `Form` — loses validity aggregation,
  label alignment and keyboard flow.
- Centered or right-aligned value columns in detail panes.
- The form saving to the server itself — mutation stays in the form, persistence with the
  caller.
