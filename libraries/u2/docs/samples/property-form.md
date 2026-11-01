---
title: A form from DG.Property metadata
group: Inside Datagrok
order: 40
summary: Describe the settings once as properties; the form, its editors, validation and write-back come from them.
live: no
---

When an object carries metadata — `DG.Property` lists, a viewer's properties, a domain row's schema — don't
hand-build an input per field: `propertyForm(props, target)` generates the editors by type (choices become a
choice, `nullable: false` refuses an empty value, `min`/`max` bound a number, `category` groups), writes each valid
edit into the target, and is a `Form` with the same `problem` and `validate()`. `overrides` replaces one field's
editor; `editors: 'auto'` asks the platform for its own (molecule, column, color). Persisting stays with you: these
are the user's own settings, so each valid edit is saved as it is made (`onChanged`), with no Save button. Settings
shared with others or set by an admin stage their edits behind "SAVE N CHANGES" instead.

## What it takes

{lines} lines for a validated settings editor that saves every valid edit. u2 handles choosing the editors, the
validation, writing valid values through (never an invalid one), and `refresh()` when the target changed elsewhere.
