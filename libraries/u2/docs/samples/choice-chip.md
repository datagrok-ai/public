---
title: Small choices as chips
group: Inputs & validation
order: 35
summary: Three R groups as chips reading "R1 H ▾": a press opens the caps, each with what it stands for; X turns its chip red.
live: yes
height: 160
---

A `ChoiceChip` is a choice that lives inline among others, its label inside it: "R1 H ▾". A press — or ArrowDown
— opens a menu of the choices, each with its `description` in grey and the current one checked, under the chip
when it fits and above it otherwise; a pick is `onInput`, two clicks from start to end. `invalid` (a reason, or a
signal of one) turns the chip the error color at once and says why as its tooltip; `disabled` takes the reason
the user reads. The value is `chip.value`, `chip.signal` wires it and `chip.follow(read)` keeps it on a model.
`actions` add entries under the choices ("Edit definition…"). A field of a form, in its label column, is a
`ChoiceInput` instead.

## What it takes

{lines} lines and no CSS. u2 handles the menu, its keyboard (arrows, Enter, Escape back to the chip), the
placement, `aria-haspopup`/`aria-expanded`/`aria-invalid`, the description and the reason as the tooltip.
