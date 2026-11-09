---
title: Buttons that say why
group: Basics
order: 30
summary: A disabled button names the reason; a split drop-down, an icon button and a segmented choice.
live: yes
---

`button(text, onClick, {disabled})` takes the *reason* it is disabled — a string or a signal of one, `null` while
it can run — and shows it as the tooltip, on hover and on keyboard focus: a disabled u2 button keeps its tab stop,
and a press on it does nothing. Every off in u2 takes a reason the same way (inputs, actions, menu items, a
dialog's `okDisabled`); there is no `enabled`. `dropDownButton` with `split` runs its main action and opens the
menu from the arrow; the menu is built on every open, so it is never stale. `iconButton` needs its tooltip, which
is also its accessible name. In a `Form`, `addButtons` puts the row at the form's end, on the right.

## What it takes

{lines} lines. u2 handles the disabled state (`aria-disabled`, the reason as tooltip and `aria-description`),
the menu's keyboard (arrows, Enter, Escape), its placement, and the segmented choice's radio semantics.
