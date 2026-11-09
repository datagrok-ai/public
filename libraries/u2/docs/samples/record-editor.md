---
title: Edit and save a record
group: Forms
order: 10
summary: Load a record, edit it, Save only when valid and changed — the reason always on the button.
live: yes
old: ui.* today
---

The everyday editor: inputs in a `Form`, loaded with `reset(value)` (no change event, no stale "touched" state),
Save disabled with the reason — the first failing field, or "No changes to save" — Revert back to the saved
record, and a status line that says what Save would do. Save's reason is `form.problem` (the first failing field,
touched or not); the status line shows `form.shownProblem` (the message the user already sees on a field they
edited or left), so a fresh form is never red. Every rule is a `computed` over the inputs' values, so there is no
refresh function to call and nothing to unsubscribe.

## What it takes

{lines} lines against {oldLines} with `ui.*` and `DG.InputBase` (second tab), which subscribes to every input,
recomputes the button state by hand, and never unsubscribes. u2 handles the touched state, the messages, the
reason on the disabled button (also on keyboard focus) and the cleanup.
