---
title: Pickers with suggestions
group: Inputs & validation
order: 30
summary: Tags from an async search, a suggesting text box, chips for a small set — the same input contract.
live: yes
height: 300
---

The pickers are inputs like any other: a `value`, a label, validation. `TagsInput` shows the picked items as
chips and suggests more as you type, from a static list or an async `source(query, abort)`; `SuggestInput` is a
text box that offers completions; `ChipsInput` toggles a small set of choices in place. An async source is
debounced, and a stale answer is dropped when a newer query is on its way. A choice of a table's column is not one
of these: it is u2/dg's `columnInput`, which lists what the platform's own pickers list.

## What it takes

{lines} lines. u2 handles the debounce, the cancellation (the `AbortSignal`), the popup, its placement and keyboard
(arrows, Enter, Escape, Backspace on an empty box removes the last chip) and the "no matches" row.
