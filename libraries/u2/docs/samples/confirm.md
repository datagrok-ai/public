---
title: Ask before what cannot be undone
group: Dialogs
order: 20
summary: One call — the question, what goes with it, CANCEL and the action as a verb; busy until the work is done.
live: yes
height: 300
old: ui.dialog today
---

`confirmDialog()` asks before an action that cannot be taken back and answers whether it ran. The footer is CANCEL
and the action as a verb (DELETE), the message names what it acts on, `destructive` adds the trash icon. With
`onConfirm` the work runs inside the dialog: it stays open and busy until the work is done, and a failure stays
in its error row instead of a balloon nobody connects to the click. Delete KHL-13 to see one.

## What it takes

{lines} lines against {oldLines} with `ui.dialog` (second tab). u2 handles the busy footer, the spinner, the
error row, Enter and Escape, the focus and its return, and the answer as a promise.
