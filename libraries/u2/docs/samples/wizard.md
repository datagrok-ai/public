---
title: A wizard in a dialog
group: Dialogs
order: 30
summary: Steps known up front, each gated by its form; the last one reviews and creates — busy until it is done.
live: yes
height: 440
---

`Wizard` takes its steps up front: a title, the content (a function is built on first visit), a `canProceed`
gate — a form's `problem` composes directly, and the reason shows beside NEXT — and an optional `commit` that runs
before leaving the step. `onFinish` may return a promise: the footer waits, Escape or ✕ abort its `signal`, and a
failure keeps the step with the message. `openInDialog` puts the whole flow in a dialog whose footer is the
wizard's own; reopening returns to the step the user left. The review shows each input's live value with
`span(input.signal)` — `span(input.value)` would show the value it had then.

## What it takes

{lines} lines. u2 handles the step rail (visited steps are clickable), BACK / NEXT / FINISH and CANCEL, the gate
and its reason, Enter for NEXT, panels kept (not rebuilt) between visits, and the abort on Escape.
