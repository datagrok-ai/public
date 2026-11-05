---
title: Signals, computed and effects
group: Reactive state
order: 10
summary: State in signals, inputs bound to them two-way, a computed total, an effect, and batch for one update.
live: yes
---

A `signal` holds state; a `computed` derives from signals and recomputes only when one it read changed; an
`effect` runs a side effect now and after every such change, and stops with whoever built it — here the
component. An input given `bind: signal` edits that signal and follows it, so the model is the signals and the
form is a view of them. `batch` makes several writes one change: Reset writes three signals and the effect runs
once. Read with `peek()` where a read must not subscribe.

## What it takes

{lines} lines. u2 handles the dependency tracking, the two-way binding (without echo loops), glitch-free
recomputation, and stopping the effect and the bindings with their owner.
