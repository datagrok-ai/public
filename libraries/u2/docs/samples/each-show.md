---
title: Keyed rows with each() and show()
group: Reactive state
order: 20
summary: Rows that follow a list signal — added, removed, edited in place with focus and caret kept.
live: yes
---

`each(host, items, keyOf, build)` keeps one element per key: built once in its own scope, kept while its key
stays — so typing in a row's input keeps focus and caret while the list changes around it — and disposed with
everything it owns when its key leaves. `update` patches a kept row whose item changed (here nothing to do: the
input already shows the value); `empty` is shown while there are no rows. `show(host, when, build)` is the same
idea for a condition: the warning is built when the total passes the limit and disposed when it drops back.
Both need an owner building — here the component. Build a row's inputs, listeners and effects inside the row,
never on the control around it.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (a width and two colours; `divH(row, {gap: 'm'})` spaces the
rows). u2 handles the keyed reconciliation, the per-row scopes and their disposal, the empty state and the focus
kept across edits.
