---
title: Numbers that load
group: Async & data
order: 20
summary: KPI cards over async values — a skeleton while loading, the message on failure, refreshed together.
live: yes
---

`AsyncValue` is the single-value twin of the search source: one abortable run at a time, a stale answer dropped,
`refresh()` to run again. `StatCard` takes one as its `source` and shows a skeleton while it loads, the message
if it fails, the formatted value when ready. An `AsyncValue` built inside an owner is disposed with it, so a run
in flight is aborted when the view closes.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the grid of cards). u2 handles the states, the abort, the
number formatting, the delta's arrow and colour, and the disposal.
