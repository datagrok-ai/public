---
title: Search a server as you type
group: Async & data
order: 10
summary: Loading, results, nothing found and failed-with-Retry — one AsyncView over a fetch that can be aborted.
live: yes
height: 320
old: ui.* today
---

`AsyncView.owned(fetch, render)` renders the states every async area has — loading (a skeleton), the results, an
empty answer, and a failure with Retry — over a `fetch(query, abort)`. Each new query aborts the one in flight and
a late answer is dropped, so results never arrive out of order; while a newer answer loads, the old results
stay, marked busy. Search for "error" to see a failure, "zzz" for nothing found.

## What it takes

{lines} lines against {oldLines} with `ui.*` (second tab), which debounces with a timer, discards late answers
with a counter and cannot cancel the request. u2 handles the debounce, the abort, the stale answers, all four
states with Retry, and `aria-busy`.
