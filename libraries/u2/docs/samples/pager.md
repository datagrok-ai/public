---
title: A list that loads as you scroll
group: Async & data
order: 30
summary: Pages accumulate into one list as the user nears the end; the total, the loading state, Retry.
live: yes
---

`AsyncPager` keeps the pages in one `items` signal: `reset()` reads the total and loads the first page,
`loadMore()` the next one — a no-op while one is in flight or once a short page said it was the last — and a
`reset()` drops answers it outran. Hand `items` to a `VirtualList` and call `loadMore()` when the user nears the
bottom. A failed page keeps what was loaded and is retried by the next `loadMore()`. The list takes the height
its host gives it: here a class of the app's own on the list's root.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the height). u2 handles the page counting, the end of data,
the race with `reset()`, the virtualization and the states, and disposes the pager with the control that built
it (a page in flight is dropped). The near-the-bottom listener is the app's own.
