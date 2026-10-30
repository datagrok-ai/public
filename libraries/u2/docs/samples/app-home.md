---
title: An app home with a catalog
group: Inside Datagrok
order: 50
summary: The header every app opens with, then the records the user comes back to — loaded, refreshable, remembered.
live: no
---

An app's home: `AppHeader` (what the app does, Learn more), then what the user starts from — here the campaigns as
a `CatalogTable` over an `AsyncValue`. The first load shows a skeleton, a refresh keeps the rows (the refresh icon
spins), a failed first load shows the reason with Retry, and the user's sort, groups and columns are kept under
one key. `ViewHost.create` builds it all under the view's host, so the view's close disposes all of it — the header
follows its own width, which is why it needs that owner.

## What it takes

{lines} lines for the whole home view. u2 handles the header, the loading states and Retry, the catalog with its
menus and memory, and the hosting.
