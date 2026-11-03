---
title: An app view in the shell
group: Inside Datagrok
order: 10
summary: A package app whose search and refresh live in the ribbon, its count in the status bar, its URL in sync.
live: no
---

A u2 app runs inside the Datagrok shell and uses its chrome instead of rebuilding it. `ViewHost.create({name,
content})` runs the builder under the new view's host and makes a platform view of what it answers: the content as
the view's root, the `ribbon` groups in the view's ribbon, a `status` signal in the shell's status bar, and `path`
in the address bar. Closing the view disposes the host and everything built under it — the async value, the inputs
in the ribbon, the table. A view's main controls go to the ribbon; a toolbar inside the content is for a panel
within it.

## What it takes

{lines} lines for a searchable, sortable, groupable catalog of the server's projects with a refresh and a live
count. u2 handles the hosting and its disposal, the ribbon and status placement, the debounced server search
(`deps`), the catalog and the URL.
