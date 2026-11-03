---
title: An app that lives in a TableView
group: Inside Datagrok
order: 30
summary: A session owned by the view, its state as signals, its actions pinned at the ribbon's end, its dialogs pinned to it.
live: no
---

When the grid is the platform's, the app's session belongs to the view's `ViewHost`: `ViewHost.of(view)` is one
host per view, whichever wrapper you hold, and it is disposed when the shell removes the view — and with it every
subscription, timer and dialog built under it (`host.build(…)`). `host.signals` gives the view's current row,
selection, columns and grid as signals: an effect over `grid` re-subscribes when a layout swaps the grid, and what
it subscribed for the old grid goes with that run. `host.ribbon` puts the app's actions at the end of the ribbon
(they fold into one "…" menu when the ribbon runs out of room; an action's `hint` pulses it) and `host.show` hides
a dialog while another view is current. The session is a plain object, not a control: what it builds is the host's.

## What it takes

{lines} lines for the whole integration. u2 handles the lifetime, the signals over the platform's events, the
ribbon group and its folding, the hint, the status text, and the dialog that follows its view.
