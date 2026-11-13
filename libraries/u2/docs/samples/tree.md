---
title: A tree that loads as it opens
group: Lists & tables
order: 40
summary: Folders whose children arrive from the server when expanded — a loading row, an error row, a context menu.
live: yes
---

`VirtualTree` takes nodes as data; a node's `children` may be a function that answers a promise — called when
the node first opens, with a loading row meanwhile and an error row (with the message) if it fails. The tree is
virtualized like a list, so a thousand open nodes cost what is on screen. `selectedNode` is a signal;
`contextActions` gives every node a right-click menu, whose actions' `run` returns a promise or nothing. Open
"Restricted" to see a failure.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the height). u2 handles the lazy loading and its states, the
tree semantics (`role=tree`, ←/→ to fold and open, ↑/↓ to move), virtualization and the context menu.
