---
title: Rich rows with actions
group: Lists & tables
order: 10
summary: Ten thousand two-line rows; icon actions on hover, the full list on right-click — declared once.
live: yes
---

A row is a small layout, not a line of text. Declare the row's actions once as `Action`s — a `text`, an optional
`icon`, a `run` that returns a promise or nothing: `rowActions` shows the ones with an icon on the hovered row
(they keep their room, so nothing shifts, and stay reachable by Tab — a test reaches them with the keyboard, as
no script can hover), and `contextActions` opens all of them on right-click — the row is selected first. The list
renders only what is in view. With `u2-adaptive` on the list, parts marked `u2-p2` give way as the pane narrows.

## What it takes

{lines} lines of TypeScript and {cssLines} of CSS (the row's own layout). u2 handles the virtualization, the
selection and its keyboard, the hover block, the context menu and its keyboard, the timestamp's tooltip, and the
responsive hiding.
