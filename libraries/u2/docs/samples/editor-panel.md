---
title: A dense editor panel
group: Forms
order: 20
summary: Short fields two to a row, a grouped caption, the rest folded under a section that says what it holds.
live: yes
---

The shape for an editor docked beside a grid. `addColumns` puts short fields two to a row — one to a row when
the panel is narrow — with every caption and editor still in line. `addGroup` puts several inputs under one
caption (Colors: Line, Background, Text). What is rarely needed folds into a `Section` whose `hint` says what it
holds, whose `validity` opens it when something inside is wrong, and whose `persistKey` remembers the user's
open or close. A `statusLine` beside Save says what it will do, or why it cannot yet; the button row sits at the
form's end, on the right.

## What it takes

{lines} lines, no CSS. u2 handles the responsive columns, the caption column shared across the rows, the
section's open or closed state (remembered in the browser) and the status line.
