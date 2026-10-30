---
title: Menus, a menu bar and a context menu
group: Navigation
order: 10
summary: One fluent Menu for popups, submenus, checked items and right-click; a MenuBar for a panel's commands.
live: yes
height: 340
---

`Menu` is built fluently — `item`, `separator`, `group` … `endGroup` for a submenu — and shown under an element
or at a point. A `MenuBar` builds each menu on every open, so a checked item always shows the current state;
`bindContext(el)` gives an element a right-click menu that closes when whoever built it goes. An item takes a
`shortcut` (shown), `check`, `disabled` (the reason the user reads), an `icon` and a grey `detail`. A view's own
commands belong in the ribbon (see Inside Datagrok); a menu bar is for a panel inside it.

## What it takes

{lines} lines. u2 handles the menu keyboard (arrows, Enter, Escape, → into a submenu), typeahead on the items,
placement that flips at the window's edge, focus returning to the opener, and `aria-checked`.
