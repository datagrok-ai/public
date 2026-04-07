---
title: "Docking"
sidebar_position: 4
format: mdx
unlisted: true
---

Datagrok provides a flexible window management system. You can either drag windows out and position them manually or set them up automatically.

![](img/table-view-add-viewers-dock.gif)

## Top-level window docking

Use `grok.shell.dockManager` to dock, undock, or reposition windows.
See also `View.dockNode`

## Nested docking

Some views contain a nested docking manager that lets you manage
windows within that particular view. These windows cannot be undocked and docked at
the top level or within a different view. See `dockManager` property of the
View class.
