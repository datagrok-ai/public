---
title: "Docking"
sidebar_position: 4
mdx:
  format: mdx
unlisted: true
description: Dock, undock, and reposition windows manually or programmatically using Datagrok's flexible window management system.
keywords:
  - dock windows
  - undock
  - arrange panels
  - split screen
  - dockManager
  - window management
  - reposition panels
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
