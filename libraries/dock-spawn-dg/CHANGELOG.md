# dock-spawn-dg changelog

## 0.2.0 (2026-11-27)

* `<dock-spawn-ts>`: Esc no longer closes the active panel (`escClosesWindow` is false); dialogs still close on Esc
* Disabled the "New browser window" context-menu item by default (DockConfig.enableBrowserWindows is false)
* Fixed a dock mounted or re-rendered while hidden losing its panel ratios (panels split evenly once shown); panels now dock, and removed panels close, once it has a size
