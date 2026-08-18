# UI Tests changelog

## v.next

* Tests: `Viewers` (close, Reset/Set default properties), `View: Events` layout events and `lineChart.api` now add viewers to a freshly added — and therefore active — table view. A view that is not active has a zero-sized dock area, so `DockView.addViewer` queues the viewer instead of docking it: it never reaches `view.viewers`, is missing from `saveLayout()`, and default settings are not applied.
* Tests: `Viewers: Core Viewers` excludes `Scaffold Tree` — it is a Chem JS viewer, and the JsViewers filter runs before Chem is loaded, so it leaked into the core-viewer list and timed out.
* Tests: `pca` drives `ML | Analyze | PCA...` through the view ribbon menu (as `cluster` does); `grok.shell.topMenu` is the app-level menu and its ML group has no `Analyze` child.
* Tests: `Layouts.delete` no longer requires the `Upload` context-menu item — the toolbox `Save` button already persists the layout, so `ViewLayoutMeta` stops offering `Upload`.
* Tests: `UI: Sharing` scopes gallery lookups to the opened view instead of `document`, so the search filters and the assertion read the same grid.
* Tests: name the projects created in `tag.add` and `projects.api` — projects now require a non-empty name on save (datlas projects_service guard), which was throwing `Project name cannot be empty`.
* Tests: `pca` now skips gracefully when the ML menu (EDA package) isn't deployed, instead of crashing with `Cannot read properties of undefined (reading 'find')`.
* Tests: `UI: Sharing` uses a unique per-run entity name so a leftover from a prior aborted/failed run can't make the gallery search return >1 result ("more than one testing entity present").
* Tests: `View: Events` onViewLayoutApplying/onViewLayoutApplied now search the whole dock tree for the viewer — they only checked the top-level `children` (splitters), never the nested viewer leaves, so `viewerElement` was always null (chronic failure on dev too).
* Tests: `Viewers: Core Viewers` runs Tile Viewer on a small demo df — the 100-row chemistry df made it build a DOM tile per row and freeze the browser past grok test's 180s inactivity watchdog, aborting the entire puppeteer pass (0 results).

## 1.0.12 (2024-09-02)

* Test fixes update 

## 1.0.4 (WIP)
