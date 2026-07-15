# Usage Analysis — Tabs Architecture

How the tabbed Usage Analysis app is wired: app entry → tabs container → shared toolbox → per-tab viewers.

## Entry point

`UsageAnalysis:usageAnalysisApp` (`src/package.ts`) is the registered app (`meta.role: adminApp`, url `/`).
It builds a `ViewHandler`, calls `handler.init(...)`, and returns `handler.view`. URL params
(`path`, `date`, `groups`, `packages`, `tags`, `categories`, `projects`) arrive as function params.

## Core pieces

| Piece                     | File                                    | Role                                                                                     |
|---------------------------|-----------------------------------------|------------------------------------------------------------------------------------------|
| `ViewHandler`             | `view-handler.ts`                       | Owns a `DG.MultiView` (the tab strip); registers tabs, routing, per-tab toggles          |
| `UaToolbox`               | `ua-toolbox.ts`                         | The single shared left toolbox (filter accordion) used by every tab                      |
| `UaView`                  | `tabs/ua.ts`                            | Base class each tab extends                                                              |
| Tabs                      | `tabs/*.ts`                             | One file per tab (Overview, Packages, Functions, Events, Clicks, Log, Projects, Metrics) |
| `UaQueryViewer`           | `viewers/abstract/ua-query-viewer.ts`   | Runs a named UA query and builds a `DG.Viewer`                                           |
| `UaFilterableQueryViewer` | `viewers/ua-filterable-query-viewer.ts` | A `UaQueryViewer` that re-runs on every filter change                                    |

## Adding / registering tabs

Tabs are a **fixed list** in `ViewHandler.init()`:

```ts
const viewClasses = [OverviewView, PackagesView, FunctionsView, EventsView,
                     ClicksView, LogView, ProjectsView, MetricsView, VulnerabilitiesView];
```

`VulnerabilitiesView` is toolbox-independent: it loads the published VEX index
(`https://data.datagrok.ai/vex/index.json`) via `grok.dapi.fetchProxy` and drills into the
selected image's per-CVE CSV; the packages/groups filter inputs are hidden on it (like Metrics).

Each is added with `this.view.addView(name, factory, false)`. The factory is **lazy** — a tab's
`tryToInitViewers()` (→ `initViewers()`) runs only when the tab is first shown, so its queries
don't execute until clicked. To add a tab: create a `UaView` subclass in `tabs/`, then add it to this array.

## A tab (`UaView` subclass)

- `name` — tab label (must match `${urlTab}View` for URL routing).
- `rout` — optional sub-route (e.g. Packages flips `/Usage` ↔ `/InstallationTime` in `switchRout()`).
- `viewers: UaQueryViewer[]` — built in `initViewers()`, appended to `this.root`.
- Waits on `_toolboxReady` so viewers never build before the shared toolbox exists.

## Filter / data flow

```
Toolbox "Apply" → filterStream.next(UaFilter)
  → UaFilterableQueryViewer subscription → reload(filter)
    → if activated: reloadViewer()
      → grok.functions.call('UsageAnalysis:<queryName>', filter)   // queries/*.sql
      → createViewer(dataFrame)  → mounted in the tab
```

- `filterStream` is a `BehaviorSubject<UaFilter>` on the toolbox; **one stream, all tabs subscribe**.
- Viewers only re-query when `activated` — hidden tabs stay idle until visited.
- `UaQueryViewer` applies shared formatting (count format, per-user color hashing) before `createViewer`.

## Shared toolbox (`UaToolbox`)

Mounted via `this.view.toolbox = toolbox.rootAccordion.root`; every tab gets the same instance through
`setToolbox()`. It's a `DG.Accordion` with one **Filters** pane: `Date` + choice inputs
(`groups`, `packages`, `tags`, `packagesCategories`, `projects`, each a `ChoiceInput*` from `src/elements/`)
+ an **Apply** button.

`onTabChanged` (in `ViewHandler` and the toolbox) toggles which inputs are visible per tab — e.g. categories
only on Packages, tags only on Functions, projects only on Projects, packages/groups hidden on Metrics —
lazily activates the tab's viewers on first visit, and updates the URL path.

## Drilldown

The Packages context panel (`showSelectionContextPanel`) has **Details** buttons that fill the toolbox's
read-only drilldown fields, reload a target tab's viewer with a derived filter, `changeTab(...)`, and set
`uaToolbox.drilldown`. While drilled down the toolbox swaps the Filters form for `formDD` (read-only summary
+ **Close** and **🠔 back**). `exitDrilldown()` restores the form and reloads the original viewers.

## Routing

`setUrlParam` / `updatePath` keep `view.path` as `/<tab><rout>?<params>` (lowercased). On load, `init()`
parses the first path segment to pick the starting tab (default `Overview`) and applies incoming filter params.
