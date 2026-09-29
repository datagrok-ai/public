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
| Tabs                      | `tabs/*.ts`                             | One file per tab (the list in `ViewHandler.init()` below)                                |
| `UaQueryViewer`           | `viewers/abstract/ua-query-viewer.ts`   | Runs a named UA query and builds a `DG.Viewer`                                           |
| `UaFilterableQueryViewer` | `viewers/ua-filterable-query-viewer.ts` | A `UaQueryViewer` that re-runs on every filter change                                    |

## Adding / registering tabs

Tabs are a **fixed list** in `ViewHandler.init()`:

```ts
const viewClasses = [OverviewView, PackagesView, FunctionsView, EventsView, ClicksView, LogView,
                     SystemActivityView, ErrorsView, CaptureView, TimelineView, ProjectsView, MetricsView, StressView,
                     VulnerabilitiesView];
```

`CaptureView` (tab `Capture`) lists the capture rules (`CaptureRules` over `capture_rules`, the columns of
`grok s capture list`: rule, author, subject, scope, reason, active, events); active rules always show, ended ones
when created within the Date filter. **New rule...** builds the `POST /logging/capture` body and calls the server
function `CaptureRuleAdd`; its Debug flags are the server's (`LoggingPolicy`, `debugFlags` of `GET /logging/policy`)
but `credentials`, and hidden when the policy can't be read; a row's context panel and context menu offer **Stop...** (`CaptureRuleStop`, active rules
only) and **Timeline**. `TimelineView` (tab `Timeline`) shows one action, request, session, report or rule
(`cap-<n>`) in time order from the server function `Timeline`; it routes as `/timeline?<key>=<id>` (the id keeps its
case; the platform hands the query parameters to `usageAnalysisApp`, which passes them on as
`TimelineView.urlParams`; other tabs drop the parameter from their path). The Clicks tab's **Clicks** sub-tab lists
single clicks (`Clicks` query) with their `request_id` (the action id); the row's context menu **Timeline** opens that
action. The server enforces the permissions on all three functions; the packages and groups inputs are hidden on both
tabs.

`ErrorsView` (tab `Errors`) is the platform's errors as data, from the server function `ErrorStats` (the
`GET /errors` query; ViewTelemetry): without Group by, occurrences; with up to three Group by dimensions, one row of
figures per combination (first seen in, trend, and, by signature, the alert state). Its inputs are the toolbox's
**Errors** pane (`UaToolbox.addTabPane`), which replaces the Filters pane while the tab is current. A row's context
panel runs `ErrorStats` narrowed to the row and lists the occurrences (request → Timeline), then the sessions, reports
and alerts of their signatures (`ErrorSessions`, `ErrorReports`, `ErrorAlerts` in `errors_query.sql`). Those three
are package queries on `System:Datagrok`, so they follow the rest of the app's access model rather than
ViewTelemetry: whoever can use that connection (administrators by default) sees them, and anyone else gets an error in
that pane only. **Export** writes the shown table as CSV (a cell starting with `=`, `+`, `-`, `@`, a tab or a carriage
return gets a leading `'`, as in the server's CSV), JSON or Parquet (`Arrow:toParquet`, disabled without Arrow);
**Save as job...** calls
`ErrorsSaveJob` with the shown query (Since only; cron in UTC). The Clicks tab's **Followed by Error** sub-tab
(`ClicksFollowedByError`) counts clicks per element and those an error followed within 5 s whose request id is the
click's action id or `<action>.<n>`; anonymous clicks count users by `anonSession`.

`SystemActivityView` (tab `System Activity`) lists the platform-level audit records datlas writes —
`server-started`, `user-logged-in`, `user-logged-out`, `user-login-failed`, `user-impersonated`,
`impersonation-failed`, `admin-session-started`, `admin-session-ended`, `dev-key-generated`,
`settings-changed`, `log-settings-changed` (`LogAudit` in `grok_shared/lib/src/log_entities.dart`).
`SystemActivitySummary` feeds a per-type timeline and `SystemActivity` the filterable grid; the grid
resolves the record's user from its `user` parameter first (a failed login or a server start has no
session), so the groups filter applies to records with a user and passes the rest through. The
packages input is hidden on it, like Projects. Both queries are uncached: this is a security log.

`VulnerabilitiesView` is toolbox-independent: it loads the published VEX index
(`https://data.datagrok.ai/vex/index.json`) via `grok.dapi.fetchProxy` and drills into the
selected image's per-CVE CSV; the packages/groups filter inputs are hidden on it (like Metrics).

`StressView` (tab `Stress`) is likewise toolbox-independent: `StressTestsSummary` feeds a
median-duration-by-build line chart (split by thread count) and a metrics grid, and
`StressTestsRaw` (latest build) feeds a threads-vs-duration scatter plot colored by pass/fail.

Each is added with `this.view.addView(name, factory, false)`. The factory is **lazy** — a tab's
`tryToInitViewers()` (→ `initViewers()`) runs only when the tab is first shown, so its queries
don't execute until clicked. To add a tab: create a `UaView` subclass in `tabs/`, then add it to this array.

## A tab (`UaView` subclass)

- `name` — tab label; its URL segment is the name without spaces, lowercased (`ViewHandler.urlName`),
  so `System Activity` routes as `/systemactivity`.
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

`setUrlParam` / `updatePath` keep `view.path` as `/<tab><rout>?<params>` (the tab and route lowercased, the
parameters as they are). On load, `init()`
parses the first path segment to pick the starting tab (default `Overview`) and applies incoming filter params.
