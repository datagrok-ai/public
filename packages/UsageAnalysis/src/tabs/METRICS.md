# Metrics tab (Usage Analysis)

## Summary

The **Metrics** tab is the operational-health dashboard of the UsageAnalysis app. It surfaces
live PostgreSQL health, object-storage usage, data-disk free space, request latency, the
function-call queue, error/session counts, and `pg_stat_statements` query rankings as a grid of
color-banded cards plus several drill-down panels. Request latency, per-route statistics, the queue,
connections, database size and the `pg_stat_statements` rankings come from one call to
`grok.dapi.admin.getMetrics()` (the admin-only `/api/admin/metrics` endpoint); table health, largest
tables, cache-miss tables, session offenders, errors and sessions come from cached `Metrics*` SQL
queries against `System:Datagrok`; storage usage comes from a server-maintained snapshot exposed via
the JS API (`grok.dapi.info.getStorageStats()`); disk free space comes from the `DiskStats` server
function (free space on the server's data volume); data growth (cards row 3 and the Data growth panel) comes
from `grok.dapi.admin.getGrowthMetrics()` (`/api/admin/metrics/growth`, samples the server collects on a
schedule). Every loader is wrapped in `safeCall` so a single failed source degrades to an "unavailable" card
rather than breaking the view.

Implemented in [`metrics.ts`](./metrics.ts) as `class MetricsView extends UaView`.

## Architecture

### Key files

| File                                                                     | Description                                                                                                                   |
|--------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| [`metrics.ts`](./metrics.ts)                                             | `MetricsView` — the whole tab: cards, panels, refresh flow, thresholds, tooltips, snapshot email.                             |
| [`../../queries/metrics.sql`](../../queries/metrics.sql)                 | All `Metrics*` SQL queries (table health, largest tables, cache-miss tables, session offenders, errors, sessions, stats reset). |
| [`../package-api.ts`](../package-api.ts)                                 | Auto-generated typed wrappers (`queries.metrics*`) for the SQL queries; do not edit.                                          |
| [`../ua-toolbox.ts`](../ua-toolbox.ts)                                   | `UaToolbox` — supplies the date filter (`getFilter().date`) and `filterStream`.                                               |
| [`./ua.ts`](./ua.ts)                                                     | `UaView` base class.                                                                                                          |
| [`../../../../../js-api/src/dapi.ts`](../../../../../js-api/src/dapi.ts) | `AdminDataSource.getMetrics()` (server metrics, `ServerMetrics` type) and `InfoDataSource.getStorageStats()` (storage snapshot). |

### Key classes and methods

`MetricsView` (`metrics.ts`):

- `initViewers()` — builds the layout: header (`as of` label + Refresh button), two
  rows of cards (Database, Table health, Storage, Connections, Queue, Disk free / Errors, Latency,
  Sessions), the Queries panel, the Largest-tables + Table-health panels, the HTTP-routes panel.
  Subscribes `uaToolbox.filterStream → refreshWindowCards()` and calls `refresh()` once.
- `refresh()` — guards against re-entry with `refreshing`, then runs all loaders in parallel via
  `Promise.all` (table-health summary, largest tables, table health, window cards, storage, disk).
- `refreshWindowCards()` — re-runs the time-window-dependent sources (errors, sessions, and
  `loadMetrics()` — the endpoint call, which also re-renders the point-in-time cards it feeds);
  called both on full refresh and on every date-filter change.
- `fetchMetrics(limit)` / `loadMetrics()` — `fetchMetrics` calls
  `grok.dapi.admin.getMetrics({date, limit})` with the toolbox date; `loadMetrics` calls it once
  with `DASHBOARD_LIMIT`, keeps the response in `this.metrics`, and renders everything it feeds:
  Latency, HTTP routes, Queue, Connections, Queries, Database.
- `safeCall(fn, label)` — static wrapper: `try { await fn() } catch { console.warn; return null }`.
  Every loader treats `null` as "unavailable" and renders a neutral `info` card or a
  `degradedMessage`.
- `loadDbStats`, `loadTableHealthSummary`, `loadConnections`, `loadQueue`, `loadErrors`,
  `loadSessions`, `loadLatency`, `loadStorage`, `loadDisk` — one per card; each sets
  value/sub-text/color and binds a rich tooltip built by a matching `build*Tooltip` static.
- `loadQueries`, `loadLargestTables`, `loadTableHealth`, `loadHttpRoutes` — render grids into
  panel hosts (`DG.Viewer.grid`) with per-cell color banding via `onCellPrepare`.
- `frame()` / `routesFrame()` / `statementsFrame()` — turn the endpoint's arrays into
  `DG.DataFrame`s (`DG.DataFrame.fromObjects`) with the column names the grids and CSVs use
  (`route`, `count`, `p50`, `p95`, `p99`, `err_pct`; `query`, `calls`, `total_ms`, `mean_ms`,
  `hit_pct`).
- `openFullView()` — re-runs a source with `FULL_VIEW_LIMIT` (100000) — a `Metrics*` query, or
  the endpoint via `fetchMetrics` for HTTP routes and the `pg_stat_statements` rankings — and opens
  the result as a standalone `TableView` ("Add to workspace").
- `collectSnapshotAttachments()` — CSV attachments of every source, for the tab's CSV report (`exportFiles`).
- `openRouteFailures()` — a click on a route's `err_pct` opens its failed requests in the window
  (`MetricsRouteFailures`: time, status, ms, user, session and the nearest error or warning of that session within
  5 s) as a table view, with the current request's log docked below (`MetricsRequestLog`: the session's events
  within 5 s of it).
- `confirmResetPgStats()` — confirmation dialog that calls `MetricsResetPgStatStatements`.

### Color banding and thresholds

`THRESH` (top of `metrics.ts`) holds the green/orange/red cut-offs for every metric.
`thresholdBand(v, t, higherIsBetter)` returns `'green' | 'orange' | 'red'`; `thresholdColor()`
maps that to `DG.Color.success` / `0xFFFFA24A` / `DG.Color.failure` for grid cells.
`deltaColor()` colors the window cards by trend (more = red, fewer = green, equal = orange).
The card sub-text gets the band as a CSS class (`ua-metrics-green` / `-orange` / `-red`).

### Data flow

1. **Endpoint-fed cards and grids.** `loadMetrics()` calls `grok.dapi.admin.getMetrics({date, limit})`
   once per refresh (and once per date-filter change) and keeps the response in `this.metrics`. The
   server resolves the toolbox `date` pattern into `[window.start, window.end]` plus an equal-length
   previous window (`window.previousStart`) for the `±delta vs prev` comparison, and returns
   `http.now` / `http.previous` (Latency card), `http.routes` ordered by p95 (HTTP routes panel),
   `queue` (Queue card), `database.sizeBytes` / `cacheHitPct` / `statsReset` (Database card),
   `database.connections` (Connections card) and `database.statements` (Queries panel — all three
   rankings arrive in one response, so the mode toggle re-renders without a round trip).
2. **SQL-fed cards and grids.** The remaining loaders call a `queries.metrics*` wrapper, which runs
   a `Metrics*` query in `queries/metrics.sql` against `System:Datagrok` (read-only role). Most
   queries carry `--meta.cache: all` + `--meta.cache.invalidateOn: */5 * * * *`, so results are
   server-cached for ~5 minutes; the dashboard `as of` label is just the client-side fetch time.
   `MetricsErrorsCount` and `MetricsSessionsCount` take the toolbox `date` filter and internally
   derive the current window `[min_date, max_date]` and an equal-length previous window.
3. **pg_stat_statements.** The server reads the extension itself (picking the `total_exec_time` /
   `total_time` column names by its version) and reports `statements.available = false` when it is
   missing or unreadable; the Queries panel then shows a degraded message. "Reset stats" still runs
   the `MetricsResetPgStatStatements` query on `System:DatagrokAdmin` and refetches the metrics.
4. **Storage card → server snapshot.** `loadStorage()` calls `grok.dapi.info.getStorageStats()`,
   which returns a **pre-computed snapshot** maintained by the server — it does **not** scan storage
   on demand (see below). `loadStorageCoords()` separately reads bucket/region from the
   `System:AppData` connection parameters for the tooltip header.
5. **Disk card → server function.** `loadDisk()` calls `grok.functions.call('DiskStats')`, gets a
   JSON string, and `JSON.parse`s it into `DiskStat`.

### Storage snapshot (`getStorageStats`)

The **Storage** card calls `grok.dapi.info.getStorageStats()` — a public JS-API call
(`InfoDataSource` in [`dapi.ts`](../../../../../js-api/src/dapi.ts)) that returns a **server-maintained
snapshot**, not an on-demand scan. Observable contract:

- The server refreshes the snapshot on a background schedule (roughly hourly); the call is a cheap
  read of the cached result, so it is safe to invoke on every dashboard refresh.
- It returns an **empty object** (`{}`) until the first snapshot has been produced after the server
  starts — `loadStorage()` treats that as "unavailable".
- The snapshot rolls up the platform's object storage regardless of backend (local file share, S3,
  GCS, Azure Blob): total bytes, total object count, and the **top prefixes** by size.
- `truncated` is `true` if any prefix could not be fully measured; the card then prefixes a `>` and
  the tooltip shows `>?` for the affected rows so partial totals aren't mistaken for exact ones.
- The **AppData** bucket/region in the tooltip header are *not* part of the snapshot — the client
  reads them from the `System:AppData` connection's `bucket` / `region` parameters via
  `loadStorageCoords()`.

Returned shape, consumed by `buildStorageTooltip`:

```jsonc
{
  "type": "S3" | "Local" | "GCS" | ...,   // storage backend
  "root": "<root path / bucket>",
  "totalBytes": 123456789,
  "objectCount": 4567,
  "truncated": false,                       // a prefix scan failed → totals incomplete
  "collectedAt": "2026-06-10T08:00:00.000Z",// snapshot time (tooltip "as of")
  "durationMs": 8421,
  "topPrefixes": [
    {"name": "Demo", "bytes": 99999, "objectCount": 1200, "truncated": false},
    ...                                     // largest prefixes by bytes
  ]
}
```

### Data growth (`getGrowthMetrics`)

The third row of cards and the **Data growth** panel show how the platform's data grows, from
`grok.dapi.admin.getGrowthMetrics({days})` (`ViewTelemetry`; `GrowthMetrics` in
[`dapi.ts`](../../../../../js-api/src/dapi.ts)). The server samples the database (table sizes and rows
hourly, categories such as events by source daily), the storage (by location and nature every 6 hours),
the servers' disks, and what the garbage collector deleted and has pending, into `metric_samples`; the
call answers one row per key with its value now, its change over 1, 7 and 30 days, a daily trend of
`days` and the metric alerts that cover it. `loadGrowth()` fetches every area once per refresh (and when
**Days** changes); the area toggle and **Only growing** filter the response without a round trip. The
data is point-in-time: the toolbox date filter does not apply.

| Card | Shows |
|---|---|
| DB growth 7d | Δ7d of the database (`(total)` row, `bytes.d7`), its % of the size a week ago banded by `THRESH.growthPct` (10 / 25 %); sub: the size now; tooltip: the 5 tables that grew most |
| Storage growth 7d | the same of the shared storage (`shared` row); sub: backend type and size; tooltip: the 5 natures that grew most |
| GC deleted | the sum of `gc.deletedLast` of the `gc` rows; sub: the last run, or "GC disabled" in orange; tooltip: deleted in 7 days by key |
| GC pending | the sum of `gc.pending`; sub: the key with the most; tooltip: pending by key |
| Growth alerts | how many (alert, key) pairs hold; red when one of their alerts is open or acknowledged, orange when they hold without an open alert (`info` severity, muted), green when none hold. The response carries no severity, so the color reads the alert state |

Before the first collection (`collectedAt` null) every card says "collecting…" and the panel says so; a
card whose 7-day delta has no sample that old shows "—" with "collecting…".

Panel: one row per response row of the selected area, sorted by Δ7d bytes then size, columns `name`
(`table`, `table · dim=value`, `location · nature[/sub]`), `size`, `count` (`~` when estimated), `d7_bytes`,
`d7_pct` (colored by `THRESH.growthPct`), `trend` (a sparkline over the hidden `t0..t{days-1}` columns),
then `written_7d` (Δ7d count + GC deleted in 7 days: the real write rate), `gc_7d`, `gc_pending` for the
database and GC areas, or `free_bytes`, `used_pct` for storage (disk rows), and `alert`: a red dot when a
metric alert holds, its tooltip naming the alerts; a click opens the first one's problem (`/alerts/<problemId>`) —
the id is only in the response for `ManageAlerts` holders, so for others the dot is only a tooltip. **+** opens the full frame (every
column, `t*` included) as a table view. The CSV report adds `growth.csv`. Without `ViewTelemetry` the call
rejects with 403 and the cards and panel show "unavailable".

### Disk stats (`DiskStats`)

The **Disk free** card calls `grok.functions.call('DiskStats')`, a server-side cached function
(refreshed every few minutes) that returns a JSON string the tab parses into `DiskStat`
(`{path, mount, totalBytes, usedBytes, freeBytes, usedPct}`). Observable contract:

- Reports **free space on the server's data volume** — the disk where the server stores tables,
  files, and caches — not the OS root.
- Returns `null` / empty on servers without a `df`-style tool (e.g. Windows); the card then shows
  "unavailable".

The card shows `formatBytes(freeBytes)` and `"<usedPct>% used"`. Color comes from `diskColor()`:

| Band   | Condition (`THRESH`)                      |
|--------|-------------------------------------------|
| red    | `usedPct ≥ 90` **or** `freeBytes < 2 GiB` |
| orange | `usedPct ≥ 75` **or** `freeBytes < 5 GiB` |
| green  | otherwise                                 |

Because either a high used-percent or a low absolute free-bytes triggers escalation, a large disk
that is only 70% full but has < 2 GiB free still goes red.

### JS API

| Surface                                                                  | Description                                                                                                                       |
|--------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `grok.dapi.admin.getMetrics({date, limit})`                              | `AdminDataSource` in [`dapi.ts`](../../../../../js-api/src/dapi.ts); `Promise<ServerMetrics>` — `window`, `http` (now / previous / routes), `queue`, `database` (size, cache hit, stats reset, connections, statements). Admin only. |
| `grok.dapi.admin.getGrowthMetrics({area, days})`                         | `AdminDataSource` in [`dapi.ts`](../../../../../js-api/src/dapi.ts); `Promise<GrowthMetrics>` — `rows` (area, key, `bytes`/`count` series with `d1`/`d7`/`d30` and `trend`, `gc`, `alerts`), `gc`, `storage`, `schedule`. `ViewTelemetry`. |
| `grok.dapi.info.getStorageStats()`                                       | `InfoDataSource` in [`dapi.ts`](../../../../../js-api/src/dapi.ts); `Promise<{[key: string]: any}>`, the hourly storage snapshot. |
| `grok.functions.call('DiskStats')`                                       | Returns a JSON string parsed to `DiskStat`.                                                                                       |
| `queries.metrics*`                                                       | Typed wrappers for the `Metrics*` SQL queries.                                                                                    |

### Metrics* SQL queries (`queries/metrics.sql`)

| Query                                   | Returns / consumed by                                                                                                                       |
|-----------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `MetricsCacheMissTables`                | Up to 5 user tables with ≥1K disk reads and a hit ratio < 95% → **Database** tooltip table (the card itself is endpoint-fed).               |
| `MetricsTableHealthSummary`             | Count of unhealthy tables (≥10K live rows, >40% dead) + max dead % + top offenders → **Table health** card + tooltip.                       |
| `MetricsTableHealth`                    | Per-table dead % and last vacuum (≥1K live rows) → **Table health** panel grid + full view.                                                 |
| `MetricsLargestTables`                  | Top tables by total relation size (total, index, #rows, total_bytes) → **Largest tables** panel grid + full view.                           |
| `MetricsConnectionsOffenders`           | Per-session offenders (idle-in-xact, long active, lock-blockers) → **Connections** tooltip table (lazy, 500 ms after card).                 |
| `MetricsErrorsCount`                    | Error-event count for the window vs the previous window (events whose `event_types.source = 'error'`) → **Errors** card.                    |
| `MetricsSessionsCount`                  | Distinct `users_sessions` started in the window vs previous → **Sessions** card.                                                            |
| `MetricsResetPgStatStatements`          | Calls `pg_stat_statements_reset()` on `System:DatagrokAdmin` → "Reset stats" menu action.                                                   |

Read queries run on `System:Datagrok` (read-only); `MetricsResetPgStatStatements` runs on
`System:DatagrokAdmin` because the reset is privileged.

### CSV report

`collectSnapshotAttachments()` (the tab's `exportFiles`, used by the Usage Analysis CSV report) calls the
endpoint once with `EMAIL_LIMIT` (100 rows) and builds the files via `safeCall` (a failed source is skipped): a
CSV per remaining metrics query, `db_summary.csv` / `connections.csv` / `queue.csv` / `http_routes.csv` /
`latency.csv` from the endpoint response, the three `pg_stat_statements` CSVs when `statements.available`, plus
`storage.json` and `disk.json`.

## Usage

- **UI.** Open the **Usage Analysis** application → **Metrics** tab. Use the date filter in the
  toolbox to change the window for the Errors / Latency / Sessions / HTTP-routes cards (the DB /
  Connections / Queue / storage / disk / pg_stat_statements cards are point-in-time and ignore the
  filter). The Queue card shows queued + running function calls and turns orange while anything is
  waiting. The
  **Refresh** button re-runs everything. Click a route's error rate in **HTTP routes** to see its failed
  requests, each with its nearest error, and the selected request's log. The Queries panel toggles between slowest / most-called /
  worst-cache-hit; its "⋯" menu offers "Add to workspace" and "Reset stats". Each panel's
  **+** icon opens the full (100000-row) result as a standalone table view. Hover any card for a
  diagnostic tooltip with thresholds and remediation hints.

Notes:
- `DiskStats` returns `null`/empty on servers without `df` (e.g. Windows).
- `getStorageStats()` returns `{}` until the first snapshot has been produced after the server boots.
- `getMetrics()` requires an admin session (403 otherwise); every card it feeds then shows
  "unavailable" and the Queries / HTTP-routes panels a degraded message.
- The `pg_stat_statements` rankings require the extension installed and readable by the server's
  database role; otherwise `statements.available` is false and the Queries panel shows a degraded
  message.
