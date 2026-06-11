# Metrics tab (Usage Analysis)

## Summary

The **Metrics** tab is the operational-health dashboard of the UsageAnalysis app. It surfaces
live PostgreSQL health, object-storage usage, data-disk free space, request latency, error/session
counts, and `pg_stat_statements` query rankings as a grid of color-banded cards plus several
drill-down panels. Database health comes from cached `Metrics*` SQL queries against
`System:Datagrok`; storage usage comes from a server-maintained snapshot exposed via the JS API
(`grok.dapi.info.getStorageStats()`); disk free space comes from the `DiskStats` server function
(free space on the server's data volume). Every loader is wrapped in `safeCall` so a single failed
source degrades to an "unavailable" card rather than breaking the view.

Implemented in [`metrics.ts`](./metrics.ts) as `class MetricsView extends UaView`.

## Architecture

### Key files

| File                                                                     | Description                                                                                                                   |
|--------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| [`metrics.ts`](./metrics.ts)                                             | `MetricsView` — the whole tab: cards, panels, refresh flow, thresholds, tooltips, snapshot email.                             |
| [`../../queries/metrics.sql`](../../queries/metrics.sql)                 | All `Metrics*` SQL queries (db stats, table health, connections, pg_stat_statements, http routes, latency, errors, sessions). |
| [`../package-api.ts`](../package-api.ts)                                 | Auto-generated typed wrappers (`queries.metrics*`) for the SQL queries; do not edit.                                          |
| [`../ua-toolbox.ts`](../ua-toolbox.ts)                                   | `UaToolbox` — supplies the date filter (`getFilter().date`) and `filterStream`.                                               |
| [`./ua.ts`](./ua.ts)                                                     | `UaView` base class.                                                                                                          |
| [`../../../../../js-api/src/dapi.ts`](../../../../../js-api/src/dapi.ts) | `InfoDataSource.getStorageStats()` — JS API surface for the storage snapshot.                                                 |

### Key classes and methods

`MetricsView` (`metrics.ts`):

- `initViewers()` — builds the layout: header (`as of` label + Refresh / Share buttons), two
  rows of cards (Database, Table health, Storage, Connections, Disk free / Errors, Latency,
  Sessions), the Queries panel, the Largest-tables + Table-health panels, the HTTP-routes panel.
  Subscribes `uaToolbox.filterStream → refreshWindowCards()` and calls `refresh()` once.
- `refresh()` — guards against re-entry with `refreshing`, then runs all loaders in parallel via
  `Promise.all` (db, table-health summary, connections, queries, largest tables, table health,
  window cards, storage, disk).
- `refreshWindowCards()` — re-runs only the time-window-dependent cards (errors, sessions,
  latency, http routes); called both on full refresh and on every date-filter change.
- `safeCall(fn, label)` — static wrapper: `try { await fn() } catch { console.warn; return null }`.
  Every loader treats `null` as "unavailable" and renders a neutral `info` card or a
  `degradedMessage`.
- `loadDbStats`, `loadTableHealthSummary`, `loadConnections`, `loadErrors`, `loadSessions`,
  `loadLatency`, `loadStorage`, `loadDisk` — one per card; each sets value/sub-text/color and
  binds a rich tooltip built by a matching `build*Tooltip` static.
- `loadQueries`, `loadLargestTables`, `loadTableHealth`, `loadHttpRoutes` — render grids into
  panel hosts (`DG.Viewer.grid`) with per-cell color banding via `onCellPrepare`.
- `detectPssMode()` / `pssFn()` — picks modern vs legacy `pg_stat_statements` query variants
  (see Data flow).
- `openFullView()` / `openFullHttpRoutes()` — re-run a query with `FULL_VIEW_LIMIT` (100000) and
  open the result as a standalone `TableView` ("Add to workspace").
- `sendToDatagrok()` / `collectSnapshotAttachments()` — build CSV/JSON attachments for every
  source and open the e-mail composer ("Share...").
- `confirmResetPgStats()` — confirmation dialog that calls `MetricsResetPgStatStatements`.

### Color banding and thresholds

`THRESH` (top of `metrics.ts`) holds the green/orange/red cut-offs for every metric.
`thresholdBand(v, t, higherIsBetter)` returns `'green' | 'orange' | 'red'`; `thresholdColor()`
maps that to `DG.Color.success` / `0xFFFFA24A` / `DG.Color.failure` for grid cells.
`deltaColor()` colors the window cards by trend (more = red, fewer = green, equal = orange).
The card sub-text gets the band as a CSS class (`ua-metrics-green` / `-orange` / `-red`).

### Data flow

1. **Cards / grids → SQL.** Each loader calls a `queries.metrics*` wrapper (or
   `grok.data.query('UsageAnalysis:Metrics*', …)`), which runs a `Metrics*` query in
   `queries/metrics.sql` against `System:Datagrok` (read-only role). Most queries carry
   `--meta.cache: all` + `--meta.cache.invalidateOn: */5 * * * *`, so results are server-cached
   for ~5 minutes; the dashboard `as of` label is just the client-side fetch time.
2. **Time-window queries.** `MetricsErrorsCount`, `MetricsSessionsCount`, `MetricsLatency`,
   `MetricsHttpRoutes` take the toolbox `date` filter and internally derive the current window
   `[min_date, max_date]` and an equal-length previous window for the `±delta vs prev` comparison.
   These carry the same `--meta.cache: all` + `*/5 * * * *` invalidation as the rest, so they are
   also server-cached for ~5 minutes, keyed by the date window.
3. **pg_stat_statements modern/legacy split.** `detectPssMode()` queries
   `MetricsPgStatStatementsVersion` once, parses `extversion`, and chooses `'modern'`
   (≥ 1.8, columns `total_exec_time` / `mean_exec_time`) or `'legacy'` (≤ 1.7, columns
   `total_time` / `mean_time`) — the column names changed in extension v1.8. `unavailable` means
   the extension is missing or `System:Datagrok` lacks `pg_read_all_stats`, and the Queries panel
   shows a degraded message. `PSS_QUERIES` maps each `QueriesMode` (`slowest` / `most-called` /
   `worst-hit`) to its modern + legacy wrapper and full-view query name; `pssFn()` dispatches.
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
| `grok.dapi.info.getStorageStats()`                                       | `InfoDataSource` in [`dapi.ts`](../../../../../js-api/src/dapi.ts); `Promise<{[key: string]: any}>`, the hourly storage snapshot. |
| `grok.functions.call('DiskStats')`                                       | Returns a JSON string parsed to `DiskStat`.                                                                                       |
| `grok.data.query('UsageAnalysis:Metrics*', params)` / `queries.metrics*` | Typed wrappers for the `Metrics*` SQL queries.                                                                                    |
| `grok.dapi.admin.getReportEmail()` / `grok.dapi.admin.sendEmail(...)`    | Used by "Share..." to address and send the snapshot e-mail.                                                                       |

### Metrics* SQL queries (`queries/metrics.sql`)

| Query                                   | Returns / consumed by                                                                                                                       |
|-----------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| `MetricsDbStats`                        | DB size, cache-hit %, stats_reset, and up to 5 cache-miss "offender" tables → **Database** card + tooltip.                                  |
| `MetricsTableHealthSummary`             | Count of unhealthy tables (≥10K live rows, >40% dead) + max dead % + top offenders → **Table health** card + tooltip.                       |
| `MetricsTableHealth`                    | Per-table dead % and last vacuum (≥1K live rows) → **Table health** panel grid + full view.                                                 |
| `MetricsLargestTables`                  | Top tables by total relation size (total, index, #rows, total_bytes) → **Largest tables** panel grid + full view.                           |
| `MetricsConnections`                    | `pg_stat_activity` summary: total / active / idle / idle-in-xact / waiting-on-lock + oldest idle-xact age → **Connections** card + tooltip. |
| `MetricsConnectionsOffenders`           | Per-session offenders (idle-in-xact, long active, lock-blockers) → **Connections** tooltip table (lazy, 500 ms after card).                 |
| `MetricsHttpRoutes`                     | Per `method + route`: count, p50/p95/p99 ms, error % (status ≥ 400) → **HTTP routes** panel grid + full view.                               |
| `MetricsErrorsCount`                    | Error-event count for the window vs the previous window (events whose `event_types.source = 'error'`) → **Errors** card.                    |
| `MetricsLatency`                        | http-request p50/p95/p99 (now) + p95 (prev) + request counts, from `event_parameter_values.ms` → **Latency** card.                          |
| `MetricsSessionsCount`                  | Distinct `users_sessions` started in the window vs previous → **Sessions** card.                                                            |
| `MetricsTopSlowestQueries` / `…Pg12`    | pg_stat_statements ranked by mean exec time → **Queries** grid (mode "slowest").                                                            |
| `MetricsTopMostCalledQueries` / `…Pg12` | Ranked by call count → **Queries** grid (mode "most called").                                                                               |
| `MetricsWorstCacheHitQueries` / `…Pg12` | Lowest shared-buffer hit ratio (≥10 calls) → **Queries** grid (mode "worst cache hit").                                                     |
| `MetricsPgStatStatementsVersion`        | `pg_extension.extversion` for `pg_stat_statements` → drives `detectPssMode`.                                                                |
| `MetricsResetPgStatStatements`          | Calls `pg_stat_statements_reset()` on `System:DatagrokAdmin` → "Reset stats" menu action.                                                   |

**modern vs Pg12 variants.** The three pg_stat_statements rankings exist in two forms because the
extension renamed its timing columns in v1.8: `total_time → total_exec_time` and
`mean_time → mean_exec_time` (confirmed by diffing the `…` vs `…Pg12` queries in `metrics.sql`).
`detectPssMode()` picks the variant from the detected `extversion`. All ranking queries exclude
their own monitoring SQL (`pg_stat_*`, `pg_database_size`) and transaction-control statements.

Most read queries run on `System:Datagrok` (read-only); `MetricsResetPgStatStatements` runs on
`System:DatagrokAdmin` because the reset is privileged.

### Share snapshot (email)

The header **Share...** button (`sendToDatagrok()`) opens the platform e-mail composer
(`ui.composeEmail`) with one attachment per data source. Recipient is pre-filled from
`grok.dapi.admin.getReportEmail()`; `onSend` calls `grok.dapi.admin.sendEmail(...)`.
`collectSnapshotAttachments()` builds the files via `safeCall` (a failed source is skipped):
a CSV per metrics query (`EMAIL_LIMIT` = 100 rows), the three `pg_stat_statements` CSVs when
available, plus `storage.json` and `disk.json`.

## Usage

- **UI.** Open the **Usage Analysis** application → **Metrics** tab. Use the date filter in the
  toolbox to change the window for the Errors / Latency / Sessions / HTTP-routes cards (the DB /
  storage / disk / pg_stat_statements cards are point-in-time and ignore the filter). The
  **Refresh** button re-runs everything; **Share...** opens an e-mail composer with one
  CSV/JSON attachment per source. The Queries panel toggles between slowest / most-called /
  worst-cache-hit; its "⋯" menu offers "Add to workspace" and "Reset stats". Each panel's
  **+** icon opens the full (100000-row) result as a standalone table view. Hover any card for a
  diagnostic tooltip with thresholds and remediation hints.

Notes:
- `DiskStats` returns `null`/empty on servers without `df` (e.g. Windows).
- `getStorageStats()` returns `{}` until the first snapshot has been produced after the server boots.
- `pg_stat_statements` cards require the extension installed and `pg_read_all_stats` granted to the
  `System:Datagrok` role; otherwise the Queries panel shows a degraded message.
