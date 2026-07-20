# Release-Readiness Dashboard (`/release`)

A single, focused screen that answers "is the next version good to ship?". Lives in the
UsageAnalysis package as a **second app entry** next to the 10-tab Usage Analysis admin app.

## Entry point & dev-only gating

`UsageAnalysis:releaseDashboardApp` (`src/package.ts`) — `@grok.decorators.app`, url `/release`,
`meta.role: adminApp`.

It is **dev-only**: `isDevHost()` (`src/release/data.ts`) checks `grok.dapi.root` (falling back to
`window.location.hostname`) against `dev.datagrok.ai`. On any other server the app warns and returns
`null`, so it is available and usable only on dev. (Function metadata is static, so the launcher entry
still exists off-dev; the runtime guard makes it a no-op there. To hide it entirely elsewhere, publish
UsageAnalysis with this app only to dev.)

## Structure

`ReleaseHandler` (`src/release/release-handler.ts`) owns a `DG.MultiView` with five tabs. All tabs are
toolbox-independent, so — unlike the main `ViewHandler` — it builds **no shared `UaToolbox`**. Each tab
extends `UaView`; the handler calls `markToolboxReady()` (added to `src/tabs/ua.ts`) to unblock
`initViewers()` without a toolbox, then registers a lazy factory (`addView(name, factory, false)`).

| Tab            | File                        | Source                                                      |
|----------------|-----------------------------|-------------------------------------------------------------|
| Overview       | `src/release/overview.ts`   | Fetches all domains in parallel; status cards + alerts       |
| Tests          | `src/release/tests.ts`      | `UsageAnalysis:ReleaseTests` query (`queries/release_tests.sql`) — unit + integration + package (non-manual) |
| Manual         | `src/release/manual.ts`     | `UsageAnalysis:ReleaseManualTests` query — Test Track batches |
| Stress         | `src/tabs/stress-tests.ts`  | **reused**; `StressTestsSummary`/`StressTestsRaw`            |
| Vulnerabilities| `src/tabs/vulnerabilities.ts`| **reused**; VEX `index.json` — core shown as bleeding-edge only (`toDashboardImages`) |
| Tickets        | `src/release/tickets.ts`    | `JiraConnect:GetJiraTicketsByFilter` (fixVersion)           |

Shared, pure helpers live in `src/release/data.ts` (host check, `ReleaseTests` fetch + client-side
pivot, alert computation, `stressRegression`) so Overview and the tabs never duplicate logic.

## Global environment (instance) picker

`ReleaseHandler` owns a single **Environment** choice input (`ENV_CHOICES` = `dev` / `release` /
`public` / `release-ec2`, from `data.ts`) placed on the MultiView **ribbon**. Because `MultiView`
overwrites its ribbon with the active child's on every tab switch — and the child tabs keep their
controls in their own in-view header — the handler re-asserts the picker in the `onTabChanged`
handler (and once after construction), so it stays visible on every page.

The selection is shared through `ReleaseContext` (a `BehaviorSubject<string>` in `data.ts`). The two
**env-scoped** tabs — **Overview** and **Tests** — read `ctx.env.value` when fetching `ReleaseTests`
and subscribe to it to re-fetch on change (guarded by `initialized` so a lazy tab loads once, with the
current env). The other tabs have no per-instance dimension in the DB, so the picker is a no-op there:
Stress (`stress_tests` has no `instance` column — runs in a dedicated job), Manual (Test Track
batches), Vulnerabilities (VEX docker images), Tickets (Jira fixVersion).

## Tests tab

`ReleaseTests` returns one row per test per build (last N builds on a chosen instance) with
per-build `status`/`duration`/`flaking` plus `package`/`owner`. Like the correctness `TestsDashboard`,
it collapses to the **latest CI/CD run per (test, build)** (`latest_runs` CTE, `rn = 1`) — so retried-
then-passed tests aren't counted as failed and non-CI (local/ad-hoc) runs are excluded. Note the failing
count spans **all** exercised suites, including the `PlaywrightPublic` E2E suite; a single package/unit
CI job (e.g. a Jenkins `test-build`) reports fewer failures because it doesn't include E2E. `fetchReleaseTests()` pivots it
client-side (`groupBy(['package','test','owner']).pivot('build_index')`) into one row per test, then
derives `failing`/`stable`/`flaky`/`slower`/`needs_attention` and joins the **Autotests** sticky-meta
schema (`ignore?`/`ignoreReason`/`lastResolved`).

- **Split by package:** a per-package summary grid (pass/fail/skip/not-run + pass rate) over the detail grid.
- **Sparklines:** two grid summary columns via `grid.columns.add({cellType:'sparkline'})` +
  `settings.columnNames` — success (per-build 1/0) and duration (per-build ms). `build_index` is
  ascending, so sparklines read oldest→newest left-to-right.
- **Mute (global):** editing `ignore?`/`ignoreReason` persists back through
  `grok.dapi.stickyMeta.setAllValues(Autotests, testColumn, values)`.
- **Mute (version-bound):** a per-row 🔔/🔕 icon column mutes a test for the *current release version*
  (`versionFromBuild()` of the latest build, e.g. `1.28.0`). The per-test list of muted versions lives in
  its own `ReleaseMutes` sticky-meta schema (`mutedVersions`, `;`-separated) so it never touches the shared
  `Autotests` schema. `muted` (= global `ignore?` OR muted-for-this-version) gates `needs_attention` and
  every failing/flaky count; clicking the icon toggles the version, updates the counts, and persists.
- **Drilldown:** selecting a row shows the latest raw output + a Jenkins link in the context panel.

## Launch

Deep link (dev only): `https://dev.datagrok.ai/apps/usage/release` (per-tab
`/apps/usage/release/<tab>`). The tab appends its name to the app path.

## Known limitations / follow-ups

- **Tickets tab is blocked by a pre-existing JiraConnect bug.** `JiraConnect`'s `loadIssues`
  (`public/packages/JiraConnect/src/api/data.ts`) calls Atlassian's `/rest/api/3/search`, which
  Atlassian **removed** (HTTP 410 — "migrate to `/rest/api/3/search/jql`"). Until JiraConnect is
  migrated to the token-paginated `/rest/api/3/search/jql` endpoint, `GetJiraTicketsByFilter` returns
  nothing and the Tickets tab / Overview card show `n/a` (they degrade gracefully). Jira credentials
  *are* configured on dev (the 410 is reached, not an auth error). `GetJiraTicketsByFilter` also builds
  equality-only JQL, so the fix version is passed quoted (`fixVersion = "<next>"`).
- **Jenkins deep-link** is the `test-build` job page, not the specific build — the `builds` table
  stores no build URL. Proper fix: persist `build_url` in test ingestion
  (`core/server/datlas/lib/src/services/action_logger_service.dart`) and surface it here.
- **`ReleaseTests` shows only tests that actually ran** in the selected builds (`active_tests` CTE) —
  this filters out never-run/malformed `tests.package` rows (Jupyter tracebacks, `{Average time…}`,
  etc.). On a partial instance like `dev`, packages not exercised by the latest build show as
  "not run"; use the `release` instance for a full release-candidate picture.
