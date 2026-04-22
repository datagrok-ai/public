# Data Enrichment — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1.1 | Go to Databases > Postgres > Northwind | 1m 5s | AMBIGUOUS | FAILED | Scenario asks for "Northwind" but the dev server has "NorthwindTest" in `Databases > Postgres` (also under MS SQL and PostgresDart). Located it via Browse tree; `grok.dapi.connections.filter(name like '%orthwind%')` and `.list()` do NOT return it — the filter field name differs and pagination caps at 80. The Playwright replay failed earlier because the `Browse` tab click alone does not expand the Databases subtree. |
| 1.2 | Create a SQL and a visual query for the 'Orders' table | 30s | SKIP | FAILED | A pre-existing visual query (`Agolovko:OrdersVQ`) already exists on the connection — so creation was not re-attempted. Running it returns `FATAL: password authentication failed for user "datagrok"`. Direct `grok.data.db.query('Dbtests:PostgresTest', 'select 1')` fails the same way. All Postgres connections I tried (PostgresTest, PostgreSQLDBTests, PostgreSQLDBTestsCached, PostgreSQLTest, PostgresChemblTest, Chembl:ChemblSql, MolTrack:moltrack, UsageAnalysis:ua_tickets, Apitestsdb:apitest_db) error with auth failure, `SocketException: Connection refused`, or `Failed host lookup: 'database'`. |
| 1.3 | Open Orders and run the created queries | 10s | SKIP | SKIPPED | Blocked by 1.2 — no working Postgres connection. |
| 1.4 | Click the 'customerid' column | 5s | SKIP | SKIPPED | Blocked by 1.3. |
| 1.5 | Context Panel > NorthwindTest > Enrich.. | 5s | SKIP | SKIPPED | Blocked by 1.4. Verified on a CSV-loaded table that the Context Panel never shows an "Enrich" subsection (only Details/Filter/Actions/Colors/Style/Settings/Stats/Plots/Advanced/Dev), so this feature is only visible for DB-backed tables. |
| 1.6 | + Add Enrichment.., add a join on public.Customers, save | 10s | SKIP | SKIPPED | Blocked by 1.5. |
| 1.7 | Apply enrichment — verify columns added, then edit it | 5s | SKIP | SKIPPED | Blocked by 1.6. |
| 2 | Multiple enrichments per column and multiple columns | 5s | SKIP | SKIPPED | Blocked by Section 1. |
| 3 | Persistence across projects, layouts, reuse in other tables | 5s | SKIP | SKIPPED | Blocked by Section 1. |
| 4 | Visibility for different users | 5s | SKIP | SKIPPED | Blocked by Section 1 and by the single-credential environment used for scripted runs. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m 30s |
| grok-browser execution (scenario steps) | 2m 40s |
| Execute via grok-browser (total) | 7m 10s |
| Spec file generation | 1m 15s |
| Spec script execution | 43s |
| **Total scenario run (with model)** | 9m 10s |

## Summary

The scenario is blocked on infrastructure. Dev has a `NorthwindTest` Postgres connection at `db.datagrok.ai:54322 / northwind` with 13 queries attached, but every query and every direct SQL call — across every Postgres connection on the server — fails with `FATAL: password authentication failed for user "datagrok"` (or Connection refused / host lookup failure). Without a working DB-backed table, the "Enrich" Context Panel entry never appears (confirmed on a CSV-loaded table: no Enrich section in the accordion), so Sections 1.3–4 cannot be exercised at all. Total scenario run (with model): 9m 10s.

## Retrospective

### What worked well
- The Browse tree and the Connections list surface three `NorthwindTest` entries (MS SQL, Postgres, PostgresDart) correctly, with distinct `icon-ds-*` icons on each — these can be targeted from tests by filtering parent `[name="icon-ds-postgres"]`.
- Context Panel for the selected connection loads metadata (server/db/port, query count, history, sharing, activity) without error, so the connection entity itself is healthy on the server side.
- The existing `Agolovko:OrdersVQ` visual query would have been a good shortcut to step 1.2 if auth worked.

### What did not work
- All Postgres credentials on dev are broken: `FATAL: password authentication failed for user "datagrok"`. The connection records point at `db.datagrok.ai:54322` with a stored credential that the target Postgres instance no longer accepts.
- PostgresDart connections fail a level earlier: `Failed host lookup: 'database'` — they're configured against a docker-compose hostname that doesn't resolve from dev's network.
- `ChemblSql` / `PostgresChemblTest` fail with `Connection refused` to `db.datagrok.ai:54325`, suggesting that port is not exposed.
- `MolTrack:moltrack` fails with `Unable to get project asset "moltrack"` — the connection record is dangling.
- `grok.dapi.connections.filter('name like "%orthwind%"').list()` returns unrelated entries; the filter seems to match on description not name, and basic pagination caps at 80 rows even when asking for 100 per page and 10 pages. Finding connections programmatically is unreliable.

### Suggestions for the platform
- Fix the stored credential on the `datagrok` Postgres role so `NorthwindTest` (and all sibling Postgres connections on dev) can authenticate. Add a scheduled smoke test that runs `select 1` against every Postgres connection and alerts on auth failure — the connections are the backbone of half the test scenarios.
- Rework PostgresDart connections to use reachable hostnames in the dev environment, or mark them as not-for-dev so smoke tests don't try to use them.
- Let `grok.dapi.connections.filter('name like ...')` actually filter on `name`, or add `nameContains` / `byName` helpers to the dapi layer. Pagination should allow fetching the full list without manual page-walking.
- Consider adding a `[name="..."]` attribute to Browse tree connection nodes that encodes the connection id (or nqName) — today the only way to distinguish three nodes called `NorthwindTest` is by a sibling `icon-ds-*` class.
- On CSV-loaded tables, the Context Panel could still show an empty "Enrich" stub with guidance ("Enrichments require a DB-backed table") instead of silently omitting the section — discovery of the feature is currently zero for anyone who hasn't seen the docs.

### Suggestions for the scenario
- Step 1.1 hardcodes "Northwind" but the dev/sandbox connections are named `NorthwindTest`. Either rename the connection, or have the scenario say "Databases > Postgres > Northwind (or NorthwindTest on dev)".
- Section 1.2 bundles three actions ("Create a SQL and a visual query for the 'Orders' table") into one bullet. Split SQL vs visual query so failures are attributable.
- Explicitly list the enrichment prerequisites at the top: working Postgres connection, `orders` and `customers` tables, `customerid` as the FK. Today the scenario assumes infrastructure that is absent on dev.
- Section 4 ("Log in with another user") is impossible to script from a single browser context without a second account and credentials. Either mark as manual-only, or provide a named secondary credential and have the scenario log out and back in explicitly.
- Section 3's "Close and reopen the project" is a stop-the-world step that needs an explicit project name and save path; otherwise the replay cannot reopen what it created.
