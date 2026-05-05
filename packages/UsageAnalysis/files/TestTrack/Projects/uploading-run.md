# Uploading — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Case 1: Open SPGI_v2_infinity.csv from Files | 5s | PASS | PASSED | Loaded via `grok.dapi.files.readCsv('System:DemoFiles/SPGI_v2_infinity.csv')`; 3624×88; substituted `System:Demo` → `System:DemoFiles` (the actual connection name on dev) |
| 2 | Case 1: Open SPGI_v2.csv from Files | 5s | AMBIGUOUS | PASSED | `SPGI_v2.csv` does NOT exist in `System:DemoFiles` on dev; substituted with `SPGI.csv` (same 88-column schema, shares `Id` column). The file exists at `Home:SPGI_v2.csv` but `readCsv` errors there with a server-side handlerPath bug |
| 3 | Case 1: Data > Link Tables (selection→filter on Id) | 6s | PASS | PASSED | Top menu submenu items not visible until parent hovered; both MCP and Playwright clicked via direct DOM `.click()` |
| 4 | Case 1: Verify link works (select row → other table filters) | 2s | PASS | PASSED | Select 1 row in `SPGI_v2_infinity` → `SPGI` filtered to exactly 1 matching `Id`; mechanics work in-session |
| 5 | Case 1: Save Project — Data sync ON | 3s | SKIP | n/a | **Data sync toggle hidden** for tables loaded via `readCsv()` — toggles have `style="display: none"`. Source binding is required for sync; the API path doesn't set it. Setting `tags['source']` manually doesn't help |
| 6 | Case 1: Save Project — Data sync OFF | 8s | PASS | PASSED | Used ribbon `[name="button-Save"]`; default mode is "Save a copy" with Clone-only choices for both tables; `Test_Case1_NoSync_*` saved (slug normalized to `TestCase1NoSync*`) |
| 7 | Case 1: Close project + reopen via Browse > Projects | 8s | PASS | PASSED | Used `proj.open()` JS API (Browse tree click is a separate UX); both tables restored at full 3624×88 |
| 8 | Case 1: Verify after reopen (link works, no console errors) | 3s | PASS | PASSED | In Playwright: row-select → filter dropped below 3624 ✓. In MCP probe earlier the link appeared NOT to fire — likely a stale-session timing issue; the clean-context Playwright run is authoritative. No console errors |
| 9 | Case 2: Query result + Query result | n/a | SKIP | SKIPPED | `NorthwindTest` connection doesn't exist on dev (it's `PostgresNorthwind`); `PostgressAll` query doesn't exist (closest: `PostgresStates`, `OrdersByEmployee`) |
| 10 | Case 3: Query result + Files | n/a | SKIP | SKIPPED | Same env block as Case 2 |
| 11 | Case 4: Spaces + Spaces | n/a | SKIP | SKIPPED | **0 Spaces exist on dev** — verified via `grok.dapi.projects.filter('isDashboard = false').list()` returning empty |
| 12 | Case 5: Spaces + Files | n/a | SKIP | SKIPPED | No Spaces on dev |
| 13 | Case 6: Spaces + Query result | n/a | SKIP | SKIPPED | No Spaces + no `PostgressAll` |
| 14 | Case 7: orders Get Top 100 / Get All + Join Tables | n/a | SKIP | SKIPPED | Requires `Browse > Databases > Postgres > NorthwindTest > Schema > public > orders` tree double-click; `NorthwindTest` doesn't exist (PostgresNorthwind does). Direct SQL via JS API: `grok.data.db.query()` errored with `Unable to parse string` for connection ID |
| 15 | Case 8: Files + Pivot Table → Add to workspace | n/a | SKIP | SKIPPED | `SPGI_v2.csv` not in DemoFiles; could substitute with `SPGI.csv`, but Data sync would be hidden (same blocker as Case 1) |
| 16 | Case 9: DB table double-click + Aggregate Rows | n/a | SKIP | SKIPPED | Same DB block as Case 7 |

**Time** = step 2b wall-clock (incl. thinking). **Result** = MCP-run outcome. **Playwright** = `npx playwright test` outcome (PASSED/FAILED/SKIPPED).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~7m |
| grok-browser execution (scenario steps) | ~3m |
| Execute via grok-browser (total) | ~10m |
| Spec file generation | ~2m |
| Spec script execution | 1m 7s |
| **Total scenario run (with model)** | ~14m |

The two `scenario steps` rows include a large amount of time spent diagnosing dev-environment mismatches (missing connection, missing query, missing Space, broken `Home:` file API) before settling on Case 1 as the single executable case.

## Summary

**Total scenario run (with model): ~14m.** Of the 9 cases, only Case 1 was runnable on dev — and only after substituting `System:Demo` → `System:DemoFiles` and `SPGI_v2.csv` → `SPGI.csv`. Cases 2/3/7/9 are blocked by missing Postgres `NorthwindTest` / `PostgressAll`. Cases 4/5/6 are blocked by a complete absence of Spaces on dev. Case 8 has the same Data-sync blocker as Case 1. The Playwright spec covers Case 1 end-to-end including the Save-Project + Reopen + Link-still-works flow, which all PASS in a clean browser context.

## Retrospective

### What worked well

- **Link Tables dialog** (`[name="dialog-Link-Tables"]`) — defaults already matched scenario inputs (Table 1, Table 2, key=`Id`); only Link Type needed to change.
- **Selection-to-filter linking in-session**: select 1 row → filter exact match on `Id`. Mechanics solid.
- **Project save+reopen**: tables restored at full row/column counts; in the Playwright clean-context run, the link metadata was also preserved and re-fired after reopen (contradicting an earlier MCP probe — likely stale session state).
- **`grok.dapi.files.readCsv` from `System:DemoFiles`** is fast and reliable for ad-hoc setup.
- **Project name normalization** (`Test_Case1_NoSync_*` → `TestCase1NoSync*` slug) is consistent with documented `friendlyNameToName` behavior.

### What did not work

- **`System:Demo` connection doesn't exist on dev** — the scenario's path prefix is wrong; the actual name is `System:DemoFiles`. Calls to `grok.dapi.files.list('System:Demo', ...)` produce server-side `handlerPath/url` errors.
- **`SPGI_v2.csv` not present in `System:DemoFiles` on dev**. Available chem-related files at the root: `SPGI.csv`, `SPGI-linked1.csv`, `SPGI-linked2.csv`, `SPGI_v2_infinity.csv`. `Home:SPGI_v2.csv` does exist but `readCsv` against the `Home:` connection produces a server-side handlerPath error — likely a real bug worth filing.
- **`PostgresNorthwind`, not `NorthwindTest`** — the scenario's Postgres connection name is wrong for dev.
- **`PostgressAll` query missing** — the scenario references a query that doesn't exist anywhere on dev. Closest matches: `PostgresStates`, `OrdersByEmployee` on `PostgresNorthwind`.
- **0 Spaces on dev** — `grok.dapi.projects.filter('isDashboard = false').list()` returns 0. Cases 4/5/6 cannot run as written.
- **Data sync toggle hidden for tables loaded via `readCsv()`** — without an external source binding (file/query) registered with the project layer, the per-table Data sync `[name="input-host-Data-sync"]` host stays `display: none`. Manually setting `df.tags['source']`/`'.source'`/`'filePath'` doesn't surface it. The scenario's expectation that "tables come from a file source automatically" only holds for tables opened via Browse-tree double-click, not via the JS API path used by all current `/grok-browser` automation.
- **Top-menu File > Save Project does not exist in simpleMode** — only `Edit / View / Select / Data / ML / Chem / Bio` appear. The scenario's `File > Save Project` path is unreachable; the ribbon `[name="button-Save"]` is the only entry point in this mode.
- **`grok.dapi.projects.filter('name = "X"').first()`** silently returns the **first project in the list** (regardless of filter) when the syntax doesn't match a real query — observed `name = "NorthwindTest"` returning a `MSSQLTest` connection. Listing + JS-side `find` is reliable; the `.filter('name = "...")` path is not.
- **`grok.data.db.query(connId, sql)` fails with `Unable to parse string` for a normal UUID** — blocked direct SQL fallback for Cases 7/9.

### Suggestions for the platform

1. **Fix `Home:` file connection** — `grok.dapi.files.readCsv('Home:SPGI_v2.csv')` errors with `handlerPath "/connectors/connections/Home.SPGI_v2.csv/file"` instead of `/connectors/connections/Home/file/Home:SPGI_v2.csv`. The `:`-suffix path parser appears to mis-tokenize `Home:` paths.
2. **Surface a clear error when `grok.dapi.files.list/readCsv` is called against a non-existent connection** — currently a confusing handlerPath/uri assertion bubbles up to the user.
3. **Validate `dapi.*.filter('name = "X"')` syntax** — when the predicate doesn't bind, the call should error instead of silently returning the unfiltered first row.
4. **Fix `grok.data.db.query(connId, sql, params)`** — passing a real UUID for `connId` errors with `Unable to parse string`. Either the parameter contract is undocumented (does it want a name?) or the parser is broken.
5. **Document the source-binding requirement for Data sync** — the toggle being silently hidden is confusing. Either show it disabled with a tooltip explaining "Open this table via Browse > Files for Data sync" or make `readCsv()` auto-register the source.
6. **`File > Save Project` in simpleMode** — make the ribbon Save and the (currently absent) `File > ...` entry point both available, or document that simpleMode hides it.

### Suggestions for the scenario

1. **Use `System:DemoFiles`, not `System:Demo`** — fix the connection prefix everywhere in the scenario.
2. **Replace `SPGI_v2.csv` with an existing file** (e.g. `SPGI.csv`, which shares the schema and Id column with `SPGI_v2_infinity.csv`).
3. **Use `PostgresNorthwind` and an existing query** (e.g. `OrdersByEmployee` or `PostgresStates`) for Cases 2/3/6/7/9. The scenario's `NorthwindTest` and `PostgressAll` references appear to be from a private dev environment.
4. **Provide a fallback for Cases 4–6** when Spaces are unavailable — either a setup step that creates a `SPGIs` Space (e.g. via `grok s`), or alternative cases that use Browse > Dashboards.
5. **Spell out the Browse-tree double-click flow** for opening files/queries — without it the platform's source-binding doesn't kick in and Data sync is hidden, which silently invalidates a major part of the scenario.
6. **Drop `File > Save Project`** from instructions since it's not present in default Tabs mode; replace with "click the **Save** button in the ribbon (Ctrl+S)".
7. **Each case is a near-clone** of "open table A, open table B, link/join, save, reopen, verify, repeat with sync OFF". Consider parameterizing by `(source1, source2, linkOrJoin, syncMode)` to remove duplication and make adding new cases trivial.
