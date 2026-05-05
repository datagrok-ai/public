# Catalogs — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Databases | 8s | PASS | PASSED | Navigated to /browse/databases; tree showed all providers including MS SQL |
| 2 | Expand MS SQL > NorthwindTest | 18s | PASS | PASSED | MS SQL expander toggled on second try (first label-click only selected). NorthwindTest visible alongside Northwind, MSSQLDBTests; expanding NorthwindTest produced the schema-group placeholder |
| 3 | Verify Catalogs node (vs Schemas) | 12s | FAIL | FAILED | Schema-group label resolved to `Unavailable` with `fa-exclamation-circle` icon — not `Catalogs`. Root cause: MS SQL DB host refuses TCP/IP connections on dev (verified via `conn.test()`: `com.microsoft.sqlserver.jdbc.SQLServerException: ... Connection refused`). All 3 MS SQL connections (NorthwindTest, Northwind, MSSQLDBTests) fail identically |
| 4 | Expand Catalogs — list of catalogs | n/a | SKIP | SKIPPED | Blocked by step 3 — no Catalogs node exists |
| 5 | Expand a catalog — schemas load | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 6 | Expand a schema — tables load | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 7 | Click on a catalog node | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 8 | Context Panel shows catalog preview | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 9 | Catalog name displayed correctly | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 10 | Select a catalog in the tree | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 11 | Comment / LLM comment meta props assignable | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 12 | Comment persists after re-selecting | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 13 | LLM comment can be set and is displayed | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 14 | Right-click catalog > Browse | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 15 | Schema view opens with all schemas | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 16 | Right-click catalog > Open as table | n/a | SKIP | SKIPPED | Blocked by step 3 |
| 17 | Table with catalog DB objects opens | n/a | SKIP | SKIPPED | Blocked by step 3 |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 30s |
| grok-browser execution (scenario steps) | 1m 5s |
| Execute via grok-browser (total) | 2m 35s |
| Spec file generation | 2m 0s |
| Spec script execution | 45s |
| **Total scenario run (with model)** | 5m 20s |

## Summary

**Catalog browsing cannot be exercised on dev** — all three MS SQL connections in
`Browse > Databases > MS SQL` (NorthwindTest, Northwind, MSSQLDBTests) point to
`db.datagrok.ai:14331` / `:14332`, both of which currently refuse TCP/IP
connections (`SQLServerException: Connection refused`). The tree therefore
expands to a single `Unavailable` placeholder under each connection instead of
the expected `Catalogs` group, so steps 4-17 (catalog list, schema/table
drill-down, Context Panel preview, meta-property persistence, right-click
Browse / Open as table) cannot be reached. Playwright spec correctly captures
this state: steps 1-2 PASSED, step 3 FAILED on the intended `Catalogs`
assertion, remaining 14 steps SKIPPED via `test.skip()` guarded on the
schema-group label. **Total scenario run (with model): ~5m 20s**.

## Retrospective

### What worked well
- `grok.dapi.connections.filter('name = ...').first()` + `conn.test()` made it trivial to confirm the underlying DB outage and rule out a UI bug
- Tree expansion via `.d4-tree-view-tri` chevron click is reliable; `name="tree-Databases---MS-SQL---NorthwindTest"` gave a stable selector
- The `Unavailable` schema-group label (with `fa-exclamation-circle` icon) is correct, informative platform behavior — exactly what the connections.md reference documents for "Schema load failed"
- `test.skip(condition, reason)` cleanly elides downstream steps when the prerequisite is impossible, and the reason string surfaces in the Playwright report

### What did not work
- The MS SQL test database (`db.datagrok.ai:14331`/`14332`) is offline on dev — root cause is infrastructure, not the platform
- First MS SQL expand attempt (clicking the label) only selected the node — had to fall back to clicking `.d4-tree-view-tri` explicitly. Documented in tree_view.md but easy to miss
- The schema-group "loading" placeholder reuses the parent connection's `name` attribute (`name="tree-Databases---MS-SQL---NorthwindTest"`), which makes it ambiguous to target by `name` alone — had to walk the DOM via `.children[0]`
- `Unavailable` node has no `title` attribute or tooltip exposing the underlying error — user has to right-click → `Test connection` to see why

### Suggestions for the platform
- **Surface the schema-load error on the `Unavailable` node** as a `title=` attribute (or hover tooltip), so users don't have to invoke `Test connection` separately to diagnose
- **Auto-retry schema enumeration** with exponential backoff on transient connection failures, or add a `Retry` action to the right-click menu of the `Unavailable` node
- **Distinguish `Unavailable` from `Loading`** more strongly in the icon — both currently render a placeholder until the load resolves
- Give the schema-group placeholder a distinct `name=` (e.g. `tree-...---{Schemas|Catalogs|Unavailable}`) so automation can target it without DOM-walking
- Add a Browse-tree health badge on the connection node (red dot) when its last schema load failed, so the failure is visible without expanding

### Suggestions for the scenario
- **Pre-condition is missing** — the scenario assumes a working MS SQL connection with `browseCatalogs` enabled. Add an explicit pre-condition row: "An MS SQL (or other multi-catalog) connection that loads its schema successfully exists in `Browse > Databases > MS SQL`. Verify with `Test connection` first."
- Spell out which Datagrok build / fixture sets up the catalog-enabled connection — without this, the scenario silently no-ops on environments where the test DB is down
- Add a fallback / alternative provider note (e.g. "On Snowflake/BigQuery with `browseCatalogs=true`, the same UI applies"). Right now `SnowflakeDBTests` on dev shows `Schemas` not `Catalogs` because the connection doesn't enable `browseCatalogs`
- Step 3 wording "double-database icon" — the actual icon class is `svg-database-tables` for `Schemas`; document the expected icon class for `Catalogs` so verification is unambiguous
- Step 11 "Comment and LLM comment meta properties can be assigned to the catalog" — ambiguous whether these are global or per-connection meta-properties; clarify where they appear in the Context Panel (Details? Properties? a new pane?)
- Step 16 "Open as table" — clarify whether this is the same as the connection-level `Open schema as table` (documented in connections.md) but scoped to the selected catalog
