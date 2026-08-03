# Schema — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Go to Browse > Databases | 5s | PASS | PASSED | Navigated to `/connect?browse=connections`; tree populated with all DB providers (104 group labels including Postgres) |
| 2 | Expand the Postgres provider | 8s | PASS | PASSED | Clicked tri on Postgres node; 16 named child connections visible (CHEMBL, ChemblTest, Datagrok, Northwind, NorthwindTest, PostgreSQLDBTests, PostgreSQLDBTests2, PostgreSQLDBTestsCached, PostgresDocker, Projects, Starbucks×2, SureCHEMBLDocker, World, new_test_postgres, test_postgres) |
| 3 | Right-click Northwind and select Browse | 6s | PASS | PASSED | Context menu opened (Browse, New Query…, New Visual Query…, Delete…, Edit…, Rename…, Clone…, Clear cache, Browse queries, Test connection, Share…, ID, Grok name, Markup, URL, Copy, Add to favorites). Clicked `[name="div-Browse"]` → URL `/queries/Samples.PostgresNorthwind?browse=db`, view name `PostgresNorthwind`, view type `queries`. The connection's saved-queries grid is shown in the center, the Details/Queries/History/Sharing/Activity/Chats accordion on the right, and the side Browse tree expands Northwind in place to expose Schemas/queries |
| 4 | Check you can interact with structures as with DB Tables | 22s | PASS | PASSED | After Browse, side tree shows Northwind → `Schemas` (3: public, pg_catalog, information_schema) plus 9 saved queries as siblings. Drilled `Schemas` → `public` (14 tables). Right-click on `customers` → menu shows: Get All, Get Top 100, New SQL Query…, New Visual Query… — canonical DB-table operations |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 35s |
| grok-browser execution (scenario steps) | 22s |
| Execute via grok-browser (total) | 57s |
| Spec file generation | 25s |
| Spec script execution | 36s |
| **Total scenario run (with model)** | 2m 35s |

## Summary

All 4 steps pass. The previous run mis-read step 3 as "Browse schema" (a non-existent menu item) and marked it AMBIGUOUS; step 3 actually says "Browse", which exists, opens the connection's queries view, and side-effect-expands Northwind in the Browse tree so step 4 can drill into `Schemas → public → customers` and verify the DB-table context menu. The Playwright spec re-runs the corrected flow end-to-end and passes in 36 seconds. **Total scenario run (with model): 2m 35s**.

## Retrospective

### What worked well
- Tree drill-down via `:scope > .d4-tree-view-group-host` is reliable once a node is visible
- Clicking `[name="div-Browse"]` from the connection context menu both navigates to the queries view AND auto-expands Northwind in the side tree, so step 4 doesn't need to expand the connection again
- The customers-table context menu exposes Get All / Get Top 100 / New SQL Query… / New Visual Query… exactly as documented in `connections.md` (Schema Group / Table Nodes)

### What did not work
- **Initial misreading of step 3.** "select **Browse**" was misread as "select **Browse schema**" — `Browse schema` is documented in `connections.md` but does not appear in the Northwind context menu on dev. The scenario step actually asks for `Browse`, which is present. Both the prior `schema-run.md` and the first attempt this run made the same mistake; corrected after re-reading the source `.md`
- First Playwright run failed because of a fixed 2s wait between page load and the Postgres-node lookup — on a cold load the tree wasn't populated yet. Replaced with a poll loop (up to 30s) for the `Postgres` label

### Suggestions for the platform
- The `Browse` menu item opens the **queries view** for the connection, not a dedicated "schema view". This is fine for the scenario as written (because the side Browse tree exposes the schema inline), but the platform-wide gap remains: `connections.md` still documents `[name="div-Browse-schema"]` ("Opens a separate schema view"), and that item is not present on `Postgres → Northwind` on dev. Either reinstate the menu item or remove the stale entry from the reference

### Suggestions for the scenario
- Step 3 wording is fine, but step 4 leans on a vague "schema view" — clarify that the structures to interact with live in the side Browse tree (Schemas → {schema} → {table}), not in the central queries grid that opens
- Add a precondition naming `Northwind` (not `NorthwindTest` or any other Postgres connection) as the expected target
- Note that Northwind exposes 9 saved queries as siblings of `Schemas` in the tree; the wording "structures presented in schema view" should explicitly scope step 4 to schema/table nodes, not the saved queries
