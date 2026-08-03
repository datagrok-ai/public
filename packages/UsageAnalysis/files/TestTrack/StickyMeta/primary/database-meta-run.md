# Database meta — Run Results

**Date**: 2026-04-23
**URL**: https://dev.datagrok.ai
**Status**: FAIL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | DB: open Browse > Databases > Postgres > CHEMBL | 25s | PASS | FAILED | MCP: Browse tree expanded, `grok.shell.o.name = "Chembl"`. Spec click race after fresh login — `grok.shell.o` undefined at assertion time. |
| 2 | DB: Context Panel shows "Database meta" section | 10s | FAIL | FAILED | Accordion headers: Details, Queries, History, Sharing, Activity, Chats, Dev. No "Database meta". |
| 3 | DB: fill Comment=test, LLM Comment=test, Save | 5s | SKIP | FAILED | Section absent → nothing to fill. Playwright reports input not found. |
| 4 | DB: reload, verify values persist | 3s | SKIP | FAILED | Save never executed. |
| 5 | DB: clear values, save, reload, verify empty | 2s | SKIP | FAILED | Save never executed. |
| 6 | Table: open NorthwindTest > Schemas > Public | 12s | FAIL | FAILED | Browse tree under `NorthwindTest` exposes queries only; no `Schemas` child node. Tried URLs `/dbtable/...`, `/db/...`, `/connection/{id}/schemas` — all resolve to Home. |
| 7 | Table: Context Panel shows "Database meta" section | 2s | SKIP | FAILED | Cannot reach a DbTable entity. |
| 8 | Table: fill/save/reload/verify/clear (3 sub-steps) | 0s | SKIP | — | Depends on step 6/7; not attempted. |
| 9 | Column: open categories > categoryid | 8s | FAIL | FAILED | `/dbtable/Postgres.PostgresTest.public.categories?browse=db` returns 200 HTML but client routes to Home view. |
| 10 | Column: Context Panel shows "Database meta" section | 2s | SKIP | FAILED | Cannot reach a DbColumn entity. |
| 11 | Column: fill/save/reload/verify/clear + is_unique (3 sub-steps) | 0s | SKIP | — | Depends on step 9/10; not attempted. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome (plain `PASSED`/`FAILED`/`SKIPPED`).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 4m |
| grok-browser execution (scenario steps) | 3m 10s |
| Execute via grok-browser (total) | 7m 10s |
| Spec file generation | 1m 20s |
| Spec script execution | 33s |
| **Total scenario run (with model)** | 9m 3s |

## Summary

The "Database meta" context-panel section is not rendered on `dev.datagrok.ai` for a
Postgres DbInfo connection entity (Chembl was verified). The platform has all the
backing sticky-meta schemas configured (`db-db`, `db-schema`, `db-table`, `db-column`,
`db-relation` with the expected `comment` / `llmComment` / `isUnique` properties), and
`grok.dapi.stickyMeta.getSchemas()` returns them — but no `Database meta` accordion
pane is registered for DbInfo entities. The Table-level and Column-level steps could
not even be reached: the Browse tree under `Postgres > NorthwindTest` exposes queries
only (no Schemas child), and no URL route resolves to a DbSchema, DbTable, or DbColumn
entity view. All 8 Playwright softSteps failed (matches MCP outcomes). **Total scenario
run (with model): 9m 3s.**

## Retrospective

### What worked well
- `grok.dapi.stickyMeta.getSchemas()` returned the full set of metadata schemas with
  properties matching the scenario's field list — confirmed backend is seeded.
- `grok.dapi.connections.filter('dataSource = "Postgres" and name = "Chembl"')`
  located the target connection cleanly.
- Playwright spec preamble (verbatim from the template) logged in from scratch and
  reached the Browse sidebar on a fresh context.

### What did not work
- **"Database meta" accordion pane missing on DbInfo context panel** — root cause
  likely a missing platform-side `#panel` registration or the sticky-meta → accordion
  wiring isn't picking up `onTypes=DbInfo`. The schema `db-db` exists but nothing
  renders on the entity.
- **No Browse-tree path to DB schemas/tables/columns.** Under
  `Databases > Postgres > NorthwindTest` the tree lists queries only — there is no
  `Schemas` sub-node. `Databases > Postgres > NorthwindTest > Schemas > Public >
  categories > categoryid` as described in the scenario does not exist on dev.
- **No URL route for DbTable/DbColumn entities.** `/dbtable/...`, `/db/...`,
  `/connection/{id}/schemas` all fall back to the Home view. The referenced
  `CmdShowSchema` function opens an empty `Model` view.
- First Playwright softStep (`DB: open Browse > Databases > Postgres > CHEMBL`) failed
  because `grok.shell.o` was still undefined when read — the `click()` + 2.5s wait is
  not enough on a post-login cold start. Increasing that wait would make the step
  PASS but the next steps would still fail on the missing panel.

### Suggestions for the platform
- Register the `Database meta` info panel for `DbInfo` / `DbSchema` / `DbTable` /
  `DbColumn` entities whenever a matching sticky-meta schema exists (`onTypes`
  contains the entity type). Without this, the schemas configured in
  `Platform > Sticky Meta` have no user-facing surface.
- Expose a `Schemas` node under each DB connection in the Browse tree, so users can
  drill to `schemas > tables > columns` as the scenario describes. Today the tree
  only shows queries — database introspection is available server-side but
  unreachable from the UI.
- Add a stable client-side route for DB entities (e.g. `/db/{conn}/{schema}/{table}`
  and `/db/{conn}/{schema}/{table}/{column}`) so tests and deep links can reach them
  without going through the tree.
- Known-issue cross-reference: GROK-19427 (column meta save: `Value test doesn't
  match expected type string_list`) and GROK-19429 (row count deleting not working)
  should be linked off the `db-column` schema's definition so the save code path is
  auditable.

### Suggestions for the scenario
- Add a precondition line: *"Sticky-meta schemas `db-db`, `db-schema`, `db-table`,
  `db-column` must be present and their corresponding info-panels registered."*
  Without this, every step fails without a clean diagnosis.
- Replace the "Open Browse > Databases > Postgres > NorthwindTest > Schemas > Public
  ..." navigation with an explicit entity-URL or a `grok s` command so the scenario
  doesn't silently depend on a Browse-tree feature that isn't universally enabled.
- Step 2 of the Database-metadata section says *"Reload platform. Open ... CHEMBL
  again. Expected: Previously entered metadata is saved and displayed correctly."* —
  the Table section copy-pastes this line but opens CHEMBL instead of NorthwindTest;
  looks like a copy-paste bug. Should be "Open ... NorthwindTest > Schemas > Public
  again."
- The "LL M Comment" label in the scenario is probably a typo for "LLM Comment"
  (which matches the `llmComment` property). Fix the label in the .md.
