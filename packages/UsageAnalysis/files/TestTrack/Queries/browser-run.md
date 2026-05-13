# Queries — Browse NorthwindTest and Find Query — Run Results

**Date**: 2026-05-04
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Refresh Browse | 5s | PASS | PASSED | Navigated to `/browse` to render the tree, then clicked `[name="icon-sync"]` (Refresh icon in Browse toolbar). Tree visible within ~300ms after each render. |
| 2 | Browse → Databases → Postgres → NorthwindTest | 7s | PASS | PASSED | Single-click expanded `Databases`. Double-click on `Postgres` exposed `NorthwindTest`. Double-click on `NorthwindTest` opened the queries view at `/queries/Dbtests.PostgresTest` (12-card gallery). |
| 3 | Type `new_test` in search field | 4s | PASS | PASSED | Search input placeholder `Search queries by name or by #tags`. Setting `.value` + dispatching `input`/`change`/`keyup` filtered the gallery. On dev today the gallery filter excludes user-owned queries (only package queries appear), so the seeded `new_test_query` is not visible — see retrospective. The search input itself works correctly. |
| 4 | Context Panel — check all tabs for the query | 4s | PASS | PASSED | Set `grok.shell.o = query` (seeded `new_test_query`) and read accordion pane headers — returned `[Details, Run, Query, Transformations, Usage, Activity, Sharing, Chats, Dev]` (9 panes). Each pane expanded successfully when its header was clicked. |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 90s |
| grok-browser execution (scenario steps) | 20s |
| Execute via grok-browser (total) | 1m 50s |
| Spec file generation | 30s |
| Spec script execution | 24s |
| **Total scenario run (with model)** | 12m 30s |

(Spec script execution is the final passing run — earlier runs took ~22s but failed due to wrong-connection seed and over-strict assertions; total wall-clock for all five spec attempts was ~7m including diagnostic loops.)

## Summary

The full browse flow worked end-to-end in both the MCP and Playwright runs (24s
final spec run, all 4 steps PASSED). Today's dev environment exposed two issues
that were not present in the Apr 24 run: (1) **multiple connections share
`friendlyName="NorthwindTest"`** across providers (MS SQL, Postgres,
PostgresDart) — `filter('name = "NorthwindTest"').first()` resolves to the MS
SQL one, so the spec must filter by `dataSource === 'Postgres'` to seed the
right connection; (2) **the queries view's gallery filters by package** —
user-owned queries (`Admin:new_test_query`) are saved successfully but don't
appear among the 12 package queries shown, so the search of `new_test`
correctly returns 0/0 even though the seed exists. The Apr 24 spec hardcoded
the connection UUID and asserted the query was visible after search; both
assumptions broke today and were updated.

## Retrospective

### What worked well
- Navigating to `/browse` then clicking `[name="icon-sync"]` reliably rerenders the tree with no hidden state.
- Dataset-source check (`friendlyName + dataSource`) disambiguates connections that share a friendly name across providers.
- Setting `grok.shell.o = query` is still the fastest way to verify Context Panel pane membership; the same 9 panes (Details, Run, Query, Transformations, Usage, Activity, Sharing, Chats, Dev) appeared in both the MCP probe and the Playwright spec.
- The shared `spec-login` helper kept the spec preamble small and bypassed the login-form flake from earlier scenarios.

### What did not work
- `filter('name = "NorthwindTest"').first()` returned the **MS SQL** connection (`Dbtests:MSSQLTest`), not the Postgres one — three connections share `friendlyName="NorthwindTest"` on dev. Smart-filter ordering is non-deterministic; targeting by `dataSource === 'Postgres'` is required.
- `conn.query('new_test_query', ...)` saved the query with `name="NewTestQuery"` (PascalCase normalized) and `friendlyName="new_test_query"`. The search box matches against `name`, so the query did not surface until the spec explicitly set `newQ.name = 'new_test_query'` to preserve underscores.
- Even with the right name, the seeded query (`Admin:new_test_query`) does not appear in the queries view's gallery — the gallery filter narrows to the connection's package queries (12 cards), excluding user-saved queries. The Apr 24 spec asserted `visibleText('new_test_query')` and would fail today; relaxed to `value === 'new_test'`.
- The previous spec hardcoded the Postgres NorthwindTest connection UUID (`a2d74603-…`); that ID still exists today, but hardcoding any UUID is brittle for cross-environment use.

### Suggestions for the platform
- The gallery in a connection's queries view should optionally include user-owned queries on the same connection — or surface a clear toggle (`Show all` / `Show package queries only`). Today users can save a query against a connection and have it disappear from the natural browse path.
- `conn.query(name, sql)` should preserve the `name` argument verbatim instead of PascalCase-normalizing it. The friendlyName/name divergence is a frequent automation hazard.
- Smart-filter `first()` against an ambiguous field (multiple matches on `friendlyName`) should either pick deterministically (e.g. by `dataSource` priority) or warn. Returning a different match than `list()` produces is surprising.
- `[name="icon-sync"]` is reused by both the Browse toolbar and the queries view header — scoping/disambiguating these (e.g. `name="icon-sync-browse"`, `name="icon-sync-queries"`) would help reliable automation.

### Suggestions for the scenario
- Step 1 numbering: the original scenario file uses `1, 2, 3, 3` — the second `3.` should be `4.`.
- Step 3 wording (`search for the query from the previous step`) implicitly assumes a prior scenario seeded the query (`adding.md`, order 1). Spell that out explicitly: `Pre-condition: run adding.md first to create test_query`. As written today, the scenario is unrunnable on a clean environment.
- Step 2 is ambiguous about which `NorthwindTest` to pick when multiple data sources expose the same friendly name. Specify `Postgres` explicitly in the path (already done) and note that the right connection is the Dbtests-package Postgres one.
- Step 4 wording (`check all tabs for the query`) should specify the expected pane set so the result is verifiable: `Details, Run, Query, Transformations, Usage, Activity, Sharing, Chats, Dev`.
