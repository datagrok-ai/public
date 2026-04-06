# Catalogs — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Databases | PASS | PASSED | Databases view accessible |
| 2 | Expand MS SQL > NorthwindTest | FAIL | FAILED | MS SQL provider expanded shows only `Northwind` connection — `NorthwindTest` connection does not exist on public.datagrok.ai; all subsequent steps are blocked |
| 3 | Verify Catalogs node appears | SKIP | SKIP | Skipped — NorthwindTest not found |
| 4 | Expand Catalogs — list of catalogs appears | SKIP | SKIP | Skipped |
| 5 | Expand catalog — schemas load | SKIP | SKIP | Skipped |
| 6 | Expand schema — tables load | SKIP | SKIP | Skipped |
| 7 | Click catalog node — Context Panel shows preview | SKIP | SKIP | Skipped |
| 8 | Check Context Panel shows catalog name | SKIP | SKIP | Skipped |
| 9 | Verify catalog name displayed correctly | SKIP | SKIP | Skipped |
| 10 | Select catalog — check Comment/LLM comment fields | SKIP | SKIP | Skipped |
| 11-17 | Catalog meta, Browse, Open as table | SKIP | SKIP | Skipped |

## Summary

1 step passed, 1 failed, 15 skipped. The required `NorthwindTest` MS SQL connection is not present on public.datagrok.ai. Only `Northwind` (Postgres-based) is available. This scenario requires a MS SQL server with catalog support and a pre-configured NorthwindTest connection.

## Retrospective

### What worked well
- Databases tree navigation works correctly; MS SQL provider is accessible

### What did not work
- `NorthwindTest` MS SQL connection not present on public.datagrok.ai

### Suggestions for the platform
- N/A (feature not testable without proper connection)

### Suggestions for the scenario
- Add a precondition: "Requires MS SQL NorthwindTest connection on the target server"
- Specify which server/environment has the NorthwindTest MS SQL connection (dev/release only?)
- Note that public.datagrok.ai may not have this connection — use dev/release environment instead
