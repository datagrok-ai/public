# Database meta — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: FAIL

## Steps

**Database metadata display (CHEMBL)**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Browse > Databases > Postgres > CHEMBL | PASS | PASSED | Navigated to Connections view; selected CHEMBL (Postgres); Context Panel shows Details, Queries, History, Sharing, Activity, Chats, Dev |
| 2 | Observe Context Panel — expect "Database meta" section | FAIL | FAILED | No "Database meta" section present in Context Panel for CHEMBL connection. Known issues GROK-19427 and GROK-19429 may be related. Sections shown: Details, Queries27, History42, Sharing, Activity5, Chats, Dev |
| 3 | Fill Comment and LLM Comment, click Save | SKIP | SKIP | Section not found |
| 4 | Reload platform; reopen CHEMBL | SKIP | SKIP | No data to verify |
| 5 | Clear metadata, save, reload | SKIP | SKIP | No data to verify |

**Table metadata display (NorthwindTest)**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Browse > Databases > Postgres > NorthwindTest > Schemas > Public | SKIP | SKIP | NorthwindTest connection not present on public.datagrok.ai (confirmed in Connections test run) |
| 2–5 | Fill/verify/clear table metadata | SKIP | SKIP | Prerequisite missing |

**Column metadata display (NorthwindTest > categories > categoryid)**

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1–4 | Column metadata for NorthwindTest column | SKIP | SKIP | NorthwindTest not available on public.datagrok.ai |

## Summary

The "Database meta" section was not found in the Context Panel for the CHEMBL database connection on public.datagrok.ai. The CHEMBL connection is accessible and selectable, but no Database meta accordion section appears. The NorthwindTest connection is not present on public.datagrok.ai, making all table/column metadata steps untestable. The scenario is fully blocked.

## Retrospective

### What worked well
- CHEMBL connection is accessible and its context panel loads correctly
- Other connection metadata (Details, Sharing, Activity) works correctly

### What did not work
- "Database meta" section absent from CHEMBL connection context panel — may be a plugin/package not loaded on public.datagrok.ai, or related to known issues GROK-19427/GROK-19429
- NorthwindTest not available on public.datagrok.ai — all table/column metadata steps blocked

### Suggestions for the platform
- The Database meta section should be visible on public.datagrok.ai for at least one accessible connection to make this scenario testable
- Add a visible indicator when the db-db / db-schema metadata schemas are not configured for a connection

### Suggestions for the scenario
- Note GROK-19427 and GROK-19429 known issues at the top as prerequisites to check before running
- Replace NorthwindTest steps with a connection available on public.datagrok.ai (e.g., CHEMBL or Northwind)
- Add a prerequisite note that the db-schema / db-table packages must be installed and configured
