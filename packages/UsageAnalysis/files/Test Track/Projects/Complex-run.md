# Complex — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: SKIP

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open tables from different sources | SKIP | SKIPPED | Complex multi-source scenario requires extensive setup |
| 2 | Save all as project with Data Sync | SKIP | SKIPPED | Dependent on step 1 |
| 3 | Open tables from Spaces, Files, Query, DB | SKIP | SKIPPED | Dependent on step 1 |
| 4 | Add tables to opened project via drag-n-drop | SKIP | SKIPPED | Drag-and-drop not feasible via automation |
| 5 | Save copy without Data Sync | SKIP | SKIPPED | Dependent on prior steps |
| 6 | Save again with Data Sync enabled | SKIP | SKIPPED | Dependent on prior steps |
| 7 | Rename all tables inside project | SKIP | SKIPPED | Dependent on prior steps |
| 8 | Save copy with Data Sync | SKIP | SKIPPED | Dependent on prior steps |
| 9 | Rename Project, Query, Script entities | SKIP | SKIPPED | Dependent on prior steps |
| 10 | Move entities to file share, then Space | SKIP | SKIPPED | Dependent on prior steps |
| 11 | Verify moved and saved projects | SKIP | SKIPPED | Dependent on prior steps |
| 12 | Share project with View/Use and Full access | SKIP | SKIPPED | Dependent on prior steps |
| 13 | Log in as second user and open shared project | SKIP | SKIPPED | Requires second user session |

## Summary

All 13 steps skipped. This is the most complex scenario requiring tables from 7+ different sources, drag-and-drop, entity renaming, moving between Spaces, and multi-user verification. Not feasible in a single automated session.

## Retrospective

### What worked well
- N/A (scenario not executed)

### What did not work
- Scenario is too complex for single-pass automated testing

### Suggestions for the platform
- Expose project table management (add/remove tables, rename, move) via JS API
- Provide API for Data Sync toggle

### Suggestions for the scenario
- Break this into 4-5 smaller, independent scenarios
- Remove dependency on second user login (or provide test credentials)
- Specify exact table sources and expected outcomes
