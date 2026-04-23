# Delete — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Find connection `new_test_postgres` in Browse → Connections | PASS | PASSED | Found `new_test_postgres` in Postgres connections browser; also accessible via search filter "new_test" → 1/1 result |
| 2a | Right-click and select Delete from context menu | PASS | PASSED | Right-clicked → "Delete..." appeared in menu; clicked it |
| 2b | In confirmation dialog, click YES | PASS | PASSED | Confirmation dialog showed DELETE and CANCEL buttons (not YES); clicked DELETE |
| 2c | Check connection disappeared | PASS | PASSED | After refresh, `new_test_postgres` no longer appears in Postgres connections list (11/11) |
| 3 | Find connection `test_postgres_2` in Browse → Databases | PASS | PASSED | Found `test_postgres_2` in Postgres connections list |
| 4a | Right-click and select Delete from context menu | PASS | PASSED | Right-clicked → Delete... clicked |
| 4b | In confirmation dialog, click YES | PASS | PASSED | Clicked DELETE in confirmation dialog |
| 4c | Check connection is no longer present | PASS | PASSED | After refresh, `test_postgres_2` no longer appears in Postgres connections list (10/11 - one test_postgres from a previous session remains) |

## Summary

All 8 steps passed. Both connections were deleted successfully. The confirmation dialog uses a red "DELETE" button (not "YES" as described in the scenario). Refreshing the view confirms the connections are removed.

## Retrospective

### What worked well
- "Delete..." context menu item works on both tree nodes and list items
- Confirmation dialog reliably appears with DELETE/CANCEL choice
- Refresh icon reloads the connection list immediately after deletion

### What did not work
- Nothing

### Suggestions for the platform
- The confirmation dialog button is labeled "DELETE" (not "YES" as the scenario says) — should be consistent with documentation
- Add a brief success toast after deletion ("Connection deleted")

### Suggestions for the scenario
- Steps 2b and 4b say "click YES" but the actual button is "DELETE" — update scenario wording
- Specify that "Browse → Platform → Connections" and "Browse → Databases" both show the same Postgres connections
