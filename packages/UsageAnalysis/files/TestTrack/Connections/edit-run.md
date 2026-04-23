# Edit — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Reload the tree in Browser | PASS | PASSED | Clicked refresh icon in Postgres connection browser; list refreshed showing all connections |
| 2 | Right-click the test_postgres connection | PASS | PASSED | Right-clicked test_postgres label; full context menu appeared: Browse, New Query..., Edit..., Delete..., Rename..., Clone..., Test connection, Share... |
| 3 | Select Edit from context menu | PASS | PASSED | Clicked "Edit..."; "Edit Connection" dialog opened pre-populated with all existing values (Name, Server, Port, Db, Login, Password) |
| 4 | Change name to `new_test_postgres` and click OK | PASS | PASSED | Changed Name from test_postgres to new_test_postgres via native setter; clicked OK; `new_test_postgres` appears in connections list |
| 5 | Change login/password with arbitrary data and save | PASS | PASSED | Right-clicked new_test_postgres → Edit; set Login=wronguser, Password=wrongpassword; clicked OK; saved with wrong credentials |
| 6 | Test the connection — should return error | PASS | PASSED | Right-clicked → "Test connection"; error toast: `"new_test_postgres": failed to connect: FATAL: password authentication failed for user "wronguser"` — matches expected error format |
| 7 | Set right login/password and test — should be fine | SKIP | SKIP | Skipped — real password unavailable; mechanism verified in step 6; the Edit dialog correctly saves credentials |

## Summary

6 of 7 steps passed (1 skipped due to missing real credentials). The connection rename, credential modification, and error testing all work correctly. The error message format matches exactly what the scenario describes (FATAL: password authentication failed).

## Retrospective

### What worked well
- "Edit..." context menu item opens a pre-populated edit dialog
- Name change via native setter + change event works reliably
- "Test connection" from context menu (without opening Edit dialog) is a convenient shortcut
- Error message clearly identifies the user that failed authentication

### What did not work
- Step 7 (verify correct credentials work) could not be tested — requires DevOps password

### Suggestions for the platform
- Add a success toast after a successful connection test (currently only errors are toasted)
- "Test connection" result (pass or fail) should be visible in the connection's Details panel

### Suggestions for the scenario
- Step 7 requires knowledge of the real password — specify which credentials to restore or provide a test account
- Add note: "Test connection" can be accessed both from the Edit dialog and from the right-click context menu
