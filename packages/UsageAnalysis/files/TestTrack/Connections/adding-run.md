# Adding — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Browse > Databases | PASS | PASSED | Navigated to `/connect?browse=connections`; Databases view shows all DB providers: AWS, Access, ClickHouse, MS SQL, MariaDB, MySQL, Neo4j, Oracle, Postgres, etc. |
| 2 | Right-click Postgres and select Add connection | PASS | PASSED | Right-clicked Postgres tree node; context menu appeared with "New connection...", "Browse queries", "Browse connections"; clicked "New connection..." |
| 3 | Enter `test_postgres` to the Name field | PASS | PASSED | "Add new connection" dialog opened with General/Credentials/Cache/Properties/Indexing tabs; Name field filled with `test_postgres` |
| 4 | Fill Server, Port, DB, login, password | PASS | PASSED | Server=db.datagrok.ai, Port=54322, Db=northwind, Login=datagrok set via native setters; password placeholder used (real password requires DevOps) |
| 5 | Click the Test button | PARTIAL | PARTIAL | TEST button clicked; returned auth error: FATAL: password authentication failed for user "datagrok" — expected without real credentials; TEST button UI works correctly |
| 6 | Click OK | PASS | PASSED | OK clicked; dialog closed; `test_postgres` appeared in Postgres connections list |
| 7 | Create another connection `test_postgres_2` | PASS | PASSED | Clicked "New connection..." button; filled same fields with Name=test_postgres_2; clicked OK; both appear in list (12/12) |

## Summary

6 of 7 steps fully passed, Step 5 was partial (TEST button works but actual connection test fails without real credentials). Both `test_postgres` and `test_postgres_2` connections were created successfully and appear in the Postgres connections browser.

## Retrospective

### What worked well
- Context menu on Postgres tree node reliably shows "New connection..." action
- "Add new connection" dialog has clear field layout with General/Credentials tabs
- Native setter pattern works for all text inputs (Name, Server, Port, Db, Login)
- Connection saves even without passing the Test (OK button always enabled)
- "New connection..." button in the connection browser toolbar for adding a second connection

### What did not work
- Real credentials (password for db.datagrok.ai:54322) not available — test fails with auth error; scenario requires DevOps to supply the password

### Suggestions for the platform
- Add `data-testid` attributes to dialog input fields (Name, Server, Port, Db, Login, Password) for reliable automation
- Show a success toast when a connection is saved (currently no feedback)

### Suggestions for the scenario
- Specify the actual password or provide a test account with known credentials
- Step 4: Clarify which fields must be filled — "Fill other fields" is too vague for automation
- Add a note that the connection can be saved without passing the Test check
