# External Provider — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Browse > Databases > Postgres > Add connection | SKIP | SKIP | Connection creation UI works (verified in Adding.md), but this scenario requires specific credentials |
| 2 | Fill Name=PostgreSQLDBTests2, Server=db.datagrok.ai, Port=54327, Db=test, Login=superuser | SKIP | SKIP | Fields can be filled (mechanism verified in Adding.md), but Password requires DevOps/QA |
| 3 | Add and run TestCreateTable query | SKIP | SKIP | Requires working connection — skipped due to missing credentials |
| 4 | Add and run TestInsertData query | SKIP | SKIP | Skipped |
| 5 | Add and run TestUpdateData query | SKIP | SKIP | Skipped |
| 6 | Add and run TestDropTable query | SKIP | SKIP | Skipped |
| 7 | Delete PostgreSQLDBTests2 connection | SKIP | SKIP | Skipped |

## Summary

All 7 steps skipped. This scenario requires a specific Postgres connection at db.datagrok.ai:54327 with superuser credentials that are not available without DevOps/QA access. The underlying UI mechanisms (create connection, run queries, delete connection) were all verified in other scenarios.

## Retrospective

### What worked well
- N/A

### What did not work
- Credentials for db.datagrok.ai:54327 (superuser) not available

### Suggestions for the platform
- Provide a publicly accessible test Postgres instance with known credentials for integration testing

### Suggestions for the scenario
- Consolidate credential storage — move all test passwords to a shared, versioned secret store
- Consider running this scenario only in the DevOps/QA environment where credentials are available
