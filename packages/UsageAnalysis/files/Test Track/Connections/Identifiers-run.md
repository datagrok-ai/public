# Identifiers — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: FAIL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Go to Data > Databases > Postgres | PASS | PASSED | Postgres connections accessible via Browse > Databases > Postgres |
| 2 | Right-click test_postgres and select Configure Identifiers... | FAIL | FAILED | `test_postgres` connection created in Adding.md has wrong credentials (no real password) — connection exists in list but cannot connect to DB; menu does show "Configure Identifiers..." option |
| 3 | Set Schema to `public` and click OK | SKIP | SKIP | Skipped — connection not functional |
| 4 | Add identifier: CUSTOMER_ID / customers / customerid | SKIP | SKIP | Skipped |
| 5 | SAVE | SKIP | SKIP | Skipped |
| 6 | Reload and verify customerid column highlighted | SKIP | SKIP | Skipped — cannot open customers table without working connection |
| 7 | Check semantic type in Context Panel | SKIP | SKIP | Skipped |
| 8 | Remove identifiers configuration | SKIP | SKIP | Skipped |
| 9 | Reload and verify no highlighting | SKIP | SKIP | Skipped |

## Summary

1 step passed, 1 failed, 7 skipped. This scenario depends on a working Postgres connection to the Northwind database. The `test_postgres` connection exists but uses incorrect credentials. The "Configure Identifiers..." menu option is present and accessible from the context menu, confirming the UI feature exists.

## Retrospective

### What worked well
- "Configure Identifiers..." context menu item is present on connection nodes
- The Postgres connections list is accessible and `test_postgres` is found

### What did not work
- Cannot test identifier configuration without a working DB connection (real credentials unavailable)

### Suggestions for the platform
- N/A

### Suggestions for the scenario
- This scenario depends on Adding.md — note explicit prerequisite: "test_postgres must connect successfully"
- Alternatively, use the existing `Northwind` connection (which works on public.datagrok.ai) instead of test_postgres
- The order value (1) conflicts with Adding.md — should be ordered after Adding (e.g., order: 10)
