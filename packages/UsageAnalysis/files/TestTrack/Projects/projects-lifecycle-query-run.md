# Projects / Lifecycle Query: Samples:PostgresCustomers source — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 56.5s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 56.5s |
| **Total scenario run** | 56.5s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 56.5s.

Lifecycle: open Samples:PostgresCustomers query → save with `.script` provenance → reopen → re-execute query.

Note: share step on this spec emits `Share skipped: ... permissions_user_group_id_fkey` (FK violation when granting to a User without materialized .group); soft-handled (warn, not assert) per spec design.

## Retrospective

### What worked well
- `openTableFromDbQuery` produces df with `.script: 'Var = Samples:PostgresCustomers()'` — reopen re-runs the query against the live DB.

### What did not work
- Nothing notable.

### Suggestions for the platform
- The `permissions_user_group_id_fkey` FK violation when granting permissions to a User (vs Group) repeats across multiple lifecycle specs.

### Suggestions for the scenario
- None from this run.
