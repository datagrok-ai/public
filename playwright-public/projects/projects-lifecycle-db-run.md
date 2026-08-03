# Projects / Lifecycle DB: NorthwindTest orders source — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Lifecycle DB / Query: Samples:PostgresProducts source | 1m 06s | PASS | PASSED | All softStep blocks completed |
| 2 | Lifecycle DB / Table: Northwind public.products via DbQuery | 1m 24s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 2m 30s (2 tests) |
| **Total scenario run** | 2m 30s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 30s (2 tests).

Two flavors covered:
1. DB / Query — Samples:PostgresProducts (saved query function).
2. DB / Table — Northwind public.products via DbQuery (ad-hoc table reference).

Both produce `.script` provenance compatible with reopen → re-execute.

Note: share step on this spec emits `Share skipped: ... permissions_user_group_id_fkey` (FK violation when granting to a User without materialized .group); soft-handled (warn, not assert) per spec design.

## Retrospective

### What worked well
- `openTableFromDbQuery` and `openTableFromDbTable` both produce reopen-stable `.script` provenance.

### What did not work
- Nothing notable.

### Suggestions for the platform
- The `permissions_user_group_id_fkey` FK violation when granting permissions to a User (vs Group) repeats across multiple lifecycle specs.

### Suggestions for the scenario
- None from this run.
