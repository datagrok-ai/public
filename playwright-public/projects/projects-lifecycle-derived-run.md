# Projects / Lifecycle Derived: pivot/aggregate/join + GROK-19103 invariant — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 1m 30s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 30s |
| **Total scenario run** | 1m 30s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 30s.

GROK-19103 invariant verified — Aggregate Rows via menu lands the result in the active project. Save+reopen restores derived table.

Note: share step on this spec emits `Share skipped: ... permissions_user_group_id_fkey` (FK violation when granting to a User without materialized .group); soft-handled (warn, not assert) per spec design.

## Retrospective

### What worked well
- Aggregate Rows menu path produces a derived dataframe with the expected `Aggregate("...", ...)` `.script` provenance.
- Save → reopen round-trip preserves the derived table.

### What did not work
- Nothing notable.

### Suggestions for the platform
- The `permissions_user_group_id_fkey` FK violation when granting permissions to a User (vs Group) repeats across multiple lifecycle specs. Either `dapi.permissions.grant` should auto-resolve User → User.group, or the API should reject User input with a clearer error.

### Suggestions for the scenario
- None from this run.
