# Projects / Complex Integration: heterogeneous sources in one project — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 39.8s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 39.8s |
| **Total scenario run** | 39.8s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 39.8s.

## Retrospective

### What worked well
- Heterogeneous source open (File + Space + DB + Query + Script + Pivot/Aggregate/Join/Clone) saved into a single project, reopens with parent-table relations restored.
- JS API path for save+reopen (canonical `uploadProject` pattern via `saveProjectWithProvenance`).

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- None from this run.
