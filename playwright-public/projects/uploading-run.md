# Projects / Uploading: representative source-matrix subset (Cases 1, 8, 9 — Sync ON) — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Case 1: Files + Files (Sync ON) | 54.1s | PASS | PASSED | All softStep blocks completed |
| 2 | Case 8: Files + Pivot Table > Add (Sync ON) | 57.2s | PASS | PASSED | All softStep blocks completed |
| 3 | Case 9: DB + Aggregate Rows > Add (Sync ON) | 1m 06s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 2m 57s (3 tests) |
| **Total scenario run** | 2m 57s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 57s (3 tests).

Representative source-matrix subset:
- Case 1: pure Files + Files via canonical `openTableFromFile` × 2.
- Case 8: Files + Pivot Table → Add to workspace (derived).
- Case 9: DB + Aggregate Rows → Add to workspace (derived).

All saved with Data Sync ON via `saveProjectWithProvenance` and round-tripped via reopen.

## Retrospective

### What worked well
- `helpers/openers.ts:openTableFromFile` (dot-form normalization workaround for bug 2a) and `addAggregateToWorkspace` produce provenance-stable dataframes.
- `saveProjectWithProvenance` mirrors the canonical `UITests/src/gui/gui-utils.ts:uploadProject` pattern.

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- Cases 4-6 (Spaces sources) remain env-skipped due to bug 1 (Spaces addEntity UUID parse).
