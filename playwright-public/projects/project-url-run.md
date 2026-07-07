# Projects / Project URL: deep-link reopen for representative project — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 2m 30s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 2m 30s |
| **Total scenario run** | 2m 30s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 2m 30s.

Deep-link reopen via `page.goto('/p/<id>/<name>')` re-materializes the project's tables via `.script` provenance.

## Retrospective

### What worked well
- URL-based deep link reopen ground-truth verified — table re-materializes by row count, layout restored.
- `tv.dataFrame.rowCount` used as load signal (avoids the `shell.tables.length` Dart-throw post-reopen).

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- None from this run.
