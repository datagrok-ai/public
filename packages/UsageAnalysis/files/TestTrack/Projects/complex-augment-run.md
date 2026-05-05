# Projects / Complex Augment: addRelation 4 tables via JS API — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 44.3s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 44.3s |
| **Total scenario run** | 44.3s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 44.3s.

## Retrospective

### What worked well
- All softStep blocks passed; `Project.addLink(savedFile)` instance method (replacing the non-existent `dapi.projects.addRelation`) produced the expected 4-link relation count.
- JS API path stable for save → reopen → relation enumeration.

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- None from this run.
