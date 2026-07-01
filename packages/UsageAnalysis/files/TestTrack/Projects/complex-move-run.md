# Projects / Complex Move: move project across namespaces via JS API — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | (whole spec) | 1m 18s | PASS | PASSED | All softStep blocks completed |

## Timing

| Phase | Duration |
|-------|----------|
| Spec script execution | 1m 18s |
| **Total scenario run** | 1m 18s |

## Summary

Spec passed end-to-end on dev.datagrok.ai. Total run: 1m 18s.

Wave 2B / complex.md Step 10 — move project across namespaces (default → file-share → Space) via JS API path. Drag-drop and right-click "Move to" UI paths are blocked / non-existent per decision-log b2-2026-05-03-drag-drop-ui-only-reclassification.

## Retrospective

### What worked well
- JS API namespace mutation (`project.namespace = '<target>'; dapi.projects.save(p)`) round-tripped through three namespace targets without errors.

### What did not work
- Nothing notable.

### Suggestions for the platform
- None from this run.

### Suggestions for the scenario
- Drag-drop / right-click "Move to" UI paths remain undocumented in references — keep the current JS-API-only spec until UI selectors are captured.
