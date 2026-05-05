# Projects / Uploading: representative source-matrix subset (Cases 1, 8 — Sync ON) — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: FAIL
**Generated from**: `uploading.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Case 1 | n/a | FAIL | FAILED | two file-share tables + Save Sync ON: expect(received).toBe(expected) // Object.is equality |
| 2 | Case 1 reopen verification | n/a | FAIL | FAILED | expect(received).toBeGreaterThan(expected) |
| 3 | Case 8 | n/a | FAIL | FAILED | Files + Pivot Table (Add to workspace): expect(received).toBe(expected) // Object.is equality |
| 4 | Case 8 reopen verification | n/a | FAIL | FAILED | expect(received).toBeGreaterThan(expected) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 2m 26s |
| **Total scenario run (with model)** | 2m 26s |

## Summary

Spec finished with status `FAIL` on dev.datagrok.ai (2m 26s). 4 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Case 1**: two file-share tables + Save Sync ON: expect(received).toBe(expected) // Object.is equality
- **Case 1 reopen verification**: expect(received).toBeGreaterThan(expected)
- **Case 8**: Files + Pivot Table (Add to workspace): expect(received).toBe(expected) // Object.is equality
- **Case 8 reopen verification**: expect(received).toBeGreaterThan(expected)

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 4 step(s) failed:
  - Case 1: two file-share tables + Save Sync ON: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Case 1 reopen verification: expect(received).toBeGreaterThan(expected)

Expected: > 0
Received:   0
  - Case 8: Files + Pivot Table (Add to workspace): expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Case 8 reopen verification: expect(received).toBeGreaterThan(expected)

Expected: > 0
Received:   0

  121 |   if (stepErrors.length > 0) {
  122 |     const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
> 123 |     throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
      |           ^
  124 |   }
  125 | });
  126 |
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\uploading-spec.ts:123:11
```
