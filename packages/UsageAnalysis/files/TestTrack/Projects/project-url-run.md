# Projects / Project URL: deep-link reopen for representative project — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL
**Generated from**: `project-url.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup | n/a | FAIL | FAILED | create representative file-share project: expect(received).toBe(expected) // Object.is equality |
| 2 | Step 3 equivalent | n/a | FAIL | FAILED | derive deep-link URL for the project: expect(received).toBeGreaterThan(expected) |
| 3 | Step 4-5 | n/a | FAIL | FAILED | open URL in new tab and verify project loads: expect(received).toBe(expected) // Object.is equality |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 1m 56s |
| **Total scenario run (with model)** | 1m 56s |

## Summary

Spec finished with status `PARTIAL` on dev.datagrok.ai (1m 56s). 3 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Setup**: create representative file-share project: expect(received).toBe(expected) // Object.is equality
- **Step 3 equivalent**: derive deep-link URL for the project: expect(received).toBeGreaterThan(expected)
- **Step 4-5**: open URL in new tab and verify project loads: expect(received).toBe(expected) // Object.is equality

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 3 step(s) failed:
  - Setup: create representative file-share project: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Step 3 equivalent: derive deep-link URL for the project: expect(received).toBeGreaterThan(expected)

Expected: > 0
Received:   0
  - Step 4-5: open URL in new tab and verify project loads: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

  58 |   if (stepErrors.length > 0) {
  59 |     const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
> 60 |     throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
     |           ^
  61 |   }
  62 | });
  63 |
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\project-url-spec.ts:60:11
```
