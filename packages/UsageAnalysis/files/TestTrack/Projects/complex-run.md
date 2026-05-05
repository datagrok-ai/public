# Projects / Complex (smoke): save, rename, share — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL
**Generated from**: `complex.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 1-2 | n/a | FAIL | FAILED | open demog (file source) and save with Sync ON: expect(received).toBe(expected) // Object.is equality |
| 2 | Step 9 equiv | n/a | FAIL | FAILED | rename project via JS API: expect(received).toBe(expected) // Object.is equality |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 1m 17s |
| **Total scenario run (with model)** | 1m 17s |

## Summary

Spec finished with status `PARTIAL` on dev.datagrok.ai (1m 17s). 2 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Step 1-2**: open demog (file source) and save with Sync ON: expect(received).toBe(expected) // Object.is equality
- **Step 9 equiv**: rename project via JS API: expect(received).toBe(expected) // Object.is equality

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 2 step(s) failed:
  - Step 1-2: open demog (file source) and save with Sync ON: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Step 9 equiv: rename project via JS API: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

  71 |   if (stepErrors.length > 0) {
  72 |     const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
> 73 |     throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
     |           ^
  74 |   }
  75 | });
  76 |
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\complex-spec.ts:73:11
```
