# Projects / Complex Move: move project across namespaces via JS API — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL
**Generated from**: `complex-move.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 1-2 | n/a | FAIL | FAILED | open demog and save with Sync ON: expect(received).toBe(expected) // Object.is equality |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 1m 40s |
| **Total scenario run (with model)** | 1m 40s |

## Summary

Spec finished with status `PARTIAL` on dev.datagrok.ai (1m 40s). 1 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Step 1-2**: open demog and save with Sync ON: expect(received).toBe(expected) // Object.is equality

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 1 step(s) failed:
  - Step 1-2: open demog and save with Sync ON: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate

  88 |   if (stepErrors.length > 0) {
  89 |     const summary = stepErrors.map((e) => `  - ${e.step}: ${e.error}`).join('\n');
> 90 |     throw new Error(`${stepErrors.length} step(s) failed:\n${summary}`);
     |           ^
  91 |   }
  92 | });
  93 |
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\complex-move-spec.ts:90:11
```
