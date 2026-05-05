# Projects / Lifecycle Derived: pivot/aggregate/join + GROK-19103 invariant — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL
**Generated from**: `projects-lifecycle-derived.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 3 | n/a | FAIL | FAILED | save with Data Sync ON: expect(received).toBe(expected) // Object.is equality |
| 2 | Step 4 | n/a | FAIL | FAILED | share via JS API + reopen (recipient-side deferred): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open') |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 1m 13s |
| **Total scenario run (with model)** | 1m 13s |

## Summary

Spec finished with status `PARTIAL` on dev.datagrok.ai (1m 13s). 2 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Step 3**: save with Data Sync ON: expect(received).toBe(expected) // Object.is equality
- **Step 4**: share via JS API + reopen (recipient-side deferred): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 2 step(s) failed:
  - Step 3: save with Data Sync ON: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Step 4: share via JS API + reopen (recipient-side deferred): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')
    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\projects-lifecycle-derived-spec.ts:104:11
```
