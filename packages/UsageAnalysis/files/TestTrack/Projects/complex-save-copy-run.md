# Projects / Complex Save Copy: round-trip Save Copy with sync OFF/ON — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: FAIL
**Generated from**: `complex-save-copy.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Step 1-2 | n/a | FAIL | FAILED | open demog and save baseline (Sync ON): expect(received).toBe(expected) // Object.is equality |
| 2 | Step 3 | n/a | FAIL | FAILED | Save Copy without Data Sync (use Save dialog with new name): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open') |
| 3 | Step 4 | n/a | FAIL | FAILED | close and reopen NoSync copy: page.evaluate: TypeError: Cannot read properties of undefined (reading 'open') |
| 4 | Step 5-6 | n/a | FAIL | FAILED | re-save the copy under Sync name: page.click: Timeout 15000ms exceeded. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 1m 27s |
| **Total scenario run (with model)** | 1m 27s |

## Summary

Spec finished with status `FAIL` on dev.datagrok.ai (1m 27s). 4 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Step 1-2**: open demog and save baseline (Sync ON): expect(received).toBe(expected) // Object.is equality
- **Step 3**: Save Copy without Data Sync (use Save dialog with new name): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')
- **Step 4**: close and reopen NoSync copy: page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')
- **Step 5-6**: re-save the copy under Sync name: page.click: Timeout 15000ms exceeded.

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 4 step(s) failed:
  - Step 1-2: open demog and save baseline (Sync ON): expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - Step 3: Save Copy without Data Sync (use Save dialog with new name): page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')
    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
  - Step 4: close and reopen NoSync copy: page.evaluate: TypeError: Cannot read properties of undefined (reading 'open')
    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
  - Step 5-6: re-save the copy under Sync name: page.click: Timeout 15000ms exceeded.
Call log:
  - waiting for locator('button:has-text("SAVE"), .ui-btn:has-text("SAVE")')

    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
    at eval (eval at evaluate (:290:30), <anonymous>:3:17)
    at async <anonymous>:316:30
    at D:\Рабочая папка\Datagrok\Bitbucket\reddata\public\packages\UsageAnalysis\files\TestTrack\Projects\complex-save-copy-spec.ts:77:11
```
