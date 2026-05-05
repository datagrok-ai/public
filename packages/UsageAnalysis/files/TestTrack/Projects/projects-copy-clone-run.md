# Projects / Copy Clone: 3 save modes + GROK-19750 invariant — Run Results

**Date**: 2026-05-05
**URL**: https://dev.datagrok.ai
**Status**: FAIL
**Generated from**: `projects-copy-clone.md` + `grok-browser` references (no MCP run).

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Setup | n/a | FAIL | FAILED | build original demog project with viewers: expect(received).toBe(expected) // Object.is equality |
| 2 | 4b | n/a | FAIL | FAILED | open original, add viewer, Save Copy as <name>-link: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function |
| 3 | 4b verification | n/a | FAIL | FAILED | reopen <name>-link, viewers intact: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator)) |
| 4 | 4b GROK-19750 invariant | n/a | FAIL | FAILED | reopen original, viewers still present: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator)) |
| 5 | 4c | n/a | FAIL | FAILED | open original, add viewer, Save Copy as <name>-clone: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function |
| 6 | 4d | n/a | FAIL | FAILED | open original, add viewer/customization, save as PVC: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function |
| 7 | Step 5 | n/a | FAIL | FAILED | re-share each variant via JS API: expect(received).toBeGreaterThan(expected) |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | n/a |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | n/a (skipped per user request) |
| Spec file generation | (batched across all 20 specs) |
| Spec script execution | 2m 53s |
| **Total scenario run (with model)** | 2m 53s |

## Summary

Spec finished with status `FAIL` on dev.datagrok.ai (2m 53s). 7 soft-step(s) failed.

## Retrospective

### What worked well
- The spec compiled and ran end-to-end (no boilerplate-level breakage).
- Reference-driven selectors load correctly via storageState auth.

### What did not work
- **Setup**: build original demog project with viewers: expect(received).toBe(expected) // Object.is equality
- **4b**: open original, add viewer, Save Copy as <name>-link: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
- **4b verification**: reopen <name>-link, viewers intact: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator))
- **4b GROK-19750 invariant**: reopen original, viewers still present: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator))
- **4c**: open original, add viewer, Save Copy as <name>-clone: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
- **4d**: open original, add viewer/customization, save as PVC: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
- **Step 5**: re-share each variant via JS API: expect(received).toBeGreaterThan(expected)

### Suggestions for the platform
- If failures cluster around Save Project dialog or `grok.dapi.files.readCsv`, consider verifying these JS APIs on the target environment and updating `grok-browser/references/projects.md` accordingly.

### Suggestions for the scenario
- The reference docs alone weren't enough to write a self-validating spec for this scenario; the **Save Project dialog flow** in particular needs an MCP-confirmed selector trace before the spec can be authored verbatim from references.


## Raw error

```
Error: 7 step(s) failed:
  - Setup: build original demog project with viewers: expect(received).toBe(expected) // Object.is equality

Expected: true
Received: false

Call Log:
- Timeout 60000ms exceeded while waiting on the predicate
  - 4b: open original, add viewer, Save Copy as <name>-link: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
    at eval (eval at evaluate (:290:30), <anonymous>:1:22)
    at eval (eval at evaluate (:290:30), <anonymous>:1:45)
    at eval (<anonymous>)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
  - 4b verification: reopen <name>-link, viewers intact: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator))
    at Array.from (<anonymous>)
    at eval (eval at evaluate (:290:30), <anonymous>:1:18)
    at eval (<anonymous>)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
  - 4b GROK-19750 invariant: reopen original, viewers still present: page.evaluate: TypeError: undefined is not iterable (cannot read property Symbol(Symbol.iterator))
    at Array.from (<anonymous>)
    at eval (eval at evaluate (:290:30), <anonymous>:1:18)
    at eval (<anonymous>)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
  - 4c: open original, add viewer, Save Copy as <name>-clone: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
    at eval (eval at evaluate (:290:30), <anonymous>:1:22)
    at eval (eval at evaluate (:290:30), <anonymous>:1:46)
    at eval (<anonymous>)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
  - 4d: open original, add viewer/customization, save as PVC: page.evaluate: TypeError: grok.shell.tv.addViewer is not a function
    at eval (eval at evaluate (:290:30), <anonymous>:1:22)
    at eval (eval at evaluate (:290:30), <anonymous>:1:45)
    at eval (<anonymous>)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
  - Step 5: re-share each variant via JS API: expect(received).toBeGreaterThan(expected)

Expected: > 0
Received:   0
    at eval (eval at evaluate (:290:30), <anonymous>:1:22)
    at eval (eval at evaluate (:290:30), <anonymous>:1:45)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
    at eval (eval at evaluate (:290:30), <anonymous>:1:18)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
    at eval (eval at evaluate (:290:30), <anonymous>:1:18)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anonymous>:1:44)
    at eval (eval at evaluate (:290:30), <anonymous>:1:22)
    at eval (eval at evaluate (:290:30), <anonymous>:1:46)
    at UtilityScript.evaluate (<anonymous>:290:30)
    at UtilityScript.<anonymous> (<anony
```
