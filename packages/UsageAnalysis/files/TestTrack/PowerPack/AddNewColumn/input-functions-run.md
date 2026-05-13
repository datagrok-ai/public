# Input Functions — Run Results

**Date**: 2026-04-22
**URL**: http://localhost:8888 (local)
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open spgi.csv dataset | 7s | PASS | PASSED | SPGI.csv loaded (3624 rows, 88 cols; Structure has `semType=Molecule`). |
| 2 | Open Add New Column dialog | 1s | PASS | PASSED | Dialog renders, 198 functions listed alphabetically. |
| 3 | Hover function, click + → function with param types inserted | 1s | PASS | PASSED | `Abs(num)` inserted into `.cm-content` via DOM `.click()` on `[name="icon-plus"]` inside row. |
| 4 | Clear, same function via drag-n-drop | 2s | PASS | PASSED | MCP `drag` (real CDP events) from Abs row to `.cm-content` inserts `Abs(num)`. Playwright reproduces via `page.mouse.down/move*20/up`. |
| 5 | Click Structure column → getCLogP + icon → `Chem:getCLogP(${Structure})` | 3s | PASS | PASSED | Canvas synthetic pointer events DO select the column (previous run claim that they don't was wrong). Column list re-sorts — `canonicalize`, `convertMolNotation`, `convertMoleculeNotation`, `getCLogP`, `getDescriptors` rank top. Plus icon on `getCLogP` inserts formula with `${Structure}` pre-filled. |
| 6 | Clear, same getCLogP via drag-n-drop | 2s | PASS | PASSED | MCP `drag` from `span[name="span-getCLogP"]` to `.cm-content` inserts `Chem:getCLogP(${Structure})`. Initial failure was a flaky one-off (three back-to-back drag attempts returned empty `cm-content`); on rerun from a fresh dialog it works consistently — even immediately after a synthetic mouse-sequence attempt (which correctly does nothing, because Dart's `makeDraggable` needs real CDP pointer events to cross the 5-px threshold). |
| 7 | Select numeric column, add Abs via + → `Abs(${NumericCol})` | 3s | PASS | PASSED | Clicked canvas row 18 (`Chemical Space Y`). `Abs(${Chemical Space Y})` inserted. Spec asserts formula matches `Abs(${Chemical Space X/Y})`. |
| 8 | Clear the text field | 1s | PASS | PASSED | `execCommand('selectAll') + execCommand('delete')` clears `.cm-content`. |
| 9 | Click Id column → sort icon → By name → functions alphabetical | 2s | PASS | PASSED | Sort dropdown `[name="div-By-name"]` works. First 10 functions: Abs, Acos, Add, And, Asin, Atan, Atan2, Avg, BinByDateTime, BinBySpecificLimits — alphabetical. |
| 10 | + on Abs with Id (string) selected → column NOT passed | 1s | PASS | PASSED | `Abs(num)` inserted — Id (string) doesn't match Abs(num) input, so the column is correctly not substituted. |

**Time** column = step 2b wall-clock. **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | ~3m |
| grok-browser execution (scenario steps) | ~30s |
| Execute via grok-browser (total) | ~3m 30s |
| Spec file generation | ~1m |
| Spec script execution | 18.8s |
| **Total scenario run (with model)** | ~5m |

## Summary

All 10 scenario steps pass end-to-end — both in interactive MCP driving and in the Playwright replay (18.8s).
A single initial MCP `drag` attempt at step 6 returned an empty formula; reproduction on a fresh dialog
succeeded every time (even immediately after a synthetic-event attempt that correctly failed to trigger Dart's
drag). Treating that as a flaky/timing artifact rather than a platform bug. Spec replay passed on first run.

**Total scenario run (with model)**: ~5m (including diagnostic re-run of step 6).

## Retrospective

### What worked well
- `[name="icon-plus"]` scoped to a function row — reliable selector for the `+` icon.
- Synthetic `MouseEvent` on the column grid **canvas** (`mousedown` → `mouseup` → `click` at `{row_height * idx + row_height/2, row/2}`) DOES select a column. The previous run's claim that "canvas-based column list doesn't respond to synthetic pointer events" was incorrect for the column grid — the click DID re-sort the function list.
- `page.mouse.down` + 20-step `page.mouse.move` + `page.mouse.up` replays Dart's custom drag protocol. Unlike HTML5 DnD (which never worked), this native CDP pointer sequence triggers `makeDraggable`'s `onMouseDown` + distance-threshold handler.

### What did not work
- Synthetic JS `MouseEvent` sequences (`mousedown` on the draggable span + `mousemove` on `document` + `mouseup` on `.cm-content`) never trigger `makeDraggable` — Dart's handler listens through the Dart event bridge and appears to ignore `isTrusted=false` events. Only real CDP pointer events (via MCP `drag` or Playwright `page.mouse.*`) start the drag.
- One earlier MCP `drag` run returned empty on three consecutive retries before recovering. Could not reproduce from a fresh dialog. Noting as flaky; not a confirmed platform bug.

### Suggestions for the platform
- Add a semantic label (`aria-label` or a hidden `<div>` overlay) on each column-grid row so automation/accessibility tooling can target rows by column name rather than by computed y-offset.
- The `+` icon row behaviour (type-matched substitution when parameter type matches selected column) is valuable — document it in the Add New Column help page, since today users have to discover it by trial.
- `ui.makeDraggable` (`core/client/d4/lib/src/widgets/ui.dart:541`) only accepts trusted pointer events. If possible, also react to synthesised events in development/test builds so unit-test automation can drive the feature without full CDP.

### Suggestions for the scenario
- Step 3 could name a specific function (e.g., `Abs`) instead of "any function" — keeps replay deterministic.
- Step 7 could say "select *any numeric column*, e.g. `Chemical Space X`" to avoid ambiguity about which column the parameter binding is meant for.
- Steps 4 and 6 (drag-and-drop) should note: the drag is a pointer-based Dart interaction, not HTML5 DnD — test runners must use real CDP pointer sequences.

---
{"order": 6}
