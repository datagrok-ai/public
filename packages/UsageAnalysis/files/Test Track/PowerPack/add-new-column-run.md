# Add New Columns — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open the Demog Dataset via JS API | PASS | 5s | PASSED | Opened via grok.dapi.files.readCsv, 5850 rows, 11 columns |
| 2 | Press "Add new column" icon, dialog opens | PASS | 2s | PASSED | Clicked [name="icon-add-new-column"], dialog-Add-New-Column appeared with name input, CodeMirror editor, OK button |
| 3 | UI Check: no overlapping text, no overflow | PASS | 3s | PASSED | No overflow; CodeMirror editor visible; dialog rect 752x502 |
| 4 | Add column "New" with Round(${HEIGHT}+${WEIGHT}), press OK | PASS | 4s | PASSED | Column added, type auto-detected as int. Values correct: Round(160.48+73.20)=234 |
| 5 | Recent Activities: reopen dialog, click history, autofill | PASS | 8s | PASSED | History icon showed entry; clicking popup menu item autofilled name="New" and formula="Round(${HEIGHT} + ${WEIGHT})" |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~22s |
| Spec file generation | ~2s |
| Spec script execution | 5.3s |

## Summary

All 5 steps passed on dev.datagrok.ai. The Add New Column dialog opens correctly, displays cleanly with no visual issues, accepts formula input via the CodeMirror editor, creates the column with correct computed values, and the history/recent activity feature restores previous formulas via autofill.

## Retrospective

### What worked well
- Dialog layout is clean with no overlapping or overflow issues
- Formula preview updates in real time with correct computed values
- Column type auto-detection works (detected int for Round result)
- History icon stores and recalls previous formulas with full autofill (name + formula)

### What did not work
- JS `dispatchEvent(click)` on history popup menu items does not trigger autofill; only Playwright's native `.click()` (via element UID in MCP, or Playwright locator) works — likely because Dart menu items rely on pointer events that `dispatchEvent` doesn't fully replicate
- First attempt at clicking history menu item accidentally hit a top-level Edit menu item (`.d4-menu-item-vert` matched globally); scoping to `.d4-menu-popup` was required

### Suggestions for the platform
- History menu items should respond to programmatic click events (dispatchEvent) the same way as native clicks
- Consider adding `name=` attributes to history popup menu items for easier targeting

### Suggestions for the scenario
- Step 1 mentions pressing a "star" icon in TestTrack — not available in automation; opening via JS API is the practical approach
- Step 4 says "using autocomplete hints and columns drag-n-drop" — drag-n-drop is hard to test via automation; consider separating into its own step or marking as manual-only
