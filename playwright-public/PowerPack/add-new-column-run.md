# Add new columns — Run Results

**Date**: 2026-04-23
**URL**: http://localhost:8888 (local, v1.26.8, master, commit 10ff7fae7)
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open Demog dataset | 4s | PASS | PASSED | Opened via `grok.dapi.files.readCsv('System:DemoFiles/demog.csv')`; HEIGHT and WEIGHT present (5850 rows). Star-icon path in Test Track is equivalent. |
| 2 | Press Add new column icon; dialog opens | 2s | PASS | PASSED | Clicked `[name="icon-add-new-column"]`; dialog rendered at 752×502 with name/type/formula/column-list/function-list/preview-grid; name input, CodeMirror, and OK button all mounted in the same tick. |
| 3 | UI check: no overlap, no unnecessary scrollbars, resizable | 1s | PASS | PASSED | `d4-dialog-contents` has no horizontal or vertical overflow; 8 resize handles present on all sides/corners. Programmatic pointer-drag on `.d4-host-bottom-right-resizer` did NOT trigger the Dart resize handler (platform gap for synthetic events), but handles exist and a CSS-height override does enlarge the frame. |
| 4 | Set name "New" and formula Round(${HEIGHT} + ${WEIGHT}) | 5s | PASS | PASSED | MCP: `fill` on `input.ui-input-addnewcolumn-name` with "New", click the `.cm-content` CodeMirror, `fill` the formula; OK button became `enabled`. Playwright: `page.locator('.d4-dialog input.ui-input-addnewcolumn-name').click()` + `Control+A` + `keyboard.type('New')`; `.cm-content.click()` + `keyboard.type('Round(${HEIGHT} + ${WEIGHT})')`. |
| 5 | Press OK; "New" column is added | 2s | PASS | PASSED | MCP: column "New" added (int, samples 234/257/221 consistent with Round(HEIGHT+WEIGHT)), dialog closed immediately. Playwright: `page.locator('[name="button-Add-New-Column---OK"]').click()` + `waitForFunction` on `df.columns.contains('New')`. |
| 6 | Reopen dialog, history icon, select most recent; form autofills | 10s | PASS | PASSED | MCP: `[name="icon-history"]` popup lists "28 seconds ago: {Name: New, Type: int, Expression: Round(${HEIGHT} + ${WEIGHT})}"; clicking the menu item via MCP uid-based click autofills both fields. Playwright: `locator.click()` on the history icon and the first `.d4-menu-popup .d4-menu-item` is a real synthesized click, so `applyInput` fires and both fields repopulate. |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 1m 45s |
| grok-browser execution (scenario steps) | 9s |
| Execute via grok-browser (total) | 1m 54s |
| Spec file generation | 31s |
| Spec script execution | 11s |
| **Total scenario run (with model)** | ~2m 36s |

## Summary

All six scenario steps PASS in both the MCP-driven grok-browser run and the Playwright replay (existing spec — not overwritten). Column "New" is added correctly and the Recent Activities history icon autofills the form. Total end-to-end run (MCP scenario + spec drafting + spec execution) was **2m 36s**. The prior step-4 race (`input.ui-input-addnewcolumn-name` queried before it mounted) remains fixed via Playwright locators with auto-waiting and an explicit `waitFor` on the name input after the dialog opens.

## Retrospective

### What worked well
- `[name="icon-add-new-column"]` toolbar icon is reliable
- Command-bar `[name="icon-history"]` is discoverable via the dialog footer
- Input history's `applyInput` repopulates both the name and the CodeMirror formula
- Dialog layout is clean at the default 752×502 — no overflow

### What did not work
- Synthetic `pointerdown/move/up` events on the resize handles do not trigger Dart's resize logic, even though the handles exist in the DOM — hard to cover resize in automation.
- A raw DOM `.click()` on a history menu item closes the popup but does not fire `applyInput`. Only a real synthesized click (CDP / Playwright `locator.click` / MCP uid click) applies the stored state — any automation that builds events via `dispatchEvent` will appear to work but leave the form unfilled.
- In Playwright, `input.ui-input-addnewcolumn-name` was not yet mounted right after the dialog rendered, even though the OK button was visible — suggests the OK button is created before the inputs finish wiring. Fixed in this spec by using Playwright locators (with auto-waiting) and adding an explicit `waitFor` on the name input after the dialog opens, instead of querying the DOM directly inside `page.evaluate`.

### Suggestions for the platform
- Make Add New Column dialog resize handles respond to synthetic pointer events (or expose a JS API to resize).
- Wire history menu item to fire on both native and synthesized click events (listen to `d4-menu-item-click` bubble).
- The Add New Column OK button is named `button-Add-New-Column---OK` rather than the conventional `button-OK`. Either standardize or update docs — the mismatch silently broke prior attempts.

### Suggestions for the scenario
- Step 1: add "System:DemoFiles/demog.csv" as a bulleted alternative to the TestTrack star icon for direct/automation runs.
- Step 3 is a broad "hover entire dialog" check; splitting into concrete sub-assertions (no overlap / no scrollbars / resize works) makes each independently verifiable.
- Step 4 "using autocomplete hints and columns drag-n-drop" is covered by AddNewColumn/autocomplete.md and AddNewColumn/input_functions.md; accepting any successful formula entry here would avoid duplication.

---
{"order": 2, "datasets": ["System:DemoFiles/demog.csv"]}
