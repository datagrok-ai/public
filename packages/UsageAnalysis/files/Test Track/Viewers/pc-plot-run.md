# PC Plot tests (Playwright) — Run Results

**Date**: 2026-04-10
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Menu Ribbon: Add PC Plot via toolbox icon | PASS | PASSED | `[name="icon-pc-plot"]` click → viewer appeared |
| 2 | Right-click > To Script > To JavaScript | PASS | PASSED | Context menu navigation → balloon with script |
| 3 | Close viewer, re-add via Toolbox | PASS | PASSED | `pc.close()` then toolbox click |
| 4 | Normalize Each Column default (true) | PASS | PASSED | `pc.props.normalizeEachColumn === true` |
| 5 | Disable normalization | PASS | PASSED | Props set to false, verified |
| 6 | Re-enable normalization | PASS | PASSED | Props back to true |
| 7 | Context menu Y Axis > Global | PASS | PASSED | Menu navigation worked |
| 8 | Context menu Y Axis > Normalized | PASS | PASSED | Reverted to normalized |
| 9 | Add AGE to log columns | PASS | PASSED | `logColumnsColumnNames = ['AGE']` |
| 10 | Add WEIGHT to log columns | PASS | PASSED | Both columns in log list |
| 11 | Remove log columns | PASS | PASSED | Empty array restored |
| 12 | Selection: showCurrentLine toggle | PASS | PASSED | Props toggle works |
| 13 | Selection: showMouseOverLine toggle | PASS | PASSED | Props toggle works |
| 14 | Selection: showAllLines toggle | PASS | PASSED | Props toggle works |
| 15 | Context menu Selection toggles | PASS | PASSED | Menu items found and clickable |
| 16 | Style: lineWidth, currentLineWidth, mouseOverLineWidth | PASS | PASSED | All set and verified |
| 17 | Style: labelsOrientation, minMaxOrientation | PASS | PASSED | Set to Vert, verified |
| 18 | Style: horzMargin, autoLayout | PASS | PASSED | Margin=60, autoLayout=false |
| 19 | showFilteredOutLines toggle | PASS | PASSED | Property toggled |
| 20 | Reset View via context menu | PASS | PASSED | Menu item found and clicked |
| 21 | Filter > Show Filters via context menu | PASS | PASSED | Submenu item found |
| 22 | Filter panel: apply AGE range | PASS | PASSED | 2870 rows after filter |
| 23 | Filter panel: reset filter | PASS | PASSED | 5850 rows restored |
| 24 | Column management: remove HEIGHT | PASS | PASSED | Axis removed |
| 25 | Column management: add HEIGHT back | PASS | PASSED | Axis restored |
| 26 | Column management: reorder | PASS | PASSED | WEIGHT first |
| 27 | Density: default circles | PASS | PASSED | `densityStyle === 'circles'` |
| 28 | Density: box plot mode | PASS | PASSED | All box plot toggles work |
| 29 | Density: violin plot mode | PASS | PASSED | Bins and whisker adjustable |
| 30 | Color: set to AGE, log, invert | PASS | PASSED | All color props work |
| 31 | Color: RACE categorical + legend positions | PASS | PASSED | Legend moves correctly |
| 32 | Color: grid color coding (linear/conditional) | PASS | PASSED | JS API color coding applied |
| 33 | Title and description | PASS | PASSED | Set and cleared |
| 34 | Pick Up / Apply | PASS | PASSED | Second plot inherited all props |
| 35 | Layout save/restore | PASS | PASSED | Scatter removed, PC Plot preserved |
| 36 | Table switching to spgi-100 | PASS | PASSED | Table changed, transformation applied |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~8 min |
| Spec file generation | ~5s |
| Spec script execution | 38.3s |

## Summary

All 36 steps across 12 sections passed. PC Plot properties, context menus, filter interaction, column management, density styles, color coding, Pick Up/Apply, layout persistence, and table switching all work correctly on dev.

## Retrospective

### What worked well
- All PC Plot properties are accessible and settable via JS API
- Context menu navigation via DOM events works reliably
- Pick Up / Apply correctly transfers all viewer settings
- Layout save/restore preserves PC Plot state

### What did not work
- Range slider drag on axes cannot be automated (canvas-based) — used property verification instead

### Suggestions for the platform
- Add name attributes to range slider handles for better automation

### Suggestions for the scenario
- Steps involving axis range slider drag should note that programmatic verification is needed
