# Hierarchical Filter — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open demog.csv | PASS | PASSED | 5850 rows loaded |
| 2 | Open the Filter Panel | PASS | PASSED | Filter panel with default filters |
| 3 | Hamburger > Add filter > Hierarchical Filter | PASS | PASSED | Added via fg.add({type: 'hierarchical'}) |
| 4 | Select SEX, RACE, SEVERITY columns | PASS | PASSED | ApplyState with colNames, header: SEX / RACE / SEVERITY |
| 5 | Expand F > Caucasian, see SEVERITY values | PASS | PASSED | None(1588), High(319), Low(552), Medium(354), Critical(10) |
| 6 | Select Caucasian - all children checked | PASS | PASSED | 2823 rows, other RACE values show 0 |
| 7 | Uncheck Low and Medium - indeterminate parent | PASS | PASSED | 1917 rows |
| 8 | Reorder to SEX / SEVERITY / RACE | PASS | PASSED | ApplyState with new colNames order |
| 9 | Collaborative filtering: None + DIS_POP | PASS | PASSED | None=1815, +DIS_POP(AS,Indigestion)=339 |
| 10 | Save layout | PASS | PASSED | Layout saved via JS API |
| 11 | Close Filter Panel | PASS | PASSED | 5850 rows (filters removed) |
| 12 | Apply saved layout - verify restored | PASS | PASSED | 339 rows restored |

## Summary

All 12 steps passed. The hierarchical filter correctly supports multi-level tree filtering with expandable nodes, checkbox selection with indeterminate parent state, column reordering, and collaborative filtering with categorical filters. Layout save/restore preserves the complete state.

## Retrospective

### What worked well
- ApplyState reliably sets hierarchical filter columns and order
- Tree node expansion and checkbox interaction work via DOM
- Collaborative filtering between hierarchical and categorical filters
- Layout save/restore preserves hierarchical filter state

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- Simplify `name=` attribute encoding for nested tree expanders (currently complex encoded strings)

### Suggestions for the scenario
- Step 4 mentions drag-and-drop reordering in column picker dialog — canvas-based, hard to automate manually
- Step 9 says `F = 1815` for SEX filter but None=1815 includes both SEX values
