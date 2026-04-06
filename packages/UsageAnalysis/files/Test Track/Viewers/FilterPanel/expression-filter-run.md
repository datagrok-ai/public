# Expression Filter — Run Results

**Date**: 2026-04-05
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|-------|-------|
| 1 | Open demog, open Filter Panel | PASS | PASSED | 5850 rows, 9 filter cards |
| 2 | Add Expression filter | PASS | PASSED | Via fg.updateOrAdd |
| 3 | Add 5 rules (WEIGHT>50, HEIGHT<160, SEX=F, RACE contains an, STARTED after 01/01/1991) | PASS | PASSED | 288 rows with AND |
| 4 | Switch to OR | PASS | PASSED | 5850 rows |
| 5 | Switch to AND | PASS | PASSED | 288 rows |
| 6 | Remove first rule (WEIGHT>50) | PASS | PASSED | 330 rows with 4 rules |
| 7 | Switch to free-text mode | PASS | PASSED | Mode changed |
| 8 | Enter AGE > 30 and SEX = M | PASS | PASSED | Added as free-text rule |
| 9 | Enter AGE < 60 and HEIGHT > 190 | PASS | PASSED | 0 rows (contradictory AND) |
| 10 | Uncheck first 4 rules | PASS | PASSED | 73 rows |
| 11 | Save layout | PASS | PASSED | Via JS API |
| 12 | Close Filter Panel | PASS | PASSED | Panel closed |
| 13 | Apply saved layout | PASS | PASSED | 73 rows, first 4 unchecked |

## Summary

Expression filter works correctly in all tested modes. Adding rules, switching AND/OR, removing rules, free-text mode, unchecking rules, and layout persistence all pass.

## Retrospective

### What worked well
- Expression filter rules added and managed entirely via JS API (fg.updateOrAdd with gridNames/gridValues/mode/expressionMode)
- AND/OR toggle produces expected row counts
- Free-text mode rules work alongside expression mode rules
- Layout save/restore preserves unchecked rule state

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- None identified

### Suggestions for the scenario
- Step 12 says "Save the layout (Ctrl+S)" but Ctrl+S saves project, not layout. Should use Toolbox > Layouts > Save or JS API
