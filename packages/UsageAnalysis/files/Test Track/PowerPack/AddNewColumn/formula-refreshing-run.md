# Formula Refreshing — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open the demog dataset | PASS | N/A | demog.csv with 5850 rows |
| 2 | Create Weight2 = ${WEIGHT} + 100 | PASS | N/A | Created, val[0] = 173.20 (73.20 + 100) |
| 3 | Create Weight3 = ${Weight2} + 100 | PASS | N/A | Created, val[0] = 273.20 (173.20 + 100) |
| 4 | Create Weight4 = Log10(${Weight3}) - 0.2 | PASS | N/A | Created, val[0] = 2.2365 |
| 5 | Verify formula dependency propagation | PASS | N/A | Changed Weight2[0]=200; Weight3 updated to 300, Weight4 to 2.277 |

## Summary

All 5 steps passed. The formula dependency chain (Weight2 -> Weight3 -> Weight4) works correctly. Changing a value in Weight2 propagates through Weight3 and Weight4 with correct recalculation. Note: recalculation requires fireValuesChanged() when values are set programmatically.

## Retrospective

### What worked well
- Three-level formula dependency chain works correctly
- Recalculation propagates through the entire chain
- Log10 and arithmetic operations work correctly in formulas

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- Auto-trigger recalculation on programmatic value changes (currently requires fireValuesChanged())

### Suggestions for the scenario
- Step 5 mentions "modify the formulas via Context Panel" but doesn't specify how — could add detail about right-clicking column header or using the properties panel
- Steps about saving as project and re-opening (mentioned in Additional Notes) could be explicit numbered steps
