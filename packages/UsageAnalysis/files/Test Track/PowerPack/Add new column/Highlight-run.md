# Highlight — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open demog.csv dataset | PASS | N/A | Already open |
| 2 | Open Add new column dialog | PASS | N/A | Edit > Add New Column... |
| 3 | Paste Abs(${age}) — column name highlighted blue | PASS | N/A | ${age} has cm-column-name class, case-insensitive match to AGE column |
| 4 | Paste Avg($[age]) — column name highlighted blue | PASS | N/A | $[age] also gets cm-column-name class highlighting |

## Summary

All 4 steps passed. Column names in both ${} (scalar) and $[] (aggregate) syntax are highlighted with the cm-column-name CSS class (renders in blue). Column name matching is case-insensitive.

## Retrospective

### What worked well
- Both column reference syntaxes (${} and $[]) are highlighted correctly
- Case-insensitive matching works (typed "age", matched column "AGE")

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- N/A

### Suggestions for the scenario
- Step 3 has a typo: "text filed" should be "text field"
- Step numbers have duplicates (two step 3s)
