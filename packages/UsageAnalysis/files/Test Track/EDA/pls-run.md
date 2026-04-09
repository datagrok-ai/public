# PLS — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv dataset | PASS | 4s | PASSED | Opened via JS API, 30 rows, 17 columns |
| 2 | Run PLS (ML > Analyze > PLS) | PASS | 3s | PASSED | Dialog opened via UI; Predict=price, Using=15 cols, Components=3 |
| 3 | Click RUN | PASS | 2s | PASSED | PLS executed successfully |
| 4 | Verify PLS1, PLS2, PLS3 columns added | PASS | 1s | PASSED | Three new columns added (17 → 20 cols): PLS1, PLS2, PLS3 |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 10s |
| Spec file generation | 2s |
| Spec script execution | 4s |

## Summary

All 4 steps passed. PLS dialog opened with correct defaults (Predict=price, Using=15 numeric columns, Components=3). After clicking RUN, three new columns PLS1, PLS2, PLS3 were added to the dataset as expected. Note: the previous run on public.datagrok.ai (2026-03-10) reported FAIL for this step — this has been fixed on dev.

## Retrospective

### What worked well
- PLS dialog auto-populated Predict=price and selected 15 numeric columns
- Components defaulted to 3 (matching scenario requirement)
- PLS1, PLS2, PLS3 columns correctly added to the table
- Menu navigation (ML → Analyze → PLS) worked reliably

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- No issues found

### Suggestions for the scenario
- Scenario is accurate and matches actual behavior on dev
