# EDA ANOVA — Run Results

**Date**: 2026-03-10
**URL**: https://public.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open demog.csv from Demo files | PASS | 5850 rows, 11 cols (USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, ...) |
| 2 | ML > Analyze > ANOVA | PASS | Dialog opens via Eda:anova(); Category=RACE (auto), Feature=AGE (auto), Alpha=0.05 |
| 3 | Click RUN | PASS | Box plot opens immediately showing AGE distribution by RACE (Asian, Black, Caucasian, Other) |
| 4 | Verify box plot | PASS | Box plot shows 4 categories with quartiles, whiskers, and individual data points |
| 5 | Verify Analysis tab | PASS | Table shows: Source of variance / SS / DF; Between groups SS=848.48, Within groups SS=1069045.00, Total SS=1069893.38 |
| 6 | Click F-test tab | PASS | Shows: Null Hypothesis (all group means equal), Alternative Hypothesis (at least one differs), Conclusion: "No significant differences detected", Learn more link |

## Summary

All steps passed. ANOVA dialog opens with good defaults (Category=RACE, Feature=AGE, Alpha=0.05). After clicking RUN, a box plot viewer appears with AGE distribution by race category. The Analysis tab shows the ANOVA table with SS and DF values. The F-test tab correctly interprets the results: no significant age difference across race groups in this dataset.

## Retrospective

### What worked well
- Dialog auto-populates Category and Feature with reasonable defaults
- Box plot shows all 4 racial categories clearly with individual data points overlaid
- Analysis tab presents a clear ANOVA table
- F-test tab gives a human-readable conclusion with hypothesis statements
- No console errors during execution

### What did not work
- Nothing — all steps passed cleanly

### Suggestions for the platform
- No major issues found

### Suggestions for the scenario
- Scenario is concise and accurate
- Could add expected values (e.g., p-value threshold, expected conclusion) to make verification more objective
- Note that the dialog auto-selects Category=RACE and Feature=AGE — users don't need to configure manually for demog.csv
