# Add New Column (Multi-source) — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PARTIAL

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open table from local storage / Home dir / query | PASS | N/A | Opened demog.csv via grok.data.getDemoTable |
| 2 | Click Add new column icon | PASS | N/A | Edit > Add New Column... menu works |
| 3 | Add formula using columns (${HEIGHT} + ${WEIGHT}) | PASS | N/A | Column CalcCol1 created, values correct (e.g., 233.68) |
| 4 | Add another column using previous calc column | PASS | N/A | CalcCol2 = ${CalcCol1} * 2, values correct (e.g., 467.37) |
| 5 | Change column names — formulas should update | PASS | N/A | Renamed HEIGHT to HT; formula auto-updated to ${HT} |
| 6 | Change values — calc columns should recalculate | PASS | N/A | After fireValuesChanged(), CalcCol1 and CalcCol2 recalculated correctly |
| 7 | Save project with datasync | SKIP | N/A | Requires datasync configuration |
| 8 | Close All | SKIP | N/A | Depends on step 7 |
| 9 | Open saved project, rename column in formula | SKIP | N/A | Depends on step 7 |
| 10 | Change values in reopened project | SKIP | N/A | Depends on step 7 |

## Summary

Core Add New Column functionality works correctly: creating calculated columns, chaining formulas, and propagating column renames and value changes. Steps 7-10 (datasync persistence) were skipped as they require server-side datasync configuration.

## Retrospective

### What worked well
- Formula creation and preview work smoothly
- Column rename propagation to formulas is automatic
- Value recalculation through dependency chain works after fireValuesChanged()

### What did not work
- Values did not auto-recalculate immediately on set(); required explicit fireValuesChanged() call

### Suggestions for the platform
- Consider auto-triggering recalculation when column values are modified via the API

### Suggestions for the scenario
- Specify which formula to use in step 3 (currently says "TODO: specify which formula")
- Clarify that datasync steps require specific server configuration
- Split datasync persistence into a separate scenario
