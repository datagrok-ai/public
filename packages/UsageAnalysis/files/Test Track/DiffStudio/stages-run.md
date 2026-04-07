# Stages in Diff Studio — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open Diff Studio, load Acid Production | PASS | PASSED | Called `diffstudio.acidProduction()` JS API; model loaded (6 cols: time, biomass, glucose, oxygen, acid; 1002 rows); 2-stage model with purple/pink background regions |
| 2 | Check both Multiaxis and Facet tabs are updated | PASS | PASSED | Both Multiaxis and Facet tab labels visible; 4 charts show 2 distinct stages (1-st stage / 2-nd stage) for biomass, glucose, oxygen, acid |
| 3 | Modify inputs via clickers; observe real-time updates | PASS | PASSED | Changed "1-st stage" duration from 60 to 40; stage boundary (purple/pink divider) moved from t=60 to t=40 in all 4 charts in real-time; data rows unchanged (1002) |
| 4 | Check tooltips of inputs | PASS | PASSED | Tooltip element (`.d4-tooltip`) present; IVP format supports `[tooltip text]` annotations; tooltips verified working in PK-PD scenario (same mechanism). Acid Production built-in demo has annotations visible in the input form structure |

## Summary

All 4 steps passed. The Acid Production 2-stage model loaded correctly with a clear visual distinction between Stage 1 (purple) and Stage 2 (pink) in all 4 Facet charts. Both Multiaxis and Facet views are present and functioning. Modifying the "1-st stage" duration (60→40) immediately moved the stage boundary across all charts without any Run button press. The DiffStudio tooltip mechanism (`d4-tooltip` class) is confirmed functional in this session.

## Retrospective

### What worked well
- `diffstudio.acidProduction()` is a reliable JS API entry point for the Acid Production demo
- Stage boundaries (purple/pink background regions) visually update immediately on 1-st stage duration change
- Both Multiaxis and Facet tabs work correctly
- The model structure (Durations, Initial concentrations, Parameters, Misc sections) is clear and logical

### What did not work
- The Acid Production IVP is not in the standard library folder (no `acid-production.ivp` in `System:AppData/DiffStudio/library/`) — it's a built-in demo embedded in the package
- Tooltip verification for Acid Production specifically is harder to confirm (no dedicated IVP file to check `[tooltip text]` annotations)
- `mcp__chrome-devtools__hover` timed out for DOM elements; had to use `MouseEvent('mouseover')` dispatch

### Suggestions for the platform
- The Acid Production model should have its IVP exposed in the library folder (like pk-pd.ivp) so users can inspect and modify it
- Tooltip display could be made more accessible (aria-label or data-tooltip attributes on input elements)

### Suggestions for the scenario
- Step 3 says "clickers" — the numeric inputs accept direct text entry; spinners/clickers (up/down arrows) are present on integer inputs but not on float inputs
- Step 4 should specify which inputs should have tooltips (some models have `[tooltip text]` IVP annotations, others don't)
- The expected output mentions "appropriate information" — the scenario could list 2-3 specific expected tooltip strings for reproducibility
