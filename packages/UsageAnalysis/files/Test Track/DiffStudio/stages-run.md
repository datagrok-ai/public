# Stages in Diff Studio — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open DiffStudio, load Acid Production | PASS | 10s | PASSED | Library > "Acid production"; GA-production model loaded with 6 cols, 1002 rows |
| 2 | Check Multiaxis and Facet plots | PASS | 1s | PASSED | Both tabs present; 11 canvases visible; Facet shows 4 panels (biomass, glucose, oxygen, acid) with stage labels |
| 3 | Modify inputs; observe real-time changes | PASS | 4s | PASSED | Changed 1st stage 60→40 via `input[name="input-1-st-stage"]`; charts updated with shifted stage boundary |
| 4 | Check tooltips on input hover | PASS | 3s | PASSED | "1-st stage" → "Duration of the 1-st stage"; "overall" → "Overall duration" |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 25s |
| Spec file generation | 3s |
| Spec script execution | 17s |

## Summary

All 4 steps passed. The Acid Production (GA-production) model loaded from Library with 2-stage simulation (1st stage + 2nd stage). Both Multiaxis and Facet tabs work with distinct colors and stage boundary labels. Input modification updates charts in real-time. Tooltips display correctly on input label hover via Datagrok's custom tooltip system.

## Retrospective

### What worked well
- Library menu shows "Acid production" (not "Acid Production" or "GA Production") — matched correctly
- Facet view shows 4 panels with clearly labeled stage boundaries ("1st stage", "2nd stage")
- Input names use `input-{param}` pattern consistently
- Tooltip system triggered by `mouseover` event on `label.ui-label` elements; tooltip text found in `.d4-tooltip` element

### What did not work
- The scenario mentions a "Count" input but the Acid Production model uses "1-st stage" and "overall" for durations — no count input exists

### Suggestions for the platform
- Stage boundary colors on the Facet plots could be customizable
- Tooltip display delay could be shortened for automation scenarios

### Suggestions for the scenario
- Step 1 says "Acid Production" but the Library menu shows "Acid production" — minor casing difference
- Step 3 mentions "clickers" and "Count input" in expected results, but the model has no Count input — update to reference "1-st stage" or "overall" duration inputs
- Step 4 could list expected tooltip texts for verification
