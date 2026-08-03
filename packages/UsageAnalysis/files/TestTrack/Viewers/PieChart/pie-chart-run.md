# Pie chart tests (Playwright) — Run Results

**Date**: 2026-04-07
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| **Sorting** | | | | |
| 1 | Category = RACE | PASS | PASSED | |
| 2 | Pie Sort Type = by value | PASS | PASSED | |
| 3 | Sort Order = desc | PASS | PASSED | |
| 4 | Sort Order = asc | PASS | PASSED | |
| 5 | Sort Type = by category | PASS | PASSED | |
| 6 | Category sort asc | PASS | PASSED | |
| 7 | Category sort desc | PASS | PASSED | |
| **Segment angle and length** | | | | |
| 1 | Category = RACE | PASS | PASSED | |
| 2 | Segment Angle Column = AGE | PASS | PASSED | |
| 3 | Angle Aggr Type = sum | PASS | PASSED | |
| 4 | Angle Aggr Type = count | PASS | PASSED | |
| 5 | Segment Length Column = WEIGHT | PASS | PASSED | |
| 6 | Length Aggr Type = max | PASS | PASSED | |
| 7 | Clear both columns | PASS | PASSED | |
| **Appearance** | | | | |
| 1 | Start Angle = 90 | PASS | PASSED | |
| 2 | Start Angle = 180 | PASS | PASSED | |
| 3 | Start Angle = 0 | PASS | PASSED | |
| 4 | Max Radius = 100 | PASS | PASSED | |
| 5 | Max Radius = 150 | PASS | PASSED | |
| 6 | Shift = 10 (exploded) | PASS | PASSED | |
| 7 | Shift = 0 | PASS | PASSED | |
| **Labels** | | | | |
| 1 | Label Position = Inside | PASS | PASSED | |
| 2 | Label Position = Outside | PASS | PASSED | |
| 3 | Label Position = Auto | PASS | PASSED | |
| 4 | Show Label off | PASS | PASSED | |
| 5 | Show Percentage off | PASS | PASSED | |
| 6 | Show Value on | PASS | PASSED | |
| 7 | Re-enable label + percentage | PASS | PASSED | |
| **Outline** | | | | |
| 1 | Outline Line Width = 5 | PASS | PASSED | |
| 2 | Outline Line Width = 0 | PASS | PASSED | |
| 3 | Outline Line Width = 1 | PASS | PASSED | |
| **Include nulls** | | | | |
| 1 | Category = DIS_POP | PASS | PASSED | |
| 2 | includeNulls = true | PASS | PASSED | |
| 3 | includeNulls = false | PASS | PASSED | |
| **Column selector** | | | | |
| 1 | Show Column Selector off | PASS | PASSED | |
| 2 | Show Column Selector on | PASS | PASSED | |
| **Legend** | | | | |
| 1 | Legend Visibility = Always | PASS | PASSED | |
| 2 | Legend Position = LeftTop | PASS | PASSED | Uses PascalCase (no space) |
| 3 | Legend Position = RightBottom | PASS | PASSED | |
| 4 | Legend Visibility = Never | PASS | PASSED | |
| 5 | Legend Visibility = Auto | PASS | PASSED | |
| **Category map (dates)** | | | | |
| 1 | Category = STARTED | PASS | PASSED | |
| 2 | Default categoryMap = year | PASS | PASSED | |
| 3 | categoryMap = month | PASS | PASSED | |
| 4 | categoryMap = quarter | PASS | PASSED | |
| **Row source** | | | | |
| 1 | Select first 50 rows | PASS | PASSED | |
| 2 | Row Source = Selected | PASS | PASSED | |
| 3 | Row Source = Filtered | PASS | PASSED | |
| 4 | Row Source = All | PASS | PASSED | |
| **Aggregation functions** | | | | |
| 1-8 | avg/min/max/sum/med/stdev/count | PASS | PASSED | All 7 aggr types work |
| **Title and description** | | | | |
| 1 | Show Title = true | PASS | PASSED | |
| 2 | Title = "Demographics" | PASS | PASSED | |
| 3 | Description = "By race" | PASS | PASSED | |
| 4 | Description Position = Top | PASS | PASSED | |
| 5 | Description Visibility = Never | PASS | PASSED | |
| 6 | Clear title | PASS | PASSED | |
| **Layout persistence** | | | | |
| 1-4 | Set cat=RACE, angle=AGE, startAngle=45, shift=5 | PASS | PASSED | |
| 5 | Save layout | PASS | PASSED | |
| 6 | Close pie chart | PASS | PASSED | |
| 7 | Apply saved layout | PASS | PASSED | All 4 props restored correctly |
| 8 | Verify restored props | PASS | PASSED | cat=RACE, angle=AGE, startAngle=45, shift=5 |
| 9 | Delete layout | PASS | PASSED | |
| **Selection and interaction** | | | | |
| 1 | Category = RACE | PASS | PASSED | |
| 2 | Click slice → 5267 rows selected | PASS | PASSED | Caucasian slice |
| 3 | showSelectedRows toggle | PASS | PASSED | |
| 4 | showMouseOverRowGroup toggle | PASS | PASSED | |
| **On Click modes** | | | | |
| 1-2 | onClick=Select, click → selected | PASS | PASSED | 5267 selected |
| 3-4 | onClick=Filter, click → filtered | PASS | PASSED | 5267 filtered |
| 5 | Click empty → filter cleared | PASS | PASSED | 5850/5850 restored |
| **Selection between grid and pie** | | | | |
| 1 | Select 50 in grid | PASS | PASSED | |
| 2 | Click pie slice → grid updated | PASS | PASSED | 5267 selected |
| 3 | Clear selection | PASS | PASSED | |
| **Auto layout** | | | | |
| 1 | Auto Layout off | PASS | PASSED | |
| 2 | Margin Left = 50 | PASS | PASSED | |
| 3 | Margin Top = 50 | PASS | PASSED | |
| 4 | Auto Layout on | PASS | PASSED | |
| **Table switching (SPGI)** | | | | |
| 1 | Add pie on demog | PASS | PASSED | |
| 2 | Switch table to SPGI | PASS | PASSED | |
| 3 | Switch back to demog | PASS | PASSED | |
| 4 | Row Source = Selected (100 rows) | PASS | PASSED | |
| 5-6 | Row Source = Filtered + RACE=Asian | PASS | PASSED | 72 rows filtered |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 32s |
| grok-browser execution (scenario steps) | 8s |
| Execute via grok-browser (total) | 40s |
| Spec file generation | 7s |
| Spec script execution | 47s |
| **Total scenario run (with model)** | 1m 34s |

All rows are full-phase wall-clock (incl. model thinking and retries). The two `scenario steps`
rows sum to `Execute via grok-browser (total)`. Spec was not rewritten (existing spec is correct
per the user's instruction); the 7s reflects verification/sanity-pass only. Spec execution used
the existing `playwright.config.files-and-sharing.ts` (testMatch: `/-spec\.ts$/`).

## Summary

All 16 pie chart test sections passed on dev.datagrok.ai. All viewer properties (sorting, segment angle/length, appearance, labels, outline, include nulls, column selector, legend, category map, row source, aggregation, title/description, layout persistence, selection, on click modes, auto layout, table switching) behave correctly. Layout save/restore preserves all custom settings. Canvas click interaction works for selection and filtering modes.

## Retrospective

### What worked well
- All pie chart properties accessible through `viewer.props`
- Layout persistence correctly saves and restores all 4 tested properties
- Canvas click interaction reliably selects pie slices
- legendPosition uses PascalCase (LeftTop, RightBottom) — discovered and corrected during run

### What did not work
- No failures encountered

### Suggestions for the platform
- legendPosition values should accept both "Left Top" and "LeftTop" for consistency

### Suggestions for the scenario
- Scenario says "Left Top" / "Right Bottom" but JS API requires "LeftTop" / "RightBottom"
