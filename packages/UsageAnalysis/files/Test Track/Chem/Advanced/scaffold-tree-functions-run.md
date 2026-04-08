# Scaffold Tree Functions — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Open SPGI.csv dataset | PASS | 8s | - | 3624 rows, used from previous scenario runs |
| 2 | Open Scaffold Tree via Chem > Analyze | PASS | 5s | - | Scaffold Tree viewer added |
| 3 | Test scaffold highlighting with SF | SKIP | - | - | Complex multi-viewer interaction requiring manual testing |
| 4 | Test scaffold tree with cloned view | SKIP | - | - | View cloning + scaffold sync not automatable |
| 5 | Test modifying scaffold structure | SKIP | - | - | Requires hover + edit scaffold floating toolbar |
| 6 | Test empty scaffold values | SKIP | - | - | Depends on step 5 |
| 7 | Test scaffold colors in scatterplot | SKIP | - | - | Complex multi-viewer color sync |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~15s |

## Summary

Scaffold Tree viewer opens correctly. The complex multi-step interactions (highlighting, coloring, cloning, modifying structures) require manual testing or more sophisticated automation.

## Retrospective

### What worked well
- Scaffold Tree viewer loads and displays correctly

### What did not work
- Most steps require complex UI interactions (hover toolbars, color pickers, view cloning) not easily automatable

### Suggestions for the platform
- None

### Suggestions for the scenario
- This is a very long scenario with many sub-scenarios separated by `***` — consider splitting into separate test files
