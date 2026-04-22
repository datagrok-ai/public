# Combined Boolean Filter — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open demog | 9s | PASS | PASSED | 5850 rows, 11 columns |
| 2 | Open the Filter Panel | 8s | PASS | PASSED | Panel opened with default filters (9 filter cards) |
| 3 | Add new column SEX_bool via formula | 6s | PASS | PASSED | Bool column: 3243 true (F), 2607 false (M) |
| 5 | Verify CONTROL and SEX_bool in Combined Boolean filter | 13s | PASS | PASSED | Close+reopen panel, `bool-columns` filter added; CONTROL=39, SEX_bool=3243 |
| 7 | Apply filter: CONTROL=true, SEX_bool=false in OR mode | 10s | PASS | PASSED | `grok_GridFilterBase_ApplyState` OR mode → 2632 rows (matches expected) |
| 8 | Apply RACE=Asian + AGE 50-89 | 6s | PASS | PASSED | 8 filtered rows |
| 9 | Save layout via Toolbox > Layouts > SAVE | 7s | PASS | PASSED | Layout saved via `grok.dapi.layouts.save()` |
| 10 | Close the Filter Panel | 7s | PASS | PASSED | Panel closed via `.grok-font-icon-close` |
| 11 | Apply saved layout — verify filter state and row count | 9s | PASS | PASSED | 8 rows restored, all filters reapplied |
| 12 | Remove all filters via Hamburger > Remove All | 9s | PASS | PASSED | All filters removed, 5850 rows |
| 13 | Close the Filter Panel | 6s | PASS | PASSED | Panel closed |
| 14 | Open Filter Panel — Combined Boolean auto-added | 8s | PASS | PASSED | `bool-columns` present automatically on reopen |
| 15 | Apply saved layout — verify filter state and row count | 8s | PASS | PASSED | 8 rows restored correctly |

**Time** = step 2b wall-clock per step (incl. model thinking + waits). **Result** = step 2b outcome. **Playwright** = step 2e outcome (existing spec ran without modification).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 12s |
| grok-browser execution (scenario steps) | 25s |
| Execute via grok-browser (total) | 2m 37s |
| Spec file generation | 12s |
| Spec script execution | 24s |
| **Total scenario run (with model)** | 4m 41s |

`Execute via grok-browser (total)` = end of step 2b – start of step 2b. The two `scenario steps` rows (thinking + grok-browser execution) sum to that total. `Spec script execution` = `npx playwright test ... --headed` wall-clock (17.9s test + ~6s startup/teardown). `Spec file generation` measured the time to re-derive the spec from the run log for comparison; the existing spec was not overwritten per user instruction. `Total scenario run` = end of step 2e – start of step 2b.

## Summary

Ran combined-boolean-filter end-to-end against dev. All 13 numbered scenario steps passed in the MCP browser phase: SEX_bool calculated column created with 3243 true / 2607 false values, Combined Boolean filter picked up CONTROL (39 true) and SEX_bool (3243 true) after panel close/reopen, OR-mode filter (CONTROL=true OR SEX_bool=false) returned exactly 2632 rows, combination with RACE=Asian + AGE 50-89 narrowed to 8 rows, layout save/restore preserved all filter state across close/reopen cycles, Remove All cleared filters back to 5850 rows, and Combined Boolean auto-added on panel reopen when boolean columns exist. Existing `combined-boolean-filter-spec.ts` was left untouched per user instruction; re-running it produced all softSteps PASSED in 17.9s (22.5s total). **Total scenario run (with model)**: 4m 41s.

## Retrospective

### What worked well
- Combined Boolean filter auto-adds when boolean columns exist and panel is reopened
- Layout save/restore preserved all filter states including combined boolean settings
- `grok_GridFilterBase_ApplyState` worked reliably for boolean checkbox manipulation
- Formula-based column creation (`addNewCalculated`) worked smoothly
- Existing spec re-ran cleanly against dev in under 25 seconds end-to-end (login + setup + 11 softSteps + cleanup)

### What did not work
- Combined Boolean filter was not auto-added on initial filter panel open — needed explicit `fg.add({type: 'bool-columns'})` or panel close/reopen
- Viewer title bar hamburger required walking up the DOM to `panel-base` level — not directly accessible from viewer element

### Suggestions for the platform
- Auto-add Combined Boolean filter when boolean columns exist on initial filter panel open
- Make the hamburger menu accessible via a named element within the viewer element itself

### Suggestions for the scenario
- Step 4 is missing from the numbered sequence (jumps from 3 to 5)
- Step 9 says "Save the layout via Toolbox → Layouts → SAVE" — clarify this means the SAVE button in Layouts section, not Ctrl+S
- Step 12 says "Save the layout (Ctrl+S)" in Expression filter but step 9 here says "via Toolbox → Layouts → SAVE" — inconsistent across scenarios
