# Basic Operations — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

### Section 1: Filtering and resetting

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 0 | Open spgi-100 + Filter Panel | 28s | PASS | PASSED | `readCsv` + `addTableView` + `getFiltersGroup()`; 100 rows × 90 cols, 42 filter cards |
| 1.1-1.3 | Sketch c1ccccc1 → 32 rows | 17s | PASS | PASSED | `.sketch-link` click → SMILES via native value setter + Enter → `[name="button-OK"]`; 32 rows |
| 1.4-1.5 | Stereo Category R_ONE → 15 rows | 9s | PASS | PASSED | JS API only (canvas categories): `fg.updateOrAdd({type: CATEGORICAL, column: 'Stereo Category', selected: ['R_ONE']})` |
| 1.6-1.7 | Average Mass max=400 → 4 rows | 7s | PASS | PASSED | JS API: `fg.updateOrAdd({type: 'histogram', column: 'Average Mass', min: col.min, max: 400})` |
| 1.8 | Disable all filters globally → 100 | 7s | PASS | PASSED | `.d4-filter-group-header input[type="checkbox"]` click |
| 1.9 | Re-enable all filters → 4 | 7s | PASS | PASSED | Same global checkbox toggled back |
| 1.10 | Disable Stereo Category → 9 | 10s | PASS | PASSED | Per-card checkbox click; rows = 9 (Structure + AvgMass active) |
| 1.11 | Close Filter Panel → 100 | 8s | PASS | PASSED | `getFiltersGroup().close()`; all rows shown |
| 1.12 | Reopen — Stereo still disabled | 9s | PASS | PASSED | `stereoChecked === false` confirmed; rows = 9 (Structure + AvgMass survive) |
| 1.13 | Reset (↺) → 100 | 14s | PASS | PASSED | `.fa-arrow-rotate-left` click; no confirm dialog appeared this run, reset still applied |
| 1.14-1.15 | Close/reopen — default state | 11s | PASS | PASSED | 100 rows, 42 cards |
| 1.16-1.17 | Remove Structure & Core (X) | 9s | PASS | PASSED | `[name="icon-times"]` scoped to filter card; 40 cards remain |
| 1.18 | Close/reopen — removed absent | 12s | PASS | PASSED | 40 cards, no Structure/Core |

### Section 2: Adding Filters

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 2.1 | Hamburger → Remove All → 0 cards | 13s | PASS | PASSED | Found titlebar hamburger via `.panel-titlebar [name="icon-font-icon-menu"]`; menu item `Remove All`; OK confirm not always shown but Remove All applied |
| 2.2 | Add ID (drag) → histogram | 4s | PASS | PASSED | JS API only (drag is canvas): `fg.updateOrAdd({type: 'histogram', column: 'Id'})` |
| 2.3 | Add CAST Idea ID (col header menu) | 4s | PASS | PASSED | JS API equivalent |
| 2.4 | Add Structure (right-click cell) | 4s | PASS | PASSED | `fg.updateOrAdd({type: 'Chem:substructureFilter', column: 'Structure', columnName: 'Structure'})` |
| 2.5 | Add Scaffold Tree (Filter Panel ctx menu) | 4s | PASS | PASSED | `fg.updateOrAdd({type: 'Chem:scaffoldTreeFilter', ...})` |
| 2.6 | Verify order: ScaffoldTree, Structure, CAST, Id | 6s | PASS | PASSED | `fg.filters` enumeration: [Structure (scaffold), Structure (substructure), CAST Idea ID, Id] |
| 2.7 | Hamburger → Remove All → 0 cards | 12s | PASS | PASSED | Same hamburger flow |
| 2.8 | Close Filter Panel | 7s | PASS | PASSED | `getFiltersGroup().close()`; panel removed |

### Section 3: Hidden Columns

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 3.1-3.3 | Hide Structure/Core/R1 → grid 86 visible cols | 12s | PASS | PASSED | `grid.columns.byName(...).visible = false` + `grid.invalidate()` |
| 3.4 | Open Filter Panel — no hidden cards | 8s | PASS | PASSED | 39 cards, none for Structure/Core/R1 |
| 3.5 | Restore columns → 89 visible | 9s | PASS | PASSED | `visible = true` for all three |
| 3.6-3.7 | Hamburger Remove All + close | 25s | PASS | PASSED | 0 cards, panel closed |
| 3.8 | Reopen — Structure/Core/R1 cards present | 9s | PASS | PASSED | 42 cards, all three present |

**Time** = step 2b wall-clock per step (incl. model thinking + waits). **Result** = step 2b outcome. **Playwright** = step 2e outcome (existing spec ran without modification).

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 53s |
| grok-browser execution (scenario steps) | 1m 34s |
| Execute via grok-browser (total) | 4m 27s |
| Spec file generation | 9s |
| Spec script execution | 50s |
| **Total scenario run (with model)** | 6m 9s |

`Execute via grok-browser (total)` = end of step 2b – start of step 2b. The two `scenario steps` rows (thinking + grok-browser execution) sum to that total. `Spec script execution` = `npx playwright test ... --headed` wall-clock (47.1s test + ~3s startup/teardown). `Total scenario run` = end of step 2e – start of step 2b.

## Summary

Ran basic-operations end-to-end against dev. All 31 scenario steps passed in the MCP browser phase (Section 1: structure → 32, Stereo R_ONE → 15, Mass max 400 → 4, global toggle, per-filter toggle, close/reopen persistence, reset, remove via X. Section 2: hamburger Remove All, add 4 filters via JS API, verify order. Section 3: hide columns, verify exclusion from filter panel, restore, verify reappearance). Existing `basic-operations-spec.ts` was left untouched per user instruction; re-running it produced 22 PASSED softSteps in 47.1s. **Total scenario run (with model)**: 6m 9s.

## Retrospective

### What worked well
- All filter operations completed successfully on the first try this run — no retries needed
- Disabled Stereo Category survived close/reopen as expected (DOM checkbox state preserved)
- JS API additions (`fg.updateOrAdd`) reliably created all four filter types in correct top-of-list order (ScaffoldTree → Structure substructure → CAST Idea ID → Id)
- Hidden columns correctly excluded from the filter panel after `Remove All` + reopen — and reappear when columns are re-shown
- Hamburger menu flow (`.panel-titlebar [name="icon-font-icon-menu"]` → `.d4-menu-item-label` text match → optional confirm OK) is robust

### What did not work
- Reset (↺) icon click did not show a confirmation dialog this run — reset was applied directly with no prompt. Spec still works because it tries to click OK if the dialog appears (no-op otherwise)
- Drag-from-grid, column-header-menu Add Filter, right-click cell → Use as filter, and Filter Panel right-click → Add filter (steps 2.2-2.5) all rely on canvas / context menu interactions that are not automatable; spec falls back to `fg.updateOrAdd` for all four

### Suggestions for the platform
- Consider exposing a `name=` attribute on the Filter Panel header reset icon (currently found via `.fa-arrow-rotate-left` class) for more stable automation
- Reset confirmation dialog inconsistency (sometimes shown, sometimes not) — clarify the trigger condition; if intentional (e.g. only when filters are dirty), document it

### Suggestions for the scenario
- Steps 2.2-2.5 should explicitly note that drag/right-click/context-menu adds are functionally equivalent to JS API `fg.updateOrAdd` calls — useful for both manual testers (knowing what to expect) and automation authors
- Step 1.10 should specify expected row count after disabling Stereo Category (~9 rows for Structure + Average Mass intersection)
- Step 1.13 wording "confirm the prompt" implies a dialog always appears — observed behavior is that the prompt does not always show; either the wording should be "and confirm if prompted" or the platform behavior should be made consistent
