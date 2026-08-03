# Collaborative Filtering for Linked Tables — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Run JS script to link 3 tables | 25s | PASS | PASSED | JS API only: 3x openTable + 2x linkTables + 3x addTableView; SPGI(3624), linked1(3624), linked2(224) |
| 2 | Select 5 rows in SPGI | 8s | PASS | PASSED | JS API only: grok.shell.v = v; df.selection.set(i, true) for i in 0..4 |
| 3 | SPGI-linked1 has 9 filtered rows | 8s | PASS | PASSED | Selection-to-filter link propagates; df.filter.trueCount == 9 |
| 4 | Switch to SPGI-linked2 | 8s | PASS | PASSED | JS API only: grok.shell.v = v |
| 5 | Open Filter Panel on SPGI-linked2 | 10s | PASS | PASSED | getFiltersGroup(); [name="viewer-Filters"] present with 5 filter cards |
| 6 | Filter link column 3 to 'v ii' | 10s | PASS | PASSED | JS API only: fg.updateOrAdd categorical; SPGI-linked2 filtered to 148 rows |
| 7 | SPGI-linked1 has 5 filtered rows | 8s | PASS | PASSED | Filter-to-filter link propagates; df.filter.trueCount == 5 |
| 8 | Open Filter Panel on SPGI-linked1 | 10s | PASS | PASSED | getFiltersGroup(); 19 filter cards |
| 9 | PAMPA Classification = 'inconclusive' — 2 rows | 10s | PASS | PASSED | fg.updateOrAdd categorical; df.filter.trueCount == 2 |

**Time** = 2b wall-clock per step (incl. thinking). **Result** = 2b outcome. **Playwright** = 2e outcome.

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 70s |
| grok-browser execution (scenario steps) | 36s |
| Execute via grok-browser (total) | 1m 46s |
| Spec file generation | 17s |
| Spec script execution | 35s |
| **Total scenario run (with model)** | 3m 42s |

## Summary

All 9 steps passed end-to-end on dev: table linking (SELECTION_TO_FILTER and FILTER_TO_FILTER) propagated correctly between SPGI / SPGI-linked1 / SPGI-linked2, and category filters applied via `fg.updateOrAdd` gave the expected row counts (9 → 5 → 2). Playwright replay of the existing spec also passed all 9 softSteps. **Total scenario run (with model)**: 3m 42s.

## Retrospective

### What worked well
- `grok.data.linkTables` with both SELECTION_TO_FILTER and FILTER_TO_FILTER sync types worked on first attempt and propagated correctly across all three linked views.
- Category filters set through `FilterGroup.updateOrAdd({type: CATEGORICAL, column, selected})` gave exact expected counts (148, 5, 2) without needing any DOM interaction.
- Replayed Playwright spec was fully aligned with the run log and passed cleanly headed against dev.

### What did not work
- `grok.shell.tableViews.map(...)` failed because `tableViews` is an iterable, not a JS array — had to switch to `for (const v of grok.shell.tableViews)`. Not a scenario blocker, just a JS API ergonomic sharp edge.
- First `npx playwright test` invocation reported "No tests found" because there is no root `playwright.config.ts` and Playwright's default pattern does not match `*-spec.ts`. Re-running with `-c playwright.config.files-and-sharing.ts` picked the test up. Cost ~14s of wasted wall-clock.

### Suggestions for the platform
- Consider making `grok.shell.tableViews` a real Array (or expose `.map/.filter/.toArray()`) so idiomatic JS traversal works without surprise.
- Document in Test Track or in the spec template how Playwright specs under `public/packages/.../Test Track/**/*-spec.ts` are expected to be discovered (either ship a tiny root `playwright.config.ts` with `testMatch: /-spec\.ts$/` or call it out explicitly in the skill).

### Suggestions for the scenario
- The scenario asserts row counts "9", "5", "2" but does not state the expected count on SPGI-linked2 after step 6 (actual: 148). Adding it would make step 6 verifiable on its own view, not only through its downstream effect on SPGI-linked1.
- Step 1 uses an inline JS script block. Consider moving the dataset paths and link definitions into the scenario's trailing JSON metadata (as other FilterPanel scenarios do) so the setup can be shared and versioned independently of the steps.
