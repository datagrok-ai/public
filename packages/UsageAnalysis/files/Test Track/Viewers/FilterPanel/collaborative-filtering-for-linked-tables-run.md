# Collaborative Filtering for Linked Tables — Run Results

**Date**: 2026-04-08
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|-------|-------|
| 1 | Run JS script to link 3 tables | PASS | 15s | PASSED | SPGI(3624), linked1(3624), linked2(224) |
| 2 | Select 5 rows in SPGI | PASS | 2s | PASSED | JS API: selection.set(i, true) |
| 3 | SPGI-linked1 has 9 filtered rows | PASS | 2s | PASSED | Selection-to-filter link works |
| 4 | Switch to SPGI-linked2 | PASS | 1s | PASSED | JS API: grok.shell.v = v |
| 5 | Open Filter Panel on SPGI-linked2 | PASS | 3s | PASSED | getFiltersGroup(), 5 filter cards |
| 6 | Filter link column 3 to v ii | PASS | 2s | PASSED | JS API: fg.updateOrAdd categorical, 148 rows |
| 7 | SPGI-linked1 has 5 filtered rows | PASS | 2s | PASSED | Filter-to-filter link propagation works |
| 8 | Open Filter Panel on SPGI-linked1 | PASS | 3s | PASSED | getFiltersGroup() |
| 9 | PAMPA Classification = inconclusive | PASS | 2s | PASSED | 2 filtered rows as expected |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~35s |
| Spec file generation | ~3s |
| Spec script execution | 30s (PASSED) |

## Summary

All 9 steps passed. Table linking (selection-to-filter and filter-to-filter) works correctly. Filtering propagates between linked tables as expected — selecting 5 rows in SPGI filters SPGI-linked1 to 9 rows, then applying a category filter on SPGI-linked2 further narrows SPGI-linked1 to 5 rows, and adding PAMPA Classification filter results in exactly 2 rows.

## Retrospective

### What worked well
- Table linking APIs (grok.data.linkTables) set up correctly on first attempt
- Selection-to-filter and filter-to-filter sync types both work as documented
- Category filter via fg.updateOrAdd worked reliably across linked tables

### What did not work
- (nothing — all steps passed)

### Suggestions for the platform
- (none)

### Suggestions for the scenario
- The scenario uses inline JS script — consider providing dataset paths in JSON metadata like other scenarios
