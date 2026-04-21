# Structure Filter — Run Results

**Date**: 2026-04-21
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI.csv and show filter panel | 9s | PASS | PASSED | 3624 rows, ≥1 filter cards |
| 2 | `grok.chem.searchSubstructure(Structure, 'c1ccccc1')` → BitSet | 4s | PASS | PASSED | matchCount > 0 and < 3624; `df.filter.trueCount` matches |
| 3 | `df.filter.setAll(true)` → reset, close + reopen panel | 4s | PASS | PASSED | closed width ≤ 0, reopened width > 0 |
| 4 | Clone view (`grok.shell.addTableView(df)`) + apply benzene filter → both views share the filter | 4s | PASS | PASSED | `views.length ≥ 2`, `df.filter.trueCount > 0` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 30s |
| grok-browser execution (scenario steps) | n/a |
| Execute via grok-browser (total) | 30s |
| Spec file generation | 25s |
| Spec script execution | 24.0s |
| **Total scenario run (with model)** | ~1m 30s |

## Summary

Substructure filtering via `grok.chem.searchSubstructure` works on SPGI.csv (3624 rows): benzene substructure yields a bitset applied through `df.filter.and(bs)`, with consistent counts between bitset.trueCount and df.filter.trueCount. Filter panel close/reopen preserves state; cloning the view and reapplying a filter yields a shared filter across both views (single dataframe). Not automated: DnD column → filter panel, right-click "Current value → Use as filter" (context-menu text-match that lacks stable selectors), the sketcher-driven flow.

## Retrospective

### What worked well
- `grok.chem.searchSubstructure(col, smarts)` is the right automation primitive — avoids the sketcher entirely
- `.first()` locator after a cloned view prevents strict-mode violation on duplicate `.d4-grid[name="viewer-Grid"]`

### What did not work
- Scenario asks for drag-and-drop of a column header into the filter panel — Playwright `dragTo()` requires precise coordinates that are fragile

### Suggestions for the platform
- Add `grok.shell.tv.addFilter(colName)` helper so automation can seed filters without DnD
- "Use as filter" context-menu entry deserves `[name="menu-use-as-filter"]` for click-through automation

### Suggestions for the scenario
- Split "DnD column → filter" into its own scenario; it's a UX-heavy micro-behavior
- The "both cloned views share filter" expectation should be stated explicitly (filter on a dataframe is single-source, so all views reflect the same state)
