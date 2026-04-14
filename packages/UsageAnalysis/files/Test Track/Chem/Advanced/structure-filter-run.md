# Structure Filter — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps

| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open SPGI.csv, open filter panel | PASS | 17s | PASSED | 3624 rows, 90 cols; filter panel with 43 filters including Structure substructure filter |
| 2 | Set substructure filter (benzene) | PASS | 5s | PASSED | `grok.chem.searchSubstructure(col, 'c1ccccc1')` found 1356/3624 matches; filter applied |
| 3 | Disable filter, close/open panel, re-enable | PASS | 5s | PASSED | `df.filter.setAll(true)` reset to 3624; `fg.close()` hid panel; `getFiltersGroup()` reopened |
| 4 | Clone view, verify filter sync | PASS | 8s | PASSED | Cloned view "Table (2)" created; applied benzene filter synced in both views (1356 rows) |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | 40s |
| Spec file generation | 3s |
| Spec script execution | 24s |

## Summary

All 4 core steps passed. The Structure filter panel opened with 43 filters for SPGI.csv. Benzene substructure search matched 1356/3624 molecules. Filter panel close/open cycle worked correctly. Cloned view shared the same filter state via shared DataFrame.

## Retrospective

### What worked well
- `grok.chem.searchSubstructure(col, 'c1ccccc1')` API works for programmatic substructure filtering
- `df.filter.and(bs)` applies the bitset as a filter immediately
- `fg.close()` and `getFiltersGroup()` toggle the filter panel correctly
- Cloned views share the same DataFrame, so filter changes are automatically synced

### What did not work
- The structure filter's `setMolecule()` method doesn't exist — had to use `grok.chem.searchSubstructure` + `df.filter.and()` instead
- Drawing in the structure filter sketcher via automation was not attempted (canvas-based)

### Suggestions for the platform
- The structure filter should expose a `setMolecule(smiles)` API for programmatic testing
- Filter sync status could show an indicator when views share the same filter

### Suggestions for the scenario
- The scenario has multiple sub-sections (4 separate test blocks) — consider splitting into numbered steps
- Step about "Current value > Use as filter" requires right-clicking a cell — not easily automatable
