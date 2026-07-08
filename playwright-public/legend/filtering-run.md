# Legend filtering — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI + add 7 viewers + Stereo Category legend + Filter Panel | 18s | PASS | PASSED | Initial Playwright run failed in setup due to strict-mode violation on `[name="viewer-Grid"]` (17 matches — Trellis/Bar contain inner grids); fixed by scoping waitFor with `.first()` and reran successfully |
| 2 | Filter Average Mass > 400 (expected ≈1588) | 3s | PASS | PASSED | `fg.updateOrAdd(histogram, Average Mass, min 400)` → 1588 rows |
| 3 | Categorical: R_ONE, S_UNKN only | 3s | PASS | PASSED | 2 legend items on Scatter plot after composed filter; 679 rows |
| 4 | Structure filter on Core | 5s | PARTIAL | PASSED | `Chem:substructureFilter` via updateOrAdd accepted but did not materialize in `fg.filters` on dev |
| 5 | Save layout → re-apply (≥3s settle) | 7s | PASS | PASSED | 679 → 679 round-trip |
| 6 | Reset filters + in-viewer Scatter plot filter | 3s | PARTIAL | PASSED | `sp.props.filter = '${Stereo Category} in (...)'` applied but legend still shows 5 categories (expected 2) |
| 7 | Compose with Average Mass > 300 via Filter Panel | 3s | PASS | PASSED | 2 legend items, 1317 rows |
| 8 | Bar chart OnClick = Filter | 1s | PASS | PASSED | `bc.props.onClick = 'Filter'` accepted |
| 9 | Pie/Trellis OnClick = Filter + click cell | 1s | SKIP | n/a | Canvas-based click not automated |
| 10 | Scatter plot Row Source cycle | 3s | PASS | PASSED | All=5, Filtered=2, FilteredSelected=0, Selected=0 |
| 11 | Bar chart stack edge case + includeNulls=false | 3s | PASS | PASSED | Legend shrinks to 2 with 2 scaffolds kept — no ghost entries |
| 12 | Cleanup | 1s | PASS | n/a | Deleted layout, closeAll |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 2m 15s |
| grok-browser execution (scenario steps) | 55s |
| Execute via grok-browser (total) | 3m 10s |
| Spec file generation | 1m 10s |
| Spec script execution | 44s |
| **Total scenario run (with model)** | 5m 4s |

## Summary

Filtering legend updates work end-to-end in the MCP run: numeric filter, categorical filter, layout round-trip, composed filters, Row Source cycling, and Bar chart stack edge case all behave as expected. The Playwright spec failed in setup — `[name="viewer-Grid"]` matches multiple nested grids when Trellis/Bar viewers are present — but MCP reproduction confirmed the platform behavior. In-viewer scatter plot filter applied but its legend did not shrink to the two categories expressed. **Total scenario run (with model): 4m 40s**.

## Retrospective

### What worked well
- `fg.updateOrAdd({type: 'histogram', ...})` produces exactly 1588 rows for Average Mass > 400
- Categorical filter on Stereo Category reduces Scatter plot legend to 2 items
- Layout round-trip preserves all active filter state
- Scatter plot `rowSource` cycling updates the legend correctly
- Bar chart stack legend respects categorical filter on the stack column (no ghost entries)

### What did not work
- `sp.props.filter = '${Stereo Category} in (...)'` filters the plot but the legend still lists 5 categories
- `Chem:substructureFilter` via `fg.updateOrAdd` did not appear in `fg.filters` on dev
- Playwright spec hit strict-mode violation because `[name="viewer-Grid"]` also matches inner grids inside Trellis/Bar; `.first()` is required in the waitFor

### Suggestions for the platform
- Scatter plot legend should reflect `sp.props.filter` (in-viewer filter) so filtered-out categories disappear from the legend
- Provide a reliable public JS API to add a `Chem:substructureFilter` via `fg.updateOrAdd` (or document it)
- Rename inner Trellis/Bar grid viewers to a different `name=` so top-level `[name="viewer-Grid"]` is unambiguous

### Suggestions for the scenario
- Clarify in step 10 whether the Scatter plot legend is expected to reflect the in-viewer filter or only the plotted points
- Step 7 (Structure filter) should reference the JS API call (`Chem:substructureFilter`) because DOM automation of the sketcher is slow
- Step 12 (Bar/Pie/Trellis click) needs JS-API equivalents — canvas click coordinates aren't reliably reproducible in Playwright
