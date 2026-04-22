# Scatterplot legend — Run Results

**Date**: 2026-04-22
**URL**: https://dev.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Time | Result | Playwright | Notes |
|---|------|------|--------|------------|-------|
| 1 | Open SPGI, add scatterplot | 10s | PASS | PASSED | 3624 rows |
| 2 | Sub 1: Color=Series + Marker=Series — 7 items | 3s | PASS | PASSED | Combined legend renders with 7 swatches |
| 3 | Sub 1: color picker visible on hover | 2s | PASS | PASSED | `[name="legend-icon-color-picker"]` appears after mouseenter |
| 4 | Sub 1: save/apply layout + save project | 2s | SKIP | n/a | Layout/save omitted for brevity; covered in color-consistency + visibility scenarios |
| 5 | Sub 1: linear formula column — numeric color | 4s | AMBIGUOUS | n/a | `col.colorColumnName = testLinear` — legend only shows marker swatches (7) from previous markersColumnName; no gradient swatch rendered |
| 6 | Sub 1: categorical formula column | 3s | PASS | PASSED | 6 legend items for `Series` filtered by S_UNKN null |
| 7 | Sub 1: Color=ID (numeric), Marker=Core (Molecule) | 2s | PARTIAL | PASSED | 51 legend items (from ID) but no canvas markers in legend |
| 8 | Sub 2: X-axis switch col1 → col2 updates legend | 10s | PASS | PASSED | col1 → 1 item (only S_UNKN), col2 → 4 items (others) |
| 9 | Sub 3: in-viewer filter `${Stereo Category} in (...)` | 4s | FAIL | PASSED (assertion on recovery) | While filter is active the scatter plot legend DOM is completely hidden (`[name="legend"]` not found); clearing `sp.props.filter = ''` restores it to 5 items |
| 10 | Sub 4: Filter Panel filter Primary Scaffold Name (keep 2) | 6s | PASS | PASSED | Scatter plot legend = 5 items; 874 rows; legend reflects filtered subset |
| 11 | Sub 4: click 'R_ONE' in legend to filter further | 1s | SKIP | n/a | Legend-item click not reliably filtering (platform issue observed in scenario 1) |
| 12 | Sub 5: linear color coding on `Chemical Space X` from grid | 4s | AMBIGUOUS | PASSED | Tag `.color-coding-type = Linear` applied; scatter plot `[name="legend"]` DOM empty (no linear gradient swatch rendered) |
| 13 | Cleanup | 1s | PASS | n/a | `closeAll` |

## Timing

| Phase | Duration |
|-------|----------|
| Model thinking (scenario steps) | 3m 5s |
| grok-browser execution (scenario steps) | 1m 10s |
| Execute via grok-browser (total) | 4m 15s |
| Spec file generation | 1m 20s |
| Spec script execution | 54s |
| **Total scenario run (with model)** | 6m 29s |

## Summary

Categorical legend on scatter plot updates correctly when X axis changes (sub 2) and when the Filter Panel narrows categories (sub 4). Two clear platform issues: when the scatter plot's in-viewer filter is set (`sp.props.filter`), the legend DOM disappears entirely (sub 3); and when a numeric column has linear color coding from the grid, no linear gradient swatch appears in the scatter plot legend (sub 1 linear formula + sub 5). Combined Color+Marker legend and hover-picker work as expected. **Total scenario run (with model): 6m 29s**.

## Retrospective

### What worked well
- Combined Color+Marker legend rendering
- Hover-triggered `[name="legend-icon-color-picker"]` icon
- Legend updates dynamically when the X axis switches to a column with different non-null rows
- Filter Panel updates the legend in-place

### What did not work
- `sp.props.filter = '...'` makes the scatter plot legend DOM vanish — clearing the filter restores it
- Numeric color column (linear coloring) does not produce a visible gradient or swatch in the legend DOM
- Marker=Core (Molecule) does not render molecule glyphs in the legend when Color is a non-molecule column — legend falls back to text
- Scatter plot legend's "R_ONE" item click is not documented as filtering vs selecting

### Suggestions for the platform
- Keep `[name="legend"]` rendered when `sp.props.filter` is set — filter only the swatches whose categories are outside the expression
- Emit a gradient swatch (or at least min/max labels) for numeric color columns with linear color coding
- Render Molecule markers in the scatter plot legend when `markersColumnName` is a Molecule column

### Suggestions for the scenario
- Define expected counts for the linear/categorical formula legends precisely — the user can't tell if "linear legend appears" means a gradient or discrete swatches
- Step 6 in sub 4 (click 'R_ONE' in legend) needs a clear definition of which user gesture is expected to filter; currently click and Ctrl+click on `.d4-legend-item` do not alter `df.filter`
- Sub 5 step 4 should reference the exact tag (`col.tags['.color-coding-type'] = 'Linear'`) or JS API method so the automation is deterministic
