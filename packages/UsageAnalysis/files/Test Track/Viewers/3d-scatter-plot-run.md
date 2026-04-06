# 3D Scatter Plot — Run Results

**Date**: 2026-04-07
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open the demog dataset | PASS | PASSED | 5850 rows, 11 columns loaded |
| 2 | Click 3d Scatter plot in Viewers tab | PASS | PASSED | Viewer opened as rotatable 3D model with axes |
| 3a | Click on single data point highlights grid row | PASS | PASSED | Clicked point selected row 3975, tooltip shown with row data |
| 3b | Zoom in/out with mouse wheel | PASS | PASSED | Zoom in/out via wheel events worked correctly |
| 4a | Modify Color property to SEX | PASS | PASSED | setOptions({colorColumnName: 'SEX'}) applied, verified via getOptions() |
| 4b | Modify Size property to AGE | PASS | PASSED | setOptions({sizeColumnName: 'AGE'}) applied, verified via getOptions() |
| 4c | Change Color to RACE and Size to WEIGHT | PASS | PASSED | Both properties changed and verified, points colored by race with varying sizes |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~20s |
| Spec file generation | ~3s |
| Spec script execution | 10s |

## Summary

All steps passed. The 3D Scatter Plot opens correctly as a rotatable 3D model, data point clicking selects the corresponding row in the grid, mouse wheel zoom works, and property modifications (Color, Size) are correctly reflected in the visualization via `setOptions()` with `colorColumnName` and `sizeColumnName` keys.

## Retrospective

### What worked well
- Viewer opened reliably via toolbox icon click
- Data point click interaction selected the correct row in the grid
- `setOptions({colorColumnName: ..., sizeColumnName: ...})` works for changing properties
- `getOptions().look.colorColumnName` / `.sizeColumnName` for reading current values
- Zoom via mouse wheel worked correctly
- Updated spec to use CDP connection for consistency

### What did not work
- `v3d.props.color` and `v3d.setOptions({color: 'RACE'})` do not work — must use `colorColumnName` and `sizeColumnName` flat keys
- `v3d.setOptions({look: {colorColumnName: ...}})` (nested under look) does not apply — must use flat keys

### Suggestions for the platform
- 3D scatter plot should support `color` and `size` shorthand property names like other viewers do

### Suggestions for the scenario
- Step 3 says "Clicking on the single data on the viewer should highlight same data line in the grid" — clearer as "Clicking a data point selects the corresponding row in the grid"
- Consider specifying which properties to modify in step 4 for deterministic testing
