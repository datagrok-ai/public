# Linear Regression — Run Results

**Date**: 2026-04-09
**URL**: https://dev.datagrok.ai
**Status**: PASS

## Steps
| # | Step | Result | Time | Playwright | Notes |
|---|------|--------|------|------------|-------|
| 1 | Open cars.csv | PASS | 1s | PASSED | 30 rows, 17 cols loaded via JS API |
| 2 | Open Train Model via ML > Models > Train Model | PASS | 3s | PASSED | Menu navigation: ML > Models > Train Model, PredictiveModel view opened |
| 3 | Set Predict = price | PASS | 1s | PASSED | UI column selector: clicked editor, typed "price", pressed Enter |
| 4 | Set Features (all except price and model) + Train | PASS (JS API fallback) | 2s | PASSED | Canvas-based column grid checkboxes cannot be toggled via DOM; fell back to `eda:trainLinearRegression` JS API |

## Timing
| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~45s |
| Spec file generation | 2s |
| Spec script execution | 6.8s |

## Summary
Linear Regression trained successfully on cars.csv predicting price. Steps 1-3 were completed via UI (menu navigation, column selector). Step 4 (Features selection) required JS API fallback because the Select Columns dialog uses a canvas-based grid where checkboxes cannot be toggled via DOM mouse events. The "All"/"None" links work but apply globally, not per-filtered-row. Model trained via `eda:trainLinearRegression` as fallback.

## Retrospective
### What worked well
- ML > Models > Train Model menu navigation works reliably via `mouseenter`/`mousemove` events on submenus
- Predict column selector popup opens on `mousedown` and responds to keyboard type + Enter
- `eda:trainLinearRegression` trains successfully with df + predictColumn params
- Select Columns dialog "All"/"None" links and Search filter work via UI clicks

### What did not work
- Canvas-based checkboxes in Select Columns grid cannot be toggled via DOM events (mousedown/click/pointerdown all fail)
- "None" link unchecks ALL columns globally, not just the filtered subset — cannot use search + None to exclude a single column
- `DG.InputBase.forElement()` returns null for Train Model form inputs — no JS API access to set Features programmatically
- `DG.Grid.fromRoot()` finds the grid but exposes no dataFrame or cell-level API for the column selector grid

### Suggestions for the platform
- Add DOM-accessible checkboxes (not canvas-rendered) in the Select Columns dialog for automation
- Make "All"/"None" apply only to filtered columns when search is active
- Expose `DG.InputBase` references on Train Model form input host elements

### Suggestions for the scenario
- Note that the Model Engine dropdown defaults to an available engine; explicit selection of "Eda: Linear Regression" was not needed since it was auto-selected
- Consider adding a verification step to check the trained model result (e.g., view type, column count)
