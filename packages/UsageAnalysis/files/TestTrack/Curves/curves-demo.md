---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p2
related_bugs:
  - GROK-7700
  - GROK-12773
  - GROK-13784
---

# Curves demos

**Browse > Apps > Demo > Curves** has two demos: **Curve Fitting** (a grid of curve columns) and
**Assay Curves** (the Assay Curves table with a Multi curve viewer). That both open is covered by
`UsageAnalysis/bdd/features/apps/demos-data-access-compute-curves.feature`.

## Scenarios

### 1. Curve Fitting: clicking cells and the Filters pane (GROK-13784)

1. Open **Browse > Apps > Demo > Curves > Curve Fitting**.
2. Click a curve cell.
3. **Verify:** the Context Panel shows **Options**, **Chart** and **Fit**; the console has no errors.
4. Open the Filter Panel.
5. **Verify:** the console has no errors (GROK-13784).
6. Scroll the grid to the end and back.
7. **Verify (manual):** scrolling is smooth (GROK-12773).

### 2. Assay Curves

1. Open **Browse > Apps > Demo > Curves > Assay Curves**.
2. **Verify:** the view holds the Assay Curves table (43 rows) with its `Fitted Curve` column and a
   MultiCurveViewer; the console has no errors.
3. Select the rows of one assay.
4. **Verify (manual):** the viewer overlays their curves with a legend (see `multi-curve-viewer.md`,
   Scenario 6).
