---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
related_bugs:
  - GROK-20679
---

# Calculate MSR, PZFX files

**Data > Curves > Calculate MSR...** computes the minimum significant ratio (a Python script).
GraphPad Prism `.pzfx` files open with their XY tables as curve columns.

## Setup

A clean session.

## Scenarios

### 1. Calculate MSR (GROK-20679) — manual

Runs a Python script on the server.

1. Open a table that already has an `assay_name` column.
2. **Data > Curves > Calculate MSR...**, fill the inputs, **OK**.
3. **Verify:** the result is added; no `ValueError: cannot insert assay_name, already exists`.

### 2. A PZFX file (Curves 1.12.0)

1. **Browse > Files > App Data > Curves**, double-click `dose-response-4pl.pzfx`.
2. **Verify:** the view `dose-response-4pl.pzfx` opens with a grid whose curve column draws
   charts; the console has no errors.
3. Close the view, right-click `dose-response-4pl.pzfx` and pick **Open**.
4. **Verify:** the same view opens; the console has no errors.
