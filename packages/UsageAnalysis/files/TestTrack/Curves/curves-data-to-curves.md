---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
related_bugs:
  - GROK-18886
  - GROK-20679
---

# Data to Curves, Calculate MSR, PZFX files

**Data > Curves > Data to Curves** groups well-level assay data (one row per well: assay, batch,
run, compound, target, concentration, readout) into one fitted curve per group and opens the result
as a new table. **Data > Curves > Calculate MSR...** computes the minimum significant ratio (a Python
script). GraphPad Prism `.pzfx` files open with their XY tables as curve columns.

## Setup

A clean session.

## Scenarios

### 1. Data to Curves (GROK-18886)

1. Open a well-level table, e.g. a CSV with the columns `Assay Name`, `Batch ID`, `Run ID`,
   `Compound ID`, `Target Entity`, `Concentration`, `% Inhibition` and two compounds with six
   concentrations each.
2. **Data > Curves > Data to Curves**.
3. **Verify:** the dialog picks the columns by their names (`Concentration` for **Concentration**,
   `% Inhibition` for **Readout**, `Batch ID`, `Run ID`, `Target Entity`, `Assay Name`).
4. Click **OK**.
5. **Verify:** a new table opens with one row per curve (one per compound here) and a curve column;
   the group names are joined with ` | ` and have no empty parts and no `null` (GROK-18886); the
   console has no errors.

On localhost (Curves 1.13.3, 2026-10-14) **Data > Curves** lists only **Calculate MSR...**:
Data to Curves is missing from the menu although the package declares it (suspected defect). The
automated test of this scenario is deferred.

### 2. Calculate MSR (GROK-20679) — manual

Runs a Python script on the server.

1. Open a table that already has an `assay_name` column.
2. **Data > Curves > Calculate MSR...**, fill the inputs, **OK**.
3. **Verify:** the result is added; no `ValueError: cannot insert assay_name, already exists`.

### 3. A PZFX file (Curves 1.12.0)

1. **Browse > Files > App Data > Curves**, double-click `dose-response-4pl.pzfx`.
2. **Verify:** the view `dose-response-4pl.pzfx` opens with a grid whose curve column draws
   charts; the console has no errors.
3. Close the view, right-click `dose-response-4pl.pzfx` and pick **Open**.
4. **Verify:** the same view opens; the console has no errors.
