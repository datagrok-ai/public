---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p0
related_bugs:
  - GROK-14876
  - GROK-15258
  - GROK-17312
  - GROK-18073
  - GROK-20232
  - public#2748
  - public#3215
  - public#3792
---

# Unusual and broken curve data

A curve column must not break the grid, the dialogs or the rest of the table when some of its
cells cannot be drawn. Curve JSON needs only `series` with `points` to be detected; everything else
(fit function, colors) has a default.

## Setup

A clean session. Open `System:DemoFiles/demog-1000.csv`.

## Scenarios

### 1. A formula that returns JSON (public#3792, GROK-20232)

1. **Edit > Add New Column...**. Formula: `{}`. Name: `jsoncol`. Click **OK**.
2. **Verify:** the dialog closes and the column `jsoncol` is added; the console has no errors.
3. **Edit > Add New Column...** again.
4. **Verify:** the dialog opens (it used to stop opening after a JSON formula).
5. Formula: `'{"series":[{"points":[{"x":0.1,"y":1},{"x":1,"y":5},{"x":10,"y":9}]}]}'`. Name:
   `curvecol`. Click **OK**.
6. **Verify:** `curvecol` is added with the semantic type `fit` and its cells draw charts — the JSON
   has no `fitFunction`, which is not needed (GROK-18073); the console has no errors.

### 2. A broken curve cell does not break the table (GROK-15258, public#2748)

1. Open a table with a curve column in which one cell holds truncated JSON, e.g. a CSV:
   ```
   id,curve
   1,"{""series"":[{""points"":[{""x"":0.1,""y"":1},{""x"":1,""y"":5},{""x"":10,""y"":9}]}]}"
   2,"{""series"":[{""points"":[{""x"":0.1"
   3,"{""series"":[{""points"":[{""x"":0.1,""y"":2},{""x"":1,""y"":6},{""x"":10,""y"":8}]}]}"
   ```
2. **Verify:** the column is recognised as curves (`fit`); the grid shows all 3 rows and the headers;
   rows 1 and 3 draw charts, row 2 shows the "cannot render" mark (red cross, manual) instead of
   throwing.
3. Click the cell of row 2.
4. **Verify:** the console has no errors and no error balloon (a cell with no series no longer
   kills the Context Panel, Curves 1.13.1).

### 3. Export to CSV and open it again (GROK-17312, public#3215, GROK-14876)

1. Open `System:DemoFiles/curves.csv`.
2. Toolbar **Export** > **As CSV**. Close all; open the downloaded file with **Browse > Open local
   file**.
3. **Verify:** the reopened table has the same 30 rows and the same 6 columns in the same order; the
   five curve columns are recognised as `fit` again; every cell of `custom fit function` (the last
   column) still holds its JSON, nothing is cropped.

### 4. Save a project and open it again (GROK-14876)

1. Open `System:DemoFiles/curves.csv`. Save the project (**SAVE** on the ribbon).
2. Close all. Open the project.
3. **Verify:** the curve columns are detected as `fit` and draw charts; the console has no errors.
4. Delete the project.
