---
feature: curves
target_layer: playwright
pyramid_layer: ui-smoke
coverage_type: regression
priority: p1
related_bugs:
  - GROK-18138
  - GROK-18166
  - GROK-18171
  - GROK-18172
  - GROK-18389
  - GROK-18467
  - GROK-19801
  - GROK-20668
  - core#646
---

# The Multi curve viewer inside a Trellis plot

A Trellis plot takes the MultiCurveViewer as its inner viewer (**Curves** in the viewer selector,
top left of the trellis). Every cell overlays the curves of its rows. The trellis Context Panel then
has two tabs, **Trellis** and **MultiCurveViewer**; the second holds the inner viewer's properties,
and a change there reaches every cell.

## Setup

A clean session. Open `System:DemoFiles/curves.csv`.

## Scenarios

### 1. Curves as the inner viewer (GROK-18138)

1. Add a **Trellis plot**. Set its X split to `smiles`.
2. In the viewer selector, pick **Curves**.
3. **Verify:** the inner viewer is MultiCurveViewer; every cell draws; the trellis shows no error;
   the console has no errors.

### 2. The full-screen view of a cell (GROK-18172)

1. Hover a trellis cell; click the full-screen icon in its corner.
2. **Verify:** a dialog shows the cell's curves (not `No data to show.`).
3. Click the dialog's **X**.
4. **Verify:** the dialog closes on the first click; the console has no errors.

### 3. The MultiCurveViewer tab of the Context Panel (GROK-19801, GROK-18166, GROK-20668)

1. Click the trellis **Gear** icon; open the **MultiCurveViewer** tab.
2. **Verify:** the tab lists the inner viewer's properties (Data, Axes, Misc …), not an empty grid
   (GROK-19801).
3. Expand **Misc**; check **Log Y**.
4. **Verify:** the trellis inner viewer has Log Y on (every cell redraws on a log Y axis, GROK-18166).
   Also repeat this scenario right after the scenarios of `multi-curve-viewer.md` in the same
   session: the cells redraw on a log Y axis, and a layout saved and applied again keeps Log Y
   ([core#854](https://github.com/datagrok-ai/core/issues/854)).
5. In **Curves**, click `...`, then **All**, **OK**.
6. **Verify:** **Curves** reads `5 / 5` and the Context Panel stays on the **MultiCurveViewer** tab
   (GROK-20668: it used to jump to the Trellis tab).

### 4. Curves on a table without curves (GROK-18171, core#646)

1. Open `System:DemoFiles/demog-1000.csv`. Add a Trellis plot split by `SEX`.
2. Pick **Curves** in the viewer selector. Move the mouse over the cells.
3. **Verify:** the cells say `The MultiCurveViewer viewer requires a minimum of 1 curves column.`;
   the console has no errors.

core#646 (open): the `NullError: method not found: 'ab' on null` comes back when no curve cell has
been drawn in the session before — repeat the scenario right after a page reload.

### 5. Switching inner viewers on a table without curves (GROK-18389)

1. On `demog-1000`, in the trellis viewer selector pick **Scatter plot**, **Bar chart**, **Box plot**,
   then **Curves**.
2. **Verify:** each switch works; the console has no errors.

### 6. Scrolling a trellis of curves (GROK-18467)

1. On `curves.csv`, Trellis split by `smiles` with **Curves**, switch on **Tiles**.
2. Scroll the trellis back and forth; toggle **Log X** in the MultiCurveViewer tab.
3. **Verify (manual):** cells already seen redraw at once; a toggle takes one recalculation.
