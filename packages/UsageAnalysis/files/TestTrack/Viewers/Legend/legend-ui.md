---
feature: legend
realizes_atlas: [legend.cp.categorical-color-coded-show, legend.cp.conditional-vs-numerical-coloring]
realizes: [viewers.line-chart, viewers.scatter-plot]
priority: p0
target_layer: manual-only
coverage_type: regression
manual_only_reason: |
  The artifact under test is a pixel/canvas-level visual property (palette
  ordering, subplot layout, gradient swatch rendering) that automation cannot
  assert — human visual verification only.
related_bugs: []
---

# Legend — visual / undecidable steps (manual QA only)

Each section is self-contained and includes the preceding setup steps so a
manual QA can reproduce end-to-end.

## Setup (shared)

Each scenario starts by opening SPGI fresh:

1. Open `System:DemoFiles/chem/SPGI.csv` — wait for table view + semantic-type detection

## Scenarios

### 1. Line chart — distinct palette colors (no two adjacent share)

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Wait for the legend to render — it lists every distinct value of `Series`
   (typically 4–10 categories)
4. **Visual check:** each category in the legend is rendered with a distinct
   color from the categorical palette — no two adjacent categories share the
   same color.

### 2. Line chart — Multi Axis: each Y line gets its own subplot

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']`
4. Enable **Multi Axis** (`lc.props.multiAxis = true`)
5. **Visual check:** each Y line is rendered in its own dedicated subplot —
   a vertical stack of plot regions, not a single merged plot.

### 3. Line chart — corresponding legend block updates on Y-column replacement

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']`
4. Wait for both lines to render with their own per-Y legend blocks
5. Replace one Y column: `lc.props.yColumnNames = ['Average Mass', 'NIBR logP']`
6. **Visual check:** the legend block corresponding to the replaced Y column
   updates to reflect the new column's values; the other Y column's legend
   block stays put.

### 4. Scatter plot — numerical color gradient swatch

1. Add a **Scatter plot**
2. Add a derived numerical column:
   `df.columns.addNewCalculated('testLinear', "if(${Stereo Category}=='S_UNKN', null, ${Average Mass})")`
3. Set the Scatter plot's **Color** to the new numerical column
   (`sp.props.colorColumnName = 'testLinear'`)
4. **Visual check:** the legend renders a linear (numerical) color scale —
   a vertical gradient swatch with min/max axis labels — instead of
   categorical legend items.

---
{
  "order": 99,
  "datasets": ["System:DemoFiles/chem/SPGI.csv"]
}
