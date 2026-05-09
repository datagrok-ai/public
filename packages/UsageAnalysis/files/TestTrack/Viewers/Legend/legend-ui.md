---
feature: legend
sub_features_covered:
  - legend.color-scale.numerical
target_layer: ui-only
coverage_type: regression
priority: p2
pyramid_layer: manual
produced_from: split
original_paths:
  - public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/line-chart.md
  - public/packages/UsageAnalysis/files/TestTrack/Viewers/Legend/scatterplot.md
date_split: 2026-05-08
related_bugs: []
---

# Legend — visual / undecidable steps (manual QA only)

Scenario steps that require human visual verification — they cannot be
deterministically asserted via JS API or Playwright UI because the artifact
under test is a pixel/canvas-level visual property (palette ordering,
subplot SVG layout, gradient swatch rendering) that isn't exposed through
the DOM on `dev.datagrok.ai` (2026-05-08).

Each section is self-contained and includes the preceding setup steps so a
manual QA can reproduce end-to-end.

`target_layer: ui-only` — no `.ts` spec is generated for this file
(per E-LAYER-COMPLIANCE-01: ui-only scenarios must not have a `.ts` body).

## Setup (shared)

Each scenario starts by opening SPGI fresh:

1. Open `System:DemoFiles/SPGI.csv` — wait for table view + semantic-type detection

## Scenarios

### 1. Line chart — distinct palette colors (no two adjacent share)

**Source:** moved from `line-chart.md` Scenario 1 step 3 (2026-05-08).

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Wait for the legend to render — it lists every distinct value of `Series`
   (typically 4–10 categories)
4. **Visual check:** each category in the legend is rendered with a distinct
   color from the categorical palette — no two adjacent categories share the
   same color.

**Why manual:** categorical palette ordering is internal; pixel-level color
comparison across legend items is not deterministic at the JS API surface.

### 2. Line chart — Multi Axis: each Y line gets its own subplot

**Source:** moved from `line-chart.md` Scenario 2 step 2 (2026-05-08).

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']`
4. Enable **Multi Axis** (`lc.props.multiAxis = true`)
5. **Visual check:** each Y line is rendered in its own dedicated subplot —
   a vertical stack of plot regions, not a single merged plot.

**Why manual:** subplot SVG layout requires DOM measurement; the Line chart
exposes only a single `[name="legend"]` element regardless of subplot
count, so subplot-count cannot be derived from the DOM tree at the JS API
surface.

### 3. Line chart — corresponding legend block updates on Y-column replacement

**Source:** moved from `line-chart.md` Scenario 4 step 4 (2026-05-08).

1. Add a **Line chart**
2. Set **Split** = `Series`
3. Configure two Y columns: `lc.props.yColumnNames = ['Average Mass', 'TPSA']`
4. Wait for both lines to render with their own per-Y legend blocks
5. Replace one Y column: `lc.props.yColumnNames = ['Average Mass', 'NIBR logP']`
6. **Visual check:** the legend block corresponding to the replaced Y column
   updates to reflect the new column's values; the other Y column's legend
   block stays put.

**Why manual:** per-line legend blocks are not individually labeled in the
DOM (single `[name="legend"]` element for the whole chart), so per-Y
legend update cannot be observed by a DOM query.

### 4. Scatter plot — numerical color gradient swatch

**Source:** moved from `scatterplot.md` Scenario 1 step 10 first sub-bullet
(2026-05-08).

1. Add a **Scatter plot**
2. Add a derived numerical column:
   `df.columns.addNewCalculated('testLinear', "if(${Stereo Category}=='S_UNKN', null, ${Average Mass})")`
3. Set the Scatter plot's **Color** to the new numerical column
   (`sp.props.colorColumnName = 'testLinear'`)
4. **Visual check:** the legend renders a linear (numerical) color scale —
   a vertical gradient swatch with min/max axis labels — instead of
   categorical legend items.

**Why manual:** when a numerical Color column is bound, the Scatter plot's
`[name="legend"]` element on `dev.datagrok.ai` (2026-05-08) is empty in
the DOM — the gradient swatch is rendered to a canvas without DOM
children. There is nothing to query at the JS API level to verify the
gradient swatch appearance.

## Notes

- `target_layer: ui-only` — no `.ts` spec is generated for this file.
- `pyramid_layer: manual` — visual verification by human QA.
- These scenarios are NOT redundant with their source `.md` files: the
  source files retain JS-API / Playwright proxy assertions (legend item
  count, prop round-trip, etc.); this file holds the visual-only
  assertions that cannot be automated.
- When the platform exposes legend rendering details via DOM (e.g.
  per-subplot `[name="legend"]` elements, per-bin gradient swatch divs),
  the corresponding scenario can move back to its source `.md`.

---
{
  "order": 99,
  "datasets": ["System:DemoFiles/SPGI.csv"]
}
