# Scatter plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `scatterplot/CLAUDE.md` "Automation
surface", the TestTrack `scatterplot-*.md` and `scatter-plot-ui.md`. The "(no value)" legend click is
in `scatter-plot-legend-no-value.feature`; the rejected line formula (GROK-17282) in
PowerPack's `formula-lines-dialog.feature`.

## Needs a reading

- **Structures on the X axis** (`scatter-plot-ui.md`, "Structures as markers and as an axis", items
  7-8, spgi-100). With a molecule column on X the axis draws structures instead of tick labels;
  the status reports the axis column and range, not what the ticks are drawn as. Needs an
  `x axis labels drawn as` reading (`text` / `structure`, from the axis painter's renderer), claimed
  after the pick and again after a saved layout is applied.

## Edge cases (`scatter-plot-edge-cases.feature`)

Written for `scatter-plot-edge-cases.feature` (26 Oct 2026), the guard suite of edge cases: data
shapes, columns removed or renamed under the open viewer, layouts naming a missing column, property
bounds, closing the viewer. Nothing here asks for a library, binding or code change made by that
feature; it lists what the feature could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

1. **`rows shown` counts the rows an empty size hides.** With Size = HEIGHT on demog-1000 the plot
   does not draw the 128 rows whose HEIGHT is empty (help scatter-plot.md, "X, Y, Colors, Sizes,
   Markers"), yet `rows shown` reads 1000 with Show Markers With Empty Size off and on: the status
   counts `visibleMarkers`, which ignores the option (`scatterplot_status.dart` `_markers` against
   `scatterplot_renderer.dart:201`). The feature claims the option by ink. Needs `rows shown` to
   skip rows the renderer skips (or a `markers drawn` reading); then: `Then scatter plot viewer
   should show 872 rows` with the option off and 1000 with it on.

#### 2. Questions for owners

1. **Per-category lines at exactly 20 categories.** The help says per-category regression and
   moving average lines work "up to 20 categories" (help scatter-plot.md:121, :161, :448, :456);
   the code requires fewer than 20 (`features/regression_line.dart:20`,
   `features/moving_average.dart:13`). At 20 categories the help promises per-category lines and
   the code draws one common line. Which one is right: "fewer than 20" in the help, or `<= 20` in
   the code?

#### Decided (`viewer-rules-edge-cases.feature`)

1. A layout naming a missing Color column leaves Color empty, and a removed Color or Size column is
   dropped with no other column picked: rows of the layout and the removal outlines. No rows, and
   renamed X and Color columns, are rows there too.
2. Moving Average Window typed as 1 or 101 in the context panel is not taken (the out-of-bounds
   rule): written in `scatter-plot-edge-cases.feature` with the moving average switched on first, a
   known failure of core#1200 (the value is written).

#### 3. API or Dart tests to write

1. Dart regTest: a layout whose X, Color or Label column the table lacks (as
   `viewers_test.dart:1881` for the trellis and `:1728` for the pivot).
2. Dart benchmarks for core#1122: hover miss and selection drag count on 1M rows (timing is not a
   UI claim).
3. ApiTests: Moving Average Window round trip at 1 and 101, once the bound question is answered.

### Edge cases: Won't do

- Scale (1M-row hover and drag): timing belongs to the Dart benchmarks, not a feature.
