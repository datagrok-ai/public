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

## Edge cases

The edge-case round of 26 Oct 2026: data shapes, columns removed or renamed under the open viewer,
layouts naming a missing column, property bounds, closing the viewer. Its checks are now `@edge-cases`
scenarios of their own, after the master scenarios of `scatter-plot-zoom-filter.feature`, and the
known failures in `scatter-plot-known-failures.feature`. Nothing here asks for a
library, binding or code change; it lists what the features could not claim and why.

### Edge cases: Still open

#### 1. Steps requested

1. **`rows shown` counts the rows an empty size hides.** With Size = HEIGHT on demog-1000 the plot
   does not draw the 128 rows whose HEIGHT is empty (help scatter-plot.md, "X, Y, Colors, Sizes,
   Markers"), yet `rows shown` reads 1000 with Show Markers With Empty Size off and on: the status
   counts `visibleMarkers`, which ignores the option (`scatterplot_status.dart` `_markers` against
   `scatterplot_renderer.dart:201`). The ink check is requested from ApiTests below. Needs `rows shown` to
   skip rows the renderer skips (or a `markers drawn` reading); then: `Then scatter plot viewer
   should show 872 rows` with the option off and 1000 with it on.

#### 2. Questions for owners

1. **Per-category lines at exactly 20 categories.** The help says per-category regression and
   moving average lines work "up to 20 categories" (help scatter-plot.md:121, :161, :448, :456);
   the code requires fewer than 20 (`features/regression_line.dart:20`,
   `features/moving_average.dart:13`). At 20 categories the help promises per-category lines and
   the code draws one common line. Which one is right: "fewer than 20" in the help, or `<= 20` in
   the code?

#### Decided

1. A layout naming a missing Color column is applied and keeps the name in Color, and a removed
   Color or Size column is dropped with no other column picked. These checks, with no rows and
   renamed X and Color columns, are requested from ApiTests in `../MISSING.md`.
2. Moving Average Window typed as 1 or 101 in the context panel is not taken (the out-of-bounds
   rule): the scenario outline of `scatter-plot-known-failures.feature` (back into
   `scatter-plot-lines.feature` once fixed), one Examples row
   per value, with the moving average switched on first; a known failure of core#1200 (the value
   is written).

#### 3. API or Dart tests to write

1. Dart regTest: a layout whose X, Color or Label column the table lacks (as
   `viewers_test.dart:1881` for the trellis and `:1728` for the pivot).
2. Dart benchmarks for core#1122: hover miss and selection drag count on 1M rows (timing is not a
   UI claim).
3. ApiTests: Moving Average Window round trip at 1 and 101, once the bound question is answered.
4. ApiTests: with Size = HEIGHT, Show Markers With Empty Size on adds ink to the plot and off takes
   it away; nothing logged.
5. ApiTests: a one-row table plotted on x and y shows 1 row and the marker of row 1; nothing logged.
6. ApiTests: removing AGE, the X column, takes it out of X; nothing logged.

Moved out of BDD (no user gesture); requested from ApiTests.

### Edge cases: Won't do

- Scale (1M-row hover and drag): timing belongs to the Dart benchmarks, not a feature.

## Closed issues the regression round of 26 Oct 2026 could not claim

Written while turning the closed scatter plot issues into scenarios, now in the scatter plot's own
journeys, `scatter-plot-labels-regressions.feature` and `../viewer-chrome.feature`.

- **The missing-column message (public#3091).** A layout applied after its Y column was removed draws
  `y: Column "HEIGHT" does not exist.` on the plot (`DataFrameViewer._validateColumns`), but
  `_ScatterPlotStatus.of` builds its `WidgetStatus` without `error:` (the line chart passes
  `error: plot.error`), so `scatter plot viewer should report the error …` reads `""`. Removing the
  column live is no test either: the plot re-picks another column. Until then
  `scatter plot viewer should report no error` cannot fail, and no scenario uses it. Needs
  `error: p.error` in `scatterplot_status.dart`:

  ```gherkin
  Scenario: A plot whose Y column is gone says which column it misses (public#3091)
    Given user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should report no error
    When user saves the layout of the current table view
    And user removes "HEIGHT" column
    And user loads the saved layout
    Then scatter plot viewer should report the error "y: Column \"HEIGHT\" does not exist."
    And no errors should have been logged
  ```
- **Markers the renderer drew (public#2108, public#1411).** `rows shown` counts the rows the axes can
  place (`isRowDrawable`), not what the renderer drew, so it reads 872 whether the dot renderer puts the
  128 rows with no HEIGHT on the zero line (public#2108) and whether a numeric Color column's markers
  come back after a layout (public#1411, which said markers went missing). A pixel claim does not
  help: a layout round trip repaints the canvas by a few pixels with every marker in place (measured
  26 Oct 2026). Needs a `markers drawn` reading counted by `_ScatterPlotRenderer` in the last frame:

  ```gherkin
  Scenario: Dot markers draw no row whose Y is empty (public#2108)
    When user sets "Marker Type" property of scatter plot viewer to "dot"
    Then the "markers drawn" reading of scatter plot viewer should be 872
    When user sets "Marker Type" property of scatter plot viewer to "circle"
    Then the "markers drawn" reading of scatter plot viewer should be 872
    And no errors should have been logged

  Scenario: A numeric Color column comes back from a layout with every marker it drew (public#1411)
    When user sets "Color" property of scatter plot viewer to "AGE"
    And user saves the layout of the current table view
    And user loads the saved layout
    Then the "markers drawn" reading of scatter plot viewer should be 872
    When user sets "Color" property of scatter plot viewer to ""
    Then no errors should have been logged
  ```
- **A wheel zoom ends after the plot says it has settled.** After `user scrolls the mouse wheel up 3
  times over the "view" area`, `user remembers the value range` and an immediate `should show the
  remembered value range` failed: the window went on narrowing (x 79.3..127.3, then 85.5..121.1), so
  the zoom animation is not in `isRenderPending`. The regression scenarios set a window with the range
  slider instead; a claim on a wheel zoom's end needs the animation in the pending flags.

### A viewer's picture logs three 404s (public#942, public#3052)

General > Save as PNG and Control+C on a viewer make the picture (what the two issues fixed holds),
but each logs `cannot fetch resource: …/font/MaterialIcons-Regular.woff2|woff|ttf, status: 404`:
the picture embeds every `@font-face` of the page, and `core/client/xamgle/web/google_fonts.css:313-318`
still declares Material Icons from `/font/`, whose files left the repository in 2018 (`70272c5e29c`).
A product bug to file in core; once fixed, `viewer-chrome.feature` gets "General > Save as PNG
downloads the viewer as a .png file (public#942)" and "Control+C on a viewer copies its picture and
says so (public#3052)" (`user watches downloads` / `a file matching "\.png$" should have been
downloaded`; `the clipboard should hold a PNG image of at least 2000 bytes`).
