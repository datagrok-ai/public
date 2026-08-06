---
feature: boxplot
realizes_atlas:
  - boxplot.cp.render-stats-color-sync
  - boxplot.int.violin-needs-bins-and-style
  - boxplot.int.pvalue-toggle-key-and-menu
realizes:
  - viewers.box-plot
realized_as:
  - boxplot-render-stats-color-spec.ts
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: github-2966
    status: fixed
  - id: github-3066
    status: fixed
  - id: GROK-18245
    status: fixed
  - id: GROK-17506
    status: fixed
expected_results:
  - anchor: Scenario 1 Step 2
    expectation: >-
      With whiskerColor set to null (the default, no explicit color), the box
      glyphs are colored per-category with sequential categorical coloring (not
      a single uniform color across all boxes). Multiple distinct hue regions
      are present on the data canvas corresponding to the distinct SEX
      categories.
  - anchor: Scenario 1 Step 3
    expectation: >-
      After sizing the Box Plot viewer to a sufficient width, the statistics
      table is visible inside the viewer canvas (the _statsValuesFit gate
      passes). The default statistics entries (mean, median, Q1, Q3, count) are
      present in the statistics strip.
  - anchor: Scenario 1 Step 6
    expectation: >-
      After enabling Show Total Count, Show Inliers Count, Show Outliers Count,
      Show Stdev, Show Q1, and Show Q3 in Context Panel > Statistics, each
      corresponding prop reads true and the statistics table gains additional
      rows. The canvas statistics strip re-renders with no console error.
  - anchor: Scenario 1 Step 8
    expectation: >-
      After setting Statistics Format to an explicit numeric format string
      (e.g. '#,##0.00') in Context Panel > Statistics, the statisticsFormat
      prop reads the
      new value and the statistics table re-renders without a console error. The
      p-value displayed alongside the statistics is UNCHANGED — Statistics
      Format applies only to the statistics table (dependsOn showStatistics) and
      does not affect the p-value display.
  - anchor: Scenario 1 Step 9
    expectation: >-
      After setting an explicit whiskerColor value (e.g. a named color or hex),
      all box glyphs repaint to a single uniform color — the settle-gated canvas
      diff confirms the per-category sequential coloring transitions to one
      uniform color.
  - anchor: Scenario 2 Step 3
    expectation: >-
      With Show Group Comparison confirmed off, the bare p-value label is
      visible on the Box Plot canvas (the t-test overlay is active for 2
      categories). Hovering over the bare p-value label reveals the
      group-comparison reveal icon (name="show-group-stats") in the DOM — no
      test-name tooltip appears on bare p hover.
  - anchor: Scenario 2 Step 5
    expectation: >-
      After pressing the T key on the keyboard, the Show P Value property in
      Context Panel > Statistics flips from true to false (or from false to
      true) and the p-value overlay on the canvas appears or disappears
      accordingly (T-key toggle transition signal;
      boxplot.int.pvalue-toggle-key-and-menu).
  - anchor: Scenario 2 Step 7
    expectation: >-
      After switching Category 1 to RACE (3+ categories), the p-value area on
      the canvas re-renders: the Alexander-Govern test branch is active. The
      bare p-value label is present (canvas floor; hover no longer produces a
      test-name tooltip — it reveals the reveal icon instead).
  - anchor: Scenario 2 Step 9
    expectation: >-
      After right-clicking the statistics region and toggling a statistic from
      the context menu (e.g. toggling Show Total Count), the corresponding
      showTotalCount prop reads the new value — the context-menu toggle path
      drives the same property as the Context Panel checkbox
      (boxplot.int.pvalue-toggle-key-and-menu: third toggle path).
  - anchor: Scenario 3 Step 3
    expectation: >-
      After setting Plot Style to violin in Context Panel > Style, the Box Plot
      canvas re-renders with violin shapes replacing the box glyphs — a large
      settle-gated canvas diff is observed on the data canvas. All category
      value distributions are present: no category is silently dropped from the
      violin rendering (github-2966 regression guard).
  - anchor: Scenario 3 Step 5
    expectation: >-
      After changing Bins from 50 to 500, the violin render changes — a
      settle-gated canvas diff is observed showing different KDE smoothing.
      After changing Interquartile Line Width to 10, another canvas diff is
      observed on the data canvas. Each property change independently produces a
      measurable canvas diff (GROK-18245 — violin props must actually affect the
      render; boxplot.int.violin-needs-bins-and-style).
  - anchor: Scenario 3 Step 7
    expectation: >-
      After changing Violin Line Width (1-4) and Violin Whisker Color, each
      change produces a distinct settle-gated canvas diff on the data canvas.
  - anchor: Scenario 4 Step 4
    expectation: >-
      After setting Marker Color Column to WEIGHT on the Box Plot and then
      changing the linear color scheme Min and Max on the WEIGHT column in the
      grid color settings, the Box Plot canvas re-renders (a settle-gated diff
      is observed) — the shared column tag ensures the box plot follows the
      grid's color-scale change (GROK-17506 regression guard).
  - anchor: Scenario 4 Step 6
    expectation: >-
      After switching the WEIGHT column coloring type from linear to conditional
      in the grid, the Box Plot canvas re-renders (a settle-gated diff is
      observed showing the new conditional palette) and the selection marker
      pixels remain visible above the conditional-color baseline (echoes
      github-3066 — selection highlight must survive non-default coloring).
  - anchor: Scenario 4 Step 8
    expectation: >-
      After switching Marker Color Column to the categorical SEX column, the Box
      Plot categorical palette on the marker circles matches the palette shown
      in the grid column header color indicator for SEX. Changing one category's
      color in the grid's column-coloring dialog updates the corresponding
      marker/legend color in the Box Plot (categorical palette sync).
  - anchor: Scenario 5 Step 3
    expectation: >-
      After setting Value to STARTED (a datetime column) and resizing the viewer
      to sufficient width, the statistics strip re-renders without a console
      error — the datetime Statistics Format fix ensures the datetime value
      column renders statistics without errors. No new page or browser errors
      appear.
---

# Box Plot — Rendering, Statistics, and Grid Color Synchronization

## Setup

1. Close all open tables and viewers.
2. Open the demog dataset: navigate to Files > App Data and open
   System:DemoFiles/demog.csv.
3. Add a Box Plot viewer to the table view via the toolbar (Add viewer
   > Box Plot).
4. In Context Panel > Data, set Value to AGE and set Category 1 to SEX.

## Scenarios

### Scenario 1: Box coloring baseline and statistics ladder

Steps:
1. Verify that Show Group Comparison is off in Context Panel > Data (confirm
   the prop is false or the toggle is unchecked). This is required for the
   t-test overlay and statistics steps below — the t-test overlay is
   suppressed while group comparison is on.
2. Inspect the Box Plot canvas: confirm that the two SEX category boxes are
   colored in distinct hues (sequential categorical coloring, one hue per
   category) rather than a single uniform color — the default look while no
   explicit Whisker Color is set.
3. Resize the Box Plot viewer (drag its edge) so it is wide and tall enough
   that the statistics strip is fully visible. Verify that the statistics
   table area is visible in the viewer (the _statsValuesFit gate passes —
   a too-small viewer hides the table entirely).
4. In Context Panel > Statistics, confirm that Show Statistics is enabled
   (on by default). Observe the statistics strip in the viewer.
5. In Context Panel > Statistics, enable Show Mean (if not already on).
6. In Context Panel > Statistics, enable Show Total Count, Show Inliers Count,
   Show Outliers Count, Show Stdev, Show Q1, and Show Q3 one at a time.
7. Verify that each statistics flag now reads true in Context Panel and the
   statistics strip in the viewer canvas has gained additional entries with no
   console error.
8. In Context Panel > Statistics, set Statistics Format to an explicit format
   string (e.g., '#,##0.00' or a date-compatible format such as
   'yyyy-MM-dd'). Confirm the prop reads the new format value.
9. Verify that the statistics table re-renders without a console error and
   that the p-value displayed (if visible below the statistics strip) is
   UNCHANGED — Statistics Format applies only to the statistics table and
   does not affect the p-value display.
10. In Context Panel > Style, set Whisker Color to an explicit color (e.g.,
    select any named color from the color picker that is visibly different
    from the sequential categorical hues currently applied).
11. Verify that all box glyphs repaint to a single uniform color — a
    settle-gated canvas diff confirms the per-category sequential coloring
    transitions to one uniform color.

Expected:
- With Whisker Color unset (the default): box glyphs show per-category sequential
  distinct hues (two distinct color regions for SEX female/male on the
  canvas).
- After enabling statistics flags, each prop reads true and the statistics
  strip in the viewer gains the corresponding entries with no console error.
- After setting Statistics Format, the statistics table re-renders without
  error; the p-value is unchanged.
- After setting an explicit Whisker Color, all boxes repaint to a single
  uniform color (canvas diff observed; sequential palette disappears).

### Scenario 2: P-value toggle — three equivalent paths

Steps:
1. In Context Panel > Statistics, enable Show P Value (the property panel
   checkbox path). Confirm that the p-value label appears on the Box Plot
   canvas with Category 1=SEX (2 categories, Welch's t-test branch).
2. Confirm Show Group Comparison is off (the t-test overlay shows only while
   group comparison is off).
3. Move the mouse pointer over the bare p-value label on the canvas (hover
   over the p-value text). Verify that the group-comparison reveal icon
   (the element with name="show-group-stats") appears in the DOM. Confirm
   that NO test-name tooltip appears on bare p hover (the hover reveals the
   icon, not a tooltip with the test conclusion).
4. In Context Panel > Statistics, disable Show P Value using the checkbox.
   Verify that the p-value overlay disappears from the canvas.
5. Click the Box Plot canvas to give it keyboard focus. Press the T key.
   Verify that the Show P Value property in Context Panel > Statistics flips
   back to true and the p-value overlay reappears on the canvas (the T-key
   toggle path).
6. Press the T key again. Verify that Show P Value flips to false and the
   p-value overlay disappears (T off).
7. In Context Panel > Data, change Category 1 to RACE (3 or more distinct
   categories). Re-enable Show P Value via the property panel checkbox.
   Verify that the p-value area re-renders on the canvas (the
   Alexander-Govern test branch is now active for 3+ categories; canvas
   floor — bare p is present).
8. Move the mouse pointer over the bare p label. Confirm again that the
   reveal icon (name="show-group-stats") appears rather than a test-name
   tooltip (the hover behavior is the same for 2 and 3+ categories).
9. Right-click the statistics region of the Box Plot (click on the statistics
   strip area). In the context menu that appears, toggle one statistic entry
   (e.g., toggle Show Total Count on or off). Verify that the corresponding
   Show Total Count setting reads the newly toggled value — the context-menu path
   drives the same property as the Context Panel checkbox.

Expected:
- Hovering the bare p-value label reveals the name="show-group-stats" icon
  in the DOM; no test-name tooltip appears on bare p hover.
- Pressing the T key toggles Show P Value on→off and off→on (both
  directions verified).
- With Category 1 set to RACE (3+ categories), the Alexander-Govern test
  branch is active and the bare p is present on the canvas (floor).
- Right-clicking the statistics region and toggling a statistic from the
  context menu changes the corresponding prop value (third toggle path).

### Scenario 3: Violin rendering — per-category completeness and prop effects

Steps:
1. In Context Panel > Data, set Category 1 back to SEX. In Context Panel >
   Style, set Plot Style to violin.
2. Verify that a large canvas diff is observed on the data canvas — the box
   glyphs are replaced by violin shapes (settle-gated diff confirming the
   plot style switch).
3. Verify that both SEX category distributions are present in the violin
   render — no category is silently dropped from the violin rendering
   (github-2966 regression guard: check that the data canvas has two distinct
   violin shapes, one per SEX category, both with non-zero ink).
4. In Context Panel > Style, set Bins to 50. Observe the canvas.
5. In Context Panel > Style, change Bins to 500. Verify that a settle-gated
   canvas diff is observed compared to Bins=50 — KDE smoothing changes
   measurably (GROK-18245).
   In Context Panel > Style, set Interquartile Line Width to 10. Verify that
   a canvas diff is observed on the data canvas (Interquartile Line Width
   affects the violin render — GROK-18245).
6. In Context Panel > Style, set Plot Style back to box (to confirm the
   round-trip before testing additional violin props).
   In Context Panel > Style, set Plot Style to violin again.
7. In Context Panel > Style, set Violin Line Width to 4 (the maximum; range
   1-4). Verify a settle-gated canvas diff is observed.
   In Context Panel > Style, set Violin Whisker Color to an explicit color
   (pick a color clearly different from the current rendering). Verify a
   settle-gated canvas diff is observed.

Expected:
- After switching Plot Style to violin, a large canvas diff is observed and
  both SEX violin shapes are present with non-zero ink (github-2966).
- Changing Bins from 50 to 500 produces a measurable canvas diff (different
  KDE smoothing; GROK-18245). Changing Interquartile Line Width to 10
  produces its own canvas diff (GROK-18245).
- Changing Violin Line Width and Violin Whisker Color each produce their own
  settle-gated canvas diff.

### Scenario 4: Grid color synchronization

Steps:
1. In Context Panel > Style, set Plot Style back to box.
2. In Context Panel > Data, set Value back to AGE and Category 1 to SEX.
3. In the grid (the main table view), right-click the WEIGHT column header
   and open the column coloring dialog (Color > Color by > Linear or use
   the column's color icon). Apply a linear color scheme to WEIGHT.
4. In Context Panel > Color, set Marker Color Column to WEIGHT. Observe the
   Box Plot markers — they are now colored by WEIGHT with the linear scheme.
5. In the grid, change the linear color scheme's Min and Max values for the
   WEIGHT column (edit the Color Min and Color Max inputs in the column
   coloring panel). Verify that the Box Plot canvas re-renders in response
   (a settle-gated canvas diff is observed on the data canvas — GROK-17506:
   the box plot must respond to linear color scheme min/max changes in the
   grid).
6. In the grid, switch the WEIGHT column coloring type from linear to
   conditional coloring (in the column coloring dialog, pick Conditional).
   Verify that the Box Plot canvas re-renders (settle-gated diff — the
   markers adopt the conditional palette). With a non-empty selection active
   (shift-drag a band over some markers in the Box Plot first), verify that
   the selection marker pixels remain visible above the conditional-color
   baseline (echoes github-3066: selection highlight survives non-default
   coloring).
7. In Context Panel > Color, change Marker Color Column to SEX (a categorical
   column). Observe that the Box Plot markers adopt the SEX categorical
   palette.
8. In the grid, change one SEX category's color (right-click the SEX column
   header > Color > pick a different color for one category). Verify that the
   corresponding marker color in the Box Plot updates to match the new
   category color (categorical palette synchronization between grid and
   box plot markers/legend).

Expected:
- After setting Marker Color Column to WEIGHT and changing the linear
  color scheme Min/Max in the grid, the Box Plot canvas re-renders
  (GROK-17506 regression guard).
- After switching WEIGHT column coloring to conditional in the grid, the
  Box Plot canvas re-renders and selection pixels remain visible
  (github-3066 guard).
- After switching Marker Color Column to SEX and changing a category's
  color in the grid, the Box Plot marker color for that category updates
  to match.

### Scenario 5: Datetime value — statistics without errors

Steps:
1. In Context Panel > Data, set Value to STARTED (a datetime column in
   demog).
2. Confirm Show Statistics is enabled. Resize the Box Plot viewer if needed
   so the statistics strip is visible.
3. Open the browser developer console (or observe the page's error indicator)
   before the next step.
4. Observe the statistics strip in the viewer. Verify that the statistics
   render without a new console error — the datetime Statistics Format fix
   ensures that datetime value columns render statistics without errors.
5. In Context Panel > Statistics, set Statistics Format to a datetime-compatible
   format string (e.g. 'yyyy-MM-dd HH:mm'). Verify that the statistics re-render
   without a console error and that the Statistics Format setting reads the
   new format value.

Expected:
- With Value set to STARTED (datetime), the statistics strip renders without
  a new browser console error (datetime Statistics Format fix).
- After setting a datetime-compatible Statistics Format, the strip re-renders
  without error and Statistics Format reads the set value.

## Automation notes

- target_layer rationale: canvas diff capture, DOM event probing
  (name="show-group-stats"), keyboard T-key delivery, right-click context
  menu navigation, and grid column coloring all require a live browser
  session; Playwright is the appropriate layer.
- Statistics table precondition (Scenarios 1 and 5): the statistics table
  auto-hides entirely when its values do not fit (_statsValuesFit) — assert
  a minimum viewer width before any statistics assertion.
- Box coloring channels (Scenario 1): whiskerColor=null means per-category
  sequential coloring — assert two distinct dominant hues on the data
  canvas, never a single uniform hue; an explicit whiskerColor is a
  settle-gated canvas diff to one uniform color.
- P-value hover (Scenario 2): the bare-p hover signal is the DOM appearance
  of [name="show-group-stats"] — park the pointer on the p-value text area
  (~40 % of the canvas height from the top) and poll. No test-name tooltip
  exists on bare-p hover; the conclusion text lives in the group-comparison
  strip, owned by cp.group-comparison-ladder.
- T-key delivery (Scenario 2 Step 5): focus the canvas before sending the
  keypress; the signal is bp.props.showPValue before and after.
- Violin data completeness (Scenario 3 Step 3): count non-white pixels in
  two vertical strips of the data canvas, one per category — both must be
  > 0 (github-2966).
- Violin prop diffs (Scenario 3): Bins and Interquartile Line Width only
  produce visible changes while Plot Style is violin
  (boxplot.int.violin-needs-bins-and-style); each prop change is its own
  settle-gated canvas diff (GROK-18245).
- Grid color sync (Scenario 4): the shared channel is the column's color
  tag (column_coloring_mixin getVisibleColorRange) — after a grid-side
  change, wait for a settled Box Plot canvas before capturing the diff.
- Selection under conditional coloring (Scenario 4 Step 6): select rows
  before the coloring-type switch; the selection-hue pixel count
  (countSelectionHuePixels against selectedRowsColor) must stay > 0 after.
- Datetime step (Scenario 5): grok.shell.warnings has no JS-visible
  accessor — the console/page-error delta is the no-error channel.
- Derived from atlas cp boxplot.cp.render-stats-color-sync:
  core/client/d4/lib/src/viewers/box_plot/box_plot_core.dart#L1387
- Interactions realized in this scenario:
  boxplot.int.violin-needs-bins-and-style (Scenario 3 Steps 4-5)
  boxplot.int.pvalue-toggle-key-and-menu (Scenario 2 Steps 4-9)
