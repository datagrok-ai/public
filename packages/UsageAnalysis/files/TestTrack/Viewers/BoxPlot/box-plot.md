---
feature: boxplot
realizes_atlas:
  - boxplot.cp.property-surface-smoke
  - boxplot.int.inside-outside-values
  - boxplot.int.auto-layout-hides-chrome
  - boxplot.int.showmarkers-gates-marker-props
realizes: [viewers.box-plot]
priority: p2
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-18677
    status: fixed
  - id: GROK-19297
    status: fixed
  - id: GROK-18361
    status: fixed
  - id: GROK-20272
    status: fixed
realized_as:
  - box-plot-spec.ts
expected_results:
  - anchor: "Context menus as property paths"
    expectation: >-
      Clicking Show Inside Values and Show Outside Values in the main context
      menu's Misc group flips the matching viewer setting (prop read-back after
      the menu click — the menu is the actuation path, the prop is the signal)
      AND changes the drawn point population (settle-gated canvas diff per
      toggle); re-clicking both items restores the settings and the points
      (round-trip; boxplot.int.inside-outside-values — the menu drives the same
      two flags as the property panel). With a marker size column set, the
      Markers menu's Size item renders grayed. No exclusive Y-axis / X-axis
      region menus are asserted — the live build consolidates their items into
      the main menu.
  - anchor: "Statistics and group-comparison menu regions"
    expectation: >-
      With the viewer sized so the statistics table fits, right-clicking the
      statistics region opens the stats menu; while group comparison is off, the
      Group Comparison group's gated items render at opacity 0.5 (grayed, not
      disabled or hidden) and hovering one shows a hint tooltip. The p-value
      region menu offers no Statistics Format item. With group comparison
      enabled, the comparison-strip region menu contains an "Add ... Table" item
      and the Group Comparison items are no longer grayed — item presence only;
      the exported table's semantics belong to boxplot-group-comparison.md.
  - anchor: "Resize and auto layout"
    expectation: >-
      With Auto Layout on, shrinking the viewer to a small size removes the
      column selectors, axes, and statistics table from the viewer chrome (the
      selectors are DOM elements — presence/absence is asserted), and restoring
      the size brings them back (boxplot.int.auto-layout-hides-chrome — auto
      layout overrides the individual show toggles at small sizes). A very
      narrow resize with a legend-bearing coloring applied completes with a zero
      console-error delta (GROK-18677 — no error on resizing the viewer).
  - anchor: "Marker gate and size scaling"
    expectation: >-
      Disabling Show Markers removes the individual data points — only the box
      glyphs remain (settle-gated canvas diff) — AND the entire Marker settings
      group in the property panel becomes disabled (every Marker prop dependsOn
      showMarkers); re-enabling restores both. With a marker size column set,
      switching Size Scaling between linear and log repaints the markers
      (settle-gated canvas diff per switch; GROK-20272).
  - anchor: "Whisker and control-band style"
    expectation: >-
      Whisker Line Width and Whisker Width Ratio changes each yield their own
      settle-gated canvas diff. Setting Control Band Color completes without
      errors (no-error floor — the band is drawn only while control comparisons
      are active). The Whisker Color uniform-vs-sequential transition is owned
      by boxplot-render-stats-color.md and is not re-asserted here.
  - anchor: "Controls visibility"
    expectation: >-
      The size selector is absent by default (the default showSizeSelector is
      false — the chrome baseline has no size selector); enabling Show Size
      Selector adds it. Disabling Show Value Selector and Show Color Selector
      collapses each selector's DOM host as toggled; the category selector and
      the value / category axes are canvas-drawn chrome, so disabling and
      re-enabling each of those yields its own settle-gated canvas repaint.
      Re-enabling restores the DOM selectors, and disabling Show Size Selector
      returns the default baseline (round-trip).
  - anchor: "Title and description"
    expectation: >-
      After enabling Show Title and setting a title, the title text appears on
      the viewer; the description text appears while Description Visibility Mode
      is Always, moves when Description Position is set to Bottom, and
      disappears when the mode is Never.
  - anchor: "Axis font"
    expectation: >-
      Changing Axis Font repaints the axis labels with a zero console-error
      delta — specifically no Infinity.floor error (GROK-19297); restoring the
      font also completes without errors.
  - anchor: "Date category mapping"
    expectation: >-
      With Category 1 set to the datetime STARTED column, setting Category 1 Map
      to Month restructures the category axis into month bins (settle-gated
      re-layout repaint), Quarter restructures it again, and returning Category
      1 to a categorical column restores plain categories — all with no console
      errors.
  - anchor: "Custom tooltip"
    expectation: >-
      With the viewer's Row Tooltip set to AGE, SEX, WEIGHT and the tooltip mode
      switched to show custom tooltip (the column list is inert while the mode
      inherits from the table), hovering a marker shows a tooltip (.d4-tooltip)
      whose content lists exactly those three columns and no others; resetting
      the mode and the Row Tooltip to inherit from the table restores the
      default tooltip columns. Tooltip show/hide modes beyond this activation
      are owned by the Tooltips section and are not asserted here.
  - anchor: "Table switching resets Category 2"
    expectation: >-
      After switching the Table property from demog to spgi-100 with a two-level
      category configuration in place, the viewer renders the new table without
      errors and Category 2 resets to a consistent state — it is either cleared
      or reads a real spgi-100 column; no stale demog column name remains in the
      Category 2 selector (GROK-18361 — broken Category2 reset after changing
      the data table).
  - anchor: "Legend minimum under coloring"
    expectation: >-
      With a legend-bearing categorical coloring applied, the render stays
      valid: the named data canvas keeps non-zero ink (countCanvasPixels > 0),
      the shell-warnings delta is zero, and markerColorColumnName stays applied.
      The box plot's coloring legend is canvas-drawn — no legend DOM host exists
      on this viewer, so legend-DOM asserts are not performed; a vanished legend
      would surface as a coloring drop or a render fault. Legend depth belongs
      to the Viewers/Legend section.
  - anchor: "Double-click resets the view"
    expectation: >-
      After the value-axis range slider narrows the visible value range,
      double-clicking an empty area of the plot fires the viewer's
      d4-boxplot-reset-view event (the gesture path to resetView) AND returns
      the viewport to the full value range. The single-click-vs-double-click
      distinction is owned by boxplot-pointer-select.md and is not re-asserted
      here.
---

# Box plot — property surface smoke

## Purpose

Verifies the Box plot's remaining property surface: context menus exercised as
paths to properties (main-menu groups plus the statistics, p-value, and
group-comparison hit regions), auto-layout chrome hiding and resize stability,
the Show Markers gate and marker size scaling, whisker and control-band style,
controls visibility, title and description, axis font, date category mapping,
the custom row tooltip, table switching, the legend-minimum render floor, and
the double-click view reset.

Selection, persistence, filtering, statistics, and group comparison have their
own dedicated tests in this folder and are not repeated here.

All sections start with the following sequence:
1. Close all open tables and viewers.
2. Open the **demog** dataset.
3. Add a Box plot by clicking its icon in the Toolbox > Viewers section.

## Context menus as property paths

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. Right-click an empty area of the plot — the viewer's context menu opens
   with the Markers, Style, Misc, and Color groups.
3. In the Misc group, click **Show Inside Values** — the points inside the
   whiskers disappear from the plot, and the **Show Inside Values** setting in
   the Context Panel now reads off.
4. Open the menu again and click Misc > **Show Outside Values** — the outlier
   points disappear, and the setting reads off.
5. Re-enable both through the same menu items — the points return and both
   settings read on again (round-trip).
6. In the Context Panel > Markers, set **Marker Size Column** to WEIGHT.
   Right-click the plot and open the Markers group — the **Size** item is
   grayed while a marker size column is set.
7. Set **Marker Size Column** back to empty.

Expected:
- Each Misc menu toggle flips the matching setting AND changes the drawn point
  population — the menu is the actuation path, the setting plus the repaint is
  the signal (round-trip, Steps 3–5).
- With a marker size column set, the Markers menu's Size item renders grayed
  (Step 6).

## Statistics and group-comparison menu regions

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX. Resize the viewer large enough that the statistics
   table is fully visible at the bottom.
2. Right-click the statistics table region — the statistics menu opens. Locate
   the **Group Comparison** group: its gated items are shown dimmed at half
   opacity (not hidden, not disabled).
3. Hover one of the dimmed Group Comparison items — a hint tooltip appears
   explaining the item needs Group Comparison enabled.
4. In the Context Panel > Statistics, make sure **Show P Value** is on.
   Right-click the p-value text region — its menu contains no
   **Statistics Format** item.
5. In the Context Panel > Data, enable **Show Group Comparison**. Right-click
   the comparison strip (the region with the per-group statistics rows) — the
   menu now offers an **Add ... Table** item, and the Group Comparison items
   are no longer dimmed.
6. Disable **Show Group Comparison**.

Expected:
- While group comparison is off, the Group Comparison menu items render dimmed
  at half opacity with a hint tooltip on hover (Steps 2–3).
- The p-value region menu offers no Statistics Format item (Step 4).
- With group comparison on, the comparison-strip region menu contains an
  "Add ... Table" item — item presence only; the exported table itself belongs
  to the group-comparison scenario (Step 5).

## Resize and auto layout

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX. Confirm **Auto Layout** is on in the Context Panel >
   Style.
2. Note the visible chrome: the column selectors, the value and category axes,
   and the statistics table.
3. Drag the viewer's edge to shrink it to a small size — the column selectors,
   axes, and statistics table disappear from the viewer.
4. Restore a large size — the chrome returns.
5. In the Context Panel > Color, set **Marker Color Column** to SEX (a
   legend-bearing coloring). Resize the viewer to a very narrow width, then
   restore it — no errors appear during the whole resize sequence
   (GROK-18677).
6. Set **Marker Color Column** back to empty.

Expected:
- Auto layout hides the column selectors, axes, and statistics table at small
  sizes and restores them at a large size, overriding the individual show
  toggles (Steps 3–4).
- The very narrow resize with a coloring applied raises no errors
  (GROK-18677, Step 5).

## Marker gate and size scaling

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. In the Context Panel > Markers, disable **Show Markers** — the individual
   data points disappear from the plot (only the box glyphs remain), and the
   entire Markers settings group in the Context Panel becomes inactive.
3. Re-enable **Show Markers** — the points return and the Markers group is
   active again.
4. Set **Marker Size Column** to WEIGHT.
5. Switch **Size Scaling** from linear to log — the marker sizes repaint.
6. Switch **Size Scaling** back to linear and set **Marker Size Column** back
   to empty.

Expected:
- Disabling Show Markers removes the data points (only boxes remain) AND
  disables the whole Markers settings group; re-enabling restores both
  (Steps 2–3).
- Switching Size Scaling between linear and log repaints the markers
  (Step 5).

## Whisker and control-band style

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX. The plot style is box (the default).
2. Set **Whisker Line Width** to 4 — the boxes repaint. Set
   **Whisker Width Ratio** to 0.3 — the boxes repaint again.
3. Set **Control Band Color** to a different color — the change completes
   without errors.

Expected:
- Whisker Line Width and Whisker Width Ratio each produce their own repaint
  (Step 2).
- Setting Control Band Color completes without errors — the band is drawn only
  while control comparisons are active, so no visible change is required
  (Step 3).

## Controls visibility

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. Observe the default chrome: the value, category, and color selectors and
   both axes are visible; the size selector is NOT shown by default.
3. In the Context Panel > Size, enable **Show Size Selector** — the size
   selector appears.
4. In the Context Panel, disable **Show Value Selector**, then
   **Show Category Selector**, then **Show Color Selector** — each selector
   disappears as it is toggled.
5. Disable **Show Value Axis** and **Show Category Axis** — the axes
   disappear.
6. Re-enable the five toggles disabled in Steps 4–5, and disable
   **Show Size Selector** — the default chrome baseline is restored.

Expected:
- The size selector is absent in the default chrome (Step 2).
- Each visibility toggle adds or removes its own chrome element (Steps 3–5),
  and the round trip restores the default baseline (Step 6).

## Title and description

1. In the Context Panel > Description, enable **Show Title** and set **Title**
   to "Age by Race" — the title text appears on the viewer.
2. Set **Description** to "Box plot of patient ages" and
   **Description Visibility Mode** to Always — the description text appears.
3. Change **Description Position** to Bottom — the description moves to the
   bottom.
4. Set **Description Visibility Mode** to Never — the description disappears.

Expected:
- Title and description text appear, move, and disappear on the viewer as
  driven (Steps 1–4).

## Axis font

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. In the Context Panel > Style, change **Axis Font** to a larger size — the
   axis labels repaint and no error appears (GROK-19297).
3. Restore the original **Axis Font** — again with no error.

Expected:
- Both Axis Font changes repaint the axis labels with no errors — in
  particular no Infinity.floor error (GROK-19297, Steps 2–3).

## Date category mapping

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to STARTED (a datetime column).
2. Set **Category 1 Map** to Month — the category axis restructures into month
   bins.
3. Set **Category 1 Map** to Quarter — the axis restructures again into
   quarter bins.
4. Set **Category 1** to RACE — plain categories return.

Expected:
- Each mapping change restructures the category axis (Month, then Quarter),
  and returning to a categorical column restores plain categories, all with no
  errors (Steps 2–4).

## Custom tooltip

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to RACE. In Markers, set **Marker Size** to 10 so individual
   points are easy to hover.
2. In the Context Panel > Tooltip, set the **Row Tooltip** column list to AGE,
   SEX, WEIGHT and switch **Show Tooltip** to "show custom tooltip" — the
   column list takes effect only in this mode.
3. Hover a marker point and hold — the tooltip shows exactly the AGE, SEX, and
   WEIGHT rows and no other columns.
4. Reset **Show Tooltip** to "inherit from table" and clear the
   **Row Tooltip** list.
5. Hover a marker point again — the default tooltip columns are shown.

Expected:
- With the custom Row Tooltip set and the custom mode active, the hover
  tooltip lists exactly the three configured columns (Step 3); resetting
  restores the default tooltip (Step 5). Tooltip show/hide modes beyond this
  activation are covered in the Tooltips section of TestTrack.

## Table switching resets Category 2

1. Open the spgi-100 dataset alongside demog (Files > App Data >
   Chem/tests/spgi-100.csv).
2. On the demog Box plot, set **Value** to AGE, **Category 1** to SEX, and
   **Category 2** to RACE (a two-level category axis is in place).
3. In the Context Panel > Data, switch **Table** to spgi-100.
4. Verify the viewer renders the new table without errors and **Category 2**
   resets to a consistent state — it is either cleared or reads a real
   spgi-100 column; no stale demog column name remains in the Category 2
   selector (GROK-18361).
5. Set **Value** to Average Mass and **Category 1** to Series — the plot
   renders the spgi-100 data.

Expected:
- After the table switch, Category 2 is reset consistently (cleared or a valid
  spgi-100 column, never a stale demog name) and the viewer renders the new
  table without errors (GROK-18361, Steps 3–5).

## Legend minimum under coloring

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. In the Context Panel > Color, set **Marker Color Column** to RACE — a
   legend-bearing categorical coloring is applied.
3. Verify the render stays valid: the plot keeps drawing its data, no new
   warnings appear, and the coloring column stays applied. The coloring legend
   on this viewer is drawn on the canvas itself — there are no separate legend
   elements to inspect.
4. Set **Marker Color Column** back to empty.

Expected:
- With the coloring applied, the plot still draws its data, the warning count
  does not grow, and the coloring column remains applied (Step 3) — a vanished
  legend would surface as a coloring drop or a render fault. Legend behavior in
  depth is covered in the Viewers > Legend section of TestTrack.

## Double-click resets the view

1. In the Context Panel > Value, set **Value** to AGE; in Data, set
   **Category 1** to SEX.
2. Zoom in with the range slider so the visible value range is narrowed.
3. Double-click an empty area of the plot — the view resets (the viewer fires
   its reset-view event).

Expected:
- The double-click on empty space fires the reset-view event and the view
  returns to the full range (Step 3). The single-click-vs-double-click
  distinction is covered by the pointer-selection test.

## Automation notes

- target_layer rationale: context-menu navigation, property-panel chrome
  presence checks, settle-gated canvas diffs, viewer resize drags, and
  tooltip DOM reads all require a live browser session; Playwright is the
  appropriate layer.
- Menu-path signal discipline: the menu click is the actuation; the signal is
  the flipped prop read back through the props channel PLUS the structural
  repaint. Never a set-then-read of a prop that was just written via props.
- Context-menu items carry name= attributes (name="div-<Group>---<Item>"
  plus a d4-name raw-label attribute) — resolve by name= with a label-text
  fallback. The popup menu DOM builds asynchronously after the contextmenu
  event: poll for the target item (or for the item count to settle) before
  reading or clicking — a synchronous query sees a partial menu.
- Region-menu surface: no exclusive Y-axis / X-axis region menus exist —
  their items are consolidated into the main menu. Author menu paths against
  the main menu and the statistics / p-value / group-comparison hit regions
  only.
- Group Comparison gate state: the gated menu items are grayed via an
  opacity-0.5 style, NOT the disabled attribute — assert the computed opacity
  plus the hint tooltip on hover.
- Statistics-region precondition: the statistics table auto-hides entirely
  when its values do not fit the available width — assert a minimum viewer
  width before any right-click on the statistics or p-value regions.
- Fresh geometry: re-read the viewer root's bounding rect immediately before
  every hover/click group that follows a viewer add, close, or resize — dock
  relayout shifts the root and cached coordinates go stale.
- Polling over fixed pauses: after each menu click or property change, poll
  for the settled condition instead of fixed timeouts.
- Auto-layout chrome channel: the column selectors are DOM elements — assert
  their presence/absence across the resize; the axes and statistics table hide
  together with them at small sizes.
- Controls-visibility channels: the value / color / size selectors hide via
  their div-column-combobox-* hosts (display:none); the category selector and
  the value / category axes are canvas-drawn — the category1 combobox host
  does NOT hide with showCategorySelector — so those three toggles are
  asserted as settle-gated canvas repaints in both directions.
- Axis font (GROK-19297): the failure signature is an Infinity.floor error —
  watch the error collector for that message specifically, plus a zero
  overall error delta across both font changes.
- Marker gate: the Marker group's disabled state is a property-panel DOM
  signal (the rows render dimmed while showMarkers is off) plus the
  settle-gated canvas diff as the points disappear and return.
- controlBandColor: no control band is drawn while control comparisons are
  off, so the honest signal for the set-path smoke is a no-error floor.
- Custom tooltip channel: the Row Tooltip column list is inert while Show
  Tooltip inherits from the table — switch the mode to "show custom tooltip"
  before the hover. The tooltip renders a .d4-row-tooltip-table with one row
  per column and the column name in the first cell — assert the set of those
  first-cell labels (order- and value-independent); poll for the tooltip after
  the hover rather than using fixed waits.
- Legend floor channel: countCanvasPixels (TestTrack/helpers/viewers.ts) on
  the named data canvas must stay > 0 with the coloring applied; the box plot
  renders no legend DOM (no [name="legend"] host, no .d4-legend-item elements
  under marker coloring) — do not author legend-DOM asserts.
- Double-click reset: zoom first with a trusted drag of the value-axis range
  slider — bp.viewport (a Rect getter) narrows under the drag and returns to
  the full range on reset, so BOTH the d4-boxplot-reset-view event (listener
  registered before the gesture) and the viewport round-trip are the signals.
  Plain marker/label clicks do not populate the viewport — only the slider
  zoom makes the viewport channel observable.
- Derived from atlas cp boxplot.cp.property-surface-smoke:
  core/client/d4/lib/src/viewers/box_plot/features/context_menu_feature.dart#L22
- Interactions realized in this scenario:
  boxplot.int.inside-outside-values (Context menus as property paths,
  Steps 3–5); boxplot.int.auto-layout-hides-chrome (Resize and auto layout,
  Steps 3–4); boxplot.int.showmarkers-gates-marker-props (Marker gate and
  size scaling).

---
{
  "order": 8,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
