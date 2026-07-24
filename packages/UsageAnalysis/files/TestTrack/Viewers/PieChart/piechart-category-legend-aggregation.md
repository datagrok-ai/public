---
feature: piechart
realizes_atlas:
  - piechart.cp.setup-aggregation-legend-persistence
  - piechart.category-column-drives-legend
realizes:
  - viewers.pie-chart
priority: p0
target_layer: playwright
coverage_type: smoke
realized_as:
  - piechart-category-legend-aggregation-spec.ts
related_bugs:
  - id: GROK-17231
    status: fixed
expected_results:
  - anchor: "Configuration ladder"
    expectation: "each configuration step builds on the previous ones and leaves
      its setting ON: with Category = RACE the legend item labels equal the RACE
      column's category list exactly (set equality both ways); turning Show
      Value on repaints the pie by a per-color canvas pixel delta above 500;
      switching the angle aggregation count → avg(AGE) and then avg(AGE) →
      sum(AGE), with Show Value staying on, each repaint the pie by a per-color
      delta above 10000; after right-clicking the Asian legend item and
      confirming the red swatch, the Asian legend item's DOM color equals
      rgb(214, 39, 40) and the RACE column carries the categorical color tag
      with Asian mapped to #d62728; no new page or console errors and the
      viewer root stays attached throughout"
  - anchor: "Layout round-trip"
    expectation: "after saving the layout of the fully configured pie (Category
      RACE, sum(AGE), Show Value on, custom Asian color, title 'Pie
      Persistence Probe'), adding a Scatter plot, and re-applying the saved
      layout, the view's viewer set equals the SAVED set (Pie chart present,
      the later-added Scatter plot absent) AND the restored pie carries the
      category, angle column, aggregation, Show Value, title, and the red
      Asian legend swatch; the probe layout is deleted afterwards even on
      failure"
  - anchor: "Project save / Close All / reopen"
    expectation: "after saving the view as a project through the ribbon Save
      button, Close All, and reopening the saved project, a Pie chart viewer is
      present AND its category, angle column, aggregation, Show Value, and
      title equal the persisted values AND the Asian legend swatch is still
      rgb(214, 39, 40) with the column color tag intact (a cross-session
      round-trip); the configured settings are NOT reverted afterwards — the
      probe layout and project are deleted, nothing else; even on failure they
      never leak"
---

# Pie Chart — Category, Legend, Aggregation, and Persistence

## Purpose

Verify the pie chart's configuration ladder without clicking any slice: settings
are stacked one on top of another — the category, the value labels, the angle
aggregation, a custom category color picked through the legend's own dialog, and
a title — with each addition confirmed by a real product signal, and none of
them ever switched back off. The fully configured state is then saved twice —
as a layout and as a project — and both round-trips must bring every setting
back.

## Setup

1. Close all open views.
2. Open the demog dataset: `System:DemoFiles/demog.csv` — wait for the table view to load.
3. Add a Pie chart to the current table view via the Toolbox viewer icon.

## Scenarios

### Scenario 1: Configuration ladder ending in layout and project persistence

Settings accumulate: nothing configured in an earlier step is undone by a later
one, and the final state is exactly what gets saved and restored.

Steps:
1. Set Category to RACE and make the legend always visible. Read the legend item labels and the RACE column's own category list from the data — they must contain exactly the same entries, nothing more and nothing less.
2. Turn Show Value on — value labels appear on the slices and the pie visibly repaints. Show Value stays on for the rest of the test.
3. Set Segment Angle Column (Context Panel > Data) to AGE with aggregation count and let the pie finish drawing.
4. Switch the aggregation to avg — the slice geometry and the value labels change: the repaint between settled frames is large.
5. Switch the aggregation to sum — the pie repaints again by a large amount. The aggregation stays at sum for the rest of the test.
6. Right-click the Asian item in the legend — a color dialog titled Asian opens. Pick the red swatch and confirm with OK — the Asian legend item turns red and the RACE column now carries Asian's custom color. The color stays.
7. Set the title to "Pie Persistence Probe".
8. Throughout steps 1-7 the viewer stays attached and no errors appear.
9. Layout round-trip: save the view layout, add a Scatter plot viewer, then re-apply the saved layout. Verify the view's viewer set equals the SAVED set — the Pie chart is present AND the later-added Scatter plot is absent — and the restored pie still carries the category, angle column, aggregation sum, Show Value, title, and the red Asian legend item. Delete the probe layout.
10. Project round-trip: save the view as a project through the ribbon Save button, dismiss the Share dialog, Close All, then reopen the saved project. Verify the Pie chart is present with the same category, angle column, aggregation, Show Value, title, and the red Asian swatch with the custom color still recorded on the column. Delete the probe project.

Expected:
- Each setting stacks on the previous ones and stays on: the legend equals the RACE categories, Show Value repaints the pie, both aggregation switches repaint it again, the legend-dialog color change turns the Asian swatch red and records the color on the column
- Re-applying the saved layout restores the SAVED viewer set (Pie chart present, Scatter plot absent) with the full configuration intact
- Reopening the saved project restores the same full configuration, including the custom color
- The configured settings are never switched back off; only the probe layout and project are deleted — even when a verification fails, they never leak

## Automation notes

- The viewer handle is `grok.shell.tv.viewers.find(v => v.type === 'Pie chart')`;
  the container is `[name="viewer-Pie-chart"]`. Category and legend visibility
  are set via `categoryColumnName` / `legendVisibility` props.
- Step 1: legend labels come from
  `[name="legend"] .d4-legend-item .d4-legend-value` inside the viewer root;
  the reference list is `df.col('RACE').categories` — sorted set equality both
  ways, nothing hard-coded.
- Steps 2, 4, 5 measure repaints as a per-color canvas histogram delta between
  SETTLED frames (a frame is settled when two consecutive snapshots are
  identical). The total non-white pixel count is NOT a usable signal for this
  viewer: slice angles redistribute inside a constant disc and labels paint
  over already-inked slices, so the total stays flat while colors move.
  Thresholds: Show Value toggle > 500 px, aggregation transitions > 10000 px.
  In demog, sum(AGE) per race is nearly proportional to the row count, so the
  count/sum pair barely changes the geometry — the asserted transitions are
  count → avg and avg → sum. All measured deltas are logged with
  `console.log`.
- Step 6 drives the REAL legend color flow with no substitution fallback:
  `contextmenu` on the Asian `.d4-legend-item` opens
  `.d4-dialog[name="dialog-Asian"]`; the swatch is the `.d4-color-bar` with
  background rgb(214, 39, 40); commit via `[name="button-OK"]`. The signals
  are the legend item's inline color style and the RACE column's
  `.color-coding-categorical` tag carrying `"Asian":"#d62728"`. (Hovering the
  legend item also shows a `[name="legend-icon-color-picker"]` icon that opens
  the same dialog — the right-click path is the scripted one.)
- Steps 9-10: the layout is saved and re-applied via the JS API
  (`tv.saveLayout()` / `grok.dapi.layouts.save` / `tv.loadLayout`). The
  project is saved through the real ribbon Save button (`[name="button-Save"]`,
  dialog `dialog-Save-project`) because only the UI Save captures the view
  layout; the "Share" dialog that follows is dismissed via its CANCEL button.
  The probe project name carries a `Date.now()` suffix. The custom color
  survives both round-trips because it lives in the dataframe column tag,
  which the layout leaves untouched and the project serializes. The probe
  layout and project are deleted in `finally` teardowns, so they are removed
  even when an assertion fails; the viewer's configured settings are
  deliberately NOT reverted.
- The single scenario is realized as three consecutive automated steps whose
  names carry the anchors: "Configuration ladder", "Layout round-trip",
  "Project save / Close All / reopen".

---
{
  "order": 5,
  "datasets": ["System:DemoFiles/demog.csv"]
}
