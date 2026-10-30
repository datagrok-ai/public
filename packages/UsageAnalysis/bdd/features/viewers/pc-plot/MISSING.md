# PC plot — what the features cannot say yet

Written during the viewer gap round of 23 Oct 2026 against `pc_plot/CLAUDE.md` "Automation surface"
and `pc-plot-ui.md`. What could be said is in `pc-plot-hover-and-layout.feature`.

## Needs a reading

- **Show Mouse Over Row Group** ("Mouse-over and row group", item 2). With the property on, hovering a
  row that belongs to a row group (a category hovered in another viewer, `dataFrame.mouseOverRowGroup`)
  draws every line of the group in the mouse-over colour on the overlay. `lines drawn` counts only
  `_renderAllLines`; nothing counts the group lines. `pc-plot-hover-and-layout.feature` claims the
  group pass as a repaint (and none with the flag off); the exact count needs a `mouse-over group lines drawn` reading
  (the rows the overlay's group pass stroked, `pc_plot_core.dart` where `showMouseOverRowGroup` is
  read):

  ```gherkin
  Scenario: Hovering a category in another viewer draws its whole group on the plot
    Given user adds a bar chart viewer with:
      | Split | RACE |
    When user hovers over the "bar Asian" area of bar chart viewer
    Then the "mouse-over group lines drawn" reading of pc plot viewer should be <Asian rows>
    When user sets "Show Mouse Over Row Group" property of pc plot viewer to "false"
    Then the "mouse-over group lines drawn" reading of pc plot viewer should be 0
    And no errors should have been logged
  ```

## Needs a step

- **Undock and dock again** ("Viewer basics", items 1-2): no step floats a viewer out of the dock
  manager (`user undocks {widget}` and `{widget} should be floating`, through the panel's title-bar
  menu, in `bindings/tiers/viewers/widgets.ts`). Docking along an edge is in `viewer-chrome.feature`.
- **General > Save as PNG / Save to Gallery** ("Context menu: General"): both leave the page — a
  download and a server entity. Needs `the download should be a PNG named {string}` (Playwright's
  `download` event) and a gallery entry swept at feature end; `viewer-chrome.feature` claims only that
  the items are offered.

## Not a gap

- **A grid column header dropped on the plot.** `pc-plot-ui.md` lists the grid-to-plot drag among the
  manual checks, but the plot takes no foreign drop: its overlay accepts a drag only when
  `AppEvents.currentlyDragging.dragSource == this` (`pc_plot_core.dart:125`), i.e. its own axis
  labels. Dropping `header STARTED` onto `view` leaves `axes` at 3 (measured). If adding an axis that
  way is wanted, it is a feature request, not a missing step.

## Edge cases (`pc-plot-edge-cases.feature`)

Written for `pc-plot-edge-cases.feature`: the edge cases of this viewer that the feature could not state with the current step library and the expected behaviour documented today. No library, binding or code change was made for them.

### Still open

**Decided**: a Log Columns axis over a column holding zero or negative values says it cannot be logarithmic (core#1203, a defect): written in `pc-plot-edge-cases.feature` as a known failure on a non-empty `error` reading (today the axis is drawn linear without a message, `pc_plot_core.dart`, `col.min > 0`; the wording of the message is not set).

**Decided**: a transformation turned on releases the filter of the axis sliders it replaces (core#967, behaviour): claimed in `pc-plot-edge-cases.feature`, green on localhost (core#967 is closed).

**Decided** (`viewer-rules-edge-cases.feature`): an axis column removed takes its axis away and the others stay; a layout naming a missing column is applied without that axis (a row of the layout outline); a renamed axis column is followed (a row of the rename outline); every row filtered out is a row of the no-rows outline.

**API or Dart tests to write**

1. Dart regTest: a constant column as an axis under Normalized with density box or violin (normalization divides by max − min).

### Won't do

- `bins`, `horzMargin` and the line-width sliders past their bounds: property-grid clamping, the same for every viewer.
