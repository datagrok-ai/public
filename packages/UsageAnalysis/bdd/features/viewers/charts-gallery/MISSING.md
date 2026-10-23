# Charts gallery (Surface plot, Multiplot) — what the features cannot say yet

Written during the Charts gap round of 23 Oct 2026, against `public/packages/Charts/src/viewers/`
and the library's rules (`public/libraries/bdd/CLAUDE.md`).

## Still open

- **Surface plot has no status.** `SurfacePlot` (`Charts/src/viewers/surface-plot/surface-plot.ts`,
  an `EChartViewer` drawing with echarts-gl / WebGL) implements neither `getWidgetStatus` nor
  `isRenderPending` / `onRendered`. In `charts-gallery.feature` the Projection and Wireframe claims are
  the property read back — an echo: a viewer that ignored both would pass. Pixels of a WebGL canvas are
  not a claim. Needs, in a `surface-plot-status.ts` next to the viewer: readings `x column`,
  `y column`, `z column` (the columns of the series the chart holds, not the properties), `points`
  (the surface's data length after `plotFilter`), `projection` (`option.grid3D.viewControl.projection`
  of the applied option) and `wireframe` (`series[0].wireframe.show`); `isRenderPending` while
  `chart.setOption` has not finished its frame (echarts `finished` event). Then:

  ```gherkin
  Scenario: Projection and Wireframe reach the drawn surface, and the filter narrows it
    Given user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Surface plot" card in "Add Viewer" dialog
    Then the "points" reading of surface plot viewer should be 1000
    When user clicks on grid
    And user clicks on settings icon of surface plot viewer
    Given "Misc" category in context panel is expanded
    When user selects "orthographic" in "Projection" property in context panel
    Then the "projection" reading of surface plot viewer should be "orthographic"
    When user unchecks "Wireframe" property in context panel
    Then the "wireframe" reading of surface plot viewer should be "false"
    When user opens the filter panel
    And user clicks on the "category F of SEX" area of filter panel
    Then the "points" reading of surface plot viewer should be 553
    And no errors should have been logged
  ```

- **Multiplot has no feature and no status.** `MultiPlotViewer` (`Charts/src/viewers/multiplot/`)
  is a `DG.JsViewer` hosting inner viewers with no `getWidgetStatus`; nothing reports which inner
  plots it laid out. Needs a status listing the inner viewers (`plots`: their types joined, and per
  plot the bound columns) and exposing them as inner widgets, so the existing inner-viewer readings
  (`grid in group analysis viewer` style) reach them. Not written until the status exists.
- **Group Analysis** has no status of its own either; `charts-gallery.feature` claims it through its
  inner grid's readings, which is sufficient for what it shows.
