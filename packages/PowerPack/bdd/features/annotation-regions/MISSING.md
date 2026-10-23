# Annotation regions and axis annotations — what the features cannot say yet

Written during the Grid / Filter panel / Legend gap round of 23 Oct 2026. Added in that round:
`annotation-regions-scatter-lasso.feature` (the scatter plot's own lasso becomes an area region) and
`../formula-lines/formula-lines-axis-annotations-box-plot.feature` (the box plot's value axis has
the Annotations group, its category axis does not).

## Still open

- **A formula line drawn on the box plot.** The box plot reports `region <n>` areas for annotation
  regions but nothing for formula lines (`box_plot/CLAUDE.md`, and a probe with a line set shows no
  `formula line` area and no `formula lines` reading), so Add Line and Add Band on its value axis are
  claimed by what they wrote into the look, not by what was drawn. Needs the scatter plot's
  convention on the box plot — `formula line "<title>"` / `formula band "<title>"` areas and a
  `formula lines` reading of what the last pass drew:

  ```gherkin
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of box plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then box plot viewer should have a "formula line AGE = 45.0" area
    And the "formula lines" reading of box plot viewer should be 1
  ```
- **A formula line's label.** The viewers that report a line report the band the line sweeps, not
  the box its title is drawn in, so "the label is shown at the line's end" or "the label moved with
  Label Position" cannot be claimed. Needs a `formula line "<title>" label` area next to the line's
  own, in each viewer's status (`scatterplot_status.dart`, the line chart's and the histogram's).
- **A lasso region keeps the free-form shape it was drawn with.** `annotation-regions-scatter-lasso.feature`
  claims an area region keyed to WEIGHT and HEIGHT (`"type":"area"` in the look and a `region 1`
  area), but `region <n>` is the polygon's bounding rectangle (`regionRect` in
  `annotation_regions_mixin.dart`), so a lasso stored or drawn as that rectangle would pass too. The
  stored points are the look's, an API test's subject; what a person sees is the outline. Needs a
  `region <n> corner` area — a small box inside the bounding rectangle but outside the drawn
  polygon, reported by `addStatus` for area regions — or a `region <n> vertices` reading of the
  polygon the last pass drew. Then:

  ```gherkin
    When user drags a lasso over the "view" area of scatter plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then scatter plot viewer should have a "region 1" area
    And the "region 1 vertices" reading of scatter plot viewer should be at least 8
    And the "region 1 corner" area of scatter plot viewer should be painted in no color
  ```
- **A box plot region spans Q1..Q3 on the value axis.** `formula-lines-axis-annotations-box-plot.feature`
  claims the region Add Region draws is as wide as the plot and shorter than it, not that its top
  and bottom sit at AGE 56 and 36. The status reports the region's rectangle in pixels and no
  value-axis scale to map them back. Needs `region <n> min` / `region <n> max` readings (the region's
  Y bounds in data units, as drawn) from `addStatus` in `annotation_regions_mixin.dart`, or a
  `y = <value>` point area on the box plot's value axis. Then:

  ```gherkin
    When user picks "Annotations > Add Region" from the context menu of the "bottom edge of y axis" area of box plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then the "region 1 min" reading of box plot viewer should be 36
    And the "region 1 max" reading of box plot viewer should be 56
    And no errors should have been logged
  ```
