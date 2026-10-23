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
