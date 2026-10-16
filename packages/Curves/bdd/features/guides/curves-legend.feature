@guide
Feature: The legend of a dose-response chart names the curves without hiding them
  A guide: the answer to "how does the new curve legend work?" (Curves 1.13.1+, GROK-20112), filmed for
  a community post. A Multi curve viewer overlays the curves of the selected rows; its Legend property
  (Context Panel > Data > Legend, a column picker) names each curve by a column. The legend takes the
  corner of the plot the curves leave freest - here the bottom left, since the curves run high on the
  left and fall only on the right - shortens a name too long for that corner (hovering the name gives the
  whole of it) and folds the rows past the corner into "+N more" (hovering it lists the curves, each in
  its own color).
  Show Legend (Context Panel > Misc > Show Legend) takes the legend off the plot and puts it back.
  Demo: the Assay Curves table of the Curves package, its 30 Collapse Assay curves, named by
  Compound|Assay|Target: more than a full-screen plot has rows for. A click on an AUC cell makes a plain
  number the current cell, so the molecule of the first row does not take the Context Panel over from
  the viewer's properties.

  Scenario: Name the curves, read a shortened name and the folded rows, then hide and show the legend
    Given user is logged in
    And simple mode is off
    And user opens assay-curves dataset
    And the grid of the current view does not show the "Molecule" column
    And user selects rows where "Assay Name" is "Collapse Assay"
    And user adds a MultiCurveViewer viewer of the selected rows taking 45 percent of the view width
    When user clicks on a "AUC" cell of grid other than the current one
    And user picks "Properties..." from the viewer menu of MultiCurveViewer viewer
    And user selects "Compound|Assay|Target" in "Legend" property in context panel
    Then "Legend" property of MultiCurveViewer viewer should be "Compound|Assay|Target"
    And the curves legend of MultiCurveViewer viewer should be in the bottom left corner
    When user hovers over a shortened name in the curves legend of MultiCurveViewer viewer
    Then the tooltip should show the whole of the shortened name
    When user hovers over the "+N more" row of the curves legend of MultiCurveViewer viewer
    Then the tooltip should list the curves of the legend, each in its own color
    When user moves the pointer off the curves legend of MultiCurveViewer viewer
    And user expands "Misc" category in context panel
    And user unchecks "Show Legend" property in context panel
    Then the curves legend of MultiCurveViewer viewer should not be drawn
    When user checks "Show Legend" property in context panel
    Then the curves legend of MultiCurveViewer viewer should be drawn
