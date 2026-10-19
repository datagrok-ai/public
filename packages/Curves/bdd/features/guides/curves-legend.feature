@guide
Feature: The legend of a dose-response chart names the curves without hiding them
  A guide: the answer to "how does the new curve legend work?" (Curves 1.13.1+, GROK-20112), filmed for
  a community post. A Multi curve viewer overlays the curves of the selected rows; its Legend property
  (Context Panel > Data > Legend, a column picker) names each curve by a column. The legend takes the
  corner of the plot the curves leave freest - here the bottom left, since the curves run high on the
  left and fall only on the right - shortens a name too long for that corner (hovering the name gives the
  whole of it) and folds the rows past the corner into "+N more" (hovering it lists the curves).
  Show Legend (Context Panel > Misc > Show Legend) takes the legend off the plot and puts it back.
  Demo: the Assay Curves table of the Curves package, its 30 Collapse Assay curves, named by
  Compound|Assay|Target ("compound_id 8 | Collapse Assay | target_10"): more than a full-screen plot has
  rows for, the last of them compound_id 20 and 21. A click on an AUC cell makes a plain number the
  current cell, so the molecule of the first row does not take the Context Panel over from the viewer's
  properties.

  The legend is read through the viewer's own status: a "legend row N" area per row it drew ("shortened
  legend row N" for a name it cut), "legend more", and the "legend corner" and "legend shown" readings.
  Which rows are shortened depends on the width of the plot, so the hover takes the first shortened one.
  Film it with BDD_GUIDE_VIEWPORT=1280x720: in the 800x500 help window the plot is below the legend's minimum.

  Scenario: Name the curves, read a shortened name and the folded rows, then hide and show the legend
    Given user is logged in
    And simple mode is off
    And user opens assay-curves dataset
    And user removes "Molecule" column
    And user selects rows where "Assay Name" is "Collapse Assay"
    And user adds a MultiCurveViewer viewer of the selected rows taking 45 percent of the view width
    When user clicks on a "AUC" cell of grid other than the current one
    And user picks "Properties..." from the viewer menu of MultiCurveViewer viewer
    And user selects "Compound|Assay|Target" in "Legend" property in context panel
    Then the "legend corner" reading of MultiCurveViewer viewer should be "bottom left"
    When user hovers over the first "shortened legend row" area of MultiCurveViewer viewer
    Then tooltip should contain text "| Collapse Assay | target_10"
    When user hovers over the "legend more" area of MultiCurveViewer viewer
    Then tooltip should contain text "compound_id 20 | Collapse Assay | target_10"
    When user moves the pointer away from MultiCurveViewer viewer
    And user expands "Misc" category in context panel
    Then "Show Legend" property in context panel should be checked
    When user unchecks "Show Legend" property in context panel
    Then the "legend shown" reading of MultiCurveViewer viewer should be "false"
    When user checks "Show Legend" property in context panel
    Then the "legend shown" reading of MultiCurveViewer viewer should be "true"
