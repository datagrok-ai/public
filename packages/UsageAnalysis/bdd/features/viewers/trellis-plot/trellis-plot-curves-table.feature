@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot on another table, with the Multi Curve viewer inside
  A trellis set to another table works on that table — the table it reports — and a Multi Curve
  viewer (from the Curves package; "Curves" in the viewer selector, MultiCurveViewer as the type the
  trellis reports) fills its cells there; set back to the first table,
  it draws its canonical grid again. Translated from the "Multi Curve inner viewer (and table
  switching)" section of TestTrack Viewers/TrellisPlot/trellis-plot.md, steps 1-4, 8 and 9, on
  demog-1000 and curves.

  Not translated: steps 5-7 (setting the curves' X and Y, paging with +/- and moving the zoom slider
  inside Multi Curve cells), which the case itself leaves manual for want of a recon of the curve
  viewer's own controls; and the curves table's row count read from the trellis (it reports the table
  it is bound to, not its rows — MISSING.md).

  Background:
    Given user is logged in
    And the "Curves" package is installed
    And user opens curves dataset
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then the "cells" reading of trellis plot viewer should be 8
    And trellis plot viewer should be bound to table "demog-1000"

  Scenario: Set to the curves table, the trellis works on it and takes the Multi Curve viewer
    When user sets "Table" property of trellis plot viewer to "curves"
    Then trellis plot viewer should be bound to table "curves"
    When user picks "Curves" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    And the "cells drawn" reading of trellis plot viewer should be 1
    And the "blank cells" reading of trellis plot viewer should be 0
    When user clicks on settings icon of trellis plot viewer
    Then context panel should be visible
    And "Table" property in context panel should be visible
    And trellis plot viewer should report no error
    And no errors should have been logged

  Scenario: Set back to demog-1000, the canonical grid comes back
    When user sets "Table" property of trellis plot viewer to "demog-1000"
    And user sets properties of trellis plot viewer:
      | X Column Names | SEX          |
      | Y Column Names | RACE         |
      | Viewer Type    | Scatter plot |
    Then trellis plot viewer should be bound to table "demog-1000"
    And the "cells" reading of trellis plot viewer should be 8
    And the "rows shown" reading of trellis plot viewer should be 1000
    And no errors should have been logged
