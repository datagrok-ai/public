@viewers @realizes:viewers.trellis-plot
Feature: Trellis plot through a layout and a project
  The persistence tail of `trellis-plot-split-and-pick-inner.md`: two trellises with different
  settings saved in one layout come back each with its own, and a change to one leaves the other
  alone (GROK-15494); a trellis on Row Source Selected goes through a project with a selection made
  by a cell click — the selection does not come back, the trellis draws nothing until rows are
  selected and draws again when they are — and once more with nothing selected at all, the saved
  state the reopen crash lived on (GROK-19902). Each scenario starts on a fresh demog-1000 view with
  a trellis split by SEX and CONTROL over RACE with pie charts inside: 16 cells; F, false, Caucasian
  holds 475 rows. The second trellis is added bare and set up through its own properties as "second
  trellis plot viewer": the set-up table of "adds a viewer with" goes to the first viewer of the
  type on the page. The layout and the projects are saved through the server's API (the md uses the
  ribbon's Save), so the lost selection is what that round trip keeps, and the trellis "drawing" is
  read as the rows it shows and the cells it lays out, not as pixels. The layout and the projects are
  deleted at the end.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX, CONTROL |
      | Y Column Names | RACE         |
      | Viewer Type    | Pie chart    |
    Then the "cells" reading of trellis plot viewer should be 16

  Scenario: Two trellises keep their own settings through one layout
    Given user adds a trellis plot viewer
    When user sets properties of second trellis plot viewer:
      | X Column Names | RACE      |
      | Y Column Names | SEX       |
      | Viewer Type    | Bar chart |
    Then the open tableview should have 2 trellis plot viewers
    And properties of first trellis plot viewer should be:
      | X Column Names | SEX, CONTROL |
      | Viewer Type    | Pie chart    |
    And properties of second trellis plot viewer should be:
      | X Column Names | RACE      |
      | Viewer Type    | Bar chart |
    When user saves the layout of the current table view to the server
    And user adds a scatter plot viewer
    Then scatter plot viewer should be visible
    When user loads the saved layout
    Then scatter plot viewer should be absent
    And the open tableview should have 2 trellis plot viewers
    And properties of first trellis plot viewer should be:
      | X Column Names | SEX, CONTROL |
      | Y Column Names | RACE         |
      | Viewer Type    | Pie chart    |
    And properties of second trellis plot viewer should be:
      | X Column Names | RACE      |
      | Y Column Names | SEX       |
      | Viewer Type    | Bar chart |
    When user sets "Viewer Type" property of first trellis plot viewer to "Histogram"
    Then the "inner viewer type" reading of first trellis plot viewer should be "Histogram"
    And the "inner viewer type" reading of second trellis plot viewer should be "Bar chart"
    When user sets "Viewer Type" property of first trellis plot viewer to "Pie chart"
    And user clicks on close icon of second trellis plot viewer
    Then the open tableview should have 1 trellis plot viewer
    And properties of trellis plot viewer should be:
      | X Column Names | SEX, CONTROL |
      | Y Column Names | RACE         |
      | Viewer Type    | Pie chart    |
    And the "cells" reading of trellis plot viewer should be 16
    And no errors should have been logged

  Scenario: A selection behind Row Source Selected does not come back with the project, and the trellis still draws
    When user sets "On Click" property of trellis plot viewer to "Select"
    And user clicks on the "cell F, false | Caucasian" area of trellis plot viewer
    Then 475 rows should be selected
    When user sets "Row Source" property of trellis plot viewer to "Selected"
    Then trellis plot viewer should show 475 rows
    When user saves the current view as project "bdd-trellis-selected"
    And user closes all views
    And user opens the "bdd-trellis-selected" project
    Then the open tableview should have 1 trellis plot viewer
    And properties of trellis plot viewer should be:
      | X Column Names | SEX, CONTROL |
      | Y Column Names | RACE         |
      | Viewer Type    | Pie chart    |
      | On Click       | Select       |
      | Row Source     | Selected     |
    And no rows should be selected
    And trellis plot viewer should show 0 rows
    When user selects all rows
    Then trellis plot viewer should show 1000 rows
    And the "cells" reading of trellis plot viewer should be 16
    When user clears the row selection
    Then trellis plot viewer should show 0 rows
    When user sets "Row Source" property of trellis plot viewer to "All"
    Then trellis plot viewer should show 1000 rows
    And the "cells" reading of trellis plot viewer should be 16
    And no errors should have been logged

  Scenario: With nothing selected a Selected trellis comes back from a project alive
    When user sets properties of trellis plot viewer:
      | On Click   | Select   |
      | Row Source | Selected |
    Then no rows should be selected
    And trellis plot viewer should show 0 rows
    When user saves the current view as project "bdd-trellis-empty"
    And user closes all views
    And user opens the "bdd-trellis-empty" project
    Then the open tableview should have 1 trellis plot viewer
    And "Row Source" property of trellis plot viewer should be "Selected"
    And no rows should be selected
    And trellis plot viewer should show 0 rows
    When user selects all rows
    Then trellis plot viewer should show 1000 rows
    And the "cells" reading of trellis plot viewer should be 16
    And no error or warning balloon should have been shown
    And no errors should have been logged
