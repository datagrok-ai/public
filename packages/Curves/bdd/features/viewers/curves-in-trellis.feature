@journey
Feature: The Multi curve viewer inside a Trellis plot
  Translated from TestTrack Curves/curves-trellis.md, scenarios 1-5. Picked as "Curves" in the trellis
  viewer selector, the MultiCurveViewer becomes the inner viewer of every cell (GROK-18138). A cell's
  full-screen view shows its rows and closes with its X (GROK-18172). The trellis Context Panel has a
  MultiCurveViewer tab with the inner viewer's properties (GROK-19801); Log Y set there reaches the
  trellis' inner look (GROK-18166), and setting the Curves columns keeps the panel on that tab
  (GROK-20668). On a table without curves every cell says the viewer needs a curve column, and moving
  over the cells or switching inner viewers logs no error (GROK-18171, GROK-18389).

  The journey adds a Multi curve viewer first, and Log Y set on the MultiCurveViewer tab after that is
  kept by a layout saved to the server and applied again (scenario 3, step 4; core#854).

  core#646 is open: its NullError needs a session where no curve cell has been drawn yet, which a
  journey that opens the curves table first does not give; the demog scenarios do not claim it away.
  Not claimed: scrolling speed of a tiled trellis (scenario 6), a judgement by eye.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the "Curves" package is installed

  Scenario: A Multi curve viewer is created in the session before the trellis
    Given user opens curves dataset
    And user adds MultiCurveViewer viewer
    Then no errors should have been logged

  Scenario: Curves becomes the inner viewer of a trellis split by molecule
    Given user opens curves dataset
    And user adds a trellis plot viewer with:
      | X Column Names | smiles |
    When user picks "Curves" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    And trellis plot viewer should report no error
    And no errors should have been logged

  Scenario: The full-screen view of a cell shows its row and closes with one click on X
    When user hovers over the "cell body 1,1" area of trellis plot viewer
    And user clicks on the "full screen icon" area of trellis plot viewer
    Then the "full screen rows" reading of trellis plot viewer should be 1
    When user clicks on close icon in dialog
    Then trellis plot viewer should not report a "full screen rows" reading
    And no errors should have been logged

  Scenario: The MultiCurveViewer tab lists the inner properties, Log Y reaches the trellis, Curves keeps the tab
    When user clicks on settings icon of trellis plot viewer
    And user clicks on "MultiCurveViewer" tab in context panel
    And user expands "Misc" category in context panel
    Then "Merge Column Series" property in context panel should be visible
    When user checks "Log Y" property in context panel
    Then "logY" inner property of trellis plot viewer should be "true"
    When user clicks on "..." button in "Curves" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "All" link in "Select columns..." dialog
    Then "5 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Curves" property in context panel should contain text "5 / 5"
    And "Merge Column Series" property in context panel should be visible
    And no errors should have been logged

  Scenario: A layout saved and applied again keeps Log Y of the inner viewer
    When user saves the layout of the current table view to the server
    And user loads the saved layout
    Then the open tableview should have 1 trellis plot viewer
    And "logY" inner property of trellis plot viewer should be "true"
    And no errors should have been logged

  Scenario: On a table without curves every cell says it needs a curve column
    Given user closes all views
    And user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX |
    When user picks "Curves" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    When user hovers over the "cell body F | High" area of trellis plot viewer
    And user hovers over the "cell body M | Low" area of trellis plot viewer
    Then "The MultiCurveViewer viewer requires a minimum of 1 curves column." text in trellis plot viewer should be visible
    And no errors should have been logged

  Scenario: Switching the inner viewer of a trellis on a table without curves logs no error
    When user picks "Scatter plot" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Scatter plot"
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    When user picks "Box plot" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Box plot"
    When user picks "Curves" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "MultiCurveViewer"
    And no errors should have been logged
