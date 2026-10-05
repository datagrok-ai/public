@viewers @realizes:charts.viewer.surface-plot @realizes:charts.viewer.globe @realizes:charts.viewer.group-analysis
Feature: The Add viewer gallery for Charts viewers, Surface plot, Globe and Group Analysis
  The ribbon's Add viewer gallery disables a Charts viewer when the table cannot feed it and says
  why in the card's tooltip; a card is disabled through aria-disabled, so its state is read as it
  is. Surface plot and Globe draw 3D views and Group Analysis a grid of groups with analysed
  columns; they report no readings of their own, so the claims are their properties, the painted
  canvas, the inner grid's own readings and the console.
  Translated from the TestTrack case Charts/charts-other-viewers. Kept without (see the request
  document): that the Globe painted — it draws with WebGL and neither pixel step reads its canvas, so
  it is claimed visible and silent, and so is the Surface plot (WebGL too); the Group Analysis layout is saved and applied through the API,
  not through View > Layout. The Surface plot is checked on demog: on
  earthquakes.csv it is not shown.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: The gallery disables the Charts viewers a table cannot feed and says why
    Given user opens demog dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    And the following elements should be enabled:
      | first "Sankey" card in "Add Viewer" dialog         |
      | first "Chord" card in "Add Viewer" dialog          |
      | first "Timelines" card in "Add Viewer" dialog      |
      | first "Radar" card in "Add Viewer" dialog          |
      | first "Sunburst" card in "Add Viewer" dialog       |
      | first "Tree" card in "Add Viewer" dialog           |
      | first "Surface plot" card in "Add Viewer" dialog   |
      | first "Globe" card in "Add Viewer" dialog          |
      | first "Group Analysis" card in "Add Viewer" dialog |
      | first "Word cloud" card in "Add Viewer" dialog     |
    When user presses Escape
    Then "Add Viewer" dialog should be absent
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Chem tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Chem---chem_standards.csv tree node inside browse tree
    Then the "chem_standards" view should be current
    And the table should have 2 columns
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    And first "Radar" card in "Add Viewer" dialog should be disabled
    And first "Sankey" card in "Add Viewer" dialog should be disabled
    And first "Sunburst" card in "Add Viewer" dialog should be enabled
    And first "Tree" card in "Add Viewer" dialog should be enabled
    When user hovers over first "Radar" card in "Add Viewer" dialog
    Then tooltip should contain text "Radar viewer needs at least 1 numerical column"
    When user hovers over first "Sankey" card in "Add Viewer" dialog
    Then tooltip should contain text "Sankey viewer needs at least 2 string columns with less than 50 categories and 1 numerical column"
    When user presses Escape
    Then "Add Viewer" dialog should be absent
    And no errors should have been logged

  Scenario: Globe draws on earthquakes
    Given user opens earthquakes dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Globe" card in "Add Viewer" dialog
    Then globe viewer should be visible
    And no errors should have been logged
    When user clicks on close icon of globe viewer
    Then globe viewer should be absent
    And no errors should have been logged

  Scenario: Surface plot on demog takes its columns, Projection and Wireframe
    Given user opens demog dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Surface plot" card in "Add Viewer" dialog
    Then surface plot viewer should be visible
    And "XColumnName" property of surface plot viewer should not be ""
    And "YColumnName" property of surface plot viewer should not be ""
    And "ZColumnName" property of surface plot viewer should not be ""
    And no errors should have been logged
    When user clicks on grid
    And user clicks on settings icon of surface plot viewer
    Given "Misc" category in context panel is expanded
    Then "Projection" property in context panel should be visible
    When user selects "orthographic" in "Projection" property in context panel
    Then "Projection" property of surface plot viewer should be "orthographic"
    When user unchecks "Wireframe" property in context panel
    Then "Wireframe" property of surface plot viewer should be "false"
    And surface plot viewer should be visible
    And no errors should have been logged

  Scenario: Group Analysis adds an analysed column and keeps it in a layout (GROK-19039, GROK-19047)
    Given user opens demog dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Group Analysis" card in "Add Viewer" dialog
    Then group analysis viewer should be visible
    When user clicks on first grid
    And user clicks on settings icon of group analysis viewer
    Then "Group By" property in context panel should be visible
    When user clicks on "..." button in "Group By" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user toggles the "SEX" column in the column list of "Select columns..." dialog
    Then "1 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Group By" property of group analysis viewer should be "SEX"
    And the "rows" reading of grid in group analysis viewer should be 2
    And the "text of cell 1 of SEX" reading of grid in group analysis viewer should be "F"
    And the "text of cell 2 of SEX" reading of grid in group analysis viewer should be "M"
    And grid in group analysis viewer should not have a "header min(AGE)" area
    When user clicks on "Add column to analyze" icon in group analysis viewer
    Then "Add column" dialog should be visible
    When user selects "AGE" in Column input in "Add column" dialog
    And user clicks on OK button in "Add column" dialog
    Then "Add column" dialog should be absent
    And grid in group analysis viewer should have a "header min(AGE)" area
    And the "text of cell 1 of min(AGE)" reading of grid in group analysis viewer should be "18.00"
    And no errors should have been logged
    When user saves the layout of the current table view to the server
    And user clicks on close icon of group analysis viewer
    Then group analysis viewer should be absent
    When user loads the saved layout
    Then group analysis viewer should be visible
    And "Group By" property of group analysis viewer should be "SEX"
    And the "rows" reading of grid in group analysis viewer should be 2
    And grid in group analysis viewer should have a "header min(AGE)" area
    And no errors should have been logged
