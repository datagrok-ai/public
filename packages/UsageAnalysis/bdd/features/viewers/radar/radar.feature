@viewers @realizes:charts.viewer.radar
Feature: Radar table switch, Values, Title, Normalization, Color legend and the project round trip
  The Radar draws one line per row over up to ten numeric axes (Values), at most the first 1000 lines
  of the rows its legend lets through. It reports the axes it laid out (`axes`), the lines it drew
  (`rows shown`) and the message it shows in place of or above the chart (`message`). The viewer is
  added from the ribbon's Add viewer gallery and set up in the Context Panel. A legend click narrows
  the lines on demog-1000, where the selected category has fewer than 1000 rows; on demog the cap
  keeps 1000 either way. The project round trip is made twice: through the ribbon's Save dialog and
  with a project saved without it, both reopened from Browse > Dashboards with no error logged
  (GROK-18085, GROK-19376).
  Translated from the TestTrack case Charts/radar.

  Background:
    Given user is logged in
    And the package autostarts have completed

  Scenario: The Radar switches its table and back (GROK-18576, GROK-18935)
    earthquakes is opened from Browse, as a user opens a file: a table the test harness opens through
    the API after the Context Panel has shown the Radar is not offered in its Table list.
    Given user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    Then "Add Viewer" dialog should be visible
    When user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be visible
    And radar viewer should be bound to table "demog-1000"
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Given "Value" category in context panel is expanded
    Then "Values" property in context panel should be visible
    And "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And radar viewer should be painted
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And "Only first 1000 shown" text in radar viewer should be absent
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    And Files---Demo---geo tree node inside browse tree is expanded
    When user double-clicks on Files---Demo---geo---earthquakes.csv tree node inside browse tree
    Then the "earthquakes" view should be current
    Given user switches to the "demog-1000" table view
    Then the "demog-1000" view should be current
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Then "Table" property in context panel should be visible
    When user selects "earthquakes" in "Table" property in context panel
    Then radar viewer should be bound to table "earthquakes"
    And "Values" property of radar viewer should be "Latitude, Longitude, Depth, Magnitude, NbStations, Gap, Distance, RMS, EventID"
    And no errors should have been logged
    When user selects "demog-1000" in "Table" property in context panel
    Then radar viewer should be bound to table "demog-1000"
    And "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And no errors should have been logged

  Scenario: Values, Title, Normalization, Color and the 1000-row notice (GROK-18408, GROK-17999)
    Given user opens demog dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be bound to table "demog"
    And "Only first 1000 shown" text in radar viewer should be visible
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Given "Value" category in context panel is expanded
    Then "Values" property in context panel should be visible
    When user clicks on "..." button in "Values" property in context panel
    Then "Select columns..." dialog should be visible
    When user clicks on "None" link in "Select columns..." dialog
    Then "0 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be visible
    And no errors should have been logged
    When user clicks on "..." button in "Values" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "AGE" column in the column list of "Select columns..." dialog
    And user toggles the "HEIGHT" column in the column list of "Select columns..." dialog
    Then "2 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Values" property of radar viewer should be "AGE, HEIGHT"
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And radar viewer should be painted
    Given "Description" category in context panel is expanded
    When user enters "Body measures" into "Title" property in context panel
    Then title of radar viewer should have text "Body measures"
    When user takes a snapshot of radar viewer
    And user selects "Global" in "Normalization" property in context panel
    Then radar viewer should have repainted
    When user takes a snapshot of radar viewer
    And user selects "Column" in "Normalization" property in context panel
    Then radar viewer should have repainted
    And no errors should have been logged
    When user clicks on "..." button in "Values" property in context panel
    Then "Select columns..." dialog should be visible
    When user toggles the "WEIGHT" column in the column list of "Select columns..." dialog
    Then "3 checked" text in "Select columns..." dialog should be visible
    When user clicks on OK button in "Select columns..." dialog
    Then "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    Given "Color" category in context panel is expanded
    When user selects "SEX" in "Color" property in context panel
    Then the legend of radar viewer should list 2 items
    And "F" legend item in legend of radar viewer should be visible
    And "M" legend item in legend of radar viewer should be visible
    And "M" legend item in legend of radar viewer should not be selected
    When user clicks on "M" item in the legend of radar viewer
    Then "M" legend item in legend of radar viewer should be selected
    And "F" legend item in legend of radar viewer should not be selected
    And 5850 rows should pass the filter
    When user clicks on "M" item in the legend of radar viewer
    Then "M" legend item in legend of radar viewer should not be selected
    And 5850 rows should pass the filter
    And no errors should have been logged

  Scenario: A legend click narrows the lines the Radar draws
    Given user opens demog-1000 dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be bound to table "demog-1000"
    And the "axes" reading of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And the "rows shown" reading of radar viewer should be 1000
    And the "message" reading of radar viewer should be ""
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Given "Color" category in context panel is expanded
    When user selects "SEX" in "Color" property in context panel
    Then the legend of radar viewer should list 2 items
    When user clicks on "M" item in the legend of radar viewer
    Then "M" legend item in legend of radar viewer should be selected
    And the "rows shown" reading of radar viewer should be 447
    And 1000 rows should pass the filter
    When user clicks on "M" item in the legend of radar viewer
    Then "M" legend item in legend of radar viewer should not be selected
    And the "rows shown" reading of radar viewer should be 1000
    And no errors should have been logged

  Scenario: A Radar rebound to another table survives a project save through the ribbon and reopen
    Given no project named "RadarRebind{time}" is on the server
    And user opens demog-1000 dataset
    And user opens spgi dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be bound to table "spgi-100"
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Given "Value" category in context panel is expanded
    Then "Values" property in context panel should be visible
    When user selects "demog-1000" in "Table" property in context panel
    Then radar viewer should be bound to table "demog-1000"
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "RadarRebind{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share RadarRebind{time}" dialog should be visible
    When user clicks on CANCEL button in "Share RadarRebind{time}" dialog
    Then the "Share RadarRebind{time}" dialog should close
    And 1 project named "RadarRebind{time}" should be on the server
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "RadarRebind{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on RadarRebind{time} gallery card
    Then the task bar should have finished "Opening project"
    And the open table views should be exactly "demog-1000, spgi-100"
    Given user switches to the "spgi-100" table view
    Then the "spgi-100" view should be current
    And radar viewer should be bound to table "demog-1000"
    And "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And radar viewer should be painted
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And no errors should have been logged

  Scenario: Reopening a project saved without the dialog brings the rebound Radar back (GROK-18085, GROK-19376)
    Given no project named "RadarReopen{time}" is on the server
    And user opens demog-1000 dataset
    And user opens spgi dataset
    When user clicks on "Add viewer" icon
    And user clicks on first "Radar" card in "Add Viewer" dialog
    Then radar viewer should be bound to table "spgi-100"
    When user clicks on grid
    And user clicks on settings icon of radar viewer
    Given "Value" category in context panel is expanded
    Then "Values" property in context panel should be visible
    When user selects "demog-1000" in "Table" property in context panel
    Then radar viewer should be bound to table "demog-1000"
    And no errors should have been logged
    When user saves all open table views as project "RadarReopen{time}"
    Then 1 project named "RadarReopen{time}" should be on the server
    When user picks "Close All" from the context menu of browse tab
    Then the "Home" view should be current
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "RadarReopen{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Given user watches the task bar
    When user double-clicks on RadarReopen{time} gallery card
    Then the task bar should have finished "Opening project"
    And the open table views should be exactly "demog-1000, spgi-100"
    Given user switches to the "spgi-100" table view
    Then the "spgi-100" view should be current
    And radar viewer should be bound to table "demog-1000"
    And "Values" property of radar viewer should be "AGE, HEIGHT, WEIGHT"
    And radar viewer should be painted
    And "The Radar viewer requires a minimum of 1 numerical column." text in radar viewer should be absent
    And no errors should have been logged
