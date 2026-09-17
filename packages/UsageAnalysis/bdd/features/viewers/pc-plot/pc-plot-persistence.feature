@journey @viewers @realizes:viewers.pc-plot
Feature: PC plot persistence and Pick Up / Apply
  A configured plot survives the three ways the platform carries a view: a layout saved to the
  server and loaded back (which also restores the viewer SET, so the scatter plot added meanwhile is
  gone), a project saved, closed and reopened, and the Pick Up / Apply pair that copies one plot's
  settings onto another without linking them — a later change to the source must not follow.
  demog-1000 with AGE, HEIGHT and WEIGHT, coloured by RACE and titled. The project and the layout
  are deleted when the feature ends.
  The Pick Up scenario configures whichever of the two plots the ordinal resolves to at that moment
  rather than assuming the one added second is the second in the DOM: a dock does not place a new
  viewer after the one already there in a fixed order, and the earlier version of this scenario —
  which configured the plot before adding the other — copied the settings the wrong way round in
  about half the runs.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a pc plot viewer with:
      | Column Names | AGE, HEIGHT, WEIGHT  |
      | Color        | RACE                 |
      | Show Title   | true                 |
      | Title        | PC Persistence Probe |
    Then pc plot viewer should show 1000 rows
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And legend of pc plot viewer should be visible

  Scenario: A saved layout restores the viewer set and the configuration
    When user saves the layout of the current table view to the server
    And user adds a scatter plot viewer
    Then the open tableview should have 1 scatter plot viewer
    When user loads the saved layout
    Then the open tableview should have 0 scatter plot viewers
    And the open tableview should have 1 pc plot viewer
    And properties of pc plot viewer should be:
      | Column Names | AGE, HEIGHT, WEIGHT  |
      | Color        | RACE                 |
      | Title        | PC Persistence Probe |
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And legend of pc plot viewer should be visible
    And the legend of pc plot viewer should list 4 items
    And pc plot viewer should show 1000 rows
    And no errors should have been logged

  Scenario: A saved project survives Close All and a reopen
    When user saves the current view as project "zz-pcplot-bdd-probe"
    And user closes all views
    And user opens the "zz-pcplot-bdd-probe" project
    Then pc plot viewer should be visible
    And properties of pc plot viewer should be:
      | Column Names | AGE, HEIGHT, WEIGHT  |
      | Color        | RACE                 |
      | Title        | PC Persistence Probe |
    And the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And pc plot viewer should show 1000 rows
    And legend of pc plot viewer should be visible
    And no errors should have been logged

  Scenario: Pick Up on one plot and Apply on another copies its settings, once
    When user adds a pc plot viewer
    Then the open tableview should have 2 pc plot viewers
    When user sets properties of first pc plot viewer:
      | Column Names    | AGE, WEIGHT, STARTED |
      | Log Columns     | AGE                  |
      | Color           | RACE                 |
      | Legend Position | Left                 |
      | Title           | Source Plot          |
    Then the axes of first pc plot viewer should be "AGE, WEIGHT, STARTED"
    And "Title" property of second pc plot viewer should not be "Source Plot"
    And "Column Names" property of second pc plot viewer should not be "AGE, WEIGHT, STARTED"
    When user picks "Pick Up / Apply > Pick Up" from the context menu of first pc plot viewer
    And user picks "Pick Up / Apply > Apply" from the context menu of second pc plot viewer
    Then properties of second pc plot viewer should be:
      | Title           | Source Plot          |
      | Color           | RACE                 |
      | Legend Position | Left                 |
      | Log Columns     | AGE                  |
      | Column Names    | AGE, WEIGHT, STARTED |
    And the axes of second pc plot viewer should be "AGE, WEIGHT, STARTED"
    When user sets "Column Names" property of first pc plot viewer to "AGE, HEIGHT, WEIGHT, STARTED"
    Then "Column Names" property of first pc plot viewer should be "AGE, HEIGHT, WEIGHT, STARTED"
    And "Column Names" property of second pc plot viewer should be "AGE, WEIGHT, STARTED"
    When user clicks on close icon of second pc plot viewer
    Then the open tableview should have 1 pc plot viewer
    When user sets properties of pc plot viewer:
      | Column Names    | AGE, HEIGHT, WEIGHT |
      | Log Columns     |                     |
      | Legend Position | Auto                |
      | Title           |                     |
      | Show Title      | false               |
      | Color           |                     |
    Then the axes of pc plot viewer should be "AGE, HEIGHT, WEIGHT"
    And legend of pc plot viewer should be hidden
    And no errors should have been logged
