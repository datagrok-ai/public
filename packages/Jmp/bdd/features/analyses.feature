@realizes:jmp.analyses
Feature: Saved JMP analyses open as Datagrok viewers, and row states come across
  A JMP table keeps its analyses as saved scripts: Graph Builder, Distribution, Fit Y by X and the
  like. When the table opens, each analysis the package can convert opens as a viewer of the table,
  titled after the script; the rest are named in a warning balloon, and the Browse preview lists
  every saved analysis with what it opens as, or why it doesn't, and its script behind a Script
  link. A row-state column comes across as two columns, "<name> marker" and "<name> color", drawn
  with the JMP markers and colours.

  No file saved by JMP itself with an analysis is in the project yet, so the analysis fixtures are
  synthetic: fixtures/saved-bar-chart.jmp and fixtures/saved-fit-model.jmp are JMPReader.jl's
  compact_UInt8.jmp (130 rows, one text column "A") with its saved "Source" script replaced by a
  Graph Builder bar chart of A, or by a Fit Model, which has no Datagrok counterpart
  (tools/make-script-fixtures.ts). fixtures/rowstate.jmp is JMPReader.jl's own: three row-state
  columns of 3 rows. The same files are published with the package under App Data/JmpReader. How each
  JMP platform maps to a viewer is the package's own test (src/tests, "JMP: analyses").

  Background:
    Given user is logged in

  Scenario: A saved Graph Builder opens as a bar chart with the table
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/saved-bar-chart.jmp" file of the project onto browse tree
    Then the "saved-bar-chart" table view should open with 130 rows
    And the open tableview should have 1 bar chart viewer
    And "Split" property of bar chart viewer should be "A"
    And "Title" property of bar chart viewer should be "Bar chart of A"
    And no viewer of the current view should report an error
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: An analysis with no Datagrok counterpart is named, and the table still opens
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/saved-fit-model.jmp" file of the project onto browse tree
    Then the "saved-fit-model" table view should open with 130 rows
    And the open tableview should have 0 bar chart viewers
    And a warning balloon containing "1 saved analysis not converted: Fit A" should have been shown
    And no errors should have been logged

  Scenario: The preview lists the saved analyses, and Open opens them
    Given no table view is open
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---JmpReader tree node inside browse tree is expanded
    And Files---App-Data---JmpReader---tests tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---JmpReader---tests---saved-bar-chart.jmp tree node inside browse tree
    Then the "saved-bar-chart.jmp" view should be current
    And "Saved analyses" section should be visible
    And "Saved analyses" section should contain text "1 of 1 open as viewers"
    And "Saved analyses" section should contain text "Bar chart of A"
    And "Saved analyses" section should be above grid
    When user clicks on Open button
    Then the "saved-bar-chart" table view should open with 130 rows
    And the open tableview should have 1 bar chart viewer
    And no errors should have been logged

  Scenario: The preview says why an analysis doesn't open
    Given the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---JmpReader tree node inside browse tree is expanded
    And Files---App-Data---JmpReader---tests tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---JmpReader---tests---saved-fit-model.jmp tree node inside browse tree
    Then the "saved-fit-model.jmp" view should be current
    And "Saved analyses" section should contain text "0 of 1 open as viewers"
    And "Saved analyses" section should contain text "Not converted: JMP's Fit Model isn't converted"
    And "Saved analyses" section should not contain text "Fit Model(Y(:A)"
    When user clicks on Script link
    Then "Saved analyses" section should contain text "Fit Model(Y(:A), Effects(:A)"
    And no errors should have been logged

  Scenario: Row states come across as marker and color columns
    Given no table view is open
    And the browse panel is open
    When user drops the "fixtures/rowstate.jmp" file of the project onto browse tree
    Then the "rowstate" table view should open with 3 rows
    And the table should have the columns "rowstate1 marker, rowstate1 color, rowstate2 marker, rowstate2 color, rowstate3 marker, rowstate3 color"
    And no errors should have been logged
    And no error or warning balloon should have been shown
