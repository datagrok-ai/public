@journey @realizes:powerpack.dialogs.add-new-column
Feature: Add New Column: a script of the user's in a formula
  A JavaScript script with a double input and a double output, saved through the platform's API,
  is called from a formula in the Add New Column dialog on a three-row table written in the feature.
  On a column of fractional numbers the preview computes it and OK adds a double column. Translated
  from the QA round of 13 Oct 2026 on calculated columns. The script runs in the browser, so the
  feature needs no server-side kernel.

  Not translated: the script's own result for each argument is ApiTests "DataFrame: Calculated
  columns" (with script formula); a double input given whole numbers (core#806) is the ddt test
  `calculated_column_test.dart` that came with its fix.

  Background:
    Given user is logged in
    And a script "BddDouble{time}" is on the server:
      """
      //language: javascript
      //input: double x
      //output: double y
      y = x * 2;
      """
    And user opens a table "calc" with:
      | f   |
      | 1.5 |
      | 2.5 |
      | 4.5 |

  Scenario: On a column of fractional numbers the script computes every row
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "twice" into column name input
    And user types "BddDouble{time}(${f})" into formula editor
    Then the preview grid should show "twice" computed as numbers
    When user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And "twice" column should have type "double"
    And the value of "twice" column in row 3 should be "9"
    And no errors should have been logged
