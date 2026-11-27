@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Moving rows between the selection, the filter and a
  column: Select > Missing Values selects the rows empty in the ticked column, Select > Selection to Filter shows
  only those rows, and Select > Filter to Column keeps the filter as a true/false column. demog (5850 rows, 751 of
  them with no HEIGHT) is opened off camera, the grid showing six of its columns; the column the clip adds lives
  in the page's table only and goes with it.

  Scenario: Turn a selection into a filter
    Select the rows with missing values, show only them, and save that as a column.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the grid shows only the columns "USUBJID, AGE, SEX, RACE, HEIGHT, WEIGHT"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    Then no rows should be selected
    # caption: Select > Missing Values finds the empty cells
    When user picks "Select > Missing Values..." from the top menu
    Then Select Missing Values dialog should be visible
    When user toggles the "HEIGHT" column in the column list of Select Missing Values dialog
    And user clicks on OK button in Select Missing Values dialog
    Then 751 rows should be selected
    And all rows should pass the filter
    When user hovers over "Selected: 751" text in status bar
    And the video says "751 rows have no HEIGHT; they are selected"
    # caption: Selection to Filter shows only the selected rows
    When user picks "Select > Selection to Filter" from the top menu
    Then 751 rows should pass the filter
    And every selected row should pass the filter
    # caption: Filter to Column keeps the filter as a column
    When user picks "Select > Filter to Column..." from the top menu
    Then Name input in "Filter to Column" dialog should have the value "Filter"
    When user enters "No HEIGHT" into Name input in "Filter to Column" dialog
    And user clicks on OK button in "Filter to Column" dialog
    Then the "Filter to Column" dialog should close
    And "No HEIGHT" column should have type "bool"
    And the filter should pass exactly the rows where "No HEIGHT" is "true"
    When the video says "The new column is true in exactly these 751 rows"
    Then no errors should have been logged
