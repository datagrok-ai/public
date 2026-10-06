@demo @realizes:u2.dg.pickers
Feature: Dataframe pickers
  The table, column and columns pickers need an open table: the platform base opens it and the
  demo page then opens without a reload, so the table stays.

  Scenario: Table, column and columns pickers over an open table
    Given user is logged in
    And user opens demog dataset
    And user opens the "Dataframes" demo page
    When user selects "demog" in "Open table" input
    Then value of table readout should have text "demog"
    When user selects "age" in "Demog column" input
    Then value of "demog column" readout should have text "age"
    When user types "hei" into Column input
    And user selects "height" in Column input
    Then value of column readout should have text "height"
    When user checks demog checkbox in Tables input
    Then value of tables readout should have text "demog"
    When user clicks on editor of Columns input
    And user clicks on All link
    And user clicks on OK button
    Then value of columns readout should contain text "age, sex"
    When user selects "age" in x input in Mapping input
    And user selects "height" in y input in Mapping input
    Then value of mapping readout should contain text "\"x\":\"age\""
