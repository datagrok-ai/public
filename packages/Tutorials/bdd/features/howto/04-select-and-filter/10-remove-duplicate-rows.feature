@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Duplicate rows: Select > Duplicates selects every row
  whose ticked columns repeat another row's, Select > Selection to Top sorts the grid so the selected rows come
  first, a Ctrl+click keeps one copy out of the selection, and Edit > Remove > Selected Rows deletes the rest.
  iris (150 flowers; rows 102 and 143 hold the same four measurements and species) is opened off camera; the rows
  are removed from the page's table only and nothing is made on the server.

  Scenario: Find and remove duplicate rows
    Select rows that repeat each other, look at them side by side, and delete the extra copy.

    Given user is logged in
    And simple mode is off
    And user opens iris dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    Then no rows should be selected
    # caption: Select > Duplicates finds rows that repeat
    When user picks "Select > Duplicates..." from the top menu
    Then Select Duplicates dialog should be visible
    # caption: Tick the columns that must match
    When user toggles the "Sepal.Length" column in the column list of Select Duplicates dialog
    And user toggles the "Sepal.Width" column in the column list of Select Duplicates dialog
    And user toggles the "Petal.Length" column in the column list of Select Duplicates dialog
    And user toggles the "Petal.Width" column in the column list of Select Duplicates dialog
    And user toggles the "Species" column in the column list of Select Duplicates dialog
    And user clicks on OK button in Select Duplicates dialog
    Then Select Duplicates dialog should be hidden
    And 2 rows should be selected
    And rows 102 to 102 should all be selected
    And rows 143 to 143 should all be selected
    # caption: Selection to Top brings them to the top of the grid
    When user picks "Select > Selection to Top" from the top menu
    And the video says "The same flower, entered twice"
    # caption: Ctrl+click one copy to keep it
    And user clicks on the "cell 102 of Species" area of grid holding Control
    Then 1 row should be selected
    # caption: Then remove the selected rows
    When user picks "Edit > Remove > Selected Rows" from the top menu
    Then the table should have 149 rows
    And no rows should be selected
    And no errors should have been logged
