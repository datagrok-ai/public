@guide @tutorial
Feature: How-to · Data
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). Two versions of the cars demo table, made off camera:
  "cars v1" holds the first 28 cars, "cars v2" all 30 with two prices changed (audi, mazda) and mercury removed.
  Data > Compare Tables proposes the two tables, model as the key (unique in both) and every shared column as a
  value; OK opens a "diff" table of the 5 rows that differ, sorted by the key: the two changed rows with the
  changed price cell in orange (its tooltip gives the old value), mercury in red (gone from v2), volvo2 and volvo3
  in green (new in v2). Everything is a table in the page; nothing is made on the server.

  Scenario: Compare two versions of a table
    See which rows were added, removed or changed between two versions of a table.

    Given user is logged in
    And simple mode is off
    And user opens cars dataset keeping the first 28 rows as "cars v1"
    And user opens cars dataset keeping the first 30 rows as "cars v2"
    And user sets "price" column in row 2 to "18950"
    And user sets "price" column in row 11 to "7650"
    And user selects rows where "model" is "mercury"
    And user deletes the selected rows
    And user clears the row selection
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Two versions of one table: cars v1, the older one
    When user clicks on the tab of the "cars v1" view
    # caption: and cars v2, edited since
    And user clicks on the tab of the "cars v2" view
    # caption: Compare Tables: what changed from one to the other
    And user picks "Data > Compare Tables..." from the top menu
    Then Compare Tables dialog should be visible
    And "Tables" input in Compare Tables dialog should have value "cars v1"
    # caption: Key: the column that tells which row is which
    When user hovers over Key Columns input in Compare Tables dialog
    Then Key Columns input in Compare Tables dialog should contain text "model"
    And Compare Tables dialog should not contain text "is not unique"
    # caption: Values: the columns to compare, all shared ones to start with
    When user hovers over "Value Columns" text in Compare Tables dialog
    And user clicks on OK button in Compare Tables dialog
    Then Compare Tables dialog should be hidden
    And the table should have 5 rows
    And the value of "model" column in row 3 should be "mercury"
    And the value of "model" column in row 4 should be "volvo2"
    # the current row (audi, changed) is tinted green like a new row, and a current red row turns olive: the click
    # moves the current row onto a green one
    When the video says "diff: only the rows that differ" at the bottom
    # caption: Green: rows new in cars v2
    And user clicks on the "cell 4 of model" area of grid
    Then row 4 should be current
    # caption: Red: a row cars v2 no longer has
    When user hovers over the "cell 3 of model" area of grid
    # caption: Orange: a changed value; point at it for the old one
    And user hovers over the "cell 1 of price" area of grid
    Then no errors should have been logged
