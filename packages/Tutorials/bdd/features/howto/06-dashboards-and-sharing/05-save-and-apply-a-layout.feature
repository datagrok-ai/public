@guide @tutorial
Feature: How-to · Dashboards
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). demog, opened as "Demographics" (the stand holds many layouts
  named "demog", and a new one would be saved as "demog (2)"), with a scatter plot of WEIGHT by HEIGHT colored by SEX
  and a histogram of AGE is arranged off camera, with the toolbox's Layouts section open (a click on its header
  holds the current object for 2 s, so a column header is clicked after it, off camera). On camera Save in the
  Layouts section saves the arrangement, named after the table; the name typed into the context panel renames it
  "Height and weight". Then demog-1000, a table with the same columns, is opened on camera from the sidebar's
  Browse, Files > Demo (a double click; the toolbox then takes Browse's place by itself), its Layouts section offers
  the layout, and a click applies it.

  The layout (under either name), with its picture, is deleted when the feature ends; one an earlier run left under
  "Height and weight" is swept at the start (the project's own step, bindings/howto-dashboards.ts).

  Scenario: Save a layout and apply it to another table
    A layout keeps the arrangement of viewers without the data, so it fits any table with the same columns.

    Given user is logged in
    And simple mode is off
    And the layouts named "Demographics, Height and weight" are deleted when the feature ends
    And the user's own layouts named "Height and weight" are removed now and at feature end
    And user opens demog dataset keeping the first 5850 rows as "Demographics"
    And user adds a scatter plot viewer with:
      | xColumnName     | HEIGHT |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And user adds a histogram viewer with:
      | valueColumnName | AGE |
    And the browse panel is closed
    And the toolbox pane is shown
    And the context panel is open
    And Layouts accordion header in toolbox is expanded
    And user clicks on the "header AGE" area of grid
    And the context panel should show "AGE"
    When the video says "A layout is the arrangement of viewers, without the data"
    # caption: Toolbox > Layouts > Save keeps this arrangement
    And user clicks on Save button in layouts pane
    Then "Demographics" layout card should be visible
    And the context panel should show "Demographics"
    # caption: It is named after the table. Give it a better name
    When user types "Height and weight" into Name details field in context panel
    And user presses Enter
    Then "Height and weight" layout card should be visible
    When the video says "Now another table with the same columns"
    And user clicks on browse tab
    Then browse panel should be visible
    # caption: Browse > Files > Demo: open demog-1000.csv
    When user expands Files tree node inside browse tree
    And user expands Files---Demo tree node inside browse tree
    And user double-clicks on Files---Demo---demog-1000.csv tree node inside browse tree
    Then the "demog-1000" view should be current
    And the open tableview should have 0 scatter plot viewers
    And toolbox should be visible
    # caption: Its Layouts section offers the saved layout
    When user expands Layouts accordion header in toolbox
    And user hovers over "Height and weight" layout card
    # caption: Click it to arrange this table the same way
    And user clicks on "Height and weight" layout card
    Then the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 histogram viewer
    And "xColumnName" property of scatter plot viewer should be "HEIGHT"
    When the video says "Same viewers, now showing this table's 1,000 rows"
    And user hovers over scatter plot viewer
    Then no errors should have been logged
