@guide @tutorial
Feature: How-to · End to end
  An end-to-end how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films
  it as a demo clip that outlines what each step presses). One task from start to finish: average BMI per race, as
  a table and a chart. demog (5850 rows, 751 of them without HEIGHT) is opened off camera. Select > Missing Values
  selects the rows with no HEIGHT and Edit > Remove > Selected Rows drops them (5099 left); Add New Column makes BMI
  from WEIGHT and HEIGHT; Data > Aggregate Rows groups by RACE with avg(BMI) (the DIS_POP group, the SEVERITY pivot
  and the avg(AGE) aggregate it starts with are taken out), and ADD turns the 4 rows into the table "demog
  aggregation", where a bar chart of avg(BMI) per RACE comes from the toolbox. The saved pivot parameters are cleared
  first. Nothing is made on the server.

  Scenario: Clean and reshape a table
    Drop incomplete rows, add a calculated column, summarize by group, and chart the summary.

    Given user is logged in
    And simple mode is off
    And user clears the saved pivot table parameters
    And user opens demog dataset
    And the toolbox pane is shown
    And the context panel is closed
    And the browse panel is closed
    When the video says "The task: the average BMI of each race, as a table and a chart"
    # caption: Clean first: Select > Missing Values finds the rows with no height
    And user picks "Select > Missing Values..." from the top menu
    Then Select Missing Values dialog should be visible
    When user toggles the "HEIGHT" column in the column list of Select Missing Values dialog
    And user clicks on OK button in Select Missing Values dialog
    Then 751 rows should be selected
    # caption: Edit > Remove > Selected Rows drops them
    When user picks "Edit > Remove > Selected Rows" from the top menu
    Then the table should have 5099 rows
    # caption: Add New Column: BMI from weight and height
    When user clicks on "Add New Column..." icon
    Then "Add New Column" dialog should be visible
    When user types "BMI" into column name input
    And user replaces the code of code editor in "Add New Column" dialog with "${WEIGHT} / Pow(${HEIGHT} / 100, 2)"
    And user clicks on OK button in "Add New Column" dialog
    Then "Add New Column" dialog should be hidden
    And the table should have a column "BMI"
    # 73.2 / 1.60484^2
    And the "BMI" cell of row 1 should be displayed as "28.42"
    # caption: Reshape: Data > Aggregate Rows summarizes group by group
    When user picks "Data > Aggregate Rows..." from the top menu
    Then the open tableview should have 1 pivot table viewer
    # caption: One row per race
    When user clicks on the "remove group by chip DIS_POP" area of pivot table viewer
    And user picks "RACE" with the plus of the "group by" row of pivot table viewer
    And user clicks on the "remove pivot chip SEVERITY" area of pivot table viewer
    Then the "group by" reading of pivot table viewer should be "RACE"
    # caption: with the average BMI of each
    When user clicks on the "remove aggregate chip avg(AGE)" area of pivot table viewer
    And user picks "BMI" with the plus of the "aggregate" row of pivot table viewer
    Then the "aggregate" reading of pivot table viewer should be "avg(BMI)"
    And the "aggregated rows" reading of pivot table viewer should be 4
    # caption: ADD makes the summary a table of its own
    When user hovers over the "add to workspace" area of pivot table viewer
    And user clicks on the "add to workspace" area of pivot table viewer
    Then the "demog aggregation" view should be current
    And the table should have 4 rows
    # caption: Chart it: a bar chart of the average BMI per race
    When user clicks on "bar chart" icon in toolbox
    Then bar chart viewer should be visible
    When user picks "avg(BMI)" in the "value" column selector of bar chart viewer
    And user picks the "avg" aggregation in the "value" column selector of bar chart viewer
    Then bar chart viewer should be painted
    When the video says "From raw rows to an answer: four numbers and a chart"
    And user hovers over bar chart viewer
    Then no errors should have been logged
