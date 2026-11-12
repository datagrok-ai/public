@serial
Feature: The AI panel works on the open table
  A user's analysis requests in the AI panel, carried out by the Hermes runtime with its tools in the
  page. Every claim reads the outcome from the table, its viewers and the grid; a number the user asked
  for is checked against the value computed from the table in the page.

  These scenarios reach a Docker container and a paid model, which the library's hard rule keeps out of
  the regular suite: they run only with GROKKY_BDD_LIVE=1 and skip otherwise. On a stand without a queue
  worker every turn first waits out the runtime's 60 s admission window. @serial: one chat at a time, so
  turns do not compete for the runtime's admission queue. demog-1000 has 128 rows without a HEIGHT.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens demog-1000 dataset
    And the AI panel is open on a new chat

  Scenario: A small dashboard is built from one request, then changed by a follow-up
    When user asks the AI panel "Help me build a small dashboard on this table. First remove the rows where HEIGHT or WEIGHT is missing. Then add a column named BMI with each patient's body-mass index (WEIGHT is in kg, HEIGHT in cm). Then add a scatter plot with HEIGHT on the X axis and WEIGHT on the Y axis, colored by SEX, and a box plot of BMI for each RACE. Finally, filter the table to patients aged 50 or older."
    Then the table should have 872 rows
    And "HEIGHT" column should have no missing values
    And every value of "BMI" column should be the BMI of "WEIGHT" in kg and "HEIGHT" in cm
    And the open tableview should have 1 Scatter plot viewer
    And properties of scatter plot viewer should be:
      | xColumnName     | HEIGHT |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    And the open tableview should have 1 Box plot viewer
    And properties of box plot viewer should be:
      | valueColumnName     | BMI  |
      | category1ColumnName | RACE |
    And the filter should pass exactly the rows where "AGE" is between 50 and 200
    And scatter plot viewer should show every row that passes the filter of its table
    When user asks the AI panel "On the scatter plot you just made, color the points by RACE instead of SEX."
    Then the open tableview should have 1 Scatter plot viewer
    And properties of scatter plot viewer should be:
      | xColumnName     | HEIGHT |
      | yColumnName     | WEIGHT |
      | colorColumnName | RACE   |
    And the open tableview should have 1 Box plot viewer
    And every turn of the chat should have used the same AI session
    And a screenshot "table-dashboard" of the AI panel is saved

  Scenario: A question about the table is answered with the value computed from it
    When user asks the AI panel "Look at the women in this table: which race is heaviest on average? Reply with just the race and its mean weight, rounded to one decimal."
    Then the AI answer should name the "RACE" with the highest mean "WEIGHT" where "SEX" is "F"
    And the AI should have used one of the tools "datagrok_exec, terminal"
    And a screenshot "table-question" of the AI panel is saved

  Scenario: The table view's own commands color-code a column and hide another
    When user asks the AI panel "Using this table view's own commands, apply linear color coding to the AGE column and hide the CONTROL column from the grid."
    Then "AGE" column should be color-coded linearly
    And the "CONTROL" column should be hidden in the grid
    And the table should have a column "CONTROL"
    And the AI should have used the "list_view_functions" tool
    And the AI should have used the "call_view_function" tool
    And a screenshot "table-view-functions" of the AI panel is saved

  Scenario: The AI asks which chart to add and adds the one picked
    When user asks the AI panel "I'd like a chart of AGE in this table, but I'm not sure which kind. Before doing anything, ask me which chart type I want and offer me a few options." and picks the option matching "histogram" in its question
    Then the AI should have asked a question offering 2 or more options
    And the AI should have used the "AskUserQuestion" tool
    And the open tableview should have 1 Histogram viewer
    And "valueColumnName" property of histogram viewer should be "AGE"
    And a screenshot "table-ask-user" of the AI panel is saved

  Scenario: A request for data the table does not have changes nothing
    When user asks the AI panel "I'd like to see how cholesterol changes with age in these patients. If this table has no cholesterol data, don't make anything up: tell me which numeric columns it has instead."
    Then the AI answer should name every numeric column of the table
    And the current table view should show nothing but its grid
    And the table should have 1000 rows
    And the table should have 11 columns
    And a screenshot "table-missing-column" of the AI panel is saved

  Scenario: A per-group summary goes into a new table view
    When user asks the AI panel "I need a per-RACE summary of this study: how many subjects are in each group, and their average AGE and WEIGHT. Put that summary into a new table and open it in its own view."
    Then the current table should hold the count and the mean "AGE, WEIGHT" of every "RACE" group of table "demog-1000"
    And a screenshot "table-summary" of the AI panel is saved

  # TODO: core xamgle/lib/src/features/ai/table_view_prompts.dart answers "which X has the highest <agg> Y" with a bar chart of all rows
  @known-failure:87
  Scenario: A question naming an average and a column is answered, not turned into a chart
    When user submits "Among the female patients in this table, which RACE has the highest average WEIGHT, and what is that average? Reply with just the race and the value rounded to one decimal." in the AI panel
    Then a screenshot "table-question-intercepted" of the AI panel is saved
    And the AI answer should name the "RACE" with the highest mean "WEIGHT" where "SEX" is "F"

  # TODO: core xamgle/lib/src/features/ai/table_view_prompts.dart turns "recolor that scatter plot by X" into a new scatter plot of X
  @known-failure:98
  Scenario: A request to recolor the scatter plot changes it rather than adding another
    Given user adds a Scatter plot viewer with:
      | xColumnName     | HEIGHT |
      | yColumnName     | WEIGHT |
      | colorColumnName | SEX    |
    When user submits "Now recolor that scatter plot by RACE." in the AI panel
    Then a screenshot "table-recolor-intercepted" of the AI panel is saved
    And the open tableview should have 1 Scatter plot viewer
    And "colorColumnName" property of scatter plot viewer should be "RACE"
