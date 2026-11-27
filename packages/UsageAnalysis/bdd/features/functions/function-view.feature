@journey @realizes:views.functions
Feature: The function view shows a run's results and keeps them in step with the form
  A script's function view (/func/<Namespace>.<Name>) puts the form above the results: RUN shows the
  outputs with the viewers the script declares (`viewer:`), a scalar output as a read-only field, and a
  re-run updates the results in place, keeping the viewers as the user left them. Changing an input marks
  the results as outdated until the next RUN; Enter runs only a form whose inputs are valid; the copy icon
  gives a link with the current inputs, and a link with values opens the form with them. A built-in
  function opens by its bare name. The ribbon adds the result table to the workspace and offers Save.

  One JavaScript script (rows, slope -> a line of `rows` points and their sum) is made through the API
  and deleted, with every run the view records, at feature start and end.

  Background:
    Given user is logged in
    And simple mode is off
    And the runs of the "BDDFuncView{time}" function are deleted at the end
    And a script "BDDFuncView{time}" is on the server:
      """
      //language: javascript
      //input: int rows = 5 {caption: Rows; nullable: false}
      //input: double slope = 2 {caption: Slope}
      //output: dataframe curve {caption: Curve; viewer: Line chart(x: "x", y: "y", title: "Response") | Grid}
      //output: double total {caption: Total}
      const xs = Array.from({length: rows}, (_, i) => i);
      curve = DG.DataFrame.fromColumns([DG.Column.fromList('int', 'x', xs), DG.Column.fromList('double', 'y', xs.map((x) => slope * x))]);
      total = xs.reduce((s, x) => s + slope * x, 0);
      """

  Scenario: RUN shows the results with the viewers the script declares
    When user opens the address "/func/Admin.BDDFuncView{time}"
    Then "Rows" input should have value "5"
    When user clicks on RUN button
    Then run header should contain text "Run #1"
    And "Line chart" viewer should be visible
    And function results should contain text "Response"
    And "Grid" viewer should be visible
    And "Total" input should have value "20"

  Scenario: A changed input marks the results as outdated until the next RUN
    When user enters "7" into "Rows" input
    Then function results should contain text "Inputs changed"
    When user clicks on RUN button
    Then function results should not contain text "Inputs changed"
    And "Total" input should have value "42"
    And run header should contain text "Run #2"

  Scenario: A re-run keeps a viewer the user closed closed, and Ctrl+Z brings it back
    When user hovers over line chart viewer
    And user clicks on close icon of line chart viewer
    Then "Line chart" viewer should not be visible
    When user clicks on RUN button
    Then run header should contain text "Run #3"
    And "Line chart" viewer should not be visible
    When user presses Control+Z
    Then "Line chart" viewer should be visible

  Scenario: Enter runs only a form whose inputs are valid
    When user clears "Rows" input
    And user presses Enter in "Rows" input
    Then run header should contain text "Run #3"
    When user enters "3" into "Rows" input
    And user presses Enter in "Rows" input
    Then run header should contain text "Run #4"
    And "Total" input should have value "6"

  Scenario: The copy icon gives a link with the current inputs, and a link opens the form with its values
    When user clicks on copy icon
    Then the clipboard should contain text "/func/Admin.BDDFuncView{time}?rows=3"
    When user opens the address "/func/Admin.BDDFuncView{time}?rows=9&slope=1"
    Then "Rows" input should have value "9"
    And "Slope" input should have value "1"

  Scenario: A built-in function opens by its bare name
    When user opens the address "/func/Sin?x=1"
    Then "X" input should have value "1"
    And RUN button should be visible

  Scenario: The ribbon adds the result table to the workspace and offers to save the run as a project
    When user opens the address "/func/Admin.BDDFuncView{time}?rows=5"
    And user clicks on RUN button
    Then run header should contain text "Run #"
    When user clicks on "Save" button
    Then "Save project" dialog should be visible
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user clicks on "Add table to workspace" icon
    Then the "curve" table view should open with 5 rows
