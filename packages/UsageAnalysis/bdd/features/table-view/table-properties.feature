@journey
Feature: A table made current shows its properties in the context panel, and its Rename... renames it
  There is no Table Properties dialog: a table becomes the current object from the Table Manager (Alt+T, a
  click on its row), and the context panel then shows it — General with its row and column counts, Actions
  with the table's commands. Rename... there opens the Rename table dialog (New name:), and OK renames the
  table, its view and its row in the Table Manager. Follows datagrok/concepts/table.md ("each table can be
  annotated with attributes") and help/datagrok/navigation/panels/panels.md (the context panel shows the
  current object).

  Left out on purpose: the same rename from the view tab's Table > Rename... (projects/projects-complex);
  editing a table's tags or description, which no dialog or pane offers today (MISSING.md, a question);
  the Table Manager's own listing (general/table-manager). Nothing is put on the server.

  Background:
    Given user is logged in
    And the context panel is open
    And user opens a table "tp" with:
      | name | city | weight |
      | Ann  | Oslo | 1.5    |
      | Bob  | Rome | 2.25   |
      | Cy   | Oslo | 3      |
    When user presses Alt+T
    And user clicks on the "cell 1 of name" area of Grid viewer in "Tables" dock panel
    Then the context panel should show "tp"

  Scenario: General shows the table's row and column counts, and Actions its commands
    When user expands General pane in context panel
    Then General pane in context panel should contain the text "Rows3"
    And General pane in context panel should contain the text "Columns3"
    When user expands Actions pane in context panel
    Then "Rename..." label in Actions pane in context panel should be visible
    And "Clone" label in Actions pane in context panel should be visible
    And no errors should have been logged

  Scenario: Rename... in Actions renames the table, its view and its row in the Table Manager
    When user expands Actions pane in context panel
    And user clicks on "Rename..." label in Actions pane in context panel
    Then "Rename table" dialog should be visible
    And "New name:" input in "Rename table" dialog should have value "tp"
    When user enters "tp2" into "New name:" input in "Rename table" dialog
    And user clicks on OK button in "Rename table" dialog
    Then the "Rename table" dialog should close
    And table "tp2" should be open
    And the "tp2" view should be current
    And the "text of cell 1 of name" reading of Grid viewer in "Tables" dock panel should be "tp2"
    When user presses Alt+T
    Then "Tables" dock panel should be absent
    And no errors should have been logged
