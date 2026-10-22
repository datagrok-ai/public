@serial
Feature: The AI panel carries out longer analyses on the open table
  Multi-step and multi-turn work in the AI panel, each step proven from the table, its views and viewers:
  one viewer refined over several turns, rows filtered, selected and extracted into a table of their own, a
  pivot table checked against an independent aggregation, and a summary table charted in its own view.
  demog-1000 has 494 patients aged 30 to 50, 225 of them Caucasian women.

  The table view's own prompt handler takes a request that names one viewer, and a viewer's quick commands take
  short edits ("split it by SEX", "20 bins"), before any model sees them; these prompts are worded so the platform
  passes them on.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1 and skip
  otherwise. @serial: one chat at a time.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens demog-1000 dataset
    And the AI panel is open on a new chat

  Scenario: A histogram is refined over three turns, every change landing on the same viewer
    When user asks the AI panel "Show the patients' ages as a histogram."
    Then the open tableview should have 1 Histogram viewer
    And properties of histogram viewer should be:
      | valueColumnName | AGE |
    When user asks the AI panel "Break that histogram down by SEX, so I can compare men and women."
    Then the open tableview should have 1 Histogram viewer
    And properties of histogram viewer should be:
      | valueColumnName | AGE |
      | splitColumnName | SEX |
    When user asks the AI panel "I'd like finer detail on that chart: twenty bins rather than the default."
    Then the open tableview should have 1 Histogram viewer
    And properties of histogram viewer should be:
      | valueColumnName | AGE |
      | splitColumnName | SEX |
      | bins            | 20  |
    And every turn of the chat should have used the same AI session
    And a screenshot "analysis-histogram" of the AI panel is saved

  Scenario: Rows are filtered, a subgroup selected, and the selection extracted into a table of its own
    When user asks the AI panel "In this table, show only the patients aged 30 to 50, hiding the others rather than deleting them. Among the patients shown, select the Caucasian women. Then put the selected rows into a new table and open it."
    Then the table should have 225 rows
    And every value of "AGE" column should lie between 30 and 50
    And every value of "SEX" column should match "F"
    And every value of "RACE" column should match "Caucasian"
    And a screenshot "analysis-extract" of the AI panel is saved

  Scenario: A pivot table gives the mean weight of every race, one column per sex
    When user asks the AI panel "Lay out the average WEIGHT as a cross-tab of RACE against SEX in a pivot table."
    Then the open tableview should have 1 Pivot table viewer
    And the aggregated values of pivot table viewer should match "avg(WEIGHT)" grouped by "RACE" pivoted on "SEX"
    And a screenshot "analysis-pivot" of the AI panel is saved

  Scenario: A per-race summary is charted in its own view
    When user asks the AI panel "Create a per-RACE summary table (patient count and average WEIGHT) as a new table in its own view, then chart that summary's average WEIGHT per RACE as bars in the new view."
    Then the current table should hold the count and the mean "WEIGHT" of every "RACE" group of table "demog-1000"
    And the open tableview should have 1 Bar chart viewer
    And properties of bar chart viewer should be:
      | splitColumnName | RACE |
    And a screenshot "analysis-summary-chart" of the AI panel is saved
