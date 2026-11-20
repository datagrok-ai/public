@serial
Feature: A bar chart of averages "colored by" a category draws bars, not an error or a silent no-op
  pilot-1120 F-009: asked for average height by race in a bar chart colored by disease, the assistant set
  colorColumnName to the string column (the chart drew the error 'avg(disease): "avg" is not supported for
  string') or stacked by it with avg (bars stay uniform), and said the bars were colored by disease.
  datagrok_verify passed on the props alone. The bar chart reports no error in getWidgetStatus, so the
  checks read its settings: no string colour column, and no stack under avg.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens the demog demo table
    And the AI panel is open on a new chat

  Scenario: The bar chart of average height by race draws, and is not stacked by disease
    When user asks the AI panel "Which race has the highest average height? Show the answer in a bar chart colored by disease, then tell me the top race and its value."
    Then the open tableview should have 1 Bar chart viewer
    And "valueAggrType" property of bar chart viewer should be "avg"
    And "colorColumnName" property of bar chart viewer should not be "disease"
    And "stackColumnName" property of bar chart viewer should not be "disease"
    And a screenshot "F-009-bar-chart" of the AI panel is saved
