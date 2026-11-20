@serial
Feature: Asked to fix uniform bar charts of averages stacked by a category, the assistant removes the stack
  pilot-1120 F-016: two bar charts of averages stacked by sex and by disease draw uniform bars, because a stack
  draws segments only under count, values or sum. Asked to fix them, the model blamed the missing column color
  coding, turned it on (which colors the grid, not the bars), kept the stacks and said the segments were now
  colored. The fix removes the stacks that draw nothing and leaves the columns' color coding alone.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens a table "patients" with:
      | study | site | sex | disease  | age | weight |
      | S1    | A    | F   | Asthma   | 34  | 61     |
      | S1    | A    | M   | Diabetes | 52  | 88     |
      | S1    | B    | F   | None     | 45  | 70     |
      | S1    | B    | M   | Asthma   | 29  | 79     |
      | S2    | A    | F   | Diabetes | 61  | 74     |
      | S2    | A    | M   | None     | 38  | 92     |
      | S2    | B    | F   | Asthma   | 57  | 66     |
      | S2    | B    | M   | Diabetes | 44  | 84     |
      | S3    | C    | F   | None     | 23  | 58     |
      | S3    | C    | M   | Asthma   | 66  | 81     |
    And user adds a bar chart viewer with:
      | splitColumnName | study  |
      | valueColumnName | weight |
      | valueAggrType   | avg    |
      | stackColumnName | sex    |
    And user adds a bar chart viewer with:
      | splitColumnName | site    |
      | valueColumnName | age     |
      | valueAggrType   | avg     |
      | stackColumnName | disease |
    And the AI panel is open on a new chat

  Scenario: The stacks under avg are removed, and the columns keep no color coding
    When user asks the AI panel "Every bar in both charts looks the same green to me. Can you check what you see and fix it?"
    Then "stackColumnName" property of first bar chart viewer should be ""
    And "stackColumnName" property of second bar chart viewer should be ""
    And "sex" column should have no color coding
    And "disease" column should have no color coding
    And no viewer of the current view should report an error
    And a screenshot "F-016-bar-charts" of the AI panel is saved
