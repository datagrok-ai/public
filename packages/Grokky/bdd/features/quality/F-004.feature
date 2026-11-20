@serial
Feature: Asked to order a bar chart's categories, the assistant changes the order the bars are drawn in
  pilot-1120 F-004: asked to put the age groups of a bar chart youngest to oldest, the assistant set the
  column's .category-order tag under "by category", or barSortType custom with categoryOrder as a plain
  array. The viewer drops both silently (categoryOrder is a JSON string {"manual": [...]}) and echoes the
  array back, so the read-back passed while the bars stayed in the old order. The check reads the bars
  as drawn, from the bar chart's hit areas.

  These scenarios reach a Docker container and a paid model: they run only with GROKKY_BDD_LIVE=1.

  Background:
    Given live AI model calls are allowed
    And the AI runtime traffic is recorded
    And user is logged in
    And the AI runtime container is healthy
    And user opens the demog demo table
    And the AI panel is open on a new chat

  Scenario: The age groups of a bar chart are drawn youngest to oldest
    Given user adds a calculated column "age group" with formula 'If(${age} < 30, "under 30", If(${age} <= 50, "30 to 50", "over 50"))'
    And user adds a bar chart viewer with:
      | splitColumnName | age group |
      | stackColumnName | sex       |
      | valueAggrType   | count     |
    When user asks the AI panel "The order of the bars looks odd, can you put the age groups youngest to oldest?"
    Then the bar chart should draw its bars in the order "under 30, 30 to 50, over 50"
    And a screenshot "F-004-bar-order" of the AI panel is saved
