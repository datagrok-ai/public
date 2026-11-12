@demo @journey @realizes:u2.dg.property-form
Feature: All inputs
  One property set and one target object, edited in one grid row per property by the platform's
  ui.input.forProperty in the Platform section and by u2's propertyForm in the u2 section. An edit
  on either side reaches the other, and the u2 editors keep the platform's row height in one column.

  Background:
    Given user opens the "All inputs" demo page

  Scenario: The u2 editors stand in one column at the platform's row height
    Then Name input in "u2" section should be as tall as Name input in "Platform" section
    And Stage input in "u2" section should be as tall as Stage input in "Platform" section
    And Replicates input in "u2" section should be as tall as Replicates input in "Platform" section
    And editor of Name input in "u2" section should be left-aligned with editor of Stage input in "u2" section
    And editor of Stage input in "u2" section should be left-aligned with editor of Replicates input in "u2" section
    And editor of Replicates input in "u2" section should be left-aligned with editor of Active checkbox in "u2" section
    And editor of Targets input in "u2" section should be at least 150 px wide

  Scenario: The u2 checkbox row is as tall as the platform's
    Then Active checkbox in "u2" section should be as tall as Active checkbox in "Platform" section

  Scenario: The two sides never lie on one another, side by side or stacked
    Then Name input in "u2" section should not overlap Name input in "Platform" section
    And "API key" input in "u2" section should not overlap "API key" input in "Platform" section

  Scenario: The platform editors take their cell as in a form, as wide as the u2 ones
    Then editor of Name input in "Platform" section should be as wide as editor of Name input in "u2" section
    And editor of "API key" input in "Platform" section should be as wide as editor of "API key" input in "u2" section

  Scenario: The page shows no layout problems
    Then demo page should have no layout problems

  Scenario: A name typed on either side shows on the other
    When user enters "Ibuprofen" into Name input in "u2" section
    Then Name input in "Platform" section should have value "Ibuprofen"
    When user enters "Aspirin" into Name input in "Platform" section
    Then Name input in "u2" section should have value "Aspirin"

  Scenario: A stage picked on either side shows on the other
    When user selects "Phase II" in Stage input in "u2" section
    Then Stage input in "Platform" section should have value "Phase II"
    When user selects "Phase I" in Stage input in "Platform" section
    Then Stage input in "u2" section should have value "Phase I"

  Scenario: Replicates entered on either side show on the other
    When user enters "7" into Replicates input in "u2" section
    Then Replicates input in "Platform" section should have value "7"
    When user enters "3" into Replicates input in "Platform" section
    Then Replicates input in "u2" section should have value "3"

  Scenario: Active toggled on either side shows on the other
    When user unchecks Active checkbox in "u2" section
    Then Active checkbox in "Platform" section should be unchecked
    When user checks Active checkbox in "Platform" section
    Then Active checkbox in "u2" section should be checked
