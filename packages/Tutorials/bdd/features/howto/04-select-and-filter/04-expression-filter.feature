@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The expression filter: Add Filter > Expression in the
  filter panel's menu puts an Expression card on top; its form (column, operation, value, Add filter) adds one
  rule at a time, two rules combine with OR until the switch says AND, and a rule's own checkbox turns it off
  and keeps it in the list. demog (5850 rows: AGE over 50 in 2176, HEIGHT under 160 in 951, both in 499,
  either in 2628) is opened off camera with its default filter panel; nothing is made on the server.

  Scenario: Filter with an expression
    Build rules like AGE > 50 from a column, an operation and a value, and combine them with OR or AND.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens the filter panel
    Then all rows should pass the filter
    When user hovers over "Filters" dock panel
    # caption: The panel's menu adds an expression filter
    And user clicks on "font-icon-menu" icon in "Filters" dock panel
    And user picks "Add Filter | Expression" from the open menu
    Then "Expression" filter card should be visible
    And all rows should pass the filter
    # caption: Pick a column, an operation and a value
    When user selects "AGE" in Column input in "Expression" filter card
    And user selects ">" in Operation input in "Expression" filter card
    And user types "50" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 2176 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 51 and 200
    # caption: Add another rule
    When user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "<" in Operation input in "Expression" filter card
    And user types "160" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "${AGE} > 50, ${HEIGHT} < 160"
    And 2628 rows should pass the filter
    # caption: OR keeps rows that match either rule; click it for AND
    When user clicks on mode of "Expression" filter card
    Then mode of "Expression" filter card should have text "AND"
    And 499 rows should pass the filter
    # caption: A rule's checkbox turns it off
    When user clicks on the "checkbox ${HEIGHT} < 160 of Expression" area of filter panel
    Then 2176 rows should pass the filter
    And the "selected categories of Expression" reading of filter panel should be "${AGE} > 50"
    And no errors should have been logged
