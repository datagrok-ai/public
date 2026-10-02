@journey @viewers @realizes:viewers.filters
Feature: Expression and text filter cards
  The expression card builds its rules from a form — a column, an operation, a value and the Add
  filter button — and keeps each rule as a row of its own: two rules combine with OR until the
  card's switch says AND, a rule's own checkbox suspends it and keeps it in the list, Remove Query
  in the rule's menu takes it out, free-text mode takes a rule typed in the card's own syntax and
  the form comes back with it, the regex operation turns a pasted comma list into an alternation
  and drops a trailing comma (GROK-20242), and the string and date operations keep the rows they
  name. The text card on beer keeps the rows whose Aroma holds any of its terms and, switched to
  AND, the rows holding all of them; at fuzziness 0 a near miss keeps nothing and a higher
  fuzziness brings rows in. The expression half runs on demog-1000 (AGE over 50 in 367 rows, HEIGHT
  under 160 in 167, both in 88, either in 446, HEIGHT under 150 in 16; SEX F 553; RACE containing
  "an" 911; STARTED after 1991 in 483; USUBJID matching 5, 15 or 25 in 357); the text half on beer,
  118 beers whose Aroma is a long text — demog has no text column (operator ruling D4); 92 of them
  hold "malt", 106 "malt" or "hop" and 88 both, whatever the case.
  Not translated: the row-by-row read of which Aroma values survive — the library's "contains" is
  case-sensitive while the card matches regardless of case, so the claim is the count computed
  independently without case; fuzziness is read as 0 only once a term shows the slider, and the
  count rises from 0.5 to 1 (the slider takes tenths: 0 keeps none, 0.5 some, 1 more); the date
  operation is claimed by its count and its rule, there being no row check for dates; the
  free-text rule is a different predicate from the form's rather than an equivalent one, since
  Enter commits it as a further rule and an equal one could not move the count.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user picks "Add Filter | Expression" from the filter panel menu
    Then "Expression" filter card should be visible
    And the "type of Expression" reading of filter panel should be "expression"
    And all rows should pass the filter

  Scenario: Two rules from the form combine with OR, and the switch makes them AND
    When user selects "AGE" in Column input in "Expression" filter card
    And user selects ">" in Operation input in "Expression" filter card
    And user types "50" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 367 rows should pass the filter
    And the filter should pass exactly the rows where "AGE" is between 51 and 89
    And the "categories of Expression" reading of filter panel should be "${AGE} > 50"
    And counter of filter panel should have text "1"
    When user selects "HEIGHT" in Column input in "Expression" filter card
    And user selects "<" in Operation input in "Expression" filter card
    And user types "160" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "${AGE} > 50, ${HEIGHT} < 160"
    And mode of "Expression" filter card should have text "OR"
    And 446 rows should pass the filter
    When user clicks on mode of "Expression" filter card
    Then mode of "Expression" filter card should have text "AND"
    And 88 rows should pass the filter
    And no errors should have been logged

  Scenario: A rule's own checkbox suspends it and keeps it in the list
    When user clicks on the "checkbox ${HEIGHT} < 160 of Expression" area of filter panel
    Then 367 rows should pass the filter
    And the "categories of Expression" reading of filter panel should be "${AGE} > 50, ${HEIGHT} < 160"
    And the "selected categories of Expression" reading of filter panel should be "${AGE} > 50"
    When user clicks on the "checkbox ${HEIGHT} < 160 of Expression" area of filter panel
    Then 88 rows should pass the filter
    And the "selected categories of Expression" reading of filter panel should be "${AGE} > 50, ${HEIGHT} < 160"
    And no errors should have been logged

  Scenario: Remove Query in a rule's menu takes the rule out
    When user picks "Remove Query" from the context menu of the "category ${AGE} > 50 of Expression" area of filter panel
    Then 167 rows should pass the filter
    And the "categories of Expression" reading of filter panel should be "${HEIGHT} < 160"
    And the filter should pass exactly the rows where "HEIGHT" is between 1 and 159.999
    And no errors should have been logged

  Scenario: Free-text mode takes a rule in the card's own syntax and the form comes back
    When user hovers over "Expression" filter card
    And user clicks on "Switch to free-text mode" icon in "Expression" filter card
    Then Column input in "Expression" filter card should be hidden
    When user types "${HEIGHT} < 150" into the search of the "Expression" filter card
    And user presses Enter
    Then the "categories of Expression" reading of filter panel should be "${HEIGHT} < 160, ${HEIGHT} < 150"
    And 16 rows should pass the filter
    When user hovers over "Expression" filter card
    And user clicks on "Switch to free-text mode" icon in "Expression" filter card
    Then Column input in "Expression" filter card should be visible
    And 16 rows should pass the filter
    And no errors should have been logged

  Scenario: The regex operation turns a pasted list into an alternation and drops a trailing comma
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "USUBJID" in Column input in "Expression" filter card
    And user selects "regex" in Operation input in "Expression" filter card
    And user pastes "5,15,25" into Value input in "Expression" filter card
    Then Value input in "Expression" filter card should have the value "5|15|25"
    When user clears Value input in "Expression" filter card
    Then Value input in "Expression" filter card should have the value ""
    When user pastes "5,15,25," into Value input in "Expression" filter card
    Then Value input in "Expression" filter card should have the value "5|15|25"
    When user clicks on "Add filter" button in "Expression" filter card
    Then the "categories of Expression" reading of filter panel should be "${USUBJID} regex 5|15|25"
    And 357 rows should pass the filter
    And no errors should have been logged

  Scenario: The string and date operations keep the rows they name
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "SEX" in Column input in "Expression" filter card
    And user selects "equals" in Operation input in "Expression" filter card
    And user types "F" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 553 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    And the "categories of Expression" reading of filter panel should be "${SEX} equals F"
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "RACE" in Column input in "Expression" filter card
    And user selects "contains" in Operation input in "Expression" filter card
    And user types "an" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 911 rows should pass the filter
    And the filter should pass exactly the rows where "RACE" contains "an"
    And the "categories of Expression" reading of filter panel should be "${RACE} contains an"
    When user picks "Remove All" from the filter panel menu
    And user picks "Add Filter | Expression" from the filter panel menu
    And user selects "STARTED" in Column input in "Expression" filter card
    And user selects "after" in Operation input in "Expression" filter card
    And user types "01/01/1991" into Value input in "Expression" filter card
    And user clicks on "Add filter" button in "Expression" filter card
    Then 483 rows should pass the filter
    And the "categories of Expression" reading of filter panel should be "${STARTED} after 01/01/1991"
    When user picks "Remove All" from the filter panel menu
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: The text card keeps the rows holding any of its terms, or all of them
    When user opens beer dataset
    And user clicks on filter icon in toolbar
    Then "Aroma" filter card should be visible
    And the "type of Aroma" reading of filter panel should be "text"
    And all rows should pass the filter
    When user types "malt" into the search of the "Aroma" filter card
    And user presses Enter
    Then 92 rows should pass the filter
    And the "categories of Aroma" reading of filter panel should be "malt"
    And counter of filter panel should have text "1"
    When user types "hop" into the search of the "Aroma" filter card
    And user presses Enter
    Then the "categories of Aroma" reading of filter panel should be "malt, hop"
    And mode of "Aroma" filter card should have text "OR"
    And 106 rows should pass the filter
    When user clicks on mode of "Aroma" filter card
    Then mode of "Aroma" filter card should have text "AND"
    And 88 rows should pass the filter
    And no errors should have been logged

  Scenario: At fuzziness 0 a near miss keeps nothing, and a higher fuzziness brings rows in
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When user types "maltx" into the search of the "Aroma" filter card
    And user presses Enter
    Then 0 rows should pass the filter
    And Fuzzyness input in "Aroma" filter card should have the value "0"
    When user drags the slider of Fuzzyness input in "Aroma" filter card to 0.5
    Then the "rows shown" reading of filter panel should be at least 1
    When user remembers the "rows shown" reading of filter panel
    And user drags the slider of Fuzzyness input in "Aroma" filter card to 1
    Then the "rows shown" reading of filter panel should be higher than remembered
    And no errors should have been logged
