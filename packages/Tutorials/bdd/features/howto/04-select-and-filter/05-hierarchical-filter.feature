@guide @tutorial @viewport:1280x720
Feature: How-to · Select and filter
  A short how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it
  as a demo clip that outlines what each step presses). The hierarchical filter: Add Filter > Hierarchical in the
  filter panel's menu puts a tree card on top (on demog it guesses SEX / DEMOG / USUBJID); its tree icon opens
  the column list, where the levels are ticked in order (SEX, then RACE); a click on a row's label keeps that
  branch alone, its triangle opens the next level with its own counts, and the box beside a row leaves that
  branch out. The rows are the "hierarchy node" kind of bindings/howto-filters.ts. demog (5850 rows: F 3243, of
  them Caucasian 2823) is opened off camera with its default filter panel; nothing is made on the server.

  Scenario: Filter with a hierarchy
    Pick columns as levels of a tree, then keep or leave out whole branches.

    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    And user opens the filter panel
    Then all rows should pass the filter
    # caption: The panel's menu adds a hierarchical filter
    When user picks "Add Filter | Hierarchical" from the filter panel menu
    Then "SEX / DEMOG / USUBJID" filter card should be visible
    When user hovers over "SEX / DEMOG / USUBJID" filter card
    # caption: The tree icon picks the columns, one per level
    And user clicks on "Select columns" icon in "SEX / DEMOG / USUBJID" filter card
    Then "Select columns..." dialog should be visible
    When user toggles the "DEMOG" column in the column list of "Select columns..." dialog
    And user toggles the "USUBJID" column in the column list of "Select columns..." dialog
    And user toggles the "RACE" column in the column list of "Select columns..." dialog
    And user clicks on OK button in "Select columns..." dialog
    Then "SEX / RACE" filter card should be visible
    And all rows should pass the filter
    # caption: Click a branch to keep only its rows
    When user clicks on label of "F" hierarchy node
    Then 3243 rows should pass the filter
    And the filter should pass exactly the rows where "SEX" is "F"
    # caption: The triangle opens the next level
    When user clicks on twistie of "F" hierarchy node
    Then "Caucasian" hierarchy node should be visible
    # caption: Untick a branch to leave it out
    When user clicks on checkbox of "Caucasian" hierarchy node
    Then 420 rows should pass the filter
    And no rows where "RACE" is "Caucasian" should pass the filter
    And no rows where "SEX" is "M" should pass the filter
    And no errors should have been logged
