@journey
Feature: The filter panel's Filter Builder filters the table by conditions built in rows or typed as a query
  PowerPack registers Filter Builder as a columnless filter, so the filter panel's menu offers it under Add
  Filter and adds its card with no column dialog. The card is u2's FilterBuilder over the table's columns: Add
  condition adds a row, the row's property, operator and value make a condition, and the complete conditions
  filter the table; removing the row lifts the filter. The icon under the rows ("Switch between the builder and
  the query text") switches the card to the query text, where a typed query filters the table the same way and
  comes back as a row. No help page covers the card; the claims
  follow PowerPack `src/filter/filter-builder-filter.ts` and the filter panel's Add Filter menu
  (`filters_core.dart`).

  Every scenario runs on a four-row table written here (v = 1..4). Left out on purpose:
  - the builder control itself (its states, read-only and invalid looks) is U2Demo's
    `features/demo/catalog.feature` and `catalog-looks.feature`; the operators and the query grammar are u2's
    unit tests (`libraries/u2`), and the molecule operators Chem's;
  - groups, and/or and not: the same tree evaluated to a bit set, one rule (u2's tests);
  - the card in a saved layout or project: what a layout stores is an API test's claim.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And user opens a table "people" with:
      | name | v |
      | Ann  | 1 |
      | Bob  | 2 |
      | Cid  | 3 |
      | Dan  | 4 |
    And user opens the filter panel
    And user picks "Add Filter | Filter Builder" from the filter panel menu

  Scenario: A condition built in a row filters the table, and removing the row lifts the filter
    Then "Add condition" link in filter builder should be visible
    And all rows should pass the filter
    When user clicks on "Add condition" link in filter builder
    And user selects "v" in property of first filter row
    And user selects "greater than" in operator of first filter row
    And user enters "2" into value of first filter row
    Then 2 rows should pass the filter
    And the filter should pass exactly the rows where "v" is between 3 and 4
    When user clicks on remove of first filter row
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A query typed into the card's query text filters the table, and the rows show it
    When user clicks on "Switch between the builder and the query text" icon in filter panel
    And user enters "v < 3" into filter query in filter panel
    Then 2 rows should pass the filter
    And the filter should pass exactly the rows where "v" is between 1 and 2
    When user clicks on "Switch between the builder and the query text" icon in filter panel
    Then operator of first filter row should have value "less than"
    And property of first filter row should have value "v"
    When user clicks on remove of first filter row
    Then all rows should pass the filter
    And no errors should have been logged
