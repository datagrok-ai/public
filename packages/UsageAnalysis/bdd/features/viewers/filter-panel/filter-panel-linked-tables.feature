@journey @viewers @realizes:viewers.filters
Feature: Filter panel over linked tables
  What travels over links between tables: rows selected in spgi-100 become the filter of
  SPGI-linked1 through a selection-to-filter link, a card of SPGI-linked2 narrows SPGI-linked1
  through a filter-to-filter link and nothing comes back the other way, a card of SPGI-linked1
  composes with both, and a link whose type is changed in Data > Link Tables stops narrowing the
  way the old type did (GROK-19137). One journey on the md's own tables (operator ruling D4):
  spgi-100 (100 rows), SPGI-linked1 (3624, of which the 175 keyed by an SPGI-linked2 row pass once
  the filter-to-filter link is made) and SPGI-linked2 (224). The first five spgi-100 rows
  key 9 SPGI-linked1 rows; "v ii" in link column 3 keeps 148 SPGI-linked2 rows, which leave 5 of
  those 9; "inconclusive" PAMPA leaves 2, and 51 without the selection link.
  Not translated: the four-key filter-to-filter link is made through the API — the dialog takes a
  key pair at a time, each through a column picker, and the scenarios are about what the link
  carries; the one-key selection-to-filter link is made in Link Tables itself (its keys, Id and
  Concept Id, are the pair the dialog proposes) and its type changed there. The first five rows are
  selected, and the selection cleared, through the table's API (the md allows it); the views are
  switched through the shell and each panel is opened empty through its API; the reads of the
  tables that must not move are single reads, not two a settle apart.

  Background:
    Given user is logged in
    And user opens spgi-100 dataset
    And user opens spgi-linked1 dataset
    And user opens spgi-linked2 dataset
    And the "SPGI-linked2" table is linked to the "SPGI-linked1" table by "Sample Name, link column 1, link column 2, link column 3" to "Sample Name, link column 1, link column 2, link column 3" as "filter to filter"
    And user switches to the "spgi-100" table view
    And user picks "Data > Link Tables..." from the top menu
    And user clicks on "New Link" text in "Link Tables" dialog
    And user selects "selection to filter" in Link Type input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "spgi-100 -> SPGI-linked1" text in "Link Tables" dialog should be visible
    When user clicks on CLOSE button in "Link Tables" dialog
    Then 100 rows of table "spgi-100" should pass the filter
    And 175 rows of table "SPGI-linked1" should pass the filter
    And 224 rows of table "SPGI-linked2" should pass the filter

  Scenario: Rows selected in spgi-100 become the filter of SPGI-linked1
    When user switches to the "spgi-100" table view
    And user selects the first 5 rows
    Then 5 rows should be selected
    And 100 rows of table "spgi-100" should pass the filter
    And 9 rows of table "SPGI-linked1" should pass the filter
    And 224 rows of table "SPGI-linked2" should pass the filter
    And no errors should have been logged

  Scenario: A card of SPGI-linked2 narrows SPGI-linked1 and nothing comes back
    When user switches to the "SPGI-linked2" table view
    And user opens an empty filter panel
    And user adds a card for "link column 3" to the filter panel
    And user clicks on the "category v ii of link column 3" area of filter panel
    Then 148 rows should pass the filter
    And the filter should pass exactly the rows where "link column 3" is "v ii"
    And 5 rows of table "SPGI-linked1" should pass the filter
    And 148 rows of table "SPGI-linked2" should pass the filter
    And 100 rows of table "spgi-100" should pass the filter
    And no errors should have been logged

  Scenario: A card of SPGI-linked1 composes with what both links bring
    When user switches to the "SPGI-linked1" table view
    And user opens an empty filter panel
    And user adds a card for "PAMPA Classification" to the filter panel
    And user clicks on the "category inconclusive of PAMPA Classification" area of filter panel
    Then 2 rows should pass the filter
    And no rows where "PAMPA Classification" is "> -4.5 cm/s" should pass the filter
    And no rows where "PAMPA Classification" is "<= -5.3 cm/s" should pass the filter
    And 148 rows of table "SPGI-linked2" should pass the filter
    And 100 rows of table "spgi-100" should pass the filter
    And no errors should have been logged

  Scenario: A link changed to selection-to-selection no longer narrows SPGI-linked1
    When user switches to the "spgi-100" table view
    And user clears the row selection
    Then 51 rows of table "SPGI-linked1" should pass the filter
    When user selects the first 5 rows
    Then 2 rows of table "SPGI-linked1" should pass the filter
    When user picks "Data > Link Tables..." from the top menu
    And user clicks on "spgi-100 -> SPGI-linked1" text in "Link Tables" dialog
    And user selects "selection to selection" in Link Type input in "Link Tables" dialog
    Then Link Type input in "Link Tables" dialog should have the value "selection to selection"
    When user clicks on CLOSE button in "Link Tables" dialog
    Then 51 rows of table "SPGI-linked1" should pass the filter
    And 148 rows of table "SPGI-linked2" should pass the filter
    When user switches to the "SPGI-linked1" table view
    Then 9 rows should be selected
    And 51 rows should pass the filter
    And no errors should have been logged
