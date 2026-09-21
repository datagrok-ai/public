@guide @help:transform
Feature: Filter one table by the rows selected in another
  A guide: the answer to "I have two tables; when I select several rows in the first one, I want
  the second one filtered to the rows with the same ID". Data > Link Tables... makes that link:
  the table you select in is the source, the one that follows is the target, the key columns are
  the ID on each side (the dialog proposes Id / Concept Id for these two tables), and the link
  type "selection to filter" says what travels. Demo pair: spgi-100 and SPGI-linked1 (3624 rows,
  of which the first five spgi-100 rows key 9).

  Scenario: Link two tables so that selecting rows in one filters the other
    Given user is logged in
    And user opens spgi-100 dataset
    And user opens spgi-linked1 dataset
    When user clicks on spgi-100 tab
    And user picks "Data > Link Tables..." from the top menu
    And user clicks on "New Link" text in "Link Tables" dialog
    And user selects "selection to filter" in Link Type input in "Link Tables" dialog
    And user clicks on LINK button in "Link Tables" dialog
    Then "spgi-100 -> SPGI-linked1" text in "Link Tables" dialog should be visible
    When user clicks on CLOSE button in "Link Tables" dialog
    And user drags the "row header 1" area of grid to the "row header 5" area
    Then 5 rows should be selected
    And 9 rows of table "SPGI-linked1" should pass the filter
    When user clicks on SPGI-linked1 tab
    Then grid should show 9 rows
