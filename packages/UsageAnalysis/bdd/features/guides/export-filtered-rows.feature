@guide @help:transform
Feature: Export only the filtered rows
  A guide: the answer to "how do you export a filtered list, without it exporting the unfiltered
  version?". An export saves the whole table, filter or not, so the filtered subset is first
  extracted into a table of its own: the "Filtered: N" counter on the status bar opens the row
  actions, and Extract Rows puts those rows into a new view. The Export (download) icon on that
  view's ribbon then saves just them. Demo: demog filtered to DIS_POP = RA (2,550 of 5,850 rows).

  Scenario: Extract the filtered rows into their own table, then export that table
    Given user is logged in
    And simple mode is off
    And user opens demog dataset
    When user clicks on filter icon in toolbar
    And user clicks on the "category RA of DIS_POP" area of filter panel
    Then 2550 rows should pass the filter
    When user clicks on "Filtered: 2,550" text in status bar
    And user picks "Extract Rows" from the open menu
    Then grid should show 2550 rows
    When user clicks on Export icon in toolbar
    And user downloads a file through "As CSV" text in toolbar
    Then the downloaded file should contain "RA"
    And the downloaded file should not contain "Psoriasis"
