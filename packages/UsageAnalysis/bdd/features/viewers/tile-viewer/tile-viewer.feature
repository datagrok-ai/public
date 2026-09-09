@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer property surface
  Every card is one row, built from the sketch form the viewer keeps in its look, and every claim
  below is what the viewer reports about the cards it laid out: the composition of a card, the
  text a field shows, the rows the cards are built from (Row Source, the viewer's own Filter
  formula, the table's filter), the lane headers Tiles Font sizes, the title and the description,
  and the table the cards are built for. Lanes, selection, the form designer, the mirroring of the
  frame and persistence have features of their own. One journey on demog-1000 — 1000 rows over 11
  columns, SEX F 553 / M 447, RACE Caucasian 896 / Other 62 / Black 27 / Asian 15, 367 rows over
  50 — and every scenario puts back what it changed.

  The card holds ten fields (`SketchForm.defaultForm` takes ten columns), so on this table exactly
  one column is left off; the form orders its fields by a relevance score and the reading that
  lists them comes back in a different order every refresh, so the claims name members and count
  them, never the order. SEVERITY is the column the score leaves over here, verified on the stand.

  Clearing the title leaves the header with no title at all rather than falling back to the viewer
  type, so the cleared titlebar text is empty and therefore hidden; rebinding the Table to
  spgi-100 builds a three-field card rather than the ten the cap allows, so the rebind is claimed
  on what the card holds and not on a count.

  Not translated: the viewer's own context menu, which needs a region no card covers and so lives
  in the lanes journey; the four items of `tile-viewer-ui.md`, which are human judgements (colour
  legibility, menu-opening cost) or have no addressable surface (the lanes column grid, the
  Edit-Form-on-another-table probe); and closing the filter panel, which belongs to
  `features/viewers/filter-panel/`.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer
    Then tile viewer should be visible
    And the "lanes" reading of tile viewer should be 1
    And the "single lane" reading of tile viewer should be "true"
    And tile viewer should show 1000 rows
    And the "fields shown" reading of tile viewer should be 10

  Scenario: A card per row, showing the grid's display string for every field it holds
    Then the table should have 11 columns
    And the "table" reading of tile viewer should be "demog-1000"
    And the "auto generate" reading of tile viewer should be "true"
    And the "form designed" reading of tile viewer should be "false"
    And the "fields" reading of tile viewer should contain "AGE"
    And the "fields" reading of tile viewer should contain "USUBJID"
    And the "fields" reading of tile viewer should contain "WEIGHT"
    And the "fields" reading of tile viewer should not contain "SEVERITY"
    And the "tiles" reading of tile viewer should be at least 3
    And the "tiles" and "tiles in lane All rows" readings of tile viewer should be the same
    And tile viewer should have a "tile of row 1" area
    And tile viewer should have a "field AGE of row 1" area
    And tile viewer should have a "label AGE of row 1" area
    And tile viewer should not have a "field SEVERITY of row 1" area
    And the "current row" reading of tile viewer should be 1
    And the "lane of row 1" reading of tile viewer should be "All rows"
    And the "USUBJID of row 1" reading of tile viewer should be "X0273T21000300003"
    And the "AGE of row 1" reading of tile viewer should be "26"
    And the "SEX of row 1" reading of tile viewer should be "F"
    And the "RACE of row 1" reading of tile viewer should be "Caucasian"
    And the "DIS_POP of row 1" reading of tile viewer should be "Indigestion"
    And the "AGE of row 2" reading of tile viewer should be "30"
    And the "RACE of row 2" reading of tile viewer should be "Other"
    And the "WEIGHT of row 1" reading of tile viewer should be "74.10"
    And the "WEIGHT of row 1" reading of tile viewer should not be "74.1"
    And no errors should have been logged

  Scenario: A right-click on a card opens the column's menu, and the viewer has no menu region here
    Then tile viewer should not have a "viewer menu" area
    When user right-clicks on the "field AGE of row 1" area of tile viewer
    Then the open menu should list "Remove"
    And the open menu should list "Rename..."
    And the open menu should not list "Edit Form..."
    When user closes the context menu
    And user right-clicks on the "label AGE of row 1" area of tile viewer
    Then the open menu should list "Rename..."
    And the open menu should not list "Show Empty Lanes"
    When user closes the context menu
    Then no errors should have been logged

  Scenario: Tiles Font sizes the lane headers
    Given user sets "Lanes Column Name" property of tile viewer to "SEX"
    Then the "lane names" reading of tile viewer should be "F, M"
    And tile viewer should have a "lane header F" area
    And the "tiles font" reading of tile viewer should be 'normal normal 13px "Roboto"'
    When user sets "Tiles Font" property of tile viewer to 'normal normal 18px "Roboto"'
    Then the "tiles font" reading of tile viewer should be 'normal normal 18px "Roboto"'
    And the "lane header F" area of tile viewer should be taller than before
    When user sets "Tiles Font" property of tile viewer to 'normal normal 13px "Roboto"'
    Then the "tiles font" reading of tile viewer should be 'normal normal 13px "Roboto"'
    And the "lane header F" area of tile viewer should be shorter than before
    When user sets "Lanes Column Name" property of tile viewer to ""
    Then the "single lane" reading of tile viewer should be "true"
    And no errors should have been logged

  Scenario: Title and description
    When user sets properties of tile viewer:
      | Show Title | true  |
      | Title      | Cards |
    Then title of tile viewer should be visible
    And title of tile viewer should have text "Cards"
    When user sets properties of tile viewer:
      | Description                 | A card per patient |
      | Description Visibility Mode | Always             |
    Then description of tile viewer should have text "A card per patient"
    And the description of tile viewer should be above its content
    When user sets "Description Position" property of tile viewer to "Bottom"
    Then description of tile viewer should be visible
    And the description of tile viewer should be below its content
    When user sets "Description Position" property of tile viewer to "Top"
    Then the description of tile viewer should be above its content
    When user sets "Description Visibility Mode" property of tile viewer to "Never"
    Then description of tile viewer should be absent
    When user sets "Title" property of tile viewer to ""
    Then title of tile viewer should be hidden
    When user sets properties of tile viewer:
      | Show Title                  | false |
      | Description                 |       |
      | Description Visibility Mode | Auto  |
    Then no errors should have been logged

  Scenario: Row Source picks the rows the cards are built from
    Then "Row Source" property of tile viewer should be "Filtered"
    And tile viewer should show 1000 rows
    When user selects rows where "RACE" is "Asian"
    Then 15 rows should be selected
    When user sets "Row Source" property of tile viewer to "Selected"
    Then tile viewer should show 15 rows
    And tile viewer should show fewer rows than before
    And every tile of tile viewer should show "Asian" in "RACE"
    And the "selected rows shown" reading of tile viewer should be "false"
    When user sets "Row Source" property of tile viewer to "All"
    Then tile viewer should show 1000 rows
    And tile viewer should show more rows than before
    And the "selected rows shown" reading of tile viewer should be "true"
    And the "RACE of row 1" reading of tile viewer should be "Caucasian"
    When user sets "Row Source" property of tile viewer to "Filtered"
    And user clears the row selection
    Then tile viewer should show 1000 rows
    And no rows should be selected
    And no errors should have been logged

  Scenario: The viewer's own Filter formula narrows the cards and leaves the table's filter alone
    Then tile viewer should show 1000 rows
    And all rows should pass the filter
    When user sets "Filter" property of tile viewer to "${AGE} > 50"
    Then tile viewer should show 367 rows
    And tile viewer should show fewer rows than before
    And every tile of tile viewer should show a value between 51 and 89 in "AGE"
    And 1000 rows should pass the filter
    And all rows should pass the filter
    When user sets "Filter" property of tile viewer to ""
    Then tile viewer should show 1000 rows
    And tile viewer should show more rows than before
    And no errors should have been logged

  Scenario: A filter on the table reaches the cards
    Then tile viewer should show 1000 rows
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And tile viewer should show 447 rows
    And tile viewer should show fewer rows than before
    And every tile of tile viewer should show "M" in "SEX"
    When user resets the filter
    Then all rows should pass the filter
    And tile viewer should show 1000 rows
    And tile viewer should show more rows than before
    And the "SEX of row 1" reading of tile viewer should be "F"
    And no errors should have been logged

  Scenario: The Table property rebinds the cards to another table
    Given user opens spgi dataset
    And user switches to the "demog-1000" table view
    When user sets "Table" property of tile viewer to "spgi-100"
    Then tile viewer should be bound to table "spgi-100"
    And the "table" reading of tile viewer should be "spgi-100"
    And tile viewer should show 100 rows
    And the "fields" reading of tile viewer should contain "Primary Series Name"
    And the "fields" reading of tile viewer should not contain "SEX"
    And tile viewer should not have a "field SEX of row 1" area
    And tile viewer should have a "tile of row 1" area
    When user sets "Table" property of tile viewer to "demog-1000"
    Then tile viewer should be bound to table "demog-1000"
    And tile viewer should show 1000 rows
    And the "fields shown" reading of tile viewer should be 10
    And the "fields" reading of tile viewer should contain "SEX"
    And the "SEX of row 1" reading of tile viewer should be "F"
    And no errors should have been logged
