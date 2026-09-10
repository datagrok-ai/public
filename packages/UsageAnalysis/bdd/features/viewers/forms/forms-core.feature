@journey @viewers @realizes:viewers.forms
Feature: Forms viewer field set, row binding, sort mirroring and pinning
  The Forms viewer lays a table out as a row of cards, one field per column, with a shared label
  column down the left. It lives in `@datagrok-libraries/utils` and is registered by PowerGrid as
  `Forms`, so it is `forms viewer` here — never `form viewer`, which is the Dart `Form`, a
  different viewer.
  Card positions are fixed and are NOT the DOM order (`virtualView.refreshItem` re-appends what it
  rebuilt): position 1 is the current row while Show Current Row is on, position 2 the mouse-over
  row while Show Mouse Over Row is on, and the records follow. So `card 2` is a real, blank card
  whenever nothing is hovered — `record of card 2` is `''` and `card kind of card 2` is
  `mouse-over`. Every claim below about which rows are on the cards is `record of card <n>`, an
  expected row order; the four specs this replaces re-derived `getSortedOrder(...)` inside the page
  and compared it with a DOM scrape, so the test computed the answer it was checking.

  Two findings that shape this feature.
  GROK-20380 (Use Grid Sort OFF still mirrored the grid sort) was open when those specs were
  written and their Step 5c was wrapped in `knownOpenBug`. It does not reproduce: with the grid
  sorted and `useGridSort` off the cards return to table order and the viewer reports no sort
  column, and `forms-viewer.ts:539` gates the mirror on the flag. That scenario is a plain positive
  claim here and carries no `@known-failure` — a tag on a fixed bug is itself a failure.
  A sort from the grid's own header clears the table's current row, and the Forms viewer then
  draws NO cards at all — not the record cards of the selected rows either. The last scenario
  states the desired behaviour and is `@known-failure`; it is last because a failing scenario does
  not put back what it changed. Every other sort scenario therefore makes a row current again after
  each header double-click, which is what a user does anyway, and claims the mirroring from there.
  The double-click cycle runs BEFORE any grid sorting for a third reason: once a scenario has sorted
  the grid and reset it from the grid's context menu, the Forms viewer keeps reporting that sort
  column while the grid itself reports none, so a claim that the viewer has no sort column reads
  the residue instead. Ordered first, the scenario reads only its own state.

  Fixture: demog-1000, 11 columns, none of them `~`-prefixed. SEVERITY has exactly five Critical
  rows — 215, 304, 428, 430, 512, all SEX M, AGE 29 / 59 / 44 / 31 / 46, WEIGHT 71.5 / 123.9 / 88 /
  83 / 85 — which is what makes an expected card order writable.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer
    Then forms viewer should be visible
    And 1000 rows should pass the filter
    And the "fields shown" reading of forms viewer should be 11
    And the "record of card 1" reading of forms viewer should be 1

  Scenario: The default field set is every visible column, and the header draws all of them
    Then the table should have 11 columns
    And the "fields" reading of forms viewer should be "USUBJID, AGE, SEX, RACE, DIS_POP, HEIGHT, WEIGHT, DEMOG, CONTROL, STARTED, SEVERITY"
    And the "fields" and "header labels" readings of forms viewer should be the same
    And the "header labels" reading of forms viewer should not contain "~"
    And forms viewer should have a "label USUBJID" area
    And forms viewer should have a "remove SEVERITY" area
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The leading card is the current row and follows it
    Then the "card kind of card 1" reading of forms viewer should be "current"
    And the "card kind of card 2" reading of forms viewer should be "mouse-over"
    And the "record of card 2" reading of forms viewer should be ""
    And the "AGE of card 1" reading of forms viewer should be "26"
    When user makes row 13 current
    Then the "record of card 1" reading of forms viewer should be 13
    And the "current record" reading of forms viewer should be 13
    And the "AGE of card 1" reading of forms viewer should be "43"
    And the "USUBJID of card 1" reading of forms viewer should be "X0273T21000900008"
    When user makes row 78 current
    Then the "record of card 1" reading of forms viewer should be 78
    And the "AGE of card 1" reading of forms viewer should be "60"
    When user makes row 1 current
    Then the "AGE of card 1" reading of forms viewer should be "26"
    And no errors should have been logged

  Scenario: Show Selected Rows gives every selected row a card beyond the two leading ones
    Then "showSelectedRows" property of forms viewer should be "true"
    And the "cards" reading of forms viewer should be 2
    When user selects rows where "SEVERITY" is "Critical"
    Then 5 rows should be selected
    And the "cards" reading of forms viewer should be 7
    And the "records shown" reading of forms viewer should be 6
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    And every record card of forms viewer should show "Critical" in "SEVERITY"
    And every record card of forms viewer should show "M" in "SEX"
    And the record cards of forms viewer should be exactly the selected rows that pass the filter
    When user clears the row selection
    Then the "cards" reading of forms viewer should be 2
    And no errors should have been logged

  Scenario: A filter that hides selected rows takes their cards with it
    Given user selects rows where "SEVERITY" is "Critical"
    Then the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user adds a range filter on "AGE" from 40 to 60
    Then fewer than 1000 rows should pass the filter
    And the record cards of forms viewer should show rows "304, 428, 512"
    And the record cards of forms viewer should be exactly the selected rows that pass the filter
    And 5 rows should be selected
    And the "records shown" reading of forms viewer should be lower than before
    When user hovers over "AGE" filter card
    And user clicks on close of "AGE" filter card
    Then 1000 rows should pass the filter
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user clears the row selection
    Then no errors should have been logged

  Scenario: A double-click on the sort label takes the next step of the cycle
    The handler branches on the state it finds: no sort column becomes that column descending, a
    descending column becomes ascending, and an ascending column clears the sort. Each of the three
    is entered deliberately through the properties and claimed on its own rather than chained —
    chaining reads the cycle through whichever state the previous double-click left, which is why
    the spec this replaces refused to claim the order of the three at all.
    Given user selects rows where "SEVERITY" is "Critical"
    When user sets "sortByColumnName" property of forms viewer to "AGE"
    Then the "sort direction" reading of forms viewer should be "↓"
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user double-clicks on the "label AGE" area of forms viewer
    Then the "sort column" reading of forms viewer should be "AGE"
    And the "sort direction" reading of forms viewer should be "↑"
    And the record cards of forms viewer should show rows "215, 430, 428, 512, 304"
    When user sets properties of forms viewer:
      | sortByColumnName | WEIGHT |
      | sortAscending    | true   |
    Then the "sort direction" reading of forms viewer should be "↑"
    And the record cards of forms viewer should show rows "215, 430, 512, 428, 304"
    When user double-clicks on the "label WEIGHT" area of forms viewer
    And the "sort column" reading of forms viewer should be ""
    And forms viewer should not have a "sort indicator WEIGHT" area
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user double-clicks on the "label USUBJID" area of forms viewer
    Then the "sort column" reading of forms viewer should be "USUBJID"
    And the "sort direction" reading of forms viewer should be "↓"
    And forms viewer should have a "sort indicator USUBJID" area
    And the record cards of forms viewer should show rows "512, 430, 428, 304, 215"
    When user sets properties of forms viewer:
      | sortByColumnName |       |
      | sortAscending    | false |
    And user clears the row selection
    Then the "sort column" reading of forms viewer should be ""
    And no errors should have been logged

  Scenario: Sorting the grid mirrors the card order and marks the sorted label
    Given user selects rows where "SEVERITY" is "Critical"
    Then the "sort column" reading of forms viewer should be ""
    And forms viewer should not have a "sort indicator AGE" area
    When user double-clicks on the "header AGE" area of grid
    And user makes row 1 current
    Then the "sort column" reading of grid should be "AGE"
    And the "sort column" reading of forms viewer should be "AGE"
    And the "sort direction" reading of forms viewer should be "↓"
    And forms viewer should have a "sort indicator AGE" area
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user double-clicks on the "header AGE" area of grid
    And user makes row 1 current
    Then the "sort direction" reading of forms viewer should be "↑"
    And the record cards of forms viewer should show rows "215, 430, 428, 512, 304"
    When user double-clicks on the "header AGE" area of grid
    And user makes row 1 current
    Then the "sort column" reading of grid should be ""
    And the "sort column" reading of forms viewer should be ""
    And forms viewer should not have a "sort indicator AGE" area
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user clears the row selection
    Then no errors should have been logged

  Scenario: Sort By overrides the grid's own sort and moves the indicator with it
    Given user selects rows where "SEVERITY" is "Critical"
    When user double-clicks on the "header AGE" area of grid
    And user makes row 1 current
    Then the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user sets "sortByColumnName" property of forms viewer to "WEIGHT"
    Then the "sort column" reading of forms viewer should be "WEIGHT"
    And forms viewer should have a "sort indicator WEIGHT" area
    And forms viewer should not have a "sort indicator AGE" area
    And the record cards of forms viewer should show rows "304, 428, 512, 430, 215"
    And the "sort column" reading of grid should be "AGE"
    When user sets "sortByColumnName" property of forms viewer to ""
    Then the "sort column" reading of forms viewer should be "AGE"
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
    And user makes row 1 current
    And user clears the row selection
    Then the "sort column" reading of grid should be ""
    And no errors should have been logged

  Scenario: Use Grid Sort OFF stops the mirroring — GROK-20380 no longer reproduces
    Given user selects rows where "SEVERITY" is "Critical"
    When user double-clicks on the "header AGE" area of grid
    And user makes row 1 current
    Then the "sort column" reading of forms viewer should be "AGE"
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user sets "useGridSort" property of forms viewer to "false"
    Then the "sort column" reading of forms viewer should be ""
    And the "sort column" reading of grid should be "AGE"
    And forms viewer should not have a "sort indicator AGE" area
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user sets "useGridSort" property of forms viewer to "true"
    Then the "sort column" reading of forms viewer should be "AGE"
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
    When user picks "Sort > Reset" from the context menu of the "header AGE" area of grid
    And user makes row 1 current
    And user clears the row selection
    Then the "sort column" reading of grid should be ""
    And no errors should have been logged

  Scenario: Pin Row moves a card into the pinned pane and takes it out of the record set
    Given user selects rows where "SEVERITY" is "Critical"
    And user sets "showMouseOverRow" property of forms viewer to "false"
    Then the "pinned pane shown" reading of forms viewer should be "false"
    And the "pinned records" reading of forms viewer should be 0
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user picks "Pin Row" from the context menu of the "field USUBJID of card 3" area of forms viewer
    Then the "pinned pane shown" reading of forms viewer should be "true"
    And the "pinned records" reading of forms viewer should be 1
    And the pinned cards of forms viewer should show rows "304"
    And the record cards of forms viewer should show rows "215, 428, 430, 512"
    And the "pinned by" reading of forms viewer should be "USUBJID"
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    And the "records shown" reading of forms viewer should be 6
    And 5 rows should be selected
    And no error or warning balloon should have been shown
    When user picks "Unpin Row" from the context menu of the "field USUBJID of pinned card 1" area of forms viewer
    Then the "pinned pane shown" reading of forms viewer should be "false"
    And the "pinned records" reading of forms viewer should be 0
    And the "pinned values" reading of forms viewer should be ""
    And the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user sets "showMouseOverRow" property of forms viewer to "true"
    And user clears the row selection
    Then no errors should have been logged

  Scenario: Pinning through a non-unique field warns that the layout will not carry it
    Given user selects rows where "SEVERITY" is "Critical"
    And user sets "showMouseOverRow" property of forms viewer to "false"
    When user picks "Pin Row" from the context menu of the "field SEX of card 3" area of forms viewer
    Then a warning balloon containing "You have pinned a non-unique value. It won't be applied from the layout." should have been shown
    And the "pinned by" reading of forms viewer should be "SEX"
    And the "pinned values" reading of forms viewer should be "M"
    And the "pinned records" reading of forms viewer should be 1
    And the "pinned pane shown" reading of forms viewer should be "true"
    When user picks "Unpin Row" from the context menu of the "field SEX of pinned card 1" area of forms viewer
    Then the "pinned records" reading of forms viewer should be 0
    And the "pinned values" reading of forms viewer should be ""
    When user sets "showMouseOverRow" property of forms viewer to "true"
    And user clears the row selection
    Then no errors should have been logged

  @known-failure
  Scenario: A sort from the grid header keeps the cards it was showing
    A grid-header sort clears the table's current row, and `render` then hands the virtual view a
    leading card built for row -1 — a stack of empty divs, zero pixels tall. The view measures its
    first item to size its rows, gets nothing, and lays out no card at all: the five selected rows
    lose their cards too, and they do not come back until some row is made current again. The claim
    below is what the viewer should do. It is the last scenario of the feature because it leaves
    the grid sorted and the viewer blank.
    Given user selects rows where "SEVERITY" is "Critical"
    Then the record cards of forms viewer should show rows "215, 304, 428, 430, 512"
    When user double-clicks on the "header AGE" area of grid
    Then the "sort column" reading of forms viewer should be "AGE"
    And the "cards" reading of forms viewer should be 7
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"
