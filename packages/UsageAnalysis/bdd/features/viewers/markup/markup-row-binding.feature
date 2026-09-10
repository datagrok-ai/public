@journey @viewers @realizes:viewers.markup
Feature: Markup column references, table expressions and the current row
  The one binding this viewer has: `${COLUMN}` is replaced with the current row's value, and
  `#{...}` is handed to the markup engine before that. Both are read back off the rendered `text`,
  which is the only place the substituted values exist — the `content` property still holds the
  placeholders.
  The old spec reached the current row by clicking the grid 150 by 140 pixels in and then polling
  `currentRowIdx` until it stopped being -1, and read the expected value out of the data frame in
  the same breath, so it compared the viewer with the frame rather than with a number. The rows are
  named here and their values are the fixture's: row 1 of demog-1000 is AGE 26, SEX F and row 4 is
  AGE 45, SEX M.
  `#{t.rowCount}` renders **1000** on demog-1000 — the older specs expect 5850, which was demog.
  The substitution walks the frame's real columns (`markup_viewer_core.dart:117-118`), so a
  placeholder naming a column that does not exist is not "handled": it is simply never matched, and
  stays on screen exactly as written. That is the shape of the loop, and it is asserted next to a
  placeholder that does resolve, so the scenario cannot pass because nothing was substituted at all.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a markup viewer with:
      | content | Age: ${AGE} Sex: ${SEX} |
    And user makes row 1 current
    Then 1000 rows should pass the filter
    And the "current row" reading of markup viewer should be 1
    And markup viewer should report no error

  Scenario: The references render the current row's values and follow it
    Then "AGE" of the current row should be "26"
    And the "text" reading of markup viewer should be "Age: 26 Sex: F"
    When user makes row 4 current
    Then the "current row" reading of markup viewer should be 4
    And "AGE" of the current row should be "45"
    And the "text" reading of markup viewer should be "Age: 45 Sex: M"
    When user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: 26 Sex: F"
    And no errors should have been logged

  Scenario: A reference to a column that does not exist is left exactly as written
    When user sets "content" property of markup viewer to "Known: ${AGE} Unknown: ${NO_SUCH_COLUMN}"
    Then the "text" reading of markup viewer should be "Known: 26 Unknown: ${NO_SUCH_COLUMN}"
    And the "text" reading of markup viewer should not include the text "${AGE}"
    When user sets "content" property of markup viewer to "Age: ${AGE} Sex: ${SEX}"
    Then the "text" reading of markup viewer should be "Age: 26 Sex: F"
    And no errors should have been logged

  Scenario: The table expressions count the rows and the selection, and follow both
    When user sets "content" property of markup viewer to "Rows: #{t.rowCount} Selected: #{t.selection.trueCount}"
    Then the "text" reading of markup viewer should be "Rows: 1000 Selected: 0"
    When user selects all rows
    Then 1000 rows should be selected
    And the "text" reading of markup viewer should be "Rows: 1000 Selected: 1000"
    When user selects rows where "SEX" is "F"
    Then 553 rows should be selected
    And the "text" reading of markup viewer should be "Rows: 1000 Selected: 553"
    When user selects no rows
    Then the "text" reading of markup viewer should be "Rows: 1000 Selected: 0"
    And no errors should have been logged

  Scenario: A heading and a list are rendered around the substituted values
    When user sets "content" property of markup viewer to "# Demographics\n\n* Age: ${AGE}\n* Sex: ${SEX}"
    Then the "mode" reading of markup viewer should be "Markup"
    And the "heading 1" reading of markup viewer should be "Demographics"
    And the "list items" reading of markup viewer should be 2
    And the "text" reading of markup viewer should include the text "Age: 26"
    And the "text" reading of markup viewer should include the text "Sex: F"
    When user sets "content" property of markup viewer to "Age: ${AGE} Sex: ${SEX}"
    Then the "heading 1" reading of markup viewer should be ""
    And no errors should have been logged
