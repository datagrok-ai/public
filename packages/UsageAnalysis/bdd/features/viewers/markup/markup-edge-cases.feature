@viewers @edge-cases @realizes:viewers.markup
Feature: Markup edge cases — HTML in a cell value, a referenced column removed, every row deleted
  The markup viewer on data its other features do not reach. A cell value that holds HTML tags is
  substituted into `${COLUMN}` as text: the tags are shown as written and nothing is rendered from
  them, since table data is never HTML (`d4/CLAUDE.md` rule 8; the value is `htmlEscape`d in
  `markup_viewer_core.dart`). A column removed from the table while the content references it
  leaves the reference exactly as written, as a reference to a column that does not exist is left
  (`markup-row-binding.feature`). A table whose rows are all deleted leaves the viewer with no
  current row and no error. `text` is what the viewer shows after substitution, `bold` the number
  of bold elements in it.

  Background:
    Given user is logged in

  Scenario: HTML in a substituted cell value is shown as text and renders nothing
    Given user opens a table "markup-html" with:
      | NAME                  | SCORE |
      | <b>bold</b> <i>it</i> | 1     |
    When user adds a markup viewer with:
      | content | Name: ${NAME} |
    And user makes row 1 current
    Then the "current row" reading of markup viewer should be 1
    And the "text" reading of markup viewer should be "Name: <b>bold</b> <i>it</i>"
    And the "bold" reading of markup viewer should be 0
    And markup viewer should report no error
    And no errors should have been logged

  Scenario: A column removed while the content references it is left as written
    Given user opens demog-1000 dataset
    When user adds a markup viewer with:
      | content | Age: ${AGE} Sex: ${SEX} |
    And user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: 26 Sex: F"
    When user removes "SEX" column
    Then the table should not have a column "SEX"
    And the "text" reading of markup viewer should be "Age: 26 Sex: ${SEX}"
    And markup viewer should report no error
    And no errors should have been logged

  Scenario: Once every row of the table is deleted the viewer has no current row and logs no error
    Given user opens demog-1000 dataset
    When user adds a markup viewer with:
      | content | Rows: #{t.rowCount} |
    Then the "text" reading of markup viewer should be "Rows: 1000"
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "text" reading of markup viewer should be "Rows: 0"
    And the "current row" reading of markup viewer should be 0
    And markup viewer should report no error
    And no errors should have been logged

  Scenario Outline: Edit content with new text typed and then <dismissal> leaves the content as it was
    Given user opens demog-1000 dataset
    When user adds a markup viewer with:
      | content | Age: ${AGE} |
    And user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: 26"
    When user picks "Edit content..." from the context menu of markup viewer
    Then Edit dialog should be visible
    When user enters "Weight: ${WEIGHT}" into input in Edit dialog
    And <dismiss>
    Then Edit dialog should be absent
    And "content" property of markup viewer should be "Age: ${AGE}"
    And the "text" reading of markup viewer should be "Age: 26"
    And no errors should have been logged

    Examples:
      | dismissal | dismiss                                       |
      | CANCEL    | user clicks on "CANCEL" button in Edit dialog |
      | Escape    | user presses Escape                           |

  @known-failure:82 @ticket-pending
  Scenario: With the current row filtered out the placeholder shows empty
    Given user opens demog-1000 dataset
    When user adds a markup viewer with:
      | content | Age: [${AGE}] |
    And user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: [26]"
    When user filters rows where "SEX" is "M"
    Then 447 rows should pass the filter
    And the "text" reading of markup viewer should be "Age: []"
    And markup viewer should report no error
    And no errors should have been logged

  @known-failure:100 @ticket-pending
  Scenario: Once the current row is deleted the table has no current row and the placeholder shows empty
    Given user opens a table "markup-no-current" with:
      | NAME  | SCORE |
      | first | 1     |
      | last  | 2     |
    When user adds a markup viewer with:
      | content | Name: [${NAME}] |
    And user makes the last row current
    Then the "text" reading of markup viewer should be "Name: [last]"
    When user selects rows where "NAME" is "last"
    And user deletes the selected rows
    Then the table should have 1 row
    And the "current row" reading of markup viewer should be 0
    And the "text" reading of markup viewer should be "Name: []"
    And markup viewer should report no error
    And no errors should have been logged

  @known-failure:113 @ticket-pending
  Scenario: A referenced column renamed is followed by the content
    Given user opens demog-1000 dataset
    When user adds a markup viewer with:
      | content | Age: ${AGE} |
    And user makes row 1 current
    Then the "text" reading of markup viewer should be "Age: 26"
    When user renames "AGE" column to "YEARS"
    And user makes row 4 current
    Then the "text" reading of markup viewer should be "Age: 45"
    And "content" property of markup viewer should be "Age: ${YEARS}"
    And markup viewer should report no error
    And no errors should have been logged
