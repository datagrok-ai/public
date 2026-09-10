@journey @viewers @realizes:viewers.forms
Feature: Forms viewer colour coding, alignment and font
  What a field looks like: the background it takes from the column's colour coding, the alignment it
  takes from the grid column's style, and the font the viewer gives it. All three are readings —
  `background of <COL> of <label>` as `#RRGGBB`, `align of <COL> of <label>` and
  `font of <COL> of <label>` — where the spec this replaces read `getComputedStyle` in the page and
  normalized two colour spellings against each other to compare them.
  The colour claims are made between two fields rather than against a colour constant: AGE is
  colour-coded and USUBJID is not, so with Color Code on their backgrounds must differ and with it
  off they must agree. That says what the property does without hard-coding a scheme's arithmetic.
  Alignment is read as the stylesheet gives it: the viewer writes `textAlign` only for a `center`
  or `right` grid column style, so a plain column reads `start`, not `left`, and reading `left`
  would mean the viewer had started writing it.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user colors "AGE" column linearly from "#FF0000" to "#00FF00"
    And user adds a forms viewer
    And user makes row 1 current
    Then forms viewer should be visible
    And "AGE" column should be color-coded linearly
    And the "AGE of card 1" reading of forms viewer should be "26"

  Scenario: A colour-coded column paints the field's background, an uncoded one does not
    Then "colorCode" property of forms viewer should be "true"
    And the "background of USUBJID of card 1" reading of forms viewer should be "#FFFFFF"
    And the "background of AGE of card 1" reading of forms viewer should not be "#FFFFFF"
    And the "background of AGE of card 1" and "background of USUBJID of card 1" readings of forms viewer should differ
    And no errors should have been logged

  Scenario: Color Code off drops the background to the uncoded one, and on restores it
    When user sets "colorCode" property of forms viewer to "false"
    Then the "background of AGE of card 1" and "background of USUBJID of card 1" readings of forms viewer should be the same
    And the "background of AGE of card 1" reading of forms viewer should be "#FFFFFF"
    When user sets "colorCode" property of forms viewer to "true"
    Then the "background of AGE of card 1" and "background of USUBJID of card 1" readings of forms viewer should differ
    And the "background of AGE of card 1" reading of forms viewer should not be "#FFFFFF"
    And no errors should have been logged

  Scenario: A row of another colour gives its field another background
    Given user remembers the "background of AGE of card 1" reading of forms viewer
    When user makes row 78 current
    Then the "AGE of card 1" reading of forms viewer should be "60"
    And the "background of AGE of card 1" reading of forms viewer should not be as remembered
    And the "background of AGE of card 1" reading of forms viewer should not be "#FFFFFF"
    When user makes row 1 current
    Then the "background of AGE of card 1" reading of forms viewer should be as remembered
    And no errors should have been logged

  Scenario: Removing the colouring leaves every field the same background
    When user removes the coloring of "AGE" column
    Then "AGE" column should have no color coding
    And the "background of AGE of card 1" and "background of USUBJID of card 1" readings of forms viewer should be the same
    When user colors "AGE" column linearly from "#FF0000" to "#00FF00"
    Then the "background of AGE of card 1" and "background of USUBJID of card 1" readings of forms viewer should differ
    And no errors should have been logged

  Scenario: The Font property reaches every field, and the alignment is the stylesheet's
    Then the "align of AGE of card 1" reading of forms viewer should be "start"
    And the "align of SEX of card 1" reading of forms viewer should be "start"
    And the "font of AGE of card 1" reading of forms viewer should be "13px Roboto"
    When user sets "font" property of forms viewer to "italic bold 15px \"Times New Roman\""
    Then the "font of AGE of card 1" reading of forms viewer should be "italic bold 15px \"Times New Roman\""
    And the "font of SEX of card 1" reading of forms viewer should be "italic bold 15px \"Times New Roman\""
    And the "align of AGE of card 1" reading of forms viewer should be "start"
    When user sets "font" property of forms viewer to "normal normal 13px \"Roboto\""
    Then the "font of AGE of card 1" reading of forms viewer should be "13px Roboto"
    And no errors should have been logged
