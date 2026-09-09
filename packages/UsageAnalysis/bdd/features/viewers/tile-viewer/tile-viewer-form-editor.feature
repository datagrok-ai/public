@journey @viewers @realizes:viewers.tile-viewer
Feature: Tile viewer form designer
  The card is a sketch form, and `Edit Form...` opens it in a designer of its own: a value host and
  a caption host per column, an EDIT button that opens the column chooser, a RESET that reverts the
  session's edits, and CLOSE AND APPLY that writes the form back and moves the viewer out of its
  auto-generated state for good. The two host sets are read apart, because a summed host count
  would hide one channel changing while the other did not — RESET putting a deleted caption back is
  the whole point of the scenario that checks it, and it reverts to the state the designer OPENED
  in, not to the factory default. Then the contrast the state pair exists for: an auto-generated
  card is rebuilt when a column leaves the table and a column that had no field takes the freed
  slot, while a designed card keeps exactly the fields it was given — the departed column's field
  stays on it, empty.
  One journey on demog-1000 — 11 columns, a card of ten fields, SEVERITY the column the relevance
  score leaves over.

  The lanes column is set in the Background because the viewer's own menu — the only way to
  `Edit Form...` — needs a region no card covers, and in single-lane mode the cards cover the
  whole lane. Not translated: `sketchState['table']` as the designer's target, replaced by the
  `table` reading; and the synthetic six-column fixture the old spec built for the refill
  contrast, which demog-1000 makes unnecessary (ten fields over eleven columns leaves exactly one
  spare). A caption whose value field was deleted stays on the card but the viewer stops reporting
  a `label` region for it, so the caption channel is claimed on the designer, where both host sets
  are visible, and not on the card.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a tile viewer with:
      | Lanes Column Name | SEX |
    Then tile viewer should be visible
    And the "lane names" reading of tile viewer should be "F, M"
    And the "fields shown" reading of tile viewer should be 10

  Scenario: The card the designer opens on is the auto-generated one
    Then the "auto generate" reading of tile viewer should be "true"
    And the "form designed" reading of tile viewer should be "false"
    And the "table" reading of tile viewer should be "demog-1000"
    And the table should have 11 columns
    And the "fields" reading of tile viewer should contain "AGE"
    And the "fields" reading of tile viewer should contain "DEMOG"
    And the "fields" reading of tile viewer should not contain "SEVERITY"
    And form designer should be absent
    And no errors should have been logged

  Scenario: Edit Form opens the designer on the viewer's own table
    Given user listens for "d4-tile-viewer-form-edit-request" event on tile viewer
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    And the form designer should show value fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    And the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And "d4-tile-viewer-form-edit-request" event should have fired on tile viewer
    And the "fields shown" reading of tile viewer should be 10
    And the "AGE of row 1" reading of tile viewer should be "26"
    And no errors should have been logged

  Scenario: The column chooser counts the fields the card shows, and CANCEL keeps them
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    When user clicks on "EDIT" button
    Then "Select columns..." dialog should be visible
    And "All" link in "Select columns..." dialog should be visible
    And "None" link in "Select columns..." dialog should be visible
    And "10 checked" text in "Select columns..." dialog should be visible
    When user clicks on "CANCEL" button in "Select columns..." dialog
    Then "Select columns..." dialog should be absent
    And the form designer should show value fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And the "fields shown" reading of tile viewer should be 10
    And no errors should have been logged

  Scenario: A deleted caption comes back with RESET, and goes with CLOSE AND APPLY
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    And the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user deletes the "SEX" label field in the form designer
    Then the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"
    And the form designer should show value fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user clicks on "RESET" button
    Then the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    And the form designer should show value fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user deletes the "SEX" label field in the form designer
    And user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And the "fields shown" reading of tile viewer should be 10
    And tile viewer should have a "field SEX of row 1" area
    And tile viewer should not have a "label SEX of row 1" area
    And the "SEX of row 1" reading of tile viewer should be "F"
    And the "auto generate" reading of tile viewer should be "false"
    And the "form designed" reading of tile viewer should be "true"
    And no errors should have been logged

  Scenario: RESET reverts to the state the designer opened in, not to the factory card
    When user picks "Edit Form..." from the viewer menu of tile viewer
    Then form designer should be visible
    And the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"
    When user deletes the "AGE" value field in the form designer
    Then the form designer should show value fields "CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    When user clicks on "RESET" button
    Then the form designer should show value fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, SEX, STARTED, USUBJID, WEIGHT"
    And the form designer should show label fields "AGE, CONTROL, DEMOG, DIS_POP, HEIGHT, RACE, STARTED, USUBJID, WEIGHT"
    When user deletes the "AGE" value field in the form designer
    And user clicks on "CLOSE AND APPLY" button
    Then form designer should be absent
    And the "fields shown" reading of tile viewer should be 9
    And the "fields" reading of tile viewer should not contain "AGE"
    And tile viewer should not have a "field AGE of row 1" area
    And the "form designed" reading of tile viewer should be "true"
    And no errors should have been logged

  Scenario: On a designed card a departing column empties its field and no column takes the slot
    Then the "form designed" reading of tile viewer should be "true"
    And the "fields shown" reading of tile viewer should be 9
    And the "fields" reading of tile viewer should contain "DEMOG"
    And the "fields" reading of tile viewer should not contain "SEVERITY"
    And the "DEMOG of row 1" reading of tile viewer should be "26 C F"
    And the table should have 11 columns
    When user remembers the fields of tile viewer
    And user picks "Remove" from the context menu of the "field DEMOG of row 1" area of tile viewer
    Then the table should not have a column "DEMOG"
    And the table should have 10 columns
    And the "fields shown" reading of tile viewer should be 9
    And the fields of tile viewer should be as remembered
    And the "fields" reading of tile viewer should not contain "SEVERITY"
    And tile viewer should have a "field DEMOG of row 1" area
    And the "DEMOG of row 1" reading of tile viewer should be ""
    And the "auto generate" reading of tile viewer should be "false"
    And no errors should have been logged

  Scenario: Auto Generate rebuilds the card, and then a departing column's slot is refilled
    Given user adds a calculated column "DEMOG" with formula "${AGE}"
    Then the table should have 11 columns
    When user sets "Auto Generate" property of tile viewer to "true"
    And user adds a calculated column "TMP" with formula "1"
    And user removes "TMP" column
    Then the table should have 11 columns
    And the "auto generate" reading of tile viewer should be "true"
    And the "form designed" reading of tile viewer should be "false"
    And the "fields shown" reading of tile viewer should be 10
    And the "fields" reading of tile viewer should contain "AGE"
    When user remembers the fields of tile viewer
    And user removes "AGE" column
    Then the table should not have a column "AGE"
    And the table should have 10 columns
    And the "fields shown" reading of tile viewer should be 10
    And the fields of tile viewer should have refilled the freed slot
    And the "fields" reading of tile viewer should not contain "AGE"
    And tile viewer should not have a "field AGE of row 1" area
    And the "auto generate" reading of tile viewer should be "true"
    And no errors should have been logged
    When user adds a calculated column "AGE" with formula "${HEIGHT}"
    Then the table should have 11 columns
