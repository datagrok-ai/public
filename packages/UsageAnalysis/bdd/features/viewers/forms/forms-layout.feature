@journey @viewers @realizes:viewers.forms
Feature: Forms viewer layout and project round-trips
  Steps 7a, 7b and 7c of the forms-core scenario, whose subject is a configured Forms viewer
  surviving the server. The state is a reversed field subset, a viewer-set sort on AGE and one row
  pinned by value — and what a round-trip has to bring back is exactly that: the field set IN ORDER,
  the sort label the indicator sits on, and the pinned row identified BY VALUE, not by row index.
  `pinned by` and `pinned values` are the pair a layout persists (`resolvePinnedRows` maps them back
  to a row), so a claim about pinning across a round-trip compares those, never a row number.
  A layout applied over a corrupted view has a second job: the histogram added after the layout was
  saved is not in it and must go.
  The pin is made by value on USUBJID, which is unique on demog-1000 — a non-unique pin warns that
  it will not survive, and that warning is the forms-core feature's business.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer with:
      | fieldsColumnNames | SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID |
      | sortByColumnName  | AGE                                                                                 |
      | showMouseOverRow  | false                                                                               |
    And user selects rows where "SEVERITY" is "Critical"
    Then forms viewer should be visible
    And the "fields" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "sort column" reading of forms viewer should be "AGE"
    And the record cards of forms viewer should show rows "304, 512, 428, 430, 215"

  Scenario: The field set is drawn in the order it was given, not in table order
    Then the "header labels" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "fields" and "header labels" readings of forms viewer should be the same
    And forms viewer should have a "sort indicator AGE" area
    And no errors should have been logged

  Scenario: A row is pinned by value, and the pinned pane holds it
    Then the "pinned pane shown" reading of forms viewer should be "false"
    When user picks "Pin Row" from the context menu of the "field USUBJID of card 2" area of forms viewer
    Then the "pinned pane shown" reading of forms viewer should be "true"
    And the "pinned by" reading of forms viewer should be "USUBJID"
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    And the pinned cards of forms viewer should show rows "304"
    And the record cards of forms viewer should show rows "512, 428, 430, 215"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Re-applying the saved layout over a corrupted view restores the viewer and drops a foreign one
    Given user remembers the fields of forms viewer
    When user saves the layout of the current table view to the server
    And user clicks on close icon of forms viewer
    Then forms viewer should be absent
    When user adds a histogram viewer
    Then histogram viewer should be visible
    When user loads the saved layout
    Then forms viewer should be visible
    And histogram viewer should be absent
    And the open tableview should have 0 histogram viewers
    And the fields of forms viewer should be as remembered
    And the "fields" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "header labels" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "sort column" reading of forms viewer should be "AGE"
    And forms viewer should have a "sort indicator AGE" area
    And the "pinned by" reading of forms viewer should be "USUBJID"
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    And the "pinned pane shown" reading of forms viewer should be "true"
    And the pinned cards of forms viewer should show rows "304"
    And no errors should have been logged

  Scenario: A project round-trip brings the field set and the pinned row back across a session
    When user saves the current view as project "bdd forms layout"
    And user closes all views
    And user opens the "bdd forms layout" project
    Then forms viewer should be visible
    And the "fields" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "header labels" reading of forms viewer should be "SEVERITY, STARTED, CONTROL, DEMOG, WEIGHT, HEIGHT, DIS_POP, RACE, SEX, AGE, USUBJID"
    And the "sort column" reading of forms viewer should be "AGE"
    And forms viewer should have a "sort indicator AGE" area
    And the "pinned by" reading of forms viewer should be "USUBJID"
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    And the "pinned pane shown" reading of forms viewer should be "true"
    And the pinned cards of forms viewer should show rows "304"
    And no errors should have been logged
