@journey @viewers @realizes:viewers.forms
Feature: Forms viewer layout and project persistence
  What the viewer carries across a round-trip: a field set that is neither the default set nor the
  table's order, the Sort By column, and a pinned row remembered by its value — re-applied over a
  view that had the Forms viewer closed and a histogram added in its place, and again after the
  view has been saved as a project and reopened. One journey on demog-1000. The sort direction is
  not claimed: GROK-20666 is open.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer with:
      | Show Mouse Over Row | false |
    And user sets properties of forms viewer:
      | Fields  | AGE, USUBJID, HEIGHT, SEX |
      | Sort By | AGE                       |
    And user selects rows where "USUBJID" is one of "X0273T21000400001, X0273T21000500008, X0273T21001500015"
    And user picks "Pin Row" from the context menu of the "field USUBJID of card 2" area of forms viewer
    Then the "fields" reading of forms viewer should be "AGE, USUBJID, HEIGHT, SEX"
    And the "sort column" reading of forms viewer should be "AGE"
    And the "pinned records" reading of forms viewer should be 1
    And the "USUBJID of pinned card 1" reading of forms viewer should be "X0273T21000400001"

  Scenario: A saved layout restores the viewer over a changed view
    When user saves the layout of the current table view to the server
    And user clicks on close icon of forms viewer
    Then forms viewer should be absent
    When user adds a histogram viewer
    And user loads the saved layout
    Then histogram viewer should be absent
    And forms viewer should be visible
    And the "fields" reading of forms viewer should be "AGE, USUBJID, HEIGHT, SEX"
    And "Fields" property of forms viewer should be "AGE, USUBJID, HEIGHT, SEX"
    And the "sort column" reading of forms viewer should be "AGE"
    And the "pinned records" reading of forms viewer should be 1
    And the "USUBJID of pinned card 1" reading of forms viewer should be "X0273T21000400001"
    And "pinnedRowValues" property of forms viewer should be "X0273T21000400001"
    And no errors should have been logged

  Scenario: A project round-trip keeps the field set, the sort column and the pinned row
    When user saves the current view as project "zz-forms-persistence"
    And user closes all views
    And user opens the "zz-forms-persistence" project
    Then forms viewer should be visible
    And the "fields" reading of forms viewer should be "AGE, USUBJID, HEIGHT, SEX"
    And the "sort column" reading of forms viewer should be "AGE"
    And the "pinned records" reading of forms viewer should be 1
    And the "USUBJID of pinned card 1" reading of forms viewer should be "X0273T21000400001"
    And no errors should have been logged
