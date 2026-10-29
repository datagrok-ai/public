@viewers @realizes:viewers.forms
Feature: Forms viewer edge cases
  Guards for the Forms viewer (the card viewer PowerGrid registers, not the Dart Form) at the edges of
  its data, on demog-1000. A table with every row deleted leaves no record card and logs no error
  (the viewer has no error state of its own, `libraries/utils/src/viewers/CLAUDE.md`). A pinned row
  that is deleted leaves the pinned pane: pinned rows are found again by their value whenever rows
  are removed, and a value no row holds any more is dropped (`forms-viewer.ts` `resolvePinnedRows`).
  A pinned row the filter hides stays pinned, as in the grid.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a forms viewer
    Then the "fields shown" reading of forms viewer should be 11

  @edge-cases
  Scenario: A table with every row deleted leaves no record card and logs no error
    Given user selects rows where "SEVERITY" is "Critical"
    Then the "records shown" reading of forms viewer should be at least 5
    When user selects all rows
    And user deletes the selected rows
    Then the table should have 0 rows
    And the "records shown" reading of forms viewer should be 0
    And no errors should have been logged

  @edge-cases
  Scenario: A pinned row that is deleted leaves the pinned pane
    Given user selects rows where "SEVERITY" is "Critical"
    And user sets "showMouseOverRow" property of forms viewer to "false"
    When user picks "Pin Row" from the context menu of the "field USUBJID of card 3" area of forms viewer
    Then the "pinned records" reading of forms viewer should be 1
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    When user clears the row selection
    And user selects rows where "USUBJID" is "X0273T29012500105"
    And user deletes the selected rows
    Then the table should have 999 rows
    And the "pinned records" reading of forms viewer should be 0
    And the "pinned pane shown" reading of forms viewer should be "false"
    And the "pinned values" reading of forms viewer should be ""
    And no errors should have been logged

  @edge-cases
  Scenario: A pinned row the filter hides stays pinned
    Given user selects rows where "SEVERITY" is "Critical"
    And user sets "showMouseOverRow" property of forms viewer to "false"
    When user picks "Pin Row" from the context menu of the "field USUBJID of card 3" area of forms viewer
    Then the "pinned values" reading of forms viewer should be "X0273T29012500105"
    When user clears the row selection
    And user filters rows where "SEVERITY" is "Low"
    Then fewer than 1000 rows should pass the filter
    And the "pinned records" reading of forms viewer should be 1
    And the "pinned values" reading of forms viewer should be "X0273T29012500105"
    When user resets the filter
    Then the "pinned records" reading of forms viewer should be 1
    And no errors should have been logged
