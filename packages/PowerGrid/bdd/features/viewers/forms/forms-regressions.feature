@viewers @realizes:viewers.forms
Feature: Forms viewer — regressions of closed issues
  Closed issues of the Forms viewer, each a claim that fails if the bug comes back. An integer cell
  with no value is stored as the int null (-2147483648), and the forms viewer used to print that
  number in the card instead of an empty field (public#3561). The table is a three-row inline one
  whose `n` column the inline table types as int, so the empty value is a real int null, not a
  string; the card of the row before it shows its value first, so the empty field is read on a card
  that has been rebuilt for the empty row. The same claim on the core Form viewer is that viewer's
  own (UsageAnalysis `viewers/form`).

  Scenario: An empty integer cell shows an empty field on the card (public#3561)
    Given user is logged in
    And user opens a table "ints" with:
      | id | n |
      | 1  | 5 |
      | 2  |   |
      | 3  | 7 |
    Then "n" column should have type "int"
    When user adds a forms viewer
    Then forms viewer should be visible
    When user makes row 1 current
    Then the "n of current card" reading of forms viewer should be "5"
    When user makes row 2 current
    Then the "current record" reading of forms viewer should be 2
    And the "id of current card" reading of forms viewer should be "2"
    And the "n of current card" reading of forms viewer should be ""
    When user makes row 3 current
    Then the "n of current card" reading of forms viewer should be "7"
    And no errors should have been logged
