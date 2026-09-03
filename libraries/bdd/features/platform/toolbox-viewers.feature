@platform @viewers @realizes:viewers.scatter-plot
Feature: Adding viewers from the toolbox
  The original example: platform vocabulary over the Dart shell's `name=` conventions. Runs in
  the `platform` Playwright project against DATAGROK_URL (a stand with admin/admin by default).

  Background:
    Given user is logged in
    And user opens spgi dataset

  Scenario: Scatter plot from the toolbox icon
    When user opens toolbox
    And user clicks on scatter plot icon on toolbox
    Then scatter plot viewer should be added to the open tableview

  Scenario Outline: Other viewers the same way
    When user opens toolbox
    And user clicks on <viewer> icon on toolbox
    Then <viewer> viewer should be added to the open tableview

    Examples:
      | viewer    |
      | histogram |
      | bar chart |
