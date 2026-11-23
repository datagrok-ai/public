@viewers @realizes:viewers.box-plot @realizes:core#1774
Feature: Box plot p-value under the viewer's own filter
  The p-value and the group comparison follow the viewer's Filter formula, as they follow the
  table's filter (core#1774). A box plot of AGE by RACE on demog-1000, filtered by SEX.

  Scenario: The p-value follows the viewer's Filter
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE  |
      | Category 1 | RACE |
    Then box plot viewer should have a "p value" area
    When user remembers the "p value" reading of box plot viewer
    And user sets "Filter" property of box plot viewer to "${SEX} == \"F\""
    Then the "p value" reading of box plot viewer should not be as remembered
    When user sets "Show Group Comparison" property of box plot viewer to "true"
    Then box plot viewer should have a "group comparison" area
    When user remembers the "p value" reading of box plot viewer
    And user sets "Filter" property of box plot viewer to "${SEX} == \"M\""
    Then the "p value" reading of box plot viewer should not be as remembered
    When user remembers the "p value" reading of box plot viewer
    And user sets "Filter" property of box plot viewer to ""
    Then the "p value" reading of box plot viewer should not be as remembered
    And no errors should have been logged
