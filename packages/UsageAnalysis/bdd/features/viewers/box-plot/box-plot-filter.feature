@journey @viewers @realizes:viewers.box-plot
Feature: Box plot filter semantics
  How the box plot answers the table's filter and its own: the value range follows the filter
  (Zoom Values By Filter), the viewer's formula filter leaves the table's filter alone, Show
  Empty Categories drops and restores an empty-valued category, and a coloring under a filter
  keeps its color scale. One journey on spgi-100 with a box plot of Average Mass by Series.

  Background:
    Given user is logged in
    And user opens spgi dataset
    And user adds a box plot viewer with:
      | Value      | Average Mass |
      | Category 1 | Series       |
    Then the table should have 100 rows

  Scenario: The value range follows the filter
    Then all rows should pass the filter
    And "Zoom Values By Filter" property of box plot viewer should be "true"
    And box plot viewer should show 100 rows
    When user filters rows where "Average Mass" is between 300 and 400
    Then fewer than 100 rows should pass the filter
    And box plot viewer should show fewer rows than before
    And box plot viewer should show a narrower value range than before
    When user resets the filter
    Then all rows should pass the filter
    And box plot viewer should show a wider value range than before

  Scenario: Zoom Values By Filter off keeps the range
    When user sets "Zoom Values By Filter" property of box plot viewer to "false"
    And user filters rows where "Average Mass" is between 300 and 400
    Then fewer than 100 rows should pass the filter
    And box plot viewer should show the same value range as before
    When user sets "Zoom Values By Filter" property of box plot viewer to "true"
    Then box plot viewer should show a narrower value range than before
    When user resets the filter

  Scenario: The viewer's own filter leaves the table's alone
    When user filters rows where "Average Mass" is between 300 and 400
    And user sets "Filter" property of box plot viewer to "${Average Mass} > 350"
    Then box plot viewer should show fewer rows than before
    And the filter should pass exactly the rows where "Average Mass" is between 300 and 400
    When user sets "Filter" property of box plot viewer to ""
    And user resets the filter

  Scenario: Show Empty Categories drops and restores an empty-valued category
    When user adds a calculated column "AverageMassFixture" with formula "if(${Series} == \"Triazoles\", null, ${Average Mass})"
    Then table "spgi-100" should have missing values in "AverageMassFixture" column
    When user sets "Value" property of box plot viewer to "AverageMassFixture"
    And user sets "Show Empty Categories" property of box plot viewer to "true"
    Then box plot viewer should have a "category Triazoles" area
    When user sets "Show Empty Categories" property of box plot viewer to "false"
    Then box plot viewer should have repainted
    And box plot viewer should not have a "category Triazoles" area
    When user sets "Show Empty Categories" property of box plot viewer to "true"
    Then box plot viewer should have repainted
    And box plot viewer should have a "category Triazoles" area
    When user sets "Value" property of box plot viewer to "Average Mass"
    And user removes "AverageMassFixture" column

  Scenario: A coloring under a filter keeps its color scale
    When user filters rows where "Average Mass" is between 300 and 400
    And user sets "Marker Color Column" property of box plot viewer to "TPSA"
    Then "Marker Color Column" property of box plot viewer should be "TPSA"
    And box plot viewer should have a "color scale" area
    And no errors should have been logged
    When user resets the filter
    Then the color scale of box plot viewer should cover a wider range than before
    When user filters rows where "Average Mass" is between 300 and 400
    Then the color scale of box plot viewer should cover a narrower range than before
    When user resets the filter
    And user sets "Marker Color Column" property of box plot viewer to ""
