@journey @viewers @realizes:viewers.box-plot
Feature: Box plot statistics and coloring
  The box coloring baseline and an explicit whisker color, the statistics strip and its ladder,
  the statistics format, the p-value toggle by key and by menu, the three-group test branch, the
  violin style with its bins and line widths, column color coding driving the marker colors, and
  a datetime value. One journey on demog-1000 with a box plot of AGE by SEX.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value      | AGE |
      | Category 1 | SEX |

  Scenario: Box coloring
    Then "Whisker Color" property of box plot viewer should be ""
    And box plot viewer should be painted in at least 2 colors
    When user sets "Whisker Color" property of box plot viewer to "#1F77B4"
    Then "Whisker Color" property of box plot viewer should be "#1F77B4"
    And box plot viewer should have repainted
    When user sets "Whisker Color" property of box plot viewer to ""
    Then box plot viewer should have repainted

  Scenario: The statistics strip and its ladder
    Then "Show Statistics" property of box plot viewer should be "true"
    And box plot viewer should have a "stats" area
    And close-stats icon in box plot viewer should be present
    And "Show Avg" property of box plot viewer should be "true"
    And "Show Med" property of box plot viewer should be "true"
    When user sets properties of box plot viewer:
      | Show Total Count    | true |
      | Show Inliers Count  | true |
      | Show Outliers Count | true |
      | Show Stdev          | true |
      | Show Q1             | true |
      | Show Q3             | true |
    Then box plot viewer should have repainted
    And properties of box plot viewer should be:
      | Show Total Count    | true |
      | Show Inliers Count  | true |
      | Show Outliers Count | true |
      | Show Stdev          | true |
      | Show Q1             | true |
      | Show Q3             | true |
    And no errors should have been logged
    When user sets "Statistics Format" property of box plot viewer to "#,##0.00"
    Then "Statistics Format" property of box plot viewer should be "#,##0.00"
    And "Show P Value" property of box plot viewer should be "true"
    And no errors should have been logged
    When user sets "Statistics Format" property of box plot viewer to "auto"
    And user picks "Show Total Count" from the context menu of the "stats" area of box plot viewer
    Then "Show Total Count" property of box plot viewer should be "false"
    When user picks "Show Total Count" from the context menu of the "stats" area of box plot viewer
    Then "Show Total Count" property of box plot viewer should be "true"
    When user sets properties of box plot viewer:
      | Show Total Count    | false |
      | Show Inliers Count  | false |
      | Show Outliers Count | false |
      | Show Stdev          | false |
      | Show Q1             | false |
      | Show Q3             | false |

  Scenario: The T key toggles the p-value
    When user sets "Show P Value" property of box plot viewer to "false"
    And user clicks on empty plot space of box plot viewer
    And user presses t
    Then "Show P Value" property of box plot viewer should be "true"
    When user presses t
    Then "Show P Value" property of box plot viewer should be "false"
    When user sets "Show P Value" property of box plot viewer to "true"

  Scenario: Three groups take the Alexander-Govern branch
    When user sets "Category 1" property of box plot viewer to "RACE"
    And user sets "Marker Color Column" property of box plot viewer to "SEX"
    Then box plot viewer should have a "p value" area
    And box plot viewer should have a "category Caucasian" area
    When user hovers over the "p value" area of box plot viewer
    Then show-group-stats icon in box plot viewer should be visible
    When user moves the pointer away from box plot viewer
    Then show-group-stats icon in box plot viewer should be hidden
    When user sets properties of box plot viewer:
      | Category 1          | SEX |
      | Marker Color Column |     |

  Scenario: The violin style
    When user sets "Plot Style" property of box plot viewer to "violin"
    Then box plot viewer should have repainted
    And the "M values" area of box plot viewer should be painted
    And the "F values" area of box plot viewer should be painted
    When user sets "Bins" property of box plot viewer to "50"
    And user sets "Bins" property of box plot viewer to "500"
    Then box plot viewer should have repainted
    When user sets "Interquartile Line Width" property of box plot viewer to "10"
    Then box plot viewer should have repainted
    When user sets "Violin Line Width" property of box plot viewer to "4"
    Then box plot viewer should have repainted
    When user sets "Violin Whisker Color" property of box plot viewer to "#00AA00"
    Then box plot viewer should have repainted
    When user sets properties of box plot viewer:
      | Plot Style               | box |
      | Bins                     | 100 |
      | Interquartile Line Width | 6   |
      | Violin Line Width        | 2   |
    Then box plot viewer should have repainted

  Scenario: Column color coding drives the marker colors
    When user colors "WEIGHT" column linearly from "#0000FF" to "#FF0000"
    And user sets "Marker Color Column" property of box plot viewer to "WEIGHT"
    Then "Marker Color Column" property of box plot viewer should be "WEIGHT"
    When user colors "WEIGHT" column linearly from "#0000FF" to "#FF0000" over 60 to 120
    Then box plot viewer should have repainted
    When user clears the row selection
    And user drags a selection box over the "M values" area of box plot viewer
    Then some rows should be selected
    And box plot viewer should show more selection highlight than before
    When user colors "WEIGHT" column conditionally:
      | 50-90  | #00FF00 |
      | 90-150 | #800080 |
    Then box plot viewer should have repainted
    And box plot viewer should show a selection highlight
    When user clears the row selection
    And user removes the coloring of "WEIGHT" column
    And user sets "Marker Color Column" property of box plot viewer to "SEX"
    And user colors "SEX" column categorically:
      | M | #E41A1C |
    Then box plot viewer should have repainted
    When user removes the coloring of "SEX" column
    And user sets "Marker Color Column" property of box plot viewer to ""

  Scenario: A datetime value
    When user sets "Value" property of box plot viewer to "STARTED"
    Then "Value" property of box plot viewer should be "STARTED"
    And box plot viewer should be painted
    And no errors should have been logged
    When user sets "Statistics Format" property of box plot viewer to "yyyy-MM-dd HH:mm"
    Then "Statistics Format" property of box plot viewer should be "yyyy-MM-dd HH:mm"
    And no errors should have been logged
    When user sets properties of box plot viewer:
      | Value             | AGE  |
      | Statistics Format | auto |
