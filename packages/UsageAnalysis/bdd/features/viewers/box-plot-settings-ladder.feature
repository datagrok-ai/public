@journey @viewers @realizes:viewers.box-plot
Feature: Box plot settings ladder
  Every setting of the ladder takes and survives what comes after it: a datetime value gates
  Axis Type, Category 1 sets the marker color, an explicit coloring survives a value change,
  value min and max, the log axis, a zoom that survives a coloring change, group comparison
  with a control and a covariate; then the whole ladder through a layout round-trip and a
  project round-trip. One journey on demog-1000 with a box plot of AGE.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | Value | AGE |

  Scenario: A datetime value disables Axis Type
    When user sets "Value" property of box plot viewer to "STARTED"
    And user clicks on settings icon of box plot viewer
    Then "Axis Type" property in context panel should be disabled
    When user sets "Value" property of box plot viewer to "AGE"
    Then "Axis Type" property in context panel should be enabled

  Scenario: Category 1 sets the marker color
    Then "Category 1" property of box plot viewer should be "DIS_POP"
    And "Marker Color Column" property of box plot viewer should be "DIS_POP"
    When user sets "Category 1" property of box plot viewer to "SEX"
    Then "Category 1" property of box plot viewer should be "SEX"
    And "Marker Color Column" property of box plot viewer should be "SEX"

  Scenario: An explicit coloring survives a value change
    When user sets "Marker Color Column" property of box plot viewer to "HEIGHT"
    And user sets properties of box plot viewer:
      | Invert Color Scheme | true |
      | Color Min           | 20   |
      | Color Max           | 80   |
    And user sets "Category 2" property of box plot viewer to "RACE"
    And user sets properties of box plot viewer:
      | Show Minor Categories | true |
      | Show All Categories   | true |
    And user sets "Value" property of box plot viewer to "WEIGHT"
    Then properties of box plot viewer should be:
      | Value                 | WEIGHT |
      | Category 1            | SEX    |
      | Category 2            | RACE   |
      | Marker Color Column   | HEIGHT |
      | Invert Color Scheme   | true   |
      | Color Min             | 20     |
      | Color Max             | 80     |
      | Show Minor Categories | true   |
      | Show All Categories   | true   |

  Scenario: Value limits and the log axis
    When user sets properties of box plot viewer:
      | Value Min | 20 |
      | Value Max | 60 |
    Then properties of box plot viewer should be:
      | Value Min | 20 |
      | Value Max | 60 |
    When user sets properties of box plot viewer:
      | Value Min | |
      | Value Max | |
    And user sets "Axis Type" property of box plot viewer to "logarithmic"
    Then "Axis Type" property of box plot viewer should be "logarithmic"
    And no errors should have been logged
    And the value range of box plot viewer should lie within "WEIGHT" column
    When user sets properties of box plot viewer:
      | Invert Y Axis | true   |
      | Plot Style    | violin |
    Then properties of box plot viewer should be:
      | Invert Y Axis | true   |
      | Plot Style    | violin |

  Scenario: A zoom survives a coloring change
    When user zooms into the value axis of box plot viewer
    Then box plot viewer should show a narrower value range than before
    When user sets "Marker Color Column" property of box plot viewer to "SEX"
    Then "Marker Color Column" property of box plot viewer should be "SEX"
    And box plot viewer should show the same value range as before

  Scenario: Group comparison with a control and a covariate
    When user sets "Show Group Comparison" property of box plot viewer to "true"
    And user sets properties of box plot viewer:
      | Control Comparisons | true |
      | Control Group       | F    |
    Then "Control Group" property of box plot viewer should be "F"
    When user sets "Adjust By" property of box plot viewer to "HEIGHT"
    Then "Adjust By" property of box plot viewer should be "HEIGHT"
    And "Adjust by" column input in box plot viewer should contain text "HEIGHT"

  Scenario: The ladder survives a layout round-trip
    When user sets properties of box plot viewer:
      | Adjust By             |       |
      | Control Comparisons   | false |
      | Show Group Comparison | false |
    And user saves the layout of the current table view
    And user clicks on close icon of box plot viewer
    Then box plot viewer should be absent
    When user adds a scatter plot viewer
    Then scatter plot viewer should be visible
    When user loads the saved layout
    Then box plot viewer should be visible
    And scatter plot viewer should be absent
    And properties of box plot viewer should be:
      | Value                 | WEIGHT      |
      | Category 1            | SEX         |
      | Category 2            | RACE        |
      | Show Minor Categories | true        |
      | Show All Categories   | true        |
      | Marker Color Column   | SEX         |
      | Invert Color Scheme   | true        |
      | Color Min             | 20          |
      | Color Max             | 80          |
      | Axis Type             | logarithmic |
      | Invert Y Axis         | true        |
      | Plot Style            | violin      |

  Scenario: The ladder survives a project round-trip
    When user sets properties of box plot viewer:
      | Show Group Comparison | true   |
      | Control Comparisons   | true   |
      | Control Group         | F      |
      | Adjust By             | HEIGHT |
    And user zooms into the value axis of box plot viewer
    Then box plot viewer should show a narrower value range than before
    When user remembers the value range of box plot viewer
    And user saves the current view as project "bdd box plot ladder"
    And user closes all views
    And user opens the "bdd box plot ladder" project
    Then box plot viewer should be visible
    And properties of box plot viewer should be:
      | Value                 | WEIGHT      |
      | Category 1            | SEX         |
      | Category 2            | RACE        |
      | Marker Color Column   | SEX         |
      | Invert Color Scheme   | true        |
      | Color Min             | 20          |
      | Color Max             | 80          |
      | Axis Type             | logarithmic |
      | Invert Y Axis         | true        |
      | Plot Style            | violin      |
      | Show Group Comparison | true        |
      | Control Group         | F           |
      | Adjust By             | HEIGHT      |
    And box plot viewer should show the remembered value range
