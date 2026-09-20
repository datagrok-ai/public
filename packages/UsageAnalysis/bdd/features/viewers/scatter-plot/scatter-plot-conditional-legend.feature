@viewers @realizes:viewers.scatter-plot
Feature: Scatter plot legend of a conditionally coloured column
  `scatterplot-legend.md`, scenario 5: a column whose colour coding is a set of rules, set as the
  plot's Color, gives a legend just as a categorical column does, with one entry per rule, named as
  the rule. The md's spgi-100 identifier column becomes demog-1000's AGE (18 to 89, no blanks), cut
  into two ranges that both occur; the colour column is set as a property, not in the on-chart
  selector. A file of its own, next to `scatter-plot-legend.feature`, which
  another round owns.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then legend of scatter plot viewer should be hidden

  Scenario: A conditionally coloured column set as Color lists one legend entry per rule
    When user colors "AGE" column conditionally:
      | 18-45 | #00FF00 |
      | 45-89 | #FF0000 |
    Then "AGE" column should be color-coded conditionally
    When user sets "Color" property of scatter plot viewer to "AGE"
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 2 items
    And "18-45" legend item in legend of scatter plot viewer should be visible
    And "45-89" legend item in legend of scatter plot viewer should be visible
    And the "18-45" item in the legend of scatter plot viewer should be colored "#00FF00"
    And the "45-89" item in the legend of scatter plot viewer should be colored "#FF0000"
    When user sets "Color" property of scatter plot viewer to ""
    And user removes the coloring of "AGE" column
    Then legend of scatter plot viewer should be hidden
    And no errors should have been logged
