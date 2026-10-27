@viewers @realizes:viewers.trellis-plot
Feature: Trellis plot regressions of closed issues
  Closed trellis issues, each claimed by what the cells draw. public#2144: an inner viewer's own
  setting survives switching the trellis to another inner viewer and back; the bar charts set to
  avg come back with the setting and with the very picture they had (the cell's signature as remembered).
  public#3617: on spgi-100 split by Stereo Category, the inner bar charts split by Series redraw
  their cells when Include Nulls goes off (the cell's signature changes), and on is the picture they
  started from; the trellis reports no area for an inner chart's bar, so the (no value) bar itself is
  not claimed. public#1547: bar charts inside a trellis (Scaffold Names by Series, Value = CAST Idea ID,
  split by Stereo Category) used to show part of their data until a bar chart was clicked; every
  drawn cell now has a picture at once, and a full redraw (Bar Sort Order descending and back)
  gives the same picture, so the first one was complete. The inner viewer is changed in the trellis's own viewer
  selector, as a person does.

  Background:
    Given user is logged in

  @realizes:public#2144
  Scenario: An inner bar chart keeps its own setting through a pie chart and back (public#2144)
    Given user opens demog-1000 dataset
    And user adds a trellis plot viewer with:
      | X Column Names | SEX  |
      | Y Column Names | RACE |
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user sets "valueAggrType" inner property of trellis plot viewer to "avg"
    Then the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    When user remembers the "cell signature F | Caucasian" reading of trellis plot viewer
    And user picks "Pie chart" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Pie chart"
    And the "cell signature F | Caucasian" reading of trellis plot viewer should not be as remembered
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    Then the "inner viewer type" reading of trellis plot viewer should be "Bar chart"
    And "valueAggrType" inner property of trellis plot viewer should be "avg"
    And the "cell signature F | Caucasian" reading of trellis plot viewer should be as remembered
    And no errors should have been logged

  @realizes:public#3617
  Scenario: Include Nulls off redraws the inner bar charts, and on brings their picture back (public#3617)
    Given user opens spgi dataset
    And user adds a trellis plot viewer with:
      | Y Column Names | Stereo Category |
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    And user sets "splitColumnName" inner property of trellis plot viewer to "Series"
    Then the "cells drawn" reading of trellis plot viewer should be 15
    When user remembers the "cell signature Chemist 13 | R_ONE" reading of trellis plot viewer
    And user sets "includeNulls" inner property of trellis plot viewer to "false"
    Then the "cell signature Chemist 13 | R_ONE" reading of trellis plot viewer should not be as remembered
    And the "cells drawn" reading of trellis plot viewer should be 15
    When user sets "includeNulls" inner property of trellis plot viewer to "true"
    Then the "cell signature Chemist 13 | R_ONE" reading of trellis plot viewer should be as remembered
    And no errors should have been logged

  @realizes:public#1547
  Scenario: Inner bar charts draw complete on the first pass (public#1547)
    Given user opens spgi dataset
    And user adds a trellis plot viewer with:
      | X Column Names | Scaffold Names |
      | Y Column Names | Series         |
    When user picks "Bar chart" in the viewer selector of trellis plot viewer
    And user sets "valueColumnName" inner property of trellis plot viewer to "CAST Idea ID"
    And user sets "splitColumnName" inner property of trellis plot viewer to "Stereo Category"
    Then the "cells drawn" reading of trellis plot viewer should be 6
    And the "blank cells" reading of trellis plot viewer should be 0
    And the "distinct cell signatures" reading of trellis plot viewer should be 6
    When user remembers the "cell signature TRISUBSTITUTED | Triazoles" reading of trellis plot viewer
    And user sets "barSortOrder" inner property of trellis plot viewer to "desc"
    Then the "cell signature TRISUBSTITUTED | Triazoles" reading of trellis plot viewer should not be as remembered
    When user sets "barSortOrder" inner property of trellis plot viewer to "asc"
    Then the "cell signature TRISUBSTITUTED | Triazoles" reading of trellis plot viewer should be as remembered
    And the "blank cells" reading of trellis plot viewer should be 0
    And no errors should have been logged
