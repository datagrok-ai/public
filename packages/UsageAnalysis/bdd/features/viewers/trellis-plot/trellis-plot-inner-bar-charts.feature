@journey @viewers @realizes:viewers.trellis-plot
Feature: Trellis plot — inner bar charts on spgi-100
  Two closed issues about bar charts inside a trellis, each claimed by what the cells draw. On
  spgi-100 split by Stereo Category, the inner bar charts split by Series redraw their cells when
  Include Nulls goes off (the cell's signature changes), and on is the picture they started from
  (public#3617); the trellis reports no area for an inner chart's bar, so the (no value) bar itself
  is not claimed. Bar charts inside a trellis (Scaffold Names by Series, Value = CAST Idea ID, split by
  Stereo Category) used to show part of their data until a bar chart was clicked; every drawn cell now
  has a picture at once, and a full redraw (Bar Sort Order descending and back) gives the same
  picture, so the first one was complete (public#1547). The inner viewer is changed in the trellis's
  own viewer selector, as a person does. One journey on spgi-100; each scenario adds its own trellis,
  and the first closes it.

  Background:
    Given user is logged in
    And user opens spgi dataset

  @realizes:public#3617
  Scenario: Include Nulls off redraws the inner bar charts, and on brings their picture back (public#3617)
    Given user adds a trellis plot viewer with:
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
    When user clicks on close icon of trellis plot viewer
    Then no errors should have been logged

  @realizes:public#1547
  Scenario: Inner bar charts draw complete on the first pass (public#1547)
    Given user adds a trellis plot viewer with:
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
