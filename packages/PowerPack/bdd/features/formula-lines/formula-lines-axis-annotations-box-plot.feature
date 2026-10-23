@viewers @realizes:viewers.box-plot @realizes:powerpack.dialogs.formula-lines
Feature: The Annotations group on the box plot's value axis
  The box plot is the one viewer with a numeric axis that `formula-lines-axis-annotations.feature`
  does not walk. Its value axis (Y) carries the Annotations group — Add Line at the column's median,
  Add Band and Add Region between its first and third quartile — and each item opens PowerPack's
  Formula Lines dialog on the new item; its category axis (X) has no such group
  (`context_menu_feature.dart` adds it for the Y axis only). On demog-1000 by RACE: AGE's median is
  45 and its quartiles 36 and 56. The axis is right-clicked at its bottom end, away from the value
  selector. Every scenario empties what it added.
  The line and the band are claimed by what the gesture wrote into the look and by the dialog's
  inputs: the box plot reports no area for a formula line (see `MISSING.md`), while a region is
  claimed by the area it is drawn over.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a box plot viewer with:
      | categoryColumnNames | RACE |
      | valueColumnName     | AGE  |
    And user resizes box plot viewer to 800 by 500

  Scenario: The value axis adds a line, a band and a region on AGE
    When user picks "Annotations > Add Line" from the context menu of the "bottom edge of y axis" area of box plot viewer
    Then "Formula Lines" dialog should be visible
    And editor of Column input in "Formula Lines" dialog should have text "AGE"
    And Value input in "Formula Lines" dialog should have value "45"
    When user clicks OK button in "Formula Lines" dialog
    Then the "Formula Lines" dialog should close
    And "formulaLines" property of box plot viewer should contain "${AGE} = 45.0"
    When user sets "formulaLines" property of box plot viewer to ""
    And user picks "Annotations > Add Band" from the context menu of the "bottom edge of y axis" area of box plot viewer
    And user clicks OK button in "Formula Lines" dialog
    Then "formulaLines" property of box plot viewer should contain "${AGE} in (36.0, 56.0)"
    When user sets "formulaLines" property of box plot viewer to ""
    Then the "viewer regions" reading of box plot viewer should be 0
    When user picks "Annotations > Add Region" from the context menu of the "bottom edge of y axis" area of box plot viewer
    Then "Formula Lines" dialog should be visible
    When user clicks OK button in "Formula Lines" dialog
    Then "annotationRegions" property of box plot viewer should contain "${AGE} = 36.0"
    And "annotationRegions" property of box plot viewer should contain "${AGE} = 56.0"
    And the "viewer regions" reading of box plot viewer should be 1
    And box plot viewer should have a "region 1" area
    And the "region 1" and "view" areas of box plot viewer should be the same width
    When user sets "annotationRegions" property of box plot viewer to "[]"
    Then box plot viewer should not have a "region 1" area
    And no errors should have been logged

  Scenario: The category axis has no Annotations group
    When user right-clicks on the "x axis" area of box plot viewer
    Then the open menu should list "Show Category Axis"
    And the open menu should not list "Annotations"
    And the open menu should not list "Add Line"
    When user closes the context menu
    And user right-clicks on the "y axis" area of box plot viewer
    Then the open menu should list "Annotations"
    When user closes the context menu
    Then no errors should have been logged
