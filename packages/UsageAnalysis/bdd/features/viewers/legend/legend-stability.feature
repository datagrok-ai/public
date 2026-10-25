@journey @viewers @realizes:viewers.legend
Feature: A legend that stays where it is put
  What must not move a legend: the scatter plot's viewport. A wheel zoom and a plain-drag pan
  change what the plot shows (the x axis span narrows, the x axis start moves) while the legend,
  docked or in a corner, keeps its place, its size and its items. A corner legend folded into the
  mini icon by its chevron stays folded while Position Auto and three resizes make the placement
  re-resolve, and only the icon's click opens it again. A legend column with hundreds of categories
  is docked on the right and lists every category of the rows the plot draws, and Position Auto
  docks it too, where the four RACE categories go into a corner.
  From the `playwright-public/legend` specs "panning and zooming never move the legend",
  "a collapsed legend stays collapsed through any auto fallback", "hundreds of categories stay
  virtualised and scrollable" and "a large category count never enters corner mode". One journey
  on demog-1000 with a scatter plot of HEIGHT over WEIGHT coloured by RACE, Zoom and Filter at
  "no action" so a zoom leaves the legend's rows alone; 872 of the 1000 rows have a HEIGHT (499 of
  them F), and USUBJID holds a different value in every row.
  Not translated: that the long list keeps only a window of its items in the page and moves it on a
  scroll — no step scrolls a list to an item that is not in the page yet, and how many items the
  page holds is the list's implementation, not something a person sees.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X                 | WEIGHT    |
      | Y                 | HEIGHT    |
      | Color             | RACE      |
      | Zoom and Filter   | no action |
      | Legend Visibility | Always    |
      | Legend Position   | Right     |
    Then the legend of scatter plot viewer should be docked
    And the legend of scatter plot viewer should list 4 items

  Scenario: A wheel zoom and a pan leave the docked legend where it was
    When user takes a snapshot of scatter plot viewer
    And user scrolls the mouse wheel up 3 times over the "view" area of scatter plot viewer
    Then the "x axis span" reading of scatter plot viewer should be lower than before
    And the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should list the same items as before
    And the legend of scatter plot viewer should be in the "right" slot
    When user takes a snapshot of scatter plot viewer
    And user drags across the "view" area of scatter plot viewer
    Then the "x axis min" reading of scatter plot viewer should differ from before
    And the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should list the same items as before
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then no errors should have been logged

  Scenario: A wheel zoom and a pan leave a corner legend in its corner
    When user sets "Legend Position" property of scatter plot viewer to "RightBottom"
    Then the legend of scatter plot viewer should be in a corner
    When user takes a snapshot of scatter plot viewer
    And user scrolls the mouse wheel up 3 times over the "view" area of scatter plot viewer
    Then the "x axis span" reading of scatter plot viewer should be lower than before
    And the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should be in the "rightBottom" slot
    When user takes a snapshot of scatter plot viewer
    And user drags across the "view" area of scatter plot viewer
    Then the "x axis min" reading of scatter plot viewer should differ from before
    And the legend of scatter plot viewer should be placed as before
    And the legend of scatter plot viewer should list the same items as before
    When user picks "Reset View" from the context menu of scatter plot viewer
    And user sets "Legend Position" property of scatter plot viewer to "Right"
    Then the legend of scatter plot viewer should be docked
    And no errors should have been logged

  Scenario: A corner legend folded by its chevron stays folded through Position Auto and resizes
    When user sets "Legend Position" property of scatter plot viewer to "RightBottom"
    Then the legend of scatter plot viewer should be in a corner
    When user moves the pointer away from scatter plot viewer
    And user hovers over legend of scatter plot viewer
    And user clicks on legend close chevron of scatter plot viewer
    Then the legend of scatter plot viewer should be collapsed to the mini icon
    When user sets "Legend Position" property of scatter plot viewer to "Auto"
    And user resizes scatter plot viewer to 380 by 300
    Then the legend of scatter plot viewer should be collapsed to the mini icon
    When user resizes scatter plot viewer to 950 by 700
    Then the legend of scatter plot viewer should be collapsed to the mini icon
    When user resizes scatter plot viewer to 700 by 620
    Then the legend of scatter plot viewer should be collapsed to the mini icon
    And mini legend icon of scatter plot viewer should be visible
    And legend of scatter plot viewer should be hidden
    When user clicks on mini legend icon of scatter plot viewer
    Then legend of scatter plot viewer should be visible
    And mini legend icon of scatter plot viewer should be hidden
    When user restores the size of scatter plot viewer
    And user sets "Legend Position" property of scatter plot viewer to "Right"
    Then the legend of scatter plot viewer should be docked
    And no errors should have been logged

  Scenario: A legend column of a thousand categories docks on the right and lists the rows drawn
    When user sets "Color" property of scatter plot viewer to "USUBJID"
    Then the legend of scatter plot viewer should be docked
    And the legend of scatter plot viewer should be in the "right" slot
    And the legend of scatter plot viewer should list 872 items
    And scatter plot viewer should show 872 rows
    When user sets "Legend Position" property of scatter plot viewer to "Auto"
    Then the legend of scatter plot viewer should be docked
    And the legend of scatter plot viewer should list 872 items
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should be in a corner
    When user sets "Color" property of scatter plot viewer to "USUBJID"
    And user sets "Legend Position" property of scatter plot viewer to "Right"
    And user filters rows where "SEX" is "F"
    Then the legend of scatter plot viewer should list 499 items
    And scatter plot viewer should show 499 rows
    When user resets the filter
    And user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    And no errors should have been logged
