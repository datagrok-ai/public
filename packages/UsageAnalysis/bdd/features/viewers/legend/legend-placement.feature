@journey @viewers @realizes:viewers.legend
Feature: Where a legend sits, how big it is and when it folds away
  The placement every viewer's legend shares, walked on seven viewers at once: the splitter between
  a docked legend and the plot resizes it; Visibility Always picked from the legend's own menu keeps
  it, and Position Auto moves it to the side the viewer's shape leaves room on; Visibility Auto folds
  it into the mini icon on a small viewer and puts it away on a tiny one; the four corner positions
  lay it over the plot, and the corner legend's chevron folds it into the mini icon, whose click
  opens it again; and the settings come back from a saved layout and a saved project.
  One journey on demog-1000 with RACE as the legend column of all seven viewers, set up as in the
  legend-across-viewers feature, each legend docked on the right with Visibility Always to start
  with. The legend publishes its mode (docked, corner, mini icon, hidden) and its slot, which is
  what every claim reads; a viewer is made small or large by the resize step, which holds the size
  until it is restored. The slots Position Auto picks are the ones the platform's placement policy
  gives each viewer: the five viewers walked there dock the legend on top of a tall narrow box, and
  on a wide flat one the scatter plot, the histogram and the bar chart lay it in a free corner
  while the line chart and the trellis plot dock it on the right.
  The pie chart and the box plot draw their categories themselves, and under Visibility Auto they
  put away a legend that repeats them, so they are left out of the Visibility Auto scenarios (the
  legend's menu cannot be opened on a legend that is not shown); they keep Always from the
  Background and have Position Auto walked in a scenario of their own, where the box plot takes a
  corner on the tall box. Visibility Always is also claimed on a 220-pixel viewer, where Auto
  would fold the legend away. The mini-icon scenario docks the legend on the right first, as the
  manual case turns auto-positioning off. A corner is claimed by the mode and slot the legend
  publishes and by where its box is drawn. The layout scenario changes every viewer's Visibility
  and Position after the save, so every one of them has to come back.
  Not translated: that the mini-legend mode survives a layout and a project — the chevron's fold is
  not a viewer setting, the platform keeps it only until the legend is reopened or its position or
  visibility changes, so the layout and project claims are about the corner positions and the
  Visibility and Position settings.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | xColumnName       | WEIGHT |
      | yColumnName       | HEIGHT |
      | colorColumnName   | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a histogram viewer with:
      | valueColumnName   | AGE    |
      | splitColumnName   | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a line chart viewer with:
      | xColumnName       | AGE    |
      | yColumnNames      | WEIGHT |
      | splitColumnNames  | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a bar chart viewer with:
      | splitColumnName   | SEX    |
      | stackColumnName   | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a pie chart viewer with:
      | categoryColumnName | RACE   |
      | Legend Visibility  | Always |
      | Legend Position    | Right  |
    And user adds a trellis plot viewer with:
      | xColumnNames      | RACE         |
      | Viewer Type       | Scatter plot |
      | Legend Visibility | Always       |
      | Legend Position   | Right        |
    And user sets "colorColumnName" inner property of trellis plot viewer to "RACE"
    And user adds a box plot viewer with:
      | categoryColumnNames   | RACE   |
      | valueColumnName       | AGE    |
      | markerColorColumnName | RACE   |
      | Legend Visibility     | Always |
      | Legend Position       | Right  |
    Then the legend of scatter plot viewer should list 4 items
    And the legend of box plot viewer should list 4 items

  Scenario Outline: The <viewer> legend grows and shrinks with its splitter
    Then the legend of <viewer> viewer should be docked
    And the legend of <viewer> viewer should be in the "right" slot
    When user drags the legend splitter of <viewer> viewer by 40 pixels to the left
    Then the legend of <viewer> viewer should be wider than before
    And the legend of <viewer> viewer should be in the "right" slot
    And the legend of <viewer> viewer should list 4 items
    When user drags the legend splitter of <viewer> viewer by 30 pixels to the right
    Then the legend of <viewer> viewer should be narrower than before
    And no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | histogram    |
      | line chart   |
      | bar chart    |
      | pie chart    |
      | trellis plot |
      | box plot     |

  Scenario Outline: Visibility Always from the <viewer> legend's menu holds, and Position Auto follows the viewer's shape
    When user sets "Legend Visibility" property of <viewer> viewer to "Auto"
    And user right-clicks on legend of <viewer> viewer
    And user picks "Always" from the open menu
    And user closes the context menu
    Then "Legend Visibility" property of <viewer> viewer should be "Always"
    When user sets "Legend Position" property of <viewer> viewer to "Auto"
    And user resizes <viewer> viewer to 320 by 640
    Then the legend of <viewer> viewer should be docked
    And the legend of <viewer> viewer should be in the "top" slot
    When user resizes <viewer> viewer to 800 by 320
    Then the legend of <viewer> viewer should be in the "<wide slot>" slot
    And legend of <viewer> viewer should be visible
    When user resizes <viewer> viewer to 220 by 220
    Then legend of <viewer> viewer should be visible
    And mini legend icon of <viewer> viewer should be hidden
    When user restores the size of <viewer> viewer
    Then no errors should have been logged

    Examples:
      | viewer       | wide slot   |
      | scatter plot | rightBottom |
      | histogram    | rightBottom |
      | line chart   | right       |
      | bar chart    | rightBottom |
      | trellis plot | right       |

  Scenario Outline: Position Auto follows the <viewer>'s shape under Visibility Always
    When user sets "Legend Position" property of <viewer> viewer to "Auto"
    And user resizes <viewer> viewer to 320 by 640
    Then the legend of <viewer> viewer should be in the "<tall slot>" slot
    When user resizes <viewer> viewer to 800 by 320
    Then the legend of <viewer> viewer should be in the "<wide slot>" slot
    And the legend of <viewer> viewer should be docked
    When user restores the size of <viewer> viewer
    Then no errors should have been logged

    Examples:
      | viewer    | tall slot | wide slot |
      | pie chart | top       | right     |
      | box plot  | leftTop   | right     |

  Scenario: Visibility Always and Position Auto come back from a saved layout
    When user saves the layout of the current table view to the server
    And user sets properties of scatter plot viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of histogram viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of line chart viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of bar chart viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of pie chart viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of trellis plot viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    And user sets properties of box plot viewer:
      | Legend Visibility | Never |
      | Legend Position   | Left  |
    Then legend of scatter plot viewer should be hidden
    And "Legend Position" property of box plot viewer should be "Left"
    When user loads the saved layout
    Then properties of scatter plot viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of histogram viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of line chart viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of bar chart viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of pie chart viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of trellis plot viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And properties of box plot viewer should be:
      | Legend Visibility | Always |
      | Legend Position   | Auto   |
    And legend of scatter plot viewer should be visible
    And legend of line chart viewer should be visible
    And the legend of scatter plot viewer should list 4 items
    And no errors should have been logged

  Scenario Outline: Visibility Auto folds the <viewer> legend into the mini icon on a small viewer and puts it away on a tiny one
    When user sets properties of <viewer> viewer:
      | Legend Position   | Right |
      | Legend Visibility | Auto  |
    And user resizes <viewer> viewer to 220 by 220
    Then the legend of <viewer> viewer should be collapsed to the mini icon
    And mini legend icon of <viewer> viewer should be visible
    When user resizes <viewer> viewer to 60 by 60
    Then the legend of <viewer> viewer should be placed nowhere
    And mini legend icon of <viewer> viewer should be hidden
    When user resizes <viewer> viewer to 480 by 480
    Then legend of <viewer> viewer should be visible
    And the legend of <viewer> viewer should be docked
    And the legend of <viewer> viewer should be in the "right" slot
    When user restores the size of <viewer> viewer
    And user sets "Legend Visibility" property of <viewer> viewer to "Always"
    Then no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | histogram    |
      | line chart   |
      | bar chart    |
      | trellis plot |

  Scenario Outline: The <viewer> legend takes each corner, and its chevron folds it into the mini icon
    When user sets "Legend Position" property of <viewer> viewer to "LeftTop"
    Then the legend of <viewer> viewer should be in a corner
    And the legend of <viewer> viewer should be in the "leftTop" slot
    When user sets "Legend Position" property of <viewer> viewer to "LeftBottom"
    Then the legend of <viewer> viewer should be in the "leftBottom" slot
    When user sets "Legend Position" property of <viewer> viewer to "RightTop"
    Then the legend of <viewer> viewer should be in the "rightTop" slot
    When user sets "Legend Position" property of <viewer> viewer to "RightBottom"
    Then the legend of <viewer> viewer should be in a corner
    And the legend of <viewer> viewer should be in the "rightBottom" slot
    And the legend of <viewer> viewer should list 4 items
    And mini legend icon of <viewer> viewer should be hidden
    When user hovers over legend of <viewer> viewer
    And user clicks on legend close chevron of <viewer> viewer
    Then the legend of <viewer> viewer should be collapsed to the mini icon
    And mini legend icon of <viewer> viewer should be visible
    And legend of <viewer> viewer should be hidden
    When user clicks on mini legend icon of <viewer> viewer
    Then the legend of <viewer> viewer should be in a corner
    And the legend of <viewer> viewer should be in the "rightBottom" slot
    And no errors should have been logged

    Examples:
      | viewer       |
      | scatter plot |
      | histogram    |
      | line chart   |
      | bar chart    |
      | pie chart    |
      | trellis plot |
      | box plot     |

  Scenario: Corner positions come back from a saved layout
    When user sets "Legend Position" property of scatter plot viewer to "LeftTop"
    And user sets "Legend Position" property of pie chart viewer to "LeftBottom"
    And user sets "Legend Position" property of box plot viewer to "RightTop"
    And user saves the layout of the current table view to the server
    And user sets "Legend Position" property of scatter plot viewer to "Right"
    And user sets "Legend Position" property of pie chart viewer to "Right"
    And user sets "Legend Position" property of box plot viewer to "Right"
    Then the legend of scatter plot viewer should be docked
    When user loads the saved layout
    Then the legend of scatter plot viewer should be in the "leftTop" slot
    And the legend of scatter plot viewer should be in a corner
    And the legend of pie chart viewer should be in the "leftBottom" slot
    And the legend of box plot viewer should be in the "rightTop" slot
    And the legend of bar chart viewer should be in the "rightBottom" slot
    And the legend of trellis plot viewer should be in a corner
    And no errors should have been logged

  Scenario: Corner positions and the visibility come back from a saved project
    When user sets "Legend Visibility" property of histogram viewer to "Never"
    And user sets "Legend Position" property of bar chart viewer to "Auto"
    And user saves the current view as project "bdd-legend-placement"
    And user closes all views
    And user opens the "bdd-legend-placement" project
    Then the open tableview should have 1 scatter plot viewer
    And the legend of scatter plot viewer should be in the "leftTop" slot
    And the legend of pie chart viewer should be in the "leftBottom" slot
    And the legend of box plot viewer should be in the "rightTop" slot
    And the legend of line chart viewer should be in the "rightBottom" slot
    And the legend of line chart viewer should be in a corner
    And "Legend Visibility" property of histogram viewer should be "Never"
    And legend of histogram viewer should be hidden
    And "Legend Visibility" property of scatter plot viewer should be "Always"
    And "Legend Position" property of bar chart viewer should be "Auto"
    And no errors should have been logged
