@journey @viewers @realizes:viewers.legend
Feature: One legend column across seven viewers
  Seven viewers on one table view share a legend column, and the column's category colors are the
  one source every legend and every canvas paints from: a coloring switched on in the grid reaches
  all seven legends, a color picked in one viewer's legend reaches the grid and the other six, the
  picker's Cancel leaves everything as it was, an empty category gets an item of its own, and the
  palette, the legend column and the legend's visibility come back from a saved layout and a
  saved project.
  One journey on demog-1000 with RACE as the legend column (Caucasian 896, Other 62, Black 27,
  Asian 15 rows; the default palette gives Caucasian #2CA02C and Other #D62728). Each viewer takes
  RACE through the property that feeds its legend: the scatter plot's Color, the histogram's and
  the line chart's Split, the pie chart's Category, the box plot's Category with its Marker Color
  (the box plot's legend is its marker color column, and left alone that column stays on the one
  the plot picked before the category was set), the bar chart's Stack (its legend lists the
  stack, so the chart is split by SEX and stacked by RACE) and the trellis plot's X with an inner
  scatter plot colored by RACE (the trellis legend is its inner viewer's color column). Every
  legend is docked on the right with Visibility Always: on a pane this small the Auto policy would
  fold some of them into the mini icon, which is the placement feature's subject. The grid shares
  the view with seven viewers and shows three columns, so it is scrolled to RACE before its header
  menu is opened.
  The scenarios build on each other the way the manual case does — the palette the grid sets is the
  one the picker changes and the one the layout and the project have to bring back; the last
  scenario puts the column back as it was found.
  The colors set are away from the hues the four RACE categories are painted in by default
  (#FF00FF, #FFFF00, #9467BD — the palette's fifth color, which a four-category column does not
  use), so a canvas that still paints the old palette cannot pass a claim about the new one; the
  "(no value)" item, on a five-category column, gets #E377C2, which none of its categories has.
  The picker is driven on the bar chart, as the color-consistency case says (the
  visibility-and-positioning case uses the scatter plot, whose picker the "(no value)" scenario
  drives). The legend-equals-canvas check of the first scenario reads the two default colors that
  matter here, Caucasian's green and Other's red, on every canvas (the trellis plot's "view" is its
  own small canvas, its cells are inner viewers, so only Caucasian's green is read there); the grid
  scenario then moves them away from the defaults and reads them on every legend and canvas again.
  Not translated: the SPGI columns of the manual case are replaced by demog-1000's (RACE for Stereo
  Category); demog-1000 has no categorical column with empty values, so the "(no value)" item is
  checked on a calculated column that empties RACE for the 153 rows with AGE over 60. The two
  category colors of the grid step are written through the column's categorical map, not through
  the color-coding editor's swatches (the coding itself is switched on from the grid's header menu).

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
    Then "RACE" column should have no color coding
    And no errors should have been logged

  Scenario Outline: The <viewer> legend lists the four races in the colors its canvas paints them
    Then the legend of <viewer> viewer should list 4 items
    And the legend of <viewer> viewer should be docked
    And the "Caucasian" item in the legend of <viewer> viewer should be colored "#2CA02C"
    And the "Other" item in the legend of <viewer> viewer should be colored "#D62728"
    And the "Caucasian" and "Asian" items in the legend of <viewer> viewer should be colored differently
    And the "view" area of <viewer> viewer should contain the color "#2CA02C"
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

  Scenario: Other's red is on every canvas that draws the rows itself
    Then the "view" area of scatter plot viewer should contain the color "#D62728"
    And the "view" area of histogram viewer should contain the color "#D62728"
    And the "view" area of line chart viewer should contain the color "#D62728"
    And the "view" area of bar chart viewer should contain the color "#D62728"
    And the "view" area of pie chart viewer should contain the color "#D62728"
    And the "view" area of box plot viewer should contain the color "#D62728"
    And no errors should have been logged

  Scenario: A categorical coloring switched on in the grid reaches every legend and every canvas
    When user drags the "x scroll handle" area of grid by 60 pixels to the right
    Then grid should have a "header RACE" area
    When user picks "Color Coding > Categorical" from the context menu of the "header RACE" area of grid
    Then "RACE" column should be color-coded categorically
    When user colors "RACE" column categorically:
      | Caucasian | #FF00FF |
      | Other     | #FFFF00 |
    Then the "Caucasian" item in the legend of scatter plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of histogram viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of line chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of bar chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of pie chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of trellis plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of box plot viewer should be colored "#FF00FF"
    And the "Other" item in the legend of scatter plot viewer should be colored "#FFFF00"
    And the "Other" item in the legend of histogram viewer should be colored "#FFFF00"
    And the "Other" item in the legend of line chart viewer should be colored "#FFFF00"
    And the "Other" item in the legend of bar chart viewer should be colored "#FFFF00"
    And the "Other" item in the legend of pie chart viewer should be colored "#FFFF00"
    And the "Other" item in the legend of trellis plot viewer should be colored "#FFFF00"
    And the "Other" item in the legend of box plot viewer should be colored "#FFFF00"
    And the "view" area of scatter plot viewer should contain the color "#FF00FF"
    And the "view" area of histogram viewer should contain the color "#FF00FF"
    And the "view" area of line chart viewer should contain the color "#FF00FF"
    And the "view" area of bar chart viewer should contain the color "#FF00FF"
    And the "view" area of pie chart viewer should contain the color "#FF00FF"
    And the "view" area of trellis plot viewer should contain the color "#FF00FF"
    And the "view" area of box plot viewer should contain the color "#FF00FF"
    And the "view" area of pie chart viewer should not contain the color "#2CA02C"
    And no errors should have been logged

  Scenario: The picker's Cancel puts the color back everywhere
    When user hovers over "Asian" legend item in legend of bar chart viewer
    And user clicks on color picker icon
    Then "Asian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    Then the "Asian" item in the legend of pie chart viewer should be colored "#9467BD"
    When user clicks on CANCEL button in "Asian" dialog
    Then "Asian" dialog should be absent
    And the categorical color of "Asian" in "RACE" column should be "#1F77B4"
    And the "Asian" item in the legend of bar chart viewer should be colored "#1F77B4"
    And the "Asian" item in the legend of pie chart viewer should be colored "#1F77B4"
    And the "Asian" item in the legend of scatter plot viewer should be colored "#1F77B4"
    And no errors should have been logged

  Scenario: A color picked in the bar chart legend reaches the grid and the other six legends
    When user hovers over "Asian" legend item in legend of bar chart viewer
    And user clicks on color picker icon
    Then "Asian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Asian" dialog
    Then "Asian" dialog should be absent
    And the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And "RACE" column should be color-coded categorically
    And the "color of cell 10 of RACE" reading of grid should be "#9467bd"
    And the "Asian" item in the legend of bar chart viewer should be colored "#9467BD"
    And the "Asian" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "Asian" item in the legend of histogram viewer should be colored "#9467BD"
    And the "Asian" item in the legend of line chart viewer should be colored "#9467BD"
    And the "Asian" item in the legend of pie chart viewer should be colored "#9467BD"
    And the "Asian" item in the legend of trellis plot viewer should be colored "#9467BD"
    And the "Asian" item in the legend of box plot viewer should be colored "#9467BD"
    And the "Caucasian" item in the legend of pie chart viewer should be colored "#FF00FF"
    And the "pie" area of pie chart viewer should contain the color "#9467BD"
    And no errors should have been logged

  Scenario: An empty category gets a "(no value)" item whose color can be picked
    When user adds a calculated column "RACE under 61" with formula "if(${AGE} > 60, null, ${RACE})"
    And user sets "colorColumnName" property of scatter plot viewer to "RACE under 61"
    Then the legend of scatter plot viewer should list 5 items
    And "(no value)" legend item in legend of scatter plot viewer should be visible
    When user hovers over "(no value)" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "(no value)" dialog should be visible
    And the "view" area of scatter plot viewer should not contain the color "#E377C2"
    When user picks the color "#E377C2" in the color picker dialog
    And user clicks on OK button in "(no value)" dialog
    Then the "(no value)" item in the legend of scatter plot viewer should be colored "#E377C2"
    And the "view" area of scatter plot viewer should contain the color "#E377C2"
    When user sets "colorColumnName" property of scatter plot viewer to "RACE"
    And user removes "RACE under 61" column
    Then the legend of scatter plot viewer should list 4 items
    And no errors should have been logged

  Scenario: The palette, the legend column and the visibility come back from a saved layout
    When user saves the layout of the current table view to the server
    And user colors "RACE" column categorically:
      | Caucasian | #2CA02C |
      | Asian     | #1F77B4 |
    And user sets "colorColumnName" property of scatter plot viewer to "SEX"
    And user sets "Legend Visibility" property of pie chart viewer to "Never"
    Then the legend of scatter plot viewer should list 2 items
    And legend of pie chart viewer should be hidden
    And the "Caucasian" item in the legend of box plot viewer should be colored "#2CA02C"
    When user loads the saved layout
    Then "colorColumnName" property of scatter plot viewer should be "RACE"
    And the legend of scatter plot viewer should list 4 items
    And legend of pie chart viewer should be visible
    And "Legend Visibility" property of pie chart viewer should be "Always"
    And the categorical color of "Caucasian" in "RACE" column should be "#FF00FF"
    And the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And the "Caucasian" item in the legend of scatter plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of histogram viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of line chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of bar chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of pie chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of trellis plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of box plot viewer should be colored "#FF00FF"
    And the "Asian" item in the legend of box plot viewer should be colored "#9467BD"
    And the "view" area of box plot viewer should contain the color "#FF00FF"
    And no errors should have been logged

  Scenario: The palette comes back on every viewer from a saved project
    When user saves the current view as project "bdd-legend-across-viewers"
    And user closes all views
    And user opens the "bdd-legend-across-viewers" project
    Then the open tableview should have 1 scatter plot viewer
    And the open tableview should have 1 box plot viewer
    And the categorical color of "Caucasian" in "RACE" column should be "#FF00FF"
    And the legend of scatter plot viewer should list 4 items
    And the legend of trellis plot viewer should list 4 items
    And the "Caucasian" item in the legend of scatter plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of histogram viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of line chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of bar chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of pie chart viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of trellis plot viewer should be colored "#FF00FF"
    And the "Caucasian" item in the legend of box plot viewer should be colored "#FF00FF"
    And the "Asian" item in the legend of line chart viewer should be colored "#9467BD"
    And the "view" area of pie chart viewer should contain the color "#FF00FF"
    And no errors should have been logged

  Scenario: The coloring goes back where it was found
    When user removes the coloring of "RACE" column
    Then "RACE" column should have no color coding
    And no errors should have been logged
