@journey @viewers @realizes:viewers.scatter-plot
Feature: Scatter plot legend
  What the scatter plot puts in its legend and what the legend does back to the plot: a categorical
  Color builds it, Markers adds a second section to the same legend, a numeric Color leaves only
  the marker section, clearing Markers takes it away again, a table filter shrinks the marker
  section and giving the filter back restores it, a click on an entry hides that category on the
  canvas without touching the table filter, and Legend Visibility and Legend Position decide
  whether and where it shows. The item counts are what the legend lists (every section, rendered
  or not — a corner legend renders only the rows that fit). The item-level rendering — glyphs,
  colors, the cross — belongs to the legend's own feature. One journey on demog-1000, X = WEIGHT,
  Y = HEIGHT (872 of the 1000 rows have a HEIGHT); every scenario puts back what it changed.
  The later scenarios translate the Legend section's scatter plot cases on the same data: colors
  picked in the legend (the palette icon of a hovered item) reach the column and come back from a
  layout and a project; a categorical Color formula that empties RACE Other keeps a "(no value)"
  item and an ID column (USUBJID) lists one item per drawn row; a numerical Color with empty values
  draws the color scale and no item; X columns that empty one race or the others leave only the
  races they draw, and a zoom-filter on them narrows the legend to the box (rows with no X are not
  placed, so they stay in the table's filter); two scatter plots with the same viewer Filter; a
  legend click on top of a filter panel criterion narrows the drawn rows further and brings none
  back; and a coloring from the grid — Linear from the header menu, a new scheme and its inversion,
  Apply to text in the Color-coding dialog, then a Conditional coding — reaches the scatter plot's
  legend or scale and the box plot's markers, through a layout and a project.
  The manual case's SPGI columns are demog-1000's: Series is RACE, Stereo Category is RACE, Primary
  scaffold name is DIS_POP, Chemical Space X is WEIGHT. demog-1000's numeric columns offer no
  Categorical coding (Off, Conditional, Linear, Linked), so the step from linear to categorical
  goes to Conditional, the coding that gives a numeric column items with colors of their own; the
  new scheme's stops and the inversion are written through the column (the dialog's scheme picker
  has no named swatches).
  Known failure, GROK-20896: a layout that holds two scatter plots with the same Filter brings the
  second one back with its Filter set but not applied — it draws all 872 rows and its legend lists
  all four races, while the first one comes back drawing its 38. The operator reproduced it by hand
  and filed the ticket; reproduced in the suite too, and outside it with plain JS: the second
  viewer's rows shown goes from 38 to 872 across saveLayout/loadLayout whether the first viewer's
  Filter, the second's or neither was changed in between. With one scatter plot the layout is
  restored correctly. Take the mark off when the ticket is closed.
  The numerical Color with empty values is legend-ui's gradient case: its scale is drawn in one
  hue, so it is claimed as painted with finite ends rather than as a gradient of several colors.
  The box plot's markers under a numeric Marker Color are already painted in every hue, so a new
  scheme and its inversion are claimed on them as repaints, and on the scatter plot's scale by the
  scheme's own colors. After the layout the linear scheme, its inversion and Apply to text are read
  back from the column and the scatter plot's scale; before the layout the scheme is overwritten,
  so only the layout can bring it back.
  Not translated here: Marker = Core with Color = ID — the molecule markers are the subject of
  legend/legend-structures.feature (spgi-100); the PC plot of the color-coding case (this feature
  claims the scatter plot and the box plot); modifying a few of the categorical colors after the
  switch — the Conditional coding is written with its two bins' colors at once.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user adds a scatter plot viewer with:
      | X | WEIGHT |
      | Y | HEIGHT |
    Then scatter plot viewer should show 872 rows
    And legend of scatter plot viewer should be hidden

  Scenario: Color and Markers build one legend of two sections
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 4 items
    And legend of scatter plot viewer should contain text "Asian"
    And legend of scatter plot viewer should contain text "Caucasian"
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 6 items
    When user sets "Color" property of scatter plot viewer to "AGE"
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Color" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 6 items
    And legend of scatter plot viewer should contain text "Asian"
    And no errors should have been logged

  Scenario: Clearing Markers leaves the color entries alone
    When user sets properties of scatter plot viewer:
      | Color   | SEX |
      | Markers | SEX |
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to ""
    Then "Markers" property of scatter plot viewer should be ""
    And the legend of scatter plot viewer should list 2 items
    And legend of scatter plot viewer should contain text "F"
    When user sets properties of scatter plot viewer:
      | Color   | RACE |
      | Markers | SEX  |
    Then the legend of scatter plot viewer should list 6 items
    And no errors should have been logged

  Scenario: A table filter drops the filtered-out categories and Markers on the color column adds no section
    When user sets "Markers" property of scatter plot viewer to ""
    Then the legend of scatter plot viewer should list 4 items
    When user opens the filter panel
    And user adds a categorical filter on "RACE" keeping "Asian, Caucasian"
    Then 911 rows should pass the filter
    And the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 2 items
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 4 items
    When user resets the filter
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 6 items
    When user sets "Markers" property of scatter plot viewer to "RACE"
    Then the legend of scatter plot viewer should list 4 items
    When user sets "Markers" property of scatter plot viewer to "SEX"
    Then the legend of scatter plot viewer should list 6 items
    And no errors should have been logged

  Scenario: A click on a legend entry hides that category on the canvas, not in the table
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    When user clicks on "Asian" item in the legend of scatter plot viewer
    Then scatter plot viewer should show fewer rows than before
    And scatter plot viewer should have less ink than before
    And all rows should pass the filter
    When user clicks on "Asian" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 872 rows
    And all rows should pass the filter
    And no errors should have been logged

  Scenario: Legend Visibility and Legend Position decide whether and where it shows
    When user sets "Legend Visibility" property of scatter plot viewer to "Never"
    Then legend of scatter plot viewer should be hidden
    When user sets "Legend Visibility" property of scatter plot viewer to "Always"
    Then legend of scatter plot viewer should be visible
    And the legend of scatter plot viewer should list 6 items
    When user sets "Legend Position" property of scatter plot viewer to "Left"
    Then the legend of scatter plot viewer should be on the left
    When user sets "Legend Position" property of scatter plot viewer to "Right"
    Then the legend of scatter plot viewer should be on the right
    When user sets properties of scatter plot viewer:
      | Legend Position   | Auto |
      | Legend Visibility | Auto |
      | Color             |      |
      | Markers           |      |
    Then legend of scatter plot viewer should be hidden
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: Colors picked in the legend reach the column and come back from a layout and a project
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then the "filters" reading of filter panel should be 0
    And all rows should pass the filter
    When user sets properties of scatter plot viewer:
      | Color   | RACE |
      | Markers | RACE |
    Then the legend of scatter plot viewer should list 4 items
    And marker of "Asian" legend item in legend of scatter plot viewer should be visible
    When user hovers over "Asian" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "Asian" dialog should be visible
    When user picks the color "#9467BD" in the color picker dialog
    And user clicks on OK button in "Asian" dialog
    And user hovers over "Black" legend item in legend of scatter plot viewer
    And user clicks on color picker icon
    Then "Black" dialog should be visible
    When user picks the color "#BCBD22" in the color picker dialog
    And user clicks on OK button in "Black" dialog
    Then the categorical color of "Black" in "RACE" column should be "#BCBD22"
    And the categorical color of "Asian" in "RACE" column should be "#9467BD"
    And the "Asian" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "Black" item in the legend of scatter plot viewer should be colored "#BCBD22"
    And the "view" area of scatter plot viewer should contain the color "#BCBD22"
    When user saves the layout of the current table view to the server
    And user colors "RACE" column categorically:
      | Asian | #1F77B4 |
      | Black | #FFBB78 |
    And user sets "Markers" property of scatter plot viewer to ""
    Then the "Black" item in the legend of scatter plot viewer should be colored "#FFBB78"
    When user loads the saved layout
    Then "Markers" property of scatter plot viewer should be "RACE"
    And the legend of scatter plot viewer should list 4 items
    And the "Asian" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "Black" item in the legend of scatter plot viewer should be colored "#BCBD22"
    And the "view" area of scatter plot viewer should contain the color "#BCBD22"
    When user saves the current view as project "bdd-scatter-plot-legend-colors"
    And user closes all views
    And user opens the "bdd-scatter-plot-legend-colors" project
    Then the open tableview should have 1 scatter plot viewer
    And the legend of scatter plot viewer should list 4 items
    And the "Asian" item in the legend of scatter plot viewer should be colored "#9467BD"
    And the "Black" item in the legend of scatter plot viewer should be colored "#BCBD22"
    And the "view" area of scatter plot viewer should contain the color "#9467BD"
    When user colors "RACE" column categorically:
      | Asian | #1F77B4 |
      | Black | #FFBB78 |
    And user removes the coloring of "RACE" column
    And user sets properties of scatter plot viewer:
      | Color   |  |
      | Markers |  |
    Then legend of scatter plot viewer should be hidden
    And no errors should have been logged

  Scenario: A categorical Color formula with empty values lists them as an item of their own
    When user adds a calculated column "DIS_POP but Other" with formula "if(${RACE}=='Other', null, ${DIS_POP})"
    And user sets "Color" property of scatter plot viewer to "DIS_POP but Other"
    Then the legend of scatter plot viewer should list 6 items
    And "(no value)" legend item in legend of scatter plot viewer should be visible
    And "RA" legend item in legend of scatter plot viewer should be visible
    And the "(no value)" and "RA" items in the legend of scatter plot viewer should be colored differently
    When user sets "Color" property of scatter plot viewer to "USUBJID"
    Then the legend of scatter plot viewer should list 872 items
    And the legend of scatter plot viewer should be docked
    When user sets "Color" property of scatter plot viewer to ""
    And user removes "DIS_POP but Other" column
    Then legend of scatter plot viewer should be hidden
    And no errors should have been logged

  Scenario: A numerical Color column with empty values draws a scale, not items
    When user adds a calculated column "WEIGHT but Other" with formula "if(${RACE}=='Other', null, ${WEIGHT})"
    And user sets "Color" property of scatter plot viewer to "WEIGHT but Other"
    Then the legend of scatter plot viewer should list 0 items
    And scatter plot viewer should have a "color scale" area
    And the "color scale" area of scatter plot viewer should be painted
    And the "color scale min" reading of scatter plot viewer should be a finite number
    And the "color scale max" reading of scatter plot viewer should be 165
    When user sets "Color" property of scatter plot viewer to ""
    And user removes "WEIGHT but Other" column
    Then scatter plot viewer should not have a "color scale" area
    And no errors should have been logged

  Scenario: The legend follows an X column that empties other races, and a zoom-filter on it
    When user adds a calculated column "WEIGHT of Caucasians" with formula "if(${RACE}!='Caucasian', null, ${WEIGHT})"
    And user adds a calculated column "WEIGHT of the rest" with formula "if(${RACE}=='Caucasian', null, ${WEIGHT})"
    And user sets properties of scatter plot viewer:
      | X     | WEIGHT of Caucasians |
      | Color | RACE                 |
    Then the legend of scatter plot viewer should list 1 item
    And "Caucasian" legend item in legend of scatter plot viewer should be visible
    When user sets "X" property of scatter plot viewer to "WEIGHT of the rest"
    Then the legend of scatter plot viewer should list 3 items
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    And scatter plot viewer should show 97 rows
    When user sets "Zoom And Filter" property of scatter plot viewer to "filter by zoom"
    And user drags from the "marker of row 133" area to the "marker of row 178" area of scatter plot viewer holding Alt
    Then the legend of scatter plot viewer should list 1 item
    And "Other" legend item in legend of scatter plot viewer should be visible
    And scatter plot viewer should show fewer rows than before
    And fewer than 1000 rows should pass the filter
    And all rows where "RACE" is "Caucasian" should pass the filter
    When user double-clicks on empty plot space of scatter plot viewer
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 3 items
    When user sets properties of scatter plot viewer:
      | Zoom And Filter | no action |
      | X               | WEIGHT    |
      | Color           |           |
    And user removes "WEIGHT of Caucasians" column
    And user removes "WEIGHT of the rest" column
    Then scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: Two scatter plots with the same viewer filter both narrow their legends, and a layout brings the first back
    When user sets properties of scatter plot viewer:
      | Markers           | RACE                          |
      | Filter            | ${RACE} in ["Asian", "Black"] |
      | Legend Visibility | Always                        |
      | Legend Position   | Right                         |
    And user adds a scatter plot viewer
    Then the open tableview should have 2 scatter plot viewers
    When user sets properties of last scatter plot viewer:
      | X                 | WEIGHT                        |
      | Y                 | HEIGHT                        |
      | Markers           | RACE                          |
      | Filter            | ${RACE} in ["Asian", "Black"] |
      | Legend Visibility | Always                        |
      | Legend Position   | Right                         |
    Then the legend of first scatter plot viewer should list 2 items
    And the legend of last scatter plot viewer should list 2 items
    And "Caucasian" legend item in legend of last scatter plot viewer should be absent
    And the "rows shown" reading of last scatter plot viewer should be 38
    When user saves the layout of the current table view to the server
    And user sets "Filter" property of first scatter plot viewer to ""
    Then the legend of first scatter plot viewer should list 4 items
    When user loads the saved layout
    Then the open tableview should have 2 scatter plot viewers
    And "Filter" property of first scatter plot viewer should be '${RACE} in ["Asian", "Black"]'
    And "Filter" property of last scatter plot viewer should be '${RACE} in ["Asian", "Black"]'
    And the legend of first scatter plot viewer should list 2 items
    And the "rows shown" reading of first scatter plot viewer should be 38
    And no errors should have been logged

  Scenario: The second scatter plot draws its Filter again after the layout (GROK-20896)
    Then the "rows shown" reading of last scatter plot viewer should be 38
    And the legend of last scatter plot viewer should list 2 items

  Scenario: The second scatter plot goes and the first gives its settings back
    When user clicks on close icon of last scatter plot viewer
    Then the open tableview should have 1 scatter plot viewer
    When user sets properties of scatter plot viewer:
      | Markers           |      |
      | Filter            |      |
      | Legend Visibility | Auto |
      | Legend Position   | Auto |
    Then legend of scatter plot viewer should be hidden
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A legend click on top of a panel filter narrows further and brings nothing back
    When user sets properties of scatter plot viewer:
      | Color   | DIS_POP |
      | Markers | RACE    |
    Then the legend of scatter plot viewer should list 9 items
    When user opens an empty filter panel
    And user adds a card for "DIS_POP" to the filter panel
    And user clicks on the "category RA of DIS_POP" area of filter panel
    Then 434 rows should pass the filter
    When user clicks on the "checkbox Psoriasis of DIS_POP" area of filter panel
    Then 638 rows should pass the filter
    And scatter plot viewer should show 633 rows
    And the legend of scatter plot viewer should list 6 items
    And "Indigestion" legend item in legend of scatter plot viewer should be absent
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 19 rows
    And 638 rows should pass the filter
    And "Black" legend item in legend of scatter plot viewer should be selected
    When user clicks on "Black" item in the legend of scatter plot viewer
    Then scatter plot viewer should show 633 rows
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    And user sets properties of scatter plot viewer:
      | Color   |  |
      | Markers |  |
    Then all rows should pass the filter
    And scatter plot viewer should show 872 rows
    And no errors should have been logged

  Scenario: A linear coloring from the grid gives the scatter plot a scale and colors the box plot's markers
    When user sets "Color" property of scatter plot viewer to "WEIGHT"
    And user adds a box plot viewer with:
      | categoryColumnNames   | RACE   |
      | valueColumnName       | AGE    |
      | markerColorColumnName | WEIGHT |
    And user drags the "x scroll handle" area of grid by 100 pixels to the right
    Then grid should have a "header WEIGHT" area
    When user picks "Color Coding > Linear" from the context menu of the "header WEIGHT" area of grid
    Then "WEIGHT" column should be color-coded linearly
    And the legend of scatter plot viewer should list 0 items
    And scatter plot viewer should have a "color scale" area
    And the "color scale" area of scatter plot viewer should not contain the color "#FF00FF"
    When user colors "WEIGHT" column linearly from "#00FF00" to "#FF00FF"
    Then the "color scale" area of scatter plot viewer should contain the color "#FF00FF"
    And box plot viewer should have repainted
    When user inverts the color scheme of "WEIGHT" column
    Then the color scheme of "WEIGHT" column should be "#FF00FF, #00FF00"
    And the "color scale" area of scatter plot viewer should have repainted
    And box plot viewer should have repainted
    When user picks "Color Coding > Edit..." from the context menu of the "header WEIGHT" area of grid
    Then "Color-coding: WEIGHT" dialog should be visible
    When user selects "text" in "Apply to" input in "Color-coding: WEIGHT" dialog
    Then "Apply to" input in "Color-coding: WEIGHT" dialog should have value "text"
    When user clicks on CLOSE button in "Color-coding: WEIGHT" dialog
    Then the text of "WEIGHT" column should be color-coded
    And the legend of scatter plot viewer should list 0 items
    And the "color scale" area of scatter plot viewer should contain the color "#00FF00"
    When user saves the layout of the current table view to the server
    And user colors "WEIGHT" column linearly from "#FFFF00" to "#00FFFF"
    And user removes the coloring of "WEIGHT" column
    Then "WEIGHT" column should have no color coding
    And the color scheme of "WEIGHT" column should be "#FFFF00, #00FFFF"
    When user loads the saved layout
    Then "WEIGHT" column should be color-coded linearly
    And the color scheme of "WEIGHT" column should be "#FF00FF, #00FF00"
    And the text of "WEIGHT" column should be color-coded
    And the "color scale" area of scatter plot viewer should contain the color "#00FF00"
    And no errors should have been logged

  Scenario: A conditional coloring turns the scale into items, and a project brings them back
    When user colors "WEIGHT" column conditionally:
      | 40-100  | #00FFFF |
      | 100-170 | #FFA500 |
    Then "WEIGHT" column should be color-coded conditionally
    And the legend of scatter plot viewer should list 2 items
    And the "40-100" item in the legend of scatter plot viewer should be colored "#00FFFF"
    And the "100-170" item in the legend of scatter plot viewer should be colored "#FFA500"
    And scatter plot viewer should not have a "color scale" area
    And the "view" area of box plot viewer should contain the color "#00FFFF"
    When user saves the current view as project "bdd-scatter-plot-legend-coding"
    And user closes all views
    And user opens the "bdd-scatter-plot-legend-coding" project
    Then "WEIGHT" column should be color-coded conditionally
    And the legend of scatter plot viewer should list 2 items
    And the "40-100" item in the legend of scatter plot viewer should be colored "#00FFFF"
    And the "view" area of box plot viewer should contain the color "#FFA500"
    When user removes the coloring of "WEIGHT" column
    And user clicks on close icon of box plot viewer
    Then the open tableview should have 0 box plot viewers
    When user sets "Color" property of scatter plot viewer to ""
    Then legend of scatter plot viewer should be hidden
    And no errors should have been logged
