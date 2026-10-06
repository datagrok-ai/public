@journey @viewers @realizes:viewers.filters @serial
Feature: The panel's criterion composes with the viewers
  A card's criterion and a viewer's own filtering intersect and neither loses the other: a scatter
  plot zoom narrows the rows the card left and Reset View gives them back, a histogram range and a
  histogram handle narrow further and let go with the viewer still open, closing a viewer releases
  only its share, and the grid shows exactly the rows the card keeps; a zoom box keeps biting when
  the card's criterion changes under it (github-2642); a layout saved with a zoom, a bar click and
  the card comes back onto a closed panel without a balloon (GROK-18281); Select > Invert and
  Select > Selection to Filter replace the filter while the card keeps its criterion and redraws
  (GROK-16713); and a click-filter of a bar chart, a pie chart, a trellis plot and a PC plot slider
  narrows the rows two cards leave, the cards keep their categories, the viewer's own release
  gesture gives the rows back with the viewer open and the panel's Reset filters clears the rest.
  One journey on demog-1000 with a RACE card keeping Caucasian — 896 rows, 416 of them M (93 of
  those RA, 60 of these with SEVERITY None), 633 of them aged 30 to 60; Other holds 62 rows;
  Caucasian and Asian together are 911, the other 89 rows Black or Other. The trellis row splits by
  DIS_POP and SEVERITY, columns no card filters, so its cell cannot hide a card's share.
  Not translated: the zoom coming back with the layout — a layout does not store the scatter plot
  zoom (operator ruling 2026-08-18), so the claim is the card and a narrowed table; the selection
  Invert works on is made through the table's API, the two menu commands are the subject; the
  layout is saved through the API rather than View > Layout > Save to Gallery; the viewers' On Click
  and Show Filters are set as properties and read back, not picked from their context menus; the
  histogram lets go by a double-click on its slider (dev restores the full range that way, which the
  md says it does not) and the PC plot by Reset View; the GROK-16713 redraw is the panel's canvases
  changing across the command after they held still across the moment before it.

  Tagged serial: the "should not have repainted" claims of this journey read the panel while another
  feature's page is working the same stand; it passes alone and failed twice in a full viewers run.

  Background:
    Given user is logged in
    And user opens demog-1000 dataset
    And user opens an empty filter panel
    And user adds a categorical filter on "RACE" keeping "Caucasian"
    Then 896 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"

  Scenario: A scatter plot zoom narrows the rows the card left and Reset View gives them back
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    Then "Zoom and Filter" property of scatter plot viewer should be "filter by zoom"
    When user scrolls the mouse wheel up over the "view" area of scatter plot viewer
    Then fewer than 896 rows should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user picks "Reset View" from the context menu of scatter plot viewer
    Then 896 rows should pass the filter
    When user clicks on close icon of scatter plot viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: A histogram range narrows the card's rows further
    When user adds a histogram viewer with:
      | Value             | AGE  |
      | Show Range Inputs | true |
      | Filtering Enabled | true |
    And user resizes histogram viewer to 500 by 400
    And user enters "30" into the "range min input" area of histogram viewer
    Then 771 rows should pass the filter
    When user enters "60" into the "range max input" area of histogram viewer
    Then 633 rows should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user clicks on close icon of histogram viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: The grid shows exactly the rows the card keeps
    Then the "rows shown" reading of grid should be 896
    When user clicks on the "category Asian of RACE" area of filter panel
    Then 15 rows should pass the filter
    And the "rows shown" reading of grid should be 15
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And the "rows shown" reading of grid should be 896
    And no errors should have been logged

  Scenario: A zoom box narrows the card's rows and a new card criterion keeps the zoom
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    And user drags a zoom box over the "view" area of scatter plot viewer
    Then fewer than 896 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And no rows where "RACE" is "Black" should pass the filter
    And counter of filter panel should have text "1"
    When user clicks on the "category Other of RACE" area of filter panel
    Then fewer than 62 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And no rows where "RACE" is "Caucasian" should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Other"
    And "Zoom and Filter" property of scatter plot viewer should be "filter by zoom"
    And counter of filter panel should have text "1"
    When user clicks on the "category Caucasian of RACE" area of filter panel
    Then no rows where "RACE" is "Other" should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And fewer than 896 rows should pass the filter
    When user clicks on close icon of scatter plot viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: A layout saved with three filtering sources comes back without a balloon
    When user adds a scatter plot viewer with:
      | X | AGE    |
      | Y | HEIGHT |
    And user drags a zoom box over the "view" area of scatter plot viewer
    And user adds a bar chart viewer with:
      | Split | SEX |
    And user sets "On Click" property of bar chart viewer to "Filter"
    And user clicks on the "bar M" area of bar chart viewer
    Then fewer than 416 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    When user saves the layout of the current table view to the server
    And user clicks on close icon of filters viewer
    Then filter panel should be hidden
    When user loads the saved layout
    Then filter panel should be visible
    And "RACE" filter card should be visible
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And fewer than 1000 rows should pass the filter
    And no rows where "RACE" is "Black" should pass the filter
    When user clicks on close icon of scatter plot viewer
    And user clicks on close icon of bar chart viewer
    Then 896 rows should pass the filter
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Invert and Selection to Filter replace the filter and the card keeps its criterion
    When user selects rows where "RACE" is one of "Caucasian, Asian"
    Then 911 rows should be selected
    When user picks "Select > Invert" from the top menu
    Then 89 rows should be selected
    And 896 rows should pass the filter
    When user takes a snapshot of filter panel
    Then filter panel should not have repainted
    When user picks "Select > Selection to Filter" from the top menu
    Then 89 rows should pass the filter
    And no rows where "RACE" is "Caucasian" should pass the filter
    And no rows where "RACE" is "Asian" should pass the filter
    And all rows where "RACE" is "Black" should pass the filter
    And all rows where "RACE" is "Other" should pass the filter
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And filter panel should have repainted
    When user clears the row selection
    And user clicks on the "category Caucasian of RACE" area of filter panel
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario: A histogram handle narrows the card's rows and lets go with the viewer open
    When user adds a histogram viewer with:
      | Value             | AGE  |
      | Filtering Enabled | true |
    And user resizes histogram viewer to 500 by 400
    And user hovers over the "view" area of histogram viewer
    And user drags the "range max handle" area of histogram viewer to the "bin 12" area
    Then fewer than 896 rows should pass the filter
    And the "rows shown" reading of filter panel should be at least 1
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "1"
    When user double-clicks on the "range slider" area of histogram viewer
    Then 896 rows should pass the filter
    When user clicks on close icon of histogram viewer
    Then 896 rows should pass the filter
    And no errors should have been logged

  Scenario Outline: A <viewer> click-filter narrows the rows two cards leave and lets go with the viewer open
    When user clicks on the "category Caucasian of RACE" area of filter panel
    And user adds a card for "SEX" to the filter panel
    And user clicks on the "category M of SEX" area of filter panel
    Then 416 rows should pass the filter
    And counter of filter panel should have text "2"
    When user adds a <viewer> viewer with:
      | <property>  | <value>  |
      | <property2> | <value2> |
      | <property3> | <value3> |
    And user sets "<switch>" property of <viewer> viewer to "<armed>"
    Then "<switch>" property of <viewer> viewer should be "<armed>"
    And 416 rows should pass the filter
    When user <gesture>
    Then <narrowed>
    And the "rows shown" reading of filter panel should be at least 1
    And the "selected categories of SEX" reading of filter panel should be "M"
    And the "selected categories of RACE" reading of filter panel should be "Caucasian"
    And counter of filter panel should have text "2"
    When user <release>
    Then 416 rows should pass the filter
    And the "selected categories of SEX" reading of filter panel should be "M"
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And counter of filter panel should be hidden
    When user clicks on close icon of <viewer> viewer
    And user hovers over "SEX" filter card
    And user clicks on close of "SEX" filter card
    Then all rows should pass the filter
    And no errors should have been logged

    Examples:
      | viewer       | property       | value               | property2      | value2       | property3    | value3       | switch       | armed  | gesture                                                                         | narrowed                                   | release                                                    |
      | bar chart    | Split          | DIS_POP             | Split          | DIS_POP      | Split        | DIS_POP      | On Click     | Filter | clicks on the "bar RA" area of bar chart viewer                                 | 93 rows should pass the filter             | clicks on empty plot space of bar chart viewer             |
      | pie chart    | Category       | DIS_POP             | Category       | DIS_POP      | Category     | DIS_POP      | On Click     | Filter | clicks on the "slice RA" area of pie chart viewer                               | 93 rows should pass the filter             | clicks on the "empty space" area of pie chart viewer       |
      | trellis plot | X Column Names | DIS_POP             | Y Column Names | SEVERITY     | Viewer Type  | Scatter plot | On Click     | Filter | clicks on the "cell RA \| None" area of trellis plot viewer                     | 60 rows should pass the filter             | presses Escape in trellis plot viewer                      |
      | pc plot      | Column Names   | AGE, HEIGHT, WEIGHT | Show Filters   | true         | Show Filters | true         | Show Filters | true   | drags the max handle of the "AGE" range slider of pc plot viewer by 120 pixels | fewer than 416 rows should pass the filter | picks "Reset View" from the context menu of pc plot viewer |
