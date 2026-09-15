@journey @viewers @realizes:viewers.legend
Feature: Legends follow the rows that are left
  Seven legends on one table view list only the categories whose rows are still there, whatever
  took the others away: a range and a category filter on the filter panel, the filter panel's
  reset, a viewer's own Filter formula alone and together with the panel, a scatter plot's
  zoom-filter, a click on a bar, a wedge or a trellis cell with On Click set to Filter, and the
  scatter plot's Row Source. A stacked bar chart with Include Nulls off lists no stack category the
  panel took away, and a saved layout brings the panel's filter back with the legends it left.
  One journey on demog-1000 with RACE as the legend column of all seven viewers, set up as in the
  legend-across-viewers feature (scatter plot Color, histogram and line chart Split, bar chart
  Stack, pie chart Category, trellis plot X with an inner scatter plot colored by RACE, box plot
  Category and Marker Color), every legend docked on the right with Visibility Always. The numbers
  are demog-1000's: WEIGHT 100..170 keeps 143 rows of three races (Caucasian 131, Black 6, Other 6
  — no Asian weighs over 91.8), RACE Black or Other then keeps 12 of them, and every row with AGE
  over 70 is Caucasian. Each scenario from the third on puts back the filters it set; the first two
  are one case in the manual one — filter, then reset. The trellis plot draws scatter plots in its
  cells, so a scenario that reads "scatter plot viewer" after the view was rebuilt checks first that
  the phrase reached the scatter plot (its X column), not a cell of the trellis.
  The categories are clicked on the panel's RACE and DIS_POP cards and the panel's own reset icon
  clears them; the WEIGHT range goes in through the filter group's API — the state a handle drag
  leaves. A range card draws a histogram of its own, earlier in the page than the histogram viewer,
  so the scenarios name that one "last histogram viewer".
  The layout scenarios come last and close the histogram viewer first: with the immediate rendering
  the harness puts every viewer in, a layout that holds a histogram viewer comes back with an empty
  filter panel, so the claims are made on the other six viewers.
  On Click = Filter comes with Row Source = All, and that is the core's own doing, confirmed by the
  operator: `viewer_base.dart`'s On Click menu sets the two together, and the trellis plot couples
  them in `onLookChanged` as well. A viewer that filters on a click therefore keeps drawing every
  row, and its own legend rightly lists all four races while the other six narrow to the clicked
  one. The bar chart's and the pie chart's claims are made through the viewer's context menu for
  that reason — the property set on its own is a different path, and outside the trellis it leaves
  Row Source alone. The bar chart makes both claims in one scenario, the menu's pair and the click
  that filters under it. The pie chart's menu is driven in a scenario of its own, near the end of
  the journey: under Row Source All its wedge click leaves the table's filter untouched, and once
  the pie chart has been through the menu the click stays inert for the rest of the run, so the
  wedge-click scenario comes earlier and sets On Click by itself, keeping Row Source on Filtered.
  A layout does not bring a click-filter back, and that too is the right behaviour by the operator's
  decision: the scenario that saves a layout while a bar click is filtering loads it again and
  claims the panel's own range filter returned while the click-filter did not — every row of the
  range passing and every legend back to the range's three races.
  Not translated: the structure filter on Core (a Chem filter, not the legend's). The manual case's numeric filter (Average Mass >
  400) is WEIGHT 100..170 here, its Stereo Category filter is RACE. Its composition case (the
  viewer's Filter with Average Mass > 300, legend unchanged at two) becomes the viewer's Filter with
  WEIGHT 100..170, whose intersection leaves one race — a composition the legend can only show by
  narrowing further. Its bar chart edge case (count of CAST Idea ID stacked by Primary Scaffold
  Name) is a count of USUBJID stacked by a calculated DIS_POP that is empty for RACE Other —
  demog-1000 has no empty category of its own, and Include Nulls has nothing to take away without
  one.
  Known failure: applying a layout that holds a histogram viewer races the filter panel the layout
  re-creates. The new panel builds its cards from the layout on a 10 ms timer (`filters_core.dart`),
  while the histogram's filter request makes the panel save its state, and the save rewrites the
  panel's filters from the cards it already has, which before the first build is none. The save
  normally waits 50 ms and lands after the build; under `immediateRendering`, which
  `libraries/bdd/src/runtime/viewer-runtime.ts`'s `arm` sets on every viewer, it lands first, and the
  layout comes back with an empty filter panel and every row passing. The saved layout is the same
  with the flag and without it. A layout holding no histogram viewer, or one whose Filtering Enabled
  is off, restores the panel, so the round-trips before the last scenario are made with the
  histogram closed. Take the mark off once the save no longer runs before the panel's first build.

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
      | splitColumnName   | RACE   |
      | stackColumnName   | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    And user adds a pie chart viewer with:
      | categoryColumnName | RACE   |
      | Legend Visibility  | Always |
      | Legend Position    | Right  |
    And user adds a trellis plot viewer with:
      | xColumnNames      | RACE         |
      | yColumnNames      |              |
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
    And user opens an empty filter panel
    And user adds a card for "RACE" to the filter panel
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 4 items
    And the legend of box plot viewer should list 4 items

  Scenario: A range filter and a category filter on the panel narrow every legend
    When user adds a range filter on "WEIGHT" from 100 to 170
    Then 143 rows should pass the filter
    And the legend of scatter plot viewer should list 3 items
    And the legend of last histogram viewer should list 3 items
    And the legend of line chart viewer should list 3 items
    And the legend of bar chart viewer should list 3 items
    And the legend of pie chart viewer should list 3 items
    And the legend of trellis plot viewer should list 3 items
    And the legend of box plot viewer should list 3 items
    And "Asian" legend item in legend of pie chart viewer should be absent
    And "Black" legend item in legend of pie chart viewer should be visible
    When user clicks on the "category Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    Then 12 rows should pass the filter
    And the legend of scatter plot viewer should list 2 items
    And the legend of last histogram viewer should list 2 items
    And the legend of line chart viewer should list 2 items
    And the legend of bar chart viewer should list 2 items
    And the legend of pie chart viewer should list 2 items
    And the legend of trellis plot viewer should list 2 items
    And the legend of box plot viewer should list 2 items
    And "Caucasian" legend item in legend of box plot viewer should be absent
    And no errors should have been logged

  Scenario: The panel's reset gives every legend all four races back
    Then 12 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 4 items
    And the legend of last histogram viewer should list 4 items
    And the legend of line chart viewer should list 4 items
    And the legend of bar chart viewer should list 4 items
    And the legend of pie chart viewer should list 4 items
    And the legend of trellis plot viewer should list 4 items
    And the legend of box plot viewer should list 4 items
    And no errors should have been logged

  Scenario: The scatter plot's own Filter narrows its legend and no other
    When user sets "Filter" property of scatter plot viewer to '${RACE} in ["Asian", "Black"]'
    Then the legend of scatter plot viewer should list 2 items
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    And all rows should pass the filter
    And the legend of last histogram viewer should list 4 items
    And the legend of line chart viewer should list 4 items
    And the legend of bar chart viewer should list 4 items
    And the legend of pie chart viewer should list 4 items
    And the legend of trellis plot viewer should list 4 items
    And the legend of box plot viewer should list 4 items
    And no errors should have been logged

  Scenario: The viewer's Filter and the panel's filter compose
    When user adds a range filter on "WEIGHT" from 100 to 170
    Then 143 rows should pass the filter
    And the legend of scatter plot viewer should list 1 item
    And "Black" legend item in legend of scatter plot viewer should be visible
    And "Asian" legend item in legend of scatter plot viewer should be absent
    And the legend of last histogram viewer should list 3 items
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    And user sets "Filter" property of scatter plot viewer to ""
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 4 items
    And no errors should have been logged

  Scenario: A zoom-filter on the scatter plot drops the races with no point in the box
    When user sets properties of scatter plot viewer:
      | xColumnName   | AGE            |
      | yColumnName   | WEIGHT         |
      | Zoom And Filter | filter by zoom |
    And user drags from the "marker of row 5" area to the "marker of row 20" area of scatter plot viewer holding Alt
    Then the legend of scatter plot viewer should list 2 items
    And "Caucasian" legend item in legend of scatter plot viewer should be visible
    And "Asian" legend item in legend of scatter plot viewer should be absent
    And fewer than 30 rows should pass the filter
    And the legend of pie chart viewer should list 2 items
    When user double-clicks on empty plot space of scatter plot viewer
    Then all rows should pass the filter
    And the legend of scatter plot viewer should list 4 items
    When user sets properties of scatter plot viewer:
      | Zoom And Filter | no action |
      | xColumnName     | WEIGHT    |
      | yColumnName     | HEIGHT    |
    Then no errors should have been logged

  Scenario: A click on a bar with On Click Filter leaves only its race in the other legends
    When user picks "On Click > Filter" from the context menu of bar chart viewer
    Then "On Click" property of bar chart viewer should be "Filter"
    And "Row Source" property of bar chart viewer should be "All"
    And the legend of bar chart viewer should list 4 items
    When user clicks on the "bar Black" area of bar chart viewer
    Then 27 rows should pass the filter
    And the legend of bar chart viewer should list 4 items
    And "Black" legend item in legend of bar chart viewer should be visible
    And the legend of scatter plot viewer should list 1 item
    And the legend of box plot viewer should list 1 item
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of bar chart viewer should list 4 items
    When user sets properties of bar chart viewer:
      | On Click   | None     |
      | Row Source | Filtered |
    Then no errors should have been logged

  Scenario: A click on a wedge with On Click Filter leaves only its race in the other legends
    When user sets "On Click" property of pie chart viewer to "Filter"
    And user clicks on the 'slice "Other"' area of pie chart viewer
    Then 62 rows should pass the filter
    And the legend of pie chart viewer should list 1 item
    And "Other" legend item in legend of pie chart viewer should be visible
    And the legend of last histogram viewer should list 1 item
    And the legend of line chart viewer should list 1 item
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of pie chart viewer should list 4 items
    When user sets properties of pie chart viewer:
      | On Click   | None     |
      | Row Source | Filtered |
    Then no errors should have been logged

  Scenario: A click on a trellis cell with On Click Filter leaves only its race in the other legends
    When user sets "On Click" property of trellis plot viewer to "Filter"
    And user clicks on the "cell body Asian" area of trellis plot viewer
    Then 15 rows should pass the filter
    And "Row Source" property of trellis plot viewer should be "All"
    And the legend of trellis plot viewer should list 4 items
    And the legend of pie chart viewer should list 1 item
    And "Asian" legend item in legend of pie chart viewer should be visible
    And the legend of line chart viewer should list 1 item
    And the legend of box plot viewer should list 1 item
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of pie chart viewer should list 4 items
    When user sets properties of trellis plot viewer:
      | On Click   | None     |
      | Row Source | Filtered |
    Then the legend of trellis plot viewer should list 4 items
    And no errors should have been logged

  Scenario: The scatter plot's Row Source decides which rows its legend lists
    Given "xColumnName" property of scatter plot viewer should be "WEIGHT"
    When user clicks on the "category Asian of RACE" area of filter panel
    And user clicks on the "checkbox Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    And user selects rows where "RACE" is one of "Black, Caucasian"
    Then 104 rows should pass the filter
    And 923 rows should be selected
    When user sets "Row Source" property of scatter plot viewer to "All"
    Then the legend of scatter plot viewer should list 4 items
    When user sets "Row Source" property of scatter plot viewer to "Filtered"
    Then the legend of scatter plot viewer should list 3 items
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    When user sets "Row Source" property of scatter plot viewer to "FilteredSelected"
    Then the legend of scatter plot viewer should list 1 item
    And "Black" legend item in legend of scatter plot viewer should be visible
    When user sets "Row Source" property of scatter plot viewer to "Selected"
    Then the legend of scatter plot viewer should list 2 items
    And "Caucasian" legend item in legend of scatter plot viewer should be visible
    And "Asian" legend item in legend of scatter plot viewer should be absent
    When user sets "Row Source" property of scatter plot viewer to "Filtered"
    And user clears the row selection
    And user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And no errors should have been logged

  Scenario: A stacked bar chart with Include Nulls off lists no category the panel took away
    When user adds a calculated column "DIS_POP but Other" with formula "if(${RACE}=='Other', null, ${DIS_POP})"
    And user sets properties of bar chart viewer:
      | valueColumnName | USUBJID           |
      | valueAggrType   | count             |
      | splitColumnName | RACE              |
      | stackColumnName | DIS_POP but Other |
      | Include Nulls   | true              |
    Then the legend of bar chart viewer should list 7 items
    And "(no value)" legend item in legend of bar chart viewer should be visible
    When user sets "Include Nulls" property of bar chart viewer to "false"
    Then the legend of bar chart viewer should list 6 items
    And "(no value)" legend item in legend of bar chart viewer should be absent
    When user adds a card for "DIS_POP but Other" to the filter panel
    And user clicks on the "category RA of DIS_POP but Other" area of filter panel
    And user clicks on the "checkbox Psoriasis of DIS_POP but Other" area of filter panel
    And user clicks on the "checkbox UC of DIS_POP but Other" area of filter panel
    Then fewer than 761 rows should pass the filter
    And the legend of bar chart viewer should list 3 items
    And "Indigestion" legend item in legend of bar chart viewer should be absent
    And "(no value)" legend item in legend of bar chart viewer should be absent
    And "RA" legend item in legend of bar chart viewer should be visible
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of bar chart viewer should list 6 items
    When user sets properties of bar chart viewer:
      | stackColumnName | RACE |
      | Include Nulls   | true |
    Then the legend of bar chart viewer should list 4 items
    And no errors should have been logged

  Scenario: The panel's filter and the legends it left come back from a saved layout
    When user clicks on close icon of last histogram viewer
    And user adds a range filter on "WEIGHT" from 100 to 170
    And user clicks on the "category Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    Then 12 rows should pass the filter
    And the legend of pie chart viewer should list 2 items
    When user saves the layout of the current table view to the server
    And user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of pie chart viewer should list 4 items
    When user loads the saved layout
    Then 12 rows should pass the filter
    And the legend of scatter plot viewer should list 2 items
    And the legend of line chart viewer should list 2 items
    And the legend of bar chart viewer should list 2 items
    And the legend of pie chart viewer should list 2 items
    And the legend of trellis plot viewer should list 2 items
    And the legend of box plot viewer should list 2 items
    And "Caucasian" legend item in legend of scatter plot viewer should be absent
    And "xColumnName" property of scatter plot viewer should be "WEIGHT"
    And no errors should have been logged

  Scenario: A layout brings the panel's filter back and the filter a click set does not come back
    Then 12 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    And user adds a range filter on "WEIGHT" from 100 to 170
    Then 143 rows should pass the filter
    When user picks "On Click > Filter" from the context menu of bar chart viewer
    And user clicks on the "bar Black" area of bar chart viewer
    Then fewer than 143 rows should pass the filter
    And the legend of box plot viewer should list 1 item
    And "Black" legend item in legend of box plot viewer should be visible
    When user saves the layout of the current table view to the server
    And user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the legend of box plot viewer should list 4 items
    When user loads the saved layout
    Then 143 rows should pass the filter
    And the legend of scatter plot viewer should list 3 items
    And the legend of line chart viewer should list 3 items
    And the legend of pie chart viewer should list 3 items
    And the legend of trellis plot viewer should list 3 items
    And the legend of box plot viewer should list 3 items
    And "Black" legend item in legend of box plot viewer should be visible
    And "Caucasian" legend item in legend of box plot viewer should be visible
    When user clicks on the "category Black of RACE" area of filter panel
    And user clicks on the "checkbox Other of RACE" area of filter panel
    Then 12 rows should pass the filter
    When user sets properties of bar chart viewer:
      | On Click   | None     |
      | Row Source | Filtered |
    Then no errors should have been logged

  Scenario: On Click Filter from the pie chart's own menu comes with Row Source All
    Then 12 rows should pass the filter
    And the legend of pie chart viewer should list 2 items
    When user picks "On Click > Filter" from the context menu of pie chart viewer
    Then "On Click" property of pie chart viewer should be "Filter"
    And "Row Source" property of pie chart viewer should be "All"
    And the legend of pie chart viewer should list 4 items
    And 12 rows should pass the filter
    When user sets properties of pie chart viewer:
      | On Click   | None     |
      | Row Source | Filtered |
    Then the legend of pie chart viewer should list 2 items
    And 12 rows should pass the filter
    And no errors should have been logged

  Scenario: A histogram viewer joins the view and a layout is saved with it
    When user adds a histogram viewer
    And user sets properties of last histogram viewer:
      | valueColumnName   | AGE    |
      | splitColumnName   | RACE   |
      | Legend Visibility | Always |
      | Legend Position   | Right  |
    Then 12 rows should pass the filter
    When user saves the layout of the current table view to the server
    And user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And no errors should have been logged

  @known-failure
  Scenario: The layout that holds the histogram viewer brings the panel's filter back
    When user loads the saved layout
    Then 12 rows should pass the filter
