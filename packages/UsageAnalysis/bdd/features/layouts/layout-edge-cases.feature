@realizes:views.layouts
Feature: Layout edge cases: the Layouts pane suggests a layout only for a table that has its columns
  A layout is suggested for a table when its "columns of interest" match a subset of the table's
  columns, and for a layout with non-grid viewers those are the data columns its viewers use
  (`help/visualize/view-layout.md`, "Layout suggestions"). A layout whose scatter plot plots HEIGHT
  against WEIGHT is saved from Toolbox > Layouts on a demog table, and a layout whose scatter plot
  plots peak.rpm against city.mpg on a cars table. Opened again, the cars table is suggested its
  own layout and not the demog one, which needs columns cars does not have; the demog table opened
  again is suggested the demog layout. The cars card shown in the same pane is what makes the
  missing demog card a finished reading rather than a pane still loading.

  A layout of a grid alone takes the grid's columns as its columns of interest
  (`TableView.getViewerColumns`), so it is suggested only for a table that has them (fixed in
  core#1179; a known failure until the stand has it). A layout
  renamed to an empty name in the context panel keeps its name.

  The tables are named with the run's time, so each layout is named after its table and is unique to
  the run; both layouts are deleted when the feature ends.

  Background:
    Given user is logged in
    And the layouts named "bdd-edge-layout-{time}" are deleted when the feature ends
    And the layouts named "bdd-edge-cars-{time}" are deleted when the feature ends
    And the layouts named "bdd-edge-grid-{time}" are deleted when the feature ends
    And the layouts named "bdd-edge-carsb-{time}" are deleted when the feature ends
    And the layouts named "bdd-edge-name-{time}" are deleted when the feature ends
    And the toolbox pane is shown
    And the context panel is open

  @edge-cases
  Scenario: A layout plotting HEIGHT and WEIGHT is not suggested for a table without those columns
    Given user opens demog dataset keeping the first 100 rows as "bdd-edge-layout-{time}"
    And user adds a scatter plot viewer with:
      | xColumnName | HEIGHT |
      | yColumnName | WEIGHT |
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-layout-{time}" layout card should be visible
    When user closes all views
    Given user opens cars dataset keeping the first 30 rows as "bdd-edge-cars-{time}"
    And user adds a scatter plot viewer with:
      | xColumnName | peak.rpm |
      | yColumnName | city.mpg |
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header model" area of grid
    Then the context panel should show "model"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-cars-{time}" layout card should be visible
    When user closes all views
    Given user opens cars dataset keeping the first 30 rows as "bdd-edge-cars-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "bdd-edge-cars-{time}" layout card should be visible
    And "bdd-edge-layout-{time}" layout card should be absent
    When user closes all views
    Given user opens demog dataset keeping the first 100 rows as "bdd-edge-layout-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "bdd-edge-layout-{time}" layout card should be visible
    And "bdd-edge-cars-{time}" layout card should be absent
    And no errors should have been logged

  @edge-cases @known-failure:85 @realizes:core#1179
  Scenario: A layout of a grid alone is suggested only for a table that has the grid's columns
    Before core#1179 the grid-only layout is suggested for the cars table too: `getViewerColumns`
    looked the grid's columns up by `name` instead of `colName`. Remove the tag once the stand has
    core#1179 and the cars table is not offered it.
    Given user opens demog dataset keeping the first 100 rows as "bdd-edge-grid-{time}"
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-grid-{time}" layout card should be visible
    When user closes all views
    Given user opens cars dataset keeping the first 30 rows as "bdd-edge-carsb-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header model" area of grid
    Then the context panel should show "model"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-carsb-{time}" layout card should be visible
    And "bdd-edge-grid-{time}" layout card should be absent
    When user closes all views
    Given user opens demog dataset keeping the first 100 rows as "bdd-edge-grid-{time}"
    And the toolbox pane is shown
    And Layouts accordion header in toolbox is expanded
    Then "bdd-edge-grid-{time}" layout card should be visible
    And no errors should have been logged

  @edge-cases
  Scenario: A layout renamed to an empty name in the context panel keeps its name
    Given user opens demog dataset keeping the first 100 rows as "bdd-edge-name-{time}"
    And Layouts accordion header in toolbox is expanded
    When user clicks on the "header AGE" area of grid
    Then the context panel should show "AGE"
    When user clicks on Save button in layouts pane
    Then "bdd-edge-name-{time}" layout card should be visible
    And the context panel should show "bdd-edge-name-{time}"
    When user clears Name details field in context panel
    And user presses Enter
    Then "bdd-edge-name-{time}" layout card should be visible
    And value of Name details field in context panel should have text "bdd-edge-name-{time}"
    And no errors should have been logged
