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

  What the pane suggests for a layout of a grid alone is not settled (see MISSING.md), so it is not
  claimed here.

  The tables are named with the run's time, so each layout is named after its table and is unique to
  the run; both layouts are deleted when the feature ends.

  Background:
    Given user is logged in
    And the layouts named "bdd-edge-layout-{time}" are deleted when the feature ends
    And the layouts named "bdd-edge-cars-{time}" are deleted when the feature ends
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
