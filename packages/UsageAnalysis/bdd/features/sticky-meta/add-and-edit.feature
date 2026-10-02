@journey @serial @sticky-meta
Feature: Sticky metadata on a molecule cell, as sticky columns and on several rows at once
  A schema that matches molecules gives every Structure cell of SPGI its own metadata: written for one
  cell in the cell's Sticky meta dialog and read back from the cell's marker tooltip and from the
  dialog reopened; materialized as sticky columns from the column's Sticky meta pane, which sort like
  any column and come back with their values after one is removed; and written for several selected
  rows at once from Edit for all properties. Translated from TestTrack StickyMeta/02-add-and-edit.md
  (Tests 2.1-2.3), its primary add-and-edit case, and
  playwright-tests/e2e/stickymeta/02-add-and-edit.test.ts.

  The marker itself is a dot the grid paints on its canvas and reports nowhere; the tooltip shown
  over the cell's top right corner is built from the same values the dot is drawn for, so the values
  are claimed there. The blue circle on a sticky column's header is painted the same way and is not
  claimed. The batch edit is claimed in the sticky columns and again after the table is reopened, so
  what reached the server is read, not what the view keeps.

  The cell's dialog and pane list one section per schema that matches molecules, and a section has no
  container of its own, so its fields are named across all of them: the feature needs a stand where its
  schema is the only one matching molecules (dev has four more — MISSING.md, section 8).

  The entity type and the schema are made through the UI in the Background and deleted through the UI
  in the last scenarios, which a journey runs even after a failure; the server-side sweep the library
  does not have yet is in sticky-meta/MISSING.md. Serial: a schema matching molecules adds a section
  to the Sticky meta pane of every molecule column on the stand.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Sticky Meta" tree node inside browse tree
    And user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    And user clicks on "New Entity Type..." button
    And user enters "bdd-sm-cells-type-{time}" into "Name" input in "Create a new entity type" dialog
    And user enters "semtype=molecule" into "Matching expression" input in "Create a new entity type" dialog
    And user clicks on OK button in "Create a new entity type" dialog
    Then the "Create a new entity type" dialog should close
    When user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    And user clicks on "New Schema..." button
    And user enters "bdd-sm-cells-{time}" into "Name" input in "Create a new schema" dialog
    And user clicks on "select entities" action in "Create a new schema" dialog
    And user checks "bdd-sm-cells-type-{time}" property in "Select types for bdd-sm-cells-{time}" dialog
    And user clicks on OK button in "Select types for bdd-sm-cells-{time}" dialog
    And user enters "rating" into second "Name" input in "Create a new schema" dialog
    And user selects "int" in "Property Type" input in "Create a new schema" dialog
    And user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "notes" into third "Name" input in "Create a new schema" dialog
    And user selects "string" in second "Property Type" input in "Create a new schema" dialog
    And user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "verified" into fourth "Name" input in "Create a new schema" dialog
    And user selects "bool" in third "Property Type" input in "Create a new schema" dialog
    And user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "review_date" into fifth "Name" input in "Create a new schema" dialog
    And user selects "datetime" in fourth "Property Type" input in "Create a new schema" dialog
    And user clicks on OK button in "Create a new schema" dialog
    Then the "Create a new schema" dialog should close
    Given user opens spgi dataset

  Scenario: Metadata written for one cell shows in its tooltip and stays in its dialog (2.1)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    Then "Sticky meta" dialog should be visible
    When user enters "5" into Rating input in "Sticky meta" dialog
    And user enters "test note" into Notes input in "Sticky meta" dialog
    And user checks Verified input in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "test note"
    And the tooltip should show "verified" as "true"
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    Then Rating input in "Sticky meta" dialog should have value "5"
    And Notes input in "Sticky meta" dialog should have value "test note"
    And Verified input in "Sticky meta" dialog should be checked
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user hovers over the "top right corner of cell 2 of Structure" area of grid
    Then tooltip should contain text "No sticky meta for this cell"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The column's Sticky meta pane adds the properties as sticky columns (2.2)
    When user clicks on the "header Structure" area of grid
    Then the context panel should show "Structure"
    Given "Sticky meta" pane in context panel is expanded
    When user clicks on "Add bdd-sm-cells-{time}'s properties as columns" button in context panel
    Then the table should have a column "rating"
    And the table should have a column "notes"
    And the table should have a column "verified"
    And the table should have a column "review_date"
    And the value of "rating" column in row 1 should be "5"
    And the value of "notes" column in row 1 should be "test note"
    And the value of "verified" column in row 1 should be "true"
    And the value of "rating" column in row 2 should be ""
    And the value of "notes" column in row 2 should be ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A sticky column sorts like any column (2.2)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 2 of Structure" area of grid
    And user enters "9" into Rating input in "Sticky meta" dialog
    And user clicks on Save button in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    And the value of "rating" column in row 2 should be "9"
    When user double-clicks on the "header rating" area of grid
    Then the "sort column" reading of grid should be "rating"
    And the "sort direction" reading of grid should be "descending"
    And the "row order" reading of grid should be "2, 1, 3, 4, 5, 6, 7, 8, 9, 10"
    And no errors should have been logged

  Scenario: A removed sticky column comes back with its values (2.2)
    When user picks "Remove" from the context menu of the "header rating" area of grid
    Then the table should not have a column "rating"
    When user clicks on the "header Structure" area of grid
    Then the context panel should show "Structure"
    Given "Sticky meta" pane in context panel is expanded
    When user clicks on "Add rating as a column" button in context panel
    Then the table should have a column "rating"
    And the value of "rating" column in row 1 should be "5"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Edit for all properties writes the same values into every selected row (2.3)
    When user selects the first 3 rows
    Then rows 1 to 3 should be selected
    When user picks "Sticky meta > Edit for all properties" from the context menu of the "cell 2 of Structure" area of grid
    Then Rows input should have value "Selected"
    When user checks "verified" property
    And user clicks on editor of "notes" property
    And user types "batch note" at the caret
    And user presses Enter
    And user switches to the "spgi-100" table view
    Then the value of "notes" column in row 1 should be "batch note"
    And the value of "notes" column in row 2 should be "batch note"
    And the value of "notes" column in row 3 should be "batch note"
    And the value of "verified" column in row 2 should be "true"
    And the value of "verified" column in row 3 should be "true"
    And the value of "notes" column in row 4 should be ""
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: What was written is read back from the server after the table is reopened
    When user closes all views
    And user opens spgi dataset
    And user clicks on the "header Structure" area of grid
    Then the context panel should show "Structure"
    Given "Sticky meta" pane in context panel is expanded
    When user clicks on "Add bdd-sm-cells-{time}'s properties as columns" button in context panel
    Then the table should have a column "rating"
    And the table should have a column "notes"
    And the table should have a column "verified"
    And the value of "rating" column in row 1 should be "5"
    And the value of "notes" column in row 1 should be "batch note"
    And the value of "notes" column in row 2 should be "batch note"
    And the value of "verified" column in row 2 should be "true"
    And the value of "notes" column in row 3 should be "batch note"
    And the value of "verified" column in row 3 should be "true"
    And the value of "notes" column in row 4 should be ""
    And no errors should have been logged

  Scenario: The schema is deleted
    When user closes all views
    And user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    And user types "bdd-sm-cells-{time}" into gallery search
    Then "bdd-sm-cells-{time}" link in gallery should be visible
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-cells-{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-cells-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged

  Scenario: The entity type is deleted
    When user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    And user types "bdd-sm-cells-type-{time}" into gallery search
    Then "bdd-sm-cells-type-{time}" link in gallery should be visible
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-cells-type-{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-cells-type-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
