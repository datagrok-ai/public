@journey @serial @sticky-meta
Feature: Sticky metadata belongs to the molecule, not to the table that shows it
  Metadata written for two SPGI molecules is found again wherever those molecules are shown: in a
  clone of the table and a new view of it, in a CSV file the table was exported to and opened again,
  in a project saved and reopened, in the same project moved into a space and opened from there,
  after the page is reloaded and after another account has signed in on the page and the running
  one has come back. Removed with Clear values for the schema, it is gone, and stays gone after a
  reload; emptied in the cell's dialog and saved, it stays. Translated from TestTrack
  StickyMeta/03-persistence-copy-delete.md (Tests 3.1-3.4), its primary copy-clone-delete case, and
  playwright-tests/e2e/stickymeta/03-persistence-copy-delete.test.ts.

  On every surface the claim names the row first (the grid's text of its Id cell), then reads the
  cell's marker tooltip — the values the grid draws the marker for. The tooltip is built from a cache
  the page keeps per molecule, so after a fresh load the values are read from the cell's Sticky meta
  pane first, whose fields appear only once the values have come from the server; an absence ("No
  sticky meta for this cell") is read only after such a read. The old spec drove these operations
  through the JS API and read the values back with getAllValues; here each one is made the way a user
  makes it. The session the running account comes back with is the one it had: the library signs an
  account in with its token, and the old spec's logout through the login form cannot run where the
  suite is given a token, not a password (a fresh session of the same account is in MISSING.md). The
  case's export and import is a CSV downloaded from the table and opened again from the file (the old
  spec's d42 round trip ran in the page, not through a file). Not translated: a server restart, which
  a feature cannot cause.

  An emptied field saved in the cell's dialog keeps its value: GROK-15602 ("cannot remove values") was
  resolved with Clear values for the schema, so the case's "delete the fields, save" is claimed as the
  platform does it — the values stay, and Clear values removes them.

  The entity type and the schema are made through the UI in the Background and deleted through the UI
  in the last scenarios; the project and the space are removed by their steps when the feature ends.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And no project named "bdd-sm-keep-{time}" is on the server
    And no space named "bdd-sm-space-{time}" is on the server
    And the browse panel is open
    And the context panel is open
    When user expands "Platform" tree node inside browse tree
    And user expands "Platform > Sticky Meta" tree node inside browse tree
    And user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    And user clicks on "New Entity Type..." button
    And user enters "bdd-sm-keep-type-{time}" into "Name" input in "Create a new entity type" dialog
    And user enters "semtype=molecule" into "Matching expression" input in "Create a new entity type" dialog
    And user clicks on OK button in "Create a new entity type" dialog
    Then the "Create a new entity type" dialog should close
    When user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    And user clicks on "New Schema..." button
    And user enters "bdd-sm-keep-{time}" into "Name" input in "Create a new schema" dialog
    And user clicks on "select entities" action in "Create a new schema" dialog
    And user checks "bdd-sm-keep-type-{time}" property in "Select types for bdd-sm-keep-{time}" dialog
    And user clicks on OK button in "Select types for bdd-sm-keep-{time}" dialog
    And user enters "rating" into second "Name" input in "Create a new schema" dialog
    And user selects "int" in "Property Type" input in "Create a new schema" dialog
    And user clicks on "Add new property to schema" button in "Create a new schema" dialog
    And user enters "notes" into third "Name" input in "Create a new schema" dialog
    And user selects "string" in second "Property Type" input in "Create a new schema" dialog
    And user clicks on OK button in "Create a new schema" dialog
    Then the "Create a new schema" dialog should close
    Given user opens spgi dataset

  Scenario: Metadata is written for two molecules
    Then the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    And user enters "5" into Rating input in "Sticky meta" dialog
    And user enters "excellent" into Notes input in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 2 of Structure" area of grid
    And user enters "4" into Rating input in "Sticky meta" dialog
    And user enters "good" into Notes input in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    When user hovers over the "top right corner of cell 2 of Structure" area of grid
    Then the tooltip should show "rating" as "4"
    And the tooltip should show "notes" as "good"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A clone of the table and a new view of it show the same metadata (3.1)
    Given simple mode is off
    When user picks "Table > Clone" from the context menu of the current view tab
    Then table "spgi-100 (2)" should be open
    And the "spgi-100 (2)" view should be current
    And the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    When user picks "View > Layout > Clone View" from the top menu
    Then the current view should hold at least 1 viewer
    And the "spgi-100 (2) copy" view should be current
    When user hovers over the "top right corner of cell 2 of Structure" area of grid
    Then the tooltip should show "rating" as "4"
    And the tooltip should show "notes" as "good"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The table exported to a file and opened from it shows the same metadata (3.2)
    Given user watches downloads
    When user closes all views
    And user opens spgi dataset
    And user clicks on Export icon in toolbar
    And user downloads a file through "As CSV" text in toolbar
    Then the downloaded file should contain "Structure"
    When user closes all views
    And user uploads the downloaded file through "Open local file" icon inside browse toolbar
    Then the table should have 100 rows
    And "Structure" column should have semantic type "Molecule"
    And the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A project saved with the table and reopened shows the same metadata (3.2)
    When user closes all views
    And user opens spgi dataset
    And user opens the Save project dialog from the ribbon
    And user types "bdd-sm-keep-{time}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "bdd-sm-keep-{time}" should be on the server
    And the "spgi-100" table of the "bdd-sm-keep-{time}" project should be saved as a snapshot
    When user closes all views
    And user opens the "bdd-sm-keep-{time}" project and waits for its table
    Then the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The project moved into a space and opened from there shows the same metadata (3.2)
    When user closes all views
    And user picks "Create Space..." from the context menu of Spaces tree node inside browse tree
    And user enters "bdd-sm-space-{time}" into Name input in Create Space dialog
    And user clicks on OK button in Create Space dialog
    Then 1 space named "bdd-sm-space-{time}" should be on the server
    And the "Create Space" dialog should close
    When user refreshes the browse tree
    And "My stuff" tree node inside browse tree is expanded
    And "My stuff > My dashboards" tree node inside browse tree is expanded
    And user drags "My stuff > My dashboards > bdd-sm-keep-{time}" tree node inside browse tree to bdd-sm-space-{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on bdd-sm-space-{time} tree node inside browse tree
    Then the "bdd-sm-space-{time}" view should be current
    When user double-clicks on bdd-sm-keep-{time} link in gallery
    Then the project "bdd-sm-keep-{time}" should be open
    And the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The metadata is there after the page is reloaded (3.3)
    When user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in context panel should have value "4"
    And Notes input in context panel should have value "good"
    When user hovers over the "top right corner of cell 2 of Structure" area of grid
    Then the tooltip should show "rating" as "4"
    And the tooltip should show "notes" as "good"
    And no errors should have been logged

  Scenario: The metadata is there after another account has signed in on the page and the running one has come back (3.3)
    When user signs in as the sharing user
    Then the sharing user should be signed in
    When user signs in as themselves again
    Then the running account should be signed in
    When user opens spgi dataset
    And user clicks on the "cell 1 of Structure" area of grid
    Then the context panel should show the current cell
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in context panel should have value "5"
    And Notes input in context panel should have value "excellent"
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then the tooltip should show "rating" as "5"
    And the tooltip should show "notes" as "excellent"
    And no errors should have been logged

  Scenario: Clear values for the schema removes the cell's metadata, also after a reload (3.4)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    Then Rating input in "Sticky meta" dialog should have value "5"
    When user hovers over the text "sm-keep-{time}" in "Sticky meta" dialog
    And user clicks on "Clear values for the schema" icon in "Sticky meta" dialog
    Then "Are you sure?" dialog should be visible
    When user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in context panel should have value "4"
    When user clicks on the "cell 1 of Structure" area of grid
    Then the context panel should show the current cell
    And Rating input in context panel should have value ""
    And Notes input in context panel should have value ""
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then tooltip should contain text "No sticky meta for this cell"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # The cell's dialog does not delete a field that is emptied: GROK-15602 ("cannot remove values") was
  # resolved with Clear values for the schema, which the scenario before this one claims.
  Scenario: Fields emptied in the cell's dialog and saved keep their values; Clear values is the way to remove them (3.4)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 2 of Structure" area of grid
    Then Rating input in "Sticky meta" dialog should have value "4"
    When user clears Rating input in "Sticky meta" dialog
    And user clears Notes input in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "Sticky meta" dialog
    Then Save button in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in context panel should have value "4"
    And Notes input in context panel should have value "good"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: The schema is deleted
    When user presses Escape
    And user closes all views
    And the browse panel is open
    And user expands "Platform" tree node inside browse tree
    And user expands "Platform > Sticky Meta" tree node inside browse tree
    And user clicks on "Platform > Sticky Meta > Schemas" tree node inside browse tree
    And user types "bdd-sm-keep-{time}" into gallery search
    Then "bdd-sm-keep-{time}" link in gallery should be visible
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-keep-{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-keep-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged

  # GROK-18980 (Open) reported the values of a deleted schema still shown on the molecule until its
  # entity type was deleted too; on 2026-10-01 the pane showed none — this claims it stays that way.
  # The pane's "Semantic Type" header is built in the same pass as its schema sections, once every
  # value has been read (stickyMetaEditorForCell), so it is the sign the absence is read after the load.
  Scenario: The values of a deleted schema are no longer shown on the molecule
    When user closes all views
    And user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    Given "Sticky meta" pane in context panel is expanded
    Then "Sticky meta" pane in context panel should contain text "Semantic Type"
    And Rating input in context panel should be absent
    And Notes input in context panel should be absent
    And "Sticky meta" pane in context panel should not contain text "sm-keep-{time}"
    And no errors should have been logged

  Scenario: The entity type is deleted
    When user closes all views
    And the browse panel is open
    And "Platform" tree node inside browse tree is expanded
    And "Platform > Sticky Meta" tree node inside browse tree is expanded
    And user clicks on "Platform > Sticky Meta > Types" tree node inside browse tree
    And user types "bdd-sm-keep-type-{time}" into gallery search
    Then "bdd-sm-keep-type-{time}" link in gallery should be visible
    When user remembers the gallery counter
    And user picks "Delete" from the context menu of "bdd-sm-keep-type-{time}" link in gallery
    And user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And the gallery counter should be lower than remembered
    And "bdd-sm-keep-type-{time}" link in gallery should be absent
    When user clears gallery search
    Then no errors should have been logged
