@journey @serial @sticky-meta
Feature: Sticky metadata belongs to the molecule, not to the table that shows it
  Metadata written for two SPGI molecules is found again wherever those molecules are shown: in a
  clone of the table and a new view of it, in a project saved and reopened, in the same project moved
  into a space and opened from there, and after the page is reloaded. Removed with Clear values for the
  schema, it is gone; emptied in the cell's dialog and saved, it stays. A CSV round trip, a
  non-canonical SMILES, and what Clear values and a deleted schema leave on the server are ApiTests'
  and datlas' claims (ApiTests `Dapi: sticky meta: molecules`, datlas `test/sticky_meta`). Translated from TestTrack
  StickyMeta/03-persistence-copy-delete.md (Tests 3.1-3.4), its primary copy-clone-delete case, and
  playwright-tests/e2e/stickymeta/03-persistence-copy-delete.test.ts.

  On every surface the claim names the row first (the grid's text of its Id cell), then reads the
  cell's marker tooltip — the values the grid draws the marker for. The tooltip is built from a cache
  the page keeps per molecule, so after a fresh load the values are read from the cell's Sticky meta
  pane first, whose fields appear only once the values have come from the server; an absence ("No
  sticky meta for this cell") is read only after such a read. The old spec drove these operations
  through the JS API and read the values back with getAllValues; here each one is made the way a user
  makes it. Another account signing in and the running one coming back (3.3) is not claimed: the
  library signs the account back in with the session it had, which proves nothing a reload does not
  (a fresh session of the same account is in MISSING.md). Not translated: a server restart, which a
  feature cannot cause.

  An emptied field saved in the cell's dialog keeps its value: GROK-15602 ("cannot remove values") was
  resolved with Clear values for the schema, so the case's "delete the fields, save" is claimed as the
  platform does it — the values stay, and Clear values removes them.

  The cell's dialog and pane list one section per schema that matches molecules, each named after its
  schema, so the fields are named in the feature's own section and another schema matching molecules on
  the stand (dev has four) is no concern.

  The entity type and the schema are made through the API in the Background and removed at feature end
  (the New Entity Type and New Schema dialogs and the galleries' Delete are schema-and-type's subject); the
  schema alone is deleted before the last scenario, its entity type kept, as GROK-18980 describes. The
  project and the space are removed by their steps when the feature ends.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And no project named "bdd-sm-keep-{time}" is on the server
    And no space named "bdd-sm-space-{time}" is on the server
    And a Sticky Meta schema "bdd-sm-keep-{time}" for molecules, with the entity type "bdd-sm-keep-type-{time}", is on the server with the properties:
      | rating | int    |
      | notes  | string |
    And the browse panel is open
    And the context panel is open
    Given user opens spgi dataset

  Scenario: Metadata is written for two molecules
    Then the "text of cell 1 of Id" reading of grid should be "CAST-634783"
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    And user enters "5" into Rating input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    And user enters "excellent" into Notes input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 2 of Structure" area of grid
    And user enters "4" into Rating input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    And user enters "good" into Notes input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be disabled
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
    And the "current row" reading of grid should be 2
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in "bdd-sm-keep-{time}" section in context panel should have value "4"
    And Notes input in "bdd-sm-keep-{time}" section in context panel should have value "good"
    When user hovers over the "top right corner of cell 2 of Structure" area of grid
    Then the tooltip should show "rating" as "4"
    And the tooltip should show "notes" as "good"
    And no errors should have been logged

  Scenario: Clear values for the schema removes the cell's metadata (3.4)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 1 of Structure" area of grid
    Then Rating input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should have value "5"
    When user hovers over the text "sm-keep-{time}" in "Sticky meta" dialog
    And user clicks on "Clear values for the schema" icon in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then "Are you sure?" dialog should be visible
    When user clicks on YES button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    And the "current row" reading of grid should be 2
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in "bdd-sm-keep-{time}" section in context panel should have value "4"
    When user clicks on the "cell 1 of Structure" area of grid
    Then the context panel should show the current cell
    And the "current row" reading of grid should be 1
    And Rating input in "bdd-sm-keep-{time}" section in context panel should have value ""
    And Notes input in "bdd-sm-keep-{time}" section in context panel should have value ""
    When user hovers over the "top right corner of cell 1 of Structure" area of grid
    Then tooltip should contain text "No sticky meta for this cell"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # The cell's dialog does not delete a field that is emptied: GROK-15602 ("cannot remove values") was
  # resolved with Clear values for the schema, which the scenario before this one claims.
  Scenario: Fields emptied in the cell's dialog and saved keep their values; Clear values is the way to remove them (3.4)
    When user picks "Sticky meta > Edit for current cell..." from the context menu of the "cell 2 of Structure" area of grid
    Then Rating input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should have value "4"
    When user clears Rating input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    And user clears Notes input in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be enabled
    When user clicks on Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog
    Then Save button in "bdd-sm-keep-{time}" section in "Sticky meta" dialog should be disabled
    When user clicks on CANCEL button in "Sticky meta" dialog
    Then the "Sticky meta" dialog should close
    When user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    And the "current row" reading of grid should be 2
    Given "Sticky meta" pane in context panel is expanded
    Then Rating input in "bdd-sm-keep-{time}" section in context panel should have value "4"
    And Notes input in "bdd-sm-keep-{time}" section in context panel should have value "good"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  # GROK-18980 (Open) reported the values of a deleted schema still shown on the molecule until its
  # entity type was deleted too; on 2026-10-01 the pane showed none — this claims it stays that way.
  # The pane's "Semantic Type" header is built in the same pass as its schema sections, once every
  # value has been read (stickyMetaEditorForCell), so it is the sign the absence is read after the load.
  Scenario: The values of a deleted schema are no longer shown on the molecule
    Given the Sticky Meta schema "bdd-sm-keep-{time}" is deleted on the server
    When user closes all views
    And user reloads the page
    And user opens spgi dataset
    And user clicks on the "cell 2 of Structure" area of grid
    Then the context panel should show the current cell
    And the "current row" reading of grid should be 2
    Given "Sticky meta" pane in context panel is expanded
    Then "Sticky meta" pane in context panel should contain text "Semantic Type"
    And Rating input in context panel should be absent
    And Notes input in context panel should be absent
    And "Sticky meta" pane in context panel should not contain text "sm-keep-{time}"
    And no errors should have been logged
