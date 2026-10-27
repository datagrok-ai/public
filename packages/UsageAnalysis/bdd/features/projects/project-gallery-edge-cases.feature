@realizes:views.projects
Feature: Project gallery edge cases: Rename with no letter, CANCEL on Rename and Delete, the address of a deleted project, a recipient's card menu
  Guards for a project's card in Browse > Dashboards. cars.csv from Browse > Files > Demo is saved
  through the ribbon's Save dialog, the address the page shows after the save is remembered, and the
  card is found by the gallery search.

  - The Rename dialog refuses a name without a letter (its validator, `DbEntityMeta.getRenameDialog`
    in `db_entity_meta.dart`): the Name box is marked invalid and OK is disabled; a letter typed in
    enables OK again. The validator's text, "Name must contain at least one letter", is not shown in
    the dialog's text, so it is not claimed.
  - CANCEL on Rename and on the "Are you sure?" dialog of Delete Project changes nothing: the server
    keeps the project under its name, and the gallery still shows its card.
  - The address of a deleted project names nothing, and opening it says so and leaves the shell as it
    was: an error balloon "Unable to get project asset", the current view kept (the contract
    `browse/browse-routing.feature` states for an address that names nothing).
  - A recipient with View and use may not rename or delete the project, so its card menu shows
    Rename... and Delete Project disabled (core#1162: today both are enabled, a known failure).

  The project names carry the run's time (letters and digits only: the Dashboards search misses "-"
  and "_"); every project of these names is removed with its table and view when a scenario starts
  and when the feature ends.

  Background:
    Given user is logged in
    And the context panel is open
    And the browse panel is open
    And no project named "BDDEdgeGal{time}" is on the server
    And no project named "BDDEdgeGalNew{time}" is on the server
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDEdgeGal{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeGal{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdgeGal{time}" dialog
    Then the "Share BDDEdgeGal{time}" dialog should close
    And the page address should contain "/p/"
    When user remembers the page address
    And user picks "Close All" from the context menu of browse tab
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeGal{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDEdgeGal{time} gallery card should be visible

  @edge-cases
  Scenario: Rename with the name emptied says a letter is needed and disables OK
    When user picks "Rename..." from the context menu of BDDEdgeGal{time} gallery card
    Then Rename project dialog should be visible
    And OK button in Rename project dialog should be enabled
    When user clears Name input in Rename project dialog
    Then Name input in Rename project dialog should be invalid
    And OK button in Rename project dialog should be disabled
    When user types "x" into Name input in Rename project dialog
    Then Name input in Rename project dialog should be valid
    And OK button in Rename project dialog should be enabled
    When user clicks on CANCEL button in Rename project dialog
    Then the "Rename project" dialog should close
    And 1 project named "BDDEdgeGal{time}" should be on the server

  @edge-cases
  Scenario: CANCEL on Rename keeps the project's name
    When user picks "Rename..." from the context menu of BDDEdgeGal{time} gallery card
    Then Rename project dialog should be visible
    When user enters "BDDEdgeGalNew{time}" into Name input in Rename project dialog
    And user clicks on CANCEL button in Rename project dialog
    Then the "Rename project" dialog should close
    And 0 projects named "BDDEdgeGalNew{time}" should be on the server
    And 1 project named "BDDEdgeGal{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDEdgeGal{time} gallery card should be visible

  @edge-cases
  Scenario: CANCEL on the Delete Project confirmation keeps the project
    When user picks "Delete Project" from the context menu of BDDEdgeGal{time} gallery card
    Then "Are you sure?" dialog should contain text 'Delete project "BDDEdgeGal{time}"?'
    When user clicks on CANCEL button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 1 project named "BDDEdgeGal{time}" should be on the server
    When user clicks on "Refresh" icon inside gallery toolbar
    Then BDDEdgeGal{time} gallery card should be visible

  @edge-cases
  Scenario: The address of a deleted project says it names nothing and keeps the current view
    When user picks "Delete Project" from the context menu of BDDEdgeGal{time} gallery card
    Then "Are you sure?" dialog should be visible
    When user clicks on DELETE button in "Are you sure?" dialog
    Then the "Are you sure?" dialog should close
    And 0 projects named "BDDEdgeGal{time}" should be on the server
    Given user opens demog-1000 dataset
    Then the "demog-1000" view should be current
    When user opens the remembered address
    Then an error balloon containing "Unable to get project asset" should have been shown
    And the "demog-1000" view should be current
    And no errors should have been logged

  @edge-cases @realizes:sharing.share-dialog @known-failure:119 @realizes:core#1162
  Scenario: A View and use recipient sees Rename and Delete Project disabled in the card menu
    core#1162: both items are enabled for a View and use recipient. Remove the tag when they are
    disabled.
    When user picks "Share..." from the context menu of BDDEdgeGal{time} gallery card
    Then "Share BDDEdgeGal{time}" dialog should be visible
    When user picks the sharing user in "User, group, or email" input in "Share BDDEdgeGal{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDEdgeGal{time}" dialog
    And user clicks on OK button in "Share BDDEdgeGal{time}" dialog
    Then the "Share BDDEdgeGal{time}" dialog should close
    When user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeGal{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    Then BDDEdgeGal{time} gallery card should be visible
    When user opens the context menu of BDDEdgeGal{time} gallery card
    Then the open menu should list "Rename..." as disabled
    And the open menu should list "Delete Project" as disabled
    When user presses Escape
    And user signs in as themselves again
    Then the running account should be signed in
