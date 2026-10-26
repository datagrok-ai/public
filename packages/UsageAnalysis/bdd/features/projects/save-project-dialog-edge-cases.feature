@realizes:views.projects
Feature: Save project dialog edge cases: CANCEL and Escape, names without a letter, a local file, a View and use recipient
  Guards for the ribbon's Save project dialog. CANCEL and Escape on the dialog of a new dashboard
  close it and save nothing: no Share dialog follows and the server has no project of the typed
  name. A name must contain at least one letter (`isValidEntityName` in the Save dialog,
  `project_meta.dart`): an emptied name, spaces only, digits only or a "/" alone mark the Name box
  invalid and disable OK, and a letter added to the same text enables OK again. A "/" inside a name
  with letters is taken: the project is saved under that name and opens from the address the page
  shows after the save (every dashboard has a URL, `dashboard.md`, "Sharing by link").

  A table opened with Open local file has no creation script, so the dialog offers it no Data sync
  switch and no Creation script button: "the only table that can't be saved with Data sync is one
  opened from a local file" (`help/datagrok/concepts/project/dashboard.md`, "Data sync").

  A recipient with View and use cannot save the original project and can still save a copy
  (`dashboard.md`, "Saving a dashboard"). The owner's own dialog for the saved project opens on
  the original, its Name box holding the project's name; the recipient's opens in its Save a copy
  mode, the Name box holding "Copy of <name>" (the name the dialog gives only in that mode). Which
  option of this Dart radio group is checked or disabled has no reading, and a click on the
  disabled option is refused by the browser before it lands (see MISSING.md), so the mode is read
  off the name.

  Data: cars.csv from Browse > Files > Demo; `fixtures/browse-import.csv` for the local file. The
  project names carry the run's time (letters and digits, the "/" one aside); every project of
  these names is removed with its table, view and grant when a scenario starts and when the
  feature ends. The recipient is the library's sharing user.

  Background:
    Given user is logged in
    And the context panel is open
    And the browse panel is open
    And no project named "BDDEdgeSave{time}" is on the server
    And no project named "BDDEdge/Slash{time}" is on the server

  @edge-cases
  Scenario: CANCEL on the Save dialog of a new dashboard saves nothing
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDEdgeSave{time}" into Name text input in "Save project" dialog
    And user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSave{time}" dialog should be absent
    And 0 projects named "BDDEdgeSave{time}" should be on the server
    And no errors should have been logged

  @edge-cases
  Scenario: Escape on the Save dialog of a new dashboard saves nothing
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDEdgeSave{time}" into Name text input in "Save project" dialog
    And user presses Escape
    Then the "Save project" dialog should close
    And "Share BDDEdgeSave{time}" dialog should be absent
    And 0 projects named "BDDEdgeSave{time}" should be on the server
    And no errors should have been logged

  @edge-cases
  Scenario: An emptied name is marked invalid and OK is disabled
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And OK button in "Save project" dialog should be enabled
    When user clears Name text input in "Save project" dialog
    Then Name text input in "Save project" dialog should be invalid
    And OK button in "Save project" dialog should be disabled
    When user types "x" into Name text input in "Save project" dialog
    Then Name text input in "Save project" dialog should be valid
    And OK button in "Save project" dialog should be enabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  @edge-cases
  Scenario Outline: A name without a letter ("<name>") is marked invalid and OK is disabled
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user types "<name>" into Name text input in "Save project" dialog
    Then Name text input in "Save project" dialog should be invalid
    And OK button in "Save project" dialog should be disabled
    When user types "<name>x" into Name text input in "Save project" dialog
    Then Name text input in "Save project" dialog should be valid
    And OK button in "Save project" dialog should be enabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

    Examples:
      | name |
      | 123  |
      | /    |
      | 1/2  |

  @edge-cases
  Scenario: A name of spaces only is marked invalid and OK is disabled
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user types "   " into Name text input in "Save project" dialog
    Then Name text input in "Save project" dialog should be invalid
    And OK button in "Save project" dialog should be disabled
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close

  @edge-cases
  Scenario: A name with a "/" between letters is saved and the project opens from its address
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDEdge/Slash{time}" into Name text input in "Save project" dialog
    Then OK button in "Save project" dialog should be enabled
    When user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And 1 project named "BDDEdge/Slash{time}" should be on the server
    And "Share BDDEdge/Slash{time}" dialog should be visible
    When user clicks on CANCEL button in "Share BDDEdge/Slash{time}" dialog
    Then the "Share BDDEdge/Slash{time}" dialog should close
    And the page address should contain "/p/"
    When user remembers the page address
    And user closes all views
    And user opens the remembered address
    Then the "cars" view should be current
    And the table should have 30 rows
    And no errors should have been logged

  @edge-cases
  Scenario: A table opened from a local file is offered no Data sync
    When user uploads "fixtures/browse-import.csv" through "Open local file" icon inside browse toolbar
    Then the "browse-import" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And "browse-import" project table in "Save project" dialog should be visible
    And Data sync switch in "browse-import" project table in "Save project" dialog should be hidden
    And "Creation script" button in "browse-import" project table in "Save project" dialog should be hidden
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    And no errors should have been logged

  @edge-cases @realizes:sharing.share-dialog
  Scenario: A View and use recipient's Save dialog opens in its Save a copy mode
    Given Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    When user double-clicks Files---Demo---cars.csv tree node inside browse tree
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    When user enters "BDDEdgeSave{time}" into Name text input in "Save project" dialog
    And user clicks on OK button in "Save project" dialog
    Then the "Save project" dialog should close
    And "Share BDDEdgeSave{time}" dialog should be visible
    And share access selector should contain text "View and use"
    When user picks the sharing user in "User, group, or email" input in "Share BDDEdgeSave{time}" dialog
    And user unchecks "Send notifications" input in "Share BDDEdgeSave{time}" dialog
    And user clicks on OK button in "Share BDDEdgeSave{time}" dialog
    Then the "Share BDDEdgeSave{time}" dialog should close
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Name text input in "Save project" dialog should have value "BDDEdgeSave{time}"
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user closes all views
    And user signs in as the sharing user
    Then the sharing user should be signed in
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeSave{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user double-clicks on BDDEdgeSave{time} gallery card and the project opens
    Then the "cars" view should be current
    When user clicks on Save ribbon item
    Then "Save project" dialog should be visible
    And Name text input in "Save project" dialog should have value "Copy of BDDEdgeSave{time}"
    When user clicks on CANCEL button in "Save project" dialog
    Then the "Save project" dialog should close
    When user closes all views
    And user signs in as themselves again
    Then the running account should be signed in
