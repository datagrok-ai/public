@journey @serial @spaces @realizes:views.space
Feature: Putting files into a space by dragging them
  A file dragged from the demo files onto a space's node in the Browse tree, and what the
  "Move entity" dialog does with it. Translated from
  files/TestTrack/Spaces/spaces-general.test.ts (tests 9, 10, 11 and 16).

  Every drag out of Demo uses Copy, never Link — the rule the old suite arrived at, because a
  linked entity that is later moved is removed from the demo files for everyone. Move is exercised
  between two spaces, on a copy, where it costs nothing.

  The old spec reached the files with page.goto(BASE + '/files/System.DemoFiles/?browse=files'),
  which reloads the whole client; here the tree gets there, the way a person does. The two spaces are
  made through the API: the Create Space dialog is spaces-create's subject.

  When the drag starts inside a space view rather than in the files, the dialog is shown before its
  Link/Copy/Move chooser is in it, so the chooser is waited for as its own claim. That lateness is
  what the old spec's "add waitForTimeout(600) before clicking YES (let the dialog initialize)" was
  covering up.

  A link follows its original: a project moved into one space, linked into a second and then moved
  from the first into a third still shows its details when its link is clicked, with no error and no
  refresh of the tree (core#382, claimed on a project: a file's Link is core#1042). A function of a
  package (here a script of BiostructureViewer) can only be linked: dropped on a space it is linked at
  once, with no "Move entity" dialog (`ProjectMeta.doDropEx`). A file linked from one space into
  another is expected to show in the target space and the dialog to close; today the Link fails with
  "FormatException: SyntaxError: Unexpected end of JSON input" and the dialog stays open (core#1042),
  a known failure in the journey's last scenario, whose tag goes when core#1042 is fixed. These use
  spaces, a file (the project's fixtures/cars-small.csv) and a project of their own, named with the
  run's time and swept with what they hold when the feature starts and ends; the package script stays
  where it is.

  Serial: the package script is found by the Scripts view's search and the project by the Dashboards
  search, both the account's own settings, and the feature drags onto the Spaces tree with the other
  serial features of Spaces.

  Background:
    Given user is logged in
    And the browse panel is open
    And a space named "BDD-DnD" is on the server
    And a space named "BDD-DnD-Src" is on the server
    And a space named "BDDEdgeSrc{time}" is on the server
    And a space named "BDDEdgeDst{time}" is on the server
    And the "fixtures/cars-small.csv" file of the project is in the space "BDDEdgeSrc{time}" as "edgecars{time}.csv"
    And Spaces tree node inside browse tree is expanded

  Scenario: The demo files
    Given Files tree node inside browse tree is expanded
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The dialog offers Link, Copy and Move, and names the target
    When user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree
    Then Move entity dialog should be visible
    And Move entity dialog should contain text "BDDDnD"
    And choice input in Move entity dialog should have value "Link"

  Scenario: Cancelling leaves the space empty
    When user clicks on CANCEL button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on Spaces---BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And demog.csv link in gallery should be absent

  Scenario: A copied file lands in the space
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags demog.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on Spaces---BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And demog.csv link in gallery should be visible

  Scenario: The file's details show in the context panel
    Given the context panel is open
    When user clicks on demog.csv link in gallery
    Then the context panel should show "demog.csv"
    And context panel should contain text "demog"

  Scenario: A second file joins the first
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags TSLA.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on Spaces---BDD-DnD tree node inside browse tree
    Then demog.csv link in gallery should be visible
    And TSLA.csv link in gallery should be visible

  Scenario: A copy is made in the source space to move later
    When user clicks on "Files > Demo" tree node inside browse tree
    And user drags beer.csv link in gallery to Spaces---BDD-DnD-Src tree node inside browse tree
    And user selects "Copy" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree
    Then the "BDD-DnD-Src" view should be current
    And beer.csv link in gallery should be visible

  Scenario: Moving takes the file out of the source space, and a linked project follows its original moved to a third space
    When user drags beer.csv link in gallery to Spaces---BDD-DnD tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on Spaces---BDD-DnD tree node inside browse tree
    Then the "BDD-DnD" view should be current
    And beer.csv link in gallery should be visible
    When user double-clicks on Spaces---BDD-DnD-Src tree node inside browse tree
    Then the "BDD-DnD-Src" view should be current
    And beer.csv link in gallery should be absent
    Given no project named "BDDEdgeMoved{time}" is on the server
    And a space named "BDDEdgeMoveTo{time}" is on the server
    And user opens demog-1000 dataset
    When user saves the current view as project "BDDEdgeMoved{time}"
    And user closes all views
    Given the browse panel is open
    When user clicks on Dashboards tree node inside browse tree
    And user enters "BDDEdgeMoved{time}" into gallery search
    And user clicks on "Refresh" icon inside gallery toolbar
    And user picks "Move to Space..." from the context menu of BDDEdgeMoved{time} gallery card
    Then "Move to space" dialog should be visible
    When user selects "BDDEdgeSrc{time}" in Space input in "Move to space" dialog
    And user clicks on OK button in "Move to space" dialog
    Then the "Move to space" dialog should close
    When user clears gallery search
    Given Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree
    Then "BDDEdgeMoved{time}" link in gallery should be visible
    When user drags "BDDEdgeMoved{time}" link in gallery to BDDEdgeDst{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user drags "BDDEdgeMoved{time}" link in gallery to BDDEdgeMoveTo{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should be visible
    When user selects "Move" in Move entity dialog
    And user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    When user double-clicks on BDDEdgeMoveTo{time} tree node inside browse tree
    Then the "BDDEdgeMoveTo{time}" view should be current
    And "BDDEdgeMoved{time}" link in gallery should be visible
    When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree
    Then the "BDDEdgeSrc{time}" view should be current
    And "edgecars{time}.csv" link in gallery should be visible
    And "BDDEdgeMoved{time}" link in gallery should be absent
    Given the context panel is open
    When user double-clicks on BDDEdgeDst{time} tree node inside browse tree
    Then the "BDDEdgeDst{time}" view should be current
    And "BDDEdgeMoved{time}" link in gallery should be visible
    When user clicks on "BDDEdgeMoved{time}" link in gallery
    Then the context panel should show "BDDEdgeMoved{time}"
    And no errors should have been logged
    And no error or warning balloon should have been shown
    Given Files tree node inside browse tree is expanded

  Scenario: The demo files kept their originals, and a package script is linked with no dialog
    When user clicks on "Files > Demo" tree node inside browse tree
    Then the "Demo" view should be current
    And demog.csv link in gallery should be visible
    And TSLA.csv link in gallery should be visible
    And beer.csv link in gallery should be visible
    Given the "BiostructureViewer" package is installed
    And Platform tree node inside browse tree is expanded
    And Platform---Functions tree node inside browse tree is expanded
    When user clicks on Platform---Functions---Scripts tree node inside browse tree
    And user enters "Protein-Ligand" into gallery search
    Then "Protein-Ligand Interaction Diagram" gallery card should be visible
    When user drags "Protein-Ligand Interaction Diagram" gallery card to BDDEdgeDst{time} tree node inside browse tree
    Then Move entity dialog should be absent
    When user clears gallery search
    And user double-clicks on BDDEdgeDst{time} tree node inside browse tree
    Then the "BDDEdgeDst{time}" view should be current
    And "Protein-Ligand Interaction Diagram" link in gallery should be visible
    And no errors should have been logged
    And no error or warning balloon should have been shown

  @known-failure:192 @realizes:core#1042
  Scenario: A file linked from one space into another shows in the target space
    core#1042: the Link fails with "FormatException: SyntaxError: Unexpected end of JSON input" and the
    "Move entity" dialog stays open. Remove the tag when core#1042 is fixed.
    When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree
    Then the "BDDEdgeSrc{time}" view should be current
    And "edgecars{time}.csv" link in gallery should be visible
    When user drags "edgecars{time}.csv" link in gallery to BDDEdgeDst{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
