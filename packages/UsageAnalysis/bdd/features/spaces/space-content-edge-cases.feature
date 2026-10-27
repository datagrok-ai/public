@serial @spaces @edge-cases @realizes:views.space
Feature: Space content edge cases: a file linked between spaces, a link whose original moved, a package script linked
  Guards for putting entities into a space by dragging them onto the space's node in the Browse tree.
  A file dropped on a space offers Link, Copy and Move in the "Move entity" dialog
  (`ProjectEntityMoveWidget.getActions`). A file linked from one space into another is expected to
  show in the target space and the dialog to close (core#1042). Today the Link fails with
  "FormatException: SyntaxError: Unexpected end of JSON input" and the dialog stays open; that
  scenario is a known failure, and the tag goes when core#1042 is fixed.

  A link follows its original: a project linked into one space and then moved from its own space into
  a third one still shows its details when its link is clicked, with no error and no refresh of the
  tree (core#382; the ticket's entity is a file, whose Link is core#1042, so the guard is claimed on a
  project, which Move to Space... puts into the source space first).

  A function of a package (here a script of BiostructureViewer) can only be linked: dropped on a space
  it is linked at once, with no "Move entity" dialog and no choice of action (`ProjectMeta.doDropEx`).

  The file is the project's fixtures/cars-small.csv, written into the source space through the API.
  The spaces and the project are this feature's own, named with the run's time, made or saved through
  the API and swept with what they hold when the feature starts and ends; the package script stays
  where it is, only its link goes with the space. Every drag of a file starts from the source space's
  own view, where the dialog is shown before its Link/Copy/Move chooser, so the chooser is waited for.

  Serial with the other features that drag onto the Spaces tree.

  Background:
    Given user is logged in
    And the browse panel is open
    And a space named "BDDEdgeSrc{time}" is on the server
    And a space named "BDDEdgeDst{time}" is on the server
    And the "fixtures/cars-small.csv" file of the project is in the space "BDDEdgeSrc{time}" as "edgecars{time}.csv"
    And Spaces tree node inside browse tree is expanded
    When user double-clicks on BDDEdgeSrc{time} tree node inside browse tree
    Then the "BDDEdgeSrc{time}" view should be current
    And "edgecars{time}.csv" link in gallery should be visible

  @known-failure:45 @realizes:core#1042
  Scenario: A file linked from one space into another shows in the target space
    core#1042: the Link fails with "FormatException: SyntaxError: Unexpected end of JSON input" and the
    "Move entity" dialog stays open. Remove the tag when core#1042 is fixed.
    When user drags "edgecars{time}.csv" link in gallery to BDDEdgeDst{time} tree node inside browse tree
    Then Move entity dialog should be visible
    And choice input in Move entity dialog should have value "Link"
    When user clicks on YES button in Move entity dialog
    Then Move entity dialog should be hidden
    And no error or warning balloon should have been shown
    When user double-clicks on BDDEdgeDst{time} tree node inside browse tree
    Then the "BDDEdgeDst{time}" view should be current
    And "edgecars{time}.csv" link in gallery should be visible

  Scenario: A package script dropped on a space is linked with no dialog
    Given Platform tree node inside browse tree is expanded
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

  Scenario: A linked project still opens its details after its original moved to another space
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
