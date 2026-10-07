@tutorials @serial @realizes:tutorials.sticky-meta
Feature: The Sticky Meta tutorial
  Walks Data access > Sticky Meta from its card to the end: an entity type made under Browse >
  Platform > Sticky Meta > Types, a schema for it under Schemas with one string property, the first
  molecule of the tutorial's table annotated through the Sticky meta pane, and the cell hovered. Each step is claimed as
  ticked and as done — the dialogs, the fields, the type and the schema on the server, the pane and
  its saved state.
  Not claimed: the annotation read back. The tutorial's last step asks for the cell's tooltip, which
  shows the molecule; the annotation has a tooltip of its own on the dot drawn in the cell's corner,
  and the grid reports no area for that dot.
  Translated from playwright-tests/e2e/stickymeta (the TestTrack Sticky Meta cases walked the same
  screens outside the tutorial).

  Everything the tutorial saves has a fixed name, so it is swept before the walk and removed after it
  (the annotated value goes with its schema). Serial, for the fixed names and because a finished
  tutorial writes its completion record into the account's settings, which every page syncs whole.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the package autostarts have completed
    And the "tutorials" user settings are put back at feature end
    And the "achievement-badges" user settings are put back at feature end
    And the Sticky Meta schema "schema for tutorial" and entity type "molecule-tutorial" are removed now and at feature end
    And the "Sticky Meta" tutorial is not completed yet
    And the Tutorials app is open

  Scenario: A learner completes the Sticky Meta tutorial
    When user starts the "Sticky Meta" tutorial
    Then the tutorial progress should be 1 of 20
    Given the tutorial step "Open Types node" should not be done yet
    When user expands Platform tree node inside browse tree
    And user expands Platform---Sticky-Meta tree node inside browse tree
    And user clicks on Platform---Sticky-Meta---Types tree node inside browse tree
    Then the tutorial step "Open Types node" should be done
    Given the tutorial step "Create a new entity type" should not be done yet
    When user clicks on "New Entity Type..." button
    Then the tutorial step "Create a new entity type" should be done
    Given the tutorial step "Explore entity type dialog" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore entity type dialog" should be done
    When user enters "molecule-tutorial" into "Name" input in "Create a new entity type" dialog
    Then the tutorial step "Set \"Name\" to \"molecule-tutorial\"" should be done
    When user enters "semtype=Molecule" into "Matching expression" input in "Create a new entity type" dialog
    Then the tutorial step "Set \"Matching expression\" to \"semtype=Molecule\"" should be done
    When user clicks on OK button in "Create a new entity type" dialog
    Then the tutorial step "Save entity type" should be done
    Then the entity type "molecule-tutorial" should exist

    Given the tutorial step "Open schemas node" should not be done yet
    When user clicks on Platform---Sticky-Meta---Schemas tree node inside browse tree
    Then the tutorial step "Open schemas node" should be done
    Given the tutorial step "Create a new schema" should not be done yet
    When user clicks on "New Schema..." button
    Then the tutorial step "Create a new schema" should be done
    Given the tutorial step "Explore schema dialog" should not be done yet
    When user goes through the tour to its end
    Then the tutorial step "Explore schema dialog" should be done
    When user enters "schema for tutorial" into "Name" input in "Create a new schema" dialog
    Then the tutorial step "Set \"Name\" to \"schema for tutorial\"" should be done
    When user clicks on "select entities" action in "Create a new schema" dialog
    Then the tutorial step "Select associated entity" should be done
    And "Select types for schema for tutorial" dialog should be visible
    When user checks "molecule-tutorial" property in "Select types for schema for tutorial" dialog
    Then the tutorial step "Select molecule-tutorial" should be done
    When user clicks on OK button in "Select types for schema for tutorial" dialog
    Then the tutorial step "Confirm entity selection" should be done
    When user enters "project name" into second "Name" input in "Create a new schema" dialog
    Then the tutorial step "Set property \"Name\" to \"project name\"" should be done
    When user selects "string" in "Property Type" input in "Create a new schema" dialog
    Then the tutorial step "Set property \"Type\" to \"string\"" should be done
    When user clicks on OK button in "Create a new schema" dialog
    Then the tutorial step "Save schema" should be done
    And the Sticky Meta schema "schema for tutorial" should exist

    Given the tutorial step "In the Sticky Meta molecules table, click the first cell in the smiles column" should not be done yet
    When user clicks on the "cell 1 of smiles" area of grid
    Then the tutorial step "In the Sticky Meta molecules table, click the first cell in the smiles column" should be done
    And "Sticky meta" pane in context panel should be visible
    And the tutorial step "Set \"project name\" to \"Tutorial\"" should not be done yet
    When user types "Tutorial" into "project name" input in context panel
    Then the tutorial step "Set \"project name\" to \"Tutorial\"" should be done
    When user clicks on Save button in "Sticky meta" pane in context panel
    Then the tutorial step "Click SAVE under \"schema for tutorial\"" should be done
    And Save button in "Sticky meta" pane in context panel should be disabled
    When user hovers over the "cell 1 of smiles" area of grid
    Then the tutorial step "Hover a cell to verify metadata tooltip." should be done

    And the "Sticky Meta" tutorial should be completed
    And the tutorial should have listed 20 steps
    And the tutorial progress should be 20 of 20
    And no hint should be shown
    And no errors should have been logged
