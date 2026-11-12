@journey @serial
Feature: Chem's settings in the Plugins view
  Platform > Plugins > Chem > Settings lists Chem's package settings under a bold header per category, in the order
  package.json lists them: Chemical engine (Sketcher, Chemical Engine), Rendering (Highlight color, Mark unspecified
  stereo, Simplified stereo label), then Enumeration (Default Reactions, Default BBs, Default Reagents).
  With the context panel narrowed, the form stacks its labels over the inputs and SAVE stays in sight (the buttons
  row's spacer label used to take the whole row and push SAVE out of the panel).
  A Rendering setting saved there reaches the open page at once: a molecule column's Chemistry | Rendering pane,
  opened after the save and without a reload, shows the new Highlight color as the column's.

  The settings are the stand's, for all its users, so the feature runs alone (@serial) and puts them back at its end.
  The settings' storage per group and `_package.settings` after a save are ApiTests' "Packages: Properties"; what the
  Rendering settings draw, in both engines, is Chem's package tests ('chemical engine', 'rendering').

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the "Chem" package settings are put back at feature end
    And the context panel is open
    And Platform tree node inside browse tree is expanded

  Scenario: The settings are grouped under their categories, in package.json's order
    When user clicks on Platform---Plugins tree node inside browse tree
    And user clicks on "Chem" gallery card
    Then the context panel should show "Chem"
    When user expands Settings accordion header in context panel
    Then first "Chemical engine" label in "Settings" pane in context panel should be above "Sketcher" input in "Settings" pane in context panel
    And "Chemical Engine" input in "Settings" pane in context panel should be above "Rendering" label in "Settings" pane in context panel
    And "Rendering" label in "Settings" pane in context panel should be above "Highlight color" input in "Settings" pane in context panel
    And "Simplified stereo label" checkbox in "Settings" pane in context panel should be above "Enumeration" label in "Settings" pane in context panel
    And "Enumeration" label in "Settings" pane in context panel should be above "Default Reactions" input in "Settings" pane in context panel
    And "Mark unspecified stereo" checkbox in "Settings" pane in context panel should be checked
    And "Simplified stereo label" checkbox in "Settings" pane in context panel should be checked
    And no errors should have been logged

  Scenario: SAVE stays in sight in a narrow context panel
    When user drags context panel splitter by 250 and 0 pixels
    Then context panel should be at most 300 px wide
    When user scrolls SAVE button in "Settings" pane in context panel into view
    Then SAVE button in "Settings" pane in context panel should be inside context panel
    And SAVE button in "Settings" pane in context panel should not be clipped
    When user drags context panel splitter by -250 and 0 pixels
    Then no errors should have been logged

  Scenario: A saved Highlight color is every molecule column's at once
    When user enters "#2080ff" into "Highlight color" input in "Settings" pane in context panel
    And user clicks on SAVE button in "Settings" pane in context panel
    Then SAVE button in "Settings" pane in context panel should be disabled
    When user opens chembl-scaffolds dataset
    And user clicks on the "header Smiles" area of grid
    And user expands Chemistry accordion header in context panel
    And user expands Rendering accordion header in context panel
    Then "Highlight color" input in "Rendering" pane in context panel should have value "#2080ff"
    And no errors should have been logged
