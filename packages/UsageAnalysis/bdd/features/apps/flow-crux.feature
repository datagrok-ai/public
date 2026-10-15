@crux @sketcher-controls
Feature: Flow's Sketcher Input node, drawn in Crux
  The Sketcher Input node holds a molecule for the flow. Its preview expands into an inline sketcher (at the canvas's
  zoom 1), and the node reports a molecule the user changed there as one edit of its parameters, the flow's
  "parameter edits" reading; a molecule put into the sketcher by the node itself (opening it again, the node's own sync
  after an edit) is not the user's and is not reported. Here the sketcher is Crux, Chem's own, pinned as the session's
  sketcher, and the user's changes are strokes on Crux's own controls.

  Needs Chem on the stand. The feature drives Crux's own controls (@sketcher-controls): a run that pins another
  sketcher skips it. Nothing is saved: the flow is never saved.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-069
  Scenario: Each stroke in the node's Crux is one edit of the flow, and opening the node's sketcher again is none
    Given user opens the "Flow" app
    When user clicks on flow blank canvas card
    And user clicks on "Inputs" pane header in toolbox
    And user double-clicks on flow sketcher input item
    Then the "parameter edits" reading of flow editor widget should be 0
    When user clicks on flow sketcher preview
    Then the "ready" reading of crux sketcher widget should be "true"
    And the "parameter edits" reading of flow editor widget should be 0
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    And the "parameter edits" reading of flow editor widget should be 1
    When user clicks on crux single bond tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccccc1"
    And the "parameter edits" reading of flow editor widget should be 2
    When user clicks on flow sketcher Done button
    And user clicks on flow sketcher preview
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccccc1"
    When user clicks on crux nitrogen tool
    And user clicks on the "atom 1" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccccn1"
    And the "parameter edits" reading of flow editor widget should be 3
