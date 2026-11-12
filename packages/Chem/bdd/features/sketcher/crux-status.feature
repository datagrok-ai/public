@sketcher-controls
Feature: Crux in the platform's automation
  Crux Sketch reports itself to Datagrok's automation (crux-sketch's spike datagrok-platform; HOST-043 to 045): its root
  is named cruxSketch (data-u2-name), which crux sketcher widget is found by; its status names its parts as areas,
  canvas, actions (the top toolbar), tools (the left one), elements (the element palette), templates (the ring bar) and
  label editor while it is open, which are also parts of the widget here; each control a toolbar shows is a "tool
  <name>" area, named as Crux names it; and its readings say what it holds: the selected atoms, the tool, whether it
  holds a query, whether it is empty, and the change events it has said. Every gesture here settles on the sketcher's
  isRenderPending and onRendered, never on a timeout. The sketcher is the cell editor of a molecule cell.
  The feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And user opens a table "molecules" with:
      | molecule |
      | CCO      |
      | c1ccccc1 |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule "CCO"

  # HOST-044: the parts, as parts of the widget and as areas of its status
  Scenario: Crux's canvas, toolbars and label editor are parts of its widget and areas of its status
    Then canvas of crux sketcher widget should be visible
    And actions of crux sketcher widget should be visible
    And tools of crux sketcher widget should be visible
    And elements of crux sketcher widget should be visible
    And templates of crux sketcher widget should be visible
    And the "canvas" area of crux sketcher widget should lie below the "actions" area
    And the "canvas" area of crux sketcher widget should lie to the right of the "tools" area
    And the "elements" area of crux sketcher widget should lie to the right of the "canvas" area
    And the "templates" area of crux sketcher widget should lie below the "canvas" area
    And label editor of crux sketcher widget should be absent
    When user double-clicks on the "atom 1" area of crux sketcher widget
    Then label editor of crux sketcher widget should be visible
    And the "label editor" area of crux sketcher widget should be at least 10 pixels tall

  # HOST-045: a tool's area chooses it, and the canvas's areas take its edit; the readings follow
  Scenario: A click on a tool's area chooses that tool, and its edit on an atom's area is read back
    Then the "tool" reading of crux sketcher widget should be "bond.single"
    And the "empty" reading of crux sketcher widget should be "false"
    And the "changes" reading of crux sketcher widget should be 1
    When user clicks on the "tool element.n" area of crux sketcher widget
    Then the "tool" reading of crux sketcher widget should be "element.n"
    When user clicks on the "atom 2" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "CCN"
    And the "changes" reading of crux sketcher widget should be 2
    And the "query" reading of crux sketcher widget should be "false"

  # HOST-045: what is selected, in molfile order, and the canvas cleared
  Scenario: The selected atoms are read in molfile order, and a cleared canvas reads as empty
    Then the "selected atoms" reading of crux sketcher widget should be ""
    When user presses Control+A in crux canvas
    Then the "selected atoms" reading of crux sketcher widget should be "0, 1, 2"
    And the "selected bonds" reading of crux sketcher widget should be "0, 1"
    And the "tool" reading of crux sketcher widget should be "select.rect"
    When user clicks on the "tool clear" area of crux sketcher widget
    Then the "empty" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0
    And the "selected atoms" reading of crux sketcher widget should be ""
