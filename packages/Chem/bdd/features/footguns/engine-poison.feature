@footguns @engine-poison
Feature: A Crux sketcher whose engine stops, or that cannot write its drawing, never writes over the cell
  The engine-poison foot-gun: after an edit, a getter of Crux Sketch throws when its engine has stopped (a
  WebAssembly trap, or it never loaded) or cannot write the drawing in the column's notation. Chem's sketcher
  answered the platform '' and the molecule cell editor's OK wrote '' over the cell, with nothing said. Crux
  mode: a drawing that cannot be written blocks OK in the footer's middle slot. An engine stop (the user is
  never asked to reload: the engine recovers inside the call) is a quiet warning in that slot: OK stays on,
  writes nothing and keeps the cell, and no balloon is shown. No build stops the engine on demand: the steps
  make the dialog's Crux getters throw in the page (own properties over its handle's four getters; nothing
  else of it changes). The package tests "Foot-guns: Crux" pin the same with three throws (a raw RuntimeError
  too).

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And the molecule sketcher is "Crux"

  Scenario: The engine stopping under the molecule cell editor warns quietly, and neither Enter nor OK empties the cell
    Given user opens a table "molecules" with a SMILES column "molecule" holding "CCO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And the Crux engine of the sketcher dialog stops
    Then the middle slot of the footer shows a warn flag mentioning "could not be read back"
    And no error or warning balloon should have been shown
    When user presses Enter while the sketcher dialog has the focus
    Then the value in row 1 of "molecule" column is exactly what it was
    When user clicks the OK button of the sketcher dialog, disabled or not
    Then the value in row 1 of "molecule" column is exactly what it was
