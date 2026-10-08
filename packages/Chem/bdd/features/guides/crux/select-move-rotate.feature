@guide @sketcher-controls
Feature: Select, move, turn, flip and clean up in Crux
  A guide: how do I move a structure in Crux, Chem's molecule sketcher? With the selection tool, a double-click on a
  structure selects all of it; a drag moves it, the round handle above it turns it, the flip button beside it
  mirrors it, and Clean Up tidies its bonds and angles. Only the drawing changes: the molecule the sketcher holds stays
  naphthalene and ethanol (read by RDKit), while the atoms' coordinates in its molblock change. The feature drives
  Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Select the naphthalene, move it, turn it, flip it and clean up
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "c1ccc2ccccc2c1.OCC"
    And user remembers the coordinates of the Crux sketcher's molblock
    # caption: Pick the selection tool
    When user clicks on Crux select tool
    # caption: Double-click the naphthalene: all of it is selected
    And user double-clicks on the "atom 0" area of Crux sketcher widget
    # caption: Drag it to the left
    And user drags the "atom 0" area of Crux sketcher widget by 120 pixels to the left
    # caption: Drag the round handle to turn it
    And user drags Crux rotate handle by 90 and 0 pixels
    # caption: Flip it
    And user clicks on Crux flip horizontal button
    # caption: Clean Up tidies the drawing
    And user clicks on Crux clean up button
    # caption: Only the drawing moved: still naphthalene and ethanol
    Then the sketcher in sketcher dialog should hold the molecule "c1ccc2ccccc2c1.CCO"
    And the Crux sketcher's molblock should have new coordinates
