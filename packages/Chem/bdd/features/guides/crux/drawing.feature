@guide @sketcher-controls
Feature: Draw a molecule in Crux with bonds, a chain and hotkeys
  A guide: how do I draw a molecule in Crux, Chem's molecule sketcher? A click with the bond tool puts a bond on the
  empty canvas, a click on an atom grows the chain from it, and the chain tool draws a zig-zag chain as far as it is
  dragged. Over an atom or a bond, Crux's hotkeys change it at once: a letter makes the atom that element (o oxygen,
  n nitrogen), a digit sets the bond order (2 double), and over an atom a digit sprouts a group (3 a phenyl ring).
  Each result is the molecule the sketcher holds, compared by RDKit. The feature drives Crux's own controls
  (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Draw bonds and a chain, then change atoms and bonds with keys
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open
    # caption: Pick the bond tool
    When user clicks on Crux single bond tool
    # caption: Click the empty canvas: a bond appears
    And user clicks on Crux canvas
    # caption: Click an end of the bond to add another
    And user clicks on the "atom 1" area of Crux sketcher widget
    # caption: Pick the chain tool
    And user clicks on Crux chain tool
    # caption: Drag from the end of the chain to draw more carbons
    And user drags the "atom 2" area of Crux sketcher widget by 175 pixels to the right
    # caption: An octane chain
    Then the sketcher in sketcher dialog should hold the molecule "CCCCCCCC"
    # caption: Point at the last carbon and press o: it becomes an oxygen
    When user presses the "o" key over the "atom 7" area of Crux sketcher widget
    # caption: Point at the first bond and press 2: it becomes a double bond
    And user presses the "2" key over the "bond 0" area of Crux sketcher widget
    # caption: Point at a carbon and press 3: a phenyl ring sprouts from it
    And user presses the "3" key over the "atom 3" area of Crux sketcher widget
    # caption: 4-Phenylhept-6-en-1-ol, drawn with three keys
    Then the sketcher in sketcher dialog should hold the molecule "C=CCC(CCCO)c1ccccc1"
