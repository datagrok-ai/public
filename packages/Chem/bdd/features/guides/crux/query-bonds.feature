@guide @sketcher-controls
Feature: Query bonds in Crux
  A guide: how do I draw bonds that match more than one bond type, or only in a chain, in Crux, Chem's molecule
  sketcher? In query mode the query bond tool (Any bond at first) makes a bond it clicks match any bond, and a bond's
  menu sets its topology: Ring or Chain. The query is the SMARTS the sketcher writes, read as written. The feature
  drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Make one bond any bond and another a chain bond
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "NCCc1ccccc1"
    # caption: Open the gear menu
    When user clicks on Crux settings button
    # caption: Switch on Query mode
    And user clicks on Crux query mode item
    # caption: Pick the query bond tool: Any bond
    And user clicks on Crux query bond tool
    # caption: Click the bond to the ring: it may now be any bond
    And user clicks on the "bond 2" area of Crux sketcher widget
    # caption: The query takes any bond there (~)
    Then the Crux sketcher should hold the query "[#7]-[#6]-[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    # caption: Right-click the C–C bond of the side chain
    When user opens the Crux context menu on the "bond 1" area
    # caption: Open Topology
    And user clicks on Crux topology item
    # caption: Choose Chain: the bond must not be in a ring
    And user clicks on Crux chain topology item
    # caption: The query asks for a chain bond there (!@)
    Then the Crux sketcher should hold the query "[#7]-[#6]-&!@[#6]~[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
