@guide @sketcher-controls
Feature: Query atoms in Crux
  A guide: how do I draw a substructure query with atom lists and generic atoms in Crux, Chem's molecule sketcher?
  Query mode (the gear's menu) turns the sketcher into a query editor; there the periodic table picks a list of
  elements (List), a list it must not be (Not list), or a generic atom (Q: any atom but carbon and hydrogen), which a
  click on an atom puts there. The query is the SMARTS the sketcher writes, read as written. The feature drives Crux's
  own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Switch to query mode, make a ring atom "N or O" and the acid's OH any heteroatom
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "OC(=O)c1ccccc1"
    # caption: Open the gear menu
    When user clicks on Crux settings button
    # caption: Switch on Query mode
    And user clicks on Crux query mode item
    # caption: Open the periodic table
    And user clicks on Crux periodic table button
    # caption: Choose List
    And user clicks on Crux periodic table list button
    # caption: Choose nitrogen
    And user clicks on "Nitrogen" button inside Crux periodic table
    # caption: … and oxygen
    And user clicks on "Oxygen" button inside Crux periodic table
    # caption: Add the list
    And user clicks on Crux periodic table Add button
    # caption: Click a ring atom: it may now be N or O
    And user clicks on the "atom 4" area of Crux sketcher widget
    # caption: The query asks for N or O at that ring position
    Then the Crux sketcher should hold the query "[#8]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1"
    # caption: Open the periodic table again
    When user clicks on Crux periodic table button
    # caption: Choose Q: any atom but carbon and hydrogen
    And user clicks on Crux periodic table Q button
    # caption: Add it
    And user clicks on Crux periodic table Add button
    # caption: Click the acid's OH oxygen: it becomes Q
    And user clicks on the "atom 0" area of Crux sketcher widget
    # caption: The query now takes any heteroatom there
    Then the Crux sketcher should hold the query "[!#6&!#1]-[#6](=[#8])-[#6]1:[#7,#8]:[#6]:[#6]:[#6]:[#6]:1"
