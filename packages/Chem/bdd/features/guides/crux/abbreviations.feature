@guide @sketcher-controls
Feature: Abbreviations in Crux
  A guide: how do I use abbreviations such as CO2Me in Crux, Chem's molecule sketcher? A double-click on an atom opens
  its label editor; a known abbreviation typed there is drawn as its label, the group behind it. Its menu expands it
  into its atoms (and Contract Abbreviation folds it back), the molecule unchanged (the molecule the sketcher holds, read by RDKit, and
  its atoms drawn, read from the status). The feature drives Crux's own controls (@sketcher-controls): a run that pins
  another sketcher skips it.

  Scenario: Type CO2Me on toluene's methyl, then expand the ester into its atoms
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "Cc1ccccc1"
    # caption: Double-click the methyl carbon to edit its label
    When user double-clicks on the "atom 0" area of Crux sketcher widget
    # caption: Type CO2Me
    And user types "CO2Me" into label editor of Crux sketcher widget
    # caption: Press Enter: a methyl ester, drawn as CO2Me
    And user presses Enter in label editor of Crux sketcher widget
    # caption: Methyl benzoate
    Then the sketcher in sketcher dialog should hold the molecule "COC(=O)c1ccccc1"
    # caption: Right-click the CO2Me label
    When user opens the Crux context menu on the "atom 0" area
    # caption: Choose Expand Abbreviation: the ester's atoms are drawn
    And user clicks on Crux expand abbreviation item
    # caption: The same molecule, drawn atom by atom
    Then the sketcher in sketcher dialog should hold the molecule "COC(=O)c1ccccc1"
