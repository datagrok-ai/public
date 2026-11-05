@guide @sketcher-controls
Feature: Charges, isotopes and radicals in Crux
  A guide: how do I set an atom's charge, isotope or radical in Crux, Chem's molecule sketcher? The charge tools raise
  or lower the charge of an atom they click; Atom Properties, from an atom's menu, sets its isotope and its radical.
  Each result is the molecule the sketcher holds, read by RDKit. The feature drives Crux's own controls
  (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Make propanoate with the charge tool, label its methyl 13C and make the middle carbon a radical
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "CCC(=O)O"
    # caption: Pick the charge minus tool
    When user clicks on Crux charge minus tool
    # caption: Click the OH oxygen: it becomes O⁻
    And user clicks on the "atom 4" area of Crux sketcher widget
    # caption: Propanoate
    Then the sketcher in sketcher dialog should hold the molecule "CCC(=O)[O-]"
    # caption: Right-click the methyl carbon
    When user opens the Crux context menu on the "atom 0" area
    # caption: Choose Atom Properties…
    And user clicks on Crux Atom Properties item
    # caption: Set the isotope to 13
    And user types "13" into Crux isotope field
    # caption: Apply
    And user clicks on Crux Atom Properties Apply button
    # caption: The methyl is now ¹³C
    Then the sketcher in sketcher dialog should hold the molecule "[13CH3]CC(=O)[O-]"
    # caption: Right-click the middle carbon
    When user opens the Crux context menu on the "atom 1" area
    # caption: Choose Atom Properties…
    And user clicks on Crux Atom Properties item
    # caption: Make it a monoradical: one unpaired electron
    And user selects "Monoradical" in Crux radical list
    # caption: Apply
    And user clicks on Crux Atom Properties Apply button
    # caption: The middle carbon now carries a radical
    Then the sketcher in sketcher dialog should hold the molecule "[13CH3][CH]C(=O)[O-]"
