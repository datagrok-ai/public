@guide @sketcher-controls
Feature: Explicit hydrogens in Crux
  A guide: how do I show every hydrogen of a molecule in Crux, Chem's molecule sketcher, and hide them again? The
  Hydrogens menu on the top toolbar draws each implicit hydrogen as an atom of its own (Add explicit hydrogens) and
  folds them back (Remove explicit hydrogens); the molecule stays L-alanine, its stereocentre kept. The atoms drawn are
  counted from the sketcher's status, the molecule read by RDKit. The feature drives Crux's own controls
  (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Draw L-alanine's hydrogens, then fold them back
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "C[C@H](N)C(=O)O"
    # caption: Open the Hydrogens menu
    When user clicks on Crux hydrogens button
    # caption: Choose Add explicit hydrogens
    And user clicks on Crux add hydrogens item
    # caption: All 13 atoms are drawn, the 7 hydrogens included
    Then the "atoms" reading of Crux sketcher widget should be 13
    And the "smiles" reading of Crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"
    # caption: Open the Hydrogens menu again
    When user clicks on Crux hydrogens button
    # caption: Choose Remove explicit hydrogens
    And user clicks on Crux remove hydrogens item
    # caption: The hydrogens fold back into their atoms: still L-alanine
    Then the "atoms" reading of Crux sketcher widget should be 6
    And the "smiles" reading of Crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"
