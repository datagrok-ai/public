@guide @sketcher-controls
Feature: E and Z double bonds in Crux
  A guide: how do I set a double bond's geometry in Crux, Chem's molecule sketcher? Crotonic acid drawn from a SMILES
  with no geometry has none; a right-click on its C=C bond offers E and Z, and the label beside the bond follows the
  choice (R, S, E and Z labels shown, as Settings shows them). The molecule the sketcher holds is read by RDKit:
  (E)-crotonic acid, then (Z)-, isocrotonic acid. Crux's settings are put back when the feature ends. The feature
  drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Make a double bond E, then Z, from its menu
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "CC=CC(=O)O", showing R, S, E and Z labels
    # caption: Right-click the C=C double bond
    When user opens the Crux context menu on the "bond 1" area
    # caption: Choose E
    And user clicks on Crux E item
    # caption: (E)-crotonic acid: the methyl and the acid on opposite sides
    Then the "smiles" reading of Crux sketcher widget should be the molecule "C/C=C/C(=O)O"
    # caption: Right-click the double bond again
    When user opens the Crux context menu on the "bond 1" area
    # caption: Choose Z
    And user clicks on Crux Z item
    # caption: The bond flips to Z: isocrotonic acid
    Then the "smiles" reading of Crux sketcher widget should be the molecule "C/C=C\C(=O)O"
