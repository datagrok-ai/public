@guide @sketcher-controls
Feature: Copy and paste in Crux
  A guide: how do I copy a molecule out of Crux, Chem's molecule sketcher, in a format I choose, and paste one in?
  Copy As, beside Copy on the top toolbar, lists the formats (SMILES, CXSMILES, SMARTS, MOL V2000 and V3000, CDXML,
  InChI, InChIKey); SMILES puts aspirin's SMILES on the clipboard. Ctrl+V on the canvas pastes it back. The clipboard
  is read as text; the molecule the sketcher holds is read by RDKit. The feature drives Crux's own controls
  (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Copy aspirin as SMILES, clear the canvas and paste it back
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "CC(=O)Oc1ccccc1C(=O)O"
    # caption: Open Copy As
    When user clicks on Crux copy as button
    # caption: Choose SMILES
    And user clicks on Crux SMILES item
    # caption: The clipboard holds aspirin's SMILES
    Then the clipboard should contain text "CC(=O)Oc1ccccc1C(=O)O"
    # caption: Clear the canvas
    When user clicks on Crux clear button
    # caption: Press Ctrl+V on the canvas
    And user presses Control+V in Crux canvas
    # caption: Aspirin is back
    Then the sketcher in sketcher dialog should hold the molecule "CC(=O)Oc1ccccc1C(=O)O"
