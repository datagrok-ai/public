@guide @sketcher-controls
Feature: Wedges and R/S labels in Crux
  A guide: how do I set a stereocentre in Crux, Chem's molecule sketcher, and see whether it is R or S? Settings (the
  gear) turns on the R, S, E and Z labels; the wedge tool, clicked on a bond, wedges it from its first atom, which
  sets the centre; a bond's own menu turns the wedge into a hashed wedge, which inverts it. The molecule the sketcher
  holds is read by RDKit after each step: alanine with no configuration, then D-alanine (R), then L-alanine (S).
  Crux's settings, which a user's change keeps in the browser, are put back when the feature ends. The feature drives
  Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Wedge a bond, read the centre's R or S, and invert it with a hashed wedge
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "CC(N)C(=O)O"
    # caption: Open the gear menu
    When user clicks on Crux settings button
    # caption: Choose Settings…
    And user clicks on Crux settings item
    # caption: Turn on "Show R, S, E and Z labels"
    And user checks Crux stereo labels checkbox
    # caption: Apply
    And user clicks on Crux settings Apply button
    # caption: Pick the wedge bond tool
    And user clicks on Crux stereo bond tool
    # caption: Click the bond from the centre to the nitrogen: it becomes a wedge
    And user clicks on the "bond 1" area of Crux sketcher widget
    # caption: The centre is labelled (R): this is D-alanine
    Then the "smiles" reading of Crux sketcher widget should be the molecule "C[C@@H](N)C(=O)O"
    # caption: Right-click the wedge
    When user opens the Crux context menu on the "bond 1" area
    # caption: Choose Hashed wedge bond
    And user clicks on Crux hash item
    # caption: The configuration flips from R to S: L-alanine
    Then the "smiles" reading of Crux sketcher widget should be the molecule "C[C@H](N)C(=O)O"
