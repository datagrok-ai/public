@guide @sketcher-controls
Feature: Enhanced stereochemistry in Crux
  A guide: how do I say that a drawing's stereocentres are relative, or that it stands for a mixture, in Crux, Chem's
  molecule sketcher? Selected stereocentres go into an enhanced stereo group from the menu: AND (&1, the drawing and
  its mirror image, a racemate) or OR (or1, one of the two, which is unknown). Crux writes the group in the CXSMILES
  the sketcher holds (`|&1:1,3|`, `|o1:1,3|`), read here as written. Crux's settings are put back when the feature
  ends. The feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Put two stereocentres in an AND group, then an OR group
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "C[C@H](O)[C@@H](C)C(=O)O", showing R, S, E and Z labels
    # caption: Select everything (Ctrl+A)
    When user presses Control+A in Crux canvas
    # caption: Right-click a stereocentre
    And user opens the Crux context menu on the "atom 1" area
    # caption: Choose Enhanced Stereochemistry…
    And user clicks on Crux enhanced stereo item
    # caption: Choose AND: the drawing and its mirror image, a racemic mixture
    And user clicks on Crux AND button
    # caption: Apply
    And user clicks on Crux enhanced stereo Apply button
    # caption: Both centres are marked &1
    Then the "smiles" reading of Crux sketcher widget should be "C[C@H](O)[C@@H](C)C(=O)O |&1:1,3|"
    # caption: Right-click a stereocentre again
    When user opens the Crux context menu on the "atom 1" area
    # caption: Choose Enhanced Stereochemistry…
    And user clicks on Crux enhanced stereo item
    # caption: Choose OR: one of the two, not known which
    And user clicks on Crux OR button
    # caption: Apply
    And user clicks on Crux enhanced stereo Apply button
    # caption: Both centres are marked or1
    Then the "smiles" reading of Crux sketcher widget should be "C[C@H](O)[C@@H](C)C(=O)O |o1:1,3|"
