@guide @sketcher-controls
Feature: R-groups and attachment points in Crux
  A guide: how do I mark R-groups and attachment points in Crux, Chem's molecule sketcher? The R-group tool, clicked on
  an atom, opens the R-Group dialog, whose R1 makes the atom R1; the attachment point tool, in the same palette,
  marks an atom where the fragment attaches. Crux writes both in the CXSMILES the sketcher holds, read as written. The
  feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Scenario: Turn the phenol's OH into R1 and mark the acid's carbon as the attachment point
    Given user is logged in
    And simple mode is off
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And the Crux sketcher is open on "Oc1ccc(cc1)C(=O)O"
    # caption: Pick the R-group tool
    When user clicks on Crux R-group tool
    # caption: Click the OH oxygen
    And user clicks on the "atom 0" area of Crux sketcher widget
    # caption: Choose R1
    And user clicks on Crux R1 button
    # caption: OK
    And user clicks on Crux R-Group OK button
    # caption: The oxygen is now R1
    Then the "smiles" reading of Crux sketcher widget should be the molecule "O=C(O)c1ccc([*:1])cc1"
    # caption: Open the R-group palette
    When user clicks on Crux R-group tool
    # caption: Pick the attachment point tool
    And user clicks on Crux attachment point tool
    # caption: Click the acid's carbon
    And user clicks on the "atom 7" area of Crux sketcher widget
    # caption: Mark it as the primary attachment point
    And user checks Crux primary attachment point checkbox
    # caption: OK
    And user clicks on Crux Attachment Points OK button
    # caption: The acid's carbon is marked as where the fragment attaches
    Then the "smiles" reading of Crux sketcher widget should be "O=C(O)c1ccc([*:1])cc1 |atomProp:1.molAttchpt.1:7.dummyLabel.R1|"
