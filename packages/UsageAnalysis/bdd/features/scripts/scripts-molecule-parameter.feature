@serial @crux @sketcher-controls
Feature: A script's molecule parameter, sketched in Crux
  A function parameter of semantic type Molecule is the platform's molecule input in the function's dialog: a click on
  it opens a sketcher dialog, and its OK writes the sketcher's SMILES into the parameter. Here the sketcher is Crux,
  Chem's own, pinned as the session's sketcher, and the molecule is drawn on Crux's own controls, one bond at a time:
  what the function receives is what Crux drew, as SMILES.

  The script is this feature's own ({time} in its name), a JavaScript script that hands its parameter back as its
  output, run from its editor (Run script, whose results show the output), and is deleted with its chats at the
  end. Needs Chem on the stand. The
  feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Serial: it works in the Scripts view, whose search text is the account's own setting (as scripts-run.feature).

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And a script "BddCruxMolecule{time}" is on the server:
      """
      //language: javascript
      //input: string mol {semType: Molecule}
      //output: string got
      got = mol;
      """

  # HOST-065
  Scenario: A pentane drawn in Crux for the Molecule parameter reaches the script as CCCCC
    Given user opens the Scripts view
    When user clears gallery search
    And user types "BddCruxMolecule{time}" into gallery search
    And user double-clicks on "BddCruxMolecule{time}" link in gallery
    Then the "BddCruxMolecule{time}" view should be current
    When user clicks on "Run script (F5)" icon
    Then "BddCruxMolecule{time}" dialog should be visible
    When user clicks on editor of mol input in "BddCruxMolecule{time}" dialog
    Then sketcher dialog should be visible
    And the "ready" reading of crux sketcher widget should be "true"
    When user clicks on crux single bond tool
    And user clicks on crux canvas
    And user clicks on the "atom 1" area of crux sketcher widget
    And user clicks on the "atom 2" area of crux sketcher widget
    And user clicks on the "atom 3" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "CCCCC"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    When user clicks on OK button in "BddCruxMolecule{time}" dialog
    Then the "BddCruxMolecule{time}" dialog should close
    And the script results should show "got" as "CCCCC"
