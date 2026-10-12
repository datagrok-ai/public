@demo @crux @sketcher-controls
Feature: A u2 form's molecule field, edited in Crux
  u2's object form picks the platform's molecule input for a property of semantic type Molecule (fromDartInput over
  ui.input.molecule): the field is a tab stop, and Enter or Space on it opens the sketcher. The Molecules demo page binds
  such a form (Name, Structure, MW) to a compound and shows the compound as its "compound" readout, so what the form
  wrote back through the bound property is what the readout shows. Here the sketcher is Crux, Chem's own, pinned as the
  session's sketcher; the edit is drawn on Crux's own controls. Enter opens the sketcher only on the field that has the
  keyboard: the Tab from Name lands on Structure, the next field.

  Needs Chem on the stand. The feature drives Crux's own controls (@sketcher-controls): a run that pins another
  sketcher skips it.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-068
  Scenario: Tab to the form's Structure field and Enter opens Crux, and the drawn molecule writes back through the bound property
    Given user opens the "Molecules" demo page
    When user clicks on Name input in object form
    And user presses Tab
    And user presses Enter
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)OC1=CC=CC=C1C(=O)O"
    When user clicks on crux clear button
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And value of compound readout should contain text "\"smiles\":\"c1ccccc1\""
    And value of compound readout should contain text "\"name\":\"Aspirin\""
