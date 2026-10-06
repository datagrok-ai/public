@demo @realizes:u2.dg.molecules
Feature: Bridged chemistry inputs
  Needs the Chem package on the stand: the sketcher input, the semType-picked editor and the
  structure typeahead are all bridged to it.

  Scenario: The structure input, the property form and the structure typeahead
    Given user opens the "Molecules" demo page
    Then value of smiles readout should have text "CC(=O)OC1=CC=CC=C1C(=O)O"
    And Structure input should be visible
    When user types "Aspirin acetate" into Name input in object form
    Then value of compound readout should contain text "Aspirin acetate"
    When user enters "181" into "MW, Da" input in object form
    Then value of compound readout should contain text "\"mw\":181"
    When user types "Caf" into compound picker
    And user selects "Caffeine" in compound picker
    Then value of picked readout should have text "Caffeine"
