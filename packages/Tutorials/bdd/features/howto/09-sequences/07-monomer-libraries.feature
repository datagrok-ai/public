@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). What a monomer is and where Datagrok keeps them: Bio > Manage >
  Monomers opens a library's monomers as a table (structure, symbol, name, the cap of each attachment point) beside
  the monomer editor, which follows the current row: the structure in the sketcher, its symbol and name, and a chip
  per attachment point (R1, R2, R3) whose menu lists the common caps; Bio > Manage > Monomer Libraries lists the
  library files Datagrok loads together. Browsing only: nothing is saved, switched or picked, so no library and no
  setting on the server changes. The molecule sketcher is pinned to Crux for the clip and restored after it; the
  sample's Activity is shown to two decimals.

  Scenario: Monomers and monomer libraries
    Every sequence is built from monomers; the libraries hold each one's structure and attachment points.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And user opens HELM_sample dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And user sets the format of "Activity" column in the grid to "0.00"
    When the video says "Monomers: the building blocks of every sequence"
    # caption: Bio > Manage > Monomers
    And user picks "Bio > Manage > Monomers" from the top menu
    Then the "Manage Monomers" view should be current
    And the monomer sketcher of the Manage Monomers view should be ready
    And "R1" column should have semantic type "Molecule"
    When the video says "A library as a table: structure, symbol, name, attachment points"
    # caption: Click one: its structure and details on the right
    And user clicks on the "cell 4 of Name" area of grid
    Then "Name" input in monomer form should have value "Glutamic acid"
    And "R3" R-group chip in monomer form should contain the text "R3 OH"
    # caption: R1, R2, R3: where it bonds; the cap when it does not
    When user clicks on "R3" R-group chip in monomer form
    Then "OH" menu item should be checked
    When user presses Escape
    # caption: Bio > Manage > Monomer Libraries: the library files in use
    And user picks "Bio > Manage > Monomer Libraries" from the top menu
    Then the "Manage Monomer Libraries" view should be current
    And "HELMCoreLibrary.json" checkbox should be checked
    When the video says "Library files load together; the ticked ones are in use"
    And the video says "A symbol two files define differently: you pick one on the right"
    Then no error or warning balloon should have been shown
    And no errors should have been logged
