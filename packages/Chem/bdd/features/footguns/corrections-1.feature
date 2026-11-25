@footguns @batch-6
Feature: Corrections: RDKit mode on this page only, and Chem's own urea formation once fixed
  RDKit mode: `Given the Chemical Engine is "RDKit"` is Chem's existing gate step: on a Crux stand it skips the
  scenario, so RDKit-mode scenarios written with it prove nothing. These use `Given the Chemical Engine of this
  page is "RDKit"`, bound page-local (this page reads Chem's engine setting as RDKit; nothing is written to the
  stand; the next scenario's page starts fresh). Each scenario below replaces the one of the same name in the
  other footguns features (named in its comment); the other scenarios of those features stand as they are.
  Urea formation: Chem's own default "Urea Formation" (dr014), whose template wrote CO2 as a third reactant, is
  fixed, so it has nothing to warn about. So the three-reactant warn uses a reaction the user saves for it, and
  Urea Formation gets a scenario of its own that it works (the fixed template gives CCCNC(=O)NCC on ethylamine +
  propylamine in both engines; the saved three-reactant one fails every pair in both). These replace
  templates.feature's "The Two-Component dialog warns before Run that Chem's urea formation takes three
  reactants".

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the Chemical Engine is "Crux"
    And the molecule sketcher is "Crux"

  # replaces flags-footer.feature's "RDKit mode shows no flag in the editor and lets OK write as before"
  @catalogue-UI-03
  Scenario: RDKit mode shows no flag in the editor and lets OK write as before
    Given the Chemical Engine of this page is "RDKit"
    And user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "C(C)(C)(C)(C)C" into the Crux sketcher of the dialog
    Then the footer of the sketcher dialog shows no flag
    And the OK button of the sketcher dialog is enabled
    When user clicks on OK button in sketcher dialog
    Then the value in row 1 of "molecule" column is no longer "CCO"

  # replaces templates.feature's "RDKit mode: a template Crux mode flags has no Checks section and the reaction editor
  # no flag"
  @catalogue-UI-03
  Scenario: RDKit mode: a template Crux mode flags has no Checks section and the reaction editor no flag
    Given the Chemical Engine of this page is "RDKit"
    And user opens a table "templates" with a reaction column "template" of units "smarts":
      | template                                            |
      | [C:1](=[O:2])[OH].[NH2:1]>>[C:1](=[O:2])[NH:3]      |
    When user clicks on the "cell 1 of template" area of grid
    Then the context panel has no Checks section
    # Chem's top menu shows only for a table with a molecule column, and this one has none
    Given user opens a table "molecules" with:
      | molecule |
      | CCO      |
    And the semantic types of the current table are detected
    When user picks "Chem | Transform | Reactions | Transformation..." from the top menu
    And user clicks the "New Reaction" button of the dialog
    And user enters "[#6:1][C:2](=O)O>>[#6:1][C:2](=O)[*${n}]" as the reaction SMARTS
    Then the footer of the "New Reaction" dialog shows no flag

  # replaces sketcher-warn.feature's "RDKit mode: the warns above show nowhere"
  @catalogue-UI-03 @catalogue-MOL-04 @catalogue-MOL-05
  Scenario: RDKit mode: the warns above show nowhere
    Given the Chemical Engine of this page is "RDKit"
    And user opens a table "molecules" with a SMILES column "molecule" holding "CO"
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user puts "[CH3:1][OH:2]" into the Crux sketcher of the dialog
    Then the footer of the sketcher dialog shows no flag
    When user clicks on OK button in sketcher dialog
    And user clicks on the "cell 1 of molecule" area of grid
    Then the context panel has no Checks section

  # replaces reactions-panel.feature's "RDKit mode has no Reaction panel and no Checks section for these cells"
  @catalogue-UI-03 @catalogue-MAP-01
  Scenario: RDKit mode has no Reaction panel and no Checks section for these cells
    Given the Chemical Engine of this page is "RDKit"
    And user opens a table "reactions" with a reaction column "reaction" of units "smiles":
      | reaction              |
      | C(C)(C)(C)(C)C>>CO    |
    When user clicks on the "cell 1 of reaction" area of grid
    Then the context panel has no Checks section
    And the context panel has no Reaction panel

  # replaced reaction-editor.feature's "RDKit mode: the reaction cell editor (where it opens) shows no flag". A
  # double-click on any reaction cell opens the reaction editor in both modes, and it is the Crux editor whatever the
  # engine, its footer's checks too, so erasing the arrow blocks OK in RDKit mode as well
  @catalogue-UI-03 @catalogue-RXN-16
  Scenario: RDKit mode: a double-click opens the reaction editor, and erasing the arrow blocks OK
    Given the Chemical Engine of this page is "RDKit"
    And user opens a table "reactions" with a reaction column "reaction" of units "smiles":
      | reaction                                  |
      | CC(=O)O.NCc1ccccc1>>CC(=O)NCc1ccccc1      |
    When user double-clicks on the "cell 1 of reaction" area of grid
    Then reaction editor dialog should be visible
    When user erases the arrow in the reaction editor
    Then the middle slot of the footer shows a block flag mentioning "arrow"
    And the OK button of the reaction editor is disabled
    When user clicks on CANCEL button of the reaction editor
    Then the value in row 1 of "reaction" column is exactly what it was

  # replaces templates.feature's "The Two-Component dialog warns before Run that Chem's urea formation takes three
  # reactants": Chem's own Urea Formation is fixed, so the three-reactant warn uses a reaction saved for it
  @catalogue-TPL-10 @catalogue-TRF-01
  Scenario: The Two-Component dialog warns before Run that a saved reaction takes three reactants
    Given user opens a table "pairs" with:
      | acid    | amine |
      | CC(=O)O | CCN   |
    And the semantic types of the current table are detected
    # Chem's detector takes no column of one short SMILES under a name like these, so the types are set
    And user sets the semantic type of "acid" column to "Molecule"
    And user sets the semantic type of "amine" column to "Molecule"
    And the user's saved reactions are put back when the feature ends
    And the user has a saved two-component reaction "Three reactants" "[C:1](=[O:2])[OH].[NH2:3].CN(C)C=O>>[C:1](=[O:2])[NH:3]"
    When user picks "Chem | Transform | Reactions | Two-Component Reaction..." from the top menu
    And user selects the reaction "Three reactants" in the dialog
    Then the middle slot of the footer shows a warn flag mentioning "3 reactants" and "2"
    And the OK button of the dialog is enabled
