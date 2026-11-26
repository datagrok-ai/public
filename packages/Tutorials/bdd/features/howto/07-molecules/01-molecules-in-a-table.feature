@guide @tutorial
Feature: How-to · Molecules
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). Molecules in a table, from opening to editing: smiles.csv opened
  from Files > Demo > chem has its SMILES text recognized as molecules (semantic type Molecule) and drawn as
  structures; a click on a structure makes it current,
  and the context panel's Chemistry > Properties pane (OpenChemLib, in the browser) lists its weight, LogP, PSA and
  formula, and follows the next molecule clicked. A double-click on a structure opens the sketcher (OpenChemLib,
  pinned for the feature and put back after it) on that molecule; a SMILES typed into its top field replaces the
  drawing, and OK writes it into the cell: row 1 becomes aspirin, and the panel then shows aspirin's formula C9H8O4.
  The Descriptors, Gasteiger, 3D Structure and Identifiers panes are never expanded (a container, server-side Python
  and outside services). The left pane stays as opening the file leaves it. The edit lives in the page's table only;
  nothing is made on the server. Merged from the former clips "Molecules in a table" and "Draw or edit a molecule".

  Scenario: View and edit molecules in a table
    Open a file of SMILES, read a molecule's properties, and draw or paste a new structure into a cell.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---Demo tree node inside browse tree is expanded
    # caption: Files > Demo > chem holds tables of molecules
    When user expands Files---Demo---chem tree node inside browse tree
    And user double-clicks on Files---Demo---chem---smiles.csv tree node inside browse tree
    Then the "smiles" view should be current
    And the table should have 1000 rows
    And "canonical_smiles" column should have semantic type "Molecule"
    When the video says "SMILES text is recognized and drawn as structures"
    # caption: Click a molecule to make it current
    And user clicks on the "cell 2 of canonical_smiles" area of grid
    Then row 2 should be current
    # caption: The context panel shows what is known about it
    When user clicks on "Context Panel" status bar toggle
    Then the context panel should show the current cell
    # caption: Chemistry > Properties: weight, LogP, PSA, formula...
    When user expands Chemistry accordion header in context panel
    And user expands Properties accordion header in context panel
    Then context panel should contain text "C22H22O4"
    # caption: Click another molecule: the panel follows
    When user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And context panel should not contain text "C22H22O4"
    # caption: Double-click a structure to open the sketcher
    When user double-clicks on the "cell 1 of canonical_smiles" area of grid
    Then sketcher dialog should be visible
    When the video says "Draw with the tools on the left, or type a SMILES above" at the top
    And the video says "Type or paste a SMILES, then press Enter" at the top
    And user clicks on molecule input of sketcher dialog
    And user types "CC(=O)Oc1ccccc1C(=O)O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And the video says "OK puts the molecule in the cell" at the top
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the value of "canonical_smiles" column in row 1 should be "CC(=O)Oc1ccccc1C(=O)O"
    And context panel should contain text "C9H8O4"
    When the video says "Row 1 is now aspirin, and the panel shows its properties"
    Then the table should have 1000 rows
    And no errors should have been logged
