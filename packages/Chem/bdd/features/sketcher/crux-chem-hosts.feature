@sketcher-controls
Feature: Crux in Chem's hosts
  Crux Sketch, pinned as the session's sketcher, in the places Chem opens a sketcher beyond the cell editor: the options
  menu's Recent and Favorites, a molblock column's cell editor, the Sketch item of a molecule's menu, the substructure
  filter's card (its dialog, its thumbnail, a column popup's inline filter, a query kept through Clone View, a layout,
  a project and the panel reopened), Chem > Search > Substructure Search, the Scaffold Tree's sketcher, the Similarity
  Search viewer's reference, the column's Rendering and Highlight panes, the R-Groups Analysis dialog and the Deprotect
  editor. What Crux holds is drawn on Crux's own controls (its tools, its atoms as hit areas) and read through its
  status: its "smiles" reading compared with a molecule by RDKit, its "atoms" count. What Datagrok does with it is read
  where Datagrok shows it: the rows the filter passes, the card's readings ("molecule of <column>" is the very string
  the card filters by), the cells, the viewers' readings, a balloon. Typing into the host's molecule field is not
  drawing (the host then keeps the typed string): every claim about Crux's own output draws.
  The feature drives Crux's own controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the molecule sketcher is "Crux"
    And the package autostarts have completed

  # HOST-050: the cell editor's OK adds the drawing to Recent; Add to Favorites adds the one shown
  Scenario: A molecule picked from the options menu's Recent or Favorites is drawn in Crux
    Given the sketcher's Recent and Favorites are empty, and come back when the feature ends
    And user opens a table "molecules" with:
      | molecule |
      | CCO      |
      | c1ccccc1 |
      | CC(=O)O  |
    And the semantic types of the current table are detected
    When user double-clicks on the "cell 1 of molecule" area of grid
    And user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then the molecule in row 1 of "molecule" column should be "C1CCCCC1"
    When user double-clicks on the "cell 2 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    When user clicks on "Options" icon in sketcher dialog
    And user picks "Favorites > Add to Favorites" from the open menu
    And user clicks on "Options" icon in sketcher dialog
    And user picks molecule 1 of the "Recent" group from the open menu
    Then the "smiles" reading of crux sketcher widget should be the molecule "C1CCCCC1"
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 2 of "molecule" column should be "c1ccccc1"
    When user double-clicks on the "cell 3 of molecule" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)O"
    When user clicks on "Options" icon in sketcher dialog
    And user picks molecule 1 of the "Favorites" group from the open menu
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"

  # HOST-051: a molblock column (a SMILES column is sketcher/cell-editor.feature's, with Crux pinned)
  Scenario: In a molblock column's cell editor, a typed C1CCCCC1 is written as a molblock and the table keeps its rows
    Given user opens spgi-100 dataset
    When user double-clicks on the "cell 1 of Structure" area of grid
    Then the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "Structure" column
    When user types "C1CCCCC1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the molecule in row 1 of "Structure" column should be "C1CCCCC1"
    And every value of "Structure" column should contain "M  END"
    And the table should have 100 rows

  # HOST-053: Sketch opens a sketcher in a dialog of its own, with no OK and no Cancel: only its close icon. A table of
  # its own name: the 2D Structure pane is kept from an earlier table of the same name, column and row
  Scenario: Sketch in the menu of a molecule's drawing shows it in Crux, and an edit there changes nothing when it closes
    Given user opens a table "sketch-menu" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | CCO                |
    And the semantic types of the current table are detected
    And the context panel is open
    And user remembers the value of "molecule" column in row 1
    When user clicks on the "cell 1 of molecule" area of grid
    Then the context panel should show the current cell
    When user expands Structure accordion header in context panel
    And user expands "2D Structure" accordion header in context panel
    And user clicks on "More" icon in "2D Structure" pane in context panel
    And user picks "Sketch" from the open menu
    Then the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "molecule" column
    When user clicks on crux clear button
    Then the "atoms" reading of crux sketcher widget should be 0
    When user closes sketcher dialog
    Then sketcher dialog should be absent
    And the value of "molecule" column in row 1 should be, byte for byte, as remembered

  # HOST-054: spgi-100 has 17 molecules with a pyridine; the card's sketcher opens Crux in query mode
  Scenario: A pyridine drawn in Crux on a filter card keeps the 17 rows that contain it
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user clicks on "Sketch" text in "Structure" filter card
    Then sketcher dialog should be visible
    And the "mode" reading of crux sketcher widget should be "query"
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"

  # HOST-055: Filter as you draw is the account's (chem-filter): Chem's step that toggles it puts it back at feature end
  Scenario: With Filter as you draw cleared, a drawn edit reaches the grid only on OK, and Cancel brings back what the dialog opened with
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user clicks on "Sketch" text in "Structure" filter card
    And user unchecks "Filter as you draw" input in sketcher dialog
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    And the "structure of Structure" reading of filter panel should be ""
    And all rows should pass the filter
    When user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user clicks on the "card Structure" area of filter panel
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on crux single bond tool
    And user clicks on the "atom 3" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccncc1"
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And 17 rows should pass the filter
    When user clicks on CANCEL button in sketcher dialog
    Then sketcher dialog should be absent
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    And 17 rows should pass the filter
    When user clicks on the "card Structure" area of filter panel
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"

  # HOST-056: the thumbnail's Clear and the panel's reset both set an empty value on the card's sketcher
  Scenario: The thumbnail's Clear and the panel's reset empty the card, bring the Sketch link back, and its Crux opens empty
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user clicks on "Sketch" text in "Structure" filter card
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user hovers over sketcher thumbnail in "Structure" filter card
    And user clicks on "Clear" button in "Structure" filter card
    Then all rows should pass the filter
    And the "structure of Structure" reading of filter panel should be ""
    And "Sketch" text in "Structure" filter card should be visible
    When user clicks on "Sketch" text in "Structure" filter card
    Then the "ready" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0
    When user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    When user hovers over filter panel
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    And the "structure of Structure" reading of filter panel should be ""
    And "Sketch" text in "Structure" filter card should be visible
    When user clicks on "Sketch" text in "Structure" filter card
    Then the "ready" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0

  # HOST-057: smiles has 1,000 molecules, 924 of them with a benzene ring
  Scenario: Chem > Search > Substructure Search opens Crux empty and ready, and a ring drawn at once filters by it
    Given user opens smiles dataset
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then sketcher dialog should be visible
    And the "ready" reading of crux sketcher widget should be "true"
    And the "atoms" reading of crux sketcher widget should be 0
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on OK button in sketcher dialog
    Then 924 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "c1ccccc1"

  # HOST-058: the card filters by the cell's molblock (row 1 alone contains its own molecule); OK on an unedited Crux
  # gives it back as it was given
  Scenario: Use as filter shows the cell's molecule in Crux, and OK without an edit gives the card the same molblock back
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user hovers over "Structure" filter card
    And user clicks on close of "Structure" filter card
    And user picks "Current Value > Use as filter" from the context menu of the "cell 1 of Structure" area of grid
    Then there should be 1 visible "Structure" filter card
    And the "structure of Structure" reading of filter panel should be the molecule of row 1 of "Structure" column
    And 1 row should pass the filter
    When user remembers the "molecule of Structure" reading of filter panel
    And user clicks on the "card Structure" area of filter panel
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "Structure" column
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the "molecule of Structure" reading of filter panel should be as remembered
    And 1 row should pass the filter

  # HOST-059: the column popup's filter holds its sketcher inline, not in a dialog
  Scenario: In a column popup's filter, Crux draws inside the popup, filters, and keeps the keyboard
    Given user opens spgi-100 dataset
    When user opens an empty filter panel
    And user hovers over the "header Structure" area of grid
    And user clicks on "Column options" icon in grid
    Then column popup should be visible
    And crux sketcher widget in column popup should be visible
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then crux canvas should be focused
    And 17 rows should pass the filter
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccncc1"
    When user presses Control+Z
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    And column popup should be visible
    And the filter should pass exactly the molecules of "Structure" column containing "c1ccccc1"

  # HOST-060: a drawn query's card filters by Crux's own molblock; each round trip sets it back on a new sketcher
  Scenario: A drawn query's molblock comes back unchanged through Clone View, its sketcher, a layout, the panel reopened and a project
    Given user opens spgi-100 dataset
    When user clicks on filter icon in toolbar
    And user clicks on "Sketch" text in "Structure" filter card
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on OK button in sketcher dialog
    Then 17 rows should pass the filter
    And the "structure of Structure" reading of filter panel should be "c1ccncc1"
    When user remembers the "molecule of Structure" reading of filter panel
    And user picks "View > Layout > Clone View" from the top menu
    Then filter panel should be visible
    And the "molecule of Structure" reading of filter panel should be as remembered
    When user clicks on the "card Structure" area of filter panel
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And the "molecule of Structure" reading of filter panel should be as remembered
    And 17 rows should pass the filter
    When user saves the layout of the current table view
    And user loads the saved layout
    Then the "molecule of Structure" reading of filter panel should be as remembered
    And 17 rows should pass the filter
    When user clicks on close icon of filters viewer
    And user clicks on filter icon in toolbar
    Then the "molecule of Structure" reading of filter panel should be as remembered
    And 17 rows should pass the filter
    When user saves the current view as project "crux-chem-hosts-query"
    And user closes all views
    And user opens the "crux-chem-hosts-query" project
    Then "Structure" filter card should be visible
    And the "molecule of Structure" reading of filter panel should be as remembered
    And 17 rows should pass the filter

  # HOST-061: a scaffold added under a node is checked against it (the edited molecule must contain its parent)
  Scenario: In the Scaffold Tree's sketcher, a scaffold that does not contain its parent is flagged until it does
    Given user opens spgi-100 dataset
    When user picks "Chem > Analyze > Scaffold Tree" from the top menu
    Then Scaffold Tree viewer should be visible
    When user hovers over Scaffold Tree viewer
    And user clicks on "Sketch scaffolds manually" icon inside Scaffold Tree viewer
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 1
    When user hovers over the "node 1" area of Scaffold Tree viewer
    And user clicks on the "add icon of node 1" area of Scaffold Tree viewer
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on crux clear button
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccccc1"
    And sketcher dialog should contain the text "The edited molecule is not a superstructure of its parent"
    When user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on crux single bond tool
    And user clicks on the "atom 3" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "Cc1ccncc1"
    And sketcher dialog should not contain the text "The edited molecule is not a superstructure of its parent"
    When user clicks on "Add" button in sketcher dialog
    Then the "nodes" reading of Scaffold Tree viewer should be 2
    And the "scaffold of node 2" reading of Scaffold Tree viewer should be the molecule "Cc1ccncc1"

  # HOST-062: the viewer reruns the search on the sketcher's changes (search as you sketch, debounced)
  Scenario: In the Similarity Search viewer, a reference drawn in Crux reruns the search with the dialog open, and OK on an empty canvas is refused
    Given user opens smiles dataset
    When user picks "Chem > Search > Similarity Search..." from the top menu
    Then Chem Similarity Search viewer should be visible
    And the "target row" reading of Chem Similarity Search viewer should be 0
    When user remembers the "scores" reading of Chem Similarity Search viewer
    And user clicks on "Edit" icon in Chem Similarity Search viewer
    Then sketcher dialog should be visible
    And the "smiles" reading of crux sketcher widget should be the molecule in row 1 of "canonical_smiles" column
    When user clicks on crux clear button
    And user clicks on crux benzene tool
    And user clicks on crux canvas
    Then the "scores" reading of Chem Similarity Search viewer should not be as remembered
    And sketcher dialog should be visible
    When user clicks on crux clear button
    And user clicks on OK button in sketcher dialog
    Then an error balloon containing "Empty molecule cannot be used for similarity search" should have been shown

  # HOST-063: the pane writes the scaffold to the column, and the grid draws its structures aligned to it
  Scenario: A scaffold drawn in Crux in the column's Rendering pane realigns the structures, and its Clear realigns them back
    Given user opens smiles-50 dataset
    And the context panel is open
    When user clicks on the "header canonical_smiles" area of grid
    Then the context panel should show "canonical_smiles"
    When user expands Chemistry accordion header in context panel
    And user expands Rendering accordion header in context panel
    And user takes a snapshot of grid
    And user clicks on "Sketch" text in "Rendering" pane in context panel
    Then sketcher dialog should be visible
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on OK button in sketcher dialog
    Then the "cell 1 of canonical_smiles" area of grid should have repainted
    When user takes a snapshot of grid
    And user hovers over sketcher thumbnail in "Rendering" pane in context panel
    And user clicks on "Clear" button in "Rendering" pane in context panel
    Then the "cell 1 of canonical_smiles" area of grid should have repainted
    And "Sketch" text in "Rendering" pane in context panel should be visible

  # HOST-064: the pane's new fragment is green until another colour is picked
  Scenario: A benzene drawn in Crux in the Highlight pane colours the cells that contain it
    Given user opens a table "molecules" with:
      | molecule          |
      | C[C@H](N)c1ccccc1 |
      | C[C@H](N)C(=O)O   |
    And the semantic types of the current table are detected
    And the context panel is open
    When user clicks on the "header molecule" area of grid
    And user expands Chemistry accordion header in context panel
    And user expands Highlight accordion header in context panel
    And user clicks on "Sketch" text in "Highlight" pane in context panel
    Then sketcher dialog should be visible
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on OK button in sketcher dialog
    Then the "cell 1 of molecule" area of grid should contain the color "#00ff00"
    And the "cell 2 of molecule" area of grid should not contain the color "#00ff00"

  # RGROUP-009: R2 on the ring and R1 on the nitrogen, against the core's atom order: a core whose labels were lost
  # would be numbered by that order, and R1 would take the ring's substituent
  Scenario: A core drawn in Crux with R1 and R2 decomposes the table at its labels
    Given user opens a table "amides" with:
      | smiles                 |
      | Cc1ccc(cc1)C(=O)NCC    |
      | Clc1ccc(cc1)C(=O)NC    |
      | Oc1ccc(cc1)C(=O)NC(C)C |
    And the semantic types of the current table are detected
    When user picks "Chem > Analyze > R-Groups Analysis..." from the top menu
    Then "R-Groups Analysis" dialog should be visible
    When user types "Cc1ccc(cc1)C(=O)NC" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on crux R-group tool
    And user clicks on the "atom 10" area of crux sketcher widget
    And user clicks on crux R1 button
    And user clicks on crux R-Group OK button
    And user clicks on the "atom 0" area of crux sketcher widget
    And user clicks on crux R2 button
    And user clicks on crux R-Group OK button
    Then the "smiles" reading of crux sketcher widget should be the molecule "[*:2]c1ccc(cc1)C(=O)N[*:1]"
    When user clicks on OK button in "R-Groups Analysis" dialog
    Then a new column "R1" should have been added
    And the value of "R1" column in row 1 should be "CC[*:1]"
    And the value of "R2" column in row 1 should be "C[*:2]"
    And the value of "R1" column in row 3 should be "CC(C)[*:1]"

  # RGROUP-010: the editor checks the fragment for an R1 ("Fragment must contain R1-group": the input marked, OK disabled),
  # and the run removes what the fragment matches up to the R1's atom
  Scenario: A protecting group drawn in Crux with R1 in the Deprotect editor validates, and the run removes it
    Given user opens a table "amides" with:
      | molecule           |
      | CC(=O)Nc1ccc(O)cc1 |
      | CC(=O)NCCO         |
      | CCO                |
    And the semantic types of the current table are detected
    When user picks "Chem > Transform > Reactions > Deprotect..." from the top menu
    Then "Deprotect" dialog should be visible
    When user clicks on "Fragment" input in "Deprotect" dialog
    And user types "CC(=O)NC" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then OK button in "Deprotect" dialog should be disabled
    When user clicks on "Fragment" input in "Deprotect" dialog
    And user types "CC(=O)NC" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on crux R-group tool
    And user clicks on the "atom 4" area of crux sketcher widget
    And user clicks on crux R1 button
    And user clicks on crux R-Group OK button
    Then the "smiles" reading of crux sketcher widget should be the molecule "CC(=O)N[*:1]"
    When user clicks on OK button in sketcher dialog
    Then OK button in "Deprotect" dialog should be enabled
    When user clicks on OK button in "Deprotect" dialog
    Then a new column "deprotected(molecule)" should have been added
    And the molecule in row 1 of "deprotected(molecule)" column should be "Oc1ccccc1"
    And the molecule in row 2 of "deprotected(molecule)" column should be "CCO"
    And the molecule in row 3 of "deprotected(molecule)" column should be "CCO"
