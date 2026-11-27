@guide @tutorial @viewport:1280x720
Feature: How-to · Molecules
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). Finding and highlighting a substructure, the two ways side by side:
  Chem > Search > Substructure Search... opens the sketcher and puts a structure card for the molecule column on the
  filter panel; a fragment (a carboxylic acid, C(=O)O) confirmed with OK keeps the 314 of 1000 molecules that contain
  it, the match drawn in each. The filter panel's reset icon brings every row back. Then the column's own Chemistry >
  Highlight pane in the context panel takes a second fragment (pyridine, c1ccncc1) from the sketcher and colors it in
  every molecule that has it, with no row filtered out. Both run in the browser (RDKit). smiles (1000 molecules) is
  opened off camera with the grid narrowed to molregno and the structures, and the OpenChemLib sketcher is pinned for
  the feature and put back after it; the filter and the highlight live on the page's table only, nothing is made on
  the server. Merged from the former clips "Find molecules that contain a fragment" and "Highlight a fragment in
  every molecule".

  Scenario: Find and highlight a substructure
    Keep only the molecules that contain a fragment, or color a fragment in every molecule without hiding any.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the molecule sketcher is "OpenChemLib"
    And the package autostarts have completed
    And user opens smiles dataset
    And the grid shows only the columns "molregno, canonical_smiles"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Chem > Search > Substructure Search
    When user picks "Chem > Search > Substructure Search..." from the top menu
    Then sketcher dialog should be visible
    And all rows should pass the filter
    When the video says "Draw the fragment, or type its SMILES" at the top
    And user clicks on molecule input of sketcher dialog
    And user types "C(=O)O" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then 314 rows should pass the filter
    And the "structure of canonical_smiles" reading of filter panel should be "C(=O)O"
    When the video says "314 of 1000 molecules contain it; the match is drawn in each"
    # caption: The fragment stays on the filter panel; click it to change it
    And user hovers over the "card canonical_smiles" area of filter panel
    # caption: Reset brings every molecule back
    And user clicks on reset icon of filter panel
    Then all rows should pass the filter
    When the video says "To mark a fragment without hiding any rows, highlight it"
    # caption: Click the column's header, then open the context panel
    And user clicks on the "header canonical_smiles" area of grid
    And user clicks on "Context Panel" status bar toggle
    Then the context panel should show "canonical_smiles"
    # caption: Chemistry > Highlight > Sketch
    When user expands Chemistry accordion header in context panel
    And user expands Highlight accordion header in context panel
    And user clicks on "Sketch" text in "Highlight" pane in context panel
    Then sketcher dialog should be visible
    # caption: Here a pyridine ring
    When user clicks on molecule input of sketcher dialog
    And user types "c1ccncc1" into molecule input of sketcher dialog
    And user presses Enter in molecule input of sketcher dialog
    And user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    And all rows should pass the filter
    When the video says "The ring is colored wherever it occurs; all 1000 rows stay"
    Then the table should have 1000 rows
    And no errors should have been logged
