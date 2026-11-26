@guide @tutorial
Feature: How-to · Molecules
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a demo
  clip that outlines what each step presses). Finding related molecules three ways: Chem > Search > Similarity
  Search... adds a viewer of the molecules most similar to the current one (Tanimoto on Morgan fingerprints, with the
  score on each card), which follows the current row; a click beside a card's drawing makes that molecule current;
  Chem > Search > Diversity Search... adds a viewer of molecules picked to differ from each other; Chem > Analyze >
  Chemical Space... (UMAP on fingerprint distances) adds two embedding columns and a cluster column and plots every
  molecule so that similar ones land close together. The pointer on a dot shows its structure, a click on a dot makes
  that molecule current (and the similarity viewer follows it), and the plot's Color selector colors the map by any
  column. All of it runs in the browser (RDKit). smiles (1000 molecules) is opened off camera with the grid narrowed
  to molregno and the structures; the new columns live in the page's table only and nothing is made on the server.
  Merged from the former clips "Find similar and diverse molecules" and "Map the chemical space".

  Scenario: Find similar and diverse molecules, and map the chemical space
    Find the molecules most like one, a set unlike each other, and a map where similar structures sit together.

    Given user is logged in
    And simple mode is off
    And the "Chem" package is installed
    And the package autostarts have completed
    And user opens smiles dataset
    And the grid shows only the columns "molregno, canonical_smiles"
    And the toolbox pane is hidden
    And the context panel is closed
    And the browse panel is closed
    # caption: Chem > Search > Similarity Search
    When user picks "Chem > Search > Similarity Search..." from the top menu
    Then Chem Similarity Search viewer should be visible
    When the video says "The current molecule first, then the most similar, with scores"
    And user remembers the "scores" reading of Chem Similarity Search viewer
    # caption: Make another molecule current: the search follows
    And user clicks on the "cell 5 of canonical_smiles" area of grid
    Then row 5 should be current
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    # caption: Click a card to make its molecule current
    When user clicks on card 2 of Chem Similarity Search viewer beside its drawing
    Then the current row should be the row of the clicked card
    # caption: Chem > Search > Diversity Search
    When user picks "Chem > Search > Diversity Search..." from the top menu
    Then Chem Diversity Search viewer should be visible
    When the video says "Molecules picked to be as different from each other as possible"
    # caption: Chem > Analyze > Chemical Space puts every molecule on a map
    And user picks "Chem > Analyze > Chemical Space..." from the top menu
    Then "Chem Space" dialog should be visible
    And Method input in "Chem Space" dialog should have value "UMAP"
    # caption: UMAP on fingerprint similarity, then OK
    When user clicks on OK button in "Chem Space" dialog
    Then the top menu command should have completed
    And a new column matching "^Embed_X_" should have been added
    And scatter plot viewer should be visible
    And scatter plot viewer should be painted
    When the video says "Each dot is a molecule; similar ones sit close together"
    # caption: Point at a dot to see its structure
    And user hovers over the first free "marker" area of scatter plot viewer
    Then tooltip should be visible
    When user remembers the "scores" reading of Chem Similarity Search viewer
    # caption: Click a dot: the similarity search follows it
    And user clicks on the "marker of row 15" area of scatter plot viewer
    Then row 15 should be current
    And the "scores" reading of Chem Similarity Search viewer should not be as remembered
    # caption: Color the map by any column
    When user picks "NumAromaticCarbocycles" in the "color" column selector of scatter plot viewer
    Then "Color Column Name" property of scatter plot viewer should be "NumAromaticCarbocycles"
    And the table should have 1000 rows
    And no errors should have been logged
