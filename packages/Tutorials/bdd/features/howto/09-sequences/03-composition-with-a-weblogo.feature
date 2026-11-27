@guide @tutorial
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). On FASTA_PT_activity (99 peptides of 16 residues in three
  families, with an activity), Bio > Analyze > Composition docks a WebLogo over the sequence column: a stack of
  letters per position, the pointer on a letter shows how many peptides carry it, and a click on it selects them: 31
  with N at position 3, all between 4.37 and 4.49, near the top of the activity (-2.81 to 4.93), then 34 with V at
  position 14, 33 of them below 2 and the lowest of all among them; after
  each click the grid is scrolled down to the rows it selected. Off camera, before the first step: the histogram of
  the activity added beside the grid. Nothing is made on the server.

  Scenario: See the composition with a WebLogo
    A WebLogo shows which residues sit at each position; a click on a letter selects those sequences.

    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "sequence, activity"
    And the "sequence" column of the grid is 300 px wide
    And user adds a histogram viewer with:
      | valueColumnName | activity |
    # caption: Bio > Analyze > Composition
    When user picks "Bio > Analyze > Composition" from the top menu
    Then WebLogo viewer should be visible
    And WebLogo viewer should be painted
    When the video says "A WebLogo: the residues at each position, the most common on top" at the top
    And user hovers over the "monomer N at position 3" area of WebLogo viewer
    And the video says "Click a letter to select the peptides that have it there" at the top
    And user clicks on the "monomer N at position 3" area of WebLogo viewer
    Then 31 rows should be selected
    And histogram viewer should show a selection highlight
    # caption: Scroll the table to them
    When user scrolls the mouse wheel down 3 times over the "cell 10 of sequence" area of grid
    Then the "row header 40" area of grid should be at least 10 pixels tall
    When the video says "Selected everywhere: in the histogram, all near the top of the activity" at the top
    And user hovers over the "monomer V at position 14" area of WebLogo viewer
    And the video says "Another residue, another family: nearly all at the low end" at the top
    And user clicks on the "monomer V at position 14" area of WebLogo viewer
    Then 34 rows should be selected
    And histogram viewer should show a selection highlight
    When user scrolls the mouse wheel down 3 times over the "cell 45 of sequence" area of grid
    Then the "row header 70" area of grid should be at least 10 pixels tall
    When user hovers over the first "selected bin" area of histogram viewer
    Then no error or warning balloon should have been shown
    And no errors should have been logged
