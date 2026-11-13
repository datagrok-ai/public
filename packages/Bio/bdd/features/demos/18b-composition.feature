@guide @demo
Feature: Bioinformatics · Composition
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). On FASTA_PT_activity (99 peptides of 16 residues in
  three families, with an activity): Bio > Analyze > Composition docks a WebLogo over the sequence
  column, a stack of letters per position; the pointer on a letter shows how many peptides carry it,
  and a click on it selects them: 31 with N at position 3, which the histogram of the activity shows
  at the top of its range, then 34 with V at position 14, at its bottom. The same interaction as the
  guide features/guides/weblogo-residue-selection, filmed as a demo. Off camera: the histogram added beside the grid, the grid scrolled to the rows
  the click selects, and the WebLogo's letters drawn by frequency (Position Height 100%) rather than
  by information content, so every position reads at full height.

  Scenario: Pick out the peptides with one residue at one position
    A WebLogo stacks each position's residues by frequency; a click on a letter selects the peptides that carry it, and the other viewers follow.

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
    And the grid is scrolled down to row 52
    When user picks "Bio > Analyze > Composition" from the top menu
    Then WebLogo viewer should be visible
    Given user sets "positionHeight" property of WebLogo viewer to "100%"
    Then WebLogo viewer should be painted
    When the video says "Click a residue at a position: its peptides light up everywhere"
    And user hovers over the "monomer N at position 3" area of WebLogo viewer
    And user clicks on the "monomer N at position 3" area of WebLogo viewer
    Then 31 rows should be selected
    And only the rows with "N" at position 3 of "sequence" column should be selected
    And histogram viewer should show a selection highlight
    When user hovers over the "monomer V at position 14" area of WebLogo viewer
    And the video says "Another residue, another family: here the least active"
    And user clicks on the "monomer V at position 14" area of WebLogo viewer
    Then 34 rows should be selected
    And only the rows with "V" at position 14 of "sequence" column should be selected
    And histogram viewer should show a selection highlight
    And no error or warning balloon should have been shown
    And no errors should have been logged
