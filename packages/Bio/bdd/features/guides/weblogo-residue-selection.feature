@guide @help:visualize/viewers
Feature: Pick out the peptides with one residue at one position
  A guide: the answer to "how do I pick out every peptide that has a given amino acid at a given
  position, and see whether those peptides are the active ones?". Bio > Analyze > Composition
  docks a WebLogo over the sequence column: a stack of letters per position, each letter as tall
  as its share of the peptides. A letter is also a selection: a click on it selects the rows that
  have that monomer at that position, and every other viewer highlights them, so a histogram of the
  activity shows where they fall. Demo: FASTA_PT_activity, 99 peptides of 16 residues with an
  activity; the N at position 3 is carried by 31 of them, whose activity averages 4.4 against 3.0
  for the whole set.

  Scenario: Click a letter of the WebLogo to select the peptides that have it
    Given user is logged in
    And simple mode is off
    And user opens FASTA_PT_activity dataset
    And the Bio package is initialized
    When user picks "Bio > Analyze > Composition" from the top menu
    Then WebLogo viewer should be visible
    When user opens toolbox
    And user clicks on histogram icon on toolbox
    And user selects "activity" in Value column input in histogram viewer
    And user clicks on the "monomer N at position 3" area of WebLogo viewer
    Then 31 rows should be selected
    And only the rows with "N" at position 3 of "sequence" column should be selected
    And histogram viewer should show a selection highlight
