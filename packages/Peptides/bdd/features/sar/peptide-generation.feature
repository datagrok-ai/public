Feature: The Peptide Generation viewer designs new peptides and opens them as a table
  Peptide Generation is an optional viewer (help `visualize/viewers/peptide-generation.md`). It ranks generated
  peptides by predicted activity toward the Activity target and says how many it generated above its grid;
  Peptide count bounds that number, Activity target switches between the highest and the lowest activity, and its
  plus icon opens the generated peptides as a new table, one row per peptide.

  Runs on the peptides dataset (647 aligned peptides, IC50). Left out on purpose: the predictions, the beam search,
  the basis columns and their cells, Novel (Exclude existing) and the generated table's columns are the package's
  'Viewers: Peptide Generation' tests and the cell renderer's (MISSING.md). Nothing is put on the server.

  Scenario: The viewer generates peptides for the target and count set, and its plus icon opens them as a table
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    And user adds Peptide Generation viewer
    Then Peptide Generation viewer should contain text "50 peptides generated for highest activity"
    When user sets "activityColumnName" property of Peptide Generation viewer to "IC50"
    And user sets "activityTarget" property of Peptide Generation viewer to "Low"
    And user sets "peptideCount" property of Peptide Generation viewer to "10"
    Then Peptide Generation viewer should contain text "10 peptides generated for lowest activity"
    When user clicks on "Open the generated peptides as a new table" icon in Peptide Generation viewer
    Then the "Generated Peptides" table view should open with 10 rows
    And no errors should have been logged
