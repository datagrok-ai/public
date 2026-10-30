@journey
Feature: The Peptide Generation viewer designs new peptides and opens them as a table
  Peptide Generation is an optional viewer added from the Add Viewer dialog (help
  `visualize/viewers/peptide-generation.md`: "Add it like any other viewer from the Add viewer menu").
  It ranks generated peptides by predicted activity toward the Activity target and says how many it
  generated above its grid; Peptide count bounds that number, Activity target switches between the
  highest and the lowest activity, and its plus icon opens the generated peptides as a new table
  ("Peptide", "Predicted activity", ..., a basis column per position of the sequences, one row per peptide; "Novel" is true for each with Exclude
  existing on, the default).

  Runs on the peptides dataset (647 aligned peptides, IC50). Left out on purpose:
  - the predicted values, the beam search and the per-position basis numbers: the package's
    `src/tests/viewers.ts` ("Viewers: Peptide Generation") is where the algorithm matrix belongs;
  - the basis cells' circles and colours: what a cell renderer draws (bdd CLAUDE.md);
  - the viewer's default Activity column (it takes ID, not IC50, on this dataset): a question in MISSING.md;
    the feature picks IC50 itself.
  Nothing is put on the server.

  Background:
    Given user is logged in
    And simple mode is off
    And the package autostarts have completed
    And the Peptides package is initialized
    And user opens peptides dataset

  Scenario: The Add Viewer dialog adds the viewer, which generates 50 peptides for the highest activity
    Then "Add viewer" icon in toolbar should be visible
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "peptide" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Peptide Generation" card in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And the open tableview should have 1 Peptide Generation viewer
    And Peptide Generation viewer should contain text "50 peptides generated for highest activity"
    And no errors should have been logged

  Scenario: Activity, Activity target and Peptide count in the properties change what is generated
    When user picks "Properties..." from the viewer menu of Peptide Generation viewer
    And user selects "IC50" in "Activity" property in context panel
    And user expands "Misc" category in context panel
    And user selects "Low" in "Activity Target" property in context panel
    Then Peptide Generation viewer should contain text "50 peptides generated for lowest activity"
    When user enters "10" into "Peptide Count" property in context panel
    Then Peptide Generation viewer should contain text "10 peptides generated for lowest activity"
    And no errors should have been logged

  Scenario: The plus icon opens the generated peptides as a new table
    When user clicks on "Open the generated peptides as a new table" icon in Peptide Generation viewer
    Then the "Generated Peptides" table view should open with 10 rows
    And the table should have the columns "Peptide, Predicted activity, Confidence, Min support, Mean p-value, Significant positions, Novel, Basis, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17"
    And every value of "Novel" column should match "^true$"
    And table "peptides" should have 647 rows
    And no errors should have been logged
