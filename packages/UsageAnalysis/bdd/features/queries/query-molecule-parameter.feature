@crux @sketcher-controls
Feature: A query's molecule parameter, sketched in Crux
  A database query's parameter of semantic type Molecule is the platform's molecule input in the query's run dialog: a
  click on it opens a sketcher dialog, and its OK writes the sketcher's SMILES into the parameter, which the query then
  receives: SMILES, as `mol_from_smiles` and the like expect, never a molblock. Here the sketcher is Crux, Chem's own,
  pinned as the session's sketcher, and the pattern is drawn on Crux's own controls.

  The query is this feature's own, on the platform's own database (System:Datagrok), and only hands its parameter back
  as its one row; it is deleted with its chats at the end. Needs Chem on the stand. The feature drives Crux's own
  controls (@sketcher-controls): a run that pins another sketcher skips it.

  Background:
    Given user is logged in
    And the "Chem" package is installed
    And the molecule sketcher is "Crux"
    And the package autostarts have completed
    And a query "BddCruxPattern" on "System:Datagrok" is on the server:
      """
      --input: string pattern {semType: Molecule}
      select @pattern as pattern
      """
    And the browse panel is open

  # HOST-066
  Scenario: A pattern drawn in Crux for the query's Molecule parameter reaches the query as SMILES
    Given Databases tree node inside browse tree is expanded
    And Databases---Postgres tree node inside browse tree is expanded
    And Databases---Postgres---Datagrok tree node inside browse tree is expanded
    When user picks "Run" from the context menu of Databases---Postgres---Datagrok---BddCruxPattern tree node inside browse tree
    Then "BddCruxPattern" dialog should be visible
    When user clicks on editor of pattern input in "BddCruxPattern" dialog
    Then sketcher dialog should be visible
    And the "ready" reading of crux sketcher widget should be "true"
    When user clicks on crux benzene tool
    And user clicks on crux canvas
    And user clicks on crux nitrogen tool
    And user clicks on the "atom 0" area of crux sketcher widget
    Then the "smiles" reading of crux sketcher widget should be the molecule "c1ccncc1"
    When user clicks on OK button in sketcher dialog
    Then sketcher dialog should be absent
    When user clicks on OK button in "BddCruxPattern" dialog
    Then the current view should be a TableView view
    And the table should have 1 row
    And every value of "pattern" column should match "^\S+$"
    And the molecule in row 1 of "pattern" column should be "c1ccncc1"
