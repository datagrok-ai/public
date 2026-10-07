@journey @realizes:bio.menu.polytool.convert @realizes:bio.menu.polytool.combine-sequences
Feature: Bio | PolyTool on a custom-notation column: Convert and Combine Sequences
  cyclized.csv holds 14 cyclic peptides in a custom notation (`seqs`, e.g.
  R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI). PolyTool Convert with Get HELM on turns it into a
  HELM column and a molfile column; Combine Sequences refuses to run with no table chosen and, with
  seqs combined with itself, opens a new table of every pair. Translated from the TestTrack case
  SequenceTranslator/oligo-nucleotide-grid, Block G.

  The table is opened from Browse > Files > App Data, as a user does: the custom notation is set by
  SequenceTranslator's notation refiner when the file's columns are detected.

  Background:
    Given user is logged in
    And the package autostarts have completed
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator tree node inside browse tree is expanded
    And Files---App-Data---SequenceTranslator---samples tree node inside browse tree is expanded
    When user double-clicks Files---App-Data---SequenceTranslator---samples---cyclized.csv tree node inside browse tree
    Then the "cyclized" view should be current
    And the table should have 14 rows
    And "seqs" column should have semantic type "Macromolecule"

  Scenario: The PolyTool submenu offers Convert, Enumerate HELM and Combine Sequences
    Then the top menu should list:
      | Bio > PolyTool > Convert...           |
      | Bio > PolyTool > Enumerate HELM...    |
      | Bio > PolyTool > Combine Sequences... |

  Scenario: Convert with Get HELM adds a HELM column and a molfile column
    When user picks "Bio > PolyTool > Convert..." from the top menu
    Then "PolyTool Conversion" dialog should be visible
    And Column input in "PolyTool Conversion" dialog should contain the text "seqs"
    And "Get HELM" checkbox in "PolyTool Conversion" dialog should be checked
    When user clicks on OK button in "PolyTool Conversion" dialog
    Then the "PolyTool Conversion" dialog should close
    And a new column "transformed(seqs)" should have been added
    And a new column "molfile(seqs)" should have been added
    And "transformed(seqs)" column should have semantic type "Macromolecule"
    And "transformed(seqs)" column should have units "helm"
    And "molfile(seqs)" column should have semantic type "Molecule"
    And the value of "transformed(seqs)" column in row 1 should be "PEPTIDE1{R.F.C.T.G.H.F.Y.G.H.F.Y.G.H.F.Y.P.C.[meI]}$PEPTIDE1,PEPTIDE1,3:R3-18:R3$$$V2.0"
    And "transformed(seqs)" column should have no missing values
    And "molfile(seqs)" column should have no missing values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Combine Sequences with no table chosen says so and makes no table
    When user picks "Bio > PolyTool > Combine Sequences..." from the top menu
    Then "Combine Sequences" dialog should be visible
    When user clicks on OK button in "Combine Sequences" dialog
    Then an error balloon containing "Please fill all the fields" should have been shown
    And the open table views should be exactly "cyclized"

  Scenario: Combine Sequences of seqs with itself opens a table of all 196 pairs
    When user picks "Bio > PolyTool > Combine Sequences..." from the top menu
    Then "Combine Sequences" dialog should be visible
    When user selects "cyclized" in Table input in "Combine Sequences" dialog
    Then Column input in "Combine Sequences" dialog should have the value "seqs"
    When user clicks on Add icon in "Combine Sequences" dialog
    And user selects "cyclized" in second Table input in "Combine Sequences" dialog
    Then second Column input in "Combine Sequences" dialog should have the value "seqs"
    When user enters "-" into Separator input in "Combine Sequences" dialog
    And user clicks on OK button in "Combine Sequences" dialog
    Then the "Combined Sequences" view should be current
    And the table should have 196 rows
    And the table should have the columns "Combined Sequences"
    And the value of "Combined Sequences" column in row 1 should be "R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI-R-F-C(1)-T-G-H-F-Y-G-H-F-Y-G-H-F-Y-P-C(1)-meI"
    And no error or warning balloon should have been shown
    And no errors should have been logged
