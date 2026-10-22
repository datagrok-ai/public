@journey @realizes:bio.import.fasta @realizes:bio.export.fasta @realizes:GROK-18616
Feature: FASTA files: open, export, reopen
  A .fasta file opened from the computer goes through Bio's file handler and comes up as a table
  whose sequence column is already a detected Macromolecule, painted by the sequence renderer and
  offered by the analysis dialogs (GROK-18616: an entry path that skipped the detector left
  Sequence Space with no column to offer). Download > As FASTA... writes the table back as FASTA
  and that file opens again as the same sequences.

  Translated from TestTrack Bio/bio-lifecycle-fasta-file.md: scenario 1 and scenario 2 (the entry
  paths) through Open local file, a drop onto the window and a double-click in Browse > Files, each
  ending in the same detected sequence table; scenario 3 (the FASTA round trip); scenario 4 (the
  imported table saved to a project and reopened, its semantic type, units and renderer back). The
  drop is the browser's own drag of the file from disk, made through the DevTools protocol: the
  platform lays its drop layer over the window and the file is dropped on it. The double-clicked
  file is the same fixture written into My files for the feature.

  Background:
    Given user is logged in
    And simple mode is off
    And the browse panel is open
    And the Bio package is initialized

  Scenario: A FASTA file opened from the computer is a detected sequence table
    When user uploads "fixtures/bdd-sample.fasta" through "Open local file" icon inside browse toolbar
    Then the "bdd-sample" view should be current
    And the table should have 6 rows
    And the table should have 2 columns
    And the value of "description" column in row 1 should be "UPI0000000595:31"
    And the value of "sequence" column in row 1 should be "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    And the "cell type of sequence" reading of grid should be "sequence"
    And the "cell 1 of sequence" area of grid should be painted in at least 3 colors
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    And editor of Column input in "Sequence Space" dialog should have text "sequence"
    When user clicks on CANCEL button in "Sequence Space" dialog
    Then "Sequence Space" dialog should be hidden
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Download as FASTA writes the table back, and the file opens as the same sequences
    When user renames "description" column to "seq id"
    And user clicks on "arrow to bottom" icon
    And user clicks on "As FASTA..." label
    Then "Save as FASTA" dialog should be visible
    When user downloads a file through OK button in "Save as FASTA" dialog
    Then the downloaded file should contain ">UPI0000000595:31"
    And the downloaded file should contain "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"
    And the downloaded file should contain ">UPI000005175B:38"
    When user clicks on "Browse" view
    When user uploads the downloaded file through "Open local file" icon inside browse toolbar
    Then the table should have a column "description"
    And the table should not have a column "seq id"
    And the table should have 6 rows
    And the value of "description" column in row 1 should be "UPI0000000595:31"
    And the value of "sequence" column in row 1 should be "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"
    And the value of "sequence" column in row 6 should be "MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF"
    And "sequence" column should have units "fasta"
    And the "cell type of sequence" reading of grid should be "sequence"
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A FASTA file dropped onto the window is a detected sequence table
    When user closes all views
    Then no table view is open
    When user drops the "fixtures/bdd-sample.fasta" file of the project onto browse tree
    Then the open table views should be exactly "bdd-sample"
    And the table should have 6 rows
    And the value of "sequence" column in row 1 should be "MDYKETLLMPKTDFPMRGGLPNKEPQIQEKW"
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    And the "cell type of sequence" reading of grid should be "sequence"
    And the "cell 1 of sequence" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A FASTA-imported table keeps its type, units and renderer through a project
    When user saves the current view as project "bdd-bio-fasta-{run}"
    And user closes all views
    And user opens the "bdd-bio-fasta-{run}" project and waits for its table
    Then the table should have 6 rows
    And the value of "sequence" column in row 6 should be "MHAILRYFIRRLFYHIFYKIYSLISKKHQSLPSDVRQF"
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    And the "cell type of sequence" reading of grid should be "sequence"
    And the "cell 1 of sequence" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A FASTA file double-clicked in Browse > Files is a detected sequence table
    Given the "fixtures/bdd-sample.fasta" file of the project is in the user's files as "bdd-fasta-{time}.fasta"
    When user closes all views
    And user refreshes the browse tree
    And Files tree node inside browse tree is expanded
    And Files---My-files tree node inside browse tree is expanded
    And user double-clicks on Files---My-files---bdd-fasta-{time}.fasta tree node inside browse tree
    Then the "bdd-fasta-{time}" table view should open with 6 rows
    And "sequence" column should have semantic type "Macromolecule"
    And "sequence" column should have units "fasta"
    And the "cell type of sequence" reading of grid should be "sequence"
    And the "cell 1 of sequence" area of grid should be painted in at least 3 colors
    And no error or warning balloon should have been shown
    And no errors should have been logged
