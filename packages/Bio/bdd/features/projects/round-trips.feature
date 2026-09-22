@serial @realizes:GROK-19928 @realizes:bio.project.round-trip
Feature: Bio results survive a project save and reopen
  A table saved as a project with what Bio added to it — Sequence Space's embeddings and the
  scatter plot over them (GROK-19928), an antibody numbering's aligned column — comes back from
  the server with the columns, their tags and the viewers; the numbering run again on the reopened
  table loads its WASM engine afresh and numbers exactly as before. A HELM table comes back
  painted by the HELM renderer from the same monomer library. Serial: the HELM scenario reads the
  monomer library selection, which the library features toggle.

  Not translated: the ribbon Save dialog and its Data Sync toggle — the projects are saved through
  the project API with their data uploaded, because a feature's table is a clone with no file
  behind it for Data Sync to follow (the FASTA file feature does the same for an imported file);
  a project "referencing" a monomer library or collection (the lifecycle md cases) — a project
  stores neither, so what those cases can claim is that the library still colours the reopened
  column, claimed here, and that the collection files are untouched, which the collections feature
  owns.

  The save step clears what an earlier run left under the same name itself, so no scenario asks
  for that separately (each sweep lists every project on the stand, about half a minute on dev).

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Sequence Space output survives a save and reopen (GROK-19928)
    Given user opens FASTA_sample dataset keeping the first 20 rows as "FASTA_sample"
    When user picks "Bio > Analyze > Sequence Space..." from the top menu
    Then "Sequence Space" dialog should be visible
    When user clicks on OK button in "Sequence Space" dialog
    Then the top menu command should have completed
    And a new column "Embed_X_1" should have been added
    And the open tableview should have 1 scatter plot viewer
    When user saves the current view as project "bdd-bio-seqspace-{run}"
    And user closes all views
    And user opens the "bdd-bio-seqspace-{run}" project
    Then the table should have 20 rows
    And "Sequence" column should have semantic type "Macromolecule"
    And "Sequence" column should have units "fasta"
    And "Embed_X_1" column should have no missing values
    And "Embed_Y_1" column should have no missing values
    And the open tableview should have 1 scatter plot viewer
    And scatter plot viewer should be painted
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: Antibody numbering survives a save and reopen, and numbers the same again
    Given user opens antibodies dataset keeping the first 20 rows as "antibodies"
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user selects "kabat" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then the top menu command should have completed
    And a new column "AntibodyHC (aligned)" should have been added
    When user saves the current view as project "bdd-bio-numbering-{run}"
    And user closes all views
    And user opens the "bdd-bio-numbering-{run}" project
    Then the table should have 20 rows
    And "AntibodyHC" column should have units "fasta"
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "kabat"
    And the ".positionNames" tag of "AntibodyHC (aligned)" column should list at least 100 values
    And "AntibodyHC (aligned)" column should have no missing values
    And "AntibodyHC" column should carry at least 7 annotations
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user selects "kabat" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then the top menu command should have completed
    And a new column "AntibodyHC (aligned) (2)" should have been added
    And "AntibodyHC (aligned) (2)" column should hold the same values as "AntibodyHC (aligned)" column
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: A HELM table comes back painted from the same monomer library
    Given user opens filter_HELM dataset
    Then the "cell type of HELM string" reading of grid should be "helm"
    And the "cell 1 of HELM string" area of grid should be painted in at least 3 colors
    When user saves the current view as project "bdd-bio-helm-{run}"
    And user closes all views
    And user opens the "bdd-bio-helm-{run}" project
    Then the table should have 4 rows
    And "HELM string" column should have semantic type "Macromolecule"
    And "HELM string" column should have units "helm"
    And the "cell type of HELM string" reading of grid should be "helm"
    And the "cell 1 of HELM string" area of grid should be painted in at least 3 colors
    And the monomer library should be loaded from "HELMCoreLibrary.json"
    And "dV" should be a known "PEPTIDE" monomer
    And no error or warning balloon should have been shown
    And no errors should have been logged
