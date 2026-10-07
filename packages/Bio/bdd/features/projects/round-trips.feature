@realizes:bio.project.round-trip
Feature: Bio results survive a project save and reopen
  A table saved as a project with what Bio added to it — an antibody numbering's aligned column —
  comes back from the server with the columns and their tags; the numbering run again on the
  reopened table loads its WASM engine afresh and numbers exactly as before.

  Sequence Space's embeddings and scatter plot after a reopen (GROK-19928) are tested in Bio
  src/tests/projects-tests.ts (sequence_space), and a reopened HELM column's detection in
  src/tests/detectors-tests.ts (samplesHelmCsv).

  Not translated: the ribbon Save dialog and its Data Sync toggle — the projects are saved through
  the project API with their data uploaded, because a feature's table is a clone with no file
  behind it for Data Sync to follow (the FASTA file feature does the same for an imported file);
  a project "referencing" a monomer library or collection (the lifecycle md cases) — a project
  stores neither, so what those cases can claim is that the library still colours the reopened
  column, which is the renderer's and tested by the package, and that the collection files are
  untouched, which the collections feature owns.

  The save step clears what an earlier run left under the same name itself, so no scenario asks
  for that separately (each sweep lists every project on the stand, about half a minute on dev).

  Background:
    Given user is logged in
    And the Bio package is initialized

  Scenario: Antibody numbering survives a save and reopen, and numbers the same again
    Given user opens antibodies dataset keeping the first 20 rows as "antibodies"
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user selects "kabat" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
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
    Then a new column "AntibodyHC (aligned) (2)" should have been added
    And "AntibodyHC (aligned) (2)" column should hold the same values as "AntibodyHC (aligned)" column
    And no error or warning balloon should have been shown
    And no errors should have been logged
