@journey @realizes:bio.annotate.numbering-scheme
Feature: Antibody numbering with the bundled immunum engine
  Bio | Annotate | Apply Numbering Scheme... runs a numbering engine in a worker: the dialog
  offers the engines registered and the schemes the chosen one supports, and the run leaves an
  aligned column tagged with the scheme and the position names, plus the FR and CDR regions as
  annotations on the source column. Kabat, chosen over the default IMGT, must reach the result, and
  the default IMGT run, made next on the same table, numbers by IMGT.
  Engines of other packages (e.g. HMM's ANARCI) come first and are the default when installed, so
  the scenario picks the bundled engine.
  The last scenario translates TestTrack Bio/bio-lifecycle-immunum-wasm.md scenarios 2 and 3: the
  table with the IMGT numbering is saved through the ribbon's Save project dialog and reopened, the
  aligned column comes back with its scheme and position-name tags and the source column with its
  region annotations, and the engine numbers the reopened table again. That the re-run gives the
  same values is tested in Bio src/tests/projects-tests.ts ('antibody_numbering').

  Not translated: the md's timing bound on the worker call — a duration is not a claim the
  platform signals; the other schemes of the md title (Chothia, AHo) — the bundled engine offers
  IMGT and Kabat only, which the dialog claim pins; the engine function's own result table, called
  directly (no UI) — see the bdd library's CLAUDE.md, "What never becomes a feature".

  Background:
    Given user is logged in
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"
    And "AntibodyHC" column should have units "fasta"

  Scenario: The dialog offers the engine and its schemes, IMGT first
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    And editor of Sequence input in "Apply Antibody Numbering" dialog should have text "AntibodyHC"
    When user selects "Immunum" in Engine input in "Apply Antibody Numbering" dialog
    Then Engine input in "Apply Antibody Numbering" dialog should have value "Immunum"
    And Scheme input in "Apply Antibody Numbering" dialog should offer "imgt, kabat"
    And Scheme input in "Apply Antibody Numbering" dialog should have value "imgt"

  Scenario: Kabat numbering aligns the column and annotates its regions
    When user selects "kabat" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "kabat"
    And every value of "AntibodyHC (aligned)" column should have the same length
    And "AntibodyHC (aligned)" column should have no missing values
    And every value of "AntibodyHC (aligned)" column should match "[A-Z]{20}"
    And the ".positionNames" tag of "AntibodyHC (aligned)" column should list at least 100 values
    And "AntibodyHC" column should carry at least 7 annotations
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The default IMGT numbering aligns the column by IMGT
    Given user removes "AntibodyHC (aligned)" column
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    When user selects "Immunum" in Engine input in "Apply Antibody Numbering" dialog
    Then Scheme input in "Apply Antibody Numbering" dialog should have value "imgt"
    When user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    And every value of "AntibodyHC (aligned)" column should have the same length
    And "AntibodyHC (aligned)" column should have no missing values
    And the ".positionNames" tag of "AntibodyHC (aligned)" column should list at least 100 values
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The numbering comes back from a saved project, and the engine runs again on it
    Given no project named "bdd-bio-numbering-{run}" is on the server
    When user opens the Save project dialog from the ribbon
    And user types "bdd-bio-numbering-{run}" into text input in "Save project" dialog
    And user clicks on OK in the Save project dialog and the project uploads
    Then the "Save project" dialog should close
    When user presses Escape
    Then 1 project named "bdd-bio-numbering-{run}" should be on the server
    When user closes all views
    And user opens the "bdd-bio-numbering-{run}" project and waits for its table
    Then the table should have 40 rows
    And "AntibodyHC" column should have semantic type "Macromolecule"
    And "AntibodyHC" column should have units "fasta"
    And the table should have a column "AntibodyHC (aligned)"
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    And every value of "AntibodyHC (aligned)" column should have the same length
    And "AntibodyHC (aligned)" column should have no missing values
    And the ".positionNames" tag of "AntibodyHC (aligned)" column should list at least 100 values
    And "AntibodyHC" column should carry at least 7 annotations
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    When user selects "Immunum" in Engine input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned) (2)" should have been added
    And "AntibodyHC (aligned) (2)" column should have tag ".numberingScheme" equal to "imgt"
    And "AntibodyHC (aligned) (2)" column should have no missing values
    And no error or warning balloon should have been shown
    And no errors should have been logged
