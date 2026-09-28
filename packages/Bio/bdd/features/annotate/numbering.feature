@journey @realizes:bio.annotate.numbering-scheme
Feature: Antibody numbering with the bundled immunum engine
  Bio | Annotate | Apply Numbering Scheme... runs a numbering engine in a worker: the dialog
  offers the engines registered and the schemes the chosen one supports, and the run leaves an
  aligned column tagged with the scheme and the position names, plus the FR and CDR regions as
  annotations on the source column. Kabat, chosen over the default IMGT, must reach the result.
  The project round-trip of a numbering and its deterministic re-run are in projects/round-trips.

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
    And Engine input in "Apply Antibody Numbering" dialog should have value "Immunum"
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
