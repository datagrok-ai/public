@journey @realizes:bio.annotate.numbering-scheme
Feature: Antibody numbering with the bundled immunum engine
  Bio | Annotate | Apply Numbering Scheme... runs a numbering engine in a worker: the dialog
  offers the engines registered and the schemes the chosen one supports, and the run leaves an
  aligned column tagged with the scheme and the position names, plus the FR and CDR regions as
  annotations on the source column. Kabat, chosen over the default IMGT, must reach the result.

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
    And Scheme input in "Apply Antibody Numbering" dialog should have value "imgt"

  Scenario: Kabat numbering aligns the column and annotates its regions
    When user selects "kabat" in Scheme input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then the top menu command should have completed
    And a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "kabat"
    And every value of "AntibodyHC (aligned)" column should have the same length
    And "AntibodyHC (aligned)" column should have no missing values
    And "AntibodyHC" column should carry at least 7 annotations
    And no error or warning balloon should have been shown
    And no errors should have been logged

  Scenario: The engine's own result honours the five-column contract
    When user calls "Bio:immunumAntibodyNumbering" function with:
      | df     | table             |
      | seqCol | column:AntibodyHC |
      | scheme | imgt              |
    Then the result should be a table with columns "position_names, chain_type, annotations_json, numbering_detail, numbering_map"
    And every column of the result table should be filled in row 1
    And no errors should have been logged
