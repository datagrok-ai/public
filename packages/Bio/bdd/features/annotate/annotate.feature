@journey @realizes:bio.analyze.compare-sequences @realizes:bio.annotate.scan-liabilities @realizes:bio.annotate.manage-annotations
Feature: Comparing, scanning and annotating antibody sequences
  Two macromolecule columns of one notation (heavy and light chains): Compare sequences pairs
  them into a difference column; Scan Liabilities marks the motifs its rules find and reports
  them per row or as a count; Manage Annotations lists what a column carries and drops it.

  Background:
    Given user is logged in
    And user opens antibodies dataset keeping the first 40 rows as "antibodies"
    And the Bio package is initialized
    Then "AntibodyHC" column should have semantic type "Macromolecule"
    And "AntibodyLC" column should have semantic type "Macromolecule"

  Scenario: Compare sequences pairs the two chains with the defaults
    When user picks "Bio > Analyze > Compare sequences..." from the top menu
    Then "Compare Sequences" dialog should be visible
    And "Sequence column 1" input in "Compare Sequences" dialog should have value "AntibodyHC"
    And "Sequence column 2" input in "Compare Sequences" dialog should have value "AntibodyLC"
    When user clicks on OK button in "Compare Sequences" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "AntibodyHC vs AntibodyLC" should have been added
    And "AntibodyHC vs AntibodyLC" column should have tag "cell.renderer" equal to "MacromoleculeDifference"
    And every value of "AntibodyHC vs AntibodyLC" column should be "AntibodyHC" and "AntibodyLC" of the same row joined by "#"
    And no error or warning balloon should have been shown

  Scenario: A custom name and swapped columns are honoured
    When user picks "Bio > Analyze > Compare sequences..." from the top menu
    And user selects "AntibodyLC" in "Sequence column 1" input in "Compare Sequences" dialog
    And user selects "AntibodyHC" in "Sequence column 2" input in "Compare Sequences" dialog
    And user enters "LC minus HC" into "Result column name" input in "Compare Sequences" dialog
    And user clicks on OK button in "Compare Sequences" dialog
    Then the top menu command should have completed
    And 1 new column should have been added
    And a new column "LC minus HC" should have been added
    And every value of "LC minus HC" column should be "AntibodyLC" and "AntibodyHC" of the same row joined by "#"

  Scenario: The same column twice is refused and adds nothing
    When user picks "Bio > Analyze > Compare sequences..." from the top menu
    And user selects "AntibodyHC" in "Sequence column 2" input in "Compare Sequences" dialog
    And user clicks on OK button in "Compare Sequences" dialog
    Then the top menu command should have completed
    And an error balloon containing "distinct columns" should have been shown
    And no new column should have been added

  Scenario: Scan Liabilities opens with its documented rule defaults
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    Then "Scan Sequence Liabilities" dialog should be visible
    And "Deamidation (NG)" checkbox in "Scan Sequence Liabilities" dialog should be checked
    And "Free Cysteine" checkbox in "Scan Sequence Liabilities" dialog should be unchecked
    And "Highlight in cell renderer" checkbox in "Scan Sequence Liabilities" dialog should be checked
    And "Create annotation column" checkbox in "Scan Sequence Liabilities" dialog should be checked
    And "Create summary count column" checkbox in "Scan Sequence Liabilities" dialog should be unchecked

  Scenario: The default rules write the per-row annotation column and hit real motifs
    When user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then the top menu command should have completed
    And a new column "~AntibodyHC_annotations" should have been added
    And the table should not have a column "AntibodyHC_liability_count"
    And "AntibodyHC" column should carry at least 1 annotation
    And every liability hit on "AntibodyHC" column should match its motif at the position it reports
    And no error or warning balloon should have been shown

  Scenario: Manage Annotations lists one row per annotation and drops them one by one or all
    When user picks "Bio > Annotate > Manage Annotations..." from the top menu
    Then "Manage Annotations" dialog should be visible
    And annotations list in "Manage Annotations" dialog should have as many items as "AntibodyHC" column has annotations
    When user clicks on "delete annotation" icon in first item in annotations list
    Then "AntibodyHC" column should carry one annotation fewer than before
    And annotations list in "Manage Annotations" dialog should have as many items as "AntibodyHC" column has annotations
    When user clicks on "Clear All" button in "Manage Annotations" dialog
    Then "Manage Annotations" dialog should contain text "No annotations on this column."
    And "AntibodyHC" column should carry no annotations
    When user clicks on CANCEL button in "Manage Annotations" dialog
    Then "Manage Annotations" dialog should be hidden

  Scenario: Only the oxidation rules, reported as a count column
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    And user unchecks "Deamidation (NG)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Deamidation (NS)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Deamidation (NA)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Deamidation (ND)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Deamidation (NT)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Isomerization (DG)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Isomerization (DS)" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "N-glycosylation" checkbox in "Scan Sequence Liabilities" dialog
    And user unchecks "Create annotation column" checkbox in "Scan Sequence Liabilities" dialog
    And user checks "Create summary count column" checkbox in "Scan Sequence Liabilities" dialog
    And user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then the top menu command should have completed
    And a new column "AntibodyHC_liability_count" should have been added
    And "AntibodyHC_liability_count" column should have type "int"
    And "AntibodyHC_liability_count" column should have no missing values
    And the total of "AntibodyHC_liability_count" column should be fewer than the liability hits found before
    And no error or warning balloon should have been shown
    And no errors should have been logged
