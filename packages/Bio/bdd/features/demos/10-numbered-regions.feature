@guide @demo
Feature: Bioinformatics · Numbered regions
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). The heavy chains of Bio's antibodies sample, numbered in
  IMGT off camera (ANARCI, the default engine with the HMM package), give an aligned column whose cell
  menu extracts any region: CDR3 becomes a column of its own that keeps the scheme's positions (its ruler
  starts at 105, insertions named 111A, 112A), and Scan Liabilities marks the deamidation, isomerization
  and oxidation motifs in it. Needs the HMM package (the clip skips without it).

  Scenario: CDR3 extracted with its IMGT positions
    A numbered antibody column gives up any region as a column of its own, in the scheme's positions, ready to scan for liabilities.

    Given user is logged in
    And simple mode is off
    And the "HMM" package is installed
    And user opens antibodies dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And the "header tracks of AntibodyHC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    Given the notification balloons are closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned)"
    And the header of the "AntibodyHC (aligned)" column starts at position "95"
    And the "AntibodyHC (aligned)" column of the grid stays 900 px wide
    Then the "header positions of AntibodyHC (aligned)" reading of grid should be "100, 110, 120, C+10"
    When user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 12 of AntibodyHC (aligned)" area of grid
    Then a new column "AntibodyHC (aligned)(CDR3)" should have been added
    And the "header tracks of AntibodyHC (aligned)(CDR3)" reading of grid should be "Conservation, WebLogo"
    Given the notification balloons are closed
    And the "AntibodyHC (aligned)" column of the grid stays 700 px wide
    And the "AntibodyHC (aligned)(CDR3)" column of the grid stays 440 px wide
    Then the "header positions of AntibodyHC (aligned)(CDR3)" reading of grid should be "105, 110"
    When the video says "CDR3 keeps its IMGT positions: 105 to 117"
    And user hovers over the "position 111A of AntibodyHC (aligned)(CDR3) header" area of grid
    Then tooltip should contain text "Position: 111A"
    When user picks "Bio > Annotate > Scan Liabilities..." from the top menu
    And user selects "AntibodyHC (aligned)(CDR3)" in Sequence input in "Scan Sequence Liabilities" dialog
    And user clicks on OK button in "Scan Sequence Liabilities" dialog
    Then a new column "~AntibodyHC (aligned)(CDR3)_annotations" should have been added
    And "AntibodyHC (aligned)(CDR3)" column should carry at least 3 annotations
    Given the notification balloons are closed
    When the video says "Liabilities marked in every CDR3"
    And user hovers over monomer 3 of the "cell 12 of AntibodyHC (aligned)(CDR3)" area of grid
    Then tooltip should contain text "Isomerization (DG)"
    And no error or warning balloon should have been shown
    And no errors should have been logged
