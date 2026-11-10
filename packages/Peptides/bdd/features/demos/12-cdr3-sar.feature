@guide @demo @viewport:1920x1080
Feature: Bioinformatics · Peptide SAR
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). Bio's antibodies sample opens from Browse; Bio numbers the heavy
  chains in IMGT (the dialog's default engine: ANARCI of the HMM package, Bio's immunum without it); the
  aligned column's cell menu cuts CDR3 out as a column of its own; Bio > Analyze > SAR... analyzes CDR3
  against Y, the binding constant in molar, scaled by -lg (pKd 3.7 to 12.4), with MCL clusters at the
  similarity threshold 90: 11 clusters of 6 to 17 related CDR3s (the default 70 lumps 428 of the 493 into
  one). The SAR clips 13 and 14 start from this state, made off camera.

  Scenario: From antibodies to a CDR3 SAR
    Number the heavy chains, cut CDR3 out and relate its residues to binding affinity: variability map, potent residues and clusters.

    Given user is logged in
    And simple mode is off
    And the Peptides package is initialized
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Bio tree node inside browse tree is expanded
    And Files---App-Data---Bio---samples tree node inside browse tree is expanded
    When user double-clicks on Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree
    Then the "antibodies" view should be current
    And "AntibodyHC" column should have semantic type "Macromolecule"
    Given the browse panel is closed
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC, Y"
    And the "AntibodyHC" column of the grid is 1300 px wide
    And user sets the format of "Y" column in the grid to "scientific"
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    When user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    Given the grid shows only the columns "Antibody_ID, AntibodyHC (aligned), Y"
    And the "AntibodyHC (aligned)" column of the grid is 1300 px wide
    And the header of the "AntibodyHC (aligned)" column starts at position "93"
    Then the "header tracks of AntibodyHC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    When the video says "CDR3, the loop that binds the antigen, becomes a column"
    And user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 1 of AntibodyHC (aligned)" area of grid
    Then a new column "AntibodyHC (aligned)(CDR3)" should have been added
    Given the grid shows only the columns "Antibody_ID, AntibodyHC (aligned), AntibodyHC (aligned)(CDR3), Y"
    And the "AntibodyHC (aligned)" column of the grid is 820 px wide
    And the "AntibodyHC (aligned)(CDR3)" column of the grid is 560 px wide
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user selects "AntibodyHC (aligned)(CDR3)" in Sequence input in "Analyze Peptides" dialog
    And the video says "Y is a Kd in molar: -lg turns it into pKd"
    And user selects "-lg" in Scaling input in "Analyze Peptides" dialog
    And user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "90" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And the SAR setting "sequenceColumnName" should be "AntibodyHC (aligned)(CDR3)"
    And the SAR setting "activityScaling" should be "-lg"
    And the "completed threshold" reading of MCL viewer should be 90
    And the "positions" reading of Sequence Variability Map viewer should be 26
    And the "clusters shown" reading of Logo Summary Table viewer should be 11
    And Sequence Variability Map viewer should be painted
    And Most Potent Residues viewer should be painted
    And Logo Summary Table viewer should be painted
    When the video says "Variability map, most potent residues, MCL clusters and their CDR3 logos"
    Then no errors should have been logged
