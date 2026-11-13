@guide @demo
Feature: Bioinformatics · Antibody numbering
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). The heavy chains of Bio's antibodies sample are numbered
  in IMGT by the default engine, ANARCI from the HMM package (HMMER in the browser): every chain shows its
  frameworks and CDRs in color, a monomer's tooltip names its region, and the aligned column the
  numbering adds has a header of conservation, WebLogo and region tracks over a ruler of IMGT positions,
  dragged along to CDR3, whose positions and region show under the pointer. The light chains are
  numbered the same way, and the clip ends on both CDR3s. Needs the HMM package (without it Bio's own
  immunum engine would be the default; the clip skips). The aligned column starts in the N-terminal
  flank of the chains that have one; the clip shows it from IMGT position 1.

  Scenario: Antibody numbering in the browser
    ANARCI numbers 493 antibodies in IMGT: frameworks and CDRs colored on every chain, and aligned columns with WebLogo and region tracks.

    Given user is logged in
    And simple mode is off
    And the "HMM" package is installed
    And user opens antibodies dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC"
    And the "AntibodyHC" column of the grid is 1400 px wide
    When user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    Then "Apply Antibody Numbering" dialog should be visible
    And Engine input in "Apply Antibody Numbering" dialog should have value "ANARCI (HMMER)"
    And Scheme input in "Apply Antibody Numbering" dialog should have value "imgt"
    When user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyHC (aligned)" should have been added
    And "AntibodyHC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    And "AntibodyHC" column should carry at least 7 annotations
    Given the "AntibodyHC" column of the grid is 1400 px wide
    When the video says "Frameworks and CDRs, colored on every chain"
    And user hovers over monomer 54 of the "cell 3 of AntibodyHC" area of grid
    Then tooltip should contain text "CDR2"
    When the video says "The aligned column, ruled in IMGT positions"
    Given the tooltip is hidden
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned)"
    And the header of the "AntibodyHC (aligned)" column starts at position "1"
    Then the "header tracks of AntibodyHC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    Given the "AntibodyHC (aligned)" column of the grid stays 1400 px wide
    Then the "header positions of AntibodyHC (aligned)" reading of grid should be "10, 20, 30, 40, 50, 60, 70, 80, 90"
    When user drags the ruler under the "position 85 of AntibodyHC (aligned) header" area of grid by 1100 pixels to the left
    Then the "header positions of AntibodyHC (aligned)" reading of grid should be "80, 90, 100, 110, 120, C+10, C+20"
    When user hovers over the "position 110 of AntibodyHC (aligned) header" area of grid
    Then tooltip should contain text "Position: 110"
    When user hovers over the region track under the "position 111A of AntibodyHC (aligned) header" area of grid
    Then tooltip should contain text "CDR3"
    And tooltip should contain text "Position: 111A"
    When the video says "Light chains, numbered the same way"
    Given user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user selects "AntibodyLC" in Sequence input in "Apply Antibody Numbering" dialog
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    Then a new column "AntibodyLC (aligned)" should have been added
    And "AntibodyLC (aligned)" column should have tag ".numberingScheme" equal to "imgt"
    Given the notification balloons are closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned), AntibodyLC (aligned)"
    And the header of the "AntibodyHC (aligned)" column starts at position "100"
    And the header of the "AntibodyLC (aligned)" column starts at position "100"
    Then the "header tracks of AntibodyLC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    Given the "AntibodyHC (aligned)" column of the grid stays 700 px wide
    And the "AntibodyLC (aligned)" column of the grid stays 700 px wide
    Then the "header positions of AntibodyLC (aligned)" reading of grid should be "100, 110, 120"
    And no error or warning balloon should have been shown
    And no errors should have been logged
