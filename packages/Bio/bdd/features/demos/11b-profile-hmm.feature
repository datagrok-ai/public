@guide @demo
Feature: Bioinformatics · Profile HMMs from your data
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). The HMM package builds profile HMMs in the browser
  (hmmbuild): from the selected chains of an aligned column (here the IMGT alignment of Bio's antibodies
  sample, numbered off camera by ANARCI), Build Profile HMM builds a profile, searches the whole table with
  it (an E-value per chain) and docks the profile's logo; a letter of the logo selects the chains that have
  it at that position.

  Scenario: A profile HMM built from selected chains
    A profile HMM built from a few chains searches all 493, and its logo selects the chains by residue.

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
    And the header of the "AntibodyHC (aligned)" column starts at position "1"
    And the "AntibodyHC (aligned)" column of the grid stays 1400 px wide
    When user selects the first 15 rows
    And the video says "15 chains picked to build a profile HMM"
    And user picks "Bio > Search > Build Profile HMM..." from the top menu
    Then "Build Profile HMM" dialog should be visible
    And Rows input in "Build Profile HMM" dialog should have value "Selected rows"
    When user enters "my15" into Name input in "Build Profile HMM" dialog
    And user clicks on OK button in "Build Profile HMM" dialog
    Then a new column "my15 E-value" should have been added
    And WebLogo viewer should be visible
    Given the notification balloons are closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned), my15 E-value"
    And the header of the "AntibodyHC (aligned)" column starts at position "1"
    And the "AntibodyHC (aligned)" column of the grid stays 1150 px wide
    And the "my15 E-value" column of the grid is 130 px wide
    When the video says "The profile scores every chain; its logo selects them"
    And user selects no rows
    And user hovers over the "monomer R at position 38" area of WebLogo viewer
    And user clicks on the "monomer R at position 38" area of WebLogo viewer
    Then 378 rows should be selected
    And no error or warning balloon should have been shown
    And no errors should have been logged
