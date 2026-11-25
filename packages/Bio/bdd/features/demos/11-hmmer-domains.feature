@guide @demo
Feature: 6 · Antibodies and proteins · Pfam domains
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag
  leaves the set-up and the checks off camera). The HMM package runs HMMER in the browser: Find Domains
  scans every heavy chain of Bio's antibodies sample against the package's built-in Pfam library and
  colors the V-set and C1-set immunoglobulin domains on the chains, with a column that names them; a chain's
  context panel shows its ANARCI numbering (chain, species, closest germlines, IMGT CDRs) and its Pfam
  domains; Profile HMM Logo docks Bio's WebLogo drawing the V-set family's profile, picked on camera from
  the dialog's Family list (the dialog opens on the list's first family, 7tm_1, a GPCR family, rather than
  on V-set, which it means to default to: a product bug). Nothing is fetched from InterPro: every model
  comes with the package. Two notes: one as the chain's panes open, one as the V-set logo is added (read over
  its arrival, gone before the clip ends on it).

  Scenario: Pfam domains, germlines and profile logos
    HMMER in the browser: Pfam domains on every chain, germlines, a family's profile logo.

    Given user is logged in
    And simple mode is off
    And user opens antibodies dataset
    And the Bio package is initialized
    And the toolbox pane is hidden
    And the context panel is closed
    And the grid shows only the columns "Antibody_ID, AntibodyHC"
    And the "AntibodyHC" column of the grid is 1400 px wide
    And the header of the "AntibodyHC" column starts at position "85"
    When user picks "Bio > Annotate > Find Domains (HMMER)..." from the top menu
    Then "Find Domains (HMMER)" dialog should be visible
    When user clicks on OK button in "Find Domains (HMMER)" dialog
    Then a new column "AntibodyHC domains" should have been added
    And some value of "AntibodyHC domains" column should contain "V-set"
    And some value of "AntibodyHC domains" column should contain "C1-set"
    Given the notification balloons are closed
    And the "AntibodyHC" column of the grid is 1100 px wide
    And the "AntibodyHC domains" column of the grid is 300 px wide
    When the video says "Pfam domains of every chain: V-set, then C1-set"
    And user clicks on the "cell 3 of AntibodyHC" area of grid
    Then the "current row" reading of grid should be 3
    Given the context panel is open
    And the "AntibodyHC" column of the grid is 700 px wide
    And the header of the "AntibodyHC" column starts at position "100"
    When user expands Bioinformatics accordion header in context panel
    And user expands ANARCI accordion header in context panel
    Then "ANARCI" pane in context panel should contain the text "IGHV"
    When the video says "ANARCI and Pfam panes: germlines, CDRs and domains"
    When user expands "Pfam Domains" accordion header in context panel
    Then "Pfam Domains" pane in context panel should contain the text "V-set"
    When user picks "Bio > Analyze > Profile HMM Logo..." from the top menu
    Then "Profile HMM Logo" dialog should be visible
    When user selects "V-set" in Family input in "Profile HMM Logo" dialog
    Then Family input in "Profile HMM Logo" dialog should have value "V-set"
    When the video says "The V-set profile: its conserved Cys and Trp stand tall"
    And user clicks on OK button in "Profile HMM Logo" dialog
    Then WebLogo viewer should be visible
    And "profile" property of WebLogo viewer should be "V-set"
    And the "positions shown" reading of WebLogo viewer should be at least 20
    And no error or warning balloon should have been shown
    And no errors should have been logged
