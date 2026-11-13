@guide @demo @viewport:1920x1080
Feature: Bioinformatics · Peptide SAR
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). It starts from the CDR3 SAR of clip 12, made off camera the same
  way (antibodies sample, IMGT numbering, CDR3 cut out of the aligned column, SAR on Y at -lg, MCL at 90).
  The Sequence Variability Map is docked along the bottom of the view and made taller; in Mutation Cliffs
  mode a circle marks each residue that pairs of CDR3s differ in alone (Min Activity Delta 0: every such
  pair), colored by the mean pKd difference; its tooltip counts the pairs, a click selects them, and the
  context panel lists the pairs with their pKd deltas, expanded to a dialog of its own. The Invariant Map
  mode colors every residue at every position by the mean pKd of the antibodies that carry it; a click
  selects them. Needs Peptides with the scaled cliff statistics (CHANGELOG v.next): before it, the cliffs'
  Mean activity read 0.000 for a Kd in molar.

  Scenario: Mutation cliffs and the invariant map
    One residue apart, far apart in affinity: the variability map finds the mutation cliffs of CDR3, lists the exact pairs and maps every residue.

    Given user is logged in
    And simple mode is off
    And the Peptides package is initialized
    And the browse panel is open
    And Files tree node inside browse tree is expanded
    And Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---Bio tree node inside browse tree is expanded
    And Files---App-Data---Bio---samples tree node inside browse tree is expanded
    And user double-clicks on Files---App-Data---Bio---samples---antibodies.csv tree node inside browse tree
    And the "antibodies" view should be current
    And the browse panel is closed
    And the toolbox pane is hidden
    And user picks "Bio > Annotate > Apply Numbering Scheme..." from the top menu
    And user clicks on OK button in "Apply Antibody Numbering" dialog
    And a new column "AntibodyHC (aligned)" should have been added
    And the grid shows only the columns "Antibody_ID, AntibodyHC (aligned), Y"
    And the "AntibodyHC (aligned)" column of the grid is 1300 px wide
    And the header of the "AntibodyHC (aligned)" column starts at position "93"
    And the "header tracks of AntibodyHC (aligned)" reading of grid should be "Conservation, WebLogo, Annotations"
    And user picks "Annotations > Extract CDR3 as Column" from the context menu of the "cell 1 of AntibodyHC (aligned)" area of grid
    And a new column "AntibodyHC (aligned)(CDR3)" should have been added
    And user picks "Bio > Analyze > SAR..." from the top menu
    And user selects "AntibodyHC (aligned)(CDR3)" in Sequence input in "Analyze Peptides" dialog
    And user selects "-lg" in Scaling input in "Analyze Peptides" dialog
    And user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "90" into "Similarity Threshold" input in "Analyze Peptides" dialog
    And user listens for "peptides-sar-ready" custom event
    And user clicks on OK button in "Analyze Peptides" dialog
    And the SAR analysis should be ready
    And the context panel is open
    And the notification balloons are closed
    When user docks Sequence Variability Map viewer to the bottom edge of the view
    And user drags the top border of Sequence Variability Map viewer by 200 pixels up
    And user hovers over the "cell W at 8" area of Sequence Variability Map viewer
    Then tooltip should be visible
    And "Pairs count" table row in tooltip should have text "Pairs count13"
    And "Mean activity" table row in tooltip should not contain text "0.000"
    When the video says "Mutation cliffs: CDR3 pairs one substitution apart"
    And user clicks on the "cell W at 8" area of Sequence Variability Map viewer
    Then the "selected monomer-positions" reading of Sequence Variability Map viewer should be "8:W"
    And context panel should contain text "Mutation Cliffs pairs"
    When user hovers over "Mutation Cliffs pairs" pane in context panel
    And the video says "Every pair: the substitution and its jump in pKd"
    And user clicks on "Expand to full screen" icon in context panel
    Then "Mutation Cliffs pairs" dialog should be visible
    When user presses Escape
    Then "Mutation Cliffs pairs" dialog should be absent
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    Then the "mode" reading of Sequence Variability Map viewer should be "Invariant Map"
    When user hovers over the "cell G at 2" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "23 ("
    When the video says "Invariant map: every residue, colored by its mean pKd"
    And user hovers over the "cell W at 6" area of Sequence Variability Map viewer
    Then Count table row in tooltip should contain text "6 ("
    When user clicks on the "cell W at 6" area of Sequence Variability Map viewer
    Then 15 rows should be selected
    And context panel should contain text "Invariant map: 6:W"
    And context panel should contain text "Mutation cliffs: 8:W"
    When the video says "A click selects those antibodies in every viewer"
    Then no errors should have been logged
