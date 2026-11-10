@guide @demo @viewport:1920x1080
Feature: Bioinformatics · Peptide SAR
  A demo clip, one of the Bio showcase (features/demos, filmed with `grok-bdd guide`: the @demo tag leaves
  the set-up and the checks off camera). It starts from the CDR3 SAR of clip 12, made off camera the same
  way (antibodies sample, IMGT numbering, CDR3 cut out of the aligned column, SAR on Y at -lg, MCL at 90).
  The Logo Summary Table sums up each MCL cluster (its tooltip, a click selects its members everywhere);
  then two more Peptides viewers come from the Add Viewer gallery: Sequence Position Statistics (the pKd of
  every residue at a CDR3 position, as violins; position 2, IMGT 106, edited in its header) and Sequence
  Mutation Cliffs (the mean pKd, by the residue there, of the antibodies that differ only at a position;
  position 8, IMGT 111A, where W comes out lowest). Not shown: the Mutation Cliffs viewer's Series picker
  (one line per cluster): filmed, the guide's pick left the in-viewer column selector empty (a plain run
  picks the column).

  Scenario: Clusters, residues and cliffs in more viewers
    Each MCL cluster summed up by its CDR3 logo and pKd; more viewers break pKd down by residue and by mutation cliff.

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
    And the context panel is closed
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
    And the notification balloons are closed
    When the video says "Each MCL cluster: its CDR3 logo and pKd"
    And user hovers over the "cluster 3" area of Logo Summary Table viewer
    Then the cluster statistics in tooltip should match cluster "3" of table "antibodies"
    When user clicks on the "cluster 3" area of Logo Summary Table viewer
    Then only rows where "Cluster (MCL)" is "3" should be selected
    When user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "Position" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Sequence Position Statistics" button in "Add Viewer" dialog
    Then "Add Viewer" dialog should be absent
    And Sequence Position Statistics viewer should be added to the open tableview
    Given the viewers have settled
    When user clicks on "Edit position" icon in Sequence Position Statistics viewer
    And user types "2" at the caret
    And user presses Enter
    Then the "positions" reading of Sequence Position Statistics viewer should be "2"
    When the video says "pKd of every residue at one CDR3 position"
    And user clicks on "Add viewer" icon in toolbar
    Then "Add Viewer" dialog should be visible
    When user types "Mutation" into viewer gallery search in "Add Viewer" dialog
    And user clicks on first "Sequence Mutation Cliffs" button in "Add Viewer" dialog
    Then Sequence Mutation Cliffs viewer should be added to the open tableview
    When user selects "8" in Position input in Sequence Mutation Cliffs viewer
    Then the "position" reading of Sequence Mutation Cliffs viewer should be 8
    And the "cliff rows" reading of Sequence Mutation Cliffs viewer should be 9
    And line chart viewer in Sequence Mutation Cliffs viewer should be painted
    When the video says "Cliffs at one position: W there costs affinity"
    Then no errors should have been logged
