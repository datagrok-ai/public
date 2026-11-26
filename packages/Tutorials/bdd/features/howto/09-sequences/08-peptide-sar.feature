@guide @tutorial @viewport:1920x1080
Feature: How-to · Sequences
  A how-to of the Datagrok series (features/howto, filmed with `grok-bdd guide`: the @tutorial tag films it as a
  demo clip that outlines what each step presses). The Peptides package's SAR on the demo peptides (Files > Demo >
  bio > peptides.csv: aligned peptides with their IC50; the first 200 opened off camera): Bio > Analyze > SAR picks
  the sequence and activity columns itself; Scaling -lg turns the IC50 into pIC50, and the clustering threshold is
  raised to 93 so the clusters come quickly. The analysis lays out the Sequence Variability Map (each residue at
  each position, in mutation-cliffs mode: a circle's size counts the peptides in cliffs, its color the mean change of
  activity), the most potent residues, the clusters with their logos and a Selection table; the pointer on a residue
  tells its numbers, a click selects its peptides. Filmed in a 1920x1080 window: the analysis lays out six panels.
  Nothing is made on the server.

  Scenario: Analyze peptide SAR
    The SAR analysis shows which residue at which position drives activity, with clusters and their logos.

    Given user is logged in
    And simple mode is off
    And the Peptides package is initialized
    And user opens peptides dataset keeping the first 200 rows
    And the toolbox pane is hidden
    And the context panel is closed
    When the video says "200 aligned peptides with their IC50"
    # caption: Bio > Analyze > SAR
    And user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    And editor of Sequence input in "Analyze Peptides" dialog should have text "AlignedSequence"
    And editor of Activity input in "Analyze Peptides" dialog should have text "IC50"
    # caption: Sequence and activity are picked; -lg turns IC50 into pIC50
    When user selects "-lg" in Scaling input in "Analyze Peptides" dialog
    # caption: Clusters: how alike peptides must be to group
    And user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "93" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    When the video says "Bottom left: every residue at every position. Top right: the clusters"
    # caption: Circle size: peptides in cliffs; color: how activity changes
    And user hovers over the "cell A at 2" area of Sequence Variability Map viewer
    Then tooltip should be visible
    # caption: Click a residue: its peptides are selected and listed
    When user clicks on the "cell A at 2" area of Sequence Variability Map viewer
    Then some rows should be selected
    And no errors should have been logged
