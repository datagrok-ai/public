Feature: Launch SAR at a raised similarity threshold
  The MCL similarity threshold 90, the value of GROK-19145, runs a fresh analysis of all 647 peptides
  and keeps the WebLogo headers and monomer selection usable.

  The cluster count a threshold leaves (1 at 10, 2 at 75, 5 or 6 at 96 on the first 200 peptides)
  and the Logo Summary Table's statistics are tested in Peptides src/tests/model.ts
  'Similarity threshold 10', '75' and '96' (category 'Model: MCL clusters').

  Not translated: the manual case's "turn the optional viewers off in the SAR dialog" — the dialog
  has no viewer toggles any more, only Generate clusters. Its "Activity" column is IC50 here.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: Similarity threshold 90 on all peptides preserves the analysis and its selection behavior
    Then the table should have 647 rows
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    And user enters "90" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And the SAR setting "mclSettings.threshold" should be "90"
    And the "completed threshold" reading of MCL viewer should be 90
    And the "members total" reading of Logo Summary Table viewer should be 647
    And the "positions" reading of Sequence Variability Map viewer should be 17
    And the "header 2" area of grid should be at least 100 pixels tall
    When user clicks on the "A at 2" area of grid
    Then 299 rows should be selected
    And only rows where "2" is "A" should be selected
    And no errors should have been logged
    And no error or warning balloon should have been shown
