Feature: Launch SAR at different similarity thresholds
  Several MCL similarity thresholds keep the WebLogo headers and monomer selections usable.
  Each threshold launches a fresh analysis of the first 200 peptides; threshold 90, the value of
  GROK-19145, also runs on all 647.

  Not translated: the manual case's "turn the optional viewers off in the SAR dialog" — the dialog
  has no viewer toggles any more, only Generate clusters. Its "Activity" column is IC50 here.

  Background:
    Given user is logged in
    And the Peptides package is initialized
    And user opens peptides dataset keeping the first 200 rows
    Then no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: Similarity threshold <threshold> preserves the analysis and its selection behavior
    When user picks "Bio > Analyze > SAR..." from the top menu
    Then "Analyze Peptides" dialog should be visible
    And "Generate clusters" checkbox in "Analyze Peptides" dialog should be checked
    When user clicks on "Adjust clustering parameters" icon in "Analyze Peptides" dialog
    Then "Similarity Threshold" input in "Analyze Peptides" dialog should be visible
    When user enters "<threshold>" into "Similarity Threshold" input in "Analyze Peptides" dialog
    Given user listens for "peptides-sar-ready" custom event
    When user clicks on OK button in "Analyze Peptides" dialog
    Then the SAR analysis should be ready
    And the SAR setting "mclSettings.threshold" should be "<threshold>"
    And Sequence Variability Map viewer should be added to the open tableview
    And Most Potent Residues viewer should be added to the open tableview
    And MCL viewer should be added to the open tableview
    And Logo Summary Table viewer should be added to the open tableview
    And the "positions" reading of Sequence Variability Map viewer should be 17
    And the "count of cell A at 2" reading of Sequence Variability Map viewer should be 59
    And the "header 2" area of grid should be at least 100 pixels tall
    And the "A at 2" area of grid should be painted
    And scatter plot viewer in MCL viewer should be painted
    And the "completed threshold" reading of MCL viewer should be <threshold>
    And no rows should be selected
    When user clicks on the "A at 2" area of grid
    Then 59 rows should be selected
    And only rows where "2" is "A" should be selected
    When user clears the row selection
    Then no rows should be selected
    When user clicks on "Invariant Map" checkbox in Sequence Variability Map viewer
    And user clicks on the "cell A at 2" area of Sequence Variability Map viewer
    Then 59 rows should be selected
    And only rows where "2" is "A" should be selected
    And the "selected monomer-positions" reading of Sequence Variability Map viewer should be "2:A"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | threshold |
      | 10        |
      | 50        |
      | 75        |
      | 90        |
      | 93        |
      | 96        |

  Scenario: Similarity threshold 90 on all peptides preserves the analysis and its selection behavior
    When user closes all views
    And user opens peptides dataset
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
