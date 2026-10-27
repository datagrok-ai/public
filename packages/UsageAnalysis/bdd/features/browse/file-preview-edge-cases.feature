@browse @edge-cases @realizes:views.browse
Feature: File preview edge cases: a file of a format the platform does not read, an empty file
  Guards for a file Browse cannot preview. The help says: "Unsupported file formats cannot be
  previewed, but you can download them" (help/access/files/files.md, the note under File browsing
  and preview). browse-file-formats claims the formats the platform reads; here a file with an
  extension no reader claims (".qqq") is clicked: nothing fails, its details show in the context
  panel, and its menu downloads it with its content. An empty (0-byte) .csv and .txt are clicked and
  previewed with no error, and their details show in the context panel.

  The file is this feature's own, written through the files API into the package's App Data folder
  (System:AppData/UsageAnalysis), named with the run's time, and deleted when the feature ends.

  Background:
    Given user is logged in
    And the browse panel is open
    And the context panel is open
    And a file "System:AppData/UsageAnalysis/bdd-edge-odd-{time}.qqq" with text "an unreadable format" is on the server
    And a file "System:AppData/UsageAnalysis/bdd-edge-empty-{time}.csv" with text "" is on the server
    And a file "System:AppData/UsageAnalysis/bdd-edge-empty-{time}.txt" with text "" is on the server
    And Files tree node inside browse tree is expanded
    When user refreshes the browse tree
    Given Files---App-Data tree node inside browse tree is expanded
    And Files---App-Data---UsageAnalysis tree node inside browse tree is expanded
    Then Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree should be visible

  Scenario: A file of an unsupported format is clicked with no error and shows its details
    When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree
    Then the context panel should show "bdd-edge-odd-{time}.qqq"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario: A file of an unsupported format is downloaded from its menu
    Given user watches downloads
    When user opens the context menu of Files---App-Data---UsageAnalysis---bdd-edge-odd-{time}.qqq tree node inside browse tree
    And user downloads a file through "Download" menu item in context menu
    Then a file "bdd-edge-odd-{time}.qqq" should have been downloaded
    And the downloaded file "bdd-edge-odd-{time}.qqq" should contain text "an unreadable format"
    And no errors should have been logged
    And no error or warning balloon should have been shown

  Scenario Outline: An empty <ext> file is previewed with no error
    When user clicks on Files---App-Data---UsageAnalysis---bdd-edge-empty-{time}.<ext> tree node inside browse tree
    Then the context panel should show "bdd-edge-empty-{time}.<ext>"
    And no errors should have been logged
    And no error or warning balloon should have been shown

    Examples:
      | ext |
      | csv |
      | txt |
