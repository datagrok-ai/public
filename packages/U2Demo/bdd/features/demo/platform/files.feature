@demo @realizes:u2.dg.file-input
Feature: File pickers

  Scenario: A typed path resolves against the file shares
    Given user opens the "Files" demo page
    Then value of file readout should have text "(none)"
    When user enters "System:DemoFiles/demog.csv" into File input
    Then value of file readout should have text "demog.csv"
