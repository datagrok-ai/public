@portal
Feature: A portfolio exported as its file
  A loaded portfolio can be read back into the portfolio file format, to diff against a file or to
  load on another instance. exportPortfolio, diffPortfolios, exportProfile and importProfile do the
  same as functions.

  Background:
    Given user opens the Project Portal app

  Scenario: The board exports the open portfolio
    When user clicks on small molecule discovery portfolio
    And user clicks on export button
    Then export dialog should be visible
    And exported portfolio file should contain text "EX-SMD-V1"
    And exported portfolio file should contain text "Hit-to-Lead"
