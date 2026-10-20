@portal
Feature: Find a program from the landing page
  The top of the landing page searches every portfolio the user can see, by a program's code,
  name, target or an alias its registry delivered. A hit opens the program in its own portfolio.
  The programs the user opened last are listed under Recent, kept in their user settings.

  Background:
    Given user opens the Project Portal app

  Scenario: Find a program of another portfolio by its target, then open it again from Recent
    Given the user's recent programs are cleared
    When user enters "ALP-1" into program search input
    Then search results should contain text "Antibody Alpha"
    And search results should contain text "target: ALP-1"
    When user clicks on EXB-201 search hit
    Then board header should contain text "EX-BIOC-V1"
    And program detail should contain text "Antibody Alpha"
    And the address should end with "/EX-BIOC/EXB-201"
    When user clicks on portfolios link
    Then recent programs should contain text "Antibody Alpha"
    When user clicks on EXB-201 recent
    Then board header should contain text "EX-BIOC-V1"
    And program detail should contain text "Antibody Alpha"
