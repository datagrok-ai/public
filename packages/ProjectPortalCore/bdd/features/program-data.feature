@portal
Feature: A program's governed data one click away
  Every portfolio declares its program dataset, a query with @program, and may name further
  datasets beside it. The program page opens each one as a table named after the program, with
  the viewers the portfolio declares, or the defaults for its setting: a discovery portfolio plots
  two measures against each other. The example small-molecule portfolio calls its program dataset
  the SAR table and declares the assay results beside it.

  Background:
    Given user opens the Project Portal app

  Scenario: The SAR table of a discovery project opens as a table named after it
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    Then open data button should be visible
    When user clicks on open data button
    Then the open table should have 11 rows and a name starting with "EXS-KIN1"
    And the open table should show the viewer types "Scatter plot, Histogram"

  Scenario: The assay results open beside it
    When user clicks on small molecule discovery portfolio
    And user clicks on EXS-KIN1 card
    And user clicks on open assay results button
    Then the open table should have 8 rows and a name ending with "Assay results"
    And the open table should show the viewer types "Box plot"
