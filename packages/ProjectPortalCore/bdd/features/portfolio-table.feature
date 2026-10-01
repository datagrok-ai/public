@portal
Feature: The portfolio as a table, with filters and analytics
  A portfolio's programs are a table: one row per program, with its stage in the model's order,
  the days it has spent there, its readiness for the next gate, its compounds by role and its work
  units. The table opens with the filter panel and a default set of viewers; the readiness matrix
  and the compound roles are tables of their own, linked to the programs, so the filters drive them
  too. The example portfolios run on synthetic registries.

  Background:
    Given user opens the Project Portal app

  Scenario: The development portfolio as a table, filtered to one stage
    When user clicks on biologics development portfolio
    And user clicks on table button
    Then the "Example Biologics Development" table view should open with 7 rows
    And the address should end with "/Portal/EX-BIOC/table"
    And the portfolio table should show the viewers "Programs by stage, Days in stage, Readiness for the next gate, Compounds by role"
    When user adds a categorical filter on "stage" keeping "D4"
    Then the filter should pass exactly the rows where "stage" is "D4"
    And the readiness matrix should show the programs that pass the filter
